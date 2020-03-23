#include "GSam.h"
#include "GArgs.h"
#include "GStr.h"

const char* USAGE="Usage:\n  bamread [--sam|--S|--bam|-B|--fasta|-F|--fastq|-Q|--gff|-G] \\\n"
"	[--ref|-r <ref.fa>] [-A|--all] \\\n"
"	[-o <outfile>] <in.bam>|<in.sam>\n";

enum OutType {
  outFASTQ,
  outFASTA,
  outGFF,
  outSAM,
};

OutType out_type=outFASTQ;
bool all_reads=false; //including unmapped
GSamWriter* samwriter=NULL;

void showfastq(GSamRecord& rec, FILE* fout) {
  if (rec.isUnmapped() && !all_reads) return;
  char* qseq=rec.sequence();
  fprintf(fout, "@%s\n%s\n", rec.name(), qseq);
  GFREE(qseq);
  qseq=rec.qualities();
  fprintf(fout, "+\n%s\n",qseq);
  GFREE(qseq);
}

void showfasta(GSamRecord& rec, FILE* fout) {
  if (!all_reads && rec.isUnmapped()) return;
  char* qseq=rec.sequence();
  char* alt_name=rec.tag_str("XA");
  if (alt_name)
	fprintf(fout, ">%s alt_name=%s\n%s\n", rec.name(), alt_name, qseq);
  else
    fprintf(fout, ">%s\n%s\n", rec.name(), qseq);
  GFREE(qseq);
}

void showgff(GSamRecord& rec, FILE* fout) {
  if (rec.isUnmapped()) return;
  char tstrand=rec.spliceStrand();
  fprintf(fout, "%s\tbam\tmRNA\t%d\t%d\t.\t%c\t.\tID=%s\n", rec.refName(),
         rec.start, rec.end, tstrand, rec.name());
  for (int i=0;i<rec.exons.Count();i++) {
    fprintf(fout, "%s\tbam\texon\t%d\t%d\t.\t%c\t.\tParent=%s\n", rec.refName(),
           rec.exons[i].start, rec.exons[i].end, tstrand, rec.name());
     }
}

void showSAM(GSamRecord* rec) {
  if (!all_reads && rec->isUnmapped()) return;
   samwriter->write(rec);
}

int main(int argc, char *argv[])  {
    GArgs args(argc, argv, "fasta;fastq;sam;bam;gff;all;help;ref="
        "hBAFSGaqo:r:");
    args.printError(USAGE, true);
    bool outBAM=false;
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
      GMessage(USAGE);
      return 1;
      }
    //args.printCmdLine(stderr);

    all_reads=(args.getOpt('A') || args.getOpt("all"));

    //mapped_only=(args.getOpt('M') || args.getOpt("mapped-only"));

    if (args.getOpt('F') || args.getOpt("fasta"))
       out_type=outFASTA;
    else if (args.getOpt('G') || args.getOpt("gff")) {
       out_type=outGFF;
       all_reads=false;
       }
    else if (args.getOpt('S') || args.getOpt("sam") ||
    		args.getOpt('B') || args.getOpt("bam")) {
        out_type=outSAM;
        if (args.getOpt('B') || args.getOpt("bam"))
        	outBAM=true;
    }
    char* cram_ref=NULL;
    cram_ref=args.getOpt('r');
    if (cram_ref==NULL) cram_ref=args.getOpt("ref");
    char* fname=args.nextNonOpt();
    if (fname==NULL || fname[0]==0) {
        GMessage(USAGE);
        return 1;
    }
    GSamReader samreader(fname, SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX,
    		cram_ref);
    //GSamReader samreader(fname);
    FILE* fout=stdout;
    const char* outfname=args.getOpt('o');
    if (outfname) {
     if (out_type != outSAM) {
       fout=fopen(outfname, "w");
       if (fout==NULL) {
           fprintf(stderr, "Error creating output file %s\n", outfname);
           return 2;
       }
     }
    }
    if (out_type==outSAM) {
    	if (outfname==NULL) outfname="-";
    	GSamFileType st=outBAM ? GSamFile_BAM : GSamFile_SAM;
    	samwriter=new GSamWriter(outfname, samreader.header(), st);
    }

    GSamRecord *aln=NULL;
    if (out_type==outFASTA) {
        while ((aln=samreader.next())!=NULL) {
           showfasta(*aln, fout);
           delete aln;
        }
    }
    else if (out_type==outGFF) {
        while ((aln=samreader.next())!=NULL) {
           showgff(*aln, fout);
           delete aln;
        }
    }
    else if (out_type==outSAM) {
        while ((aln=samreader.next())!=NULL) {
           showSAM(aln);
           delete aln;
        }
    }
    else { //default: FASTQ output
        while ((aln=samreader.next())!=NULL) {
          showfastq(*aln, fout);
          delete aln;
        }
    }
    if (out_type==outSAM) {
    	delete samwriter;
    } else  if (fout!=stdout) fclose(fout);
    return 0;
}
