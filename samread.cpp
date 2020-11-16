#include "GSam.h"
#include "GArgs.h"
#include "GStr.h"
#include "GHashMap.hh"

const char* USAGE="Usage:\n samread [--sam|--S|--bam|-B|--fasta|-F|--fastq|-Q|--gff|-G] \n\
   [--ref|-r <ref.fa>] [-A|--all] [--table|-T] [-Y] [--nstrand] \n\
   [-o <outfile>] <in.bam>|<in.sam> ..\n";
/*
		"
 Recognized fields for the --table output option:\n\
  SAM columns: @qname, @flag, @rname, @pos, @mapg, @cigar, @rnext, @pnext,\n\
               @tlen ,@seq, @qual, @aux\n\
  SAM tag names: use the 2-letter string of the tag name directly; its value\n\
                 will be shown instead (or . if tag is not found)\n\
  Special constructs: @end   :the coordinate of the right-most aligned base\n\
                      @exons :comma delimited list of exon segments\n\
                      @strand :+/- alignment strand from SAM flag data\n\
                      @tstrand :+/-/. transcription strand as found in\n\
                                'XS' or 'ts' SAM tags\n\
";
*/
enum OutType {
  outFASTQ,
  outFASTA,
  outGFF,
  outTable,
  outSAM
};

OutType out_type=outFASTQ;
bool all_reads=false; //including unmapped
bool addYC=false;
bool nstrand=false; //--nstrand enforce '.' strand for unspliced alignments
GSamWriter* samwriter=NULL;

GHash<int> rnames;
int last_refid=-1;


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

int getAlnId(GSamRecord& rec) {
	if (rec.refId()!=last_refid) {
		last_refid=rec.refId();
		rnames.Clear();
	}
	if (rec.isPrimary())
		return 1;
    int* c=rnames.Find(rec.name());
    if (c) {
    	(*c)++;
    	return *c;
    }
    rnames.Add(rec.name(), 2);
    return 2;
}

void showgff(GSamRecord& rec, FILE* fout) {
  if (rec.isUnmapped()) return;
  char tstrand=rec.spliceStrand();
  char alnstrand=rec.alnStrand();
  int alnid=getAlnId(rec);
  fprintf(fout, "%s\tbam\tmRNA\t%d\t%d\t.\t%c\t.\tID=%s.%d;aln_strand=%c;num_exons=%d\n", rec.refName(),
         rec.start, rec.end, tstrand, rec.name(), alnid, alnstrand, rec.exons.Count());
  for (int i=0;i<rec.exons.Count();i++) {
    fprintf(fout, "%s\tbam\texon\t%d\t%d\t.\t%c\t.\tParent=%s.%d\n", rec.refName(),
           rec.exons[i].start, rec.exons[i].end, tstrand, rec.name(), alnid);
     }
}

void showTable(GSamRecord& rec, FILE* fout) {
	static const char* dot=".";
	if (rec.isUnmapped()) return;
	char tstrand=rec.spliceStrand();
	const char* md=rec.tag_str("MD");
	if (md==NULL) md=dot;
	// readName, refName, start, tstrand, cigar, exons, MD
	GStr exons;
	for (int i=0;i<rec.exons.Count();i++) {
		exons+=rec.exons[i].start;exons+='-';
		exons+=rec.exons[i].end;
		if (i+1<rec.exons.Count()) exons+=',';
	}
	if (nstrand && rec.exons.Count()==1) tstrand='.';
	fprintf(fout, "%s\t%s\t%d\t%c\t%s\t%s\t%s",rec.name(), rec.refName(),
			rec.start, tstrand, rec.cigar(), exons.chars(), md);
	if (addYC) {
		int v=rec.tag_int("YC");
		if (v==0) v=1;
		fprintf(fout, "\t%d", v);
	}
    fprintf(fout, "\n");
}

void showSAM(GSamRecord& rec) {
  if (!all_reads && rec.isUnmapped()) return;
   samwriter->write(&rec);
}

int main(int argc, char *argv[])  {
    GArgs args(argc, argv, "fasta;fastq;sam;nstrand;bam;gff;all;table;help;ref="
        "hBAFTSGYaqo:r:");
    args.printError(USAGE, true);
    bool outBAM=false;
    if (args.getOpt('h') || args.getOpt("help") || args.startNonOpt()==0) {
      GMessage(USAGE);
      return 1;
    }
    //args.printCmdLine(stderr);
    all_reads=(args.getOpt('A') || args.getOpt("all"));
    addYC=args.getOpt('Y');
    nstrand=args.getOpt("nstrand");
    //mapped_only=(args.getOpt('M') || args.getOpt("mapped-only"));

    if (args.getOpt('F') || args.getOpt("fasta"))
       out_type=outFASTA;
    else if (args.getOpt('G') || args.getOpt("gff")) {
       out_type=outGFF;
       all_reads=false;
    }
    else if (args.getOpt('T') || args.getOpt("table")) {
       out_type=outTable;
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


    const char* fname=NULL;
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
	}

    while ((fname=args.nextNonOpt())) {
    	GStr infn(fname);
    	if (infn!="-")
    	  if (fileExists(fname)<2) {
    		GError("Error: %s is not a valid file\n",fname);
    	  }
    }

	bool writerCreated=false;
	args.startNonOpt(); //start parsing again the non-option arguments
    while ((fname=args.nextNonOpt())) {
		GSamReader samreader(fname, SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX,
				cram_ref);
		if (out_type==outSAM && !writerCreated) {
		   GSamFileType st=outBAM ? GSamFile_BAM : GSamFile_SAM;
		   samwriter=new GSamWriter(outfname, samreader.header(), st);
		   writerCreated=true;
		}

		GSamRecord aln;
		if (out_type==outFASTA) {
			while (samreader.next(aln)) {
			   showfasta(aln, fout);
			}
		}
		else if (out_type==outGFF) {
			while (samreader.next(aln)) {
			   showgff(aln, fout);
			}
		}
		else if (out_type==outTable) {
			while (samreader.next(aln)) {
			   showTable(aln, fout);
			}
		}
		else if (out_type==outSAM) {
			while (samreader.next(aln)) {
			   showSAM(aln);
			}
		}
		else { //default: FASTQ output
			while (samreader.next(aln)) {
			  showfastq(aln, fout);
			}
		}
    }
	if (out_type==outSAM) {
		delete samwriter;
	} else  if (fout!=stdout) fclose(fout);
    return 0;
}
