#include "GSam.h"
#include "GArgs.h"
#include "GStr.h"
#include "GHashMap.hh"

/* const char* USAGE="Usage:\n samread [--bundle] [--sam|-S|-B|-U|--fasta|-F|--fastq|-Q|--gff|-G]\n\
   [--ref|-r <ref.fa>] [--best] [-A|--all] [--table|-T] [-Y] [--human-rat] [--nstrand] \n\   
   [-o <outfile>] <in.bam>|<in.sam> ..\n"; */

const char* USAGE = R"TXT(
Usage:
 samread [--sam|-S|--bam|-B|--fasta|-F|--fastq|-Q|--gff|-G] [--ref|-r <ref.fa>]
   [--best] [-A|--all] [--table|-T] [-Y] [--human-rat] [--nstrand]
   [-o <outfile>] <in.bam>|<in.sam> ..

Options:
  --sam,-S          Output alignments in SAM format
  --bam,-B          Output alignments in BAM format
  --fasta,-F        Output sequences in FASTA format
  --fastq,-Q        Output sequences in FASTQ format
  --gff,-G          Output annotations in GFF format
  --ref,-r <ref.fa> Use the specified FASTA file as reference for CRAM input
  --best            outputs only the 'best' alignments for each read that have
                    the highest alignment score (AS tag). Also adds adds a YH
                    tag to each output alignment indicating the number
                    of 'best' alignments (only if it differs from NH tag value)
  -A|--all          include unmapped reads in the output (ignored in some formats)
  --table|-T        Output alignments in a tab-delimited table format
  -Y                add the TieBrush YC tag
  --human-rat       collect per-genome alignment summary for dual genome human+rat
                    alignments (human iPSC cultures on rat astrocytes)
  --nstrand         enforce neutral ('.') strand for unspliced alignments in 
                    table output (override the transcription strand information)
                    processed or filtered based on strand-specific information
  -o <outfile>      specify the output file, instead of default stdout

Description:
  'samread' is a utility for processing and converting alignment files. It supports
  multiple input and output formats, including SAM, BAM, FASTA, FASTQ, and GFF. Users
  can filter alignments, extract sequences, and perform format conversions. The tool
  is designed to work with genomic data, offering various options to customize the
  output according to research needs.

Examples:
  samread --bam -o filtered.bam --best in.bam
    Filters the input 'in.bam' to output only the best alignments for each read into
    'filtered.bam' in BAM format.

  samread --fasta -A -o sequences.fasta in.bam
    Extracts all sequences from 'in.bam' and outputs them in FASTA format to
    'sequences.fasta', including non-primary alignments.

)TXT";

/*
 
 TODO: Recognized fields for the --table output option:\n\
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
  outSAM,
  outHumanRat // alignments vs dual-genome, output per-genome stats
};

OutType out_type=outFASTQ;
bool all_reads=false; //including unmapped
bool addYC=false;
bool nstrand=false; //--nstrand enforce '.' strand for unspliced alignments
bool bestAlns=false; // --best, only keep highest score

GSamWriter* samwriter=NULL;

GHash<int> rnames;
int last_refid=-1;
bool headerPrinted=false;

//aln stats per file:
struct TAlnStats { // [2] : per mate in paired reads
   int totalReads=0; //total reads seen
   int totalAlignments=0; //total alignments found
   int numAligned[2]={0}; //how many reads were aligned at least once
   int numAlignedPairs=0; //how many aligned pairs (no matter how)
   int unpaired=0; //how many unpaired reads were found
   //^6
   int unaligned[2]={0};
   int uniqAligned[2]={0}; //how many reads aligned exactly once
   int uniqAlignedPairs=0;
   int multiMapped[2]={0}; //how many reads have at least 2 mappings (NH>1)
   int multiMappedPairs=0; //how many pairs have at least 2 mappings (NH>1)
   //^8
   int haveSecAln[2]={0}; //how many reads have seconday alignments
   int mmover5=0; //how many pairs have NH>5
   int mmover10=0; //how many pairs have NH>10
   int mmover20=0; //how many pairs have NH>20
   int mmover40=0; //how many pairs have NH>40
   int maxNH=0;
   //^7
   int concPairs=0; //concordantly aligned pairs
   int discPairs=0; //discordantly aligned pairs
   int discPairsHuman=0; //discordant pairs, both reads on human genome
   int discPairsRat=0; //discordant pairs, both reads on rat genome
   int discPairsHumanRat=0; //discordant pairs, split between genomes
       // -- human-rat stats:
   int humanOnly[2]={0}; //how many reads have no mappings reported on rat at all
   int humanOnlyMM[2]={0}; //of these, how many are multi-mappings within the human genome
   int ratOnly[2]={0}; //how many reads have no mappings reported on human at all
   int ratOnlyMM[2]={0}; //of these, how many are multi-mappings within the rat genome
   int humanOnlyPairs=0; //same as above, but for the whole pair
   int humanOnlyPairsMM=0; // multimapped human-only pairs
   int ratOnlyPairs=0;
   int ratOnlyPairsMM=0; // multimapped rat-only pairs
   //^17
   int HumanRat[2]={0}; //how many reads have primary human and rat alignments with same score
   int HumanRatPairs=0; //how many pairs have primary human and rat alignments with same score
   int betterHuman[2]={0}; //how many reads map to both human and rat, but higher score in human
   int betterRat[2]={0}; //how many reads map to both rat and human, but higher score in rat
   int betterHumanPairs=0; //same as above but as pairs
   int betterRatPairs=0;
   //^9
   static void header(FILE* fout) {
	fprintf(fout,"file\ttotalReads\ttotalAlignments\tnumAligned_1\tnumAligned_2\tunpaired\tpairsAligned\t"
"unaligned_1\tunaligned_2\tuniqAligned_1\tuniqAligned_2\tuniqAlignedPairs\t"
"multiMapped_1\tmultiMapped_2\tmultiMappedPairs\thaveSecAln_1\thaveSecAln_2\t"
"mm_over5\tmm_over10\tmm_over20\tmm_over40\tmaxNH\tconcordantPairs\tdiscordantPairs\tdiscPairsHuman\tdiscPairsRat\t"
"discPairsHumanRat\thumanOnly_1\thumanOnlyMM_1\thumanOnly_2\thumanOnlyMM_2\thumanOnlyPairs\thumanOnlyPairsMM\t"
"ratOnly_1\tratOnlyMM_1\tratOnly_2\tratOnlyMM_2\tratOnlyPairs\tratOnlyPairsMM\tHumanRat_1\tHumanRat_2\tHumanRatPairs\t"
"betterHuman_1\tbetterHuman_2\tbetterHumanPairs\tbetterRat_1\tbetterRat_2\tbetterRatPairs\n");
   }
   void report(FILE* fout, const char* fname) {
	fprintf(fout,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t" //1+14+7
  "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", //26
	fname, totalReads, totalAlignments, numAligned[0], numAligned[1], unpaired, numAlignedPairs, unaligned[0], unaligned[1],
	uniqAligned[0], uniqAligned[1], uniqAlignedPairs, multiMapped[0], multiMapped[1], multiMappedPairs,
	haveSecAln[0], haveSecAln[1], mmover5, mmover10, mmover20, mmover40, maxNH, concPairs, discPairs,
	discPairsHuman, discPairsRat, discPairsHumanRat, humanOnly[0],humanOnlyMM[0], humanOnly[1], humanOnlyMM[1],
	humanOnlyPairs, humanOnlyPairsMM, ratOnly[0], ratOnlyMM[0], ratOnly[1], ratOnlyMM[1], ratOnlyPairs, ratOnlyPairsMM,
	HumanRat[0], HumanRat[1], HumanRatPairs, betterHuman[0], betterHuman[1], betterHumanPairs,
	betterRat[0], betterRat[1], betterRatPairs );
   }
};

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
	/*
	if (rec.refId()!=last_refid) {
		last_refid=rec.refId();
		rnames.Clear();
	}
	if (rec.isPrimary())
	   return 1;
	*/
    int* c=rnames.Find(rec.name());
    if (c) {
    	(*c)++;
    	return *c;
    }
    rnames.Add(rec.name(), 1);
    return 1;
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

struct TPairData {
	bool hasHuman[2]={false,false};
	int numaln[2]={0,0}; //NH for each mate
	bool hasRat[2]={false, false};
	bool hasSec[2]={false, false};
	int maxNH[2]={0,0};
	int maxNHPair=0;
	bool pair_aligned=false;
	int maxHscore[2]={INT_MIN, INT_MIN};
	int maxRscore[2]={INT_MIN, INT_MIN};
	bool newrec[2]={true,true};
	bool discordant=false;
};


void flushPairData(TPairData& rdata, TAlnStats& stats) {
	 int both[2]={0,0}; // 0 =  Human xor Rat, 1 = both Human & Rat, same score
	                    // 2 = human > rat, -1 = rat > human
	 for (int m=0;m<2;m++) {
	    if (!rdata.hasHuman[m]) {
	    	stats.ratOnly[m]++;
	    	if (rdata.maxNH[m]>1) stats.ratOnlyMM[m]++;
	    }
	    if (!rdata.hasRat[m]) {
	    	stats.humanOnly[m]++;
	    	if (rdata.maxNH[m]>1) stats.humanOnlyMM[m]++;
	    }
	    if (rdata.hasHuman[m] && rdata.hasRat[m]) {
	    	both[m]=1;
	    	if (rdata.maxHscore[m]==rdata.maxRscore[m]) {
	    		stats.HumanRat[m]++;
	    	}
	    	else { // scores are not equal
	    		if (rdata.maxHscore[m]>rdata.maxRscore[m]) {
	    		  both[m]=2;
	    		  stats.betterHuman[m]++;
	    	    } else { //rat score is higher
	    	      both[m]=3;
	    		  stats.betterRat[m]++;
	    	    }
	    	}
	    }
	 }
	 rdata.maxNHPair=GMAX(rdata.maxNH[0], rdata.maxNH[1]);
	 if (both[0]==both[1]) {
		if (both[0]==0) {
		 				if (rdata.hasHuman[0]) {
		 					stats.humanOnlyPairs++;
		 					if (rdata.maxNHPair>1) stats.humanOnlyPairsMM++;
		 				} else {
		 					stats.ratOnlyPairs++;
		 					if (rdata.maxNHPair>1) stats.ratOnlyPairsMM++;
		 				}
		} else if (both[0]==1) {
			stats.HumanRatPairs++;
		} else if (both[0]==2) stats.betterHumanPairs++;
		                  else stats.betterRatPairs++; //both = -1

	 }
	 if (rdata.maxNH[0]==1) stats.uniqAligned[0]++;
	 else if (rdata.maxNH[0]>1)  stats.multiMapped[0]++;
	 if (rdata.maxNH[1]==1) stats.uniqAligned[1]++;
	 else if (rdata.maxNH[1]>1)  stats.multiMapped[1]++;
	 if (rdata.pair_aligned) {
	     if (rdata.maxNHPair==1) stats.uniqAlignedPairs++;
	       else if (rdata.maxNHPair>1) stats.multiMappedPairs++;
	 }
	 if (rdata.maxNHPair>stats.maxNH) stats.maxNH=rdata.maxNHPair;
	 if (rdata.maxNHPair>5) stats.mmover5++;
	 if (rdata.maxNHPair>10) stats.mmover10++;
	 if (rdata.maxNHPair>20) stats.mmover20++;
	 if (rdata.maxNHPair>40) stats.mmover40++;
	 rdata={};
}

void statsHumanRat(GSamReader& samreader, FILE* fout) {
  GSamRecord rec;
  TAlnStats stats;
  kstring_t last_qname=KS_INITIALIZE;
  ks_resize(&last_qname, 80);
  TPairData rdata={};
  while (samreader.next(rec)) {
	 const char* qname=rec.name();
	 if (strcmp(qname, ks_c_str(&last_qname))!=0) {
		 //--flush rdata into stats
		 flushPairData(rdata, stats);
		 ks_clear(&last_qname);
		 kputs(qname, &last_qname);
	 }
	 int mate=rec.pairOrder();
	 if (mate>0) mate--;
	 else {
		 if (rdata.newrec[0]) {
		   stats.unpaired++;
		 }
	 }
     if (rec.isUnmapped()) {
		 stats.totalReads++;
    	 stats.unaligned[mate]++;
    	 continue;
     }
     int nh=rec.tag_int("NH", 1);
     if (nh>rdata.maxNH[mate]) rdata.maxNH[mate]=nh;
     int score=rec.tag_int("AS", 100);
     if (!rec.isPrimary()) rdata.hasSec[mate]=true;
     const char* yt=rec.tag_str("YT");
     //pairing: UU=not in a pair, CP=concordantly aligned pair, DP=discordantly aligned pair, UP=pair unaligned
     if (yt==NULL) yt=".";

     const char* refname=rec.refName();
     if (startsWith(refname, "rat_")) {
    	 rdata.hasRat[mate]=true;
    	 if (rdata.maxRscore[mate]<score)
    		 rdata.maxRscore[mate]=score;
     } else { //human mapping
    	 rdata.hasHuman[mate]=true;
    	 if (rdata.maxHscore[mate]<score)
    		 rdata.maxHscore[mate]=score;
     }
     rdata.numaln[mate]++;
     stats.totalAlignments++;
     if (rdata.numaln[0]==1 && rdata.numaln[1]==1) {
    	//this is the first paired alignment for this read pair
    	rdata.pair_aligned=true;
	    stats.numAlignedPairs++;
        if (strcmp(yt, "CP")==0) stats.concPairs++;
        else if (strcmp(yt, "DP")==0) {
        	stats.discPairs++;
        	rdata.discordant=true;
        	if (rdata.hasHuman[0] && rdata.hasHuman[1]) {
        		stats.discPairsHuman++;
        	} else if (rdata.hasRat[0] && rdata.hasRat[1]) {
        		stats.discPairsRat++;
        	} else stats.discPairsHumanRat++;
        }
     }
	 if (rdata.newrec[mate]) {
		 stats.totalReads++;
		 stats.numAligned[mate]++;
		 rdata.newrec[mate]=false;
	 }

  }
  flushPairData(rdata, stats);
  ks_free(&last_qname);
  stats.report(fout, samreader.fileName());
}

void showTable(GSamRecord& rec, FILE* fout) {
	static const char* dot=".";
	if (rec.isUnmapped()) return;
	char tstrand=rec.spliceStrand();
	char alnstrand=rec.alnStrand();
	int isPrimary=rec.isPrimary() ? 1 : 0;
	const char* md=rec.tag_str("MD");
	if (md==NULL) md=dot;
	int as=rec.tag_int("AS",0);
	int nh=rec.tag_int("NH", -1);
	//if (nh==NULL) nh=dot;
	int nm=rec.tag_int("NM", -1);
	const char* yt=rec.tag_str("YT"); //pairing: UU=not in a pair, CP=concordantly aligned pair, DP=discordantly aligned pair, UP=pair unaligned
	if (yt==NULL) yt=dot;
	if (!headerPrinted) {
	  fprintf(fout, "qname\trefname\tstart\talnstrand\ttstrand\tprim\tcigar\texons\tMD\tAS\tNH\tNM\tYT");
	  if (addYC) fprintf(fout,"\tYC");
	  headerPrinted=true;
	}
	fprintf(fout, "\n");
	GStr exons;
	for (int i=0;i<rec.exons.Count();i++) {
		exons+=rec.exons[i].start;exons+='-';
		exons+=rec.exons[i].end;
		if (i+1<rec.exons.Count()) exons+=',';
	}
	if (nstrand && rec.exons.Count()==1) tstrand='.';
	int mate=rec.pairOrder();
	if (mate>0) fprintf(fout, "%s/%d", rec.name(), mate);
	       else fprintf(fout, "%s", rec.name());
	const char* cigar=rec.cigar();
	fprintf(fout, "\t%s\t%d\t%c\t%c\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s", rec.refName(),
			rec.start, alnstrand, tstrand, isPrimary, cigar, exons.chars(), md, as, nh, nm, yt);
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
    GArgs args(argc, argv, "fasta;fastq;sam;nstrand;bam;gff;all;table;human-rat;help;best;ref="
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
    } else if (args.getOpt("human-rat")) {
    	out_type=outHumanRat;
    }
    if (args.getOpt("best")) {
      bestAlns=true;
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
	headerPrinted=false;
	if (out_type==outHumanRat)
	    TAlnStats::header(fout);
    while ((fname=args.nextNonOpt())) {
		GSamReader samreader(fname, cram_ref,
				SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_CIGAR|SAM_AUX);
		if (out_type==outSAM && !writerCreated) {
		   GSamFileType st=outBAM ? GSamFile_BAM : GSamFile_SAM;
		   samwriter=new GSamWriter(outfname, samreader.header(), st);
		   writerCreated=true;
		}
        //bool isSorted=samreader.isSortedByCoordinate();

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
		else if (out_type==outHumanRat) {
			statsHumanRat(samreader, fout);
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
