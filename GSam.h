#ifndef _G_SAM_H
#define _G_SAM_H
#include "GBase.h"
#include "GList.hh"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/cram.h"

class GSamReader;
class GSamWriter;

enum GSamFileType {
   GSamFile_SAM=1,
   GSamFile_UBAM,
   GSamFile_BAM,
   GSamFile_CRAM
};

class GSamRecord: public GSeg {
   friend class GSamReader;
   friend class GSamWriter;
   bam1_t* b;
   // b->data has the following strings concatenated:
   //  qname (including the terminal \0)
   //  +cigar (each event encoded on 32 bits)
   //   +seq  (4bit-encoded)
   //    +qual
   //     +aux

   union {
	  uint16_t iflags;
      struct {
    	  bool novel         :1; //if set, the destructor must free b
    	  bool hard_Clipped  :1;
    	  bool soft_Clipped  :1;
    	  bool has_Introns   :1;
      };
   };
   sam_hdr_t* b_hdr;
 public:
   GVec<GSeg> exons; //coordinates will be 1-based
   int clipL; //soft clipping data, as seen in the CIGAR string
   int clipR;
   int mapped_len; //sum of exon lengths
   bool isHardClipped() { return hard_Clipped; }
   bool isSoftClipped() { return soft_Clipped; }
   bool hasIntrons() { return has_Introns; }
   //created from a reader:
   void bfree_on_delete(bool b_free=true) { novel=b_free; }
   GSamRecord(bam1_t* from_b=NULL, sam_hdr_t* b_header=NULL, bool b_free=true):b(NULL),
		   iflags(0), b_hdr(b_header), exons(1),  clipL(0), clipR(0), mapped_len(0) {
      if (from_b==NULL) {
           b=bam_init1();
           novel=true;
      }
      else {
           b=from_b; //it'll take over from_b
           novel=b_free;
      }

      b_hdr=b_header;
      setupCoordinates();//set 1-based coordinates (start, end and exons array)
   }

   //deep copy constructor:
   GSamRecord(GSamRecord& r):GSeg(r.start, r.end), iflags(r.iflags), b_hdr(r.b_hdr),
		   exons(r.exons), clipL(r.clipL), clipR(r.clipR), mapped_len(r.mapped_len) {
	      //makes a new copy of the bam1_t record etc.
	      b=bam_dup1(r.b);
	      novel=true; //will also free b when destroyed
   }

   const GSamRecord& operator=(GSamRecord& r) {
	  //copy operator
      //makes a new copy of the bam1_t record etc.
      clear();
      b=bam_dup1(r.b);
      iflags=r.iflags;
      novel=true; //will also free b when destroyed
      start=r.start;
      end=r.end;
      exons = r.exons;
      clipL = r.clipL;
      clipR = r.clipR;
      mapped_len=r.mapped_len;
      return *this;
      }

     void setupCoordinates();
     void clear() {
        if (novel) {
           bam_destroy1(b);
           //novel=false;
        }
        b=NULL;
        exons.Clear();
        mapped_len=0;
        b_hdr=NULL;
        iflags=0;
    }

    ~GSamRecord() {
       clear();
    }

    void parse_error(const char* s) {
      GError("SAM parsing error: %s\n", s);
    }

    bam1_t* get_b() { return b; }

    void set_mdata(int32_t mtid, int32_t m0pos, //0-based coordinate, -1 if not available
                     int32_t isize=0) { //mate info for current record
      b->core.mtid=mtid;
      b->core.mpos=m0pos; // should be -1 if '*'
      b->core.isize=isize; //should be 0 if not available
    }

    void set_flags(uint16_t samflags) {
      b->core.flag=samflags;
    }

    //creates a new record from 1-based alignment coordinate
    //quals should be given as Phred33
    //Warning: pos and mate_pos must be given 1-based!
    GSamRecord(const char* qname, int32_t gseq_tid,
                    int pos, bool reverse, const char* qseq, const char* cigar=NULL, const char* quals=NULL);
    GSamRecord(const char* qname, int32_t samflags, int32_t g_tid,
             int pos, int map_qual, const char* cigar, int32_t mg_tid, int mate_pos,
             int insert_size, const char* qseq, const char* quals=NULL,
             GVec<char*>* aux_strings=NULL);
             //const std::vector<std::string>* aux_strings=NULL);
    void set_cigar(const char* str); //converts and adds CIGAR string given in plain SAM text format
    void add_sequence(const char* qseq, int slen=-1); //adds the DNA sequence given in plain text format
    void add_quals(const char* quals); //quality values string in Phred33 format
    void add_aux(const char* str); //adds one aux field in plain SAM text format (e.g. "NM:i:1")
    int  add_aux(const char tag[2], char atype, int len, uint8_t *data) {
      //IMPORTANT:  strings (Z,H) should include the terminal \0
     return bam_aux_append(b, tag, atype, len, data);
    }

    int add_tag(const char tag[2], char atype, int len, uint8_t *data) {
      //same with add_aux()
      //IMPORTANT:  strings type (Z,H) should include the terminal \0
      return bam_aux_append(b, tag, atype, len, data);
    }

    int add_int_tag(const char tag[2], int64_t val) { //add or update int tag
    	return bam_aux_update_int(b, tag, val);
    }
 //--query methods:
 uint32_t flags() { return b->core.flag; } //return SAM flags
 bool isUnmapped() { return ((b->core.flag & BAM_FUNMAP) != 0); }
 bool isMapped() { return ((b->core.flag & BAM_FUNMAP) == 0); }
 bool isPaired() { return ((b->core.flag & BAM_FPAIRED) != 0); }
 const char* name() { return bam_get_qname(b); }
 int pairOrder() {
    //which read in the pair: 0 = unpaired, 1=first read, 2=second read
    int r=0;
    if ((b->core.flag & BAM_FREAD1) != 0) r=1;
    else if ((b->core.flag & BAM_FREAD2) != 0) r=2;
    return r;
    }
 bool revStrand() {
   //this is the raw alignment strand, NOT the transcription/splice strand
   return ((b->core.flag & BAM_FREVERSE) != 0);
 }
 const char* refName() {
   return (b_hdr!=NULL) ?
         ((b->core.tid<0) ? "*" : b_hdr->target_name[b->core.tid]) : NULL;
   }
 int32_t refId() { return b->core.tid; }
 int32_t mate_refId() { return b->core.mtid; }
 const char* mate_refName() {
    return (b_hdr!=NULL) ?
       ((b->core.mtid<0) ? "*" : b_hdr->target_name[b->core.mtid]) : NULL;
    }
 int32_t insertSize() { return b->core.isize; }
 int32_t mate_start() { return b->core.mpos<0? 0 : b->core.mpos+1; }

 //int find_tag(const char tag[2], uint8_t* & s, char& tag_type);
 uint8_t* find_tag(const char tag[2]);
 //position s at the beginning of tag data, tag_type is set to the found tag type
 //returns length of tag data, or 0 if tag not found

 char* tag_str(const char tag[2]); //return tag value for tag type 'Z'
 int64_t tag_int(const char tag[2]); //return numeric value of tag (for numeric types)
 double tag_float(const char tag[2]); //return float value of tag (for float types)
 char tag_char(const char tag[2]); //return char value of tag (for type 'A')
 char tag_char1(const char tag[2]);
 char spliceStrand(); // '+', '-' from the XS tag, or '.' if no XS tag
 char* sequence(); //user should free after use
 char* qualities();//user should free after use
 char* cigar(); //returns text version of the CIGAR string; user must free
};

// from sam.c:
#define FTYPE_BAM  1
#define FTYPE_READ 2

class GSamReader {
   htsFile* hts_file;
   char* fname;
   sam_hdr_t* hdr;
 public:
   void bopen(const char* filename, int32_t required_fields,
		   const char* cram_refseq=NULL) {
	      hts_file=hts_open(filename, "r");
	      if (hts_file==NULL)
	         GError("Error: could not open alignment file %s \n",filename);
	      if (hts_file->is_cram && cram_refseq!=NULL) {
	              hts_set_opt(hts_file, CRAM_OPT_REFERENCE, cram_refseq);
    	  }

          hts_set_opt(hts_file, CRAM_OPT_REQUIRED_FIELDS,
	    			  required_fields);
	      fname=Gstrdup(filename);
	      hdr=sam_hdr_read(hts_file);
   }

   void bopen(const char* filename, const char* cram_refseq=NULL) {
      hts_file=hts_open(filename, "r");
      if (hts_file==NULL)
         GError("Error: could not open alignment file %s \n",filename);
      if (hts_file->is_cram) {
    	  if (cram_refseq!=NULL) {
              hts_set_opt(hts_file, CRAM_OPT_REFERENCE, cram_refseq);
    	  }
    	  else hts_set_opt(hts_file, CRAM_OPT_REQUIRED_FIELDS,
    		SAM_QNAME|SAM_FLAG|SAM_RNAME|SAM_POS|SAM_MAPQ|SAM_CIGAR|
			SAM_RNEXT|SAM_PNEXT|SAM_TLEN|SAM_AUX );

      }
      fname=Gstrdup(filename);
      hdr=sam_hdr_read(hts_file);
   }

   GSamReader(const char* fn, int32_t required_fields,
		   const char* cram_ref=NULL):hts_file(NULL),fname(NULL), hdr(NULL) {
      bopen(fn, required_fields, cram_ref);
   }

   GSamReader(const char* fn, const char* cram_ref=NULL):hts_file(NULL),fname(NULL), hdr(NULL) {
      bopen(fn, cram_ref);
   }

   sam_hdr_t* header() {
      return hts_file ? hdr : NULL;
   }
   const char* fileName() {
      return fname;
   }

   void bclose() {
      if (hts_file) {
   	    if (hdr!=NULL) sam_hdr_destroy(hdr);
   	    hdr=NULL;
        hts_close(hts_file);
        hts_file=NULL;
        }
    }

   ~GSamReader() {
      bclose();
      GFREE(fname);
    }
   /*
   int64_t fpos() { //ftell
     if (hts_file->is_bgzf) { // bam_ptell() from sam.c
    	    if (hts_file->fp.bgzf==NULL)
    	        return -1L;
    	    return bgzf_tell(hts_file->fp.bgzf);
     }
     else if (hts_file->is_cram) { // cram_ptell() from sam.c
    	    cram_container *c;
    	    cram_slice *s;
    	    int64_t ret = -1L;
    	    if (hts_file->fp.cram) {
    	        if ((c = hts_file->fp.cram->ctr) != NULL) {
    	            if ((s = c->slice) != NULL && s->max_rec) {
    	                if ((c->curr_slice + s->curr_rec/s->max_rec) >= (c->max_slice + 1))
    	                	hts_file->fp.cram->curr_position += c->offset + c->length;
    	            }
    	        }
    	        ret = hts_file->fp.cram->curr_position;
    	    }

    	    return ret;
     }
     else {
    	 return htell(hts_file->fp.hfile);
     }

   }

   int64_t fseek(int64_t offs) {
    if (hts_file->is_bgzf) { //bam_pseek() from sam.c
    	 return bgzf_seek(hts_file->fp.bgzf, offs, SEEK_SET);
       }
    else if (hts_file->is_cram) { //cram_pseek() from sam.c
        if ((0 != cram_seek(hts_file->fp.cram, offs, SEEK_SET))
         && (0 != cram_seek(hts_file->fp.cram, offs - hts_file->fp.cram->first_container, SEEK_CUR)))
            return -1;

        hts_file->fp.cram->curr_position = offs;

        if (hts_file->fp.cram->ctr) {
            cram_free_container(hts_file->fp.cram->ctr);
            if (hts_file->fp.cram->ctr_mt && hts_file->fp.cram->ctr_mt != hts_file->fp.cram->ctr)
                cram_free_container(hts_file->fp.cram->ctr_mt);

            hts_file->fp.cram->ctr = NULL;
            hts_file->fp.cram->ctr_mt = NULL;
            hts_file->fp.cram->ooc = 0;
        }
        return offs;
      }
    else
        return hseek(hts_file->fp.hfile, offs, SEEK_SET);
   }
   */
   void rewind() {
     if (fname==NULL) {
       GMessage("Warning: GSamReader::rewind() called without a file name.\n");
       return;
     }
     bclose();
     char* ifname=fname;
     bopen(ifname);
     GFREE(ifname);
  }

  GSamRecord* next() {
      if (hts_file==NULL)
        GError("Warning: GSamReader::next() called with no open file.\n");
      bam1_t* b = bam_init1();
      if (sam_read1(hts_file, hdr, b) >= 0) {
        GSamRecord* bamrec=new GSamRecord(b, hdr, true);
        return bamrec;
        }
      else {
        bam_destroy1(b);
        return NULL;
        }
      }
};

//basic BAM/SAM/CRAM writer class
// limitations: cannot add new reference sequences, just new alignments to
//  existing reference sequences;
class GSamWriter {
   htsFile* bam_file;
   sam_hdr_t* hdr;
 public:
   void create(const char* fname, sam_hdr_t* bh, GSamFileType ftype=GSamFile_BAM) {
     hdr=sam_hdr_dup(bh);
     create(fname, ftype);
   }

   GSamWriter(const char* fname, sam_hdr_t* bh, GSamFileType ftype=GSamFile_BAM):
	                                    bam_file(NULL),hdr(NULL) {
      create(fname, bh, ftype);
   }

   void create(const char* fname, GSamFileType ftype=GSamFile_BAM) {
      if (hdr==NULL)
         GError("Error: no header data provided for GSamWriter::create()!\n");
	  kstring_t mode=KS_INITIALIZE;
      kputc('w', &mode);
      switch (ftype) {
         case GSamFile_BAM:
        	kputc('b', &mode);
        	break;
         case GSamFile_UBAM:
        	kputs("bu", &mode);
        	break;
         case GSamFile_CRAM:
        	kputc('c', &mode);
        	break;
         case GSamFile_SAM:
         	break;
         default:
      	   GError("Error: unrecognized output file type!\n");
      }
      bam_file = hts_open(fname, mode.s);
      if (bam_file==NULL)
         GError("Error: could not create output file %s\n", fname);
      if (sam_hdr_write(bam_file, hdr)<0)
    	  GError("Error writing header data to file %s\n", fname);
   }
   sam_hdr_t* header() { return hdr; }
   GSamWriter(const char* fname, const char* hdr_file, GSamFileType ftype=GSamFile_BAM):
	                                             bam_file(NULL),hdr(NULL) {
	  //create an output file fname with the SAM header copied from hdr_file
      htsFile* samf=hts_open(hdr_file, "r");
      if (samf==NULL)
    	  GError("Error: could not open SAM file %s\n", hdr_file);
      hdr=sam_hdr_read(samf);
      if (hdr==NULL)
    	  GError("Error: could not read header data from %s\n", hdr_file);
      hts_close(samf);
      create(fname, ftype);
   }

   ~GSamWriter() {
      hts_close(bam_file);
      sam_hdr_destroy(hdr);
   }

   sam_hdr_t* get_header() { return hdr; }

   int32_t get_tid(const char *seq_name) {
      if (hdr==NULL)
         GError("Error: missing SAM header (get_tid())\n");
      return sam_hdr_name2tid(hdr, seq_name);
      }

   //just a convenience function for creating a new record, but it's NOT written
   //given pos must be 1-based (so it'll be stored as pos-1 because BAM is 0-based)
   GSamRecord* new_record(const char* qname, const char* gseqname,
            int pos, bool reverse, const char* qseq, const char* cigar=NULL, const char* qual=NULL) {
      int32_t gseq_tid=get_tid(gseqname);
      if (gseq_tid < 0 && strcmp(gseqname, "*")) {
            if (hdr->n_targets == 0) {
               GError("Error: missing/invalid SAM header\n");
               } else
                   GMessage("Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   gseqname);
            }

      return (new GSamRecord(qname, gseq_tid, pos, reverse, qseq, cigar, qual));
      }

   GSamRecord* new_record(const char* qname, int32_t samflags, const char* gseqname,
         int pos, int map_qual, const char* cigar, const char* mgseqname, int mate_pos,
         int insert_size, const char* qseq, const char* quals=NULL,
                          GVec<char*>* aux_strings=NULL) {
      int32_t gseq_tid=get_tid(gseqname);
      if (gseq_tid < 0 && strcmp(gseqname, "*")) {
            if (hdr->n_targets == 0) {
               GError("Error: missing/invalid SAM header\n");
               } else
                   GMessage("Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   gseqname);
            }
      int32_t mgseq_tid=-1;
      if (mgseqname!=NULL) {
         if (strcmp(mgseqname, "=")==0) {
            mgseq_tid=gseq_tid;
            }
          else {
            mgseq_tid=get_tid(mgseqname);
            if (mgseq_tid < 0 && strcmp(mgseqname, "*")) {
                GMessage("Warning: reference '%s' not found in header, will consider it '*'.\n",
                                   mgseqname);
                }
            }
          }
      return (new GSamRecord(qname, samflags, gseq_tid, pos, map_qual, cigar,
              mgseq_tid, mate_pos, insert_size, qseq, quals, aux_strings));
      }

   void write(GSamRecord* brec) {
      if (brec!=NULL) {
          if (sam_write1(this->bam_file,this->hdr, brec->b)<0)
        	  GError("Error writing SAM record!\n");
      }
   }

   void write(bam1_t* xb) {
     if (sam_write1(this->bam_file, this->hdr, xb)<0)
    	 GError("Error writing SAM record!\n");
   }
};

#endif
