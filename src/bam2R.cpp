/**********************************************************************
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015-2018 drjsanger@github
 ***********************************************************************/

#include <stdio.h>
#include <string.h>
#include "htslib/sam.h"
#include "htslib/khash.h"
#include <map>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

using namespace std;
KHASH_MAP_INIT_STR(strh,uint8_t) //readname -> readbase map used to prevent overlapping reads calls happening twice

//const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

typedef struct {
	int beg, end, q, s, head_clip;
	int i;
	int* counts;
	map<char, int> nt_idx;
	htsFile *in;
} nttable_t;

char NUCLEOTIDES[] = {'A','T','C','G','*','N','+','-','^','$','Q'};
int N = 11;

extern "C" {

void bam2R_pileup_function(const bam_pileup1_t *pl, int pos, int n_plp, nttable_t& nttable)
{
  int i, s;
  int len = nttable.end - nttable.beg;
  map<char, int> nt_freq;
	khash_t(strh) *h;
	khiter_t k;
	h = kh_init(strh);
  if ((int)pos >= nttable.beg && (int)pos < nttable.end)
  {

    int* counts = nttable.counts + (int)pos - nttable.beg ;
    for (i=0; i<n_plp; i++)
    {
      const bam_pileup1_t *p = pl + i;
      s = bam_is_rev(p->b) * len * N;
			int absent;
	    k = kh_put(strh, h, bam_get_qname(p->b), &absent);
			uint8_t cbase = bam_seqi(bam_get_seq(p->b),p->qpos);
			uint8_t pre_b;
			if(!absent){ //Read already processed to get base processed (we only increment if base is different between overlapping read pairs)
				k = kh_get(strh, h, bam_get_qname(p->b));
				pre_b = kh_val(h,k);
			}else{
				//Add the value to the hash
				kh_value(h, k) = cbase;
			}

      {
				if(!absent && pre_b == cbase) continue;
        if (p->is_tail) counts[s + len * nttable.nt_idx['$']]++;
        else if (p->is_head) counts[s + len * nttable.nt_idx['^']]++;

        if(p->qpos < nttable.head_clip || (bam_is_rev(p->b) && ((bam1_core_t)p->b->core).l_qseq - p->qpos < nttable.head_clip)){
          counts[s + len * nttable.nt_idx['N']]++;
        }else{
          if (!p->is_del) {
            int  c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
            if( bam_get_qual(p->b)[p->qpos] > nttable.q)
              counts[s + len * nttable.nt_idx[char(c)]]++;
            else
              counts[s + len * nttable.nt_idx['N']]++;
            if (p->indel > 0)
              counts[s + len * nttable.nt_idx['+']]++;
            else if (p->indel < 0)
              counts[s + len * nttable.nt_idx['-']]++;
          } else counts[s + len * nttable.nt_idx['*']]++;
          counts[s + len * nttable.nt_idx['Q']] += p->b->core.qual;
        }
      }
    }
    nttable.i++;
  }
	kh_destroy(strh, h);
}

static inline int64_t getNM(const bam1_t *b, unsigned long long& count)
{
	const uint8_t *nm = bam_aux_get(b, "NM");
	if (nm)
		return bam_aux2i(nm);
	else {
		count++;
		return 0;  // Dummy NM value that always passes the filter
	}
}

int bam2R(char** bamfile, char** ref, int* beg, int* end, int* counts, int* q, int* mq, int* s, 
          int* head_clip, int* maxdepth, int* verbose, int* mask, int *keepflag, int *maxmismatches )
{

	bam_plp_t buf = NULL;
	bam1_t *b = NULL;
	bam_hdr_t *head = NULL;

	nttable_t nttable;
	nttable.q = *q; //Base quality cutoff
	nttable.s = *s; //Strand (2=both)
	nttable.head_clip = *head_clip;
	nttable.i = 0;
	nttable.counts = counts;

	int64_t maxNM = (*maxmismatches != -1)? *maxmismatches : INT64_MAX;
	unsigned long long no_NM_count = 0;

	int i;
	for (i=0; i<N; i++)
		nttable.nt_idx[NUCLEOTIDES[i]] = i;

	nttable.beg = *beg -1;
	nttable.end = *end;
	nttable.in = hts_open(*bamfile, "r");
	if (nttable.in == 0) {
		Rf_error("Fail to open input BAM/CRAM file %s\n", *bamfile);
		return 1;
	}

	buf = bam_plp_init(0,(void *)&nttable); // initialize pileup
	bam_plp_set_maxcnt(buf,*maxdepth);
	b = bam_init1();
	//get header
	head = sam_hdr_read(nttable.in);
	//int mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY;
  int tid, pos, n_plp = -1;
	const bam_pileup1_t *pl;

	if (strcmp(*ref, "") == 0) { // if a region is not specified
		//Replicate sampileup functionality (uses above mask without supplementary)
		int ret;
		while((ret = sam_read1(nttable.in, head, b)) >= 0){
			if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && getNM(b, no_NM_count) <= maxNM) {
					bam_plp_push(buf, b);
            };
			while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
				bam2R_pileup_function(pl,pos,n_plp,nttable);
			}
		}
	}
	else {
		int tid;
		hts_idx_t *idx;
		idx = sam_index_load(nttable.in,*bamfile); // load BAM index
		if (idx == 0) {
			Rf_error("BAM/CRAM index file is not available.\n");
			return 1;
		}
		tid = bam_name2id(head, *ref);
		if (tid < 0) {
			Rf_error("Invalid sequence %s\n", *ref);
			return 1;
		}

		if(*verbose)
			Rprintf("Reading %s, %s:%d-%d\n", *bamfile, *ref, nttable.beg+1, nttable.end);

		//Implement a fetch style iterator
		hts_itr_t *iter = sam_itr_queryi(idx, tid, nttable.beg, nttable.end);
		int result;
		while ((result = sam_itr_next(nttable.in, iter, b)) >= 0) {
			if ((b->core.flag & *mask)==0 && b->core.qual >= *mq && (b->core.flag & *keepflag)==*keepflag && getNM(b, no_NM_count) <= maxNM) {
				bam_plp_push(buf, b);
			};
			while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
				bam2R_pileup_function(pl,pos,n_plp,nttable);
			}
		}
    if(result < -1){
      Rf_error("Error code (%d) encountered reading sam iterator.\n", result);
			return 1;
    }
		sam_itr_destroy(iter);
		hts_idx_destroy(idx);
	}

	bam_plp_push(buf,0); // finalize pileup

  while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) != 0) {
    bam2R_pileup_function(pl,pos,n_plp,nttable);
  }

	if (*maxmismatches != -1 && no_NM_count > 0) {
		Rf_warning("%llu reads did not have NM tags; max.mismatches filter was not applied to them.\n", no_NM_count);
	}

	bam_destroy1(b);
	bam_hdr_destroy(head);
	bam_plp_destroy(buf);
	hts_close(nttable.in);
	return 0;
}

R_CMethodDef cMethods[] = {
		{"bam2R", (DL_FUNC) &bam2R, 12}
};

void R_init_bam2R(DllInfo *info) {
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

} // extern "C"
