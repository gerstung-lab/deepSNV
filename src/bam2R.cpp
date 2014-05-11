/**********************************************************************
* bam2R.cpp An interface for R to count nucleotides in a .bam alignment
* Copyright (C) 2011 Moritz Gerstung, ETH Zurich, gemoritz@ethz.ch
***********************************************************************/

#include <stdio.h>
#include <string.h>
#include "samtools/sam.h"
#include <map>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

using namespace std;


typedef struct {
  int beg, end, q, s, head_clip;
  int i;
  int* counts;
  map<char, int> nt_idx;
  samfile_t *in;
} nttable_t;

char NUCLEOTIDES[] = {'A','T','C','G','*','N','+','-','^','$','Q'};
int N = 11;


// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
  bam_plbuf_t *buf = (bam_plbuf_t*)data;
  bam_plbuf_push(b, buf);
  return 0;
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	nttable_t *nttable = (nttable_t*)data;
	int i, s;
	int len = nttable->end - nttable->beg;
	map<char, int> nt_freq;
	if ((int)pos >= nttable->beg && (int)pos < nttable->end)
	{
		int* counts = nttable->counts + (int)pos - nttable->beg ;
		for (i=0; i<n; i++)
		{
			const bam_pileup1_t *p = pl + i;
			//if((bam1_strand(p->b)==nttable->s || nttable->s==2)){
			s = bam1_strand(p->b) * len * N;
			{
				if (p->is_tail) counts[s + len * nttable->nt_idx['$']]++;
				else if (p->is_head) counts[s + len * nttable->nt_idx['^']]++;

				if(p->qpos < nttable->head_clip || (bam1_strand(p->b) && ((bam1_core_t)p->b->core).l_qseq - p->qpos < nttable->head_clip)){
					counts[s + len * nttable->nt_idx['N']]++;
				}else{
					if (!p->is_del) {
						int  c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
						if( bam1_qual(p->b)[p->qpos] > nttable->q)
							counts[s + len * nttable->nt_idx[char(c)]]++;
						else
							counts[s + len * nttable->nt_idx['N']]++;
						if (p->indel > 0)
							counts[s + len * nttable->nt_idx['+']]++;
						else if (p->indel < 0)
							counts[s + len * nttable->nt_idx['-']]++;
					} else counts[s + len * nttable->nt_idx['*']]++;
					counts[s + len * nttable->nt_idx['Q']] += p->b->core.qual;
				}
			}
		}
		nttable->i++;
	}
	return 0;
}

extern "C" {
void bam_init_header_hash(bam_header_t *header);

int bam2R(char** bamfile, char** ref, int* beg, int* end, int* counts, int* q, int* s, int* head_clip, int* maxdepth, int* verbose)
{
	int c = 0;
	nttable_t nttable;
	nttable.q = *q; //Base quality cutoff
	nttable.s = *s; //Strand (2=both)
	nttable.head_clip = *head_clip;
	nttable.i = 0;
	nttable.counts = counts;

	int i;
	for (i=0; i<N; i++)
		nttable.nt_idx[NUCLEOTIDES[i]] = i;

	nttable.beg = *beg -1;
	nttable.end = *end;
	nttable.in = samopen(*bamfile, "rb", 0);
	if (nttable.in == 0) {
		Rf_error("Fail to open BAM file %s\n", *bamfile);
		return 1;
	}

	if (strcmp(*ref, "") == 0) { // if a region is not specified
		sampileup(nttable.in, -1, pileup_func, &nttable);
	}
	else {
		int tid;
		bam_index_t *idx;
		bam_plbuf_t *buf;
		idx = bam_index_load(*bamfile); // load BAM index
		if (idx == 0) {
			Rf_error("BAM indexing file is not available.\n");
			return 1;
		}

		bam_init_header_hash(nttable.in->header);
		tid = bam_get_tid(nttable.in->header, *ref);
		if (tid < 0) {
			Rf_error("Invalid sequence %s\n", *ref);
			return 1;
		}
		if(*verbose)
			Rprintf("Reading %s, %s:%i-%i\n", *bamfile, *ref, nttable.beg, nttable.end);

		buf = bam_plbuf_init(pileup_func, &nttable); // initialize pileup
		bam_plp_set_maxcnt(buf->iter, *maxdepth);
		bam_fetch(nttable.in->x.bam, idx, tid, nttable.beg, nttable.end, buf, fetch_func);
		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
	}
	samclose(nttable.in);
	return 0;
}

R_CMethodDef cMethods[] = {
   {"bam2R", (DL_FUNC) &bam2R, 10}
};

void R_init_bam2R(DllInfo *info) {
   R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

} // extern "C"
