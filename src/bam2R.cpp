/**********************************************************************
 * bamcram2R.cpp An interface for R to count nucleotides in a .bam
 * or .cram alignment
 * Copyright (C) 2015 drjsanger@github
 ***********************************************************************/

#include <stdio.h>
#include <string.h>
#include "sam.h"
#include <map>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

using namespace std;

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

//
static int pileup_function(void *data, bam1_t *b){
	return 0;
}

extern "C" {

int bam2R(char** bamfile, char** ref, int* beg, int* end, int* counts, int* q, int* s, int* head_clip, int* maxdepth, int* verbose)
{

	bam_plp_t buf = NULL;
	bam1_t *b = NULL;
	hts_itr_t *iter = NULL;
	bam_hdr_t *head = NULL;

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
	nttable.in = hts_open(*bamfile, "r");
	if (nttable.in == 0) {
		Rf_error("Fail to open input BAM/CRAM file %s\n", *bamfile);
		return 1;
	}

	buf = bam_plp_init(pileup_function,(void *)&nttable); // initialize pileup
	bam_plp_set_maxcnt(buf,*maxdepth);
	b = bam_init1();
	//get header
	head = sam_hdr_read(nttable.in);
	int mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY;

	if (strcmp(*ref, "") == 0) { // if a region is not specified
		//Replicate sampileup functionality (uses above mask without supplementary)
		int ret;
		while((ret = sam_read1(nttable.in, head, b)) >= 0){
			if (b->core.flag & mask) b->core.flag |= BAM_FUNMAP;
			bam_plp_push(buf, b);
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

		char *region = NULL;
		region = (char*) malloc(sizeof(*ref)+sizeof(":")+sizeof("-")+(sizeof(char)*50));
		sprintf(region,"%s:%d-%d",*ref,nttable.beg,nttable.end);
		if(*verbose)
			Rprintf("Reading %s, %s\n", *bamfile, region);

		//Implement a fetch style iterator
		hts_itr_t *iter = sam_itr_querys(idx, head, region);
		int result;
		while ((result = sam_itr_next(nttable.in, iter, b)) >= 0) {

			bam_plp_push(buf,b); //NB in your fetch you aren't doing the funky masking done in sampileup;
		}
		free(region);
		sam_itr_destroy(iter);
		hts_idx_destroy(idx);
	}

	bam_plp_push(buf,0); // finalize pileup
	int tid, pos, n_plp = -1;
	const bam_pileup1_t *pl;

	while ( (pl=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
		int i, s;
		int len = nttable.end - nttable.beg;
		map<char, int> nt_freq;
		if ((int)pos >= nttable.beg && (int)pos < nttable.end)
		{
			int* counts = nttable.counts + (int)pos - nttable.beg ;
			for (i=0; i<n_plp; i++)
			{
				const bam_pileup1_t *p = pl + i;
				s = bam_is_rev(p->b) * len * N;
				{
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
	}
	bam_destroy1(b);
	bam_hdr_destroy(head);
	bam_plp_destroy(buf);
	hts_close(nttable.in);
	return 0;
}

R_CMethodDef cMethods[] = {
		{"bam2R", (DL_FUNC) &bam2R, 10}
};

void R_init_bam2R(DllInfo *info) {
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

} // extern "C"

