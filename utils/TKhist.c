/*
 * TKhist.c
 *
 *  Created on: 30 Nov 2016
 *      Author: pippelmn
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <sys/param.h>

#include "lib/colors.h"
#include "lib/oflags.h"
#include "lib/utils.h"
#include "lib/pass.h"
#include "lib/tracks.h"
#include "lib/stats.h"

#include "lib.ext/types.h"
#include "lib.ext/bitarr.h"

#include "db/DB.h"
#include "dalign/align.h"

#undef DEBUG

#define DEF_ARG_Q TRACK_Q
#define DEF_ARG_T TRACK_TRIM
#define DEF_ARG_R TRACK_REPEATS

typedef struct
{
	HITS_DB* db;
	HITS_TRACK *qTrack;
	HITS_TRACK *rTrack;
	HITS_TRACK *tTrack;

	// file names
	char* pcPathReadsIn;
	char* pcPathOverlaps;

	int verbose;
	int nbin;
	int64 *hist;

	bit* vReads;

	ovl_header_twidth twidth;

} TKhistContext;

static void usage()
{
	printf("[-v] [-bs <int>] [-qrt track] <db> <ovl>\n");

	printf("options: -v ... verbose\n");
	printf("         -q ... q track (%s)\n", DEF_ARG_Q);
	printf("         -r ... repeats track (%s)\n", DEF_ARG_R);
	printf("         -t ... trim track (%s)\n", DEF_ARG_T);
}

static void tkhist_pre(PassContext *pctx, TKhistContext* tctx)
{
	tctx->nbin = 100;
	tctx->hist = (int64 *) malloc(sizeof(int64) * tctx->nbin);

	bzero(tctx->hist, sizeof(int64) * tctx->nbin);

	tctx->vReads = ba_new(DB_NREADS(tctx->db));
	assert(tctx->vReads != NULL);

	tctx->twidth = pctx->twidth;
}

static void tkhist_post(TKhistContext* tctx)
{
	int i;
	printf("qValue num\n");
	for (i = tctx->nbin - 1; i >= 0; i--)
	  {
		//if (tctx->hist[i] > 0)
		  {
			printf("%d %lld\n", i, tctx->hist[i]);
		  }
	  }
}

static void addQvalues(TKhistContext *tctx, int rid)
{
	int trim_b, trim_e;

	track_anno* ranno = (track_anno*) (tctx->rTrack->anno);
	track_data* rdata = (track_data*) (tctx->rTrack->data);

	track_anno* qanno = (track_anno*) (tctx->qTrack->anno);
	track_data* qdata = (track_data*) (tctx->qTrack->data);

	// get trim track
	get_trim(tctx->db, tctx->tTrack, rid, &trim_b, &trim_e);
	//printf("%7d | trim  %5d..%5d\n", rid, trim_b, trim_e);

	if (trim_e - trim_b < 1000)
	{
		printf("trim interval too short!\n");
		return;
	}

	track_anno rob = ranno[rid] / sizeof(track_data);
	track_anno roe = ranno[rid + 1] / sizeof(track_data);

	assert(rob <= roe);

	track_anno qob = qanno[rid] / sizeof(track_data);
	track_anno qoe = qanno[rid + 1] / sizeof(track_data);

	assert(qob <= qoe);

	// no repeats, take all complete segments within trim interval

	if (rob >= roe)
	{
		int c = 0;
		while (qob < qoe)
		{
			if(c >= trim_b + tctx->twidth && c <= trim_e - tctx->twidth)
			{
				int bin = qdata[qob];
				if(bin > 100)
					bin = 100;
				if(tctx->verbose)
				if(bin > 60 || bin == 0)
					printf("Allert1! q%d in read %d in [%d, %d] trim [%d, %d]\n",bin, rid, c, c+tctx->twidth-1, trim_b, trim_e);
				tctx->hist[bin]++;
#ifdef DEBUG
				printf("%d - %d: q%d\n", rid, c, qdata[qob]);
#endif
			}
			qob++;
			c += tctx->twidth;
		}
	}
	else
	{
		int uniq_b, uniq_e;
		uniq_b=trim_b;
		int c = 0;
		while (rob < roe)
		{
#ifdef DEBUG
			printf("%d trim [%d, %d] repeat [%d, %d]\n", rid, trim_b, trim_e, rdata[rob], rdata[rob+1]);
#endif
			uniq_e = MIN(rdata[rob], trim_e);
			if(uniq_e - uniq_b > tctx->twidth*3)
			{
#ifdef DEBUG
				printf("valid range: %d - %d\n", uniq_b, uniq_e);
#endif
				//c=0;
				while (qob < qoe && c < uniq_e)
				{
					if(c >= uniq_b + tctx->twidth && c <= uniq_e - tctx->twidth)
					{
						int bin = qdata[qob];
						if(bin > 100)
							bin = 100;

						if(tctx->verbose)
						if(bin > 60 || bin == 0)
							printf("Allert2! q%d in read %d in [%d, %d] trim [%d, %d]\n",bin, rid, c, c+tctx->twidth-1, trim_b, trim_e);

						tctx->hist[bin]++;
#ifdef DEBUG
						printf("%d - %d: q%d\n", rid, c, qdata[qob]);
#endif
					}
					qob++;
					c += tctx->twidth;
				}

			}
			uniq_b=MAX(rdata[rob+1], trim_b);
			rob+=2;
		}
		//c=0;
		if(trim_e - uniq_b > tctx->twidth*3)
		{
			uniq_e = trim_e;
			while (qob < qoe && c < uniq_e)
			{
				if(c >= uniq_b + tctx->twidth && c <= uniq_e - tctx->twidth)
				{
					int bin = qdata[qob];
					if(bin > 100)
						bin = 100;

					if(tctx->verbose)
                                                if(bin > 60 || bin == 0)
                                                        printf("Allert2! q%d in read %d in [%d, %d] trim [%d, %d]\n",bin, rid, c, c+tctx->twidth-1, trim_b, trim_e);
					
					tctx->hist[bin]++;
#ifdef DEBUG
					printf("%d - %d: q%d\n", rid, c, qdata[qob]);
#endif
				}
				qob++;
				c += tctx->twidth;
			}

		}

	}
}

static int tkhist_handler(void* _ctx, Overlap *ovls, int novl)
{
	TKhistContext *tctx = (TKhistContext*) _ctx;

	int i;
	int aread, bread;

	aread = ovls->aread;

	if (ba_value(tctx->vReads, aread) == FALSE)
	{
		addQvalues(tctx, aread);
		ba_assign(tctx->vReads, aread, TRUE);
	}

	for (i = 0; i < novl; i++)
	{
		Overlap *ovl = ovls + i;

		bread = ovl->bread;

		if (ba_value(tctx->vReads, bread) == TRUE)
		{
			continue;
		}
		addQvalues(tctx, bread);

		ba_assign(tctx->vReads, bread, TRUE);
	}

	return 1;
}

int main(int argc, char* argv[])
{
	HITS_DB db;
	FILE* fileOvlIn;

	PassContext *pctx;
	TKhistContext tctx;
	bzero(&tctx, sizeof(TKhistContext));

	// set default values
	tctx.verbose = 0;
	char *qName = DEF_ARG_Q;
	char *rName = DEF_ARG_R;
	char *tName = DEF_ARG_T;

// process arguments

	int c;

	opterr = 0;

	while ((c = getopt(argc, argv, "vq:r:t:")) != -1)
	{
		switch (c)
		{
		case 'v':
			tctx.verbose++;
			break;
		case 'q':
			qName = optarg;
			break;
		case 'r':
			rName = optarg;
			break;
		case 't':
			tName = optarg;
			break;
		default:
			usage();
			exit(1);
		}
	}

	if (argc - optind == 2)
	{
		tctx.pcPathReadsIn = argv[optind++];
		tctx.pcPathOverlaps = argv[optind++];
	}
	else
	{
		usage();
		exit(1);
	}

	if ((fileOvlIn = fopen(tctx.pcPathOverlaps, "r")) == NULL)
	{
		fprintf(stderr, "could not open '%s'\n", tctx.pcPathOverlaps);
		exit(1);
	}

	if (Open_DB(tctx.pcPathReadsIn, &db))
	{
		printf("could not open '%s'\n", tctx.pcPathReadsIn);
	}

	// load tracks
	if ((tctx.qTrack = track_load(&db, qName)) == NULL)
	{
		printf("Could not load track %s\n", qName);
		exit(1);
	}
	if ((tctx.rTrack = track_load(&db, rName)) == NULL)
	{
		printf("Could not load track %s\n", rName);
		exit(1);
	}
	if ((tctx.tTrack = track_load(&db, tName)) == NULL)
	{
		printf("Could not load track %s\n", tName);
		exit(1);
	}

	tctx.db = &db;

	pctx = pass_init(fileOvlIn, NULL);

	pctx->split_b = 0;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &tctx;

	tkhist_pre(pctx, &tctx);

	pass(pctx, tkhist_handler);

	tkhist_post(&tctx);

	Close_DB(&db);

	return 0;
}

