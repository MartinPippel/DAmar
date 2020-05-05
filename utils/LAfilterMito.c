/*******************************************************************************************
 *
 *  filters overlaps by various criteria
 *
 *  Author :  Martin
 *
 *  Date   :  December 2018
 *
 *******************************************************************************************/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>

#include "lib/colors.h"
#include "lib/oflags.h"
#include "lib/pass.h"
#include "lib/read_loader.h"
#include "lib/tracks.h"
#include "lib/trim.h"
#include "lib/utils.h"

#include "dalign/align.h"
#include "db/DB.h"

#define DEF_ARG_U 1000
#define DEF_ARG_G 1000
#define DEF_ARG_O 25
#define DEF_ARG_C 75

#define VERBOSE

typedef struct
{
	// stats counters
	int nFiltReads;
	int nFiltAlns;
	int nFilteredMinLenReads;
	int nFilteredMaxLenReads;
	int nFilteredMinLenAln;

	int nMinReadLength;
	int nMaxReadLength;
	int nMaxUnalignedTips;
	int nMaxOverhangLength;
	int nMaxGapLength;
	int nMinCoveredReadLength;

	// settings
	int nVerbose;

	HITS_DB* db;

	ovl_header_twidth twidth;


} FilterContext;

extern char* optarg;
extern int optind, opterr, optopt;

static int cyclicChain(FilterContext *ctx, Overlap *ovls, int n)
{
	int aread, bread;
	int alen, blen;
	int i;

	// aread should be the mito reference
	aread = ovls->aread;
	// bread raw reads, that we want to fish out
	bread = ovls->bread;

	alen = DB_READ_LEN(ctx->db, aread);
	blen = DB_READ_LEN(ctx->db, bread);

	// only pull out reads that overlap front to end with following criteria

	// 1. all alignments map in one orientation only
	// 2. are within given length thresholds [Min, Max] --> exclude missed adapters, nucleic DNA
	// 3. fulfill maximum unaligned sequence

	Overlap *ovl1 = ovls;

	if (ovl1->path.abpos > ctx->nMaxUnalignedTips && ovl1->path.bbpos > ctx->nMaxUnalignedTips)
	{
	  printf("ctx->nMaxUnalignedTip! : %d > %d && %d > %d\n", ovl1->path.abpos, ctx->nMaxUnalignedTips, ovl1->path.bbpos, ctx->nMaxUnalignedTips);
	  return 0;
	}

	if (n == 1)
	{
		if (alen - ovl1->path.aepos > ctx->nMaxUnalignedTips && blen - ovl1->path.bepos > ctx->nMaxUnalignedTips)
			return 0;

		if ((ovl1->path.bepos - ovl1->path.bbpos)*100.0/blen < ctx->nMinCoveredReadLength)
			return 0;

		return 1;
	}

	int numCyclicBreaks = 0;
	int bCumAlnLen = ovl1->path.bepos - ovl1->path.bbpos;

	for (i = 1; i < n; i++)
	{
		Overlap *ovl2 = ovls + i;

		bCumAlnLen += ovl2->path.bepos - ovl2->path.bbpos;

		if ((ovl1->flags & OVL_COMP) != (ovl2->flags & OVL_COMP))
		{
			printf("wrong orientation\n");
			return 0;
		}
		// TODO: urgent check gaps at the very end, otherwise we loos to many reads
		if(ovl2->path.abpos - ovl1->path.aepos > ctx->nMaxGapLength)
		{
			printf("gap length to large: %d > %d\n",ovl2->path.abpos - ovl1->path.aepos, ctx->nMaxGapLength );
			return 0;
		}

		if (ovl1->path.aepos - ctx->nMaxOverhangLength < ovl2->path.abpos)
		{
			// check if we have a valid extension into proper direction
			if (ovl1->path.bepos - ctx->nMaxOverhangLength < ovl2->path.bbpos)
			{
				if(ovl2->path.bbpos - ovl1->path.bepos > ctx->nMaxGapLength)
					return 0;

				ovl1 = ovl2;
				continue;
			}
			// check if we have a cyclic alignment
			if (blen - ovl1->path.bepos < ctx->nMaxUnalignedTips && ovl2->path.bepos - ctx->nMaxUnalignedTips < ovls->path.bbpos)
			{
				ovl1 = ovl2;
				numCyclicBreaks++;
				continue;
			}
		}
		return 0;
	}

	if (numCyclicBreaks > 1)
        {
		printf("numCyclicBreaks(%d) > 1\n",numCyclicBreaks);
		return 0;
	}
	// check covered b-read percentage
	if(bCumAlnLen*100.0/blen < ctx->nMinCoveredReadLength)
	{
		printf("bCumAlnLen*100.0/blen < ctx->nMinCoveredReadLength: %.3f < %d\n", bCumAlnLen*100.0/blen, ctx->nMinCoveredReadLength);
		return 0;
	}


	return 1;
}

static void filter_pre(PassContext* pctx, FilterContext* fctx)
{
#ifdef VERBOSE
	printf( ANSI_COLOR_GREEN "PASS filtering\n" ANSI_COLOR_RESET);

	printf( ANSI_COLOR_RED "   OPTIONS\n");
	printf( "   verbose : %d\n", fctx->nVerbose);
	printf( "   nMinReadLength : %d\n", fctx->nMinReadLength);
	printf( "   nMinCoveredReadLength : %d\n", fctx->nMinCoveredReadLength);
	printf( "   nMaxUnalignedTips : %d\n", fctx->nMaxUnalignedTips);
	printf( "   nMaxReadLength : %d\n", fctx->nMaxReadLength);
	printf( "   nMaxOverhangLength : %d\n", fctx->nMaxOverhangLength);
	printf( "   nMaxGapLength : %d\n", fctx->nMaxGapLength);
	printf( "   nFilteredMinLenReads : %d\n", fctx->nFilteredMinLenReads);
	printf( "   nFilteredMinLenAln : %d\n", fctx->nFilteredMinLenAln);
	printf( "   nFilteredMaxLenReads : %d\n", fctx->nFilteredMaxLenReads);
	printf( "   DB name : %s\n", fctx->db->path);
	printf("\n" ANSI_COLOR_RESET);
#endif

	fctx->twidth = pctx->twidth;

}

static void filter_post(FilterContext* ctx)
{
	if (ctx->nVerbose)
	{
		if (ctx->nFiltAlns > 0)
		{
			printf("total discarded alignments %10d\n", ctx->nFiltAlns);
		}
	}

}

static int filter_handler(void* _ctx, Overlap* ovl, int novl)
{
	FilterContext* ctx = (FilterContext*) _ctx;
	int i, j, k;

	j = k = 0;

	while (j < novl)
	{
		int bReadLen = DB_READ_LEN(ctx->db, ovl[j].bread);
		int skipAlns = 0;

		if (ctx->nMinReadLength != -1 && (bReadLen < ctx->nMinReadLength || bReadLen > ctx->nMaxReadLength))
		{
			skipAlns = 1;
			ovl[j].flags |= OVL_DISCARD | OVL_RLEN;
		}

		while (k < novl - 1 && ovl[j].bread == ovl[k + 1].bread)
		{
			if (skipAlns)
				ovl[k + 1].flags |= OVL_DISCARD | OVL_RLEN;
			k++;
		}

		if (!skipAlns)
		{
			// discard all alignments no valid chain could be found
			if (cyclicChain(ctx, ovl + j, k - j + 1) == 0)
			{
				for (i = j; i<=k; i++)
				{
					printf("discard OVL: %d vs  %d [%d, %d] [%d, %d] comp? %c, l [%d, %d]\n", ovl[i].aread, ovl[i].bread, ovl[i].path.abpos, ovl[i].path.aepos, ovl[i].path.bbpos, ovl[i].path.bepos, (ovl[i].flags & OVL_COMP) ? 'C' : 'N', DB_READ_LEN(ctx->db, ovl[j].aread), DB_READ_LEN(ctx->db, ovl[j].bread));
					ovl[i].flags |= OVL_DISCARD;
				}
			}
		}

		j = k + 1;
	}

	return 1;
}

static void usage()
{
	fprintf(stderr, "[-vp] [-lLugoc <int>] <db> <overlaps_in> <overlaps_out>\n");

	fprintf(stderr, "options: -v ... verbose\n");
	fprintf(stderr, "         -l ... min read length\n");
	fprintf(stderr, "         -L ... max read length\n");
	fprintf(stderr, "         -p ... purge discarded overlaps\n");
	fprintf(stderr, "         -u ... maximum unaligned tips (default: %d)\n", DEF_ARG_U);
	fprintf(stderr, "         -g ... maximum gap length in overlap chain (default: %d)\n", DEF_ARG_G);
	fprintf(stderr, "         -o ... maximum overlap overhang of chained alignments (default: %d)\n", DEF_ARG_O);
	fprintf(stderr, "         -c ... minimum covered read length in percent (default: %d, must be in [50,100])\n", DEF_ARG_C);
}

int main(int argc, char* argv[])
{
	HITS_DB db;
	FilterContext fctx;
	PassContext* pctx;
	FILE* fileOvlIn;
	FILE* fileOvlOut;

	bzero(&fctx, sizeof(FilterContext));

	fctx.db = &db;

// args

	int arg_purge = 0;

	fctx.nMinReadLength = -1;
	fctx.nMaxReadLength = -1;
	fctx.nMaxUnalignedTips = DEF_ARG_U;
	fctx.nMaxOverhangLength = DEF_ARG_O;
	fctx.nMaxGapLength = DEF_ARG_G;
	fctx.nMinCoveredReadLength = DEF_ARG_C;
	fctx.nVerbose = 0;

	int c;

	opterr = 0;
	while ((c = getopt(argc, argv, "vpl:L:u:g:o:c:")) != -1)
	{
		switch (c)
		{
		case 'v':
			fctx.nVerbose += 1;
			break;

		case 'p':
			arg_purge = 1;
			break;

		case 'c':
			fctx.nMinCoveredReadLength = atoi(optarg);
			break;

		case 'u':
			fctx.nMaxUnalignedTips = atoi(optarg);
			break;

		case 'g':
			fctx.nMaxGapLength = atoi(optarg);
			break;

		case 'o':
			fctx.nMaxOverhangLength = atoi(optarg);
			break;

		case 'l':
			fctx.nMinReadLength = atoi(optarg);
			break;

		case 'L':
			fctx.nMaxReadLength = atoi(optarg);
			break;

		default:
			fprintf(stderr, "unknown option %c\n", c);
			usage();
			exit(1);
		}
	}

	if (argc - optind != 3)
	{
		usage();
		exit(1);
	}

	if(fctx.nMinCoveredReadLength < 50 || fctx.nMinCoveredReadLength > 100)
	{
		fprintf(stderr, "Minimum covered read length in percent: %d not supported! Valid range [50, 100]\n", fctx.nMinCoveredReadLength);
		exit(1);
	}

	char* pcPathReadsIn = argv[optind++];
	char* pcPathOverlapsIn = argv[optind++];
	char* pcPathOverlapsOut = argv[optind++];

	if ((fileOvlIn = fopen(pcPathOverlapsIn, "r")) == NULL)
	{
		fprintf(stderr, "could not open %s\n", pcPathOverlapsIn);
		exit(1);
	}

	if ((fileOvlOut = fopen(pcPathOverlapsOut, "w")) == NULL)
	{
		fprintf(stderr, "could not open %s\n", pcPathOverlapsOut);
		exit(1);
	}

	if (Open_DB(pcPathReadsIn, &db))
	{
		fprintf(stderr, "could not open %s\n", pcPathReadsIn);
		exit(1);
	}

// passes

	pctx = pass_init(fileOvlIn, fileOvlOut);

	pctx->split_b = 0;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &fctx;
	pctx->write_overlaps = 1;
	pctx->purge_discarded = arg_purge;

	filter_pre(pctx, &fctx);

	pass(pctx, filter_handler);

	filter_post(&fctx);

	pass_free(pctx);

// cleanup

	Close_DB(&db);

	fclose(fileOvlOut);
	fclose(fileOvlIn);

	return 0;
}
