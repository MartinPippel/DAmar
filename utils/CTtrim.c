/*******************************************************************************************
 *
 *  Trims contigs based on valid contig chain overlaps 
 *
 *  Author :  DAmar Team
 *
 *  Date   :  May 2020
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

#define MIN_BIONANO_GAP_SIZE 13
#define TRIM_OFFSET 100

typedef struct
{
	// stats counters
	int statsTrimmedContigs;
    int statsTrimmedBases;

    // db and I/O files
	HITS_DB* db;
    HITS_TRACK* trackDust;
	HITS_TRACK* trackTan;

    FILE* fileInBionanoGaps;
	FILE* fileOutTrimmedContigs;
    FILE* fileOutTrimmedStats;
    FILE* fileOutDicardedContigParts;

	ovl_header_twidth twidth;

    // matrix to store overlapping contigs 
    int ** LAStrimMatrix;

    // other options
    int verbose;
    int minBionanoGapLen;
    int maxTrimLength;
    int maxLowCompTrimPerc;
    int trimOffset;
} TrimContext;

static void trim_pre(PassContext* pctx, TrimContext* tctx)
{
    if(tctx->verbose)
    {
	    printf( ANSI_COLOR_GREEN "PASS contig trimming\n" ANSI_COLOR_RESET);

	    printf( ANSI_COLOR_RED "OPTIONS\n" ANSI_COLOR_RESET);
	    printf( ANSI_COLOR_RED "  verbose %d\n" ANSI_COLOR_RESET, tctx->verbose);
	    printf( ANSI_COLOR_RED "  minBionanoGapLen %d\n" ANSI_COLOR_RESET, tctx->minBionanoGapLen);
	    printf( ANSI_COLOR_RED "  maxTrimLength %d\n" ANSI_COLOR_RESET, tctx->maxTrimLength);
	    printf( ANSI_COLOR_RED "  maxLowCompTrimPerc %d\n" ANSI_COLOR_RESET, tctx->maxLowCompTrimPerc);
	    printf( ANSI_COLOR_RED "  trimOffset %d\n" ANSI_COLOR_RESET, tctx->trimOffset);
      
        if(tctx->trackDust)
		    printf( ANSI_COLOR_RED "  dust Track %s\n" ANSI_COLOR_RESET, tctx->trackDust->name);
	    if(tctx->trackTan)
            printf( ANSI_COLOR_RED "  tandem Track %s\n" ANSI_COLOR_RESET, tctx->trackTan->name);
    }

	tctx->twidth = pctx->twidth;

}

static void trim_post(TrimContext* ctx)
{
    if(ctx->verbose)
    {
        if (ctx->statsTrimmedContigs > 0)
        {
            printf("#trimmed contigs %d\n", ctx->statsTrimmedContigs);
        }

        if (ctx->statsTrimmedBases > 0)
        {
            printf("#trimmed bases: %d\n", ctx->statsTrimmedBases);
        }
    }
}

static int trim_handler(void* _ctx, Overlap* ovl, int novl)
{
	TrimContext* ctx = (TrimContext*) _ctx;
	int i, j, k;


    // assumption: overlaps must be chained with LAfilterChains !!!
    // there can be up to two chains possible -> one chain at the end and one chain at the beginning of the contig!!!
    
    printf("trimHander Begin\n");
    printf(" r[%d, %d]\n", ovl->aread, ovl->bread);
    for(i=0; i<novl;i++)
    {
        Overlap *o = ovl + i;
        printf("   a[%d, %d] %c b[%d, %d]\n",o->path.abpos, o->path.aepos, (o->flags & OVL_COMP) ? 'c' : 'n', o->path.bbpos, o->path.bepos);        

    }

    printf("trimHander End\n");
    return 1;
}


static void usage()
{
	fprintf(stderr, "[-v] [-GTLO <int>] [-bsp <file>] [-dt <track>] <db> <overlaps_in> <contigs_out>\n");

	fprintf(stderr, "options: -v        verbose\n");
	fprintf(stderr, "         -d <trc>  low complexity track (e.g. dust)\n");
	fprintf(stderr, "         -t <trc>  tandem repeat track  (e,f, tan)\n");
	fprintf(stderr, "         -b <file> gap csv file i.e. based on bionano agp. Must be in following format: +|-<ContigID1> <gapLen> +|-<contigID2>\n");
	fprintf(stderr, "                   where: ContigID1, and ContigID2 must refer to the DAmar seqeuence database (i.e. 0-based contig IDs)\n");
    fprintf(stderr, "                   the mandatory +|- prefix of the ContigID describes the orientation of the contig\n");
    fprintf(stderr, "                   If a bionano-gap file is given, then only gaps up the minimum gaps size of (default: %d) are trimmed\n", MIN_BIONANO_GAP_SIZE);
    fprintf(stderr, "                   --> idea behind this: If Bionano inserts a 13bp gap, then it's most probable that the adjacent contigs overlap with each other\n");
    fprintf(stderr, "         -G <int>  min Bionano gap size (default: %d)\n", MIN_BIONANO_GAP_SIZE);
    fprintf(stderr, "         -T <int>  maximum trim length (default: -1)\n");
    fprintf(stderr, "         -L <int>  do not trim contigs if trim length contains more then -S bases (in %%) of tandem repeats (default: -1, valid range: [0,100])\n");
    fprintf(stderr, "         -s <file> some stats are written to the statistics file <file>\n");
    fprintf(stderr, "         -p <file> write seqences that were trimmed off to -p file\n");
    fprintf(stderr, "         -O <int>  trim offset in bases (default %d), i.e. in best case (if we have single overlap between 2 contigs) a gap of size 2xtrim_offset is created )\n", TRIM_OFFSET);
    fprintf(stderr, "                   in case a valid alignment chain consisting of multiple alignments is present (representing heterozygous variations). The first last and the last alignment are used, (- trimOffset and + trimOffset, accordingly) \n");
    fprintf(stderr, "                   (- trimOffset and + trimOffset, accordingly) creates a larger gap size, but heopefully removes the heterozygous difference.\n");
}


int main(int argc, char* argv[])
{
	HITS_DB db;
	TrimContext tctx;
	PassContext* pctx;
	
    FILE* fileOvlIn = NULL;
    FILE* fileGapsIn = NULL;
	FILE* fileCtgOut = NULL;
    FILE* filePurgeOut = NULL;
    FILE* fileStatsOut = NULL;

	bzero(&tctx, sizeof(TrimContext));

	tctx.db = &db;

// args

	char* pcTrackDust = NULL;
	char* pcTrackTan = NULL;

    char* pathInBionanoGapCSV = NULL;
	char* pathOutStats = NULL;
    char* pathOutPurged = NULL;

	char* pcPathReadsIn = NULL;
	char* pcPathOverlapsIn = NULL;
	char* pcPathContigsOut = NULL;

	int c,tmp;

    tctx.minBionanoGapLen = MIN_BIONANO_GAP_SIZE;
    tctx.maxTrimLength = -1;
    tctx.maxLowCompTrimPerc = -1;
    tctx.trimOffset = TRIM_OFFSET;

	opterr = 0;
	while ((c = getopt(argc, argv, "vd:t:b:G:T:L:s:p:O:")) != -1)
	{
		switch (c)
		{
			case 'v':
				tctx.verbose = 1;
				break;
            case 'd':
				pcTrackDust = optarg;
				break;
            case 't':
				pcTrackTan = optarg;
				break;            
            case 'b':
				pathInBionanoGapCSV = optarg;
				break;             
            case 'G':
				tctx.minBionanoGapLen = atoi(optarg);
				break;                           
            case 'T':
				tctx.maxTrimLength = atoi(optarg);
				break;            
            case 'L':
                tmp = atoi(optarg);
                if (tmp < 0 || tmp > 100)
                {
                    fprintf(stderr, "[ERROR] Invalid range for tandem repeat fraction %d. Must be in [0,100]\n", tmp);
                    exit(1);
                }
				tctx.maxLowCompTrimPerc = tmp;
				break;                                                           
            case 's':
				pathOutStats = optarg;
				break;            
            case 'p':
				pathOutPurged = optarg;
				break;            
            case 'O':
				tctx.trimOffset = atoi(optarg);
				break;            
			default:
				fprintf(stderr, "unknown option %c\n", optopt);
				usage();
				exit(1);
		}
	}

	if (argc - optind != 3)
	{
		usage();
		exit(1);
	}

	pcPathReadsIn = argv[optind++];
	pcPathOverlapsIn = argv[optind++];
	pcPathContigsOut = argv[optind++];

	if ((fileOvlIn = fopen(pcPathOverlapsIn, "r")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open %s\n", pcPathOverlapsIn);
		exit(1);
	}

	if ((fileCtgOut = fopen(pcPathContigsOut, "w")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open %s\n", pcPathContigsOut);
		exit(1);
	}

	if (Open_DB(pcPathReadsIn, &db))
	{
		fprintf(stderr, "[ERROR] - could not open %s\n", pcPathReadsIn);
		exit(1);
	}

    if (pcTrackDust)
    {
        tctx.trackDust = track_load(&db, pcTrackDust);
	    if (! tctx.trackDust)
		{
            fprintf(stderr, "[ERROR] - could not load track %s\n", pcTrackDust);
            exit(1);
        }
    }

    if (pcTrackTan)
    {
        tctx.trackTan = track_load(&db, pcTrackTan);
	    if (! tctx.trackTan)
		{
            fprintf(stderr, "[ERROR] - could not load track %s\n", pcTrackTan);
            exit(1);
        }
    }

	if (pathOutStats)
	{
		tctx.fileOutTrimmedStats = fopen(pathOutStats, "w");

		if (tctx.fileOutTrimmedStats == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathOutStats);
			exit(1);
		}		
	}

    if (pathOutPurged)
	{
		tctx.fileOutDicardedContigParts = fopen(pathOutPurged, "w");

		if (tctx.fileOutDicardedContigParts == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathOutPurged);
			exit(1);
		}
	}

    if (pathInBionanoGapCSV)
	{
		tctx.fileInBionanoGaps = fopen(pathInBionanoGapCSV, "r");

		if (tctx.fileInBionanoGaps == NULL)
		{
			fprintf(stderr, "could not open %s\n", pathInBionanoGapCSV);
			exit(1);
		}
	}



// passes

	pctx = pass_init(fileOvlIn, NULL);

	pctx->split_b = 1;
	pctx->load_trace = 1;
	pctx->unpack_trace = 1;
	pctx->data = &tctx;
	pctx->write_overlaps = 0;
	pctx->purge_discarded = 0;

	trim_pre(pctx, &tctx);

	pass(pctx, trim_handler);

	trim_post(&tctx);

	pass_free(pctx);

    //todo cleanup

	Close_DB(&db);
	fclose(fileOvlIn);

	return 0;
}
