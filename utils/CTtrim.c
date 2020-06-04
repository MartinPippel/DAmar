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
#define FUZZY_BASES 1500

#define DEBUG_MASKING

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

    // vector to store overlapping contigs: length: #contigs x #contigs
    int *LAStrimMatrix;

    // other options
    int verbose;
    int minBionanoGapLen;
    int maxTrimLength;
    int maxLowCompTrimPerc;
    int trimOffset;
    int maxFuzzyBases;
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
        printf( ANSI_COLOR_RED "  maxFuzzyBases %d\n" ANSI_COLOR_RESET, tctx->maxFuzzyBases);
      
        if(tctx->trackDust)
		    printf( ANSI_COLOR_RED "  dust Track %s\n" ANSI_COLOR_RESET, tctx->trackDust->name);
	    if(tctx->trackTan)
            printf( ANSI_COLOR_RED "  tandem Track %s\n" ANSI_COLOR_RESET, tctx->trackTan->name);
    }

	tctx->twidth = pctx->twidth;
    tctx->LAStrimMatrix = (int*)malloc(DB_NREADS(tctx->db)*sizeof(int)*DB_NREADS(tctx->db));
    assert(tctx->LAStrimMatrix != NULL);
    bzero(tctx->LAStrimMatrix,DB_NREADS(tctx->db)*sizeof(int)*DB_NREADS(tctx->db));
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
    free(ctx->LAStrimMatrix);
}

static int getTrimPositions(TrimContext *ctx, Overlap *ovl, int pointA, int* cutA, int *cutB)
{
    int abeg = ovl->path.abpos;
	int aend = ovl->path.aepos;

	int bbeg = ovl->path.bbpos;
	int bend = ovl->path.bepos;

	int twidth = ctx->twidth;

    if(pointA < abeg || pointA > aend)
        return 1;

	//printf("getTrimPositions %d x %d, a(%d, %d) %c b(%d, %d) pointA: %d\n", ovl->aread, ovl->bread, abeg, aend, (ovl->flags & OVL_COMP) ? 'C' : 'N', bbeg, bend, pointA);
	
    int dist = pointA - abeg;
    int apos, bpos;
    if (ovl->path.tlen)
    {		        	
        ovl_trace* trace = ovl->path.trace;
        apos = abeg;
        bpos = ovl->path.bbpos;

        int j = 0;
        while (j < ovl->path.tlen)
        {
      //      printf("apos %6d, bpos %6d, oldDist %6d, newDist %6d\n", apos, bpos, dist, abs(pointA-((apos / twidth + 1) * twidth)));
            if(dist < abs(pointA - ((apos / twidth + 1) * twidth)))
                break;
            apos = (apos / twidth + 1) * twidth;
            bpos += trace[j + 1];
            //printf("apos %6d, bpos %6d\n", apos, bpos);
            dist = abs(pointA - apos); 
            j += 2;
        }
    }
    else
    {
        apos = pointA;
        bpos = bbeg + dist;  
    }
	
  //  printf("apos: %d, bpos: %d\n", apos, bpos);
    
    *cutA = apos;
    *cutB = bpos;

    if(*cutA < abeg || *cutA > aend)
        return 1;

    if(*cutB < bbeg || *cutB > bend)
        return 1;

	//printf("final range: %d, %d\n", *cutA, *cutB);
    return 0;
}

static int analyzeContigOverlaps(TrimContext *ctx, Overlap *ovl, int novl)
{
    int i;

    int aLen = DB_READ_LEN(ctx->db, ovl->aread);
    int bLen = DB_READ_LEN(ctx->db, ovl->bread);

    // assumption: input overlaps must be chained with LAfilterChains !!!
    // one chain at the end and one chain at the beginning of a contig are possible!!!
    
    // sanity check
    Overlap *o1 = ovl;
    int cutA = -1;
    int cutB = -1;

    if(novl == 1 )
    {
        if((o1->path.abpos > ctx->maxFuzzyBases && o1->path.aepos < aLen - ctx->maxFuzzyBases) || (o1->path.bbpos > ctx->maxFuzzyBases && o1->path.bepos < bLen - ctx->maxFuzzyBases))
        {
            if(ctx->verbose)
            {
                printf("[WARNGING] fuzzy base check failed! Ignore invalid chain [%d, %d] a[%d,%d] %c b[%d,%d]!\n",o1->aread, o1->bread,o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n',o1->path.bbpos, o1->path.bepos);
            }
            return 1;
        }

        int pointA =  o1->path.abpos + (o1->path.aepos - o1->path.abpos)/2;

        if(getTrimPositions(ctx, o1, pointA, &cutA, &cutB))
        {
            printf("Unable to get cutPosition for OVL [%d,%d] a[%d,%d] %c b[%d,%d] and pointA: %d\n", o1->aread, o1->bread, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n',o1->path.bbpos, o1->path.bepos, pointA);
             return 1;
        }   

        assert((cutA - ctx->trimOffset > 0) && (cutA + ctx->trimOffset < aLen));
        assert((cutB - ctx->trimOffset > 0) && (cutB + ctx->trimOffset < bLen));

        // set cut position of contig_A
        if(o1->path.abpos < aLen - o1->path.aepos) // trim off contig at begin 
        {
            ctx->LAStrimMatrix[o1->aread*DB_NREADS(ctx->db)+o1->bread] = -(cutA + ctx->trimOffset);
        }
        else if(o1->path.abpos > aLen - o1->path.aepos) // trim off contig at end 
        {
            ctx->LAStrimMatrix[o1->aread*DB_NREADS(ctx->db)+o1->bread] = cutA - ctx->trimOffset;
        }
        else // containment 
        {
            printf("Contained overlap: [%d,%d] a[%d,%d] %c b[%d,%d] and pointA: %d\n", o1->aread, o1->bread, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n',o1->path.bbpos, o1->path.bepos, pointA);
            return 1;
        }
        
        // set cut position of contig_B
        if(o1->path.bbpos < bLen - o1->path.bepos) // trim off contig at begin 
        {
            if(o1->flags & OVL_COMP)
            {
                ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = bLen - (cutB + ctx->trimOffset);  
            }
            else
            {
                ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = -(cutB + ctx->trimOffset);    
            }                        
        }
        else if(o1->path.bbpos > bLen - o1->path.bepos) // trim off contig at end 
        {
            if(o1->flags & OVL_COMP)
            {
                ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = -(bLen - (cutB - ctx->trimOffset));
            }
            else 
            {
                ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = cutB - ctx->trimOffset;
            }            
        }    
    }
    else
    {
        int validChain=1;
        Overlap *o2;
        // first: sanity check for LAS chain
        for(i=1; i<novl;i++)
        {
            o2 = ovl + i;
            if (abs(o1->path.aepos - o2->path.abpos) > ctx->maxFuzzyBases || ((o1->flags & OVL_COMP) != (o2->flags & OVL_COMP))) 
            {
                validChain = 0;
                break;
            }
            o1=o2;
        }

        if (!validChain)
        {
            printf("INVALID chain: %d vs %d\n", ovl->aread, ovl->bread);
   
            for(i=0; i<novl;i++)
            {
                printf("   a[%d,%d] %c b[%d,%d]\n", ovl[i].path.abpos, ovl[i].path.aepos, (ovl[i].flags & OVL_COMP) ? 'c' : 'n',ovl[i].path.bbpos, ovl[i].path.bepos);
   
            }
            return 1;
        }

        o1 = ovl;
        o2 = ovl + (novl-1);

        // set cut position of contig_A
        if(o1->path.abpos < aLen - o2->path.aepos) // trim off contig at begin 
        {
            if(o2->path.abpos >= o1->path.aepos)
            {
                cutA = o2->path.abpos + ctx->trimOffset;
            }       
            else 
            {
                cutA = o1->path.aepos + ctx->trimOffset;
            }
            ctx->LAStrimMatrix[o1->aread*DB_NREADS(ctx->db)+o1->bread] = -(cutA);
        }
        else if(o1->path.abpos > aLen - o1->path.aepos) // trim off contig at end 
        {
            if(o2->path.abpos >= o1->path.aepos)
            {
                cutA = o1->path.aepos - ctx->trimOffset;
            }       
            else 
            {
                cutA = o2->path.abpos - ctx->trimOffset;
            }
            ctx->LAStrimMatrix[o1->aread*DB_NREADS(ctx->db)+o1->bread] = cutA;            
        }
        else // containment 
        {
            printf("Contained overlap: [%d,%d] a[%d,%d] %c b[%d,%d]\n", o1->aread, o1->bread, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n',o1->path.bbpos, o1->path.bepos);
            return 1;
        }
        // set cut position of contig_B

        if(o1->path.bbpos < bLen - o2->path.bepos) // trim off contig at begin 
        {
            if(o2->path.bbpos >= o1->path.bepos)
            {
                cutB = o2->path.bbpos + ctx->trimOffset;
                if(o1->flags & OVL_COMP)
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = bLen - cutB;
                }
                else
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = -(cutB);
                }                
            }
            else 
            {
                cutB = o1->path.bepos + ctx->trimOffset;
                if(o1->flags & OVL_COMP)
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = bLen - cutB;
                }
                else
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = -(cutB);
                }                
            }
        }
        else if(o1->path.bbpos > bLen - o2->path.bepos) // trim off contig at end 
        {
            if(o2->path.bbpos >= o1->path.bepos)
            {
                cutB = o1->path.bepos - ctx->trimOffset;
                if(o1->flags & OVL_COMP)
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = -(bLen - cutB);
                }
                else 
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = cutB;   
                }
            }
            else
            {
                cutB = o2->path.bbpos - ctx->trimOffset;   
                if(o1->flags & OVL_COMP)
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = -(bLen - cutB);
                }
                else 
                {
                    ctx->LAStrimMatrix[o1->bread*DB_NREADS(ctx->db)+o1->aread] = cutB;   
                }

            }
        }
        else // containment 
        {
            printf("Contained overlap: [%d,%d] a[%d,%d] %c b[%d,%d]\n", o1->aread, o1->bread, o1->path.abpos, o1->path.aepos, (o1->flags & OVL_COMP) ? 'c' : 'n',o1->path.bbpos, o1->path.bepos);
            return 1;
        }
    }
    return 0;
}

static int trim_handler(void* _ctx, Overlap* ovl, int novl)
{
	TrimContext* ctx = (TrimContext*) _ctx;
	int i, j;

    // analyze overlaps and find contig trim position 
    analyzeContigOverlaps(ctx, ovl, novl);
                    
    return 1;
}

static int getMaskedBases(TrimContext * ctx, HITS_TRACK * t, int contigID, int beg, int end)
{
#ifdef DEBUG_MASKING
    printf("call getMaskedBases on track %s, contigID: %d, in: [%d, %d] ",t->name, contigID, beg, end);
#endif
    if(t == NULL)
    {
#ifdef DEBUG_MASKING        
        printf(" --> masked bases 0 (track is Null)\n");
#endif
        return 0;
    }

	track_anno* mask_anno = t->anno;
	track_data* mask_data = t->data;

    if (contigID < 0 || contigID >= DB_NREADS(ctx->db))
    {
        fprintf(stderr, "[ERROR] - getMaskedBases contigID: %d out of bounds [0, %d]\n", contigID, DB_NREADS(ctx->db) - 1);
        fflush(stderr);
        exit(1);
    }

	track_anno rb, re;

	int maskBases = 0;
	int rBeg, rEnd;

	// repeat bases in a-read
	rb = mask_anno[contigID] / sizeof(track_data);
	re = mask_anno[contigID + 1] / sizeof(track_data);

	while (rb < re)
	{
		rBeg = mask_data[rb];
		rEnd = mask_data[rb + 1];

		maskBases += intersect(beg, end, rBeg, rEnd);

		rb += 2;
	}

#ifdef DEBUG_MASKING
    printf(" --> masked bases %d\m", maskBases);
#endif

	return maskBases;
}

static void usage()
{
	fprintf(stderr, "[-v] [-GTLOF <int>] [-bsp <file>] [-dt <track>] <db> <overlaps_in> <contigs_out>\n");

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
    fprintf(stderr, "         -F <int>  number of fuzzy bases (default: %d)\n", FUZZY_BASES);    
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
    tctx.maxFuzzyBases = FUZZY_BASES;

	opterr = 0;
	while ((c = getopt(argc, argv, "vd:t:b:G:T:L:s:p:O:F:")) != -1)
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
            case 'F':
				tctx.maxFuzzyBases = atoi(optarg);
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


        // debug report trim positions
    int nContigs = DB_NREADS(&db);
    int i,j;
    for(i=0; i<nContigs; i++)
    {
        int maxBeg=0;
        int minEnd=DB_READ_LEN(&db,i);
        for(j=0; j<nContigs; j++)
        {
            int cutPos= tctx.LAStrimMatrix[i*nContigs+j];
            if(cutPos < 0 && abs(cutPos) > maxBeg)
                maxBeg = abs(cutPos);
            if(cutPos > 0 && cutPos < minEnd)
                minEnd = cutPos;
            if(cutPos != 0)
                printf("FOUND CONTIG TRIM POSITION: CONTIG %d; TRIM: %d, TRIMLEN (%d) (OVL with: %d)\n", i, cutPos, (cutPos < 0) ? abs(cutPos) : DB_READ_LEN(&db,i)-cutPos, j);
        }
        if(maxBeg > 0 || minEnd != DB_READ_LEN(&db,i))
        {
            float dustBegFract,dustEndFract,tanBegFract,tanEndFract;
            dustBegFract = dustEndFract = tanBegFract = tanEndFract = 0.0;
            if (maxBeg > 0)
            {   
                dustBegFract = getMaskedBases(&tctx, tctx.trackDust, i, 0, maxBeg)*100.0/maxBeg;
                tanBegFract = getMaskedBases(&tctx, tctx.trackTan, i, 0, maxBeg)*100.0/maxBeg;
            }
            if(minEnd != DB_READ_LEN(&db,i))
            {   
                dustEndFract = getMaskedBases(&tctx, tctx.trackDust, i, minEnd, DB_READ_LEN(&db,i))*100.0/(DB_READ_LEN(&db,i)-minEnd);
                tanEndFract = getMaskedBases(&tctx, tctx.trackTan, i, minEnd, DB_READ_LEN(&db,i))*100.0/(DB_READ_LEN(&db,i)-minEnd);                
            }

            printf(" --> final trim Interval: [%d, %d] -> trimmed [%d, %d] dustFract(in %%) [%.2f] tanFract(in %%) [%.2f]\n", maxBeg, minEnd, maxBeg, DB_READ_LEN(&db,i)-minEnd, dustBegFract, dustEndFract, tanBegFract, tanEndFract);
        }
    }

	trim_post(&tctx);

	pass_free(pctx);

    //todo cleanup

	Close_DB(&db);
	fclose(fileOvlIn);

	return 0;
}
