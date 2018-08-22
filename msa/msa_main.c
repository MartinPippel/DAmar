
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <assert.h>
#include <unistd.h>

#include "lib/oflags.h"
#include "lib/tracks.h"
#include "lib/pass.h"
#include "lib/colors.h"
#include "msa.h"

#include "dalign/DB.h"
#include "dalign/align.h"

// development only

#undef DEBUG
#undef DEBUG_OFFSETS
#define VERBOSE

typedef struct
{
    msa* m;

    int rid;
    int maxovls;
    
    int twidth;
    
    int nreads;
    char** reads;
    
    char* base_out;
    
    FILE* fileOutMsa;
    FILE* fileOutConsensus;
    HITS_DB* db;
    
    // re-alignment
    
    Work_Data* align_work_data;
    Align_Spec* align_spec;
    
} MsaContext;


extern char *optarg;
extern int optind, opterr, optopt;


static int round_down(int n, int f)
{
    return n - (n % f);
}

static int round_up(int n, int f)
{
    return (n + f - 1) - ((n-1) % f);
}

static void load_reads(MsaContext* mctx, Overlap* pOvls, int nOvls)
{
    if (nOvls >= mctx->nreads)
    {
        int nreads = mctx->nreads * 1.2 + nOvls + 1;
        mctx->reads = (char**)realloc(mctx->reads, sizeof(char*)*nreads);
         
        for (;mctx->nreads < nreads; mctx->nreads++)
        {
            mctx->reads[ mctx->nreads ] = New_Read_Buffer(mctx->db);
        }
    }

    Load_Read(mctx->db, pOvls->aread, mctx->reads[0], 0);

    int i;
    for (i = 0; i < nOvls; i++)
    {
        Load_Read(mctx->db, pOvls[i].bread, mctx->reads[i+1], 0);
        
        if (pOvls[i].flags & OVL_COMP)
        {
            Complement_Seq( mctx->reads[i+1] );
        }
    }
}

static void adjust_offsets(MsaContext* mctx, Overlap* ovl, int novl)
{
    Alignment aln;
    Path path;

    int twidth = mctx->twidth;
        
    aln.path = &path;
    aln.aseq = mctx->reads[0];
    aln.alen = ovl->alen;    

    int ntmax = 1000;
    int* ntrace = (int*)malloc(sizeof(int)*ntmax);
    int ntlen;
    
    int i;    
    for (i = 0; i < novl; i++)
    {
#ifdef DEBUG_OFFSETS
        printf("OVL %5d x %5d @ %5d..%5d x %5d..%5d\n", ovl[i].alen, ovl[i].blen, ovl[i].path.abpos, ovl[i].path.aepos, ovl[i].path.bbpos, ovl[i].path.bepos);
#endif

        aln.bseq = mctx->reads[i+1];
        aln.blen = ovl[i].blen;
        
        path = ovl[i].path;
        
        Compute_Trace_MID(&aln, mctx->align_work_data, twidth);
        
        // Print_Reference(stdout, &aln, mctx->align_work_data, 0, 100, 0, 0, 5);
        
        {
            int a = ovl[i].path.abpos;
            int b = ovl[i].path.bbpos;
            int p, t;
            int diffs = 0;
            int matches = 0;

            int ntcur = 0;
            ntlen = (round_down(ovl[i].path.aepos - 1, twidth) - round_up(ovl[i].path.abpos + 1, twidth) + twidth) / twidth * 2 + 2;
        
            if (ntlen > ntmax)
            {
                ntmax = 1.2 * ntmax + ntlen;
                ntrace = (int*)malloc(sizeof(int) * ntmax);
            }

            // printf("pts  ");

            int bprev = aln.path->bbpos;

            for (t = 0; t < aln.path->tlen; t++)
            {
                if ((p = ((int*)(aln.path->trace))[t]) < 0)
                { 
                    p = -p - 1;
                    while (a < p)
                    { 
                        if (aln.aseq[a] != aln.bseq[b]) diffs++;
                        else matches++;
                
                        a += 1;
                        b += 1;

                        if (a % twidth == 0) 
                        {
                            ntrace[ntcur++] = diffs;
                            ntrace[ntcur++] = b - bprev;
                            bprev = b;
                    
                            // printf(" %4dx%4d %3d %3d", a, b, diffs, matches);
                            diffs = matches = 0;
                        }
                    }

                    diffs++;
                    b += 1;
                }
                else
                { 
                    p--;
                
                    while (b < p)
                    { 
                        if (aln.aseq[a] != aln.bseq[b]) diffs++;
                        else matches++;

                        a += 1;
                        b += 1;

                        if (a % twidth == 0) 
                        {
                            ntrace[ntcur++] = diffs;
                            ntrace[ntcur++] = b - bprev;
                            bprev = b;

                            // printf(" %4dx%4d %3d %3d", a, b, diffs, matches);
                            diffs = matches = 0;
                        }
                    }
                
                    diffs++;
                    a += 1;

                    if (a % twidth == 0) 
                    {
                        ntrace[ntcur++] = diffs;
                        ntrace[ntcur++] = b - bprev;
                        bprev = b;

                        // printf(" %4dx%4d %3d %3d", a, b, diffs, matches);
                        diffs = matches = 0;
                    }
                }
            }

            p = aln.path->aepos;
            while (a < p)
            { 
                if (aln.aseq[a] != aln.bseq[b]) diffs++;
                else matches++;

                a += 1;
                b += 1;

                if (a % twidth == 0 && a != ovl[i].path.aepos) 
                {
                    ntrace[ntcur++] = diffs;
                    ntrace[ntcur++] = b - bprev;
                    bprev = b;

                    // printf(" %4dx%4d %3d %3d", a, b, diffs, matches);
                    diffs = matches = 0;
                }
            }

            ntrace[ntcur++] = diffs;
            ntrace[ntcur++] = b - bprev;
        
            /*
            if (ntcur != ntlen)
            {
                a = round_up(ovl[i].path.abpos+1,twidth);
                b = ovl[i].path.bbpos;
                printf("ntcur %d != ntlen %d\n", ntcur, ntlen);
                printf("%d..%d x %d..%d -> %d pts\n", ovl[i].path.abpos, ovl[i].path.aepos, ovl[i].path.bbpos, ovl[i].path.bepos, ntlen);

                printf("npts  ");
                for (t = 0; t < ntcur-1; t+=2)
                {
                    b += ntrace[t+1];
                    printf(" %d@%d/%dx%d", ntrace[t], ntrace[t+1],a,b);
                
                    a += twidth;
                }
                printf(" %d@end", ntrace[ntcur-1]);
            
                printf("\n");
            }
            */
        }


        int j;

#ifdef DEBUG_OFFSETS
        printf("(%3d)", ovl[i].path.tlen);
        for (j = 1; j < ovl[i].path.tlen; j += 2)
        {
            if ( j > ntlen || j > ovl[i].path.tlen || ((uint16*)(ovl[i].path.trace))[j] != ntrace[j] )
            {
                printf(ANSI_COLOR_RED " %3d" ANSI_COLOR_RESET, ((uint16*)(ovl[i].path.trace))[j]);
            }
            else
            {
                printf(" %3d", ((uint16*)(ovl[i].path.trace))[j]);
            }
        }
        printf("\n");

        printf("(%3d)", ntlen);
        for (j = 1; j < ntlen; j += 2)
        {
            printf(" %3d", ntrace[j]);
        }
        printf("\n\n");
#endif

        assert(ovl[i].path.tlen == ntlen);

        for (j = 0; j < ntlen; j++)
        {
            ((uint16*)ovl[i].path.trace)[j] = ntrace[j];
        }
    }

}

static void write_seq(FILE* file, char* seq)
{
    const int width = 80;
    int len = strlen(seq);
    int j;
        
    for (j = 0; j + width < len; j += width)
    {
        fprintf(file, "%.*s\n", width, seq + j);
    }
    
    if (j < len)
    {
        fprintf(file, "%s\n", seq + j);
    }
}

static void pre_msa(PassContext* pctx, MsaContext* mctx)
{
    char* fname = (char*)malloc( strlen(mctx->base_out) + 100);
    sprintf(fname, "%s.cons.fa", mctx->base_out);
    
    mctx->fileOutConsensus = fopen(fname, "w");
    
    sprintf(fname, "%s.msa", mctx->base_out);
    
    mctx->fileOutMsa = fopen(fname, "w");
    
    free(fname);
    
    mctx->twidth = pctx->twidth;
}

static void post_msa(MsaContext* mctx)
{
    fclose(mctx->fileOutConsensus);
    fclose(mctx->fileOutMsa);
}

static int handler_msa(void* _ctx, Overlap* pOvls, int nOvls)
{
    MsaContext* ctx = _ctx;

    int a = pOvls->aread;
    int alen = pOvls->alen;

    if (ctx->rid != -1 && ctx->rid != a)
    {
        return 1;
    }

    printf("read %d (%d) ... %d overlaps,", a, alen, nOvls);

    if (nOvls > ctx->maxovls) 
    {
        nOvls = ctx->maxovls;
    }
    
    printf(" using %d\n", nOvls);
    
    msa* m;
    m = msa_init();
    
    m->twidth = ctx->twidth;
    
    load_reads(ctx, pOvls, nOvls);

    msa_add(m, ctx->reads[0], 0, alen, 0, alen, NULL, 0);

    int j;
    for (j = 0; j < nOvls; j++)
    {
        msa_add(m, ctx->reads[j+1], 
                pOvls[j].path.abpos, pOvls[j].path.aepos, 
                pOvls[j].path.bbpos, pOvls[j].path.bepos, 
                pOvls[j].path.trace, pOvls[j].path.tlen);
                
        // msa_print(m, stdout, 1);
        // msa_print_profile(m, stdout, 1);
        // printf("\n");
    }
    
    msa_print_v(m, ctx->fileOutMsa);
    // msa_print_profile(m, stdout, 1);
    
    char* cons = msa_consensus(m, 0);
    write_seq(ctx->fileOutConsensus, cons);
    
    msa_free(m);
    
    return ( ctx->rid == -1 );
}

static void usage()
{
    printf("usage: <db> <overlaps> <base_out> <read.id> <max.ovls>\n");
}

int main(int argc, char* argv[])
{
    FILE* fileOvls;
    HITS_DB db;
    PassContext* pctx;
    MsaContext mctx;
    
    bzero(&mctx, sizeof(MsaContext));

    mctx.db = &db;

    if (argc != 6)
    {
        usage();
        exit(1);
    }

    char* pcPathReadsIn = argv[1];
    char* pcPathOverlaps = argv[2];
    mctx.base_out = argv[3];
    mctx.rid = atoi(argv[4]);
    mctx.maxovls = atoi(argv[5]);
    
    if (mctx.maxovls < 1)
    {
        mctx.maxovls = INT_MAX;
    }

    if ( (fileOvls = fopen(pcPathOverlaps, "r")) == NULL )
    {
        fprintf(stderr, "could not open '%s'\n", pcPathOverlaps);
        exit(1);
    }

    if (Open_DB(pcPathReadsIn, &db))
    {
        fprintf(stderr, "could not open '%s'\n", pcPathReadsIn);
        exit(1);
    }
    
    pctx = pass_init(fileOvls, NULL);
    
    pctx->split_b = 0;
    pctx->load_trace = 1;
    pctx->unpack_trace = 1;
    pctx->data = &mctx;
    
    pre_msa(pctx, &mctx);
    
    Trim_DB(&db);

    pass(pctx, handler_msa);
    
    post_msa(&mctx);
    
    Close_DB(&db);

    return 0;
}
