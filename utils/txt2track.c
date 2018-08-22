/*
 * txt2track.c
 *
 *  Created on: 11 Jan 2017
 *      Author: pippelmn
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>

#include <sys/param.h>

#include "lib/colors.h"
#include "lib/oflags.h"
#include "lib/pass.h"
#include "lib/tracks.h"
#include "lib/utils.h"
#include "lib/lasidx.h"

#include "db/DB.h"
#include "dalign/align.h"

// getopt

extern char* optarg;
extern int optind, opterr, optopt;

static void usage()
{
    printf("<db> <txt> <track>\n");
}

static int cmp_ints(const void* x, const void* y)
{
    int* a = (int*)x;
    int* b = (int*)y;

    // compare areads
    int l = a[0];
    int r = b[0];

    // compare breads
    if (l == r)
    {
        l = a[1];
        r = b[1];
    }
/*
    // compare abpos
    if (l == r)
    {
        l = a[2];
        r = b[2];
    }
*/
    return l - r;
}

int main(int argc, char* argv[])
{
    HITS_DB db;

    // process arguments

    int c;

    opterr = 0;

    while ((c = getopt(argc, argv, "")) != -1)
    {
        switch (c)
        {

            default:
                printf("Unknow option: %s\n", argv[optind - 1]);
                usage();
                exit(1);
        }
    }

    if (argc - optind < 3)
    {
        usage();
        exit(1);
    }

    char* pathDb = argv[optind++];
    char* pathTxt = argv[optind++];
    char* nameTrack = argv[optind++];

    if (Open_DB(pathDb, &db))
    {
        printf("could not open '%s'\n", pathDb);
        return 1;
    }

    FILE* fileIn;

    if (strcmp(pathTxt, "-") == 0)
    {
        fileIn = stdin;
    }
    else
    {
        fileIn = fopen(pathTxt, "r");

        if (fileIn == NULL)
        {
            fprintf(stderr, "failed to open '%s'\n", pathTxt);
            exit(1);
        }
    }

    int aRead;
    int beg;
    int end;

    int* ints = NULL;
    int nints = 0;
    int maxints = 0;

    // todo dynamically read number of fields from first lines
    // for now: 3 columns, aread, abpos, aepos

    while ( 1 )
    {
        int count = fscanf(fileIn, "%d %d\n", &aRead, &beg);//, &end);

        if (count != 2)
        {
            break;
        }

        if (nints + 2 >= maxints)
        {
            maxints = 1.2 * maxints + 1000;
            ints = realloc(ints, maxints * sizeof(int));
        }

        ints[nints + 0] = aRead;
        ints[nints + 1] = beg;
        //ints[nints + 2] = end;



/**
        int rlen = DB_READ_LEN(&db, aRead);

        if ( beg < 0 || beg > rlen || end < beg || end > rlen )
        {
            printf("interval out of bounds %d %d..%d (%d)\n", aRead, beg, end, rlen);
        }
**/

        nints += 2;
    }

    qsort(ints, nints/2, sizeof(int) * 2, cmp_ints);

    track_anno* anno = malloc( sizeof(track_anno) * (DB_NREADS(&db) + 1) );
    track_data* data = NULL;
    int ndata = 0;
    int maxdata = 0;

    bzero(anno, sizeof(track_anno) * (DB_NREADS(&db) + 1));

    int i;
    for (i = 0; i < nints; i += 2)
    {
        if (ndata + 1 >= maxdata)
        {
            maxdata = maxdata * 1.2 + 1000;
            data = realloc( data, sizeof(track_data) * maxdata );
        }

        anno[ ints[i + 0] ] += 1 * sizeof(track_data);
        data[ ndata + 0 ] = ints[ i + 1 ];
        //data[ ndata + 1 ] = ints[ i + 2 ];

        ndata += 1;
    }

    track_anno coff, off;
    off = 0;

    for (i = 0; i <= DB_NREADS(&db); i++)
    {
        coff = anno[i];
        anno[i] = off;
        off += coff;
    }

    track_write(&db, nameTrack, 0, anno, data, ndata);

    free(anno);
    free(data);

    Close_DB(&db);

    return 0;
}




