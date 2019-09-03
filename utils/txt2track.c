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

int SORT;

static void usage(FILE* out, char *name)
{
	fprintf(out, "Usage:\n         %s [-vmr] [-s <int>] <db> <txt> <track>\n", name);
	fprintf(out, "Options:	-m       track is mask track and has the format readid pos pos\n");
	fprintf(out, "	        -r       track is read pair track and has the format: readID readID\n");
	fprintf(out, "	        -s <int> sort fields from 0. to -s <int> \n");
	fprintf(out, "	        -v       verbose mode\n");
}

static int cmp_ints(const void* x, const void* y)
{
	int* a = (int*) x;
	int* b = (int*) y;

	// compare areads
	int l = a[0];
	int r = b[0];

	// compare breads or abpos
	if (SORT == 2 && l == r)
	{
		l = a[1];
		r = b[1];
	}

	// compare abepos
	if (SORT > 2 && l == r)
	{
		l = a[2];
		r = b[2];
	}

	return l - r;
}

int main(int argc, char* argv[])
{
	HITS_DB db;

	// process arguments

	int c;

	opterr = 0;

	int verbose = 1;
	int isMask = 0;
	int isReadPair = 0;
	SORT = -1;

	while ((c = getopt(argc, argv, "vmrs:")) != -1)
	{
		switch (c)
		{
			case 'm':
				isMask = 1;
				break;

			case 'r':
				isReadPair = 1;
				break;
			case 's':
				SORT = atoi(optarg);
				break;
			case 'v':
                                verbose = 1; 
                                break;
			default:
				printf("Unknow option: %s\n", argv[optind]);
				usage(stderr, argv[0]);
				exit(1);
		}
	}

	if (argc - optind < 3)
	{
		usage(stderr, argv[0]);
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

	if(SORT>=0)
	{
		if(isMask && SORT > 2)
		{
			fprintf(stderr, "[WARNING]: Sort fields cannot be greater then 2, if isMask is set!\n");
			exit(1);
		}

		if(isReadPair && SORT > 1)
		{
			fprintf(stderr, "[WARNING]: Sort fields cannot be greater then 1, if isReadPair is set!\n");
			exit(1);
		}

		if(SORT > 2)
		{
			fprintf(stderr, "[WARNING]: Unsupported number for sort fields(%d). (isMask [0,1,2], isReadPair [0,1])\n", SORT);
				exit(1);
		}
	}

	int field1;
	int field2;
	int field3;

	int* ints = NULL;
	int nints = 0;
	int maxints = 0;

	// todo dynamically read number of fields from first lines
	// for now: if isMask: 3 fields, aread, abpos, aepos
	// 				  if isReadPair: 2 fields, aread, bread

	int count;
        int nLines=0;

	if(verbose)
	  printf("read file %s:", pathTxt);
	while (1)
	{
		nLines++;
		if(isMask)
		{
			count = fscanf(fileIn, "%d%*c%d%*c%d\n", &field1, &field2, &field3);
			if (count != 3)
			{
				if(count != EOF)
 					fprintf(stderr, "[WARNING] - Line %d does not have proper mask format. Count is %d. Stop!\n", nLines, count);
				break;
			}
		}
		else
		{
			count = fscanf(fileIn, "%d%*c%d\n", &field1, &field2);
			if (count != 2)
			{
				if(count !=EOF)
					fprintf(stderr, "[WARNING] - Line %d does not have proper read pair format. Count is %d. Stop!\n", nLines, count);
				break;
			}
		}

		if (nints + count >= maxints)
		{
			maxints = 1.2 * maxints + 1000;
			ints = realloc(ints, maxints * sizeof(int));
		}

		ints[nints + 0] = field1;
		ints[nints + 1] = field2;
		if(isMask)
			ints[nints + 2] = field3;

		int nReads = DB_NREADS(&db);
		if (isMask)
		{
			int rlen = DB_READ_LEN(&db, field1);

			if (field1 < 0 || field1 >= nReads || field3 > rlen || field2 > field3 || field3 > rlen)
			{
				printf("interval out of bounds %d %d..%d (%d)\n", field1, field2, field3, rlen);
			}
		}
		else
		{
			if (field1 < 0 || field1 >= nReads || field2 < 0 || field2 >= nReads)
			{
				printf("read ids %d, %d are out of DB read range [%d, %d] \n", field1, field2, 0, nReads);
			}
		}

		nints += count;
	}
        if(verbose)
		printf(" found %d lines\n", nLines);

	int fields = (isMask) ? 3 : 2;
	if (SORT>=0)
	{
		if(verbose)
			printf("Sort by: %d columns\n", SORT);
		qsort(ints, nints / fields, sizeof(int) * fields, cmp_ints);
	}
	track_anno* anno = malloc(sizeof(track_anno) * (DB_NREADS(&db) + 1));
	track_data* data = NULL;
	int ndata = 0;
	int maxdata = 0;

	bzero(anno, sizeof(track_anno) * (DB_NREADS(&db) + 1));
	int i;
	for (i = 0; i < nints; i += fields)
	{
		if (ndata + fields >= maxdata)
		{
			maxdata = maxdata * 1.2 + 1000;
			data = realloc(data, sizeof(track_data) * maxdata);
		}

		anno[ints[i]] += (fields - 1) * sizeof(track_data);
		data[ndata] = ints[i + 1];
		if(isMask)
			data[ndata + 1] = ints[i + 2];

		ndata += fields - 1;
	}

	track_anno coff, off;
	off = 0;

	for (i = 0; i <= DB_NREADS(&db); i++)
	{
		coff = anno[i];
		anno[i] = off;
		off += coff;
	}

	if(verbose)
		printf("write track: %s\n", nameTrack);
	track_write(&db, nameTrack, 0, anno, data, ndata);

	free(anno);
	free(data);
	free(ints);
	Close_DB(&db);

	if(verbose)
		printf("DONE\n");
	return 0;
}
