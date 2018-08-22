/*********************************************************************************************\
 *
 *  Find all local self-alignment between long, noisy DNA reads:
 *    Compare sequences in each supplied blocks against themselves search for local alignments
 *    of MIN_OVERLAP or more above the diagonal (A start coord > B start coord).  An output
 *    stream of 'Overlap' records (see align.h) is written in binary to a set of files, one
 *    per thread, each encoding a given found local alignment.  The -v option turns on a verbose
 *    reporting mode that reports progress and gives statistics on each major stage.
 *
 *    The filter operates by looking for a pair of diagonal bands of width 2^'s' that contain
 *    a collection of exact matching 'k'-mers between positions of a sequence, such that the total
 *    number of bases covered by 'k'-mer hits is 'h'.  k cannot be larger than 15 in the
 *    current implementation.
 *
 *    For each subject, say XXX, the program outputs a file containing LAs of the form
 *    XXX.XXX.T#.las where # is the thread that detected and wrote out the collection of LAs.
 *    For example, if NTHREAD in the program is 4, then 4 files are output for each subject block.
 *
 *  Author:  Gene Myers
 *  Date  :  March 27, 2016
 *
 *********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>

#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#include "lib/oflags.h"
#include "lib/colors.h"
#include "lib/utils.h"
#include "lib/tracks.h"

#include "db/DB.h"
#include "tandem.h"

extern char *optarg;
extern int optind, opterr, optopt;

#define DEF_ARG_O "tan"

static void usage()
{
	fprintf(stderr, "usage:  \n\n");
	fprintf(stderr, "datander [-v] [-kwhjls<int>] [-e<double(.70)] [-o<string(tan)] <subject:db|dam>\n\n");
	fprintf(stderr, "options: -v ... verbose\n");
	fprintf(stderr, "         -k ... kmer length (defaults: raw pacbio reads: 12)\n");
	fprintf(stderr, "         -w ... diagonal band width (default: 4)\n");
	fprintf(stderr, "         -h ... hit theshold (in bp.s, default: 35\n");
	fprintf(stderr, "         -j ... number of threads (default: 4). Must be a power of 2!\n");
	fprintf(stderr, "         -e ... average correlation rate for local alignments (default: 0.7). Must be in [.5,1.)\n");
	fprintf(stderr, "         -l ... minimum alignment length (defaul: 500)\n");
	fprintf(stderr, "         -s ... trace spacing, i.e. report trace point every -s base pairs (default: raw pacbio reads: 100)\n");
	fprintf(stderr, "         -o ... output folder (default: %s)\n", DEF_ARG_O);
}

int VERBOSE;   //   Globally visible to tandem.c
int MINOVER;

static int read_DB(HITS_DB *block, char *name, int kmer)
{
	int i, isdam;

	isdam = Open_DB(name, block);
	if (isdam < 0)
	{
		exit(1);
	}

	{
		for (i = 0; i < block->nreads; i++)
			if (block->reads[i].rlen < kmer)
			{
				fprintf(stderr, "%s: Block %s contains reads < %dbp long !  Run DBsplit.\n", Prog_Name, name, kmer);
				exit(1);
			}
	}

	Read_All_Sequences(block, 0);

	return (isdam);
}

static void createSubdir(char *out)
{
	struct stat s;

	int err = stat(out, &s);

	if (err == -1)
	{
		if (errno == ENOENT)
		{
			mkdir(out, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
		}
		else
		{
			fprintf(stderr, "Cannot create output directory: %s\n", out);
			exit(1);
		}
	}
	else
	{
		if (!S_ISDIR(s.st_mode))
		{
			fprintf(stderr, "Output directory name: \"%s\" exist - but its not a directory\n", out);
			exit(1);
		}
	}
}

int main(int argc, char *argv[])
{
	HITS_DB _bblock;
	HITS_DB *bblock = &_bblock;
	char *bfile;
	char *broot;
	Align_Spec *settings;
	int isdam;

	int KMER_LEN;
	int BIN_SHIFT;
	int HIT_MIN;
	double AVE_ERROR;
	int SPACING;
	int NTHREADS;
	char *OUT_DIR;

	{
		KMER_LEN = 12;
		HIT_MIN = 35;
		BIN_SHIFT = 4;
		AVE_ERROR = .70;
		SPACING = 100;
		MINOVER = 500;    //   Globally visible to filter.c
		NTHREADS = 4;
		OUT_DIR = DEF_ARG_O;

		// parse arguments
		int c;
		opterr = 0;

		while ((c = getopt(argc, argv, "vk:w:h:e:l:s:o:j:")) != -1)
		{
			switch (c)
			{
			case 'v':
				VERBOSE = 1;
				break;
			case 'k':
				KMER_LEN = atoi(optarg);
				break;
			case 'w':
				BIN_SHIFT = atoi(optarg);
				break;
			case 'h':
				HIT_MIN = atoi(optarg);
				break;
			case 'e':
				AVE_ERROR = atof(optarg);
				break;
			case 'l':
				MINOVER = atoi(optarg);
				break;
			case 's':
				SPACING = atoi(optarg);
				break;
			case 'j':
				NTHREADS = atoi(optarg);
				break;
			case 'o':
				OUT_DIR = optarg;
				break;
			default:
				fprintf(stderr, "Unsupported option: %s\n", argv[optind - 1]);
				usage();
				exit(1);
			}
		}

		if (KMER_LEN < 0)
		{
			fprintf(stderr, "invalid kmer length of %d\n", KMER_LEN);
			exit(1);
		}
		if (HIT_MIN < 0)
		{
			fprintf(stderr, "invalid hit threshold of %d\n", HIT_MIN);
			exit(1);
		}
		if (AVE_ERROR < .5 || AVE_ERROR >= 1.)
		{
			fprintf(stderr, "Average correlation must be in [.5,1.) (%g)\n", AVE_ERROR);
			exit(1);
		}
		if (MINOVER < 0)
		{
			fprintf(stderr, "invalid minimum alignmnet length of (%d)\n", MINOVER);
			exit(1);
		}
		if (SPACING < 0)
		{
			fprintf(stderr, "invalid trace spacing of (%d)\n", SPACING);
			exit(1);
		}
		if (optind + 1 > argc)
		{
			fprintf(stderr, "[ERROR] - at least one subject block is required\n\n");
			usage();
			exit(1);
		}
	}

	MINOVER *= 2;
	if (Set_Filter_Params(KMER_LEN, BIN_SHIFT, HIT_MIN, NTHREADS))
	{
		fprintf(stderr, "Illegal combination of filter parameters\n");
		exit(1);
	}

	/* Create subdirectory */
	createSubdir(OUT_DIR);

	/* Compare each block against itself */

	{
		int i;

		broot = NULL;
		for (i = optind; i < argc; i++)
		{
			bfile = argv[i];
			isdam = read_DB(bblock, bfile, KMER_LEN);
			if (isdam)
				broot = Root(bfile, ".dam");
			else
				broot = Root(bfile, ".db");


			settings = New_Align_Spec(AVE_ERROR, SPACING, bblock->freq, NTHREADS, 1, 0, 0, 0);

			Match_Self(broot, bblock, settings);
			Write_Overlap_Buffer(settings, OUT_DIR, OUT_DIR, broot, broot, bblock->ufirst + bblock->nreads - 1);
			Reset_Overlap_Buffer(settings);

			Free_Align_Spec(settings);
			Close_DB(bblock);
		}
	}

	exit(0);
}
