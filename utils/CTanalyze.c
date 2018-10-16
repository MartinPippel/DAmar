/*
 * analyzeContigs.c
 *
 *  assign contigs into following sets: haploid, alternative, repeat
 *
 *  use:
 *  		self alignment information from contigs
 *  		union of fixed reads that map to the contigs (derived from fixedReadOverlaps)
 *
 *
 *  Created on: 3 Aug 2017
 *      Author: pippel
 */

#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#include "lib/colors.h"
#include "lib/oflags.h"
#include "lib/pass.h"
#include "lib/tracks.h"
#include "lib/utils.h"

#include "lib.ext/types.h"
#include "lib.ext/bitarr.h"
#include "utils/CTanalyze.h"

#undef DEBUG_STEP1A
#undef DEBUG_STEP1B
#define DEBUG_STEP1C

#define DEBUG_STEP2A

#define DEBUG
#define DEBUG_READ_INTERSECTION
#define DEBUG_CONTIG_COVERAGE
#define DEBUG_CONTIGOVERLAP

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#define DEF_REPEATS_TRACK "repeats"

#define MIN_OVL_LOOKAHEAD 5000
#define STRIDE_OVL_LOOKAHEAD 5000
#define MAX_OVL_LOOKAHEAD 30000    // todo actually should be adapted according to underlying repeat length

#define	DEF_SPUR_LEN 50000
#define DEF_TIP_LEN 100000

/// defines for contig vs contig alignments
#define ConVsConMinAlign 5000
#define ConVsConAnchorBases 2000
#define ConVsConContainmentThreshold 0.7 // contig is contained if 70% of the bases are covered by another larger contig
#define ConVsConShortContigLen 100000
#define ConVsConShortContigContainmentThreshold 0.5

ContigFlag2Label cflag2label[] =
{
{ CONTIG_DISCARD, "discard", 'd' },
{
CONTIG_IS_CONTAINED, "isContained", 'c' },
{ CONTIG_HAS_CONTAINED, "hasContained", 'C' },
{ CONTIG_UNCLASSIFIED, "unclassified", 'U' },
{ CONTIG_NORELATION, "norelation", 'N' },
{ CONTIG_UNIQUE, "unique", 'u' },
{ CONTIG_AL50PCTCOVERED, "atLeast50PctCoveredByOtherContigs", 'o' },

{ 0, NULL, '\0' } };

void cflags2str(char* pc, int flags)
{
	int i;
	int last = -1;
	for (i = 0; cflag2label[i].mask; i++)
	{
		if (flags & cflag2label[i].mask)
		{
			pc[i] = cflag2label[i].indicator;
			last = i;
		}
		else
		{
			pc[i] = ' ';
		}
	}

	pc[last + 1] = '\0';
}

static void printContigClassification(FILE* out, int flags, char sep)
{
	int i;
	for (i = 0; cflag2label[i].mask; i++)
	{
		if (flags & cflag2label[i].mask)
		{
			fprintf(out, "%s%c", cflag2label[i].label, sep);
		}
	}
}

//compare contigs by tmp flag, smallest first
static int cmp_contig_bytmp(const void* a, const void* b)
{
	Contig* c1 = *(Contig**) a;
	Contig* c2 = *(Contig**) b;

	return (c1->tmp - c2->tmp);
}

//compare contigs by length, longest first
static int cmp_contig_length(const void* a, const void* b)
{
	Contig* c1 = *(Contig**) a;
	Contig* c2 = *(Contig**) b;

	int cmp = c2->len - c1->len;

	if (!cmp)
	{
		cmp = (c2->avgCov - c1->avgCov);
	}

	return cmp;
}

void initAnalyzeContext(AnalyzeContext *actx)
{
	int i, j, k, l;
	int numOfContigs = DB_NREADS(actx->corContigDB);

	if (actx->VERBOSE)
		printf("numOfContigs: %d\n", numOfContigs);

	printf("allocate contig buffer");
	Contig *contigs = (Contig*) malloc(sizeof(Contig) * numOfContigs);
	actx->numContigs = numOfContigs;
	bzero(contigs, sizeof(Contig) * numOfContigs);

	// patched read IDs and Beg,End position from contig
	track_anno* patchedRead_anno = actx->corContigPatchReadsAndPos_track->anno;
	track_data* patchedRead_data = actx->corContigPatchReadsAndPos_track->data;

	// corrected read IDs and Beg,End position from corrected contig
	track_anno* correctedRead_anno = actx->corContigCorrectReadsAndPos_track->anno;
	track_data* correctedRead_data = actx->corContigCorrectReadsAndPos_track->data;

	printf(" - done\n");

	printf("allocate VreadMask");
	fflush(stdout);
	actx->maxVReadMask = 100;
	actx->vreadMask = (int**) malloc(sizeof(int*) * actx->maxVReadMask);
	assert(actx->vreadMask != NULL);
	for (i = 0; i < actx->maxVReadMask; i++)
	{
		actx->vreadMask[i] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
		bzero(actx->vreadMask[i], sizeof(int) * DB_NREADS(actx->patchedReadDB));
	}
	printf(" - done\n");

	actx->readSeq = New_Read_Buffer(actx->corContigDB);
	bzero(actx->readSeq, DB_READ_MAXLEN(actx->corContigDB));

	printf(" init contigs - START\n");
	for (i = 0; i < numOfContigs; i++)
	{
		// init contigs with reads
		Contig *contig = contigs + i;

		contig->len = DB_READ_LEN(actx->corContigDB, i);
		contig->idx = i;

		track_anno patched_rb, patched_re;
		track_anno corrected_rb, corrected_re;

		patched_rb = patchedRead_anno[i] / sizeof(track_data);
		patched_re = patchedRead_anno[i + 1] / sizeof(track_data);
		assert(patched_rb < patched_re);

		corrected_rb = correctedRead_anno[i] / sizeof(track_data);
		corrected_re = correctedRead_anno[i + 1] / sizeof(track_data);
		assert(corrected_rb < corrected_re);

		// get number of contig reads
		contig->numcReads = (int) ((patched_re - patched_rb) / 3);
		assert(contig->numcReads == (int ) ((corrected_re - corrected_rb) / 3));

		if (contig->numcReads == 0)
		{
			fprintf(stderr, "!strange! - no reads found for contig: %d. DISCARDED\n", i);
			contig->flag |= CONTIG_DISCARD;
			continue;
		}

		contig->numGClassific = 0;
		contig->maxGClassific = 1;
		contig->gClassific = (ContigGraphClassification*) malloc(sizeof(ContigGraphClassification) * contig->maxGClassific);

		contig->curContigChains = 0;
		contig->maxContigChains = 0;
		contig->contigChains = NULL;

		// init contig read buffer for overlapping reads
		contig->cReads = (ContigRead*) malloc(sizeof(ContigRead) * contig->numcReads);

		// init position buffer
		if (contig->cReads == NULL)
		{
			fprintf(stderr, "[ERROR] - ananlyzeContigs: Unable to allocate read buffers for Contig %d!\n", i);
			exit(1);
		}
		bzero(contig->cReads, sizeof(ContigRead) * contig->numcReads);

		int cov = MAX((int )(1.5 * actx->exp_cov), 30);
		for (j = 0, k = patched_rb; j < contig->numcReads; j++, k += 3)
		{
			ContigRead *cread = contig->cReads + j;

			cread->type |= CONTIG_READ_VALID;
			cread->patchedID = patchedRead_data[k];
			cread->correctedID = correctedRead_data[k];
			if (j)
			{
				cread->patchedContigPosBeg = contig->cReads[j - 1].patchedContigPosEnd;
				cread->patchedContigPosEnd = cread->patchedContigPosBeg + abs(patchedRead_data[k + 2] - patchedRead_data[k + 1]);

				cread->correctedContigPosBeg = contig->cReads[j - 1].correctedContigPosEnd;
				cread->correctedContigPosEnd = cread->correctedContigPosBeg + abs(correctedRead_data[k + 2] - correctedRead_data[k + 1]);
			}
			else
			{
				cread->patchedContigPosBeg = 0;
				cread->patchedContigPosEnd = abs(patchedRead_data[k + 2] - patchedRead_data[k + 1]);

				cread->correctedContigPosBeg = 0;
				cread->correctedContigPosEnd = abs(correctedRead_data[k + 2] - correctedRead_data[k + 1]);
			}
			cread->correctedReadPosBeg = correctedRead_data[k + 1];
			cread->correctedReadPosEnd = correctedRead_data[k + 2];

			cread->ovlReads = (OvlRead*) malloc(cov * sizeof(OvlRead));
			if (cread->ovlReads == NULL)
			{
				fprintf(stderr, "[ERROR] - ananlyzeContigs: Unable to allocate OvlRead buffer for Contig %d at patchedRead %d!\n", i, patchedRead_data[k]);
				exit(1);
			}

			cread->maxOvlReads = cov;

			// insert cRead into ovlReads, cBeg and cEnd will be later translated into absolute contig positions!
			cread->ovlReads->patchedID = -patchedRead_data[k];
			cread->ovlReads->beg = patchedRead_data[k + 1];
			cread->ovlReads->end = patchedRead_data[k + 2];
			cread->ovlReads->cBeg = cread->patchedContigPosBeg;
			cread->ovlReads->cEnd = cread->patchedContigPosEnd;

			cread->numOvlReads = 1;

			// set vmaskread
			l = 0;
			while (l < actx->maxVReadMask && actx->vreadMask[l][patchedRead_data[k]] != 0)
				l++;

			if (l == actx->maxVReadMask)
			{
				printf("TERROR read %d is part of %d contigs!!!!", patchedRead_data[k], l);
				actx->maxVReadMask++;

				actx->vreadMask = (int**) realloc(actx->vreadMask, sizeof(int*) * actx->maxVReadMask);
				actx->vreadMask[l] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
				bzero(actx->vreadMask[l], sizeof(int) * DB_NREADS(actx->patchedReadDB));
			}
			// store negative 1-based contig id

			actx->vreadMask[l][patchedRead_data[k]] = -(i + 1);
		}

		if (actx->VERBOSE)
		{
			printf("Contig_%d, reads (%d): [ unCorRead bpos epos len cBpos cEpos | corRead bpos epos len cBpos cEpos| corCumLen unCorCumLen]\n", i,
					contig->numcReads);
			int corCumLen = 0;
			int unCorCumLen = 0;
			for (j = 0; j < contig->numcReads; j++)
			{
				corCumLen += contig->cReads[j].patchedContigPosEnd - contig->cReads[j].patchedContigPosBeg;
				unCorCumLen += contig->cReads[j].correctedContigPosEnd - contig->cReads[j].correctedContigPosBeg;
				// uncorrected (patched) coordinates
				printf("%8d", contig->cReads[j].patchedID);
				printf(" %8d", contig->cReads[j].ovlReads[0].beg);
				printf(" %8d", contig->cReads[j].ovlReads[0].end);
				printf(" %5d", abs(contig->cReads[j].ovlReads[0].end - contig->cReads[j].ovlReads[0].beg));
				printf(" %8d", contig->cReads[j].patchedContigPosBeg);
				printf(" %8d", contig->cReads[j].patchedContigPosEnd);
				//  correx=cted coordinates
				printf(" |	 %8d", contig->cReads[j].correctedID);
				printf(" %8d", contig->cReads[j].correctedReadPosBeg);
				printf(" %8d", contig->cReads[j].correctedReadPosEnd);
				printf(" %5d", abs(contig->cReads[j].correctedReadPosEnd - contig->cReads[j].correctedReadPosBeg));
				printf(" %8d", contig->cReads[j].correctedContigPosBeg);
				printf(" %8d", contig->cReads[j].correctedContigPosEnd);

				printf(" | %8d", corCumLen);
				printf(" %8d", unCorCumLen);
				printf(" %5.3f%%\n", unCorCumLen * 100.0 / corCumLen);
			}
			printf(" --> len: %10d\n", contig->len);
		}
	}
	printf(" init contigs - DONE\n");

	parseDBFileNames(actx->corContigDBName, &actx->contigFileNameOffsets, &actx->ContigFileNames, NULL);

	actx->contigs = contigs;
	printf("sort contigs by length - START\n");
	// sort contig pointer according length (longest first)
	actx->contigs_sorted = malloc(sizeof(Contig*) * actx->numContigs);
	assert(actx->contigs_sorted != NULL);
	for (i = 0; i < actx->numContigs; i++)
	{
		actx->contigs_sorted[i] = actx->contigs + i;
	}
	qsort(actx->contigs_sorted, actx->numContigs, sizeof(Contig*), cmp_contig_length);
	printf("sort contigs by length - DONE\n");
}

char* getFastaFileNameFromDB(AnalyzeContext *actx, int contigID)
{
	if (contigID < 0 || contigID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getContigFileName: contigID %d out of bounds [%d, %d [\n", contigID, 0, DB_NREADS(actx->corContigDB));
		fflush(stderr);
		exit(1);
	}

	int i = 1;
	while (1)
	{
		if (contigID >= actx->contigFileNameOffsets[i - 1] && contigID < actx->contigFileNameOffsets[i])
		{
			break;
		}
		i++;

		assert(i <= 1000000);
	}

	return actx->ContigFileNames[i];
}

int getFirstDBIdxFromFastaFileName(AnalyzeContext *actx, int contigID)
{
	if (contigID < 0 || contigID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getNumberOfSequencesFromFastaFileName: contigID %d out of bounds [%d, %d [\n", contigID, 0, DB_NREADS(actx->corContigDB));
		fflush(stderr);
		exit(1);
	}

	int i = 1;
	while (1)
	{
		if (contigID >= actx->contigFileNameOffsets[i - 1] && contigID < actx->contigFileNameOffsets[i])
		{
			break;
		}
		i++;

		assert(i <= 1000000);
	}

	return actx->contigFileNameOffsets[i - 1];
}

int getNumberOfSequencesFromFastaFileName(AnalyzeContext *actx, int contigID)
{
	if (contigID < 0 || contigID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getNumberOfSequencesFromFastaFileName: contigID %d out of bounds [%d, %d [\n", contigID, 0, DB_NREADS(actx->corContigDB));
		fflush(stderr);
		exit(1);
	}

	int i = 1;
	while (1)
	{
		if (contigID >= actx->contigFileNameOffsets[i - 1] && contigID < actx->contigFileNameOffsets[i])
		{
			break;
		}
		i++;

		assert(i <= 1000000);
	}

	return actx->contigFileNameOffsets[i] - actx->contigFileNameOffsets[i - 1];
}

void parseDBFileNames(char *dbName, int **fileOffsets, char ***fileNames, char ***readNames)
{
	// parse old database file and store bax names as well as their read offsets
	FILE *f;
	int numFiles = 0;
	char name[MAX_NAME];

	f = fopen(dbName, "r");
	if (f == NULL) // try to add .db extension
	{
		sprintf(name, "%s.db", dbName);
		f = fopen(name, "r");
	}
	if (f == NULL)
	{
		fprintf(stderr, "[ERROR] parseDBFileNames: Cannot open file %s for reading\n", dbName);
		exit(1);
	}

	if (fscanf(f, DB_NFILE, &numFiles) != 1)
	{
		fprintf(stderr, "[ERROR] parseDBFileNames: Cannot parse \"%s\" from file %s.\n", DB_NFILE, dbName);
		fclose(f);
		exit(1);
	}

	int *offsets = (int *) malloc(sizeof(int) * (numFiles + 2));
	memset(offsets, 0, numFiles + 2);

	if (offsets == NULL)
	{
		fprintf(stderr, "[ERROR] parseDBFileNames: Cannot allocate buffer for read offsets.\n");
		fclose(f);
		exit(1);
	}

	char **fnames = (char **) malloc(sizeof(char*) * (numFiles + 1));
	char **rnames = (char **) malloc(sizeof(char*) * (numFiles + 1));
	if (fnames == NULL || rnames == NULL)
	{
		fprintf(stderr, "[ERROR] parseDBFileNames: Cannot allocate buffer for file names.\n");
		fclose(f);
		exit(1);
	}

	int i;
	for (i = 1; i <= numFiles; i++)
	{
		fnames[i] = (char*) malloc(sizeof(char) * MAX_NAME);
		rnames[i] = (char*) malloc(sizeof(char) * MAX_NAME);
		if (fnames[i] == NULL || rnames[i] == NULL)
		{
			fprintf(stderr, "[ERROR] parseDBFileNames: Cannot allocate %d st file name buffer of size %d.\n", i, MAX_NAME);
			fclose(f);
			exit(1);
		}

	}

	for (i = 1; i <= numFiles; i++)
	{
		if (fscanf(f, DB_FDATA, &offsets[i], fnames[i], rnames[i]) != 3)
		{
			fprintf(stderr, "[ERROR] parseDBFileNames: Cannot parse \"%s\" in line %d from file %s.\n", DB_FDATA, i, dbName);
			fclose(f);
			exit(1);
		}
	}
	// add number of db reads as last entry
	offsets[i] = offsets[i - 1];
	fclose(f);

	if (fileNames != NULL)
		*fileNames = &(*fnames);
	if (readNames != NULL)
		*readNames = &(*rnames);
	*fileOffsets = &(*offsets);
}

static void pre(PassContext* pctx, AnalyzeContext* actx)
{
	actx->twidth = pctx->twidth;
}

ContigRead* getContigRead(Contig* contig, int readID)
{
	int i;

	for (i = 0; i < contig->numcReads; i++)
	{
		if (contig->cReads[i].patchedID == readID)
		{
			return contig->cReads + i;
		}
	}
	return NULL;
}

int getFullReadPositionInContig(AnalyzeContext *actx, Contig* contig, int readID, int *cbeg, int *cend)
{
	int rbeg;
	int tbeg, tend;
	get_trim(actx->patchedReadDB, actx->patchedReadTrim_track, readID, &tbeg, &tend);

	int i, j;
	int readPos = -1;

	for (i = 0; i < contig->numcReads; i++)
	{
		if (contig->cReads[i].patchedID == readID)
		{
			readPos = i;
			break;
		}
	}

	assert(readPos >= 0 && readPos < contig->numcReads);

	ContigRead *cread = contig->cReads + readPos;
//	printf("reads %d %d %d %d %d cIdx: %d\n", reads->id, reads->cBeg, reads->cEnd, reads->beg, reads->end, readPos);

	assert(abs(cread->patchedID) == readID);

	if (cread->ovlReads->beg < cread->ovlReads->end) // forward orientation
	{
		*cend = cread->ovlReads->cEnd;
		assert(cread->ovlReads->end == tend);
		*cbeg = cread->ovlReads->cBeg;
		if (cread->ovlReads->beg != tbeg) // try to find overlaps for readID with previous cReads
		{
			i = 1;
			while (readPos - i >= 0)
			{
//				printf("try to find aread %d in overlaps with previous Cread %d at position %d\n", readID, neighborCread, readPos - i);
				cread = contig->cReads + (readPos - i);
//				printf("cur cread %d %d %d %d %d\n", reads->id, reads->cBeg, reads->cEnd, reads->beg, reads->end);
				int found = 0;
				for (j = 1; j < cread->numOvlReads; j++)
				{
//					printf("check %d c(%d,%d) r(%d, %d)\n", reads[j].id, reads[j].cBeg, reads[j].cEnd, reads[j].beg, reads[j].end);
					if (cread->ovlReads[j].patchedID == readID)
					{
//						printf("found\n");
						assert(*cbeg > cread->ovlReads[j].cBeg);
						*cbeg = cread->ovlReads[j].cBeg;
						found = 1;

						// i.e. previous cread is in forward orientation
						if (((cread->ovlReads->beg < cread->ovlReads->end)))
						{
							rbeg = cread->ovlReads[j].beg;
							if (cread->ovlReads[j].beg == tbeg)
								found = 2;
						}
						if (cread->ovlReads->beg > cread->ovlReads->end)
						{
							rbeg = cread->ovlReads[j].end;
							if (cread->ovlReads[j].end == tbeg)
								found = 2;
						}
						break; // found start of current readID in contig
					}
				}
				if (found == 2)
				{
//					printf("iteration %d: found proper start of aread. i.e. proper cut position\n", i);
					break;
				}
				if (found == 0)
				{
					if (abs(rbeg) - tbeg < 1000) // i.e. there is no overlap, due to the minimum overlap length of daligner
					{
						*cbeg = (*cbeg) - (abs(rbeg) - tbeg);
						assert((*cbeg) >= 0 && (*cbeg) < (*cend) && (*cend) <= contig->len);
					}
					else
					{
						fprintf(stderr, "TERROR: iteration %d: UNABLE to find proper start of aread. i.e. unclear cut position. Remaining bases: %d\n", i,
								abs(rbeg) - tbeg);
					}
					break;
				}
				i++;
			}
		}
	}
	else // reverse complement orientation
	{
		*cend = cread->ovlReads->cEnd;
		assert(cread->ovlReads->end == tbeg);
		*cbeg = cread->ovlReads->cBeg;
		if (cread->ovlReads->beg != tend) // try to find overlaps for readID with previous cReads
		{
			i = 1;
			while (readPos - i >= 0)
			{
//				int neighborCread = contig->contigReadIds[readPos - i];
//				printf("try to find aread %d in overlaps with previous Cread %d at position %d\n", readID, neighborCread, readPos - i);
				cread = contig->cReads + (readPos - i);
//				printf("cur cread %d %d %d %d %d\n", reads->id, reads->cBeg, reads->cEnd, reads->beg, reads->end);
				int found = 0;
				for (j = 1; j < cread->numOvlReads; j++)
				{
//					printf("check %d c(%d,%d) r(%d, %d)\n", reads[j].id, reads[j].cBeg, reads[j].cEnd, reads[j].beg, reads[j].end);
					if (cread->ovlReads[j].patchedID == readID)
					{
//						printf("found\n");
						assert(*cbeg > cread->ovlReads[j].cBeg);
						*cbeg = cread->ovlReads[j].cBeg;

						found = 1;

						// i.e. previous cread is in forward orientation
						if (cread->ovlReads->beg < cread->ovlReads->end)
						{
							rbeg = cread->ovlReads[j].beg;
							if (cread->ovlReads[j].beg == tend)
								found = 2;
						}

						if (cread->ovlReads->beg > cread->ovlReads->end)
						{
							rbeg = cread->ovlReads[j].end;
							if (cread->ovlReads[j].end == tend)
								found = 2;
						}

						break; // found start of current readID in contig
					}
				}
				if (found == 2)
				{
//					printf("iteration %d: found proper start of aread. i.e. proper cut position\n", i);
					break;
				}
				if (found == 0)
				{
					if (tend - abs(rbeg) < 1000) // i.e. there is no overlap, due to the minimum overlap length of daligner
					{
						*cbeg = (*cbeg) - (tend - abs(rbeg));
						assert((*cbeg) >= 0 && (*cbeg) < (*cend) && (*cend) <= contig->len);
					}
					else
					{
						fprintf(stderr, "TERROR: iteration %d: UNABLE to find proper start of aread. i.e. unclear cut position. Remaining bases: %d\n", i,
								tend - abs(rbeg));
					}
					break;
				}
				i++;
			}
		}
	}

	return 0;
}

//void printSplitEvent(Contig *contig, SplitEvent *split)
//{
//	printf("CONTIG %d len %d new SplitEvent [%s] duplRead %d absPos (%d, %d, %d) leftRead (%d, %d,%d) rightRead (%d, %d,%d)\n", contig->idx, contig->len,
//			getSplitType(split), split->duplicRead, split->type, split->contigPos1, split->contigPos2, split->leftProperRead, split->leftProperReadPos1,
//			split->leftProperReadPos2, split->rightProperRead, split->rightProperReadPos1, split->rightProperReadPos2);
//}

void analyzeJunction(AnalyzeContext* actx, int read, int numContigs)
{
	printf("analyzeJunction(AnalyzeContext* actx, int read %d , int numContigs %d)\n", read, numContigs);
	if (numContigs < 2)
		return;

	int i, j;

	// first do the easy cases:
	if (numContigs == 2)
	{
		// 1. check classification of both contigs
		// 1a) haploid vs contained todo check if we hav

		Contig *contig_i, *contig_j;
		contig_i = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		contig_j = actx->contigs + abs(actx->vreadMask[1][read]) - 1;

		// check if we have a normal primary and secondary contig connection
		if (contig_i->len < contig_j->len)
		{
			contig_i = actx->contigs + abs(actx->vreadMask[1][read]) - 1;
			contig_j = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		}

		// check putative heterozygous bubbles and spurs and their coverage (switch contents if necessary)
		// i.e.
		//						A1     								H1
		//    i.e.			 _______  						     ________
		//   		H1	____/	    \________ H1    OR    H1 ____/        \_____ H1
		//   				\_______/							   |____	/
		//                      H1        							A1
		// check haploid coverage
		int path_j_beg, path_j_end;
		getContigsEndIDs(actx, contig_j->idx, &path_j_beg, &path_j_end);

		int containedByContigOverlaps = 0;
//		for (i = 0; i < contig_j->numOvlGrps; i++)
//		{
//			OverlapGroup *olg = contig_j->ovlGrps[i];
//			if ((olg->bread == contig_i->idx) && (olg->flag & OVLGRP_AREAD_IS_CONTAINED))
//			{
//				containedByContigOverlaps = 1;
//				break;
//			}
//		}

		int containedByFixedReadCoverage = 0;
		for (i = 0; i < contig_j->numGClassific; i++)
		{
			ContigGraphClassification *gclass = contig_j->gClassific + i;
			if ((gclass->correspID == contig_i->idx) && (gclass->flag & CONTIG_IS_CONTAINED))
			{
				containedByFixedReadCoverage = 1;
				break;
			}
		}

		// check coverage along junction
		checkJunctionCoverage(actx, contig_i, read);
		checkJunctionCoverage(actx, contig_j, read);

	}
	return;

	int cutOut = 0;
	//////////////////////////////////////////////////////////////////////// MULTIJOIN - EXACTLY 3 CONTIGS ////////////////////////////////////////////////////////////////////////
	if (numContigs == 3)
	{
		//ignore split if, two consecutive alt contigs collapse at same read in hap contig
		//						A1      H1 									H1         H2
		//    i.e.			 _______  ______								 ________  _________
		//   		H1	____/	    \/      \________ H1    OR    H1 ____/        \/         \_____ H2
		//   				\_______/\______/							   |____	/ \________/
		//                      H1       A2 									H2        A1

		Contig *contig_i, *contig_j, *contig_k;
		contig_i = contig_j = contig_k = NULL;

		contig_i = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		contig_j = actx->contigs + abs(actx->vreadMask[1][read]) - 1;
		contig_k = actx->contigs + abs(actx->vreadMask[2][read]) - 1;

		// all 3 contigs are haploid | or none of them is haploid
		if (((contig_i->class & CONTIG_CLASS_HAPLOID) && (contig_j->class & CONTIG_CLASS_HAPLOID) && (contig_k->class & CONTIG_CLASS_HAPLOID))
				|| (!(contig_i->class & CONTIG_CLASS_HAPLOID) && !(contig_j->class & CONTIG_CLASS_HAPLOID) && !(contig_k->class & CONTIG_CLASS_HAPLOID)))
		{
			cutOut = 1;
		}
		else
		{
			Contig *contig_hap1 = NULL;
			Contig *contig_hap2 = NULL;
			Contig *contig_alt1 = NULL;
			Contig *contig_alt2 = NULL;

			if (contig_i->class & CONTIG_CLASS_HAPLOID)
			{
				contig_hap1 = contig_i;

				if (contig_j->class & CONTIG_CLASS_HAPLOID)
				{
					contig_hap2 = contig_j;
					contig_alt1 = contig_k;
				}
				else if (contig_k->class & CONTIG_CLASS_HAPLOID)
				{
					contig_hap2 = contig_k;
					contig_alt1 = contig_j;
				}
				else
				{
					contig_alt1 = contig_j;
					contig_alt2 = contig_k;
				}
			}
			else if (contig_j->class & CONTIG_CLASS_HAPLOID)
			{
				contig_hap1 = contig_j;
				contig_alt1 = contig_i;
				if (contig_k->class & CONTIG_CLASS_HAPLOID)
				{
					contig_hap2 = contig_k;
				}
				else
				{
					contig_alt2 = contig_k;
				}
			}
			else if (contig_k->class & CONTIG_CLASS_HAPLOID)
			{
				contig_hap1 = contig_k;
				contig_alt1 = contig_i;
				contig_alt2 = contig_j;
			}
			else
			{
				fprintf(stderr, "SHIT happens\n");
				exit(1);
			}

			// check for first pattern
			if (contig_alt1 != NULL && contig_alt2 != NULL && contig_hap1 != NULL)
			{
				int a1End1, a1End2, a2End1, a2End2;
				getContigsEndIDs(actx, contig_alt1->idx, &a1End1, &a1End2);
				getContigsEndIDs(actx, contig_alt2->idx, &a2End1, &a2End2);

				if (((((contig_alt1->class & CONTIG_CLASS_ALT_BUBBLE) || (contig_alt1->class & CONTIG_CLASS_ALT_SPUR))
						&& (a1End1 == contig_hap1->idx || a1End2 == contig_hap1->idx))
						|| ((contig_alt1->class & CONTIG_CLASS_ALT_HUGEDIFF) && a1End1 == contig_hap1->idx && a1End2 == contig_hap1->idx))
						&& ((((contig_alt2->class & CONTIG_CLASS_ALT_BUBBLE) || (contig_alt2->class & CONTIG_CLASS_ALT_SPUR))
								&& (a1End1 == contig_hap1->idx || a1End2 == contig_hap1->idx))
								|| ((contig_alt2->class & CONTIG_CLASS_ALT_HUGEDIFF) && a1End1 == contig_hap1->idx && a1End2 == contig_hap1->idx))
						//&& checkJunctionCoverage(actx, read, numContigs) // returns true if everything is ok, false otherwise
						)
				{
					// assign to both alt contigs a heterozygous split
					if (contig_alt1->cReads[0].patchedID == read && contig_alt2->cReads[contig_alt2->numcReads - 1].patchedID == read)
					{
//						createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
//						createSplitEvent(actx, contig_alt2, read, SPLIT_HETEROZYGRIGHT);
//						createSplitEvent(actx, contig_hap1, read, (SPLIT_HETEROZYGLEFT | SPLIT_HETEROZYGRIGHT));
					}
					else if (contig_alt2->cReads[0].patchedID == read && contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
					{
//						createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
//						createSplitEvent(actx, contig_alt2, read, SPLIT_HETEROZYGLEFT);
//						createSplitEvent(actx, contig_hap1, read, (SPLIT_HETEROZYGLEFT | SPLIT_HETEROZYGRIGHT));
					}
					else // always fall back to split behavior, i.e. cut out read from all contigs
					{
						cutOut = 1;
					}

				}
				else
				{
					cutOut = 1;
				}
			}
			// check for second pattern
			else if (contig_alt1 != NULL && contig_hap1 != NULL && contig_hap2 != NULL)
			{
				int end1, end2;
				getContigsEndIDs(actx, contig_alt1->idx, &end1, &end2);

				// check if haploid contigs have dead end
				int deadEnd = 0;
				int contigCutPos1, contigCutPos2;
				contigCutPos1 = contigCutPos2 = -1;
				if (contig_hap1->cReads[0].patchedID == read || contig_hap1->cReads[contig_hap1->numcReads - 1].patchedID == read)
				{
					getFullReadPositionInContig(actx, contig_hap1, read, &contigCutPos1, &contigCutPos2);

					if (contigCutPos1 <= actx->TIP_LEN)
					{
						deadEnd = 1;
					}
					else if (contig_hap1->len - contigCutPos2 <= actx->TIP_LEN)
					{
						deadEnd = 2;
					}
				}
				else if (contig_hap2->cReads[0].patchedID == read || contig_hap2->cReads[contig_hap2->numcReads - 1].patchedID == read)
				{
					getFullReadPositionInContig(actx, contig_hap2, read, &contigCutPos1, &contigCutPos2);

					if (contigCutPos1 <= actx->TIP_LEN)
					{
						deadEnd = 1;
					}
					else if (contig_hap1->len - contigCutPos2 <= actx->TIP_LEN)
					{
						deadEnd = 2;
					}
				}

				if (deadEnd
						&& ((((contig_alt1->class & CONTIG_CLASS_ALT_BUBBLE) || (contig_alt1->class & CONTIG_CLASS_ALT_SPUR))
								&& (end1 == contig_hap1->idx || end2 == contig_hap1->idx || end1 == contig_hap2->idx || end2 == contig_hap2->idx))
								|| ((contig_alt1->class & CONTIG_CLASS_ALT_HUGEDIFF)
										&& (((end1 == contig_hap1->idx && end2 == contig_hap1->idx) || (end1 == contig_hap2->idx && end2 == contig_hap2->idx)))))
						//&& checkJunctionCoverage(actx, read, numContigs) // returns true if everything is ok, false otherwise
						)
				{
					if (end1 == contig_hap1->idx && end2 == contig_hap1->idx) // ALT contig is a bubble of contig_hap1
					{
						if (contig_alt1->cReads[0].patchedID == read)
						{
//							createSplitEvent(actx, contig_hap1, read, SPLIT_HETEROZYGLEFT);
//							createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
						}
						else if (contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
						{
//							createsplitevent(actx, contig_hap1, read, split_heterozygright);
//							createsplitevent(actx, contig_alt1, read, split_heterozygright);
						}
						else
						{
							cutOut = 1;
						}
					}
					else if (end1 == contig_hap2->idx && end2 == contig_hap2->idx) // ALT contig is a bubble of contig_hap2
					{
						if (contig_alt1->cReads[0].patchedID == read)
						{
//							createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGLEFT);
//							createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
						}
						else if (contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
						{
//							createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGRIGHT);
//							createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
						}
						else
						{
							cutOut = 1;
						}
					}
					else if (contig_alt1->class == CONTIG_CLASS_ALT_SPUR)
					{
						if (contig_alt1->cReads[0].patchedID == read)
						{
							if (end1 == contig_hap1->idx)
							{
//								createSplitEvent(actx, contig_hap1, read, SPLIT_HETEROZYGLEFT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
							}
							else if (end1 == contig_hap2->idx)
							{
//								createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGLEFT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGLEFT);
							}
							else
							{
								cutOut = 1;
							}
						}
						else if (contig_alt1->cReads[contig_alt1->numcReads - 1].patchedID == read)
						{
							if (end2 == contig_hap1->idx)
							{
//								createSplitEvent(actx, contig_hap1, read, SPLIT_HETEROZYGRIGHT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
							}
							else if (end2 == contig_hap2->idx)
							{
//								createSplitEvent(actx, contig_hap2, read, SPLIT_HETEROZYGRIGHT);
//								createSplitEvent(actx, contig_alt1, read, SPLIT_HETEROZYGRIGHT);
							}
							else
							{
								cutOut = 1;
							}
						}
						else
							cutOut = 1;
					}
					else
						cutOut = 1;

					if (!cutOut)
					{
						if (deadEnd == 1)
						{
//							createSplitEvent(actx, contig_hap1, read, SPLIT_DEADENDLEFT);
//							createSplitEvent(actx, contig_hap2, read, SPLIT_DEADENDLEFT);
						}
						else
						{
//							createSplitEvent(actx, contig_hap1, read, SPLIT_DEADENDRIGHT);
//							createSplitEvent(actx, contig_hap2, read, SPLIT_DEADENDRIGHT);
						}
					}

				}
				else
				{
					cutOut = 1;
				}

			}
			// split otherwise
			else
			{
				cutOut = 1;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////// MULTIJOIN - MORE THAN 3 CONTIGS ////////////////////////////////////////////////////////////////////////
	if (numContigs > 3 || (numContigs == 3 && cutOut == 1)) // thats really strange and should never happen
	{
		fprintf(stderr, "[WARNING] - Read %d is part of multiple contigs: ", read);

		for (i = 0; i < numContigs; i++)
		{
			printf(" %d", actx->vreadMask[i][read]);
		}
		printf("\n");
		assert(i == numContigs);

		j = 0;
		while (j < numContigs)
		{
			Contig *contig = actx->contigs + abs(actx->vreadMask[j][read]) - 1;
//			createSplitEvent(actx, contig, read, SPLIT_MULTIREAD);

		}
	}

	//////////////////////////////////////////////////////////////////////// READ IS PART OF EXACTLY 2 CONTIGS ////////////////////////////////////////////////////////////////////////
	///
	/// can be of type
	///
	/// 1. START_END      ------contig_i------X-----contig_j----, where X is read that is part of both contigs ---> always split those things
	///
	/// 2. BUBBLE/SPUR
	///
	///  			A)	ALT BUBBLE						b) ALT SPUR							C) DEADEND (SPUR or HAP and HAP)                  D) FALSE JOIN
	///
	///                   H1									H1										H1/H2
	///                _________							 __________								 ___________
	///      H1 ______/         \______ H1       H1 _____/          \________ H1			H1  ____/    		\______ H1/H2		     CONTIG_I_________X____________CONTIG_I
	///               \_________/						\_____|									\______|								    			      |
	///                   A1								   A1									   A1/H1											      |________ CONTIG_J
	else
	{

		// check START_END
		Contig *contig_j = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
		Contig *contig_k = actx->contigs + abs(actx->vreadMask[1][read]) - 1;

		printf("else: %d vs %d\n", contig_j->idx, contig_k->idx);
		/////
		//// check the weird stuff first
		///
		// contig_k is circular
		if (contig_k->cReads[0].patchedID == read && contig_k->cReads[contig_k->numcReads - 1].patchedID == read)
		{
			fprintf(stderr, "-----------------------> circle in contig_k %d\n", contig_k->idx);
			exit(1);
		}
		// contig_j is circular
		else if (contig_j->cReads[0].patchedID == read && contig_j->cReads[contig_j->numcReads - 1].patchedID == read)
		{
			fprintf(stderr, "-----------------------> circle in contig_j %d\n", contig_j->idx);
			exit(1);
		}
		// contig_j stops at read X and contig_k starts at read X
		else if ((contig_k->cReads[0].patchedID == read || contig_k->cReads[contig_k->numcReads - 1].patchedID == read)
				&& (contig_j->cReads[contig_j->numcReads - 1].patchedID == read || contig_j->cReads[0].patchedID == read))
		{
			printf("// contig_j stops at read X and contig_k starst at read X\n");
//			createSplitEvent(actx, contig_j, read, SPLIT_STARTEND);
//			createSplitEvent(actx, contig_k, read, SPLIT_STARTEND);
		}
		else
		{
			//// first, get position of read X in both contigs
			int contigJCutPos1, contigJCutPos2;
			int contigKCutPos1, contigKCutPos2;

			contigJCutPos1 = contigJCutPos2 = contigKCutPos1 = contigKCutPos2 = -1;

			Contig *contig_j = actx->contigs + abs(actx->vreadMask[0][read]) - 1;
			Contig *contig_k = actx->contigs + abs(actx->vreadMask[1][read]) - 1;

			getFullReadPositionInContig(actx, contig_j, read, &contigJCutPos1, &contigJCutPos2);
			getFullReadPositionInContig(actx, contig_k, read, &contigKCutPos1, &contigKCutPos2);

			assert(contigJCutPos1 != -1 && contigJCutPos2 != -1);
			assert(contigKCutPos1 != -1 && contigKCutPos2 != -1);

			printf("contig %d read %d at [%d, %d] contig %d read %d at [%d, %d]\n", contig_j->idx, read, contigJCutPos1, contigJCutPos2, contig_k->idx, read,
					contigKCutPos1, contigKCutPos2);
			//int type = 0;

			//int junctionCov = checkJunctionCoverage(actx, read, numContigs); // returns true if everything is ok, false otherwise

			// both contigs are haploid
			//if (contig_j->fClass == HAPLOID && contig_k->fClass == HAPLOID)
			if (1)
			{
				// check if there is a valuable dead end, otherwise its a false join

				// todo remove this after testing
				//	assert(contig_j->len > 200000 || contig_k->idx > 200000);

				if (contigJCutPos1 == 0)
				{
					if (contigKCutPos1 > actx->TIP_LEN && contig_k->len - contigKCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("1 DEAD END ----> check number of read that follow both paths\n");
					}
				}
				else if (contigJCutPos2 == contig_j->len)
				{
					if (contigKCutPos1 > actx->TIP_LEN && contig_k->len - contigKCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("2 DEAD END ----> check number of read that follow both paths\n");
					}
				}

				else if (contigKCutPos1 == 0)
				{
					if (contigJCutPos1 > actx->TIP_LEN && contig_j->len - contigJCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("3 DEAD END ----> check number of read that follow both paths\n");
					}
				}
				else if (contigKCutPos2 == contig_k->len)
				{
					if (contigJCutPos1 > actx->TIP_LEN && contig_j->len - contigJCutPos2 > actx->TIP_LEN)
					{
//						type = SPLIT_FALSEJOINLEFT;
					}
					else
					{
						// check number of reads that follow both paths
						printf("4 DEAD END ----> check number of read that follow both paths\n");

					}
				}

			}
			// only one contig is haploid
			else if ((contig_j->class & CONTIG_CLASS_HAPLOID) || (contig_k->class & CONTIG_CLASS_HAPLOID))
			{
				// TODO
			}
			// both contigs are NOT haploid
			else
			{
				// TODO
			}

			/*

			 /////
			 //// check the "normal" problems. Distinguish between
			 ///  1) normal heterozygous bubbles/spurs (i.e. clear ALT contigs)
			 ///  2) putative false joins
			 ///  3) as well as overlapping contigs in heterozygous regions // i.e. dead end
			 ///




			 int stype = 0;

			 if (contigJCutPos1 == 0)
			 {
			 if (contigKCutPos1 > actx->TIP_LEN && contig_k->len - contigKCutPos2 > actx->TIP_LEN) // well inside contigK
			 {
			 if (contig_j->len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGRIGHT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINRIGHT;
			 }
			 }
			 else // at tips of contigK --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDRIGHT;
			 }
			 }
			 else if (contigKCutPos1 == 0)
			 {
			 if (contigJCutPos1 > actx->TIP_LEN && contig_j->len - contigJCutPos2 > actx->TIP_LEN) // well inside contigJ
			 {
			 if (contig_k->len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGRIGHT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINRIGHT;
			 }
			 }
			 else // at tips of contigJ --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDRIGHT;
			 }
			 }

			 else if (contigJCutPos2 == contig_j->len)
			 {
			 if (contigKCutPos1 > actx->TIP_LEN && contig_k->len - contigKCutPos2 > actx->TIP_LEN) // well inside contigK
			 {
			 if (contig_j->len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGLEFT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINLEFT;
			 }
			 }
			 else // at tips of contigK --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDLEFT;
			 }
			 }
			 else if (contigKCutPos2 == contig_k->len)
			 {
			 if (contigJCutPos1 > actx->TIP_LEN && contig_j->len - contigJCutPos2 > actx->TIP_LEN) // well inside contigJ
			 {
			 if (contig_k->len <= actx->SPUR_LEN)
			 {
			 stype |= SPLIT_HETEROZYGLEFT;
			 }
			 else
			 {
			 stype |= SPLIT_FALSEJOINLEFT;
			 }
			 }
			 else // at tips of contigK --> can be a fork (2 dead ends == heterozygous + repeat, or 1 dead end == heterozygous )
			 {
			 stype |= SPLIT_DEADENDLEFT;
			 }
			 }
			 else
			 {
			 printf("SHIT no classification\n");
			 exit(1);
			 }

			 printf("############## remaining duplicate Read: %d: contig_j %d len %d pos %d, %d TYPE: %d and contig_k %d len %d pos %d, %d TYPE %d\n",
			 aread, contig_j->idx, contig_j->len, contigJCutPos1, contigJCutPos2, stype, contig_k->idx, contig_k->len, contigKCutPos1,
			 contigKCutPos2, stype);

			 // add split position to contig_j
			 if (contig_j->maxSplits == contig_j->numSplits)
			 {
			 contig_j->maxSplits = (contig_j->maxSplits * 1.2) + 2;
			 contig_j->split = (SplitEvent*) realloc(contig_j->split, sizeof(SplitEvent) * contig_j->maxSplits);
			 assert(contig_j->split != NULL);
			 }
			 assert(contigJCutPos1 >= 0 && contigJCutPos1 < contigJCutPos2 && contigJCutPos2 <= contig_j->len);
			 {
			 // find contig reads and their positions for corresponding cutPositions
			 int leftProperRead, rightProperRead;
			 leftProperRead = rightProperRead = -1;
			 int leftProperReadPos1, leftProperReadPos2, rightProperReadPos1, rightProperReadPos2;
			 leftProperReadPos1 = leftProperReadPos2 = rightProperReadPos1 = rightProperReadPos2 = 0;
			 for (k = 0; k < contig_j->numContigReadIds; k++)
			 {
			 if (leftProperRead < 0 && contig_j->reads[k]->cEnd >= contigJCutPos1 && contigJCutPos1 > 0)
			 {
			 leftProperRead = abs(contig_j->reads[k]->id);
			 leftProperReadPos1 = contig_j->reads[k]->beg;
			 if (contig_j->reads[k]->beg > contig_j->reads[k]->end) // reverse complement orientation
			 {
			 leftProperReadPos2 = contig_j->reads[k]->end + (contig_j->reads[k]->cEnd - contigJCutPos1);
			 printf("RV %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 > leftProperReadPos2);
			 }
			 else // forward orientation
			 {
			 leftProperReadPos2 = contig_j->reads[k]->end - (contig_j->reads[k]->cEnd - contigJCutPos1);
			 printf("FW %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 < leftProperReadPos2);
			 }
			 }
			 if (contig_j->reads[k]->cEnd > contigJCutPos2 && contigJCutPos2 < contig_j->len)
			 {
			 rightProperRead = abs(contig_j->reads[k]->id);
			 rightProperReadPos2 = contig_j->reads[k]->end;
			 if (contig_j->reads[k]->beg > contig_j->reads[k]->end) // reverse complement orientation
			 {
			 rightProperReadPos1 = contig_j->reads[k]->end + (contig_j->reads[k]->cEnd - contigJCutPos2);
			 printf("FW %d c(%d, %d) r(%d,%d) right(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 > rightProperReadPos2);
			 }
			 else // forward orientation
			 {
			 rightProperReadPos1 = contig_j->reads[k]->end - (contig_j->reads[k]->cEnd - contigJCutPos2);
			 printf("RV %d c(%d, %d) r(%d,%d) right(%d, %d, %d)\n", contig_j->reads[k]->id, contig_j->reads[k]->cBeg,
			 contig_j->reads[k]->cEnd, contig_j->reads[k]->beg, contig_j->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 < rightProperReadPos2);
			 }
			 break;
			 }
			 }

			 SplitEvent *split = contig_j->split + contig_j->numSplits;
			 split->duplicRead = aread; // TODO do we need thise aread info
			 split->contigPos1 = contigJCutPos1;
			 split->contigPos2 = contigJCutPos2;
			 split->type = stype;

			 split->leftProperRead = leftProperRead;
			 split->leftProperReadPos1 = leftProperReadPos1;
			 split->leftProperReadPos2 = leftProperReadPos2;

			 split->rightProperRead = rightProperRead;
			 split->rightProperReadPos1 = rightProperReadPos1;
			 split->rightProperReadPos2 = rightProperReadPos2;

			 printSplitEvent(contig_j, split);
			 contig_j->numSplits++;
			 }

			 // add split position to contig
			 if (contig_k->maxSplits == contig_k->numSplits)
			 {
			 contig_k->maxSplits = (contig_k->maxSplits * 1.2) + 2;
			 contig_k->split = (SplitEvent*) realloc(contig_k->split, sizeof(SplitEvent) * contig_k->maxSplits);
			 assert(contig_k->split != NULL);
			 }
			 assert(contigKCutPos1 >= 0 && contigKCutPos1 < contigKCutPos2 && contigKCutPos2 <= contig_k->len);
			 {
			 // find contig reads and their positions for corresponding cutPositions
			 int leftProperRead, rightProperRead;
			 leftProperRead = rightProperRead = -1;
			 int leftProperReadPos1, leftProperReadPos2, rightProperReadPos1, rightProperReadPos2;
			 leftProperReadPos1 = leftProperReadPos2 = rightProperReadPos1 = rightProperReadPos2 = 0;
			 for (k = 0; k < contig_k->numContigReadIds; k++)
			 {
			 if (leftProperRead < 0 && contig_k->reads[k]->cEnd >= contigKCutPos1 && contigKCutPos1 > 0)
			 {
			 leftProperRead = abs(contig_k->reads[k]->id);
			 leftProperReadPos1 = contig_k->reads[k]->beg;
			 if (contig_k->reads[k]->beg > contig_k->reads[k]->end) // reverse complement orientation
			 {
			 leftProperReadPos2 = contig_k->reads[k]->end + (contig_k->reads[k]->cEnd - contigKCutPos1);
			 printf("RV %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 > leftProperReadPos2);
			 }
			 else // forward orientation
			 {
			 leftProperReadPos2 = contig_k->reads[k]->end - (contig_k->reads[k]->cEnd - contigKCutPos1);
			 printf("FW %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, leftProperRead, leftProperReadPos1,
			 leftProperReadPos2);
			 assert(leftProperReadPos1 < leftProperReadPos2);
			 }
			 }
			 if (contig_k->reads[k]->cEnd > contigKCutPos2 && contigKCutPos2 < contig_k->len)
			 {
			 rightProperRead = abs(contig_k->reads[k]->id);
			 rightProperReadPos2 = contig_k->reads[k]->end;
			 if (contig_k->reads[k]->beg > contig_k->reads[k]->end) // reverse complement orientation
			 {
			 rightProperReadPos1 = contig_k->reads[k]->end + (contig_k->reads[k]->cEnd - contigKCutPos2);
			 printf("RV %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 > rightProperReadPos2);
			 }
			 else // forward orientation
			 {
			 rightProperReadPos1 = contig_k->reads[k]->end - (contig_k->reads[k]->cEnd - contigKCutPos2);
			 printf("FW %d c(%d, %d) r(%d,%d) left(%d, %d, %d)\n", contig_k->reads[k]->id, contig_k->reads[k]->cBeg,
			 contig_k->reads[k]->cEnd, contig_k->reads[k]->beg, contig_k->reads[k]->end, rightProperRead, rightProperReadPos1,
			 rightProperReadPos2);
			 assert(rightProperReadPos1 < rightProperReadPos2);
			 }
			 break;
			 }
			 }

			 SplitEvent *split = contig_k->split + contig_k->numSplits;
			 split->duplicRead = aread; // TODO do we need thise aread info
			 split->contigPos1 = contigKCutPos1;
			 split->contigPos2 = contigKCutPos2;
			 split->type = stype;

			 split->leftProperRead = leftProperRead;
			 split->leftProperReadPos1 = leftProperReadPos1;
			 split->leftProperReadPos2 = leftProperReadPos2;

			 split->rightProperRead = rightProperRead;
			 split->rightProperReadPos1 = rightProperReadPos1;
			 split->rightProperReadPos2 = rightProperReadPos2;

			 printSplitEvent(contig_k, split);
			 contig_k->numSplits++;
			 }*/

		}

	}
}

void createSplitEvent(AnalyzeContext * actx, Contig *contig, int read, int sType)
{
	int contigCutPos1, contigCutPos2;

	contigCutPos1 = contigCutPos2 = -1;

	getFullReadPositionInContig(actx, contig, read, &contigCutPos1, &contigCutPos2);

	assert(contigCutPos1 != -1 && contigCutPos2 != -1);

	// add split position to contig
//	if (contig->maxSplits == contig->numSplits)
//	{
//		contig->maxSplits = (contig->maxSplits * 1.2) + 2;
//		contig->split = (SplitEvent*) realloc(contig->split, sizeof(SplitEvent) * contig->maxSplits);
//		assert(contig->split != NULL);
//	}
	assert(contigCutPos1 >= 0 && contigCutPos1 < contigCutPos2 && contigCutPos2 <= contig->len);
	{
		// find contig reads and their positions for corresponding cutPositions
		int leftProperRead, rightProperRead;
		leftProperRead = rightProperRead = -1;
		int leftProperReadPos1, leftProperReadPos2, rightProperReadPos1, rightProperReadPos2;
		leftProperReadPos1 = leftProperReadPos2 = rightProperReadPos1 = rightProperReadPos2 = 0;

		///// TODO fix this
//      int k;
//		for (k = 0; k < contig->numContigReadIds; k++)
//		{
//			if (leftProperRead < 0 && contig->reads[k]->cEnd >= contigCutPos1 && contigCutPos1 > 0)
//			{
//				leftProperRead = abs(contig->reads[k]->id);
//				leftProperReadPos1 = contig->reads[k]->beg;
//				if (contig->reads[k]->beg > contig->reads[k]->end) // reverse complement orientation
//				{
//					leftProperReadPos2 = contig->reads[k]->end + (contig->reads[k]->cEnd - contigCutPos1);
//					assert(leftProperReadPos1 > leftProperReadPos2);
//				}
//				else // forward orientation
//				{
//					leftProperReadPos2 = contig->reads[k]->end - (contig->reads[k]->cEnd - contigCutPos1);
//					assert(leftProperReadPos1 < leftProperReadPos2);
//				}
//			}
//			if (contig->reads[k]->cEnd > contigCutPos2 && contigCutPos2 < contig->len)
//			{
//				rightProperRead = abs(contig->reads[k]->id);
//				rightProperReadPos2 = contig->reads[k]->end;
//				if (contig->reads[k]->beg > contig->reads[k]->end) // reverse complement orientation
//				{
//					rightProperReadPos1 = contig->reads[k]->end + (contig->reads[k]->cEnd - contigCutPos2);
//					assert(rightProperReadPos1 > rightProperReadPos2);
//				}
//				else // forward orientation
//				{
//					rightProperReadPos1 = contig->reads[k]->end - (contig->reads[k]->cEnd - contigCutPos2);
//					assert(rightProperReadPos1 < rightProperReadPos2);
//				}
//				break;
//			}
//		}

//		SplitEvent *split = contig->split + contig->numSplits;
//		split->duplicRead = read; // TODO do we need thise aread info
//		split->contigPos1 = contigCutPos1;
//		split->contigPos2 = contigCutPos2;
//		split->type = sType;
//
//		split->leftProperRead = leftProperRead;
//		split->leftProperReadPos1 = leftProperReadPos1;
//		split->leftProperReadPos2 = leftProperReadPos2;
//
//		split->rightProperRead = rightProperRead;
//		split->rightProperReadPos1 = rightProperReadPos1;
//		split->rightProperReadPos2 = rightProperReadPos2;
//
//		contig->numSplits++;
	}
}

void finalContigValidation(AnalyzeContext *actx)
{
	// reset CONTIG_READ_CHECKED flag as it will be reused
	int i, j;
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;
		for (j = 0; j < contig->numcReads; j++)
			contig->cReads[j].type &= ~( CONTIG_READ_CHECKED);
	}

	int MIN_OUT_COV = 3;
	int MIN_COM_OUT_COV = 1;

	// first: check coverage drops in all contigs, and mark contig read as putative "BREAK"
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;
		int j;
		for (j = 0; j < contig->numcReads - 1; j++)
		{
			ContigRead *cread = contig->cReads + j;
			if (cread->nOutReads < MIN_OUT_COV || cread->nComReadsWithNextContigRead < MIN_COM_OUT_COV)
			{
				printf(
						"CONTIG %d FOUND LOW  COVERAGE CP: read %d, # in %d, # out %d, #comRead %d, cPos [%d, %d], rPos [%d, %d], avgCov %.2f, covDrop %d\n",
						contig->idx, cread->patchedID, cread->nInReads, cread->nOutReads, cread->nComReadsWithNextContigRead,
						cread->patchedContigPosBeg, cread->patchedContigPosEnd, cread->ovlReads->beg,
						cread->ovlReads->end, cread->avgCov, cread->lowCov);

				cread->type |= (CONTIG_READ_BREAK_UNKNOWN);
			}
			cread->type |= (CONTIG_READ_CHECKED);
		}
	}

	// todo really necessary?
	for (i = 0; i < actx->numContigs; i++)
		actx->contigs[i].flag &= ~( CONTIG_VISITED);

	int k;
	Contig **secondaryContigs;
	int maxSecondaryContigs = 10;
	int curSecondaryContigs = 0;
	secondaryContigs = malloc(sizeof(Contig*) * maxSecondaryContigs);
	assert(secondaryContigs != NULL);

	// second: full analysis of contigs, in each step consider all contigs that are somehow related to each other
	// 2. loop over length-sorted contigs until all contigs are visited:
	//        a) take next longest contig
	//        b) add all contigs that are related to 2a)
	//        c) analyze all relations of contigs in 2b): heterozygous, repeat induced, false joins, and split contigs if necessary
	for (i = 0; i < actx->numContigs; i++)
	{
		if (actx->contigs_sorted[i]->flag & CONTIG_VISITED)
		{
			continue;
		}

		Contig * primaryContig = actx->contigs_sorted[i];
		curSecondaryContigs = 0;

		// check all reads, that open a junction
		for (j = 0; j < primaryContig->numcReads; j++)
		{
			int aread = abs(primaryContig->cReads[j].patchedID);

			int nOcc = 0;
			// check only reads that were used in touring
			while (nOcc < actx->maxVReadMask && actx->vreadMask[nOcc][aread] < 0)
			{
				nOcc++;
			}

			if (nOcc > 1)
			{
				for (k = 0; k < nOcc; k++)
				{
					int contigIdx = abs(actx->vreadMask[k][aread]) - 1;
					if (contigIdx == primaryContig->idx)
						continue;

					int add = 1;
					int l;
					for (l = 0; l < curSecondaryContigs; l++)
					{
						if (secondaryContigs[l]->idx == contigIdx)
						{
							add = 0;
							break;
						}
					}
					if (add)
					{
						if (curSecondaryContigs == maxSecondaryContigs)
						{
							maxSecondaryContigs = (maxSecondaryContigs * 1.2 + 5);
							secondaryContigs = (Contig**) realloc(secondaryContigs, sizeof(Contig*) * maxSecondaryContigs);
							assert(secondaryContigs != NULL);
						}
						secondaryContigs[curSecondaryContigs] = actx->contigs + contigIdx;
						secondaryContigs[curSecondaryContigs]->tmp = primaryContig->cReads[j].patchedContigPosBeg;
						curSecondaryContigs++;
					}
				}
			}
		}
		// check for remaining repeat-induced contigs, that don't share a common contig read
		// check contig overlap groups
//		for (j = 0; j < primaryContig->numOvlGrps; j++)
//		{
//			if (primaryContig->ovlGrps[j]->flag & OVLGRP_DISCARD)
//				continue;
//
//			if (primaryContig->ovlGrps[j]->flag & OVLGRP_AREAD_IS_CONTAINED)
//			{
//				int add = 1;
//				for (k = 0; k < curSecondaryContigs; k++)
//				{
//					if (secondaryContigs[k]->idx == primaryContig->ovlGrps[j]->bread)
//					{
//						add = 0;
//						break;
//					}
//				}
//				if (add)
//				{
//					if (curSecondaryContigs == maxSecondaryContigs)
//					{
//						maxSecondaryContigs = (maxSecondaryContigs * 1.2 + 5);
//						secondaryContigs = (Contig**) realloc(secondaryContigs, sizeof(Contig*) * maxSecondaryContigs);
//						assert(secondaryContigs != NULL);
//					}
//					secondaryContigs[curSecondaryContigs] = actx->contigs + primaryContig->ovlGrps[j]->bread;
//					secondaryContigs[curSecondaryContigs]->tmp = primaryContig->ovlGrps[j]->first_abpos;
//					curSecondaryContigs++;
//				}
//			}
//		}

		// check contig graph classification
		for (j = 0; j < primaryContig->numGClassific; j++)
		{
			if (primaryContig->gClassific[j].flag & CONTIG_HAS_CONTAINED)
			{
				int add = 1;
				for (k = 0; k < curSecondaryContigs; k++)
				{
					if (secondaryContigs[k]->idx == primaryContig->gClassific[j].correspID)
					{
						add = 0;
						break;
					}
				}
				if (add)
				{
					if (curSecondaryContigs == maxSecondaryContigs)
					{
						maxSecondaryContigs = (maxSecondaryContigs * 1.2 + 5);
						secondaryContigs = (Contig**) realloc(secondaryContigs, sizeof(Contig*) * maxSecondaryContigs);
						assert(secondaryContigs != NULL);
					}
					secondaryContigs[curSecondaryContigs] = actx->contigs + primaryContig->gClassific[j].correspID;
					secondaryContigs[curSecondaryContigs]->tmp = primaryContig->gClassific[j].coveredIntervals[0];
					curSecondaryContigs++;
				}
			}
		}

		// sort secondary contigs according to mapping position in primary contig
		qsort(secondaryContigs, curSecondaryContigs, sizeof(Contig*), cmp_contig_bytmp);

		printf("Primary Contig: %d LEN %d class %d flag %d\n", primaryContig->idx, primaryContig->len, primaryContig->class, primaryContig->flag);

		// analyze contig relations
		for (j = 0; j < curSecondaryContigs; j++)
		{
			Contig *secondaryContig = secondaryContigs[j];
			printf("Secondary Contig: %d LEN %d start at: %d class %d flag %d\n", secondaryContig->idx, secondaryContig->len, secondaryContig->tmp,
					secondaryContig->class, secondaryContig->flag);
		}

		primaryContig->flag |= CONTIG_VISITED;
		for (j = 0; j < curSecondaryContigs; j++)
			secondaryContigs[j]->flag |= CONTIG_VISITED;
	}

	// clean up temp contig buffer
	free(secondaryContigs);

}

//void analyzeRepeatContent(AnalyzeContext *actx)
//{
//	int i;
//
//	for (i = 0; i < actx->numContigs; i++)
//	{
//		Contig *contig = actx->contigs + i;
//
////		printf("analyzeRepeatContent of contig %d (%d), len %d", i, contig->idx, contig->len);
//
//		track_anno* ranno = (track_anno*) (actx->contigRepeats_track->anno);
//		track_data* rdata = (track_data*) (actx->contigRepeats_track->data);
//
//		track_anno ob = ranno[i];
//		track_anno oe = ranno[i + 1];
//
//		// no repeats, nothing to do
//
//		if (ob >= oe)
//		{
////			printf(" repeat 0 0 0 REALLY?\n");
//			continue;
//		}
//
//		int nrepeat = 0;
//		int repMerged = 0;
//		int repBases = 0;
//		int MERGE_BASES = 50;
//		int tipL = 0;
//		int tipR = 0;
//
//		ob /= sizeof(track_data);
//		oe /= sizeof(track_data);
//
//		track_anno numIntervals = oe - ob;
//		if (actx->VERBOSE > 2)
//			printf(" number of repeats: (%llu)\n", numIntervals);
//
//		int prvb = rdata[ob];
//		int prve = rdata[ob + 1];
//		ob += 2;
//
//		while (ob < oe)
//		{
//			int b = rdata[ob];
//			int e = rdata[ob + 1];
//
//			if (actx->VERBOSE > 2)
//				printf(" [%d, %d]\n", b, e);
//
//			if (b - prve < MERGE_BASES)
//			{
//				// merge intervals
//				if (actx->VERBOSE > 2)
//					printf("merge intervals : [%d,%d] [%d,%d]", prvb, prve, b, e);
//				prve = e;
//				repMerged++;
//				if (actx->VERBOSE > 2)
//					printf(" --> [%d,%d] (merged: %d)\n", prvb, prve, repMerged);
//			}
//			else
//			{
//				nrepeat++;
//				repBases += prve - prvb;
//				if (prvb < 25000)
//				{
//					tipL += MIN(25000, prve) - prvb;
//				}
//				if (prve > contig->len - 25000)
//				{
//					tipR += prve - MAX(prvb, contig->len - 25000);
//				}
//				prvb = b;
//				prve = e;
//
//			}
//			ob += 2;
//		}
//		nrepeat++;
//		repBases += prve - prvb;
//
//		if (prvb < 25000)
//		{
//			tipL += MIN(25000, prve) - prvb;
//		}
//		if (prve > contig->len - 25000)
//		{
//			tipR += prve - MAX(prvb, contig->len - 25000);
//		}
//
//		contig->repBases = repBases;
//		contig->repBasesTipLeft = tipL;
//		contig->repBasesTipRight = tipR;
//		contig->numRepeats = nrepeat;
////		printf(" repeat %d %d %d (%.2f%%) tipL %d, tipR %d\n", nrepeat, repBases, repMerged, repBases * 100.0 / contig->len, tipL, tipR);
//	}
//}

void getCorrespondingPositionInARead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2)
{
	int cBeg = *pos1;
	int cEnd = *pos2;

	int abeg = ovl->path.abpos;
	int aend = ovl->path.aepos;

	int bbeg = ovl->path.bbpos;
	int bend = ovl->path.bepos;

	int twidth = actx->twidth;
#ifdef DEBUG_STEP1A
	printf("getCorrespondingPositionInARead %d x %d, c(%d, %d), a(%d, %d) b(%d, %d)", ovl->aread, ovl->bread, cBeg, cEnd, abeg, aend, bbeg, bend);
	if (ovl->flags & OVL_COMP)
	printf(" COMP");
	printf("\n");
#endif

#ifdef DEBUG_STEP1A
	printf("comp: bbeg: %d, bend: %d\n", bbeg, bend);
#endif
	// find corresponding begin
	if (cBeg <= bbeg)
	{
#ifdef DEBUG_STEP1A
		printf("cBeg <= bbeg, %d <= %d\n", cBeg, bbeg);
#endif
		*pos1 = abeg;
	}
	else
	{
		if (ovl->path.tlen)
		{
#ifdef DEBUG_STEP1A
			printf("ccBeg: %d\n", cBeg);
#endif
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = ovl->path.bbpos;

			int j = 0;
			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\n", apos, bpos);
#endif
				if (bpos + twidth / 2 > cBeg)
				{
					apos += cBeg - bpos;
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\n", apos, bpos);
#endif
					*pos1 = -apos;
					break;
				}
				j += 2;
			}
		}
		else
		{
			if ((bend - cBeg) < (cBeg - bbeg))
			{
				*pos1 = -(abeg + (bend - cBeg));
			}
			else
			{
				*pos1 = -(aend - (cBeg - bbeg));
			}
		}

	}
	// find corresponding end
	if (cEnd >= bend)
	{
		*pos2 = aend;
	}
	else
	{
		if (ovl->path.tlen)
		{
#ifdef DEBUG_STEP1A
			printf("ccBeg: %d\n", cEnd);
#endif
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = ovl->path.bbpos;

			int j = 0;
			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\n", apos, bpos);
#endif
				if (bpos + twidth / 2 > cEnd)
				{
					apos += cEnd - bpos;
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\n", apos, bpos);
#endif

					*pos2 = -apos;
					break;
				}
				j += 2;
			}
		}
		else
		{
			if ((bend - cEnd) < (cEnd - bbeg))
			{
				*pos2 = -(abeg + (bend - cEnd));
			}
			else
			{
				*pos2 = -(aend - (cEnd - bbeg));
			}
		}
	}

#ifdef DEBUG_STEP1A
	printf("final range: %d, %d\n", *pos1, *pos2);
#endif
}

void getCorrespondingPositionInBRead(AnalyzeContext *actx, Overlap *ovl, int* pos1, int *pos2)
{
	int cBeg = *pos1;
	int cEnd = *pos2;

	int abeg = ovl->path.abpos;
	int aend = ovl->path.aepos;

	int bbeg = ovl->path.bbpos;
	int bend = ovl->path.bepos;

	int twidth = actx->twidth;
	int transCBeg;
	int transCEnd;

#ifdef DEBUG_STEP1A
	{
		printf("getCorrespondingPositionInBRead %d x %d, c(%d, %d), a(%d, %d) b(%d, %d)", ovl->aread, ovl->bread, cBeg, cEnd, abeg, aend, bbeg, bend);
		if (ovl->flags & OVL_COMP)
		printf(" COMP");
		printf("\n");
	}
#endif
// find begin
	if (cBeg <= abeg)
	{
		transCBeg = bbeg;
#ifdef DEBUG_STEP1A
		printf("begin: %d <= %d, set transCBeg to %d\n", cBeg, abeg, transCBeg);
#endif
	}
	else
	{
		if (ovl->path.tlen)
		{
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = bbeg;

			int bb = 0;

			int j = 0;

			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];

				bb = bpos;
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif
				if (apos + twidth / 2 > cBeg)
				{
					bb += (cBeg - apos);
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif
					transCBeg = -bb;
					break;
				}
				j += 2;
			}
		}
		else
		{
#ifdef DEBUG_STEP1A
			printf("path.tlen is 0. set transCBeg\n");
#endif
			if ((cBeg - abeg) < (aend - cBeg))
			{
				transCBeg = -(bbeg + (cBeg - abeg));
			}
			else
			{
				transCBeg = -(bend - (aend - cBeg));
			}
		}

	}

	// find end
	if (cEnd >= aend)
	{
#ifdef DEBUG_STEP1A
		printf("end: %d >0 %d\n", cEnd, aend);
#endif
		transCEnd = bend;
	}
	else
	{
		if (ovl->path.tlen)
		{
			ovl_trace* trace = ovl->path.trace;
			int apos = abeg;
			int bpos = bbeg;

			int bb = 0;

			int j = 0;

			while (j < ovl->path.tlen)
			{
				apos = (apos / twidth + 1) * twidth;
				bpos += trace[j + 1];

				bb = bpos;
#ifdef DEBUG_STEP1A
				printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif

				if (apos + twidth / 2 > cEnd)
				{
					bb += (cEnd - apos);
#ifdef DEBUG_STEP1A
					printf("apos %6d, bpos %6d\t bb %6d\n", apos, bpos, bb);
#endif
					transCEnd = -bb;
					break;
				}
				j += 2;
			}
		}
		else
		{
#ifdef DEBUG_STEP1A
			printf("transCEnd path.tlen is 0\n");
#endif

			if ((aend - cEnd) < (cEnd - abeg))
			{
				transCEnd = -(bend - (aend - cEnd));
			}
			else
			{
				transCEnd = -(bbeg + (cEnd - abeg));
			}

		}
	}

	if (ovl->flags & OVL_COMP)
	{
#ifdef DEBUG_STEP1A
		printf("switch positions transCBeg: %d, transCEnd: %d\n", transCBeg, transCEnd);
#endif
		*pos1 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - abs(transCBeg);
		*pos2 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - abs(transCEnd);

		if (transCBeg < 0)
			(*pos1) = -(*pos1);

		if (transCEnd < 0)
			(*pos2) = -(*pos2);
	}
	else
	{
		*pos1 = transCBeg;
		*pos2 = transCEnd;
	}
#ifdef DEBUG_STEP1A
	printf("getCorrespondingPositionInBRead: final range: %d, %d\n", *pos1, *pos2);
#endif
}

int processFixedReadOverlap_handler1(void* _ctx, Overlap* ovls, int novl)
{
	AnalyzeContext* actx = (AnalyzeContext*) _ctx;

	int i;

	int aread = ovls->aread;
	int bread;

	for (i = 0; i < novl; i++)
	{
		Overlap *ovl = ovls + i;
		bread = ovl->bread;

		if (ovl->flags & OVL_DISCARD)
		{
#ifdef DEBUG_STEP1A
			printf("skip overlap %d x %d, DISCARD\n", aread, bread);
#endif
			continue;
		}

		if (aread == bread)
		{
#ifdef DEBUG_STEP1A
			printf("skip identity overlap %d == %d\n", aread, bread);
#endif
			continue;
		}

		/// TODO ugly hack, to ignore reads that overlap multiple times with the current a-read, actually they should be removed in LAfilter, why are they still there????
		{
			int j = i;
			while (j + 1 < novl && (ovl->bread == ovls[j + 1].bread))
			{
				printf("ovl[%d] and ovl[%d] have same bread %d? 1(%d vs %d [%d, %d][%d,%d]) and 2(%d vs %d [%d, %d][%d,%d])\n", i, j + 1, ovl->bread, ovl->aread,
						ovl->bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, ovls[j + 1].aread, ovls[j + 1].bread, ovls[j + 1].path.abpos,
						ovls[j + 1].path.aepos, ovls[j + 1].path.bbpos, ovls[j + 1].path.bepos);
				j++;
			}

			if (j > i)
			{
#ifdef DEBUG_STEP1A
				printf("FUCK skip read %d. maps multiple times (%d) with current read %d\n", aread, j-i, bread);
#endif
				i = j;
				continue;
			}
		}

		int oLen = MAX(ovl->path.aepos - ovl->path.abpos, ovl->path.bepos - ovl->path.bbpos);

		if (oLen < actx->min_olen)
		{
#ifdef DEBUG_STEP1A
			printf("skip overlap %d x %d, len: %d\n", aread, bread, oLen);
#endif
			continue;
		}

		if (DB_READ_LEN(actx->patchedReadDB, ovl->bread) < actx->min_rlen)
		{
#ifdef DEBUG_STEP1A
			printf("skip overlap %d x %d, rlen: %d\n", aread, bread, DB_READ_LEN(actx->patchedReadDB, ovl->bread));
#endif
			continue;
		}

		// check if a-read is present in any contig, if yes than add corresponding breads to the contigs
		int j, k, l;
		for (j = 0; actx->vreadMask[j][aread] < 0; j++)
		{
			int conId = abs(actx->vreadMask[j][aread]);

			Contig * contig = actx->contigs + (conId - 1);

#ifdef DEBUG_STEP1A
			printf("found aread %d in contig %d: (len: %d, rep: %d) bread %d\n", aread, conId - 1, contig->len,
					getRepeatBasesFromInterval(actx, contig->idx, 0, contig->len), bread);
#endif
			// check if corresponding b-read is already visited
			int found = 0;
			for (k = 0; k < actx->maxVReadMask; k++)
			{
				if (abs(actx->vreadMask[k][bread]) == conId)
					found = 1;
				if (actx->vreadMask[k][bread] == 0)
					break;
			}

			if (k == actx->maxVReadMask)
			{
				//printf("TERROR b-read %d is part of %d contigs!!!!\n", bread, k);
				l = 0;
				while (l < k)
				{
					//	printf("read %d is part of contig %d\n", bread, actx->vreadMask[l][bread]);
					l++;
				}
				//fflush(stdout);
				actx->maxVReadMask++;
				actx->vreadMask = (int**) realloc(actx->vreadMask, sizeof(int*) * actx->maxVReadMask);
				actx->vreadMask[k] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
				bzero(actx->vreadMask[k], sizeof(int) * DB_NREADS(actx->patchedReadDB));
			}

			// which coordinates of the a-read are used within the contig
			ContigRead *cread = getContigRead(contig, aread);

			assert(cread != NULL);

			int areadContigBeg = cread->patchedContigPosBeg;
			int areadContigEnd = cread->patchedContigPosEnd;

			int areadReadBeg = cread->ovlReads->beg;
			int areadReadEnd = cread->ovlReads->end;

#ifdef DEBUG_STEP1A
			printf("OVERLAP: %d vs %d [%8d, %8d] %c [%8d, %8d] l [%8d, %8d]\n",ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos, (ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos,
					DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread));
			printf("Aread %d in contig %d: cPos [%d, %d] rPos: [%d, %d] %c\n", aread, contig->idx, areadContigBeg, areadContigEnd, areadReadBeg, areadReadEnd, (areadReadBeg > areadReadEnd) ? 'C' : 'N');
			printf("convert overlap coordinates from %d %d [%d, %d] %c [%d, %d]\n", ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos,
					(ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos);
#endif
			// complement sequence
			int cReadInComplement = 0;
			if (areadReadBeg > areadReadEnd)
			{
				int tmp = areadReadBeg;
				areadReadBeg = areadReadEnd;
				areadReadEnd = tmp;
				cReadInComplement = 1;
			}

			// do the a-read contig coordinates have an intersection with the b-read?
			int intersection = intersect(areadReadBeg, areadReadEnd, ovl->path.abpos, ovl->path.aepos);
			if (intersection == 0)
			{
#ifdef DEBUG_STEP1A
				printf("skip overlap %d x %d, No intersection\n", aread, bread);
#endif
				continue;
			}

			// mark b-read for current contig as visited
			if (!found)
				actx->vreadMask[k][bread] = conId;

			// check if OvlRead is already included due to transitivity
			int check = 1;
			{
				int u;

				for (u = 0; u < cread->numOvlReads; u++)
				{
					OvlRead *ovlU = cread->ovlReads + u;
					if (ovlU->patchedID == bread)
					{
						check = 0;
						break;
					}
				}
			}

			if (check)
			{
				// add valid bread to set of OvlReads (--> proper read for correction)
				if (cread->numOvlReads == cread->maxOvlReads)
				{
					cread->maxOvlReads = (int) (cread->maxOvlReads * 1.2) + 10;
					cread->ovlReads = (OvlRead*) realloc(cread->ovlReads, cread->maxOvlReads * sizeof(OvlRead));

					if (cread->ovlReads == NULL)
					{
						fprintf(stderr, "[ERROR]: Unable to increase OvlRead buffer of contig %d at read: %d, from %lu to %lu!\n", conId, aread,
								cread->numOvlReads * sizeof(OvlRead), cread->maxOvlReads * sizeof(OvlRead));
						exit(1);
					}
				}

				OvlRead *curOvlRead = cread->ovlReads + cread->numOvlReads;
				curOvlRead->patchedID = bread;

				int pos1, pos2;
				int abeg, aend;

				abeg = ovl->path.abpos;
				aend = ovl->path.aepos;

				pos1 = ovl->path.bbpos;
				pos2 = ovl->path.bepos;

				if (abeg < areadReadBeg)
					abeg = areadReadBeg;

				if (aend > areadReadEnd)
					aend = areadReadEnd;

				// contig read sticks in reverse complements in contig
				if (cReadInComplement)
				{
					if (aend > areadReadEnd)
						curOvlRead->cBeg = areadContigBeg;
					else
						curOvlRead->cBeg = areadContigBeg + (areadReadEnd - aend);

					if (abeg < areadReadBeg)
						curOvlRead->cEnd = areadContigEnd;
					else
						curOvlRead->cEnd = areadContigBeg + (areadReadEnd - abeg);
				}
				else
				{
					if (abeg < areadReadBeg)
						curOvlRead->cBeg = areadContigBeg;
					else
						curOvlRead->cBeg = areadContigBeg + (abeg - areadReadBeg);

					if (aend > areadReadEnd)
						curOvlRead->cEnd = areadContigEnd;
					else
						curOvlRead->cEnd = areadContigEnd - (areadReadEnd - aend);

				}

				if (abeg != ovl->path.abpos || aend != ovl->path.aepos)
				{
					pos1 = abeg;
					pos2 = aend;
#ifdef DEBUG_STEP1A
					printf("getCorrespondingPositionInBRead(actx, ovl, &pos1, &pos2);\n");
					getCorrespondingPositionInBRead(actx, ovl, &pos1, &pos2);
#endif
				}

#ifdef DEBUG_STEP1A
				printf("set pos1 to %d and pos2 to %d\n", pos1, pos2);
#endif

				if (cReadInComplement)
				{
					curOvlRead->beg = pos2;
					curOvlRead->end = pos1;
				}
				else
				{
					curOvlRead->beg = pos1;
					curOvlRead->end = pos2;
				}
				cread->numOvlReads++;

				assert(curOvlRead->cBeg < curOvlRead->cEnd);
				assert(curOvlRead->cBeg >= cread->patchedContigPosBeg);
				assert(curOvlRead->cEnd <= cread->patchedContigPosEnd);

#ifdef DEBUG_STEP1A
				ContigRead *falk = getContigRead(contig, aread);
				assert(falk != NULL);
				printf("add overlap %d x %d, a(%d, %d), b(%d, %d), l(%d, %d), %c at ConPos [%d, %d] in Read [%d, %d]\n", aread, bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos,
						ovl->path.bepos, DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread), (ovl->flags & OVL_COMP) ? 'c' : 'n',
						falk->patchedContigPosBeg, falk->patchedContigPosEnd, falk->ovlReads->beg, falk->ovlReads->end);

				printf("=======<<<< added OVL-B-read: %d %d-%d, in c%d-%d-%d %c %c crLen %d oLen %d ocrLen %d ocLen %d\n", curOvlRead->patchedID, curOvlRead->beg, curOvlRead->end, contig->idx, curOvlRead->cBeg, curOvlRead->cEnd,
						(cReadInComplement) ? 'C': 'N', (ovl->flags & OVL_COMP) ? 'C': 'N', cread->ovlReads->cEnd - cread->ovlReads->cBeg, ovl->path.aepos - ovl->path.abpos, curOvlRead->cEnd - curOvlRead->cBeg, abs(abs(curOvlRead->beg)-abs(curOvlRead->end)));

#endif
			}

		}

		/// As there is no proper transitivity given, we have to check if one of the b-reads from all overlaps of the current
		/// aread is part of a contig as well, and was not already incorporated as valid overlap
		// check if b-read is present in any contig, if yes than add corresponding breads to the contigs
		for (j = 0; actx->vreadMask[j][bread] < 0; j++)
		{

			int conId = abs(actx->vreadMask[j][bread]);

			Contig * contig = actx->contigs + (conId - 1);

#ifdef DEBUG_STEP1A
			printf("found bread %d in contig %d: (len: %d, rep: %d)\n", bread, conId - 1, contig->len,
					getRepeatBasesFromInterval(actx, contig->idx, 0, contig->len));
#endif
			// check if corresponding a-read is already visited
			int found = 0;
			for (k = 0; k < actx->maxVReadMask; k++)
			{
				if (abs(actx->vreadMask[k][aread]) == conId)
					found = 1;
				if (actx->vreadMask[k][aread] == 0)
					break;
			}

			if (k == actx->maxVReadMask)
			{
				//printf("TERROR a-read %d is part of %d contigs!!!!\n", aread, k);
				l = 0;
				while (l < k)
				{
					//printf("read %d is part of contig %d\n", aread, actx->vreadMask[l][aread]);
					l++;
				}
				//fflush(stdout);
				actx->maxVReadMask++;
				actx->vreadMask = (int**) realloc(actx->vreadMask, sizeof(int*) * actx->maxVReadMask);
				actx->vreadMask[k] = (int *) malloc(sizeof(int) * DB_NREADS(actx->patchedReadDB));
				bzero(actx->vreadMask[k], sizeof(int) * DB_NREADS(actx->patchedReadDB));
			}

			ContigRead *cread = getContigRead(contig, bread);
			assert(cread != NULL);

			int breadContigBeg = cread->patchedContigPosBeg;
			int breadContigEnd = cread->patchedContigPosEnd;

			int breadReadBeg = cread->ovlReads->beg;
			int breadReadEnd = cread->ovlReads->end;

#ifdef DEBUG_STEP1A
			printf("OVERLAP: %d vs %d [%8d, %8d] %c [%8d, %8d] l [%8d, %8d]\n",ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos, (ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos,
					DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread));
			printf("Bread %d in contig %d: cPos [%d, %d] rPos: [%d, %d]\n", bread, contig->idx, breadContigBeg, breadContigEnd, breadReadBeg, breadReadEnd);
			printf("convert overlap coordinates from %d %d [%d, %d] %c [%d, %d]\n", ovl->aread, ovl->bread, ovl->path.abpos, ovl->path.aepos,
					(ovl->flags & OVL_COMP) ? 'C' : 'N', ovl->path.bbpos, ovl->path.bepos);
#endif
			int bBeg, bEnd;

			int cReadInComplement = 0;

			// complement sequence
			if (breadReadBeg > breadReadEnd)
			{
				int tmp = breadReadBeg;
				breadReadBeg = breadReadEnd;
				breadReadEnd = tmp;
				cReadInComplement = 1;
			}
#ifdef DEBUG_STEP1A
			printf("cBeg: %d, cEnd: %d\n", breadContigBeg, breadContigEnd);
#endif
			if (ovl->flags & OVL_COMP)
			{
				bBeg = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - ovl->path.bepos;
				bEnd = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - ovl->path.bbpos;
			}
			else
			{
				bBeg = ovl->path.bbpos;
				bEnd = ovl->path.bepos;
			}
#ifdef DEBUG_STEP1A
			printf("get intersection from %d in contig (%d, %d) and b(%d, %d)\n", bread, breadReadBeg, breadReadEnd, bBeg, bEnd);
#endif
			int intersection = intersect(breadReadBeg, breadReadEnd, bBeg, bEnd);
			if (intersection == 0)
			{
				continue;
			}

			// mark a-read for current contig as visited
			if (!found)
				actx->vreadMask[k][aread] = conId;

			// check if we added the same OvlRead twice due to symmetry
			int check = 1;
			{
				int u;

				for (u = 0; u < cread->numOvlReads; u++)
				{
					OvlRead *ovlU = cread->ovlReads + u;
					if (ovlU->patchedID == aread)
					{
						check = 0;
						break;
					}

				}
			}

			if (check)
			{
				if (cread->numOvlReads == cread->maxOvlReads)
				{
					cread->maxOvlReads = (int) (cread->maxOvlReads * 1.2) + 10;
					cread->ovlReads = (OvlRead*) realloc(cread->ovlReads, sizeof(OvlRead) * cread->maxOvlReads);

					if (cread->ovlReads == NULL)
					{
						fprintf(stderr, "[ERROR]: Unable to increase OvlRead buffer of contig %d at read: %d, from %d to %d!\n", contig->idx, bread, cread->numOvlReads,
								cread->maxOvlReads);
						exit(1);
					}
				}

				// add valid a-read to set of OvlReads (--> proper read for correction)
				OvlRead *curOvlRead = cread->ovlReads + cread->numOvlReads;
				curOvlRead->patchedID = aread;

				int pos1, pos2;

				if (ovl->flags & OVL_COMP)
				{
					if (bEnd > breadReadEnd)
						pos1 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bEnd + (bEnd - breadReadEnd);
					else
						pos1 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bEnd;

					if (bBeg < breadReadBeg)
						pos2 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bBeg - (breadReadBeg - bBeg);
					else
						pos2 = DB_READ_LEN(actx->patchedReadDB, ovl->bread) - bBeg;

				}
				else
				{
					if (bBeg < breadReadBeg)
						pos1 = breadReadBeg;
					else
						pos1 = bBeg;

					if (bEnd > breadReadEnd)
						pos2 = breadReadEnd;
					else
						pos2 = bEnd;
				}
#ifdef DEBUG_STEP1A
				printf("set pos1 to %d and pos2 to %d\n", pos1, pos2);
#endif
				if (cReadInComplement)
				{
					if (bEnd > breadReadEnd)
						curOvlRead->cBeg = breadContigBeg;
					else
						curOvlRead->cBeg = breadContigBeg + (breadReadEnd - bEnd);

					if (bBeg < breadReadBeg)
						curOvlRead->cEnd = breadContigEnd;
					else
						curOvlRead->cEnd = breadContigBeg + (breadReadEnd - bBeg);
				}
				else
				{
					if (bBeg < breadReadBeg)
						curOvlRead->cBeg = breadContigBeg;
					else
						curOvlRead->cBeg = breadContigBeg + (bBeg - breadReadBeg);

					if (bEnd > breadReadEnd)
						curOvlRead->cEnd = breadContigEnd;
					else
						curOvlRead->cEnd = breadContigEnd - (breadReadEnd - bEnd);
				}

				getCorrespondingPositionInARead(actx, ovl, &pos1, &pos2);

				if (cReadInComplement != (ovl->flags & OVL_COMP))
				{
					curOvlRead->beg = pos2;
					curOvlRead->end = pos1;
				}
				else
				{
					curOvlRead->beg = pos1;
					curOvlRead->end = pos2;
				}

				cread->numOvlReads++;

#ifdef DEBUG_STEP1A
				ContigRead *falk = getContigRead(contig, bread);
				assert(falk != NULL);

				printf("add overlap %d x %d, a(%d, %d), b(%d, %d), l(%d, %d) %c at ConPos [%d, %d] in Read [%d, %d]\n", aread, bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos,
						ovl->path.bepos, DB_READ_LEN(actx->patchedReadDB, ovl->aread), DB_READ_LEN(actx->patchedReadDB, ovl->bread), (ovl->flags & OVL_COMP) ? 'c' : 'n',
						falk->patchedContigPosBeg, falk->patchedContigPosEnd, falk->ovlReads->beg, falk->ovlReads->end);

				printf("=======>>>> added OVL-A-read: %d %d-%d, in c%d-%d-%d %c %c crLen %d oLen %d ocrLen %d ocLen %d\n", curOvlRead->patchedID, curOvlRead->beg, curOvlRead->end, contig->idx, curOvlRead->cBeg, curOvlRead->cEnd,
						(cReadInComplement) ? 'C': 'N', (ovl->flags & OVL_COMP) ? 'C': 'N', cread->ovlReads->cEnd - cread->ovlReads->cBeg, ovl->path.aepos - ovl->path.abpos, curOvlRead->cEnd - curOvlRead->cBeg, abs(abs(curOvlRead->beg)-abs(curOvlRead->end)));
#endif
			}
		}
	}
	return 1;
}

int cmpOVLreadById(const void *a, const void *b)
{
	return ((const OvlRead*) a)->patchedID - ((const OvlRead*) b)->patchedID;
}

int getRepeatBasesFromInterval(AnalyzeContext* actx, int readID, int beg, int end)
{
	track_anno* rep_anno = actx->corContigRepeats_track->anno;
	track_data* rep_data = actx->corContigRepeats_track->data;

	track_anno rb, re;

	int repBases = 0;
	int rBeg, rEnd;

	if (readID < 0 || readID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getRepeatBasesFromInterval readID: %d out of bounds [0, %d]\n", readID, DB_NREADS(actx->corContigDB) - 1);
		fflush(stderr);
		exit(1);
	}

	// repeat bases in a-read
	rb = rep_anno[readID] / sizeof(track_data);
	re = rep_anno[readID + 1] / sizeof(track_data);

	while (rb < re)
	{
		rBeg = rep_data[rb];
		rEnd = rep_data[rb + 1];

		repBases += intersect(beg, end, rBeg, rEnd);

		rb += 2;
	}

	printf("getRepeatBasesFromInterval: %d %d %d --> %d\n", readID, beg, end, repBases);
	return repBases;
}

void updateOverlapGroup(AnalyzeContext *actx, OverlapGroup *ovlgrp, Overlap *pOvls, int *ovlIdx, int numOvlIdx)
{
	int i;
	int gap_1, gap_2;
	int its_1, its_2;
	int aread, bread;

#ifdef DEBUG_STEP2A
	int a = 30, b = 32;
#endif

	if (numOvlIdx == 0)
	{
#ifdef DEBUG_STEP2A
			printf("updateOverlapGroup: numOvlIdx == 0\n");
#endif
		return;
	}

	Overlap *ovl_1, *ovl_2;
	ovl_1 = pOvls + ovlIdx[0];

	if (ovl_1->flags & OVL_COMP)
		ovlgrp->flag |= OVLGRP_COMP;

	ovlgrp->coveredBasesInA = ovl_1->path.aepos - ovl_1->path.abpos;
	ovlgrp->coveredBasesInB = ovl_1->path.bepos - ovl_1->path.bbpos;

#ifdef DEBUG_STEP2A
		printf("ovlgrp->coveredBasesInA: %d\n", ovlgrp->coveredBasesInA);
		printf("ovlgrp->coveredBasesInB: %d\n", ovlgrp->coveredBasesInB);
#endif

	int ab1, ae1, bb1, be1;
	int ab2, ae2, bb2, be2;

	aread = ovl_1->aread;
	bread = ovl_1->bread;
	int brlen = DB_READ_LEN(actx->corContigDB, bread);

	ovlgrp->repeatBasesInA = getRepeatBasesFromInterval(actx, aread, ovl_1->path.abpos, ovl_1->path.aepos);
	if (ovl_1->flags & OVL_COMP)
		ovlgrp->repeatBasesInB = getRepeatBasesFromInterval(actx, bread, brlen - ovl_1->path.bepos, brlen - ovl_1->path.bbpos);
	else
		ovlgrp->repeatBasesInB = getRepeatBasesFromInterval(actx, bread, ovl_1->path.bbpos, ovl_1->path.bepos);

	for (i = 1; i < numOvlIdx; i++)
	{
		ovl_1 = pOvls + ovlIdx[i - 1];
		ovl_2 = pOvls + ovlIdx[i];

		assert((ovl_1->flags & OVL_COMP) == (ovl_2->flags & OVL_COMP));
		assert(ovl_1->aread == ovl_2->aread);
		assert(ovl_1->bread == ovl_2->bread);

		assert((ovlgrp->flag & OVLGRP_COMP) == (ovl_1->flags & OVL_COMP));

		ab1 = ovl_1->path.abpos;
		ae1 = ovl_1->path.aepos;
		bb1 = ovl_1->path.bbpos;
		be1 = ovl_1->path.bepos;

		ab2 = ovl_2->path.abpos;
		ae2 = ovl_2->path.aepos;
		bb2 = ovl_2->path.bbpos;
		be2 = ovl_2->path.bepos;

#ifdef DEBUG_STEP2A
			printf("ab1: %d, ae1: %d, bb1: %d, be1: %d\n", ab1, ae1, bb1, be1);
			printf("ab2: %d, ae2: %d, bb2: %d, be2: %d\n", ab2, ae2, bb2, be2);
#endif

		// update covered bases
		ovlgrp->coveredBasesInA += ae2 - ab2;
		ovlgrp->coveredBasesInB += be2 - bb2;

#ifdef DEBUG_STEP2A
			printf("ovlgrp->coveredBasesInA: %d\n", ovlgrp->coveredBasesInA);
			printf("ovlgrp->coveredBasesInB: %d\n", ovlgrp->coveredBasesInB);
#endif

		// determine overhang/gap bases between overlaps
		its_1 = intersect(ab1, ae1, ab2, ae2);
		its_2 = intersect(bb1, be1, bb2, be2);

#ifdef DEBUG_STEP2A
			printf("its_1: %d\n", its_1);
			printf("its_2: %d\n", its_2);
#endif

			// update repeat bases
			ovlgrp->repeatBasesInA += getRepeatBasesFromInterval(actx, ovl_2->aread, ae1, ae2);
			if((ovl_2->flags & OVL_COMP))
			{
				ovlgrp->repeatBasesInB += getRepeatBasesFromInterval(actx, ovl_2->bread, brlen - be2, brlen - bb2);
			}
			else
			{
				ovlgrp->repeatBasesInB += getRepeatBasesFromInterval(actx, ovl_2->bread, bb2, be2);
			}

#ifdef DEBUG_STEP2A
			printf("ovlgrp->repeatBasesInA: %d\n", ovlgrp->repeatBasesInA);
			printf("ovlgrp->repeatBasesInB: %d\n", ovlgrp->repeatBasesInB);
#endif

		gap_1 = gap_2 = 0;
		if (its_1 == 0)
		{
			gap_1 = ab2 - ae1;
			ovlgrp->gapBasesInA += gap_1;
#ifdef DEBUG_STEP2A
			printf("ovlgrp->gapBasesInA: %d\n", ovlgrp->gapBasesInA);
#endif
		}
		else    // overlap of neighboring alignments --> i.e. to avoid that overlap accounts twice remove its_1 bases and repeats from its_1 range
		{
			ovlgrp->coveredBasesInA -= its_1;
#ifdef DEBUG_STEP2A
			printf("ovlgrp->coveredBasesInA: %d (removed its1 %d)\n", ovlgrp->coveredBasesInA, its_1);
#endif

			ovlgrp->repeatBasesInA -= getRepeatBasesFromInterval(actx, ovl_2->aread, ab2, ae1);

#ifdef DEBUG_STEP2A
			printf("ovlgrp->repeatBasesInA: %d (removed its1 repeat bases %d)\n", ovlgrp->coveredBasesInA, getRepeatBasesFromInterval(actx, ovl_2->aread, ab2, ae1));
#endif
		}

		if (its_2 == 0)
		{
			gap_2 = bb2 - be1;
			ovlgrp->gapBasesInB += gap_2;

#ifdef DEBUG_STEP2A
			printf("ovlgrp->gapBasesInB: %d\n", ovlgrp->gapBasesInB);
#endif

		}
		else     // overlap of neighboring alignments --> i.e. to avoid that overlap accounts twice remove its_2 bases and repeats from its_2 range
		{
			ovlgrp->coveredBasesInB -= its_2;
#ifdef DEBUG_STEP2A
			printf("ovlgrp->coveredBasesInB: %d (removed its2 %d)\n", ovlgrp->coveredBasesInB, its_2);
#endif

			if (ovl_2->flags & OVL_COMP)
			{
				ovlgrp->repeatBasesInB -= getRepeatBasesFromInterval(actx, ovl_2->bread, brlen - be1, brlen - bb2);
#ifdef DEBUG_STEP2A
			printf("ovlgrp->repeatBasesInB: %d (removed its2 repeat bases %d)\n", ovlgrp->repeatBasesInB, getRepeatBasesFromInterval(actx, ovl_2->bread, brlen - be1, brlen - bb2));
#endif
			}
			else
			{
				ovlgrp->repeatBasesInB -= getRepeatBasesFromInterval(actx, ovl_2->bread, bb2, be1);
#ifdef DEBUG_STEP2A
			printf("ovlgrp->repeatBasesInB: %d (removed its2 repeat bases %d)\n", ovlgrp->repeatBasesInB, getRepeatBasesFromInterval(actx, ovl_2->bread, bb2, be1));
#endif
			}
		}

#ifdef DEBUG_STEP2A
		printf("%d x %d [%d, %d, %d, %d] [%d, %d, %d, %d], gap1 %d, gap_2 %d --> gapBasesInA %d, gapBasesInB %d\n", ovl_1->aread, ovl_1->bread,
				ovl_1->path.abpos, ovl_1->path.aepos, ovl_1->path.bbpos, ovl_1->path.bepos, ovl_2->path.abpos, ovl_2->path.aepos, ovl_2->path.bbpos, ovl_2->path.bepos,
				gap_1, gap_2, ovlgrp->gapBasesInA, ovlgrp->gapBasesInB);
#endif
	}

	// set positions of first and last overlaps
	ovlgrp->first_abpos = pOvls[ovlIdx[0]].path.abpos;
	ovlgrp->first_aepos = pOvls[ovlIdx[0]].path.aepos;

	ovlgrp->last_abpos = pOvls[ovlIdx[numOvlIdx - 1]].path.abpos;
	ovlgrp->last_aepos = pOvls[ovlIdx[numOvlIdx - 1]].path.aepos;

	if (ovlgrp->flag & OVLGRP_COMP)
	{
		ovlgrp->first_bbpos = brlen - pOvls[ovlIdx[0]].path.bepos;
		ovlgrp->first_bepos = brlen - pOvls[ovlIdx[0]].path.bbpos;

		ovlgrp->last_bbpos = brlen - pOvls[ovlIdx[numOvlIdx - 1]].path.bepos;
		ovlgrp->last_bepos = brlen - pOvls[ovlIdx[numOvlIdx - 1]].path.bbpos;

	}
	else
	{
		ovlgrp->first_bbpos = pOvls[ovlIdx[0]].path.bbpos;
		ovlgrp->first_bepos = pOvls[ovlIdx[0]].path.bepos;

		ovlgrp->last_bbpos = pOvls[ovlIdx[numOvlIdx - 1]].path.bbpos;
		ovlgrp->last_bepos = pOvls[ovlIdx[numOvlIdx - 1]].path.bepos;

	}
}

void analyzeContigOverlapGraph(AnalyzeContext *actx)
{
	int i, j, k;
	int tmp;
	int crComp;

	// sort b-reads according read id and remove duplicates
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;

		if (contig->flag & CONTIG_DISCARD)
			continue;

		for (j = 0; j < contig->numcReads; j++)
		{
			qsort(contig->cReads[j].ovlReads, contig->cReads[j].numOvlReads, sizeof(OvlRead), cmpOVLreadById);
		}
	}

	// todo pass variable, set by user ?
	int COV_BIN_SIZE = 500;
	int MIN_COV = 3;

	bit *inReads = ba_new(DB_NREADS(actx->patchedReadDB));
	bit *outReads = ba_new(DB_NREADS(actx->patchedReadDB));
	int *readCovHist = (int *) malloc(sizeof(int) * (DB_READ_MAXLEN(actx->patchedReadDB) / COV_BIN_SIZE + 1));

	// update cBeg, End intervals for each B-read according global contig positions
	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;

		if (contig->flag & CONTIG_DISCARD)
			continue;

#ifdef DEBUG_STEP1B

		printf("Contig %d, reads: %d\n", contig->idx, contig->numcReads);
#endif
		for (j = 0; j < contig->numcReads; j++)
		{
			ContigRead *cread = contig->cReads + j;

#ifdef DEBUG_STEP1B

			printf("prev: %d cbe %d %d, be %d, %d", cread->patchedID, cread->patchedContigPosBeg, cread->patchedContigPosEnd, cread->ovlReads->beg,cread->ovlReads->end);
#endif

#ifdef DEBUG_STEP1B
			printf(" comp? %d, crlen: %d, nOvlReads %d\n", crComp, DB_READ_LEN(actx->patchedReadDB, abs(cread->patchedID)), cread->numOvlReads);
#endif
			for (k = 1; k < cread->numOvlReads; k++)
			{
#ifdef DEBUG_STEP1B
				printf("\t %d beg,end [%d, %d], cBeg,cEnd [%d, %d]\n", cread->ovlReads[k].patchedID, cread->ovlReads[k].beg, cread->ovlReads[k].end,
						cread->ovlReads[k].cBeg, cread->ovlReads[k].cEnd);
#endif
			}

			{
				int l, m;
				// 1. analyze in-coverage and mark them as visited
				ba_all_assign(inReads, DB_NREADS(actx->patchedReadDB), FALSE);

				for (l = 1; l < cread->numOvlReads; l++)
				{
					if (cread->patchedContigPosBeg == cread->ovlReads[l].cBeg)
					{
						ba_assign(inReads, cread->ovlReads[l].patchedID, TRUE);
						cread->nInReads++;
					}
				}

				// 2. analyze intersection of out-reads from prev reads and in-reads of current read
				if (j) // skip first contig read
				{
					/////
					/////                   r_i    r_j     r_k     r_l     r_m
					/////    contigA:    -------|-------|-------|-------|-------
					/////                ---    |       |       |       |                fixed reads that overlap with each r_i, ... , r_m
					/////                 ---   |       |       |       |
					///                     --- |       |       |       |
					/////                      1 2     1 2     1 2     1 2               1 ... OUT reads (arrival rate) and 2 ... IN reads (start rate)

					int intersect = 0;
					/// loop over IN-reads of current contig read
					for (l = 1; l < cread->numOvlReads; l++)
					{
						if (ba_value(outReads, cread->ovlReads[l].patchedID))
							intersect++;
					}
					contig->cReads[j - 1].nComReadsWithNextContigRead = intersect;
				}

				// 3. analyze
				//  a) out-coverage and mark them as visited
				//  b) average coverage
				//  c) low coverage
				ba_all_assign(outReads, DB_NREADS(actx->patchedReadDB), FALSE);
				bzero(readCovHist, sizeof(int) * (DB_READ_MAXLEN(actx->patchedReadDB) / COV_BIN_SIZE + 1));
				// todo add contig read as well --> i.e. in most cases the number of common exiting/entering reads is not zero
//				ba_assign(outReads, cread->id, TRUE);
				//3a) + 3b)
				for (l = 1; l < cread->numOvlReads; l++)
				{
					if (cread->patchedContigPosEnd == cread->ovlReads[l].cEnd)
					{
						ba_assign(outReads, cread->ovlReads[l].patchedID, TRUE);
						cread->nOutReads++;
					}
					cread->avgCov += cread->ovlReads[l].cEnd - cread->ovlReads[l].cBeg;
				}
				contig->avgCov += cread->avgCov;
				cread->avgCov /= (cread->patchedContigPosEnd - cread->patchedContigPosBeg);

				// 3c) i.e. create coverage histogram  ...
				//printf("Coverage Histogram for read: %d PreOut-CurIn-Cov %d,%d\n", contig->reads[k][0].id, curCheckPoint[curNum-1].nComReadsWithNextRead, curCheckPoint[curNum].nComReadsWithprevRead);
				m = cread->patchedContigPosBeg;
				while (m < cread->patchedContigPosEnd)
				{
					int cFrom, cTo;
					cFrom = m;
					if (m + COV_BIN_SIZE > cread->patchedContigPosEnd)
						cTo = cread->patchedContigPosEnd;
					else
						cTo = cFrom + COV_BIN_SIZE;

					int bin = (cFrom - cread->patchedContigPosBeg) / COV_BIN_SIZE;

					// increase coverage count only if read fully spans histogram bin
					for (l = 1; l < cread->numOvlReads; l++)
					{
						if (cread->ovlReads[l].cBeg <= cFrom && cread->ovlReads[l].cEnd >= cTo)
							readCovHist[bin]++;
					}

					//				printf("<%9d, %9d> c %3d\n", cFrom, cTo, readCovHist[bin]);
					m = cTo;
				}

				// 3c) ... and check coverage drops
				for (l = 0, m = cread->patchedContigPosBeg; m < cread->patchedContigPosEnd; l++, m += COV_BIN_SIZE)
				{
					if (readCovHist[l] < MIN_COV)
					{
						cread->lowCov = 1;
						break;
					}
				}
			}
		}
		// calculate average coverage of contig
		contig->avgCov /= contig->len;

//		if(contig->idx == 0)
		{
			printf("STATS of contig %d\n", contig->idx);
			printf("avgCov: %.2f\n", contig->avgCov);
			for (j = 0; j < contig->numcReads; j++)
			{
				ContigRead *cread = contig->cReads + j;
				printf("%9d - %9d, id %10d, c[%9d, %9d] r[%6d,%6d] inR %10d, outR %d, allR %d, avgCov %.2f, comWNR %d\n", cread->patchedContigPosBeg,
						cread->patchedContigPosEnd, cread->patchedID, cread->ovlReads->cBeg, cread->ovlReads->cEnd, cread->ovlReads->beg, cread->ovlReads->end,
						cread->nInReads, cread->nOutReads, cread->numOvlReads, cread->avgCov, cread->nComReadsWithNextContigRead);
			}
		}
	}

	// cleanup read coverage buffers
	free(readCovHist);
	free(inReads);
	free(outReads);
}

int cmpContigLength(const void *a, const void *b)
{
	return ((const ContigIDAndLen*) b)->len - ((const ContigIDAndLen *) a)->len;
}

int getPathID(AnalyzeContext *actx, int contigID)
{
	if (contigID < 0 || contigID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getPathID: contigID %d out of bounds [%d, %d [\n", contigID, 0, DB_NREADS(actx->corContigDB));
		fflush(stderr);
		exit(1);
	}

	if (Check_Track(actx->corContigDB, "path") != 0)
		return -1;

	HITS_TRACK * pathTrack = track_load(actx->corContigDB, "path");

	assert(pathTrack != NULL);

	track_anno* path_anno = pathTrack->anno;
	track_data* path_data = pathTrack->data;

	track_anno b, e;

	b = path_anno[contigID] / sizeof(track_data);
	e = path_anno[contigID + 1] / sizeof(track_data);

	assert(b < e);

	return (int) path_data[b];
}

void getContigsEndIDs(AnalyzeContext *actx, int contigID, int* beg, int *end)
{
	if (contigID < 0 || contigID >= DB_NREADS(actx->corContigDB))
	{
		fprintf(stderr, "[ERROR] - getContigsEndIDs: contigID %d out of bounds [%d, %d [\n", contigID, 0, DB_NREADS(actx->corContigDB));
		fflush(stderr);
		exit(1);
	}

	if (Check_Track(actx->corContigDB, "ends") != 0)
		return;

	HITS_TRACK * pathTrack = track_load(actx->corContigDB, "ends");

	assert(pathTrack != NULL);

	track_anno* path_anno = pathTrack->anno;
	track_data* path_data = pathTrack->data;

	track_anno b, e;

	b = path_anno[contigID] / sizeof(track_data);
	e = path_anno[contigID + 1] / sizeof(track_data);

	assert(b < e);

	(*beg) = path_data[b];
	(*end) = path_data[b + 1];
}

void classifyContigsByBReadsAndPath(AnalyzeContext *actx)
{
	int i, j, k, l, m;

	ContigIDAndLen *contigIDAndLength = (ContigIDAndLen*) malloc(sizeof(ContigIDAndLen) * actx->numContigs);
	int numCL = 0;

	// todo define a proper struct?
	bit * contigJReads, *contigKReads, *contigIntersectionJKReads;
	int * contigJBegRange, *contigJEndRange;
	int * contigKBegRange, *contigKEndRange;

	int * contigJCovHist, *contigKCovHist;
	int binSize = 100;
	int nBins = (DB_READ_MAXLEN(actx->corContigDB) / binSize) + 1;
	contigJCovHist = (int*) malloc(sizeof(int) * nBins);
	contigKCovHist = (int*) malloc(sizeof(int) * nBins);

	int numReads = DB_NREADS(actx->patchedReadDB);
	contigJReads = ba_new(numReads);
	contigKReads = ba_new(numReads);
	contigIntersectionJKReads = ba_new(numReads);

	contigJBegRange = (int*) malloc(sizeof(int) * numReads);
	contigJEndRange = (int*) malloc(sizeof(int) * numReads);
	contigKBegRange = (int*) malloc(sizeof(int) * numReads);
	contigKEndRange = (int*) malloc(sizeof(int) * numReads);

	i = 1;
	int offset, prevOffset;

	offset = prevOffset = 0;

	int beg, end;

	int faId = 0;
	while (offset < actx->numContigs)
	{
		faId++;
		offset = actx->contigFileNameOffsets[i];

#ifdef DEBUG_STEP1C
		printf("i %d, %d %d\n", i, prevOffset, offset);
#endif
		beg = end = -1;

		getContigsEndIDs(actx, i - 1, &beg, &end);
#ifdef DEBUG_STEP1C
		printf("%s ends=%d,%d:\n", actx->ContigFileNames[i], beg, end);
#endif
		if (offset - prevOffset == 1)
		{
#ifdef DEBUG_STEP1C
			printf("fasta file %d: %s --> has only one contig\n", faId, actx->ContigFileNames[i]);
#endif
			actx->contigs[prevOffset].gClassificFlag |= CONTIG_UNIQUE;
		}
		else // at this point, current compound contains at least 2 contigs
		{
			numCL = 0;

			// add all contigs from current connected compound (same file == subgraph)
			for (j = prevOffset; j < offset; j++)
			{
				Contig *contig = actx->contigs + j;
				if (contig->flag & CONTIG_DISCARD)
					continue;
				contigIDAndLength[numCL].idx = contig->idx;
				contigIDAndLength[numCL].len = contig->len;
#ifdef DEBUG_STEP1C
				printf("%d len: %d\n", contig->idx, contig->len);
#endif
				numCL++;
			}

			// sort contigs according length
			qsort(contigIDAndLength, numCL, sizeof(contigIDAndLength), cmpContigLength);

			// check
			for (j = 0; j < numCL; j++)
			{
				Contig * contig_j = actx->contigs + contigIDAndLength[j].idx;

				if (contig_j->flag & CONTIG_DISCARD)
					continue;

#ifdef DEBUG_STEP1C
				printf("Check against contig_j: %d, len %d\n", contig_j->idx, contig_j->len);
#endif

				if (contig_j->gClassificFlag & CONTIG_IS_CONTAINED)
				{
#ifdef DEBUG_STEP1C
					printf("Skip contig %d, already contained!\n", contig_j->idx);
#endif
//					continue;
				}

				ba_all_assign(contigJReads, numReads, FALSE);
				memset(contigJBegRange, -1, sizeof(int) * numReads);
				memset(contigJEndRange, -1, sizeof(int) * numReads);

				// fill breads and positions for contig_j
				for (k = 0; k < contig_j->numcReads; k++)
				{
					ContigRead *cread = contig_j->cReads + k;

#ifdef DEBUG_STEP1C
					printf("contig_j %d add breads for cread %d\n", contig_j->idx, cread->patchedID);
#endif
					for (l = 0; l < cread->numOvlReads; l++)
					{
						int readID = abs(cread->ovlReads[l].patchedID);

						ba_assign(contigJReads, readID, TRUE);
						if (contigJBegRange[readID] < 0)
						{
							contigJBegRange[readID] = cread->ovlReads[l].cBeg;
							contigJEndRange[readID] = cread->ovlReads[l].cEnd;
						}
						else
						{
							if (contigJBegRange[readID] > cread->ovlReads[l].cBeg)
								contigJBegRange[readID] = cread->ovlReads[l].cBeg;
							if (contigJEndRange[readID] < cread->ovlReads[l].cEnd)
								contigJEndRange[readID] = cread->ovlReads[l].cEnd;
						}
					}
				}

				// compare all other contigs from current subgraph against contig_j
				for (k = j + 1; k < numCL; k++)
				{
					Contig * contig_k = actx->contigs + contigIDAndLength[k].idx;

					if (contig_k->flag & CONTIG_DISCARD)
						continue;

//					if (contig_k->gClassificFlag & CONTIG_IS_CONTAINED)
//					{
#ifdef DEBUG_STEP1C
					printf("Skip contig %d, already contained!\n", contig_k->idx);
#endif
//						continue;
//					}

					ba_all_assign(contigKReads, numReads, FALSE);
					memset(contigKBegRange, -1, sizeof(int) * numReads);
					memset(contigKEndRange, -1, sizeof(int) * numReads);

					for (l = 0; l < contig_k->numcReads; l++)
					{
						ContigRead *cread = contig_k->cReads + l;

#ifdef DEBUG_STEP1C
						printf("contig_k %d add breads for cread %d\n", contig_k->idx, cread->patchedID);
#endif
						for (m = 0; m < cread->numOvlReads; m++)
						{
							int readID = abs(cread->ovlReads[m].patchedID);

							ba_assign(contigKReads, readID, TRUE);
							if (contigKBegRange[readID] < 0)
							{
								contigKBegRange[readID] = cread->ovlReads[m].cBeg;
								contigKEndRange[readID] = cread->ovlReads[m].cEnd;
							}
							else
							{
								if (contigKBegRange[readID] > cread->ovlReads[m].cBeg)
									contigKBegRange[readID] = cread->ovlReads[m].cBeg;
								if (contigKEndRange[readID] < cread->ovlReads[m].cEnd)
									contigKEndRange[readID] = cread->ovlReads[m].cEnd;
							}
						}
					}

					// get read intersection between contig_k and contig_j
					ba_intersection(contigKReads, contigJReads, &contigIntersectionJKReads, numReads, numReads);

					int numComReads = ba_count(contigIntersectionJKReads, numReads);
					int numkReads = ba_count(contigKReads, numReads);

					if (numComReads * 100.0 / numkReads > 1.0)
					{
						memset(contigJCovHist, 0, sizeof(int) * (nBins));
						memset(contigKCovHist, 0, sizeof(int) * (nBins));

						int numCovIvlJ = 0;
						int *covIvlJ = (int*)malloc(sizeof(int)*((numCovIvlJ+1)*3));
						int numCovIvlK = 0;
						int *covIvlK = (int*)malloc(sizeof(int)*((numCovIvlK+1)*3));

						for (l = 0; l < numReads; l++)
						{
							if (ba_value(contigIntersectionJKReads, l))
							{
								int from = contigJBegRange[l] / binSize;
								int to = contigJEndRange[l] / binSize;
								if(from < 0)
								{
									printf("contigIntersectionJKReads Jrange failed for read %d: from %d to %d\n",l, from, to );
								}

								if(to >= nBins)
								{
									printf("contigIntersectionJKReads Jrange failed for read %d: from %d to %d\n",l, from, to );
								}

								while (from <= to)
								{
									if(contig_k->idx == 263)
									{
										printf("from = %d, to = %d Jcount %d nBins %d\n",from, to, contigJCovHist[from],nBins);
										fflush(stdout);
									}
									contigJCovHist[from]++;
									from++;
								}
								from = contigKBegRange[l] / binSize;
								to = contigKEndRange[l] / binSize;

								if(from < 0)
								{
									printf("contigIntersectionJKReads Krange failed for read %d: from %d to %d\n",l, from, to );
								}

								if(to >= nBins)
								{
									printf("contigIntersectionJKReads Krange failed for read %d: from %d to %d\n",l, from, to );
								}
								while (from <= to)
								{
									if(contig_k->idx == 263)
																		{
										printf("from = %d, to = %d Kcount %d nBins %d\n",from, to, contigKCovHist[from],nBins);
																				fflush(stdout);
																		}
									contigKCovHist[from]++;
									from++;
								}
							}
						}
#ifdef DEBUG_STEP1C
						printf("analyze contigs: %d p%d (cov %.2f) %d p%d (cov %.2f) len %d %d --> numOfcommonReads: %d of %d: %.3f\n", contig_j->idx,
								getPathID(actx, contig_j->idx), contig_j->avgCov, contig_k->idx, getPathID(actx, contig_k->idx), contig_k->avgCov, contig_j->len, contig_k->len,
								numComReads, numkReads, numComReads * 100.0 / numkReads);
						printf("contigKCovHist:\n");
#endif
						int firstKbin, lastKbin, firstJbin, lastJbin;
						int cumHistKreads, cumHistJreads;
						cumHistKreads = cumHistJreads = 0;
						firstKbin = lastKbin = firstJbin = lastJbin = -1;
						for (l = 0; l < nBins; l++)
						{
							if (contigKCovHist[l] > 0)
							{
								if (firstKbin == -1)
									firstKbin = l;

								lastKbin = l;
								cumHistKreads += contigKCovHist[l];
#ifdef DEBUG_STEP1C
								printf("%10d-%10d: %4d\n", l * binSize, (l * binSize) + binSize - 1, contigKCovHist[l]);
#endif
							}
							else if(firstKbin != -1)
							{
								covIvlK[numCovIvlK*3] = firstKbin * binSize;
								covIvlK[numCovIvlK*3+1] = lastKbin * binSize + (binSize-1);
								covIvlK[numCovIvlK*3+2] = cumHistKreads / (lastKbin - firstKbin + 1);
								numCovIvlK++;
								covIvlK = (int*)realloc(covIvlK, sizeof(int)*((numCovIvlK+1)*3));
								firstKbin = -1;
								cumHistKreads = 0;
							}
						}
						if(firstKbin != -1)
						{
							covIvlK[numCovIvlK*3] = firstKbin * binSize;
							covIvlK[numCovIvlK*3+1] = lastKbin * binSize + (binSize-1);
							covIvlK[numCovIvlK*3+2] = cumHistKreads / (lastKbin - firstKbin + 1);
							numCovIvlK++;
						}

#ifdef DEBUG_STEP1C
						printf("contigJCovHist:\n");
#endif
						for (l = 0; l < nBins; l++)
						{
							if (contigJCovHist[l] > 0)
							{
								if (firstJbin == -1)
									firstJbin = l;
								cumHistJreads += contigJCovHist[l];
								lastJbin = l;
#ifdef DEBUG_STEP1C
								printf("%10d-%10d: %4d\n", l * binSize, (l * binSize) + binSize - 1, contigJCovHist[l]);
#endif
							}
							else if(firstJbin != -1)
							{
								covIvlJ[numCovIvlJ*3] = firstJbin * binSize;
								covIvlJ[numCovIvlJ*3+1] = lastJbin * binSize + (binSize-1);
								covIvlJ[numCovIvlJ*3+2] = cumHistJreads / (lastJbin - firstJbin + 1);
								numCovIvlJ++;
								covIvlJ = (int*)realloc(covIvlJ, sizeof(int)*((numCovIvlJ+1)*3));
								firstJbin = -1;
								cumHistJreads = 0;
							}
						}
						if(firstJbin != -1)
						{
							covIvlJ[numCovIvlJ*3] = firstJbin * binSize;
							covIvlJ[numCovIvlJ*3+1] = lastJbin * binSize + (binSize-1);
							covIvlJ[numCovIvlJ*3+2] = cumHistJreads / (lastJbin - firstJbin + 1);
							numCovIvlJ++;
						}

#ifdef DEBUG_STEP1C
						printf("common reads for both contigs\n");

						m = 0;
						for (l = 0; l < numReads; l++)
						{
							if (ba_value(contigIntersectionJKReads, l))
							{
								printf("%4d: c %7d vs c %7d read: %7d [%7d, %7d] --> [%7d, %7d]\n", m++, contig_j->idx, contig_k->idx, l, contigJBegRange[l],
										contigJEndRange[l], contigKBegRange[l], contigKEndRange[l]);
							}
						}
#endif



						int cumKGapBases = covIvlK[0] + contig_k->len-covIvlK[3*numCovIvlK-2];
						int cumKCoveredBases = 0;
						for (l=0; l<numCovIvlK; l++)
						{
#ifdef DEBUG_STEP1C
							printf("CoveredBasesInK [%8d, %8d, cov %3d]\n", covIvlK[l*3], covIvlK[l*3+1], covIvlK[l*3+2]);
#endif
							cumKCoveredBases += (covIvlK[l*3+1] - covIvlK[l*3]);

							if(l+1<numCovIvlK)
								cumKGapBases += (covIvlK[(l+1)*3] - covIvlK[l*3+1]);
						}
#ifdef DEBUG_STEP1C
						printf("cumKGapBases %d (%.2f%%) cumKCoveredBases %d (%.2f%%)\n", cumKGapBases, cumKGapBases*100.0/contig_k->len, cumKCoveredBases, cumKCoveredBases*100.0/contig_k->len);
#endif
						int cumJGapBases = 0;
						int cumJCoveredBases = 0;
						for (l=0; l<numCovIvlJ; l++)
						{
#ifdef DEBUG_STEP1C
							printf("CoveredBasesInJ [%8d, %8d, cov %3d]\n", covIvlJ[l*3], covIvlJ[l*3+1], covIvlJ[l*3+2]);
#endif
							cumJCoveredBases += (covIvlJ[l*3+1] - covIvlJ[l*3]);

							if(l+1<numCovIvlJ)
								cumJGapBases += (covIvlJ[(l+1)*3] - covIvlJ[l*3+1]);
						}
#ifdef DEBUG_STEP1C
						printf("cumJGapBases %d (%.2f%%) cumJCoveredBases %d (%.2f%%)\n", cumJGapBases, cumJGapBases*100.0/(covIvlJ[3*numCovIvlJ-2] - covIvlJ[0]), cumJCoveredBases, cumJCoveredBases*100.0/(covIvlJ[3*numCovIvlJ-2] - covIvlJ[0]));
#endif
						// preliminary classification

						int preClassFlagK = 0;
						int preClassFlagJ = 0;

						// contig k is contained in J
						if (numComReads * 100.0 / numkReads > 70.0 || (cumKCoveredBases*100.0/contig_k->len >= 50))
						{
							preClassFlagK |= CONTIG_IS_CONTAINED;
							preClassFlagJ |= CONTIG_HAS_CONTAINED;
						}
						else
						{
							preClassFlagK |= CONTIG_UNCLASSIFIED;
							preClassFlagJ |= CONTIG_UNCLASSIFIED;
						}

#ifdef DEBUG_STEP1C
						printContigClassification(stdout, preClassFlagK, ' ');
						printf("contig %d vs contig %d [%d, %d] [%d,%d]\n", contig_k->idx, contig_j->idx, covIvlK[0], covIvlK[numCovIvlK*3-2], covIvlJ[0], covIvlJ[numCovIvlJ*3-2]);

						printContigClassification(stdout, preClassFlagJ, ' ');
						printf("contig %d vs contig %d [%d, %d] [%d,%d]\n", contig_j->idx, contig_k->idx, covIvlJ[0], covIvlJ[numCovIvlJ*3-2], covIvlK[0], covIvlK[numCovIvlK*3-2]);
#endif

						if (contig_k->numGClassific == contig_k->maxGClassific)
						{
							contig_k->maxGClassific = (int) (contig_k->maxGClassific * 1.2) + 5;
							contig_k->gClassific = (ContigGraphClassification*) realloc(contig_k->gClassific, sizeof(ContigGraphClassification) * contig_k->maxGClassific);
						}

						int num = contig_k->numGClassific;

						contig_k->gClassificFlag = preClassFlagK;

						contig_k->gClassific[num].flag = preClassFlagK;
						contig_k->gClassific[num].numCoveredIntervals = numCovIvlK;
						contig_k->gClassific[num].coveredIntervals  = covIvlK;
						contig_k->gClassific[num].correspID = contig_j->idx;
						contig_k->gClassific[num].nJointReads = numComReads;
						contig_k->numGClassific++;

						// store corresponding info for contig J

						if (contig_j->numGClassific == contig_j->maxGClassific)
						{
							contig_j->maxGClassific = (int) (contig_j->maxGClassific * 1.2) + 5;
							contig_j->gClassific = (ContigGraphClassification*) realloc(contig_j->gClassific, sizeof(ContigGraphClassification) * contig_j->maxGClassific);
						}

						num = contig_j->numGClassific;

						contig_j->gClassificFlag = preClassFlagJ;

						contig_j->gClassific[num].flag = preClassFlagJ;
						contig_j->gClassific[num].numCoveredIntervals = numCovIvlJ;
						contig_j->gClassific[num].coveredIntervals  = covIvlJ;
						contig_j->gClassific[num].correspID = contig_k->idx;
						contig_j->gClassific[num].nJointReads = numComReads;

						contig_j->numGClassific++;
					}
				}
			}
		}
		prevOffset = offset;
		i++;
	}

// cleanup
	free(contigIDAndLength);

	free(contigJReads);
	free(contigKReads);
	free(contigIntersectionJKReads);

	free(contigJBegRange);
	free(contigJEndRange);
	free(contigKEndRange);
	free(contigKBegRange);
}

int compareInt(const void * a, const void * b)
{
	return (*(int*) a - *(int*) b);
}

int contained(int ab, int ae, int bb, int be)
{
	if (ab >= bb && ae <= be)
	{
		return 1;
	}

	return 0;
}

void printOverlapGroup(FILE* out, AnalyzeContext *actx, OverlapGroup *ovlgrp)
{
	int arlen = DB_READ_LEN(actx->corContigDB, ovlgrp->aread);
	int brlen = DB_READ_LEN(actx->corContigDB, ovlgrp->bread);

	fprintf(out, "overlap group %d - %d, l(%d, %d)\n", ovlgrp->aread, ovlgrp->bread, arlen, brlen);
	fprintf(out, "coverBiA: %d (%5.2f)\n", ovlgrp->coveredBasesInA, ovlgrp->coveredBasesInA * 100.0 / arlen);
	fprintf(out, "coverBiB: %d (%5.2f)\n", ovlgrp->coveredBasesInB, ovlgrp->coveredBasesInB * 100.0 / brlen);
	fprintf(out, "gapsA: %d (%5.2f)\n", ovlgrp->gapBasesInA, ovlgrp->gapBasesInA * 100.0 / arlen);
	fprintf(out, "gapsB: %d (%5.2f)\n", ovlgrp->gapBasesInB, ovlgrp->gapBasesInB * 100.0 / brlen);
	fprintf(out, "repeatBiA: %d (%5.2f)\n", ovlgrp->repeatBasesInA, ovlgrp->repeatBasesInA * 100.0 / arlen);
	fprintf(out, "repeatBiB: %d (%5.2f)\n", ovlgrp->repeatBasesInB, ovlgrp->repeatBasesInB * 100.0 / brlen);
	fprintf(out, "uniqA: %d (%5.2f)\n", ovlgrp->coveredBasesInA - ovlgrp->repeatBasesInA, (ovlgrp->coveredBasesInA - ovlgrp->repeatBasesInA) * 100.0 / arlen);
	fprintf(out, "uniqB: %d (%5.2f)\n", ovlgrp->coveredBasesInB - ovlgrp->repeatBasesInB, (ovlgrp->coveredBasesInB - ovlgrp->repeatBasesInB) * 100.0 / brlen);
	fprintf(out, "valid: %d, comp: %c\n", !(ovlgrp->flag & OVLGRP_DISCARD), (ovlgrp->flag & OVLGRP_COMP) ? 'c' : 'n');
	fprintf(out, "first a[%d, %d] b[%d, %d] \n", ovlgrp->first_abpos, ovlgrp->first_aepos, ovlgrp->first_bbpos, ovlgrp->first_bepos);
	fprintf(out, "last a[%d, %d] b[%d, %d] \n", ovlgrp->last_abpos, ovlgrp->last_aepos, ovlgrp->last_bbpos, ovlgrp->last_bepos);
	fprintf(out, "spanning aread: [%d, %d]\n", ovlgrp->first_abpos, ovlgrp->last_aepos);
	if (ovlgrp->flag & OVLGRP_COMP)
		fprintf(out, "spanning bread: [%d, %d]\n", ovlgrp->last_bbpos, ovlgrp->first_bepos);
	else
		fprintf(out, "spanning bread: [%d, %d]\n", ovlgrp->first_bbpos, ovlgrp->last_bepos);
}

int isOverlapGroupValid(AnalyzeContext *actx, OverlapGroup *ovlgrp)
{
	if (ovlgrp->coveredBasesInA - ovlgrp->repeatBasesInA < 100 || ovlgrp->coveredBasesInB - ovlgrp->repeatBasesInB < 100)
	{
#ifdef DEBUG_CONTIGOVERLAP
		printf("isOverlapGroupValid NO: anchor bases in a %d - %d < 100, anchor bases in b %d - %d < 100\n", ovlgrp->coveredBasesInA, ovlgrp->repeatBasesInA,
				ovlgrp->coveredBasesInB, ovlgrp->repeatBasesInB);
#endif
		ovlgrp->flag |= OVLGRP_DISCARD;
		return 0;
	}

	if (ovlgrp->gapBasesInA > ovlgrp->coveredBasesInA || ovlgrp->gapBasesInB > ovlgrp->coveredBasesInB)
	{
#ifdef DEBUG_CONTIGOVERLAP
		printf("isOverlapGroupValid NO: in a %d gaps > %d coverdbases, in b %d gaps > %d coverdbases\n", ovlgrp->gapBasesInA, ovlgrp->coveredBasesInA,
				ovlgrp->gapBasesInB, ovlgrp->coveredBasesInB);
#endif
		ovlgrp->flag |= OVLGRP_DISCARD;
		return 0;
	}

	if (ovlgrp->coveredBasesInA * 100.0 / DB_READ_LEN(actx->corContigDB, ovlgrp->aread) < 10)
	{
#ifdef DEBUG_CONTIGOVERLAP
		printf("isOverlapGroupValid NO: covered based in a smaller than 10%% of A's length\n");
#endif
		ovlgrp->flag |= OVLGRP_DISCARD;
		return 0;
	}

	return 1;
}

OverlapGroup* copySymmetricOverlapGroup(OverlapGroup *ovlgrp)
{
	OverlapGroup * symOvlgrp = (OverlapGroup*) malloc(sizeof(OverlapGroup));
	bzero(symOvlgrp, sizeof(OverlapGroup));

	symOvlgrp->aread = ovlgrp->bread;
	symOvlgrp->bread = ovlgrp->aread;

	symOvlgrp->coveredBasesInA = ovlgrp->coveredBasesInB;
	symOvlgrp->coveredBasesInB = ovlgrp->coveredBasesInA;

	if (ovlgrp->flag & OVLGRP_COMP)
		symOvlgrp->flag |= OVLGRP_COMP;
	if (ovlgrp->flag & OVLGRP_AREAD_HAS_CONTAINED)
		symOvlgrp->flag |= OVLGRP_AREAD_IS_CONTAINED;
	if (ovlgrp->flag & OVLGRP_AREAD_IS_CONTAINED)
		symOvlgrp->flag |= OVLGRP_AREAD_HAS_CONTAINED;
	if (ovlgrp->flag & OVLGRP_DISCARD)
		symOvlgrp->flag |= OVLGRP_DISCARD;
	if (ovlgrp->flag & OVLGRP_TOO_SHORT)
		symOvlgrp->flag |= OVLGRP_TOO_SHORT;

	symOvlgrp->gapBasesInA = ovlgrp->gapBasesInB;
	symOvlgrp->gapBasesInB = ovlgrp->gapBasesInA;

	symOvlgrp->repeatBasesInA = ovlgrp->repeatBasesInB;
	symOvlgrp->repeatBasesInB = ovlgrp->repeatBasesInA;

	if (ovlgrp->flag & OVLGRP_COMP)
	{
		symOvlgrp->first_abpos = ovlgrp->last_bbpos;
		symOvlgrp->first_aepos = ovlgrp->last_bepos;

		symOvlgrp->last_abpos = ovlgrp->first_bbpos;
		symOvlgrp->last_aepos = ovlgrp->first_bepos;

		symOvlgrp->first_bbpos = ovlgrp->last_abpos;
		symOvlgrp->first_bepos = ovlgrp->last_aepos;

		symOvlgrp->last_bbpos = ovlgrp->first_abpos;
		symOvlgrp->last_bepos = ovlgrp->first_aepos;
	}
	else
	{
		symOvlgrp->first_abpos = ovlgrp->first_bbpos;
		symOvlgrp->first_aepos = ovlgrp->first_bepos;

		symOvlgrp->first_bbpos = ovlgrp->first_abpos;
		symOvlgrp->first_bepos = ovlgrp->first_aepos;

		symOvlgrp->last_abpos = ovlgrp->last_bbpos;
		symOvlgrp->last_aepos = ovlgrp->last_bepos;

		symOvlgrp->last_bbpos = ovlgrp->last_abpos;
		symOvlgrp->last_bepos = ovlgrp->last_aepos;
	}

	return symOvlgrp;
}
/*
int findOverlapGroups(AnalyzeContext* actx, Overlap* pOvls, int n)
{
	static int maxOverlapIdx = 0;
	static int* overlapIdx = NULL;
	int numOverlapIdx = 0;

	int aread, bread;
	int i;

	aread = pOvls->aread;
	bread = pOvls->bread;

// mark contained overlaps
	{
		int j;
		for (i = 0; i < n; i++)
		{
			Overlap *ovl_i = pOvls + i;

			if (ovl_i->flags & OVL_CONT)
				continue;

			for (j = i + 1; j < n; j++)
			{
				Overlap *ovl_j = pOvls + j;

				if (contained(ovl_j->path.abpos, ovl_j->path.aepos, ovl_i->path.abpos, ovl_i->path.aepos)
						&& contained(ovl_j->path.bbpos, ovl_j->path.bepos, ovl_i->path.bbpos, ovl_i->path.bepos))
				{
					ovl_j->flags |= OVL_CONT;

				}
			}

		}
	}

// get repeats track
	track_anno* rep_anno = actx->corContigRepeats_track->anno;
	track_data* rep_data = actx->corContigRepeats_track->data;

	track_anno rab, rae;
	track_anno rbb, rbe;

	int longestUniqBases = -1;
	int longestUniqOVLidx = -1;
	int longestOverlap = -1;
	int longestOverlapIdx = -1;

	// find longest overlap based on number of unique bases
	for (i = 0; i < n; i++)
	{
		Overlap *ovl = pOvls + i;

		if ((ovl->flags & OVL_CONT) || (ovl->flags & OVL_DISCARD) || (actx->contigs[ovl->bread].flag & CONTIG_DISCARD))
		{
			continue;
		}

		int aLen = ovl->path.aepos - ovl->path.abpos;
		int bLen = ovl->path.bepos - ovl->path.bbpos;

		int aRep = 0;
		int bRep = 0;

		rab = rep_anno[aread] / sizeof(track_data);
		rae = rep_anno[aread + 1] / sizeof(track_data);

		// loop through all repeats in a
		int rBeg, rEnd;
		while (rab < rae)
		{
			rBeg = rep_data[rab];
			rEnd = rep_data[rab + 1];

			aRep += intersect(ovl->path.abpos, ovl->path.aepos, rBeg, rEnd);
			rab += 2;
		}

		rbb = rep_anno[bread] / sizeof(track_data);
		rbe = rep_anno[bread + 1] / sizeof(track_data);

		// loop through all repeats in b
		while (rbb < rbe)
		{
			rBeg = rep_data[rbb];
			rEnd = rep_data[rbb + 1];

			if (ovl->flags & OVL_COMP)
			{
				bRep += intersect(DB_READ_LEN(actx->corContigDB, bread) - ovl->path.bepos, DB_READ_LEN(actx->corContigDB, bread) - ovl->path.bbpos, rBeg, rEnd);
			}
			else
			{
				bRep += intersect(ovl->path.bbpos, ovl->path.bepos, rBeg, rEnd);
			}
			rbb += 2;
		}
#ifdef DEBUG_STEP2A
//        printf("%d - %d [%d, %d] [%d, %d], aR %d/%d, bR %d/%d\n", aread, bread, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, aLen, aRep, bLen, bRep);
#endif
		int tmpBases = MAX(aLen - aRep, bLen - bRep);
		if (tmpBases > longestUniqBases)
		{
			longestUniqBases = tmpBases;
			longestUniqOVLidx = i;
		}

		tmpBases = MAX(aLen, bLen);
		if (tmpBases > longestOverlap)
		{
			longestOverlap = tmpBases;
			longestOverlapIdx = i;
		}
	}

// TODO min unique bases and min repeat bases hard coded
	if (longestUniqBases < 1000)
	{
		if (longestOverlap > 5000)
		{
			longestUniqBases = longestOverlap;
			longestUniqOVLidx = longestOverlapIdx;
		}
		else
		{
			return 0;
		}
	}
#ifdef DEBUG_STEP2A
	printf("longest overlap:\n");
	printf("idx: %d --> uB %d, %d - %d [%d, %d] [%d, %d]\n", longestUniqOVLidx, longestUniqBases, pOvls[longestUniqOVLidx].aread, pOvls[longestUniqOVLidx].bread,
			pOvls[longestUniqOVLidx].path.abpos, pOvls[longestUniqOVLidx].path.aepos, pOvls[longestUniqOVLidx].path.bbpos, pOvls[longestUniqOVLidx].path.bepos);
#endif
// try to "elongate" longest overlap
// 1st on the right
// 2nd on the left side
	OverlapGroup * ovlgrp = (OverlapGroup*) malloc(sizeof(OverlapGroup));
	bzero(ovlgrp, sizeof(OverlapGroup));

	ovlgrp->aread = pOvls->aread;
	ovlgrp->bread = pOvls->bread;

	if (numOverlapIdx == maxOverlapIdx)
	{
		maxOverlapIdx = maxOverlapIdx * 1.2 + 10;
		overlapIdx = (int *) realloc(overlapIdx, sizeof(int) * maxOverlapIdx);
	}

	overlapIdx[0] = longestUniqOVLidx;
	numOverlapIdx = 1;

	int ab1, ae1;
	int bb1, be1;

	int ab2, ae2;
	int bb2, be2;

	ab1 = pOvls[longestUniqOVLidx].path.abpos;
	ae1 = pOvls[longestUniqOVLidx].path.aepos;

	bb1 = pOvls[longestUniqOVLidx].path.bbpos;
	be1 = pOvls[longestUniqOVLidx].path.bepos;

#ifdef DEBUG_STEP2A
	printf("extend longest overlap in right direction\n");
#endif
// 1st right
	int cont = 1;
	int curBestOffset = 1;
	int curBestUniqBases = -1;
	int aReadLen = DB_READ_LEN(actx->corContigDB, pOvls->aread);
	int bReadLen = DB_READ_LEN(actx->corContigDB, pOvls->bread);
	while (cont && ae1 < aReadLen)
	{
		int stepSize;
		for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestUniqBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
		{
#ifdef DEBUG_STEP2A
					printf("FOR LOOP stepsize %d\n", stepSize);
#endif
			for (i = longestUniqOVLidx + curBestOffset; i < n; i++)
			{
				Overlap * ovl = pOvls + i;

				if (ovl->flags & OVL_DISCARD)
					continue;

				ab2 = pOvls[i].path.abpos;
				ae2 = pOvls[i].path.aepos;

				bb2 = pOvls[i].path.bbpos;
				be2 = pOvls[i].path.bepos;

				if ((ovl->flags & OVL_COMP) != (pOvls[longestUniqOVLidx].flags & OVL_COMP))
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> different orientations\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if (ovl->flags & OVL_CONT)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> really contained repeat\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if (ae1 >= ae2)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in A-interval\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if (ae2 - ab2 < actx->min_olen || be2 - bb2 < actx->min_olen)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> to short\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				// add some sanity checks
				// 1. overlap must go into right direction
				// 2. if overlaps are contained replace them
				if (ae2 - 100 <= ae1 || be2 - 100 <= be1) // at least 100-bases overhang to the right
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> no overhang on right side\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if (ab2 - ae1 > stepSize)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> gap size too large (stepSize %d)\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2,
							stepSize);
#endif
					continue;
				}

				// check if current overlap is better (longer/more unique bases ?) then curBest

				int curUniqBasesInAIvl = (ae2 - ab2) - getRepeatBasesFromInterval(actx, ovl->aread, ab2, ae2);
				int curUniqBasesInBIvl;

				if (ovl->flags & OVL_COMP)
				{
					curUniqBasesInBIvl = (be2 - bb2) - getRepeatBasesFromInterval(actx, ovl->bread, bReadLen - be2, bReadLen - bb2);
				}
				else
				{
					curUniqBasesInBIvl = (be2 - bb2) - getRepeatBasesFromInterval(actx, ovl->bread, bb2, be2);
				}

				if (curBestUniqBases < curUniqBasesInAIvl + curUniqBasesInBIvl && curUniqBasesInAIvl > 0 && curUniqBasesInBIvl > 0)
				{
#ifdef DEBUG_STEP2A
					printf("found right current best overlap %d %d %d: [%d, %d] [%d, %d] right side\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					curBestOffset = i - longestUniqOVLidx;
					curBestUniqBases = curUniqBasesInAIvl + curUniqBasesInBIvl;
				}
				else
				{
#ifdef DEBUG_STEP2A
					printf("NO right ovl found skip: %d vs %d [%d, %d] [%d, %d] left side curUniqBasesInAIvl %d curUniqBasesInBIvl %d\n", pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2, curUniqBasesInAIvl, curUniqBasesInBIvl);
#endif
				}
			}
		}

		// check if left best overlap can be used to extend overlap group on the right side
		if (curBestUniqBases < 0) // i.e. there was no good overlap at right side
		{
#ifdef DEBUG_STEP2A
			printf("could not extend ovlgroup on right side with proper overlap (with stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
#endif
			break;
		}

		/// todo further sanity check necessary ???
		ab2 = pOvls[longestUniqOVLidx + curBestOffset].path.abpos;
		ae2 = pOvls[longestUniqOVLidx + curBestOffset].path.aepos;

		bb2 = pOvls[longestUniqOVLidx + curBestOffset].path.bbpos;
		be2 = pOvls[longestUniqOVLidx + curBestOffset].path.bepos;
#ifdef DEBUG_STEP2A
		printf("extend ovlgroup with (right): %d %d %d: [%d, %d] [%d, %d] stepSize %d\n", longestUniqOVLidx + curBestOffset,
				pOvls[longestUniqOVLidx + curBestOffset].aread, pOvls[longestUniqOVLidx + curBestOffset].bread, ab1, ae1, ab2, ae2, stepSize - STRIDE_OVL_LOOKAHEAD);
#endif
		if (numOverlapIdx == maxOverlapIdx)
		{
			maxOverlapIdx = maxOverlapIdx * 1.2 + 10;
			overlapIdx = (int *) realloc(overlapIdx, sizeof(int) * maxOverlapIdx);
		}

		// append left side overlaps at the end of overlap group, i.e. overlap group must be sorted afterwards by overlap index
		overlapIdx[numOverlapIdx] = longestUniqOVLidx + curBestOffset;
		numOverlapIdx++;

		ab1 = ab2;
		ae1 = ae2;
		bb1 = bb2;
		be1 = be2;

		curBestOffset++;
		curBestUniqBases = -1;
		if (longestUniqOVLidx + curBestOffset >= n)
		{
			cont = 0;
		}
	}

	ab1 = pOvls[longestUniqOVLidx].path.abpos;
	ae1 = pOvls[longestUniqOVLidx].path.aepos;

	bb1 = pOvls[longestUniqOVLidx].path.bbpos;
	be1 = pOvls[longestUniqOVLidx].path.bepos;

#ifdef DEBUG_STEP2A
	printf("extend longest overlap in left direction\n");
#endif
// 2nd left side
	cont = 1;
	curBestOffset = 1;
	curBestUniqBases = -1;
	while (cont && ab1 > 0)
	{
		int stepSize;
		for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestUniqBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
		{
#ifdef DEBUG_STEP2A
					printf("FOR LOOP stepsize %d\n", stepSize);
#endif

			// try to find next best overlap with lookahead of stepSize bases
			for (i = longestUniqOVLidx - curBestOffset; i >= 0; --i)
			{
				Overlap * ovl = pOvls + i;

				if (ovl->flags & OVL_DISCARD)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> DISCARD flag set! flag %d\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2, pOvls[i].flags);
#endif
					continue;
				}

				ab2 = pOvls[i].path.abpos;
				ae2 = pOvls[i].path.aepos;

				bb2 = pOvls[i].path.bbpos;
				be2 = pOvls[i].path.bepos;

				if ((ovl->flags & OVL_COMP) != (pOvls[longestUniqOVLidx].flags & OVL_COMP))
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> different orientations\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if (ovl->flags & OVL_CONT)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> really contained repeat\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if ((ab2 >= ab1 && ae2 <= ae1) || (bb2 >= bb1 && be2 <= be1))
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in A-interval\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if (ae2 - ab2 < actx->min_olen || be2 - bb2 < actx->min_olen)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> to short\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				// add some sanity checks
				// 1. overlap must go into left direction
				// 2. if overlaps are contained replace them
				if ((ab2 + 100 >= ab1 || bb2 + 100 >= bb1))
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> no overhang on left side\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					continue;
				}

				if (ab1 - ae2 > stepSize)
				{
#ifdef DEBUG_STEP2A
					printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> gap size too large (stepSize %d)\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2,
							stepSize);
#endif
					continue;
				}

				// check if current overlap is better (longer/more unique bases ?) then curLeftBest

				int curUniqBasesInAIvl = (ae2 - ab2) - getRepeatBasesFromInterval(actx, ovl->aread, ab2, ae2);
				int curUniqBasesInBIvl;

				if (ovl->flags & OVL_COMP)
				{
					curUniqBasesInBIvl = (be2 - bb2) - getRepeatBasesFromInterval(actx, ovl->bread, bReadLen - be2, bReadLen - bb2);
				}
				else
				{
					curUniqBasesInBIvl = (be2 - bb2) - getRepeatBasesFromInterval(actx, ovl->bread, bb2, be2);
				}

				if (curBestUniqBases < curUniqBasesInAIvl + curUniqBasesInBIvl && curUniqBasesInAIvl > 0 && curUniqBasesInBIvl > 0)
				{
#ifdef DEBUG_STEP2A
					printf("found left current best overlap %d %d %d: [%d, %d] [%d, %d] left side\n", i, pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2);
#endif
					curBestOffset = longestUniqOVLidx - i;
					curBestUniqBases = curUniqBasesInAIvl + curUniqBasesInBIvl;
				}
				else
				{
#ifdef DEBUG_STEP2A
					printf("NO left ovl found skip: %d vs %d [%d, %d] [%d, %d] left side curUniqBasesInAIvl %d curUniqBasesInBIvl %d\n", pOvls[i].aread, pOvls[i].bread, ab2, ae2, bb2, be2, curUniqBasesInAIvl, curUniqBasesInBIvl);
#endif
				}
			}
		}

		// check if left best overlap can be used to extend overlap group on the left side
		if (curBestUniqBases < 0) // i.e. there was no good overlap at left side
		{
#ifdef DEBUG_STEP2A
			printf("could not extend ovlgroup on left side with proper overlap (stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
#endif
			break;
		}

		/// todo further sanity check necessary ???
		ab2 = pOvls[longestUniqOVLidx - curBestOffset].path.abpos;
		ae2 = pOvls[longestUniqOVLidx - curBestOffset].path.aepos;

		bb2 = pOvls[longestUniqOVLidx - curBestOffset].path.bbpos;
		be2 = pOvls[longestUniqOVLidx - curBestOffset].path.bepos;

#ifdef DEBUG_STEP2A
		printf("extend ovlgroup with (left): %d %d %d: [%d, %d] [%d, %d] with stepSize %d\n", longestUniqOVLidx - curBestOffset,
				pOvls[longestUniqOVLidx - curBestOffset].aread, pOvls[longestUniqOVLidx - curBestOffset].bread, ab1, ae1, ab2, ae2, stepSize - STRIDE_OVL_LOOKAHEAD);
#endif

		if (numOverlapIdx == maxOverlapIdx)
		{
			maxOverlapIdx = maxOverlapIdx * 1.2 + 10;
			overlapIdx = (int *) realloc(overlapIdx, sizeof(int) * maxOverlapIdx);
		}

		// append left side overlaps at the end of overlap group, i.e. overlap group must be sorted afterwards by overlap index
		overlapIdx[numOverlapIdx] = longestUniqOVLidx - curBestOffset;
		numOverlapIdx++;

		ab1 = ab2;
		ae1 = ae2;
		bb1 = bb2;
		be1 = be2;

		curBestOffset++;
		curBestUniqBases = -1;
		if (longestUniqOVLidx - curBestOffset < 0)
		{
			cont = 0;
		}
	}

// sort overlaps according overlap indexes
	qsort(overlapIdx, numOverlapIdx, sizeof(int), compareInt);

// update number of: covered bases, gap bases, unique bases,
	updateOverlapGroup(actx, ovlgrp, pOvls, overlapIdx, numOverlapIdx);

#ifdef DEBUG_STEP2A
	printf("final overlap group:\n");
	printOverlapGroup(stdout, actx, ovlgrp);
	printf("%% of a: %4.1f\n", ovlgrp->coveredBasesInA * 100.0 / DB_READ_LEN(actx->corContigDB, aread));
	printf("%% of b: %4.1f\n", ovlgrp->coveredBasesInB * 100.0 / DB_READ_LEN(actx->corContigDB, bread));

	printf("%% gaps in a: %4.1f\n", ovlgrp->gapBasesInA * 100.0 / DB_READ_LEN(actx->corContigDB, aread));
	printf("%% gaps in b: %4.1f\n", ovlgrp->gapBasesInB * 100.0 / DB_READ_LEN(actx->corContigDB, bread));

	printf("%% spanning a: %4.1f\n", (ovlgrp->coveredBasesInA + ovlgrp->gapBasesInA) * 100.0 / DB_READ_LEN(actx->corContigDB, aread));
	printf("%% spanning b: %4.1f\n", (ovlgrp->coveredBasesInB + ovlgrp->gapBasesInB) * 100.0 / DB_READ_LEN(actx->corContigDB, bread));
#endif
	if (isOverlapGroupValid(actx, ovlgrp))
	{
#ifdef DEBUG_STEP2A
		printOverlapGroup(stdout, actx, ovlgrp);
#endif
//		addOverlapGroupToContig(actx, ovlgrp, 1);
	}
	else
		free(ovlgrp);
	return 1;
}
*/

static int getRepeatBasesOfContigRange(AnalyzeContext *ctx, Overlap *ovl, int read)
{
	if (ctx->corContigRepeats_track == NULL)
	{
		return 0;
	}

	assert(ovl->aread == read || ovl->bread == read);

	int bLen = ovl->path.bepos - ovl->path.bbpos;

	// get repeats track
	track_anno* rep_anno = ctx->corContigRepeats_track->anno;
	track_data* rep_data = ctx->corContigRepeats_track->data;

	track_anno rb, re;
	int nrep = 0;

	rb = rep_anno[read] / sizeof(track_data);
	re = rep_anno[read + 1] / sizeof(track_data);

	// loop through all repeats in a
	int rBeg, rEnd;
	while (rb < re)
	{
		rBeg = rep_data[rb];
		rEnd = rep_data[rb + 1];

		if (ovl->aread == read)
		{
			nrep += intersect(ovl->path.abpos, ovl->path.aepos, rBeg, rEnd);
			if(rEnd > ovl->path.aepos)
				break;
		}
		else
		{
			if (ovl->flags & OVL_COMP)
			{
				nrep += intersect(bLen - ovl->path.bepos, bLen - ovl->path.bbpos, rBeg, rEnd);
				if(rEnd > bLen - ovl->path.bbpos)
					break;
			}
			else
			{
				nrep += intersect(ovl->path.bbpos, ovl->path.bepos, rBeg, rEnd);
				if(rEnd > ovl->path.bepos)
					break;
			}
		}
		rb += 2;
	}
	return nrep;
}

static int cmp_ovls_abeg(const void* a, const void* b)
{
	Overlap* o1 = *(Overlap**) a;
	Overlap* o2 = *(Overlap**) b;

	int cmp = o1->path.abpos - o2->path.abpos;

	if (!cmp)
	{
		cmp = (o1->path.aepos - o1->path.abpos) - (o2->path.aepos - o2->path.abpos);
	}

	return cmp;
}

static int cmp_chain_len(const void *a, const void *b)
{
	Chain* c1 = (Chain *) a;
	Chain* c2 = (Chain *) b;

	int i;

	int olen1 = c1->ovls[0]->path.aepos - c1->ovls[0]->path.abpos;

	for (i = 1; i < c1->novl; i++)
	{
		olen1 += c1->ovls[i]->path.aepos - c1->ovls[i]->path.abpos;

		if (c1->ovls[i - 1]->path.aepos > c1->ovls[i]->path.abpos)
			olen1 -= c1->ovls[i - 1]->path.aepos > c1->ovls[i]->path.abpos;
	}

	int olen2 = c2->ovls[0]->path.aepos - c2->ovls[0]->path.abpos;

	for (i = 1; i < c2->novl; i++)
	{
		olen2 += c2->ovls[i]->path.aepos - c2->ovls[i]->path.abpos;

		if (c2->ovls[i - 1]->path.aepos > c2->ovls[i]->path.abpos)
			olen2 -= c2->ovls[i - 1]->path.aepos > c2->ovls[i]->path.abpos;
	}

	return (olen2 - olen1);
}

void chainContigOverlaps(AnalyzeContext* ctx, Overlap* ovls, int n)
{
	int conAId = ovls->aread;
	int conBId = ovls->bread;

	assert(conAId != conBId);

	Contig *conA = ctx->contigs + conAId;
	Contig *conB = ctx->contigs + conBId;

	if (n < 2)
	{
		// add a single overlap into chain if its somehow long enough
		if((ovls->path.aepos - ovls->path.abpos) >= 0.5*conA->len || (ovls->path.bepos - ovls->path.bbpos) >= 0.5*conB->len ||
		  (	((ovls->path.abpos = 0) || (ovls->path.aepos = conA->len)) && ((ovls->path.bbpos = 0) || (ovls->path.bepos = conB->len)) && (ovls->path.aepos - ovls->path.abpos > 15000)) // putative join
			)
		{
			if (ctx->curChains == ctx->maxChains)
			{
				ctx->maxChains = ctx->maxChains * 1.2 + 5;
				ctx->ovlChains = (Chain*) realloc(ctx->ovlChains, sizeof(Chain) * ctx->maxChains);
				bzero(ctx->ovlChains + ctx->curChains, sizeof(Chain) * (ctx->maxChains - ctx->curChains));
			}

			Chain *curChain = ctx->ovlChains + ctx->curChains;

			if(curChain->novl == curChain->maxOvl)
			{
				curChain->maxOvl = curChain->novl * 1.2 + n;
				curChain->ovls = (Overlap**) realloc(curChain->ovls, sizeof(Overlap *)*curChain->maxOvl);
				bzero(curChain->ovls + curChain->novl, sizeof(Overlap*) * (curChain->maxOvl - curChain->novl));
			}
			curChain->ovls[0] = ovls;
			curChain->novl = 1;
			ovls->flags |= OVL_TEMP;
		}
		else
		{
			ovls->flags |= OVL_DISCARD;
		}
	}
	else
	{
		int i;

		// try to detect all possible chains
		int nremain = n;

		// mark contained overlaps
//		#ifdef DEBUG_CHAIN
		if(ovls->aread == 4 && ovls->bread == 839)
			printf("mark contained overlaps\n");
//		#endif
		{
			int j;
			for (i = 0; i < n; i++)
			{
				Overlap *ovl_i = ovls + i;

				if (ovl_i->flags & (OVL_CONT))
					continue;

				for (j = i + 1; j < n; j++)
				{
					Overlap *ovl_j = ovls + j;

					if (ovl_j->flags & (OVL_CONT))
						continue;

					if (contained(ovl_j->path.abpos, ovl_j->path.aepos, ovl_i->path.abpos, ovl_i->path.aepos)
							&& contained(ovl_j->path.bbpos, ovl_j->path.bepos, ovl_i->path.bbpos, ovl_i->path.bepos))
					{
						ovl_j->flags |= (OVL_CONT);
					}
				}

			}
		}

		//		#ifdef DEBUG_CHAIN
		if(ovls->aread == 4 && ovls->bread == 839)
			printf("check initial overlaps flags\n");
//		#endif
		{
			int j;
			for (i = 0; i < n; i++)
			{
				Overlap *ovl_i = ovls + i;

				if (ovl_i->flags & (OVL_CONT))
				{
					ovl_i->flags |= OVL_DISCARD;
					nremain--;
				}
			}
		}

//		#ifdef DEBUG_CHAIN
		if(ovls->aread == 4 && ovls->bread == 839)
			printf("nremain %d\n", nremain);
//		#endif
		assert(nremain >= 1);

		int whileIdx=0;
		while (nremain > 0)
		{
			int longestUniqOvlBases = -1;
			int longestUniqOvlIdx = -1;
			int longestOvlBases = -1;
			int longestOvlIdx = -1;

			// find longest overlap based on number of unique bases
			for (i = 0; i < n; i++)
			{
				Overlap *ovl = ovls + i;

				if (ovl->flags & (OVL_CONT | OVL_DISCARD | OVL_TEMP))
				{
					continue;
				}

				int aLen = ovl->path.aepos - ovl->path.abpos;
				int bLen = ovl->path.bepos - ovl->path.bbpos;

				int aRep = getRepeatBasesOfContigRange(ctx, ovl, ovl->aread);
				int bRep = getRepeatBasesOfContigRange(ctx, ovl, ovl->bread);

//#ifdef DEBUG_CHAIN
				if(ovls->aread == 4 && ovls->bread == 839)
			printf("%d - %d [%d, %d] [%d, %d], aR %d/%d, bR %d/%d\n", conAId, conBId, ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, aLen, aRep,
					bLen, bRep);
//#endif
				int tmpBases = MAX(aLen - aRep, bLen - bRep);
				if (tmpBases > longestUniqOvlBases)
				{
					longestUniqOvlBases = tmpBases;
					longestUniqOvlIdx = i;
				}

				tmpBases = MAX(aLen, bLen);
				if (tmpBases > longestOvlBases)
				{
					longestOvlBases = tmpBases;
					longestOvlIdx = i;
				}
			}

			if (longestUniqOvlBases < ctx->twidth && longestOvlBases > longestUniqOvlBases)
			{
//#ifdef DEBUG_CHAIN
				if(ovls->aread == 4 && ovls->bread == 839)
			printf("Number of unique bases to low. Use longest overlap.\n");
//#endif
				longestUniqOvlBases = longestOvlBases;
				longestUniqOvlIdx = longestOvlIdx;
			}

			// break out
			if(ctx->curChains && (longestOvlBases < 5000 && !((longestOvlBases*1.0/conB->len > 0.5)||(longestOvlBases*1.0/conA->len > 0.5))) )
			{
				for (i = 0; i < n; i++)
				{
					Overlap *ovl = ovls + i;

					if (ovl->flags & (OVL_CONT | OVL_DISCARD | OVL_TEMP))
					{
						continue;
					}

					ovl->flags |= (OVL_DISCARD | OVL_OLEN);
				}
				nremain = 0;
				break;
			}

//#ifdef DEBUG_CHAIN
			if(ovls->aread == 4 && ovls->bread == 839)
			{
		printf("longest overlap:\n");
		printf("idx: %d --> uB %d, %d - %d [%d, %d] [%d, %d]\n", longestUniqOvlIdx, longestUniqOvlBases, ovls[longestUniqOvlIdx].aread,
				ovls[longestUniqOvlIdx].bread, ovls[longestUniqOvlIdx].path.abpos, ovls[longestUniqOvlIdx].path.aepos, ovls[longestUniqOvlIdx].path.bbpos,
				ovls[longestUniqOvlIdx].path.bepos);
			}
//#endif

			// try to "elongate" longest overlap
			// 1st on the right
			// 2nd on the left side

			if (ctx->curChains == ctx->maxChains)
			{
				ctx->maxChains = ctx->maxChains * 1.2 + 5;
				ctx->ovlChains = (Chain*) realloc(ctx->ovlChains, sizeof(Chain) * ctx->maxChains);

				bzero(ctx->ovlChains + ctx->curChains, sizeof(Chain) * (ctx->maxChains - ctx->curChains));
			}

			Chain *curChain = ctx->ovlChains + ctx->curChains;

//#ifdef DEBUG_CHAIN
			if(ovls->aread == 4 && ovls->bread == 839)
		printf("chain: nOvl: %d, maxOvl %d, nremain: %d\n", curChain->novl, curChain->maxOvl, nremain);
//#endif
			if(curChain->novl == curChain->maxOvl)
			{
				curChain->maxOvl = curChain->novl * 1.2 + n;
				curChain->ovls = (Overlap**) realloc(curChain->ovls, sizeof(Overlap *)*curChain->maxOvl);
				bzero(curChain->ovls + curChain->novl, sizeof(Overlap*) * (curChain->maxOvl - curChain->novl));
//#ifdef DEBUG_CHAIN
				if(ovls->aread == 4 && ovls->bread == 839)
			printf("realloc > chain: nOvl: %d, maxOvl %d\n", curChain->novl, curChain->maxOvl);
//#endif
			}

			if (longestUniqOvlIdx < 0 || longestUniqOvlIdx >=n)
			{
				printf("TERROR: longest index out of bounds!!! longestUniqOvlIdx %d, contigs %d vs %d\n", longestUniqOvlIdx, ovls->aread, ovls->bread);
				fflush(stdout);
			}

			assert(longestUniqOvlIdx >=0 && longestUniqOvlIdx < n);
			// first: add longest overlap to chain
			curChain->ovls[0] = ovls + longestUniqOvlIdx;
			curChain->ovls[0]->flags |= OVL_TEMP;
			curChain->novl++;
			nremain--;
//	#ifdef DEBUG_CHAIN
			if(ovls->aread == 4 && ovls->bread == 839)
			printf("chain: nOvl: %d, maxOvl %d, nremain %d\n", curChain->novl, curChain->maxOvl, nremain);
//	#endif

			int ab1, ae1;
			int bb1, be1;

			int ab2, ae2;
			int bb2, be2;

			if (nremain && longestUniqOvlIdx + 1 < n)
			{
				ab1 = ovls[longestUniqOvlIdx].path.abpos;
				ae1 = ovls[longestUniqOvlIdx].path.aepos;

				bb1 = ovls[longestUniqOvlIdx].path.bbpos;
				be1 = ovls[longestUniqOvlIdx].path.bepos;

	#ifdef DEBUG_CHAIN
				printf("extend longest overlap in right direction\n");
	#endif
				// extend longest overlap into right direction
				int cont = 1;
				int curBestUniqOffset = 1;
				int curBestUniqBases = -1;
				int curBestBases = -1;
				int curBestOffset = 1;
				int curBestIntersection = MAX(conA->len, conB->len);

				while (cont)
				{
					int stepSize;
					for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestUniqBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
					{
	#ifdef DEBUG_CHAIN
						printf("FOR LOOP stepsize %d\n", stepSize);
	#endif

						for (i = longestUniqOvlIdx + curBestUniqOffset; i < n; i++)
						{
							Overlap * ovl = ovls + i;

							ab2 = ovl->path.abpos;
							ae2 = ovl->path.aepos;

							bb2 = ovl->path.bbpos;
							be2 = ovl->path.bepos;

							if ((ovl->flags & OVL_COMP) != (ovls[longestUniqOvlIdx].flags & OVL_COMP))
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> different orientations\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ovl->flags & OVL_CONT)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> really contained repeat\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ovl->flags & OVL_TEMP)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> is part of another chain\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ovl->flags & OVL_DISCARD)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> discarded\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							// todo mark those as discard ????
							if (contained(ab2, ae2, ab1, ae1))
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in A-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							// todo mark those as discard ????
							if (contained(bb2, be2, bb1, be1))
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in B-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ae2 < ae1 || be2 < be1) // also overlap must follow right direction
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> improper right extension direction (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2,
										bb2, be2, stepSize);
	#endif
								continue;
							}

							if (MAX(ab2 - ae1, bb2 - be1) > stepSize)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> gap size too large (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
										stepSize);
	#endif
								continue;
							}

							if (ae1 - ab2 > ae2 - ae1 || be1 - bb2 > be2 - be1) // at least 50% overhang
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> overhang too short (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
										stepSize);
	#endif
								continue;
							}

							// check if current overlap is better (longer/more unique bases ?) then curBest

							int curUniqBasesInAIvl = (ae2 - ab2) - getRepeatBasesOfContigRange(ctx, ovl, conAId);
							int curUniqBasesInBIvl = (be2 - bb2) - getRepeatBasesOfContigRange(ctx, ovl, conBId);

							if (curBestIntersection > MAX(intersect(ab1, ae1, ab2, ae2), intersect(bb1, be1, bb2, be2)) && curBestBases < MIN(ae2 - ab2, be2 - bb2))
							{
								curBestBases = MIN(ae2 - ab2, be2 - bb2);
								curBestOffset = i - longestUniqOvlIdx;
								curBestIntersection = MAX(intersect(ab1, ae1, ab2, ae2), intersect(bb1, be1, bb2, be2));
							}

							if (curBestUniqBases < MIN(curUniqBasesInAIvl, curUniqBasesInBIvl))
							{
	#ifdef DEBUG_CHAIN
								printf("found right current best overlap %d %d %d: [%d, %d] [%d, %d] right side\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								curBestUniqOffset = i - longestUniqOvlIdx;
								curBestUniqBases = MIN(curUniqBasesInAIvl, curUniqBasesInBIvl);
							}
							else if (curBestUniqBases == -1 && stepSize + STRIDE_OVL_LOOKAHEAD > MAX_OVL_LOOKAHEAD) // for repetitive genomes
							{
								Overlap *tmpOvl = ovls + (longestUniqOvlIdx + curBestOffset);

								if ((intersect(ab1, ae1, tmpOvl->path.abpos, tmpOvl->path.aepos) < ae1 - tmpOvl->path.abpos)
										&& (intersect(bb1, be1, tmpOvl->path.bbpos, tmpOvl->path.bepos) < be1 - tmpOvl->path.bbpos))
								{
									curBestUniqOffset = curBestOffset;
									curBestUniqBases = 1;
								}
							}
							else
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> cannot be anchored (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
										stepSize);
	#endif
							}
						}
					}

					// check if overlap can be used to extend overlap group on the right side
					if (curBestUniqBases < 0) // i.e. there was no good overlap at right side
					{
	#ifdef DEBUG_CHAIN
						printf("could not extend ovlgroup on right side with proper overlap (with stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
	#endif
						break;
					}

					/// todo further sanity check necessary ???
					ab2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.abpos;
					ae2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.aepos;

					bb2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.bbpos;
					be2 = ovls[longestUniqOvlIdx + curBestUniqOffset].path.bepos;
	#ifdef DEBUG_CHAIN
					printf("extend ovlgroup with (right): %d %d %d: [%d, %d] [%d, %d] stepSize %d\n", longestUniqOvlIdx + curBestUniqOffset,
							ovls[longestUniqOvlIdx + curBestUniqOffset].aread, ovls[longestUniqOvlIdx + curBestUniqOffset].bread, ab1, ae1, ab2, ae2,
							stepSize - STRIDE_OVL_LOOKAHEAD);
	#endif

					if (curChain->novl == curChain->maxOvl)
					{
						curChain->maxOvl = curChain->maxOvl * 1.2 + n;
						curChain->ovls = (Overlap**) realloc(curChain->ovls, sizeof(Overlap*) * curChain->maxOvl);
						bzero(curChain->ovls + curChain->novl, sizeof(Overlap*) * (curChain->maxOvl - curChain->novl));
					}

					// append right side overlaps at the end of chain, i.e. chain must be sorted afterwards by abpos
					curChain->ovls[curChain->novl] = ovls + (longestUniqOvlIdx + curBestUniqOffset);
					curChain->ovls[curChain->novl]->flags |= OVL_TEMP;
					curChain->novl++;
					nremain--;
	#ifdef DEBUG_CHAIN
					printf("chain: nOvl: %d, maxOvl %d nremain %d\n", curChain->novl, curChain->maxOvl, nremain);
	#endif

					ab1 = ab2;
					ae1 = ae2;
					bb1 = bb2;
					be1 = be2;

					curBestUniqOffset++;
					curBestOffset = curBestUniqOffset;

					curBestUniqBases = -1;
					curBestBases = -1;

					curBestIntersection = MAX(conA->len, conB->len);

					if (longestUniqOvlIdx + curBestUniqOffset >= n)
					{
						cont = 0;
					}
				}
			}

			// try to extend chain into left direction
			if (nremain && longestUniqOvlIdx > 0)
			{
				ab1 = ovls[longestUniqOvlIdx].path.abpos;
				ae1 = ovls[longestUniqOvlIdx].path.aepos;

				bb1 = ovls[longestUniqOvlIdx].path.bbpos;
				be1 = ovls[longestUniqOvlIdx].path.bepos;

	#ifdef DEBUG_CHAIN
				printf("extend longest overlap in left direction\n");
	#endif
				// 2nd left side
				int cont = 1;
				int curBestUniqOffset = 1;
				int curBestUniqBases = -1;
				int curBestBases = -1;
				int curBestOffset = 1;
				int curBestIntersection = MAX(conA->len, conB->len);

				while (cont)
				{
					int stepSize;
					for (stepSize = MIN_OVL_LOOKAHEAD; stepSize <= MAX_OVL_LOOKAHEAD && curBestUniqBases == -1; stepSize += STRIDE_OVL_LOOKAHEAD)
					{
#ifdef DEBUG_CHAIN
					printf("FOR LOOP stepsize %d\n", stepSize);
#endif

						// try to find next best overlap with lookahead of stepSize bases
						for (i = longestUniqOvlIdx - curBestUniqOffset; i >= 0; --i)
						{
							Overlap * ovl = ovls + i;
	#ifdef DEBUG_CHAIN
							printf("LEFT: Check ovl: a[%d, %d] b[%d,%d]\n", ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos);
	#endif

							ab2 = ovl->path.abpos;
							ae2 = ovl->path.aepos;

							bb2 = ovl->path.bbpos;
							be2 = ovl->path.bepos;

							if ((ovl->flags & OVL_COMP) != (ovls[longestUniqOvlIdx].flags & OVL_COMP))
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> different orientations\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ovl->flags & OVL_CONT)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> really contained repeat\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ovl->flags & OVL_TEMP)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> is part of another chain\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ovl->flags & OVL_DISCARD)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> discarded\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							// todo mark those as discard ????
							if (contained(ab2, ae2, ab1, ae1))
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in A-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							// todo mark those as discard ????
							if (contained(bb2, be2, bb1, be1))
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> contained repeat in B-interval\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								continue;
							}

							if (ab2 > ab1 || bb2 > bb1) // also overlap must follow left direction
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> improper left extension direction (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2,
										bb2, be2, stepSize);
	#endif
								continue;
							}

							if (MAX(ab1 - ae2, bb1 - be2) > stepSize)
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> gap size too large (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
										stepSize);
	#endif
								continue;
							}

							if (ae2 - ab1 > ab1 - ab2 || be2 - bb1 > bb1 - bb2) // at least 50% overhang
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> overhang too short (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
										stepSize);
	#endif
								continue;
							}

							// check if current overlap is better (longer/more unique bases ?) then curLeftBest

							int curUniqBasesInAIvl = (ae2 - ab2) - getRepeatBasesOfContigRange(ctx, ovl, conAId);
							int curUniqBasesInBIvl = (be2 - bb2) - getRepeatBasesOfContigRange(ctx, ovl, conBId);

							if (curBestIntersection > MAX(intersect(ab2, ae2, ab1, ae1), intersect(bb2, be2, bb1, be1)) && curBestBases < MIN(ae2 - ab2, be2 - bb2))
							{
								curBestBases = MIN(ae2 - ab2, be2 - bb2);
								curBestOffset = longestUniqOvlIdx - i;
								curBestIntersection = MAX(intersect(ab2, ae2, ab1, ae1), intersect(bb2, be2, bb1, be1));
							}

							if (curBestUniqBases < MIN(curUniqBasesInAIvl, curUniqBasesInBIvl))
							{
	#ifdef DEBUG_CHAIN
								printf("found left current best overlap %d %d %d: [ab2 %d, ae2 %d] [bb2 %d, be2 %d] left side\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2);
	#endif
								curBestUniqOffset = longestUniqOvlIdx - i;
								curBestUniqBases = curUniqBasesInAIvl + curUniqBasesInBIvl;
							}
							else if (curBestUniqBases == -1 && stepSize + STRIDE_OVL_LOOKAHEAD > MAX_OVL_LOOKAHEAD) // for repetitive genomes
							{
								Overlap *tmpOvl = ovls + (longestUniqOvlIdx - curBestOffset);

								if ((intersect(tmpOvl->path.abpos, tmpOvl->path.aepos, ab1, ae1) < ae1 - tmpOvl->path.abpos)
										&& (intersect(tmpOvl->path.bbpos, tmpOvl->path.bepos, bb1, be1) < be1 - tmpOvl->path.bbpos))
								{
									curBestUniqOffset = curBestOffset;
									curBestUniqBases = 1;
								}
							}
							else
							{
	#ifdef DEBUG_CHAIN
								printf("ignore overlap %d %d %d: [%d, %d] [%d, %d] --> cannot be anchored (stepSize %d)\n", i, ovl->aread, ovl->bread, ab2, ae2, bb2, be2,
										stepSize);
	#endif
							}
						}
					}

					// check if overlap can be used to extend overlap group on the left side
					if (curBestUniqBases < 0) // i.e. there was no good overlap at left side
					{
	#ifdef DEBUG_CHAIN
						printf("could not extend ovlgroup on left side with proper overlap (stepSize %d)\n", stepSize - STRIDE_OVL_LOOKAHEAD);
	#endif
						break;
					}

					/// todo further sanity check necessary ???
					ab2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.abpos;
					ae2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.aepos;

					bb2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.bbpos;
					be2 = ovls[longestUniqOvlIdx - curBestUniqOffset].path.bepos;

	#ifdef DEBUG_CHAIN
					printf("extend ovlgroup with (left): %d %d %d: [%d, %d] [%d, %d] with stepSize %d\n", longestUniqOvlIdx - curBestUniqOffset,
							ovls[longestUniqOvlIdx - curBestUniqOffset].aread, ovls[longestUniqOvlIdx - curBestUniqOffset].bread, ab1, ae1, ab2, ae2,
							stepSize - STRIDE_OVL_LOOKAHEAD);
	#endif

					if (curChain->novl == curChain->maxOvl)
					{
						curChain->maxOvl = curChain->maxOvl * 1.2 + n;
						curChain->ovls = (Overlap**) realloc(curChain->ovls, sizeof(Overlap*) * curChain->maxOvl);
						bzero(curChain->ovls + curChain->novl, sizeof(Overlap*) * (curChain->maxOvl - curChain->novl));
					}

					// append left side overlaps at the end of chain, i.e. chain must be sorted afterwards by abpos
					curChain->ovls[curChain->novl] = ovls + (longestUniqOvlIdx - curBestUniqOffset);
					curChain->ovls[curChain->novl]->flags |= OVL_TEMP;
					curChain->novl++;
					nremain--;
	#ifdef DEBUG_CHAIN
					printf("chain: nOvl: %d, maxOvl %d nremain %d\n", chain->novl, chain->maxOvl, nremain);
	#endif

					ab1 = ab2;
					ae1 = ae2;
					bb1 = bb2;
					be1 = be2;

					curBestUniqOffset++;
					curBestOffset = curBestUniqOffset;

					curBestUniqBases = -1;
					curBestBases = -1;

					curBestIntersection = MAX(conA->len, conB->len);

					if (longestUniqOvlIdx - curBestUniqOffset < 0)
					{
						cont = 0;
					}
				}

				if (curChain->novl > 1)
				{	// sort chain
					qsort(curChain->ovls, curChain->novl, sizeof(Overlap*), cmp_ovls_abeg);
				}
			}
#ifdef DEBUG_CHAIN
		printf("chain: nOvl: %d, maxOvl %d nremain: %d\n", curChain->novl, curChain->maxOvl, nremain);
#endif

			// find possible ovls that could be added to chain (i.e. fill gaps)

			if (curChain->novl > 1 && nremain > 0)
			{
	#ifdef DEBUG_CHAIN
				printf("find possible ovls that could be added to chain (i.e. fill gaps)\n");
	#endif
				int chainIdx = 0;
				int chainLastIdx = curChain->novl - 1;
				int j;
				for (i = 0; i < n; i++)
				{
					Overlap *ovl = ovls + i;
					if ((ovl->flags & (OVL_TEMP | OVL_CONT | OVL_DISCARD)) || ((ovl->flags & OVL_COMP) != (curChain->ovls[chainIdx]->flags & OVL_COMP)))
						continue;

					if (ovl->path.abpos < curChain->ovls[chainIdx]->path.abpos)
						continue;

					if (ovl->path.abpos > curChain->ovls[chainLastIdx]->path.abpos)
						break;

					for (j = chainIdx; j < chainLastIdx; j++)
					{
						if (curChain->ovls[j]->path.aepos <= ovl->path.abpos && curChain->ovls[j + 1]->path.abpos >= ovl->path.aepos
								&& curChain->ovls[j]->path.bepos <= ovl->path.bbpos && curChain->ovls[j + 1]->path.bbpos >= ovl->path.bepos)
						{
							Overlap *lastAddedOvl = curChain->ovls[curChain->novl - 1];

							if (intersect(ovl->path.abpos, ovl->path.aepos, lastAddedOvl->path.abpos, lastAddedOvl->path.aepos) > 100
									|| intersect(ovl->path.bbpos, ovl->path.bepos, lastAddedOvl->path.bbpos, lastAddedOvl->path.bepos) > 100)
								break;

							if (curChain->novl == curChain->maxOvl)
							{
								curChain->maxOvl = curChain->maxOvl * 1.2 + n;
								curChain->ovls = (Overlap**) realloc(curChain->ovls, sizeof(Overlap*) * curChain->maxOvl);
								bzero(curChain->ovls + curChain->novl, sizeof(Overlap*) * (curChain->maxOvl - curChain->novl));
							}

							// append left side overlaps at the end of chain, i.e. chain must be sorted afterwards by abpos
							ovl->flags |= OVL_TEMP;
							curChain->ovls[curChain->novl] = ovl;
							curChain->novl++;
							nremain--;
	#ifdef DEBUG_CHAIN
							printf("chain: nOvl: %d, maxOvl %d nremain %d\n", chain->novl, chain->maxOvl, nremain);
	#endif
						}

						if (ovl->path.abpos > curChain->ovls[j + 1]->path.abpos)
							chainIdx++;
					}
				}

				if (chainLastIdx < curChain->novl - 1)
				{
					qsort(curChain->ovls, curChain->novl, sizeof(Overlap*), cmp_ovls_abeg);
				}
			}

			if (nremain)
			{
				// mark remaining ovls as DISCARD if they overlap with a chain overlap !!
	#ifdef DEBUG_CHAIN
				printf("// mark remaining ovls as DISCARD if they overlap with a chain overlap !!\n");
	#endif
				int chainIdx = 0;
				int chainLastIdx = curChain->novl - 1;
				int j;
				for (i = 0; i < n; i++)
				{
					Overlap *ovl = ovls + i;
					if ((ovl->flags & (OVL_TEMP | OVL_CONT | OVL_DISCARD)) || ((ovl->flags & OVL_COMP) != (curChain->ovls[chainIdx]->flags & OVL_COMP)))
						continue;

					for (j = chainIdx; j <= chainLastIdx; j++)
					{
						if (intersect(curChain->ovls[j]->path.abpos, curChain->ovls[j]->path.aepos, ovl->path.abpos, ovl->path.aepos)
								|| intersect(curChain->ovls[j]->path.bbpos, curChain->ovls[j]->path.bepos, ovl->path.bbpos, ovl->path.bepos))
						{
							ovl->flags |= OVL_DISCARD;
							nremain--;
	#ifdef DEBUG_CHAIN
							printf("DISCARD [%d, %d] [%d, %d] nremain %d\n", ovl->path.abpos, ovl->path.aepos, ovl->path.bbpos, ovl->path.bepos, nremain);
	#endif
							break;
						}

						if (j + 1 < curChain->novl && ovl->path.abpos > curChain->ovls[j + 1]->path.abpos)
							chainIdx++;
					}
				}
			}
#ifdef DEBUG_CHAIN
		printChain(chain);
#endif
			// sanity check // there should be no intersection with other chains (with same orientation) possible
			int valid = 1;
			if (ctx->curChains)
			{
#ifdef DEBUG_CHAIN
				printf("DO SANITY CHECK\n");
#endif
				int j;
				for (i = 0; i < ctx->curChains && valid; i++)   // loop through all existing chains
				{
					Chain *validChain = ctx->ovlChains + i;

					// check contigA coordinates
					if(intersect(validChain->ovls[0]->path.abpos,validChain->ovls[validChain->novl-1]->path.aepos,curChain->ovls[0]->path.abpos,curChain->ovls[curChain->novl-1]->path.aepos))
					{
#ifdef DEBUG_CHAIN
							printf("CHAIN is invalid - DISCARD\n");
#endif
							valid = 0;
							break;
					}
					// check contigB coordinates
					if ((validChain->ovls[0]->flags & OVL_COMP) == (curChain->ovls[0]->flags && OVL_COMP))
					{
						if(intersect(validChain->ovls[0]->path.bbpos,validChain->ovls[validChain->novl-1]->path.bepos,curChain->ovls[0]->path.bbpos,curChain->ovls[curChain->novl-1]->path.bepos))
						{
	#ifdef DEBUG_CHAIN
								printf("CHAIN is invalid - DISCARD\n");
	#endif
							valid = 0;
							break;
						}
					}
					else if (validChain->ovls[0]->flags & OVL_COMP)
					{
						if(intersect(conB->len - validChain->ovls[validChain->novl-1]->path.bepos, conB->len - validChain->ovls[0]->path.bbpos,curChain->ovls[0]->path.bbpos,curChain->ovls[curChain->novl-1]->path.bepos))
						{
	#ifdef DEBUG_CHAIN
								printf("CHAIN is invalid - DISCARD\n");
	#endif
							valid = 0;
							break;
						}
					}
					else // i.e. (curChain->ovls[0]->flags && OVL_COMP)
					{
						if(intersect(validChain->ovls[0]->path.bbpos,validChain->ovls[validChain->novl-1]->path.bepos,conB->len - curChain->ovls[curChain->novl-1]->path.bepos, conB->len - curChain->ovls[0]->path.bbpos))
						{
	#ifdef DEBUG_CHAIN
								printf("CHAIN is invalid - DISCARD\n");
	#endif
							valid = 0;
							break;
						}
					}
				}
			}

			int aCoveredBases=0;
			int bCoveredBases=0;
			for (i = 0; i < curChain->novl; i++)
			{
				aCoveredBases+=(curChain->ovls[i]->path.aepos - curChain->ovls[i]->path.abpos);
				bCoveredBases+=(curChain->ovls[i]->path.aepos - curChain->ovls[i]->path.abpos);
			}

			if(
				(aCoveredBases < 0.5*conA->len && bCoveredBases < 0.5*conB->len)
				&&
				!(
						((curChain->ovls[0]->path.abpos == 0) || (curChain->ovls[curChain->novl-1]->path.aepos = conA->len)) &&
						((curChain->ovls[0]->path.bbpos == 0) || (curChain->ovls[curChain->novl-1]->path.bepos = conB->len)) &&
						MIN(aCoveredBases, bCoveredBases) > 15000
				 )
			)
			{
				valid = 0;
			}

			if (valid)
				ctx->curChains++;
			else
			{
				int j;
				for (j = 0; j < curChain->novl; j++)
				{
					curChain->ovls[j]->flags |= OVL_DISCARD;
				}
				curChain->novl = 0;
			}

#ifdef DEBUG_CHAIN
			printf("curChain: %d, remain unchained OVls: %d\n", ctx->curChains, nremain);
#endif
		}

#ifdef DEBUG_CHAIN
		printf("FINAL CHAINS: %d %7d vs %7d\n", ctx->curChains, ctx->ovlChains[0].ovls[0]->aread, ctx->ovlChains[0].ovls[0]->bread);

		for (i = 0; i < ctx->curChains; i++)
		{
			printf(" CHAIN %d/%d: #novl %d\n", i + 1, ctx->curChains, ctx->ovlChains[0].novl);

			int j;
			for (j = 0; j < ctx->ovlChains[i].novl; j++)
			{
				printf("  OVL %d/%d: a[%7d, %7d] b[%7d, %7d] %s\n", j + 1, ctx->ovlChains[i].novl, ctx->ovlChains[i].ovls[j]->path.abpos,
						ctx->ovlChains[i].ovls[j]->path.aepos, ctx->ovlChains[i].ovls[j]->path.bbpos, ctx->ovlChains[i].ovls[j]->path.bepos,
						(ctx->ovlChains[i].ovls[j]->flags & OVL_COMP) ? "COMP" : "NORM");
			}
		}
#endif

		// sort chains according to alignment lengths
		if (ctx->curChains > 1)
		{
	#ifdef DEBUG_CHAIN
			printf("SORT CHAINS (longest first):\n");
	#endif
			qsort(ctx->ovlChains, ctx->curChains, sizeof(Chain), cmp_chain_len);
	#ifdef DEBUG_CHAIN
			printf("FINAL CHAINS: %d %7d vs %7d\n", ctx->curChains, ctx->ovlChains[0].ovls[0]->aread, ctx->ovlChains[0].ovls[0]->bread);
			for (i = 0; i < ctx->curChains; i++)
			{
				printf(" CHAIN %d/%d: #novl %d\n", i + 1, ctx->curChains, ctx->ovlChains[0].novl);
				int j;
				for (j = 0; j < ctx->ovlChains[i].novl; j++)
				{
					printf("  OVL %d/%d: a[%7d, %7d] b[%7d, %7d] %s\n", j + 1, ctx->ovlChains[i].novl, ctx->ovlChains[i].ovls[j]->path.abpos,
							ctx->ovlChains[i].ovls[j]->path.aepos, ctx->ovlChains[i].ovls[j]->path.bbpos, ctx->ovlChains[i].ovls[j]->path.bepos,
							(ctx->ovlChains[i].ovls[j]->flags & OVL_COMP) ? "COMP" : "NORM");
				}
			}
	#endif
		}
	}
}

int analyzeChains(AnalyzeContext *ctx)
{
	if(ctx->curChains == 0)
		return 0;

	printf("[analyzeChains] c%d vs c%d\n", ctx->ovlChains[0].ovls[0]->aread, ctx->ovlChains[0].ovls[0]->bread);
	int i,j;
	for (i=0; i < ctx->curChains; i++)
	{
		Chain *chain = ctx->ovlChains + i;
		if(chain->novl == 0)
			return 0;

		printf("	Chain_%d [", i);
		for (j=0; j<chain->novl; j++)
		{
			printf("%d-%d", chain->ovls[j]->path.abpos,chain->ovls[j]->path.aepos);
			if(j+1 < chain->novl)
				printf(",");
		}
		printf("]\n");
	}
	return 0;
}

int processContigOverlap_handler(void* _ctx, Overlap* ovls, int novl)
{
	AnalyzeContext* actx = (AnalyzeContext*) _ctx;

	if (actx->contigs[ovls->aread].flag & CONTIG_DISCARD)
		return 1;

	if (actx->VERBOSE)
		printf("Analyze overlaps for contig: %d numOvls: %d\n", ovls->aread, novl);

	int j;

	int k;
	j = k = 0;
	while (j < novl)
	{
		while (k < novl - 1 && ovls[j].bread == ovls[k + 1].bread)
			k++;

		if (ovls[j].aread == ovls[j].bread) // ignore overlaps between same contig
		{
			// TODO discard overlaps only when writeOutFilteredContigChains is enabled
			int i;
			for (i = 0; i < (k - j + 1); i++)
			{
				ovls[i].flags |= OVL_DISCARD;
			}

		}
//		else if (DB_READ_LEN(actx->corContigDB, ovls[j].aread) > DB_READ_LEN(actx->corContigDB, ovls[j].bread)) // len aread must be < then len bread !!!
//			;
		else
		{
			chainContigOverlaps(actx, ovls + j, k - j + 1);
			if(analyzeChains(actx))
			{
				// convert valid chains into ContigChain struct
				Contig *conA = actx->contigs + ovls->aread;
				Contig *conB = actx->contigs + ovls->bread;
				int i;
				for (i=0; i < actx->curChains; i++)
				{
					Chain *chain = actx->ovlChains;
					if(chain->novl == 0)
						continue;

					if(conA->numCChains == conA->maxCChains)
					{
						conA->maxCChains = conA->maxCChains * 1.2 + 5;
						conA->cChains = (ContigChain*) realloc(conA->cChains, sizeof(ContigChain) * conA->maxCChains);
						bzero(conA->cChains + conA->numCChains, sizeof(ContigChain) * (conA->maxCChains - conA->numCChains));
					}

					ContigChain *cchain = conA->cChains + conA->numCChains;
					cchain->corContigIdx = conB->idx;
					cchain->numPos = chain->novl;
					cchain->abpos = (int*)malloc(sizeof(int)*chain->novl*2);
					cchain->aepos = cchain->abpos + chain->novl;
					int h;
					for (h=0; h<chain->novl; h++)
					{
						cchain->abpos[h] = chain->ovls[h]->path.abpos;
						cchain->aepos[h] = chain->ovls[h]->path.aepos;
					}
				}
			}
			// reset chain and ovl counter
			int i;
			for (i = 0; i < actx->curChains; i++)
			{
				actx->ovlChains[i].novl = 0;
			}
			actx->curChains = 0;
		}

		j = k + 1;
	}
	return 1;
}

void classifyContigsByOverlaps(AnalyzeContext *actx)
{
// sort contigs according length (smallest first)
	int i, j;
	int repBasesAbeg, repBasesBbeg;
	int repBasesAend, repBasesBend;

	for (i = 0; i < actx->numContigs; i++)
	{
		//printf("%5d contig: %d, len: %d", i, actx->contigs[i].idx, actx->contigs[i].len);

		Contig *contigA = actx->contigs + i;

		if (contigA->flag & CONTIG_DISCARD)
			continue;

		// no valid overlap groups --> its unique
//		if (contigA->numOvlGrps == 0)
//		{
//			contigA->ovlGrpsClassificFlag |= CONTIG_UNIQUE;
//			continue;
//		}

		// analyze all overlap groups
		int overAllNotClassified = 1;
//		for (j = 0; j < actx->contigs[i].numOvlGrps; j++)
//		{
//			OverlapGroup * olg = contigA->ovlGrps[j];
////			printOverlapGroup(stdout, actx, olg);
//
//			// get repeat bases in front and behind overlap group (aread)
//			repBasesAbeg = getRepeatBasesFromInterval(actx, contigA->idx, 0, olg->first_abpos);
//			repBasesAend = getRepeatBasesFromInterval(actx, contigA->idx, olg->last_aepos, contigA->len);
//
//			Contig *contigB = actx->contigs + olg->bread;
//
//			if (contigA->flag & CONTIG_DISCARD)
//			{
//				olg->flag |= OVLGRP_DISCARD;
//				continue;
//			}
//
//			// get overall repeat bases of bread
//			// get repeat bases in front and behind overlap group (bread)
//			if (olg->flag & OVLGRP_COMP)
//			{
//				repBasesBbeg = getRepeatBasesFromInterval(actx, contigB->idx, 0, olg->last_bbpos);
//				repBasesBend = getRepeatBasesFromInterval(actx, contigB->idx, olg->first_bepos, contigB->len);
//			}
//			else
//			{
//				repBasesBbeg = getRepeatBasesFromInterval(actx, contigB->idx, 0, olg->first_bbpos);
//				repBasesBend = getRepeatBasesFromInterval(actx, contigB->idx, olg->last_bepos, contigB->len);
//			}
//
//			int spanningA, spanningB;
//			spanningA = olg->last_aepos - olg->first_abpos;
//			if (olg->flag & OVLGRP_COMP)
//				spanningB = olg->first_bepos - olg->last_bbpos;
//			else
//				spanningB = olg->last_bepos - olg->first_bbpos;
//
////			printf("c %d vs c %d \n", olg->aread, olg->bread);
////			printf("spanning a %d , b: %d\n", spanningA, spanningB);
////			printf("anchors a %d , b: %d\n", olg->coveredBasesInA - olg->repeatBasesInA, olg->coveredBasesInB - olg->repeatBasesInB);
////			printf("repA: %.2f, repB: %.2f\n", contigA->repBases * 100.0 / contigA->len, contigB->repBases * 100.0 / contigB->len);
////			if (olg->flag & OVLGRP_COMP)
////			{
////				printf("unaligned a [%d, %d] , b: [%d, %d]\n", olg->first_abpos, contigA->len - olg->last_aepos, olg->last_bbpos,
////						contigB->len - olg->first_bepos);
////				printf("unaligned repeat bases A before [%d, %.2f] after [%d, %.2f]\n", repBasesAbeg, 100.0 * repBasesAbeg / olg->first_abpos, repBasesAend,
////						100.0 * repBasesAend / (contigA->len - olg->last_aepos));
////				printf("unaligned repeat bases B before [%d, %.2f] after [%d, %.2f]\n", repBasesBbeg, 100.0 * repBasesBbeg / olg->last_bbpos, repBasesBend,
////						100.0 * repBasesBend / (contigB->len - olg->first_bepos));
////			}
////			else
////			{
////				printf("unaligned a [%d, %d] , b: [%d, %d]\n", olg->first_abpos, contigA->len - olg->last_aepos, olg->first_bbpos,
////						contigB->len - olg->last_bepos);
////				printf("unaligned repeat bases A before [%d, %.2f] after [%d, %.2f]\n", repBasesAbeg, 100.0 * repBasesAbeg / olg->first_abpos, repBasesAend,
////						100.0 * repBasesAend / (contigA->len - olg->last_aepos));
////				printf("unaligned repeat bases B before [%d, %.2f] after [%d, %.2f]\n", repBasesBbeg, 100.0 * repBasesBbeg / olg->first_bbpos, repBasesBend,
////						100.0 * repBasesBend / (contigB->len - olg->last_bepos));
////			}
//
//			if (spanningA < ConVsConMinAlign && spanningB < ConVsConMinAlign)
//			{
//				olg->flag |= OVLGRP_TOO_SHORT;
//				olg->flag |= OVLGRP_DISCARD;
//				continue;
//			}
//
//			int notClassified = 1;
//			if (olg->coveredBasesInA - olg->repeatBasesInA > ConVsConAnchorBases || olg->coveredBasesInB - olg->repeatBasesInB > ConVsConAnchorBases)
//			{
//				// 70 % of the bread overlaps with aread
//				if (spanningB * 1.0 / contigB->len >= ConVsConContainmentThreshold)
//				{
////					printf("contig %d is contained >= %3f%% in contig %d\n", contigB->idx, 100*ConVsConContainmentThreshold, contigA->idx, );
////					contigB->ovlGrpsClassificFlag |= CONTIG_IS_CONTAINED;
////					contigA->ovlGrpsClassificFlag |= CONTIG_HAS_CONTAINED;
//
//					olg->flag |= OVLGRP_AREAD_HAS_CONTAINED;
//					notClassified = 0;
//				}
//				// 70 % of the repeat-trimmed bread overlaps with aread, end begin and end of bread are tagged repetitive
//				else if (spanningB * 1.0 / (contigB->len - repBasesBbeg - repBasesBend) >= ConVsConContainmentThreshold)
//				{
////					printf("repeat trimmed contig %d is contained >= %3f%% in contig %d\n", contigB->idx, 100*ConVsConContainmentThreshold, contigA->idx);
////					contigB->ovlGrpsClassificFlag |= CONTIG_IS_CONTAINED;
////					contigA->ovlGrpsClassificFlag |= CONTIG_HAS_CONTAINED;
//
//					olg->flag |= OVLGRP_AREAD_HAS_CONTAINED;
//					notClassified = 0;
//				}
//				// for short contigs < 100k allow 50% match in the middle of contigA to be contained
//				else if (contigB->len < ConVsConShortContigLen && spanningB * 1.0 / contigB->len >= ConVsConShortContigContainmentThreshold)
//				{
////					printf("contig %d is short (< %d) and contained >= %3f%% in contig %d\n", contigB->idx, ConVsConShortContigLen, 100*ConVsConShortContigContainmentThreshold, contigA->idx);
////					contigB->ovlGrpsClassificFlag |= CONTIG_IS_CONTAINED;
////					contigA->ovlGrpsClassificFlag |= CONTIG_HAS_CONTAINED;
//
//					olg->flag |= OVLGRP_AREAD_HAS_CONTAINED;
//					notClassified = 0;
//				}
//
//				// 70 % of the aread overlaps with bread
//				if (spanningA * 1.0 / contigA->len >= ConVsConContainmentThreshold)
//				{
////					printf("contig %d is contained > %3f%% in contig %d\n", contigA->idx, 100*ConVsConContainmentThreshold, contigB->idx);
////					contigA->ovlGrpsClassificFlag |= CONTIG_IS_CONTAINED;
////					contigB->ovlGrpsClassificFlag |= CONTIG_HAS_CONTAINED;
//
//					olg->flag |= OVLGRP_AREAD_IS_CONTAINED;
//					notClassified = 0;
//				}
//
//				// 70 % of the aread overlaps with bread, end begin and end of aread are tagged repetitive
//				else if (spanningA * 1.0 / (contigA->len - repBasesAbeg - repBasesAend) >= ConVsConContainmentThreshold)
//				{
////					printf("repeat trimmed contig %d is contained > %3f%% in contig %d\n", contigA->idx, 100*ConVsConContainmentThreshold, contigB->idx);
////					contigA->ovlGrpsClassificFlag |= CONTIG_IS_CONTAINED;
////					contigB->ovlGrpsClassificFlag |= CONTIG_HAS_CONTAINED;
//
//					olg->flag |= OVLGRP_AREAD_IS_CONTAINED;
//					notClassified = 0;
//				}
//				// for short contigs < 100k allow 50% match in the middle of contigA to be contained
//				else if (contigA->len < ConVsConShortContigLen && spanningA * 1.0 / contigA->len >= ConVsConShortContigContainmentThreshold)
//				{
////					printf("contig %d is short (< %d) and contained >= %3f%% in contig %d\n", contigB->idx, ConVsConShortContigLen, 100*ConVsConShortContigContainmentThreshold, contigA->idx);
////					contigA->ovlGrpsClassificFlag |= CONTIG_IS_CONTAINED;
////					contigB->ovlGrpsClassificFlag |= CONTIG_HAS_CONTAINED;
//
//					olg->flag |= OVLGRP_AREAD_IS_CONTAINED;
//					notClassified = 0;
//				}
//			}
//			if (notClassified)
//			{
//				olg->flag |= OVLGRP_DISCARD;
//			}
//			else
//			{
//				overAllNotClassified = 0;
//			}
//		}
		if (overAllNotClassified)
		{
//			contigA->ovlGrpsClassificFlag |= (CONTIG_UNIQUE | CONTIG_UNCLASSIFIED);
		}
	}
}

int createOutDir(char *out)
{
	int res = 0;

	if (out)
	{
		struct stat s;

		int err = stat(out, &s);

		if (err == -1)
		{
			if (errno == ENOENT)
			{
				err = mkdir(out, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
				if (err != 0)
				{
					fprintf(stderr, "Cannot create output directory: %s\n", out);
					res = 1;
				}
			}
			else
			{
				fprintf(stderr, "Cannot create output directory: %s\n", out);
				res = 1;
			}
		}
		else
		{
			if (!S_ISDIR(s.st_mode))
			{
				fprintf(stderr, "Output directory name: \"%s\" exist - but its not a directory\n", out);
				res = 1;
			}
		}
	}
	return res;
}

static void writeContigInfoFile(AnalyzeContext *actx, Contig* contig, FILE* contigFile, int cBegPos, int cEndPos, bit *allReads, int *contigCov)
{
	ba_all_assign(allReads, DB_NREADS(actx->patchedReadDB), FALSE);
	bzero(contigCov, sizeof(int) * DB_READ_MAXLEN(actx->corContigDB));

// determine all overlapping reads and coverage of contig
	int i, j, k;
	for (i = 0; i < contig->numcReads; i++)
	{
		ContigRead *cread = contig->cReads + i;

		if (intersect(cread->ovlReads->cBeg, cread->ovlReads->cEnd, cBegPos, cEndPos))
		{
			for (j = 0; j < cread->numOvlReads; j++)
			{
				if (intersect(cread->ovlReads[j].cBeg, cread->ovlReads[j].cEnd, cBegPos, cEndPos))
				{
					ba_assign(allReads, abs(cread->ovlReads[j].patchedID), TRUE);

					int covBeg = MAX(cread->ovlReads[j].cBeg, cBegPos);
					int covEnd = MIN(cread->ovlReads[j].cEnd, cEndPos);
					for (k = covBeg; k < covEnd; k++)
						contigCov[k]++;
				}
			}
		}
	}

// print all patched reads that overlap with current contig
	fprintf(contigFile, "PATCH_READS");
	for (i = 0; i < DB_NREADS(actx->patchedReadDB); i++)
		if (ba_value(allReads, i))
			fprintf(contigFile, " %d", i);
	fprintf(contigFile, "\n");

// print all source reads that should overlap with current contig
	if (actx->patchedReadSource_track)
	{
		// fixed source read ids from contig
		track_anno* sreadanno = actx->patchedReadSource_track->anno;
		track_data* sreaddata = actx->patchedReadSource_track->data;

		fprintf(contigFile, "SOURCE_READS");
		for (i = 0; i < DB_NREADS(actx->patchedReadDB); i++)
			if (ba_value(allReads, i))
			{
				track_anno rb = sreadanno[i] / sizeof(track_data);
				track_anno re = sreadanno[i + 1] / sizeof(track_data);

				assert(rb < re);

				fprintf(contigFile, " %d", sreaddata[rb]);
			}
		fprintf(contigFile, "\n");
	}

// print coverage histogram
	int binSize = 1000;
	fprintf(contigFile, "COV_HIST_B%d", binSize);

	int overallCov = 0;
	for (i = cBegPos; i + binSize < cEndPos; i++)
		overallCov += contigCov[i];

	fprintf(contigFile, " %.1f", overallCov * 1.0 / (cEndPos - cBegPos));

	int cumCov;
	for (i = cBegPos; i + binSize < cEndPos; i += binSize)
	{
		cumCov = 0;
		for (k = i; k < i + binSize; k++)
			cumCov += contigCov[k];
		fprintf(contigFile, " %.1f", cumCov * 1.0 / binSize);
	}
	if (i < cEndPos)
	{
		cumCov = 0;
		binSize = cEndPos - i;
		for (k = i; k < cEndPos; k++)
			cumCov += contigCov[k];
		fprintf(contigFile, " %.1f", cumCov * 1.0 / binSize);
	}
	fprintf(contigFile, "\n");

// repeat annotation:
	fprintf(contigFile, "REPEAT_%s", actx->corContigRepeats_track->name);
	fprintf(contigFile, " %d", getRepeatBasesFromInterval(actx, contig->idx, cBegPos, cEndPos));
	{
		track_anno* rep_anno = actx->corContigRepeats_track->anno;
		track_data* rep_data = actx->corContigRepeats_track->data;

		track_anno rb, re;
		// repeat bases in a-read
		rb = rep_anno[contig->idx] / sizeof(track_data);
		re = rep_anno[contig->idx + 1] / sizeof(track_data);

		while (rb < re)
		{
			int rBeg = rep_data[rb];
			int rEnd = rep_data[rb + 1];

			if (rEnd < cBegPos)
			{
				continue;
				rb += 2;
			}
			if (rBeg > cEndPos)
				break;

			if (rBeg < cBegPos)
				rBeg = cBegPos;

			if (rEnd > cEndPos)
				rEnd = cEndPos;

			fprintf(contigFile, " %d,%d", rBeg, rEnd);

			rb += 2;
		}
	}
	fprintf(contigFile, "\n");

// report splits for raw contigs
//	if (contig->len == cEndPos - cBegPos && contig->numSplits && !(contig->split->type & SPLIT_IGNORE))
//	{
//		fprintf(contigFile, "SPLITS");
//
//		for (i = 0; i < contig->numSplits; i++)
//		{
//			if (contig->split[i].type & SPLIT_IGNORE)
//				break;
//			fprintf(contigFile, " %d,%d,%d", contig->split[i].type, contig->split[i].contigPos1, contig->split[i].contigPos2);
//		}
//		fprintf(contigFile, "\n");
//	}
}

// todo check if cEndPos is really exclusive
void writeFasta(AnalyzeContext *actx, Contig* contig, FILE* contigFile, int splitIdx, int cBegPos, int cEndPos)
{

	int i;
	int WIDTH = 100;

	int path = getPathID(actx, contig->idx);
	int end1, end2;
	getContigsEndIDs(actx, contig->idx, &end1, &end2);
	char *InFastaName = getFastaFileNameFromDB(actx, contig->idx);

///// create header
// new name = fasta file name + Contig idx
	fprintf(contigFile, ">%s_%d", InFastaName, contig->idx);

///// add some track info
// path track
	fprintf(contigFile, " path=%d", path);
// length track
	fprintf(contigFile, " length=%d", cEndPos - cBegPos);

	if (splitIdx == 0) // otherwise the ends can be invalid
	{
		// ends track
		fprintf(contigFile, " ends=%d,%d", end1, end2);
		// repeat track
		//fprintf(contigFile, " repeat=%d,%d,%d", contig->repBases, contig->repBasesTipLeft, contig->repBasesTipRight);
	}
	else
	{
		fprintf(contigFile, " splitIdx=%d", splitIdx);
	}
// reads track
	if (contig->numcReads)
	{
		fprintf(contigFile, " reads=");
		int curPos = 0;
		for (i = 0; i < contig->numcReads; i++)
		{
			if (curPos < cBegPos)
			{
				curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);
				continue;
			}
			if (curPos > cEndPos)
				break;

			// first read of split contig
			if (cBegPos > curPos && cBegPos < curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg))
			{
				if (contig->cReads[i].ovlReads->beg < contig->cReads[i].ovlReads->end) // forward orientation
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg + (cBegPos - curPos),
							contig->cReads[i].ovlReads->end);
				}
				else
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg - (cBegPos - curPos),
							contig->cReads[i].ovlReads->end);
				}
			}
			// last read of split contig
			else if (curPos < cEndPos && cEndPos < curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg))
			{
				if (contig->cReads[i].ovlReads->beg < contig->cReads[i].ovlReads->end) // forward orientation
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg,
							contig->cReads[i].ovlReads->end - (curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg) - cEndPos));
				}
				else
				{
					fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg,
							contig->cReads[i].ovlReads->end + (curPos + abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg) - cEndPos));
				}
				break;
			}
			else
			{
				fprintf(contigFile, "%d,%d,%d", abs(contig->cReads[i].patchedID), contig->cReads[i].ovlReads->beg, contig->cReads[i].ovlReads->end);
			}

			curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);
			if (curPos < cEndPos)
				fprintf(contigFile, ",");

		}
	}
// sreads track
	if (actx->corContigRawReads_track)
	{
		// fixed source read ids from contig
		track_anno* sreadanno = actx->corContigRawReads_track->anno;
		track_data* sreaddata = actx->corContigRawReads_track->data;

		track_anno rb = sreadanno[contig->idx] / sizeof(track_data);
		track_anno re = sreadanno[contig->idx + 1] / sizeof(track_data);

		assert(rb < re && (int ) (re - rb) == contig->numcReads);

		fprintf(contigFile, " sreads=");
		int curPos = 0;
		for (i = 0; i < contig->numcReads; i++, rb++)
		{
			if (curPos < cBegPos)
			{
				curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);
				continue;
			}
			if (curPos > cEndPos)
				break;

			fprintf(contigFile, "%d", sreaddata[rb]);
			curPos += abs(contig->cReads[i].ovlReads->end - contig->cReads[i].ovlReads->beg);

			if (curPos < cEndPos)
				fprintf(contigFile, ",");
		}

	}
// classification tracks
// 1)  based on fixed read overlaps
	if (contig->gClassificFlag & CONTIG_HAS_CONTAINED)
	{
		int multi = 0;
		ContigGraphClassification *gclass;
		for (i = 0; i < contig->numGClassific; i++)
		{
			gclass = contig->gClassific + i;

			if (gclass->flag & CONTIG_HAS_CONTAINED)
			{
				if (multi)
				{
					fprintf(contigFile, ",%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0], gclass->coveredIntervals[gclass->numCoveredIntervals*3-2]);
				}
				else
				{
					fprintf(contigFile, " gContains=%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0], gclass->coveredIntervals[gclass->numCoveredIntervals*3-2]);
					multi = 1;
				}
			}
		}
	}
	if (contig->gClassificFlag & CONTIG_IS_CONTAINED)
	{
		int multi = 0;
		ContigGraphClassification *gclass;
		for (i = 0; i < contig->numGClassific; i++)
		{
			gclass = contig->gClassific + i;

			if (gclass->flag & CONTIG_IS_CONTAINED)
			{
				if (multi)
				{
					fprintf(contigFile, ",%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0], gclass->coveredIntervals[gclass->numCoveredIntervals*3-2]);
				}
				else
				{
					fprintf(contigFile, " gContainedIn=%d,%d,%d", gclass->correspID, gclass->coveredIntervals[0], gclass->coveredIntervals[gclass->numCoveredIntervals*3-2]);
					multi = 1;
				}
			}
		}
	}
// 2)  based on contig overlap groups
	if (contig->flag & CONTIG_HAS_CONTAINED)
	{
		int multi = 0;
		OverlapGroup *ogr;
//		for (i = 0; i < contig->numOvlGrps; i++)
//		{
//			ogr = contig->ovlGrps[i];
//
//			if ((ogr->flag & CONTIG_HAS_CONTAINED) && contig->len > DB_READ_LEN(actx->corContigDB, ogr->bread))
//			{
//				if (multi)
//				{
//					fprintf(contigFile, ",%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
//				}
//				else
//				{
//					fprintf(contigFile, " cContains=%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
//					multi = 1;
//				}
//			}
//		}
	}
	if (contig->flag & CONTIG_IS_CONTAINED)
	{
		int multi = 0;
		OverlapGroup *ogr;
//		for (i = 0; i < contig->numOvlGrps; i++)
//		{
//			ogr = contig->ovlGrps[i];
//
//			if ((ogr->flag & CONTIG_IS_CONTAINED) && contig->len < DB_READ_LEN(actx->corContigDB, ogr->bread))
//			{
//				if (multi)
//				{
//					fprintf(contigFile, ",%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
//				}
//				else
//				{
//					fprintf(contigFile, " cContainedIn=%d,%d,%d", ogr->bread, ogr->first_abpos, ogr->last_aepos);
//					multi = 1;
//				}
//			}
//		}
	}

	if (contig->flag & CONTIG_AL50PCTCOVERED)
	{
		int multi = 0;

		for (i = 0; i < actx->numContigs; i++)
		{
//			if ((int) contig->pctCoveredBasesInOtherContigs[i] >= 50)
//			{
//				if (multi)
//				{
//					fprintf(contigFile, ",%d,%d", i, (int) contig->pctCoveredBasesInOtherContigs[i]);
//				}
//				else
//				{
//					fprintf(contigFile, " cCoveredIn=%d,%d", i, (int) contig->pctCoveredBasesInOtherContigs[i]);
//					multi = 1;
//				}
//			}
		}
	}

// todo sort splits according positions ?
// todo write out only real splits, i.e. skip heterozygoues ones, and print out only one valid split position either left or, or right or sometimes both
//	if (contig->numSplits && !(contig->split->type & SPLIT_IGNORE))
//	{
//		int multi = 0;
//		// print absolute contig Positions
//		for (i = 0; i < contig->numSplits; i++)
//		{
//			if (multi)
//			{
//				fprintf(contigFile, ",%d,%d,%d", contig->split[i].type, contig->split[i].contigPos1, contig->split[i].contigPos2);
//			}
//			else
//			{
//				fprintf(contigFile, " splitAbsPos=%d,%d,%d", contig->split[i].type, contig->split[i].contigPos1, contig->split[i].contigPos2);
//				multi = 1;
//			}
//		}
//		multi = 0;
//		// print split reads and their positions
//		for (i = 0; i < contig->numSplits; i++)
//		{
//			if (multi)
//			{
//				fprintf(contigFile, ",%d,%d,%d,%d,%d,%d", contig->split[i].leftProperRead, contig->split[i].leftProperReadPos1,
//						contig->split[i].leftProperReadPos2, contig->split[i].rightProperRead, contig->split[i].rightProperReadPos1,
//						contig->split[i].rightProperReadPos2);
//			}
//			else
//			{
//				fprintf(contigFile, " split=%d,%d,%d,%d,%d,%d", contig->split[i].leftProperRead, contig->split[i].leftProperReadPos1,
//						contig->split[i].leftProperReadPos2, contig->split[i].rightProperRead, contig->split[i].rightProperReadPos1,
//						contig->split[i].rightProperReadPos2);
//				multi = 1;
//			}
//		}
//	}

	fprintf(contigFile, "\n");

///// add sequence
	Load_Read(actx->corContigDB, contig->idx, actx->readSeq, 1);

	for (i = cBegPos; i + WIDTH < cEndPos; i += WIDTH)
		fprintf(contigFile, "%.*s\n", WIDTH, actx->readSeq + i);
	if (i < cEndPos)
		fprintf(contigFile, "%.*s\n", cEndPos - i, actx->readSeq + i);

	fflush(contigFile);
}

static int cmpInts(const void * a, const void * b)
{
	return (*(int*) b - *(int*) a);
}

//static int cmpSplitEvents(const void * a, const void * b)
//{
//	SplitEvent *sa = (SplitEvent*) a;
//	SplitEvent *sb = (SplitEvent*) b;
//
//	if ((sa->type & SPLIT_IGNORE) && (sb->type & SPLIT_IGNORE))
//		return (sa->contigPos1 - sb->contigPos1);
//	if (sa->type & SPLIT_IGNORE)
//		return 1;
//	if (sb->type & SPLIT_IGNORE)
//		return -1;
//	return (sa->contigPos1 - sb->contigPos1);
//}

/*
 void printFinalAltContigClassification(FILE *out, FinalContigAltClass facc)
 {
 switch (facc)
 {
 case BUBBLE_ALT:
 fprintf(out, "BUBBLE_ALT");
 break;
 case BUBBLE_HUGEDIF:
 fprintf(out, "BUBBLE_HUGEDIF");
 break;
 case BUBBLE_LOWCOVEREDREPEAT:
 fprintf(out, "BUBBLE_LOWCOVEREDREPEAT");
 break;
 case BUBBLE_REPEAT:
 fprintf(out, "BUBBLE_REPEAT");
 break;
 case BUBBLE_WEIRD:
 fprintf(out, "BUBBLE_WEIRD");
 break;
 case SPUR_ALT:
 fprintf(out, "SPUR_ALT");
 break;
 case SPUR_ATCONTIGEND:
 fprintf(out, "SPUR_ATCONTIGEND");
 break;
 default:
 fprintf(out, "NONE");
 break;
 }
 }

 void printFinalContigClassification(FILE *out, FinalContigClass fcc)
 {
 switch (fcc)
 {
 case HAPLOID:
 fprintf(out, "HAPLOID");
 break;
 case CONTAINED_BY_BOTH:
 fprintf(out, "CONTAINED_BY_BOTH");
 break;
 case CONTAINED_BY_BOTHSTRANGE:
 fprintf(out, "CONTAINED_BY_BOTHSTRANGE");
 break;
 case CONTAINED_BY_CONTIG:
 fprintf(out, "CONTAINED_BY_CONTIG");
 break;
 case CONTAINED_BY_FIXREAD:
 fprintf(out, "CONTAINED_BY_FIXREAD");
 break;
 case CONTAINED_BY_COVEREDBASES:
 fprintf(out, "CONTAINED_BY_COVEREDBASES");
 break;
 default:
 fprintf(out, "UNKNOWN");
 break;
 }
 }

 char getClassification(FinalContigClass fcc)
 {
 switch (fcc)
 {
 case HAPLOID:
 return 'H';
 case CONTAINED_BY_BOTH:
 return 'A';
 case CONTAINED_BY_BOTHSTRANGE:
 return 'W';
 case CONTAINED_BY_CONTIG:
 return 'C';
 case CONTAINED_BY_FIXREAD:
 return 'R';
 case CONTAINED_BY_COVEREDBASES:
 return 'B';
 default:
 return 'U';
 }
 }*/

//char* getSplitType(SplitEvent* se)
//{
//	static char* result = NULL;
//
//	if (result == NULL)
//		result = malloc(20);
//
//	result[0] = '\0';
//
//	if (se->type & SPLIT_IGNORE)
//	{
//		strcat(result, "I");
//	}
//	if (se->type & SPLIT_UNKNOWN)
//	{
//		strcat(result, "u");
//	}
//	if (se->type & SPLIT_MULTIREAD)
//	{
//		strcat(result, "M");
//	}
//	if (se->type & SPLIT_STARTEND)
//	{
//		strcat(result, "S");
//	}
//	if (se->type & SPLIT_FALSEJOINLEFT)
//	{
//		strcat(result, "f");
//	}
//	if (se->type & SPLIT_FALSEJOINRIGHT)
//	{
//		strcat(result, "F");
//	}
//	if (se->type & SPLIT_HETEROZYGLEFT)
//	{
//		strcat(result, "h");
//	}
//	if (se->type & SPLIT_HETEROZYGRIGHT)
//	{
//		strcat(result, "H");
//	}
//	if (se->type & SPLIT_DEADENDLEFT)
//	{
//		strcat(result, "d");
//	}
//	if (se->type & SPLIT_DEADENDRIGHT)
//	{
//		strcat(result, "D");
//	}
//
//	return result;
//}
static int cmpOVLGroupByPos(const void* a, const void* b)
{
	OverlapGroup* o1 = *(OverlapGroup**) a;
	OverlapGroup* o2 = *(OverlapGroup**) b;

	if (o1->flag & OVLGRP_DISCARD)
		return 1;

	if (o2->flag & OVLGRP_DISCARD)
		return -1;

	return (o1->first_abpos - o2->first_abpos);
}

static int cmpGraphClassByPos(const void* a, const void* b)
{
	ContigGraphClassification* c1 = (ContigGraphClassification*) a;
	ContigGraphClassification* c2 = (ContigGraphClassification*) b;

	if (c1->flag & (CONTIG_UNCLASSIFIED | CONTIG_DISCARD | CONTIG_NORELATION | CONTIG_UNIQUE))
		return 1;

	if (c2->flag & (CONTIG_UNCLASSIFIED | CONTIG_DISCARD | CONTIG_NORELATION | CONTIG_UNIQUE))
		return -1;

	return (c1->coveredIntervals - c2->coveredIntervals);
}

void rawClassification(AnalyzeContext *actx)
{
	/*
	 int i, j, k;
	 int end1_hap, end2_hap; // ends from touring
	 int path_hap;

	 Contig *contig_hap;

	 bit *contigVisited = ba_new(actx->numContigs);

	 int idx = 0;
	 int curStart = 0;
	 char *InFastaName;
	 int path_other, end1_other, end2_other;

	 while (curStart < actx->numContigs)
	 {
	 int process = 0;
	 int curNum = getNumberOfSequencesFromFastaFileName(actx, curStart);

	 int compoundNumAllBases = 0;
	 int compoundNumHapContigs = 0;
	 int compoundNumHapBases = 0;
	 int compoundNumALTcontigs = 0;
	 int compoundNumALTbases = 0;
	 int compoundNumREPcontigs = 0;
	 int compoundNumREPbases = 0;
	 int compoundNumWEIRDcontigs = 0;
	 int compoundNumWEIRDbases = 0;

	 printf("#--------- Result of Connected Compound %5d ---------\n", idx);
	 while (1)
	 {
	 contig_hap = NULL;
	 for (i = curStart; i < curStart + curNum; i++)
	 {
	 if (actx->contigs[i].fClass == HAPLOID && ba_value(contigVisited, i) == FALSE)
	 {
	 contig_hap = actx->contigs + i;
	 path_hap = getPathID(actx, contig_hap->idx);
	 getContigsEndIDs(actx, contig_hap->idx, &end1_hap, &end2_hap);
	 break;
	 }
	 }

	 if (contig_hap != NULL)
	 {
	 ba_assign(contigVisited, i, TRUE);
	 process++;

	 int numALTcontigs = 0;
	 int numALTbases = 0;
	 int numREPcontigs = 0;
	 int numREPbases = 0;
	 int numWEIRDcontigs = 0;
	 int numWEIRDbases = 0;

	 // sort overlap groups and contig classification according to begin position in current haploid contig

	 qsort(contig_hap->ovlGrps, contig_hap->numOvlGrps, sizeof(OverlapGroup*), cmpOVLGroupByPos);
	 qsort(contig_hap->gClassific, contig_hap->numGClassific, sizeof(ContigGraphClassification), cmpGraphClassByPos);

	 for (j = 0, k = 0; j < contig_hap->numOvlGrps && k < contig_hap->numGClassific;)
	 {
	 OverlapGroup *og = contig_hap->ovlGrps[j];
	 ContigGraphClassification *cgc = contig_hap->gClassific + k;

	 assert(og->bread >= curStart && og->bread <= curStart + curNum);
	 assert(cgc->correspID >= curStart && cgc->correspID <= curStart + curNum);

	 if ((og->flag & OVLGRP_DISCARD) || (cgc->flag & (CONTIG_NORELATION | CONTIG_DISCARD | CONTIG_UNCLASSIFIED | CONTIG_UNIQUE)))
	 break;

	 if (og->bread == cgc->correspID)
	 {
	 assert(intersect(og->first_abpos, og->last_aepos, cgc->bpos, cgc->epos) > 0);

	 Contig * contig_other = actx->contigs + og->bread;

	 if (ba_value(contigVisited, contig_other->idx) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->idx);
	 path_other = getPathID(actx, contig_other->idx);
	 getContigsEndIDs(actx, contig_other->idx, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->idx);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_other, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->len);
	 // position in contig --> only meaningful in contained contigs
	 if (og->flag & OVLGRP_COMP)
	 fprintf(stdout, "\t%8d\t%8d", og->first_bepos, og->last_bbpos);
	 else
	 fprintf(stdout, "\t%8d\t%8d", og->first_bbpos, og->last_bepos);
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", og->first_abpos, og->last_aepos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);

	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->idx, TRUE);
	 process++;

	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->len;
	 break;
	 }
	 }
	 else
	 {
	 printf("found ovl group and gclass between contigs %d and %d more then once!!\n", contig_hap->idx, contig_other->idx);
	 }
	 j++;
	 k++;
	 continue;
	 }
	 if (og->first_abpos < cgc->bpos)
	 {
	 Contig * contig_other = actx->contigs + og->bread;

	 if (ba_value(contigVisited, contig_other->idx) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->idx);
	 path_other = getPathID(actx, contig_other->idx);
	 getContigsEndIDs(actx, contig_other->idx, &end1_other, &end2_other);

	 //
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->idx);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->len);
	 // position in contig --> only meaningful in contained contigs
	 if (og->flag & OVLGRP_COMP)
	 fprintf(stdout, "\t%8d\t%8d", og->first_bepos, og->last_bbpos);
	 else
	 fprintf(stdout, "\t%8d\t%8d", og->first_bbpos, og->last_bepos);
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", og->first_abpos, og->last_aepos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);

	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->idx, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->len;
	 break;
	 }
	 }
	 else
	 {
	 printf("found ovl group and gclass between contigs %d and %d more then once!!\n", contig_hap->idx, contig_other->idx);
	 }
	 j++;

	 continue;
	 }
	 else
	 {
	 Contig * contig_other = actx->contigs + cgc->correspID;

	 if (ba_value(contigVisited, contig_other->idx) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->idx);
	 path_other = getPathID(actx, contig_other->idx);
	 getContigsEndIDs(actx, contig_other->idx, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->idx);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", cgc->bpos, cgc->epos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);
	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->idx, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->len;
	 break;
	 }

	 }
	 else
	 {
	 printf("found ovl group and gclass between contigs %d and %d more then once!!\n", contig_hap->idx, contig_other->idx);
	 }
	 k++;
	 continue;
	 }
	 }
	 // check if there are remaining bubbles
	 for (; j < contig_hap->numOvlGrps; j++)
	 {
	 OverlapGroup *og = contig_hap->ovlGrps[j];

	 assert(og->bread >= curStart && og->bread <= curStart + curNum);
	 if (og->flag & OVLGRP_DISCARD)
	 break;
	 Contig * contig_other = actx->contigs + og->bread;

	 if (ba_value(contigVisited, contig_other->idx) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->idx);
	 path_other = getPathID(actx, contig_other->idx);
	 getContigsEndIDs(actx, contig_other->idx, &end1_other, &end2_other);

	 //
	 fprintf(stdout, "  %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->idx);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->len);
	 // position in contig --> only meaningful in contained contigs
	 if (og->flag & OVLGRP_COMP)
	 fprintf(stdout, "\t%8d\t%8d", og->first_bepos, og->last_bbpos);
	 else
	 fprintf(stdout, "\t%8d\t%8d", og->first_bbpos, og->last_bepos);
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", og->first_abpos, og->last_aepos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);
	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->idx, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->len;
	 break;
	 }
	 }
	 }
	 for (; k < contig_hap->numGClassific; k++)
	 {
	 ContigGraphClassification *cgc = contig_hap->gClassific + k;
	 assert(cgc->correspID >= curStart && cgc->correspID <= curStart + curNum);

	 if ((cgc->flag & (CONTIG_NORELATION | CONTIG_DISCARD | CONTIG_UNCLASSIFIED | CONTIG_UNIQUE)))
	 break;

	 Contig * contig_other = actx->contigs + cgc->correspID;

	 if (ba_value(contigVisited, contig_other->idx) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->idx);
	 path_other = getPathID(actx, contig_other->idx);
	 getContigsEndIDs(actx, contig_other->idx, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "   %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->idx);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t%8d\t%8d", cgc->bpos, cgc->epos);
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);
	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->idx, TRUE);
	 process++;
	 switch (contig_other->fAltClass)
	 {
	 case BUBBLE_ALT:
	 case SPUR_ALT:
	 numALTcontigs++;
	 numALTbases += contig_other->len;
	 break;
	 case BUBBLE_REPEAT:
	 case BUBBLE_HUGEDIF:
	 case BUBBLE_LOWCOVEREDREPEAT:
	 numREPcontigs++;
	 numREPbases += contig_other->len;
	 break;
	 default:
	 numWEIRDcontigs++;
	 numWEIRDbases += contig_other->len;
	 break;
	 }
	 }
	 }

	 InFastaName = getFastaFileNameFromDB(actx, contig_hap->idx);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "%c", getClassification(contig_hap->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_hap);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_hap->idx);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_hap->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_hap->len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_hap->repBases * 100.0 / contig_hap->len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_hap, end2_hap);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_hap->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_hap->fAltClass);
	 fprintf(stdout, "\n");
	 // haploid contig summary: print number of bubbles/spurs and percentage of bases according to contig length
	 fprintf(stdout, "# HAPLOID-Summary\n");
	 fprintf(stdout, "# HAPcontig  \tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", 1, contig_hap->len, 100.00);
	 fprintf(stdout, "# ALTcontig  \tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", numALTcontigs, numALTbases, numALTbases * 100.0 / contig_hap->len);
	 fprintf(stdout, "# REPcontig  \tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", numREPcontigs, numREPbases, numREPbases * 100.0 / contig_hap->len);
	 fprintf(stdout, "# WEIRDcontig\tBp\t%%Bp\t-\t%3d\t%8d\t%5.2f\n", numWEIRDcontigs, numWEIRDbases, numWEIRDbases * 100.0 / contig_hap->len);

	 // log splits
	 printf("splits: \n");
	 for (k = 0; k < contig_hap->numSplits; k++)
	 {
	 printSplitEvent(contig_hap, contig_hap->split + k);
	 }

	 printf("#-----------------------------------------------------\n");

	 // update compound numbers
	 compoundNumAllBases += (numALTbases + numREPbases + numWEIRDbases);
	 compoundNumALTcontigs += numALTcontigs;
	 compoundNumALTbases += numALTbases;
	 compoundNumREPcontigs += numREPcontigs;
	 compoundNumREPbases += numREPbases;
	 compoundNumWEIRDcontigs += numWEIRDcontigs;
	 compoundNumWEIRDbases += numWEIRDbases;
	 compoundNumHapContigs++;
	 compoundNumHapBases += contig_hap->len;
	 continue;
	 }
	 else // print out remaining contigs
	 {
	 if (process < curNum)
	 {

	 for (i = curStart; i < curStart + curNum; i++)
	 {
	 Contig * contig_other = actx->contigs + i;

	 if (ba_value(contigVisited, contig_other->idx) == FALSE)
	 {
	 InFastaName = getFastaFileNameFromDB(actx, contig_other->idx);
	 path_other = getPathID(actx, contig_other->idx);
	 getContigsEndIDs(actx, contig_other->idx, &end1_other, &end2_other);

	 // write + for haploid seqeunce, - otherwise
	 fprintf(stdout, "   %c", getClassification(contig_other->fClass));
	 // input fasta name
	 fprintf(stdout, "\t%s.fa", InFastaName);
	 // write path idx
	 fprintf(stdout, "\t%d", path_other);
	 // write db idx
	 fprintf(stdout, "\t%d", contig_other->idx);
	 // output fasta name: in + path + split + classifier
	 fprintf(stdout, "\t%s_%d_%d_%c.fa", InFastaName, path_hap, 0, getClassification(contig_other->fClass));
	 // length
	 fprintf(stdout, "\t%8d", contig_other->len);
	 // percent repeat
	 fprintf(stdout, "\t%5.2f", contig_other->repBases * 100.0 / contig_other->len);
	 // position in contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // position in corresponding contig --> only meaningful in contained contigs
	 fprintf(stdout, "\t      na\t      na");
	 // do contig tips link with other contig ids?
	 fprintf(stdout, "\t%d,%d", end1_other, end2_other);
	 fprintf(stdout, "\t");
	 // print contig classification
	 printFinalContigClassification(stdout, contig_other->fClass);
	 fprintf(stdout, "\t");
	 // print alt classification, only relevant if contig is somehow contained in another one
	 printFinalAltContigClassification(stdout, contig_other->fAltClass);

	 fprintf(stdout, "\n");

	 ba_assign(contigVisited, contig_other->idx, TRUE);
	 process++;

	 // update compound numbers
	 compoundNumAllBases += contig_other->len;

	 }
	 }
	 //printf("#-----------------------------------------------------\n");x
	 }
	 break;
	 }
	 }
	 printf("#--------- Summary of Connected Compound %5d ---------\n", idx);
	 fprintf(stdout, "# HAPcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumHapContigs, compoundNumHapBases, compoundNumHapBases * 100.0 / compoundNumAllBases);
	 fprintf(stdout, "# ALTcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumALTcontigs, compoundNumALTbases, compoundNumALTbases * 100.0 / compoundNumAllBases);
	 fprintf(stdout, "# REPcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumREPcontigs, compoundNumREPbases, compoundNumREPbases * 100.0 / compoundNumAllBases);
	 fprintf(stdout, "# WEIRDcontig\t%%Bp\t-\t%d\t%d\t%5.2f\n", compoundNumWEIRDcontigs, compoundNumWEIRDbases,
	 compoundNumWEIRDbases * 100.0 / compoundNumAllBases);

	 curStart = curStart + curNum + 1;
	 idx++;
	 break;
	 }*/
}

void classify(AnalyzeContext *actx)
{
	int i, j, k;
	int end1, end2; // ends from touring
	int path;

	for (i = 0; i < actx->numContigs; i++)
	{
		Contig *contig = actx->contigs + i;
		path = getPathID(actx, contig->idx);
		getContigsEndIDs(actx, contig->idx, &end1, &end2);

		contig->flag |= CONTIG_UNCLASSIFIED;

//		if ((contig->ovlGrpsClassificFlag & CONTIG_IS_CONTAINED))
//		{
//			// 1. check if its contained in multiple contigs?
//			OverlapGroup *ovlgr;
//			int contCount = 0;
//			for (j = 0; j < contig->numOvlGrps; j++)
//			{
//				if (contig->ovlGrps[j]->flag & CONTIG_IS_CONTAINED)
//				{
//					contCount++;
//					printf("%d contained %d vs %d", contCount, contig->ovlGrps[j]->aread, contig->ovlGrps[j]->coveredBasesInA);
//				}
//			}
//
//			assert(contCount > 0);
//
//			if (contCount > 1)
//			{
//				contig->flag &= ~( CONTIG_UNCLASSIFIED);
//				contig->flag |= CONTIG_IS_CONTAINED;
//				contig->class |= (CONTIG_CLASS_REPEAT | CONTIG_CLASS_WEIRD);
//				continue;
//			}
//
//			//	CONTINUE WITH FINDING same contained graph classification
//		}
//
//		if ((contig->ovlGrpsClassificFlag & CONTIG_IS_CONTAINED) && (contig->gClassificFlag & CONTIG_IS_CONTAINED))
//		{
//			OverlapGroup *ovlgr;
//			ContigGraphClassification *gclass;
//			for (j = 0; j < contig->numOvlGrps; j++)
//			{
//				ovlgr = contig->ovlGrps[j];
//
//			}
//
//			int match = 0;
//			for (j = 0; j < contig->numOvlGrps; j++)
//			{
//				ovlgr = contig->ovlGrps[j];
//
//				if ((ovlgr->flag & OVLGRP_AREAD_IS_CONTAINED) && DB_READ_LEN(actx->corContigDB, ovlgr->aread) < DB_READ_LEN(actx->corContigDB, ovlgr->bread))
//				{
//					for (k = 0; k < contig->numGClassific; k++)
//					{
//						gclass = contig->gClassific + k;
//
//						if ((gclass->flag & CONTIG_IS_CONTAINED) && contig->idx == 145)
//							printf("Contig 145 is contained in %d,  check by graphClassific\n", gclass->correspID);
//
//						if ((gclass->flag & CONTIG_IS_CONTAINED) && ovlgr->bread == gclass->correspID)
//							break;
//					}
//
//					if (k < contig->numGClassific) // found same classification
//					{
//
//						assert(intersect(ovlgr->first_abpos, ovlgr->last_aepos, gclass->coveredIntervals[0], gclass->coveredIntervals[gclass->numCoveredIntervals*3-2]) > 0);
//
//						// thats the clear part
//						contig->class |= CONTIG_CLASS_ALT;
//
//						if (end1 == end2)
//						{
//							contig->class |= CONTIG_CLASS_ALT_BUBBLE;
//							if (end1 != ovlgr->bread)                      // that should never happen !!!!!
//								contig->class |= CONTIG_CLASS_WEIRD;
//						}
//						else
//						{
//							if ((end1 == ovlgr->bread && end2 < 0) || (end2 == ovlgr->bread && end1 < 0))
//								contig->class = CONTIG_CLASS_ALT_SPUR;
//							else
//								contig->class |= CONTIG_CLASS_WEIRD;
//						}
//						match = 1;
//						break;
//					}
//				}
//			}

//			// todo how often does this happen, how to make a reliable classification ?
//			if (!match)
//			{
//				printf("Contig %d --> match != 1\n", contig->idx);
//				fflush(stdout);
//			}
//			assert(match == 1);
//			if (!match)    // classified contig id between gclass and ovlgrp differs !!!!
//			{
//				contig->fClass = CONTAINED_BY_BOTHSTRANGE;
//				if (end1 == end2)
//					contig->fAltClass = BUBBLE_ALT;
//				else
//					contig->fAltClass = SPUR_ALT;
//			}
//		}
//		else if (contig->flag & CONTIG_IS_CONTAINED)
//		{
//			// those are the fully repeat-induced "alternative" contigs check those guys
//			contig->class |= CONTIG_CLASS_REPEAT;
//		}
//		else if (contig->gClassificFlag & CONTIG_IS_CONTAINED)
//		{
//			// thats strange, should not happen !!! Actually it happens :) for spurs at contig tips, or ALT contigs with huge indels compared to corrersponding haploid contig
//			contig->class |= CONTIG_CLASS_ALT;
//			if (end1 < 0 || end2 < 0)
//			{
//				contig->class |= CONTIG_CLASS_ALT_SPUR;
//			}
//			else if (end1 == end2) // todo sanity check if a proper overlap group exists, that fails only due to some thresholds e.g. huge insert
//			{
//				contig->class |= CONTIG_CLASS_ALT_HUGEDIFF;
//			}
//			else
//			{
//				contig->class |= CONTIG_CLASS_WEIRD;
//			}
//		}
//		else
//		{
//			int exclude = 0;
//			if (contig->flag & CONTIG_AL50PCTCOVERED)
//			{
//				// check which one is contained
//				for (j = 0; j < actx->numContigs; j++)
//				{
//					if ((int) contig->pctCoveredBasesInOtherContigs[j] > 70 && DB_READ_LEN(actx->corContigDB, contig->idx) < DB_READ_LEN(actx->corContigDB, j))
//					{
//						exclude = 1;
//						break;
//					}
//				}
//			}
//			if (exclude)
//			{
//				contig->class |= (CONTIG_CLASS_ALT | CONTIG_CLASS_REPEAT | CONTIG_CLASS_WEIRD);
//			}
//			else
//			{
//				contig->class |= CONTIG_CLASS_HAPLOID;
//			}
//		}
	}
}

int checkJunctionCoverage(AnalyzeContext *actx, Contig *contig, int read) // returns true if everything is ok, false otherwise
{
	int i, j;

	int result = 1;

// todo pass variable, set by user ?
	int NUM_NODES_TOCHECK = 3; // number or reads to check before and after duplicated reads

	printf("checkJunctionCoverage for read %d and contig %d\n", read, contig->idx);

	j = 0;
	while (j < contig->numcReads)
	{
		if (abs(contig->cReads[j].patchedID) == read)
			break;

		j++;
	}

	assert(j < contig->numcReads);

// try to start NUM_NODES_TOCHECK before read
	int from, to;
	from = to = j;
	int k = NUM_NODES_TOCHECK;
	int k_len_before = 0;

	while (k > 0)
	{
		if (from - 1 > 0)
			from--;
		if (to + 1 < contig->numcReads)
			to++;

		k--;
	}

	for (k = from; k <= to; k++)
	{
		ContigRead *cread = contig->cReads + k;

		if (cread->nComReadsWithNextContigRead == 0 || cread->lowCov)
			result = 0;

		printf("CP: read %d, # in %d, # out %d, #comRead %d, cPos [%d, %d], rPos [%d, %d], avgCov %.2f, covDrop %d\n", cread->patchedID, cread->nInReads,
				cread->nOutReads, cread->nComReadsWithNextContigRead, cread->patchedContigPosBeg, cread->patchedContigPosEnd, cread->ovlReads->beg,
				cread->ovlReads->end, cread->avgCov, cread->lowCov);

	}
	return result;
}

static void usage()
{
	printf(" [-v] [-clxse <int>] [-d <dir>] [-rt <Track>] [-o <file>] -C <contigDB> <ContigOverlaps> -F <fixedReadDB> <fixedReadOverlaps> -D <correctedReadDB> \n");
	printf("options: -v         ... verbose\n");
	printf("         -x         ... min read length (default: 0)\n");
	printf("         -l         ... min alignment length (default: 1000)\n");
	printf("         -c         ... expected coverage\n");
	printf("         -s         ... maximum spur/bubble length (default: %d)\n", DEF_SPUR_LEN);
	printf("         -e         ... maximum contig end length (tips) to consider false joins (default: %d)\n", DEF_TIP_LEN);
	printf("         -r         ... repeats track for contigs (default: repeats)\n");
	printf("         -t         ... trim track for fixed reads (default: none)\n");
	printf("         -C DB OVL  ... contig database and contig overlaps, required to classify contigs\n");
	printf("         -F DB OVL  ... fixed read database and fixed read overlaps, required to extract all b-reads for blasr mapping\n");
	printf("         -D DB 		... corrected read database\n");
	printf("         -d         ... write classified contigs ands stats file into -d directoryName (default: cwd)\n");
	printf("         -o         ... write out filtered chained overlaps (default: cwd)\n");
}

int main(int argc, char* argv[])
{
	HITS_DB correctedContigDB;
	HITS_DB patchedReadDB;
	HITS_DB correctedReadDB;

	FILE* patchedReadLAS;
	FILE* correctedContigLAS;

	PassContext* contig_pctx;
	PassContext* patched_pctx;
	AnalyzeContext actx;

	char *patchedReadTrimTrack = NULL;
	char *corContigRepeatsTrack = DEF_REPEATS_TRACK;

	bzero(&actx, sizeof(AnalyzeContext));
	opterr = 0;

	char *filteredContigOvlsName = NULL;
	FILE *filteredContigOvlsFile = NULL;

	actx.VERBOSE = 0;
	actx.min_rlen = 0;
	actx.min_olen = 1000;
	actx.SPUR_LEN = DEF_SPUR_LEN;
	actx.TIP_LEN = DEF_TIP_LEN;

	int c;
	while ((c = getopt(argc, argv, "vx:l:c:d:t:r:C:F:s:e:D:o:")) != -1)
	{
		switch (c)
		{
		case 'o':
			filteredContigOvlsName = optarg;
			break;
		case 'v':
			actx.VERBOSE++;
			break;
		case 'x':
			actx.min_rlen = atoi(optarg);
			break;
		case 'l':
			actx.min_olen = atoi(optarg);
			break;
		case 'e':
			actx.TIP_LEN = atoi(optarg);
			break;
		case 's':
			actx.SPUR_LEN = atoi(optarg);
			break;
		case 'c':
		{
			actx.exp_cov = atoi(optarg);
			if (actx.exp_cov <= 0)
			{
				fprintf(stderr, "[ERROR] - analyzeContigs: expected coverage must be positive (>= 1)\n");
				usage();
				exit(1);
			}
		}
			break;
		case 'd':
			actx.outDir = optarg;
			break;
		case 't':
			patchedReadTrimTrack = optarg;
			break;
		case 'r':
			corContigRepeatsTrack = optarg;
			break;
		case 'C':
		{
			actx.corContigDBName = optarg;
			if (optind < argc && *argv[optind] != '-')
			{
				actx.corContigLASName = argv[optind];
				optind++;
			}
			else
			{
				fprintf(stderr, "\n-C option require TWO arguments <contigDB> <ContigOverlaps>\n\n");
				usage();
				exit(1);
			}
			break;
		}
		case 'F':
		{
			actx.patchedReadDBName = optarg;
			if (optind < argc && *argv[optind] != '-')
			{
				actx.patchedReadLASName = argv[optind];
				optind++;
			}
			else
			{
				fprintf(stderr, "\n-F option require TWO arguments <fixedReadDB> <fixedReadOverlaps>\n\n");
				usage();
				exit(1);
			}
			break;
		}
		case 'D':
			actx.corReadDBName = optarg;
			break;
		default:
			fprintf(stderr, "Unknown option: %s\n", argv[optind - 1]);
			usage();
			exit(1);
		}
	}

	if (actx.corReadDBName == NULL)
	{
		fprintf(stderr, "[ERROR] - corrected read database is required!\n");
		exit(1);
	}

	if (actx.corContigDBName == NULL && actx.corContigLASName == NULL)
	{
		fprintf(stderr, "[ERROR] - At least a contig dabase and the corresponding overlaps are required to analyze contigs!\n");
		exit(1);
	}

	if ((correctedContigLAS = fopen(actx.corContigLASName, "r")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.corContigLASName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Open database: %s\n", actx.corReadDBName);

	if (Open_DB(actx.corReadDBName, &correctedReadDB))
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.corReadDBName);
		exit(1);
	}

	if (actx.VERBOSE)
		printf("Open database: %s\n", actx.corContigDBName);

	if (Open_DB(actx.corContigDBName, &correctedContigDB))
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.corContigDBName);
		exit(1);
	}

	if ((actx.corContigRepeats_track = track_load(&correctedContigDB, corContigRepeatsTrack)) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'. Create a repeat track with LArepeat first!\n", corContigRepeatsTrack,
				actx.corContigDBName);
		exit(1);
	}

	if ((actx.corContigCorrectReadsAndPos_track = track_load(&correctedContigDB, "creads")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "creads", actx.corContigDBName);
		exit(1);
	}

	if ((actx.corContigRawReads_track = track_load(&correctedContigDB, "rreads")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "rreads", actx.corContigDBName);
		exit(1);
	}

	if ((actx.corContigPatchReadsAndPos_track = track_load(&correctedContigDB, "preads")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "preads", actx.corContigDBName);
		exit(1);
	}

	if ((patchedReadLAS = fopen(actx.patchedReadLASName, "r")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open '%s'\n", actx.patchedReadLASName);
		exit(1);
	}

///todo make this optional
	if (actx.VERBOSE)
		printf("Open database: %s\n", actx.patchedReadDBName);

	if (Open_DB(actx.patchedReadDBName, &patchedReadDB))
	{
		fprintf(stderr, "could not open '%s'\n", actx.patchedReadDBName);
		exit(1);
	}
	if (patchedReadTrimTrack != NULL)
	{
		if (actx.VERBOSE)
			printf("Load track %s for db %s\n", patchedReadTrimTrack, actx.patchedReadDBName);
		if ((actx.patchedReadTrim_track = track_load(&patchedReadDB, patchedReadTrimTrack)) == NULL)
		{
			fprintf(stderr, "could not open track '%s' for database '%s'\n", patchedReadTrimTrack, actx.patchedReadDBName);
			exit(1);
		}
	}
	else if ((actx.patchedReadTrim_track = track_load(&patchedReadDB, TRACK_TRIM)) == NULL)
	{
		fprintf(stderr, "could not open track '%s' for database '%s'\n",
		TRACK_TRIM, actx.patchedReadDBName);
		exit(1);
	}

	if ((actx.patchedReadSource_track = track_load(&patchedReadDB, "source")) == NULL)
	{
		fprintf(stderr, "[ERROR] - could not open track '%s' for database '%s'\n", "source", actx.patchedReadDBName);
		exit(1);
	}

// init

	actx.corContigDB = &correctedContigDB;
	actx.patchedReadDB = &patchedReadDB;
	actx.corReadDB = &correctedReadDB;

	char *cwd = malloc(10);
	if (actx.outDir == NULL)
	{
		sprintf(cwd, ".");
		actx.outDir = cwd;
	}

	{
		if (createOutDir(actx.outDir))
			exit(1);

		char *out = malloc(strlen(actx.outDir) + 30);
		sprintf(out, "%s/raw", actx.outDir);
		if (createOutDir(out))
			exit(1);

		sprintf(out, "%s/split", actx.outDir);
		if (createOutDir(out))
			exit(1);

		free(out);
	}

	if (filteredContigOvlsName)
	{
		char *out = malloc(strlen(actx.outDir) + strlen(filteredContigOvlsName) + 10);
		sprintf(out, "%s/%s", actx.outDir, filteredContigOvlsName);

		if((filteredContigOvlsFile = fopen(out, "w")) == NULL)
		{
			fprintf(stderr, "could not open %s\n", out);
			exit(1);
		}
		free(out);
	}


// contig pass context
	{
		contig_pctx = pass_init(correctedContigLAS, filteredContigOvlsFile);

		contig_pctx->split_b = 0;
		contig_pctx->load_trace = 0;
		contig_pctx->unpack_trace = 0;
		contig_pctx->data = &actx;

		pre(contig_pctx, &actx);
	}
// fixed read pass context
	{
		patched_pctx = pass_init(patchedReadLAS, NULL);
		patched_pctx->split_b = 0;
		patched_pctx->load_trace = 1;
		patched_pctx->unpack_trace = 1;
		patched_pctx->data = &actx;

		pre(patched_pctx, &actx);
	}
// initialize contigs
	printf("STEP0a: initialize AnalyzeContext - START\n");
	initAnalyzeContext(&actx);
	printf("STEP0a: initialize AnalyzeContext - DONE\n");
// todo for now do all steps: input format should from MARVEL assembly run
// analyze overlaps of fixed reads for each contig ( coverage, containments, duplicate reads used in different contigs )
// STEP1
	if (1)
	{
		printf("START    ---   STEP1A: analyze overlapping reads from patched database, with contig reads\n");
		pass(patched_pctx, processFixedReadOverlap_handler1);
		printf("DONE     ---   STEP1A: analyze overlapping reads from patched database, with contig reads\n");
		printf("START    ---   STEP1B: -- analyzeContigOverlapGraph\n");
		analyzeContigOverlapGraph(&actx);
		printf("DONE     ---   STEP1B: -- analyzeContigOverlapGraph\n");
		printf("START    ---   STEP1C: -- classifyContigsByBReadsAndPath\n");
//		classifyContigsByBReadsAndPath(&actx);
		printf("DONE     ---   STEP1C: -- classifyContigsByBReadsAndPath\n");
	}

// analyze overlaps between contigs
// STEP2
	if (1)
	{
		printf("START    ---   STEP2a: analyze contig alignments\n");
		pass(contig_pctx, processContigOverlap_handler);
		printf("DONE     ---   STEP2a: analyze contig alignments\n");
		printf("START    ---   STEP2b: classify contig by contig alignments\n");
//		classifyContigsByOverlaps(&actx);
		printf("DONE     ---   STEP2b: classify contig by contig alignments\n");
	}

// final classification, based on STEP1 and STEP2
// STEP3
	if (1)
	{
		printf("START    ---   STEP3: refine contig classification (based on STEP1 and STEP2)");
//		classify(&actx);
		printf("DONE     ---   STEP3: refine contig classification (based on STEP1 and STEP2)");
	}
// check reads that occur multiple times in contigs, and set proper split positions if required
// STEP4
	{
		printf("STEP4: analyze reads that occur multiple times in contigs, and set proper split positions if required\n");
		// analyze number of shared reads between contigs, trim them back if possible /// split contigs ??
//		finalContigValidation(&actx);
	}
// create output: contigs, splitted contigs, stat and bed files with i.e. coverage histogram, repeat tracks, split coordinates, ...
// STEP5
	{
		printf("STEP5: create output\n");
//		rawClassification(&actx);

	}

	pass_free(contig_pctx);
	pass_free(patched_pctx);

	fflush(stdout);
	fflush(stderr);

	free(cwd);
	return 0;
}
