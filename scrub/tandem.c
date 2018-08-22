/************************************************************************************\
*                                                                                    *
* Copyright (c) 2014, Dr. Eugene W. Myers (EWM). All rights reserved.                *
*                                                                                    *
* Redistribution and use in source and binary forms, with or without modification,   *
* are permitted provided that the following conditions are met:                      *
*                                                                                    *
*  · Redistributions of source code must retain the above copyright notice, this     *
*    list of conditions and the following disclaimer.                                *
*                                                                                    *
*  · Redistributions in binary form must reproduce the above copyright notice, this  *
*    list of conditions and the following disclaimer in the documentation and/or     *
*    other materials provided with the distribution.                                 *
*                                                                                    *
*  · The name of EWM may not be used to endorse or promote products derived from     *
*    this software without specific prior written permission.                        *
*                                                                                    *
* THIS SOFTWARE IS PROVIDED BY EWM ”AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND       *
* FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL EWM BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS  *
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY      *
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     *
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN  *
* IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                      *
*                                                                                    *
* For any issues regarding this software and its use, contact EWM at:                *
*                                                                                    *
*   Eugene W. Myers Jr.                                                              *
*   Bautzner Str. 122e                                                               *
*   01099 Dresden                                                                    *
*   GERMANY                                                                          *
*   Email: gene.myers@gmail.com                                                      *
*                                                                                    *
\************************************************************************************/

/*******************************************************************************************
 *
 *  Adaption of the fast local alignment filter of "daligner" for finding alignments between
 *     a read and itself above the diagonal.  Satellites create a "ladder" of these self alignments
 *     that can be detected in order to recognize said repetitive sequences within the read.
 *
 *  Author :  Gene Myers
 *  First  :  March 2016
 *  Current:  March 27, 2016
 *
 ********************************************************************************************/

//  A complete threaded code for the filter

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include "db/DB.h"
#include "dalign/align.h"
#include "tandem.h"

#define THREAD    pthread_t

#define PANEL_SIZE     50000   //  Size to break up very long A-reads
#define PANEL_OVERLAP  10000   //  Overlap of A-panels

#define MATCH_CHUNK    100     //  Max initial number of hits between two reads
#define TRACE_CHUNK  20000     //  Max initial trace points in hits between two reads

#undef  TEST_LSORT
#undef  TEST_KSORT
#undef  TEST_PAIRS
#undef  TEST_CSORT
#define    HOW_MANY   3000   //  Print first HOW_MANY items for each of the TEST options above

#undef  TEST_GATHER
#undef  TEST_CONTAIN
#undef  SHOW_OVERLAP          //  Show the cartoon
#undef  SHOW_ALIGNMENT        //  Show the alignment
#define   ALIGN_WIDTH    80   //     Parameters for alignment
#define   ALIGN_INDENT   20
#define   ALIGN_BORDER   10

#ifdef SHOW_OVERLAP
#define NOTHREAD
#endif

#ifdef TEST_GATHER
#define NOTHREAD
#endif

#ifdef TEST_CONTAIN
#define NOTHREAD
#endif

typedef struct
  { uint64 p1;   //  The lower half
    uint64 p2;
  } Double;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__

typedef struct
  { uint64 code;
    int    rpos;
    int    read;
  } KmerPos;

#else

typedef struct
  { uint64 code;
    int    read;
    int    rpos;
  } KmerPos;

#endif

/*******************************************************************************************
 *
 *  PARAMETER SETUP
 *
 ********************************************************************************************/

static int Kmer;
static int Hitmin;
static int Binshift;

static int    Kshift;         //  2*Kmer
static uint64 Kmask;          //  4^Kmer-1

static int NTHREADS;    //  Must be a power of 2
static int NSHIFT;      //  log_2 NTHREADS

int Set_Filter_Params(int kmer, int binshift, int hitmin, int nthreads)
{ if (kmer <= 1)
    return (1);

  Kmer     = kmer;
  Binshift = binshift;
  Hitmin   = hitmin;

  Kshift = 2*Kmer;
  if (Kmer == 32)
    Kmask = 0xffffffffffffffffllu;
  else
    Kmask = (0x1llu << Kshift) - 1;

  NTHREADS = 1;
  NSHIFT   = 0;
  while (2*NTHREADS <= nthreads)
    { NTHREADS *= 2;
      NSHIFT   += 1;
    }

  return (0);
}


/*******************************************************************************************
 *
 *  LEXICOGRAPHIC SORT
 *
 ********************************************************************************************/

#define BMER      4
#define BSHIFT    8             //  = 2*BMER
#define BPOWR   256             //  = 2^BSHIFT
#define BMASK  0xffllu          //  = BPOWR-1

static uint64  QMASK;           //  = BMASK << NSHIFT
static int     LEX_shift;
static int64   LEX_zsize;
static int     LEX_last;
static int     LEX_next;
static Double *LEX_src;
static Double *LEX_trg;

typedef struct
  { int64  beg;
    int64  end;
    int64  tptr[BPOWR];
    int64 *sptr;
  } Lex_Arg;

static void *lex_thread(void *arg)
{ Lex_Arg    *data  = (Lex_Arg *) arg;
  int64      *sptr  = data->sptr;
  int64      *tptr  = data->tptr;
  int         shift = LEX_shift;   //  Must be a multiple of 8 in [0,120]
  int        qshift = (LEX_next - LEX_shift) - NSHIFT;
  int64       zsize = LEX_zsize;
  Double     *src   = LEX_src;
  Double     *trg   = LEX_trg;
  int64       i, n, x;
  uint64      c, b;

  n = data->end;
  if (shift >= 64)
    { shift -= 64;
      if (LEX_last)
        for (i = data->beg; i < n; i++)
          { c = src[i].p2;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p2;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((b >> qshift) & QMASK) + x/zsize] += 1;
          }
    }

  else if ( ! LEX_last && LEX_next >= 64)   //  && LEX_shift < 64

    { qshift = (LEX_next - 64) - NSHIFT;
      if (qshift < 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((src[i].p2 << NSHIFT) & QMASK) + x/zsize] += 1;
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((src[i].p2 >> qshift) & QMASK) + x/zsize] += 1;
          }
    }

  else // LEX_last || LEX_next < 64
    if (LEX_last)
      if (shift == 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            x = tptr[c&BMASK]++;
            trg[x] = src[i];
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
          }
    else
      if (shift == 0)
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            x = tptr[c&BMASK]++;
            trg[x] = src[i];
            sptr[((c >> qshift) & QMASK) + x/zsize] += 1;
          }
      else
        for (i = data->beg; i < n; i++)
          { c = src[i].p1;
            b = (c >> shift);
            x = tptr[b&BMASK]++;
            trg[x] = src[i];
            sptr[((b >> qshift) & QMASK) + x/zsize] += 1;
          }

  return (NULL);
}

static Double *lex_sort(int bytes[16], Double *src, Double *trg, Lex_Arg *parmx)
{ THREAD  threads[NTHREADS];

  int64   len, x, y;
  Double *xch;
  int     i, j, k, z;
  int     b, c, fb;

  len       = parmx[NTHREADS-1].end;
  LEX_zsize = (len-1)/NTHREADS + 1;
  LEX_src   = src;
  LEX_trg   = trg;
  QMASK     = (BMASK << NSHIFT);

  for (c = 0; c < 16; c++)
    if (bytes[c])
      break;
  fb = c;
  for (b = c; b < 16; b = c)
    { for (c = b+1; c < 16; c++)
        if (bytes[c])
          break;
      LEX_last  = (c >= 16);
      LEX_shift = (b << 3);
      LEX_next  = (c << 3);

      if (b == fb)
        { for (i = 0; i < NTHREADS; i++)
            for (z = 0; z < NTHREADS*BPOWR; z++)
              parmx[i].sptr[z] = 0;
        }
      else
        { x = 0;
          for (i = 0; i < NTHREADS; i++)
            { parmx[i].beg = x;
              x = LEX_zsize*(i+1);
              if (x > len)
                x = len;
              parmx[i].end = x;
              for (j = 0; j < BPOWR; j++)
                parmx[i].tptr[j] = 0;
            }
          parmx[NTHREADS-1].end = len;

          for (j = 0; j < BPOWR; j++)
            { k = (j << NSHIFT);
              for (z = 0; z < NTHREADS; z++)
                for (i = 0; i < NTHREADS; i++)
                  { parmx[i].tptr[j] += parmx[z].sptr[k+i];
                    parmx[z].sptr[k+i] = 0;
                  }
            }
        }

      x = 0;
      for (j = 0; j < BPOWR; j++)
        for (i = 0; i < NTHREADS; i++)
          { y = parmx[i].tptr[j];
            parmx[i].tptr[j] = x;
            x += y;
          }

      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,lex_thread,parmx+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      xch     = LEX_src;
      LEX_src = LEX_trg;
      LEX_trg = xch;

#ifdef TEST_LSORT
      printf("\nLSORT %d\n",LEX_shift);
      if (LEX_shift >= 64)
        { x = (1 << ((LEX_shift-64)+BSHIFT))-1;
          for (i = 0; i < len; i++)
            { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                     i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                     LEX_src[i].p1&0xffffffffll,LEX_src[i].p2&x);
              if (i > 0 && (LEX_src[i].p1 < LEX_src[i].p1 ||
                             (LEX_src[i].p1 == LEX_src[i].p1 &&
                             (LEX_src[i].p2 & x) < (LEX_src[i-1].p2 & x))))
                printf(" OO");
              printf("\n");
            }
        }
      else
        { x = (1 << (LEX_shift+BSHIFT))-1;
          for (i = 0; i < len; i++)
            { printf("%6d: %8llx %8llx %8llx %8llx : %4llx",
                     i,LEX_src[i].p2>>32,(LEX_src[i].p2)&0xffffffffll,LEX_src[i].p1>>32,
                     LEX_src[i].p1&0xffffffffll,LEX_src[i].p1&x);
              if (i > 0 && (LEX_src[i].p1 & x) < (LEX_src[i-1].p1 & x))
                printf(" OO");
              printf("\n");
            }
        }
#endif
    }

  return (LEX_src);
}


/*******************************************************************************************
 *
 *  INDEX BUILD
 *
 ********************************************************************************************/

static HITS_DB    *TA_block;
static KmerPos    *TA_list;

typedef struct
  { int    tnum;
    int64 *kptr;
  } Tuple_Arg;

static void *tuple_thread(void *arg)
{ Tuple_Arg  *data  = (Tuple_Arg *) arg;
  int         tnum  = data->tnum;
  int64      *kptr  = data->kptr;
  KmerPos    *list  = TA_list;
  int         i, m, n, x, p;
  uint64      c;
  char       *s;

  c  = TA_block->nreads;
  i  = (c * tnum) >> NSHIFT;
  n  = TA_block->reads[i].boff;
  s  = ((char *) (TA_block->bases)) + n;
  n -= Kmer*i;

  for (m = (c * (tnum+1)) >> NSHIFT; i < m; i++)
    { c = p = 0;
      for (x = 1; x < Kmer; x++)
        c = (c << 2) | s[p++];
      while ((x = s[p]) != 4)
        { c = ((c << 2) | x) & Kmask;
          list[n].read = i;
          list[n].rpos = ++p;
          list[n].code = c;
          n += 1;
          kptr[c & BMASK] += 1;
        }
      s += (p+1);
    }

  return (NULL);
}

static KmerPos *Sort_Kmers(HITS_DB *block, int *len, KmerPos **buffer)
{ THREAD    threads[NTHREADS];
  Tuple_Arg parmt[NTHREADS];
  Lex_Arg   parmx[NTHREADS];
  int       mersort[16];

  KmerPos  *src, *trg, *rez;
  int       kmers, nreads;
  int       i, j, x;

  for (i = 0; i < NTHREADS; i++)
    parmx[i].sptr = (int64 *) alloca(NTHREADS*BPOWR*sizeof(int64));

  for (i = 0; i < 16; i++)
    mersort[i] = 0;
  for (i = 0; i < Kshift; i += 8)
    mersort[i>>3] = 1;

  nreads = block->nreads;
  kmers  = block->reads[nreads].boff - Kmer * nreads;

  if (kmers <= 0)
    goto no_mers;

  if (((Kshift-1)/BSHIFT) & 0x1)
    { trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+1),"Allocating Sort_Kmers vectors");
      src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+1),"Allocating Sort_Kmers vectors");
    }
  else
    { src = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+1),"Allocating Sort_Kmers vectors");
      trg = (KmerPos *) Malloc(sizeof(KmerPos)*(kmers+1),"Allocating Sort_Kmers vectors");
    }
  if (src == NULL || trg == NULL)
    exit (1);

  if (VERBOSE)
    { printf("\n   Kmer count = ");
      Print_Number((int64) kmers,0,stdout);
      printf("\n   Using %.2fGb of space\n",(1. * kmers) / 33554432);
      fflush(stdout);
    }

  TA_block = block;
  TA_list  = src;

  for (i = 0; i < NTHREADS; i++)
    { parmt[i].tnum = i;
      parmt[i].kptr = parmx[i].tptr;
      for (j = 0; j < BPOWR; j++)
        parmt[i].kptr[j] = 0;
    }

  for (i = 0; i < NTHREADS; i++)
    pthread_create(threads+i,NULL,tuple_thread,parmt+i);

  for (i = 0; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

  x = 0;
  for (i = 0; i < NTHREADS; i++)
    { parmx[i].beg = x;
      j = (int) ((((int64) nreads) * (i+1)) >> NSHIFT);
      parmx[i].end = x = block->reads[j].boff - j*Kmer;
    }

  rez = (KmerPos *) lex_sort(mersort,(Double *) src,(Double *) trg,parmx);

  if (rez[kmers-1].code == 0xffffffffffffffffllu)
    rez[kmers].code = 0;
  else
    rez[kmers].code = 0xffffffffffffffffllu;
  rez[kmers].read = nreads;

  if (src != rez)
    *buffer = src;
  else
    *buffer = trg;

#ifdef TEST_KSORT
  { int i;

    printf("\nKMER SORT:\n");
    for (i = 0; i < HOW_MANY && i < kmers; i++)
      { KmerPos *c = rez+i;
        printf(" %5d / %5d / %10lld\n",c->read,c->rpos,c->code);
      }
    fflush(stdout);
  }
#endif

  if (VERBOSE)
    { printf("   Index occupies %.2fGb\n",(1. * kmers) / 67108864);
      fflush(stdout);
    }

  if (kmers <= 0)
    { free(rez);
      goto no_mers;
    }

  *len = kmers;
  return (rez);

no_mers:
  *len = 0;
  return (NULL);
}


/*******************************************************************************************
 *
 *  FILTER MATCH
 *
 ********************************************************************************************/

  //  After the initial sort, all equal K-tuples are contiguous, and for a given K-tuple read,rpos
  //    is in sorted order because the sort is stable.  For a given read, record the distance between
  //    consecutive positions with the same K-mer in the .code field.  Afterwords, the array is
  //    stably resorted on read,rpos so that for each read one has effectively a "linked list" of
  //    positions with equal K-mers.

static KmerPos  *MG_alist;

typedef struct
  { int    abeg, aend;
    int64 *kptr;
  } Merge_Arg;

static void *count_thread(void *arg)
{ Merge_Arg  *data  = (Merge_Arg *) arg;
  int64      *kptr  = data->kptr;
  KmerPos    *asort = MG_alist;
  int         aend  = data->aend;

  uint64 ca, da;
  int    ia;
  int    ar, nr;
  int    ap, np;

  ia = data->abeg;
  da = asort[ia].code;
  while (ia < aend)
    { ca = da;
      ar = asort[ia].read;
      ap = asort[ia].rpos;
      kptr[ap & BMASK] += 1;
      while ((da = asort[++ia].code) == ca)
        { np = asort[ia].rpos;
          if ((nr = asort[ia].read) == ar)
            asort[ia].code = np - ap;
          else
            { asort[ia].code = 0;
              ar = nr;
            }
          ap = np;
          kptr[ap & BMASK] += 1;
        }
      asort[ia].code = 0;
    }

  return (NULL);
}

  //  Report threads: using the linked lists of consereved K-mers find likely seeds and then check
  //    for alignments as per "daligner".

static HITS_DB    *MR_ablock;
static Align_Spec *MR_spec;
static int         MR_tspace;

typedef struct
  { uint64   max;
    uint64   top;
    uint16  *trace;
  } Trace_Buffer;

  //  Determine if the minimum B-distance between the overlapping trace points of jpath and kpath
  //    Return this minimum and also the A-coordinate of a trace point where it is acheived.

static int Entwine(Path *jpath, Path *kpath, Trace_Buffer *tbuf, int *where)
{ int   ac, b2, y2, ae;
  int   i, j, k;
  int   num, den, min;
#ifdef SEE_ENTWINE
  int   strt, iflare, oflare;
#endif

  uint16 *ktrace = tbuf->trace + (uint64) (kpath->trace);
  uint16 *jtrace = tbuf->trace + (uint64) (jpath->trace);

  min   = 10000;
  num   = 0;
  den   = 0;

#ifdef SEE_ENTWINE
  strt = 1;
  printf("\n");
#endif

  y2 = jpath->bbpos;
  j  = jpath->abpos/MR_tspace;

  b2 = kpath->bbpos;
  k  = kpath->abpos/MR_tspace;

  if (jpath->abpos == kpath->abpos)
    { min = abs(y2-b2);
      if (min == 0)
        *where = kpath->abpos;
    }

  if (j < k)
    { ac = k*MR_tspace;

      j = 1 + 2*(k-j);
      k = 1;

      for (i = 1; i < j; i += 2)
        y2 += jtrace[i];
    }
  else
    { ac = j*MR_tspace;

      k = 1 + 2*(j-k);
      j = 1;

      for (i = 1; i < k; i += 2)
        b2 += ktrace[i];
    }

  ae = jpath->aepos;
  if (ae > kpath->aepos)
    ae = kpath->aepos;

  while (1)
    { ac += MR_tspace;
      if (ac >= ae)
        break;
      y2 += jtrace[j];
      b2 += ktrace[k];
      j += 2;
      k += 2;

#ifdef SEE_ENTWINE
      printf("   @ %5d : %5d %5d = %4d\n",ac,y2,b2,abs(b2-y2));
#endif

      i = abs(y2-b2);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = ac;
        }
      num += i;
      den += 1;
#ifdef SEE_ENTWINE
      if (strt)
        { strt   = 0;
          iflare = i;
        }
      oflare = i;
#endif
    }

  if (jpath->aepos == kpath->aepos)
    { i = abs(jpath->bepos-kpath->bepos);
      if (i <= min)
        { min = i;
          if (i == 0)
            *where = kpath->aepos;
        }
    }

#ifdef SEE_ENTWINE
  if (den == 0)
    printf("Nothing\n");
  else
    printf("MINIM = %d AVERAGE = %d  IFLARE = %d  OFLARE = %d\n",min,num/den,iflare,oflare);
#endif

  if (den == 0)
    return (-1);
  else
    return (min);
}

  //  Produce the concatentation of path1 and path2 where they are known to meet at
  //    the trace point with coordinate ap. Place this result in a big growing buffer,
  //    that gets reset when fusion is called with path1 = NULL

static void Fusion(Path *path1, int ap, Path *path2, Trace_Buffer *tbuf)
{ int     k, k1, k2;
  int     len, diff;
  uint16 *trace;

  k1 = 2 * ((ap/MR_tspace) - (path1->abpos/MR_tspace));
  k2 = 2 * ((ap/MR_tspace) - (path2->abpos/MR_tspace));

  len = k1+(path2->tlen-k2);

  if (tbuf->top + len >= tbuf->max)
    { tbuf->max = 1.2*(tbuf->top+len) + 1000;
      tbuf->trace = (uint16 *) Realloc(tbuf->trace,sizeof(uint16)*tbuf->max,"Allocating paths");
      if (tbuf->trace == NULL)
        exit (1);
    }

  trace = tbuf->trace + tbuf->top;
  tbuf->top += len;

  diff = 0;
  len  = 0;
  if (k1 > 0)
    { uint16 *t = tbuf->trace + (uint64) (path1->trace);
      for (k = 0; k < k1; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }
  if (k2 < path2->tlen)
    { uint16 *t = tbuf->trace + (uint64) (path2->trace);
      for (k = k2; k < path2->tlen; k += 2)
        { trace[len++] = t[k];
          trace[len++] = t[k+1];
          diff += t[k];
        }
    }

  path1->aepos = path2->aepos;
  path1->bepos = path2->bepos;
  path1->diffs = diff;
  path1->trace = (void *) (trace - tbuf->trace);
  path1->tlen  = len;
}

  //  Given all the LA's for a given read in amatch[0..novls-1], merge any overlapping LA's and
  //    remove any redundant ones.

static int Handle_Redundancies(Path *amatch, int novls, Trace_Buffer *tbuf)
{ Path *jpath, *kpath;
  int   j, k, no;
  int   dist, awhen = 0;

#ifdef TEST_CONTAIN
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].abpos,amatch[j].aepos,
                                              amatch[j].bbpos,amatch[j].bepos);
#endif

  for (j = 1; j < novls; j++)
    { jpath = amatch+j;
      for (k = j-1; k >= 0; k--)
        { kpath = amatch+k;

          if (kpath->abpos < 0)
            continue;

          if (jpath->abpos < kpath->abpos)

            { if (kpath->abpos <= jpath->aepos && kpath->bbpos <= jpath->bepos)
                { dist = Entwine(jpath,kpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->aepos > jpath->aepos)
                        { Fusion(jpath,awhen,kpath,tbuf);
#ifdef TEST_CONTAIN
                          printf("  Really 3");
#endif
                          k = j;
                        }
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! A %d %d\n",j,k);
#endif
                    }
                }
            }

          else // kpath->abpos <= jpath->abpos

            { if (jpath->abpos <= kpath->aepos && jpath->bbpos <= kpath->bepos)
                { dist = Entwine(kpath,jpath,tbuf,&awhen);
                  if (dist == 0)
                    { if (kpath->abpos == jpath->abpos)
                        { if (kpath->aepos > jpath->aepos)
                            *jpath = *kpath;
                        }
                      else if (jpath->aepos > kpath->aepos)
                        { Fusion(kpath,awhen,jpath,tbuf);
                          *jpath = *kpath;
#ifdef TEST_CONTAIN
                          printf("  Really 6");
#endif
                          k = j;
                        }
                      else
                        *jpath = *kpath;
                      kpath->abpos = -1;
#ifdef TEST_CONTAIN
                      printf("  Fuse! B %d %d\n",j,k);
#endif
                    }
                }
            }
        }
    }

  no = 0;
  for (j = 0; j < novls; j++)
    if (amatch[j].abpos >= 0)
      amatch[no++] = amatch[j];
  novls = no;

#ifdef TEST_CONTAIN
  for (j = 0; j < novls; j++)
    printf("  %3d: [%5d,%5d] x [%5d,%5d]\n",j,amatch[j].abpos,amatch[j].aepos,
                                              amatch[j].bbpos,amatch[j].bepos);
#endif

  return (novls);
}

  //  Determine the range of diagonal bins containing trace points of path

void Diagonal_Span(Path *path, int *mind, int *maxd)
{ uint16 *points;
  int     i, tlen;
  int     dd, low, hgh;

  points = path->trace;
  tlen   = path->tlen;

  dd = path->abpos - path->bbpos;
  low = hgh = dd;

  dd = path->aepos - path->bepos;
  if (dd < low)
    low = dd;
  else if (dd > hgh)
    hgh = dd;

  dd = (path->abpos/MR_tspace)*MR_tspace - path->bbpos;
  tlen -= 2;
  for (i = 1; i < tlen; i += 2)
    { dd += MR_tspace - points[i];
      if (dd < low)
        low = dd;
      else if (dd > hgh)
        hgh = dd;
    }

  *mind = (low >> Binshift)-1;
  *maxd = (hgh >> Binshift)+1;
}

typedef struct
  { int64       beg, end;
    int        *score;
    int        *lastp;
    int        *lasta;
    Work_Data  *work;
    FILE       *ofile;
    int64       nfilt;
    int64       ncheck;
    Overlap_IO_Buffer *iobuf;
  } Report_Arg;

static void *report_thread(void *arg)
{ Report_Arg  *data   = (Report_Arg *) arg;

  HITS_READ   *aread  = MR_ablock->reads;
  KmerPos     *asort  = MG_alist;

  int         *score  = data->score;
  int         *scorp  = data->score + 1;
  int         *scorm  = data->score - 1;
  int         *lastp  = data->lastp;
  int         *lasta  = data->lasta;
  int          maxdiag = ( MR_ablock->maxlen >> Binshift);
  int          mindiag = (-MR_ablock->maxlen >> Binshift);
  Overlap_IO_Buffer *obuf = data->iobuf;

//  FILE        *ofile  = data->ofile;

  char        *aseq   = (char *) (MR_ablock->bases);
  Work_Data   *work   = data->work;
  int          afirst = MR_ablock->ufirst;

  Overlap     _ovla, *ovla = &_ovla;
  Alignment   _align, *align = &_align;
  Path        *apath = &(ovla->path);
  int64        nfilt = 0;
  int64        ahits = 0;
  int          small, tbytes;

  int    novla;
  int    AOmax;
  Path  *amatch;

  Trace_Buffer _tbuf, *tbuf = &_tbuf;

  int      ar, aend;
  KmerPos *aoff;

  align->flags = ovla->flags = 0;
  align->path  = apath;

  if (MR_tspace <= TRACE_XOVR)
    { small  = 1;
      tbytes = sizeof(uint8);
    }
  else
    { small  = 0;
      tbytes = sizeof(uint16);
    }

  AOmax  = MATCH_CHUNK;
  amatch = Malloc(sizeof(Path)*AOmax,"Allocating match vector");

  tbuf->max   = 2*TRACE_CHUNK;
  tbuf->trace = Malloc(sizeof(short)*tbuf->max,"Allocating trace vector");

  if (amatch == NULL || tbuf->trace == NULL)
    exit (1);

//  fwrite(&ahits,sizeof(int64),1,ofile);
//  fwrite(&MR_tspace,sizeof(int),1,ofile);

#ifdef TEST_GATHER
  printf("\nNEW THREAD %5d(%5lld)-%5d(%5lld)\n",asort[data->beg].read,data->beg,
                                                asort[data->end-1].read,data->end);
  fflush(stdout);
#endif

  aoff = asort + (data->beg - Kmer);
  aend = asort[data->end-1].read;
  for (ar = asort[data->beg].read; ar <= aend; ar++)
    { int alen, amarkb, amarke;
      int apos, diag;
      int setaln;

#ifdef TEST_GATHER
      printf("Read %5d\n",ar);
      fflush(stdout);
#endif
      setaln = 1;
      novla  = 0;
      tbuf->top = 0;

      alen   = aread[ar].rlen;
      amarkb = Kmer;
      amarke = PANEL_SIZE;
      if (amarke >= alen)
        amarke = alen+1;
      while (1)
        {
          // Accumulate diagonal scores

          for (apos = amarkb; apos < amarke; apos++)
            { diag = aoff[apos].code;
              if (diag == 0) continue;
              diag >>= Binshift;
              if (apos - lastp[diag] >= Kmer)
                score[diag] += Kmer;
              else
                score[diag] += apos - lastp[diag];
              lastp[diag] = apos;
            }

          // Examine diagonal scores for hits to check out

          for (apos = amarkb; apos < amarke; apos++)
            { diag = aoff[apos].code;
              if (diag == 0) continue;
              diag >>= Binshift;
              if (apos > lasta[diag] &&
                  (score[diag] + scorp[diag] >= Hitmin || score[diag] + scorm[diag] >= Hitmin))
                { int bpos;

                  bpos = apos - aoff[apos].code;
                  if (setaln)
                    { setaln = 0;
                      align->aseq = align->bseq = aseq + aread[ar].boff;
                      align->alen = align->blen = alen;
                      ovla->aread = ovla->bread = ar + afirst;
                    }
#ifdef TEST_GATHER
                  if (scorm[diag] > scorp[diag])
                    printf("  %5d.. x %5d.. %5d (%3d)",
                           bpos,apos,apos-bpos,score[diag]+scorm[diag]);
                  else
                    printf("  %5d.. x %5d.. %5d (%3d)",
                           bpos,apos,apos-bpos,score[diag]+scorp[diag]);
                  fflush(stdout);
#endif

                  nfilt += 1;

                  Local_Alignment(align,work,MR_spec,apos-bpos,apos-bpos,apos+bpos,-1,-1);

                  { int low, hgh, ae;

                    Diagonal_Span(apath,&low,&hgh);
                    if (diag < low)
                      low = diag;
                    else if (diag > hgh)
                      hgh = diag;
                    ae = apath->aepos;
                    for (diag = low; diag <= hgh; diag++)
                      if (ae > lasta[diag])
                        lasta[diag] = ae;
#ifdef TEST_GATHER
                    printf(" %d - %d @ %d",low,hgh,apath->aepos);
#endif
                  }

                  if ((apath->aepos-apath->abpos) + (apath->bepos-apath->bbpos) >= MINOVER)
                    { if (novla >= AOmax)
                        { AOmax = 1.2*novla + MATCH_CHUNK;
                          amatch = Realloc(amatch,sizeof(Path)*AOmax,
                                           "Reallocating match vector");
                          if (amatch == NULL)
                            exit (1);
                        }
                      if (tbuf->top + apath->tlen > tbuf->max)
                        { tbuf->max = 1.2*(tbuf->top+apath->tlen) + TRACE_CHUNK;
                          tbuf->trace = Realloc(tbuf->trace,sizeof(short)*tbuf->max,
                                                "Reallocating trace vector");
                          if (tbuf->trace == NULL)
                            exit (1);
                        }
                      amatch[novla] = *apath;
                      amatch[novla].trace = (void *) (tbuf->top);
                      memmove(tbuf->trace+tbuf->top,apath->trace,sizeof(short)*apath->tlen);
                      novla += 1;
                      tbuf->top += apath->tlen;

#ifdef TEST_GATHER
                      printf("  [%5d,%5d] x [%5d,%5d] = %4d",
                             apath->abpos,apath->aepos,apath->bbpos,apath->bepos,apath->diffs);
#endif
#ifdef SHOW_OVERLAP
                      printf("\n\n                    %d(%d) vs %d(%d)\n\n",
                             ovla->aread,ovla->alen,ovla->bread,ovla->blen);
                      Print_ACartoon(stdout,align,ALIGN_INDENT);
#ifdef SHOW_ALIGNMENT
                      Compute_Trace_ALL(align,work);
                      printf("\n                      Diff = %d\n",align->path->diffs);
                      Print_Alignment(stdout,align,work,
                                      ALIGN_INDENT,ALIGN_WIDTH,ALIGN_BORDER,0,5);
#endif
#endif // SHOW_OVERLAP
                    }
#ifdef TEST_GATHER
                  else
                    printf("  No alignment %d",
                            ((apath->aepos-apath->abpos) + (apath->bepos-apath->bbpos))/2);
                  printf("\n");
#endif
                }
            }

          // Clear diagonal scores

          for (apos = amarkb; apos < amarke; apos++)
            { diag = aoff[apos].code;
              if (diag == 0) continue;
              diag >>= Binshift;
              score[diag] = lastp[diag] = 0;
            }

          if (amarke > alen) break;

          amarkb = amarke - PANEL_OVERLAP;
          amarke = amarkb + PANEL_SIZE;
          if (amarke > alen)
            amarke = alen+1;
        }

      // Clear diagonal last positions

      for (apos = Kmer; apos <= alen; apos++)
        { int d;

          diag = aoff[apos].code;
          if (diag == 0) continue;
          diag >>= Binshift;
          for (d = diag; d <= maxdiag; d++)
            if (lasta[d] == 0)
              break;
            else
              lasta[d] = 0;
          for (d = diag-1; d >= mindiag; d--)
            if (lasta[d] == 0)
              break;
            else
              lasta[d] = 0;
        }

      // Merge overlapping alignments and remove redundant ones

      { int i;

#ifdef TEST_CONTAIN
        if (novla > 1)
          printf("\n%5d vs %5d:\n",ar,ar);
#endif

        if (novla > 1)
          novla = Handle_Redundancies(amatch,novla,tbuf);

        for (i = 0; i < novla; i++)
          { ovla->path = amatch[i];
            ovla->path.trace = tbuf->trace + (uint64) (ovla->path.trace);
            if (small)
              Compress_TraceTo8(ovla,1);
#ifdef THREAD_OUTPUT
					if (Write_Overlap(ofile, ovla, tbytes))
					{
						fprintf(stderr, "%s: Cannot write to %s too small?\n",
								SORT_PATH, Prog_Name);
						exit(1);
					}
#endif
			  AddOverlapToBuffer(obuf, ovla, tbytes);
          }
        ahits += novla;
      }

      aoff += alen - (Kmer-1);
    }

  free(tbuf->trace);
  free(amatch);

  data->nfilt  = nfilt;
  data->ncheck = ahits;

#ifdef THREAD_OUTPUT
  rewind(ofile);
  fwrite(&ahits,sizeof(int64),1,ofile);
  fclose(ofile);
#endif

  return (NULL);
}


/*******************************************************************************************
 *
 *  THE ALGORITHM
 *
 ********************************************************************************************/

void Match_Self(char *aname, HITS_DB *ablock, Align_Spec *aspec)
{ THREAD     threads[NTHREADS];
  Merge_Arg  parmm[NTHREADS];
  Lex_Arg    parmx[NTHREADS];
  Report_Arg parmr[NTHREADS];
  int        pairsort[16];

  int64     nfilt, ncheck;

  KmerPos  *asort, *osort;
  int       alen;
  int64     atot;

  //  Setup

  atot = ablock->totlen;

  nfilt = ncheck = 0;

  { int64 powr;
    int   i, nbyte;

    for (i = 0; i < NTHREADS; i++)
      parmx[i].sptr = (int64 *) alloca(NTHREADS*BPOWR*sizeof(int64));

    for (i = 0; i < 16; i++)
      pairsort[i] = 0;

    powr = 1;
    for (nbyte = 0; powr < ablock->maxlen; nbyte += 1)
      powr <<= 8;
    for (i = 8; i < 8+nbyte; i++)
      pairsort[i] = 1;

    powr = 1;
    for (nbyte = 0; powr < ablock->nreads; nbyte += 1)
      powr <<= 8;
    for (i = 12; i < 12+nbyte; i++)
      pairsort[i] = 1;

    {
    		Overlap_IO_Buffer *buffer = OVL_IO_Buffer(aspec);
    		for (i = 0; i < NTHREADS; i++)
    		{
    			parmr[i].iobuf = &(buffer[i]);
    			if (parmr[i].iobuf == NULL)
    				exit(1);
    		}
    	}
  }

  //  Build K-mer sorted index of ablock

  if (VERBOSE)
    printf("\nIndexing %s\n",aname);

  osort = NULL;   //  Just to shut up dumb compilers
  asort = Sort_Kmers(ablock,&alen,&osort);

  if (VERBOSE)
    printf("\nComparing %s to itself\n",aname);

  if (alen == 0)
    goto zerowork;

  //  Determine equal K-mer links

  { int    i, p;
    uint64 c;

    parmm[0].abeg = 0;
    for (i = 1; i < NTHREADS; i++)
      { p = (int) ((((int64) alen) * i) >> NSHIFT);
        if (p > 0)
          { c = asort[p-1].code;
            while (asort[p].code == c)
              p += 1;
          }
        parmm[i].abeg = parmm[i-1].aend = p;
      }
    parmm[NTHREADS-1].aend = alen;

    for (i = 0; i < NTHREADS; i++)
      { parmm[i].kptr = parmx[i].tptr;
        for (p = 0; p < BPOWR; p++)
          parmm[i].kptr[p] = 0;
      }

    MG_alist = asort;

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,count_thread,parmm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#ifdef TEST_PAIRS
    printf("\nCROSS SORT %d:\n",alen);
    for (i = 0; i < HOW_MANY && i <= alen; i++)
      { KmerPos *c = asort+i;
        printf(" %5d / %5d / %10lld\n",c->read,c->rpos,c->code);
        fflush(stdout);
      }
#endif
  }

  //  Resort the list on read,rpos

  { int      i;
    KmerPos *rez;

    for (i = 0; i < NTHREADS; i++)
      { parmx[i].beg = parmm[i].abeg;
        parmx[i].end = parmm[i].aend;
      }

    rez = (KmerPos *) lex_sort(pairsort,(Double *) asort,(Double *) osort,parmx);
    if (rez != asort)
      { osort = asort;
        asort = rez;
      }

#ifdef TEST_CSORT
    printf("\nCROSS SORT %d:\n",alen);
    for (i = 0; i < HOW_MANY && i <= alen; i++)
      { KmerPos *c = asort+i;
        printf(" %5d / %5d / %10lld\n",c->read,c->rpos,c->code);
        fflush(stdout);
      }
#endif
  }

  //  Apply the diagonal filter and find local alignments about seed hits

  { int    i, w;
    int64  p;
    int    ar;
    int   *counters;

    MG_alist  = asort;
    MR_ablock = ablock;
    MR_spec   = aspec;
    MR_tspace = Trace_Spacing(aspec);

    asort[alen].read = ablock->nreads;
    parmr[0].beg = 0;
    for (i = 1; i < NTHREADS; i++)
      { p = (int) ((((int64) alen) * i) >> NSHIFT);
        if (p > 0)
          { ar = asort[p-1].read;
            while ((asort[p].read) == ar)
              p += 1;
          }
        parmr[i].beg = parmr[i-1].end = p;
      }
    parmr[NTHREADS-1].end = alen;

    w = ((ablock->maxlen >> Binshift) - ((-ablock->maxlen) >> Binshift)) + 1;
    counters = (int *) Malloc(NTHREADS*3*w*sizeof(int),"Allocating diagonal buckets");
    if (counters == NULL)
      exit (1);

    for (i = 0; i < 3*w*NTHREADS; i++)
      counters[i] = 0;
    for (i = 0; i < NTHREADS; i++)
      { if (i == 0)
          parmr[i].score = counters - ((-ablock->maxlen) >> Binshift);
        else
          parmr[i].score = parmr[i-1].lasta + w;
        parmr[i].lastp = parmr[i].score + w;
        parmr[i].lasta = parmr[i].lastp + w;
        parmr[i].work  = New_Work_Data();

#ifdef THREAD_OUTPUT
        parmr[i].ofile =
             Fopen(Catenate(SORT_PATH,"/",aname,Numbered_Suffix(".T",i+1,".las")),"w");
        if (parmr[i].ofile == NULL)
          exit (1);
#endif
      }

#ifdef NOTHREAD

    for (i = 0; i < NTHREADS; i++)
      report_thread(parmr+i);

#else

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,report_thread,parmr+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

#endif

    if (VERBOSE)
      for (i = 0; i < NTHREADS; i++)
        { nfilt  += parmr[i].nfilt;
          ncheck += parmr[i].ncheck;
        }

    for (i = 0; i < NTHREADS; i++)
      Free_Work_Data(parmr[i].work);
    free(counters);
  }

  //  Finish up

  goto epilogue;

zerowork:
  {
//	FILE *ofile;
//    int   i;
//
//    nfilt  = 0;
//    for (i = 0; i < NTHREADS; i++)
//      { ofile = Fopen(Catenate(SORT_PATH,"/",aname,Numbered_Suffix(".T",i+1,".las")),"w");
//        fwrite(&nfilt,sizeof(int64),1,ofile);
//        fwrite(&MR_tspace,sizeof(int),1,ofile);
//        fclose(ofile);
//      }
  }

epilogue:

  free(asort);
  free(osort);

  if (VERBOSE)
    { int width;

      if (nfilt <= 0)
        width = 1;
      else
        width = ((int) log10((double) nfilt)) + 1;
      width += (width-1)/3;

      printf("\n     ");
      Print_Number(nfilt,width,stdout);
      printf(" seed hits (%e of matrix)\n     ",(1.*nfilt/atot)/atot);
      Print_Number(ncheck,width,stdout);
      printf(" confirmed hits (%e of matrix)\n",(1.*ncheck/atot)/atot);
      fflush(stdout);
    }
}
