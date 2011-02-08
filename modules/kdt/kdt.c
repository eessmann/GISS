/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

#include "kdt.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* sort */

/* threshold below which in-memory qsort() is used rather than on-disk
   mergesort */
#define LENMIN 1000000

void kdt_heap_create (KdtHeap * h, int fd, long len)
{
  h->fd = fd;
  h->len = len;
  h->i = 0;
  h->p = malloc (len*sizeof (KdtPoint));
  h->end = read (fd, h->p, len*sizeof (KdtPoint));
  assert (h->end >= 0);
  h->end /= sizeof (KdtPoint);
}

int kdt_heap_get (KdtHeap * h, KdtPoint * p)
{
  if (h->i < h->end) {
    *p = h->p[h->i++];
    return 1;
  }
  if (h->end < h->len)
    return 0;
  h->end = read (h->fd, h->p, h->len*sizeof (KdtPoint));
  assert (h->end >= 0);
  h->end /= sizeof (KdtPoint);
  h->i = 0;
  return kdt_heap_get (h, p);
}

void kdt_heap_put (KdtHeap * h, KdtPoint * p)
{
  if (h->i == h->len) {
    assert (write (h->fd, h->p, h->len*sizeof (KdtPoint)) == h->len*sizeof (KdtPoint));
    h->i = 0;
  }
  h->p[h->i++] = *p;
}

void kdt_heap_flush (KdtHeap * h)
{
  if (h->i > 0)
    assert (write (h->fd, h->p, h->i*sizeof (KdtPoint)) == h->i*sizeof (KdtPoint));
}

void kdt_heap_free (KdtHeap * h)
{
  free (h->p);
}

static int put (KdtHeap * h, KdtPoint * p, KdtHeap * merged)
{
  kdt_heap_put (merged, p);
  return kdt_heap_get (h, p);
}

static int fdtemp (void)
{
  char name[] = "XXXXXX";
  int fd = mkstemp (name);
  assert (fd >= 0);
  assert (unlink (name) == 0);
  return fd;
}

static int merge (int fd1, int fd2, int (*compar) (const void *, const void *), long len)
{
  KdtHeap h1, h2, hm;
  assert (lseek (fd1, 0, SEEK_SET) == 0);
  kdt_heap_create (&h1, fd1, LENMIN/3);
  assert (lseek (fd2, 0, SEEK_SET) == 0);
  kdt_heap_create (&h2, fd2, LENMIN/3);
  int merged = fdtemp ();
  kdt_heap_create (&hm, merged, LENMIN/3);
  KdtPoint p1, p2;
  int r1 = kdt_heap_get (&h1, &p1);
  int r2 = kdt_heap_get (&h2, &p2);
  while (r1 && r2) {
    if ((* compar) (&p2, &p1))
      r1 = put (&h1, &p1, &hm);
    else
      r2 = put (&h2, &p2, &hm);
  }
  while (r1)
    r1 = put (&h1, &p1, &hm);
  while (r2)
    r2 = put (&h2, &p2, &hm);
  kdt_heap_free (&h1);
  kdt_heap_free (&h2);
  close (fd1);
  close (fd2);
  kdt_heap_flush (&hm);
  kdt_heap_free (&hm);
  return merged;
}

static int half (int fd, long len)
{
  long len1 = len/2;
  assert (lseek (fd, len1*sizeof (KdtPoint), SEEK_SET) > 0);
  int fd1 = fd, fd2 = fdtemp ();
  char a[4096];
  ssize_t size;
  while ((size = read (fd1, a, 4096)) > 0)
    assert (write (fd2, a, size) == size);
  assert (ftruncate (fd1, len1*sizeof (KdtPoint)) == 0);
  return fd2;
}

#if TIMING
static double elapsed (const struct timeval * start, const struct timeval * end)
{
  return (double) (end->tv_usec - start->tv_usec) + 1e6*(double) (end->tv_sec - start->tv_sec);
}
#endif

static int sort (int fd, long len, 
		 int  (*compar)   (const void *, const void *),
		 void (*progress) (void *), void * data)
{
  struct timeval start;
  gettimeofday (&start, NULL);
  assert (lseek (fd, 0, SEEK_SET) == 0);
  if (len <= LENMIN) {
    KdtPoint * a = malloc (len*sizeof (KdtPoint));
    assert (read (fd, a, len*sizeof (KdtPoint)) == len*sizeof (KdtPoint));
    qsort (a, len, sizeof (KdtPoint), compar);
    assert (lseek (fd, 0, SEEK_SET) == 0);
    assert (write (fd, a, len*sizeof (KdtPoint)) == len*sizeof (KdtPoint));
    free (a);
#if TIMING
    struct timeval end;
    gettimeofday (&end, NULL);
    fprintf (stderr, "sort %ld %g\n", len, elapsed (&start, &end));
#endif
    if (progress)
      (* progress) (data);
    return fd;
  }
  else {
    long len1 = len/2, len2 = len - len1;
    int fd2 = half (fd, len);
    int fd3 = merge (sort (fd, len1, compar, progress, data), 
		     sort (fd2, len2, compar, progress, data), 
		     compar, len);
#if TIMING
    struct timeval end;
    gettimeofday (&end, NULL);
    fprintf (stderr, "sort %ld %g\n", len, elapsed (&start, &end));
#endif
    return fd3;
  }
}

/* Kdt */

typedef struct {
  double split;
} Node;

#define SYSTEM_32_BITS (!defined (__LP64__) && !defined (__64BIT__) && \
                        !defined (_LP64) && !(__WORDSIZE == 64))

typedef struct {
  KdtRect bound;
  long len;
  PADDING_32_BITS;
  long np;
#if SYSTEM_32_BITS
  int padding1;
#endif
} Header;

struct _Kdt {
  Header h;
  FILE * nodes, * sums, * leaves;
  KdtPoint * buffer;
  /* progress stuff */
  void (* progress) (float complete, void * data);
  void * data;
  int i, m;
};

static int check_32_bits (const Kdt * kdt)
{
#if SYSTEM_32_BITS
  long maxlen = (1 << 31)/sizeof (KdtPoint);
  if (kdt->h.len > maxlen) {
    fprintf (stderr, "kdt: 32-bits systems are limited to %ld data points\n", maxlen);
    return 1;
  }
#endif
  return 0;
}

#define KDTSIZE(len) (((len) - 1)*sizeof (Node) + (len)*sizeof (KdtPoint))

int kdt_size (const Kdt * kdt, long len)
{
  long len1 = 1;
  while (len > kdt->h.np) {
    len /= 2; len1 *= 2;
  }
  return (len1 - 1)*sizeof (Node) + len1*sizeof (KdtPoint)*kdt->h.np;
}

void kdt_sizes (const Kdt * kdt, long len, long * nodes, long * sums, long * leaves)
{
  long len1 = 1;
  while (len > kdt->h.np) {
    len /= 2; len1 *= 2;
  }
  *nodes = (len1 - 1)*sizeof (Node);
  *sums  = (len1 - 1)*sizeof (KdtSum);
  *leaves = len1*sizeof (KdtPoint)*kdt->h.np;
}

static void relative (const KdtRect rect, double * o, double * h)
{
  o[0] = (((double)rect[0].l) + ((double)rect[0].h))/2.;
  o[1] = (((double)rect[1].l) + ((double)rect[1].h))/2.;
  *h = ((double)rect[0].h) - ((double)rect[0].l);
  if (((double)rect[1].h) - ((double)rect[1].l) > *h)
    *h = ((double)rect[1].h) - ((double)rect[1].l);
}

static void sum_add_point (const KdtRect parent, KdtSum * sum, const KdtPoint * a)
{
  double p[3], o[2], h;

  relative (parent, o, &h);

  p[0] = (a->x - o[0])/h; p[1] = (a->y - o[1])/h; p[2] = a->z;
#if AVG_TERRAIN
  sum->H0 += p[2];
  sum->n++;
  if (p[2] < sum->Hmin) sum->Hmin = p[2];
  if (p[2] > sum->Hmax) sum->Hmax = p[2];
#else
  sum->m01 += p[0];
  sum->m02 += p[1];
  sum->m03 += p[0]*p[1];
  sum->m11 += p[0]*p[0];
  sum->m13 += p[0]*p[0]*p[1];
  sum->m22 += p[1]*p[1];
  sum->m23 += p[0]*p[1]*p[1];
  sum->m33 += p[0]*p[0]*p[1]*p[1];
  sum->m44 += p[0]*p[0]*p[0];
  sum->m55 += p[1]*p[1]*p[1];
  sum->m66 += p[0]*p[0]*p[0]*p[0];
  sum->m77 += p[1]*p[1]*p[1]*p[1];
  sum->m67 += p[0]*p[0]*p[0]*p[1];
  sum->m76 += p[1]*p[1]*p[1]*p[0];
  sum->H0 += p[2];
  sum->H1 += p[0]*p[2];
  sum->H2 += p[1]*p[2];
  sum->H3 += p[0]*p[1]*p[2];
  sum->H4 += p[2]*p[2];
  sum->H5 += p[0]*p[0]*p[2];
  sum->H6 += p[1]*p[1]*p[2];
  sum->n++;
  if (p[2] < sum->Hmin) sum->Hmin = p[2];
  if (p[2] > sum->Hmax) sum->Hmax = p[2];
#endif
}

static void sum_add_sum (const KdtRect parent, KdtSum * sum, const KdtRect rect, const KdtSum * a)
{
  double op[2], oa[2], hp, ha;

  relative (parent, op, &hp);
  relative (rect, oa, &ha);
  
  double oap0 = oa[0] - op[0], oap1 = oa[1] - op[1];
  double an = a->n;
  double ha2 = ha*ha, hp2 = hp*hp;
  sum->m01 += (an*oap0 + a->m01*ha)/hp;
  sum->m02 += (an*oap1 + a->m02*ha)/hp;
  sum->m03 += (oap0*(an*oap1 + a->m02*ha) + ha*(a->m01*oap1 + a->m03*ha))/hp2;
  double m11 = (oap0*(an*oap0 + 2.*a->m01*ha) + a->m11*ha2)/hp2;
  sum->m11 += m11;
  double m13 = ha*(oap0*(a->m02*oap0 + 2.*a->m03*ha) + a->m13*ha2)/hp2;
  sum->m13 += (oap1*m11 + m13)/hp;
  double m22 = (oap1*(an*oap1 + 2.*a->m02*ha) + a->m22*ha2)/hp2;
  sum->m22 += m22;
  sum->m23 += (oap0*m22 + ha*(oap1*(oap1*a->m01 + 2.*a->m03*ha) + a->m23*ha2)/hp2)/hp;
  sum->m33 += (oap1*(oap1*m11 + 2.*m13) + 
		ha2*(oap0*(oap0*a->m22 + 2.*a->m23*ha) + ha2*a->m33)/hp2)/hp2;
  double ha3 = ha2*ha, hp3 = hp2*hp;
  sum->m44 += (oap0*(oap0*(oap0*an + 3.*ha*a->m01) + 3.*ha2*a->m11) + ha3*a->m44)/hp3;
  sum->m55 += (oap1*(oap1*(oap1*an + 3.*ha*a->m02) + 3.*ha2*a->m22) + ha3*a->m55)/hp3;
  double ha4 = ha3*ha, hp4 = hp3*hp;
  sum->m66 += (oap0*(oap0*(oap0*(oap0*an + 4.*ha*a->m01) + 6.*ha2*a->m11) 
		      + 4.*ha3*a->m44) + ha4*a->m66)/hp4;
  sum->m77 += (oap1*(oap1*(oap1*(oap1*an + 4.*ha*a->m02) + 6.*ha2*a->m22)
		      + 4.*ha3*a->m55) + ha4*a->m77)/hp4;
  sum->m67 += (oap1*(oap0*(oap0*(oap0*an + 3.*ha*a->m01) + 3.*ha2*a->m11) + ha3*a->m44)
		+ oap0*(oap0*(ha*a->m02*oap0 + 3.*ha2*a->m03) + 3.*ha3*a->m13) 
		+ ha4*a->m67)/hp4;
  sum->m76 += (oap0*(oap1*(oap1*(oap1*an + 3.*ha*a->m02) + 3.*ha2*a->m22) + ha3*a->m55)
		+ oap1*(oap1*(ha*a->m01*oap1 + 3.*ha2*a->m03) + 3.*ha3*a->m23)
		+ ha4*a->m76)/hp4;
  sum->H0 += a->H0;
  sum->H1 += (a->H0*oap0 + a->H1*ha)/hp;
  sum->H2 += (a->H0*oap1 + a->H2*ha)/hp;
  sum->H3 += (ha*(ha*a->H3 + oap0*a->H2 + oap1*a->H1) + oap0*oap1*a->H0)/hp2;
  sum->H4 += a->H4;
  sum->H5 += (oap0*(2.*ha*a->H1 + oap0*a->H0) + ha2*a->H5)/hp2;
  sum->H6 += (oap1*(2.*ha*a->H2 + oap1*a->H0) + ha2*a->H6)/hp2;
  sum->n += a->n;
  if (a->Hmin < sum->Hmin) sum->Hmin = a->Hmin;
  if (a->Hmax > sum->Hmax) sum->Hmax = a->Hmax;
}

static int sort_x (const void * p1, const void * p2)
{
  return ((KdtPoint *) p1)->x > ((KdtPoint *) p2)->x;
}

static int sort_y (const void * p1, const void * p2)
{
  return ((KdtPoint *) p1)->y > ((KdtPoint *) p2)->y;
}

static void update_sum (const KdtRect rect, KdtSum * n, KdtPoint * a, long len)
{
  kdt_sum_init (n);
  long i;
  for (i = 0; i < len; i++, a++)
    sum_add_point (rect, n, a);
}

#define PADDING 1e100

static void split (KdtPoint * a, KdtRect bound, long len, int index, Kdt * kdt)
{
  if (len > kdt->h.np) {
    //    fprintf (stderr, " splitting: %ld       \r", len);
    qsort (a, len, sizeof (KdtPoint), index == 0 ? sort_x : sort_y);
    long len1 = len/2;
    long len2 = len - len1;
    Node n;
    n.split = (&a[len1].x)[index];
    KdtSum s;
    update_sum (bound, &s, a, len);
    fwrite (&s, sizeof (KdtSum), 1, kdt->sums);
    fwrite (&n, sizeof (Node), 1, kdt->nodes);
    KdtRect bound1;
    bound1[0].l = bound[0].l; bound1[0].h = bound[0].h;
    bound1[1].l = bound[1].l; bound1[1].h = bound[1].h;
    bound1[index].h = n.split;
    split (a, bound1, len1, (index + 1) % 2, kdt);

    bound1[0].l = bound[0].l; bound1[0].h = bound[0].h;
    bound1[1].l = bound[1].l; bound1[1].h = bound[1].h;
    bound1[index].l = n.split;
    split (&(a[len1]), bound1, len2, (index + 1) % 2, kdt);
  }
  else if (len > 0) {
    fwrite (a, sizeof (KdtPoint), len, kdt->leaves);
    /* padding */
    KdtPoint dummy = { PADDING, PADDING, PADDING };
    int np = kdt->h.np;
    while (np-- > len)
      fwrite (&dummy, sizeof (KdtPoint), 1, kdt->leaves);
  }
}

static void file_update_sum (const KdtRect rect, KdtSum * n, int fd, long len)
{
  KdtHeap h;
  kdt_heap_create (&h, fd, LENMIN);
  kdt_sum_init (n);
  long i;
  for (i = 0; i < len; i++) {
    KdtPoint p;
    assert (kdt_heap_get (&h, &p));
    sum_add_point (rect, n, &p);
  }
  kdt_heap_free (&h);
}

static void progress (void * data)
{
  Kdt * kdt = data;
  if (kdt->progress)
    (* kdt->progress) (++kdt->i/(float) kdt->m, kdt->data);
}

static void splitfile (int fd, KdtRect bound, long len, int index, Kdt * kdt)
{
#if TIMING
  struct timeval start;
  gettimeofday (&start, NULL);
#endif
  assert (lseek (fd, 0, SEEK_SET) == 0);
  if (len > LENMIN) {
    //    fprintf (stderr, " splitting: %ld      \r", len);
    KdtSum s;
    file_update_sum (bound, &s, fd, len);
    fwrite (&s, sizeof (KdtSum), 1, kdt->sums);

    long len1 = len/2;
    long len2 = len - len1;
    fd = sort (fd, len,  index == 0 ? sort_x : sort_y, progress, kdt);
#if TIMING
    struct timeval s1;
    gettimeofday (&s1, NULL);
#endif
    int fd2 = half (fd, len);
#if TIMING
    struct timeval end;
    gettimeofday (&end, NULL);
    fprintf (stderr, "half %ld %g\n", len, elapsed (&s1, &end));
#endif
    KdtPoint p;
    assert (lseek (fd2, 0, SEEK_SET) == 0);
    assert (read (fd2, &p, sizeof (KdtPoint)) == sizeof (KdtPoint));

    Node n;
    n.split = (&p.x)[index];
    fwrite (&n, sizeof (Node), 1, kdt->nodes);

    KdtRect bound1;
    bound1[0].l = bound[0].l; bound1[0].h = bound[0].h;
    bound1[1].l = bound[1].l; bound1[1].h = bound[1].h;
    bound1[index].h = n.split;
    splitfile (fd, bound1, len1, (index + 1) % 2, kdt);
    close (fd);

    bound1[0].l = bound[0].l; bound1[0].h = bound[0].h;
    bound1[1].l = bound[1].l; bound1[1].h = bound[1].h;
    bound1[index].l = n.split;
    splitfile (fd2, bound1, len2, (index + 1) % 2, kdt);
    close (fd2);
  }
  else if (len > 0) {
    KdtPoint * a = malloc (len*sizeof (KdtPoint));
    assert (read (fd, a, len*sizeof (KdtPoint)) == len*sizeof (KdtPoint));
    split (a, bound, len, index, kdt);
    free (a);
  }
#if TIMING
  struct timeval end;
  gettimeofday (&end, NULL);
  fprintf (stderr, "splitfile %ld %g\n", len, elapsed (&start, &end));
#endif
}

Kdt * kdt_new (void)
{
  Kdt * kdt = calloc (1, sizeof (Kdt));
  return kdt;
}

static FILE * open_ext (const char * name, const char * ext, const char * mode)
{
  int len = strlen (name), len1 = strlen (ext);
  char * fname = malloc (sizeof(char)*(len + len1 + 1));
  strcpy (fname, name);
  strcpy (&fname[len], ext);
  FILE * fp = fopen (fname, mode);
  free (fname);
  return fp;
}

static int kdt_init (Kdt * kdt, const char * name, int npmax, long len)
{
  kdt->nodes  = open_ext (name, ".kdt", "w");
  if (!kdt->nodes)
    return -1;
  
  kdt->sums   = open_ext (name, ".sum", "w");
  if (!kdt->sums)
    return -1;

  kdt->leaves = open_ext (name, ".pts", "w");
  if (!kdt->leaves)
    return -1;

  kdt->h.len = len;
  kdt->h.np = len;
  while (kdt->h.np > npmax)
    kdt->h.np /= 2;
  kdt->h.np++;
  kdt->h.bound[0].l = kdt->h.bound[1].l =  1e30;
  kdt->h.bound[0].h = kdt->h.bound[1].h = -1e30;

  if (check_32_bits (kdt))
    return -1;
  
  return 0;
}

int kdt_create (Kdt * kdt, const char * name, int blksize,
		KdtPoint * a, long len)
{
  int npmax = blksize/sizeof (KdtPoint);
  if (kdt_init (kdt, name, npmax, len))
    return -1;
  
  long i;
  for (i = 0; i < len; i++) {
    if (a[i].x > kdt->h.bound[0].h) kdt->h.bound[0].h = a[i].x;
    if (a[i].x < kdt->h.bound[0].l) kdt->h.bound[0].l = a[i].x;
    if (a[i].y > kdt->h.bound[1].h) kdt->h.bound[1].h = a[i].y;
    if (a[i].y < kdt->h.bound[1].l) kdt->h.bound[1].l = a[i].y;
  }

  fwrite (&kdt->h, sizeof (Header), 1, kdt->nodes);
  split (a, kdt->h.bound, kdt->h.len, 0, kdt);

  return 0;
}

int kdt_create_from_file (Kdt * kdt, const char * name, int blksize,
			  int fd,
			  void (* progress) (float complete, void * data),
			  void * data)
{
  KdtHeap h;
  kdt_heap_create (&h, fd, LENMIN);
  long len = 0;
  KdtPoint p;
  KdtRect bound = {{ 1e30, -1e30}, {1e30, -1e30}};
  while (kdt_heap_get (&h, &p)) {
    if (p.x > bound[0].h) bound[0].h = p.x;
    if (p.x < bound[0].l) bound[0].l = p.x;
    if (p.y > bound[1].h) bound[1].h = p.y;
    if (p.y < bound[1].l) bound[1].l = p.y;
    len++;
  }
  kdt_heap_free (&h);

  int npmax = blksize/sizeof (KdtPoint);
  if (kdt_init (kdt, name, npmax, len))
    return -1;
  kdt->h.bound[0].l = bound[0].l; kdt->h.bound[0].h = bound[0].h;
  kdt->h.bound[1].l = bound[1].l; kdt->h.bound[1].h = bound[1].h;
  
  fwrite (&kdt->h, sizeof (Header), 1, kdt->nodes);
  kdt->m = kdt->i = 0;
  int m2 = 1;
  while (len > LENMIN) {
    kdt->m++;
    len /= 2;
    m2 *= 2;
  }
  kdt->m = kdt->m*m2;
  kdt->progress = progress;
  kdt->data = data;
  splitfile (fd, kdt->h.bound, kdt->h.len, 0, kdt);

  return 0;
}

int kdt_open (Kdt * kdt, const char * name)
{
  kdt->nodes  = open_ext (name, ".kdt", "r");
  if (!kdt->nodes)
    return -1;
  
  kdt->sums   = open_ext (name, ".sum", "r");
  if (!kdt->sums)
    return -1;

  kdt->leaves = open_ext (name, ".pts", "r");
  if (!kdt->leaves)
    return -1;

  if (fread (&kdt->h, sizeof (Header), 1, kdt->nodes) != 1)
    return -1;

  kdt->buffer = malloc (sizeof (KdtPoint)*kdt->h.np);

  if (check_32_bits (kdt))
    return -1;

  return 0;
}

void kdt_destroy (Kdt * kdt)
{
  if (kdt->nodes)
    fclose (kdt->nodes);
  if (kdt->sums)
    fclose (kdt->sums);
  if (kdt->leaves)
    fclose (kdt->leaves);
  if (kdt->buffer)
    free (kdt->buffer);
  free (kdt);
}

static long query (const Kdt * kdt, const KdtRect rect, int index, long len)
{
  if (len > kdt->h.np) {
    Node node;
    long len1 = len/2;
    if (fread (&node, sizeof (Node), 1, kdt->nodes) != 1)
      return -1;
    long pos = ftell (kdt->nodes), lpos = ftell (kdt->leaves);
    if (pos < 0 || lpos < 0)
      return -1;
    long n = 0;
    if (rect[index].l < node.split) {
      long n1 = query (kdt, rect, (index + 1) % 2, len1);
      if (n1 < 0)
	return -1;
      n += n1;
    }
    if (rect[index].h >= node.split) {
      long snodes, ssums, sleaves;
      kdt_sizes (kdt, len1, &snodes, &ssums, &sleaves);
      if (fseek (kdt->nodes, pos + snodes, SEEK_SET))
	return -1;
      if (fseek (kdt->leaves, lpos + sleaves, SEEK_SET))
	return -1;
      long n1 = query (kdt, rect, (index + 1) % 2, len - len1);
      if (n1 < 0)
	return -1;
      n += n1;
    }
    return n;
  }
  else if (len > 0) {
    if (fread (kdt->buffer, sizeof (KdtPoint), len, kdt->leaves) != len)
      return -1;
    int i, n = 0;
    for (i = 0; i < len; i++)
      if (kdt->buffer[i].x >= rect[0].l && kdt->buffer[i].x <= rect[0].h && 
	  kdt->buffer[i].y >= rect[1].l && kdt->buffer[i].y <= rect[1].h) {
       	printf ("%.8f %.8f %g\n", kdt->buffer[i].x, kdt->buffer[i].y, kdt->buffer[i].z);
	n++;
      }
    return n;
  }
  return 0;
}

long kdt_query (const Kdt * kdt, const KdtRect rect)
{
  rewind (kdt->nodes);
  rewind (kdt->leaves);
  Header h;
  if (fread (&h, sizeof (Header), 1, kdt->nodes) != 1)
    return -1;
  if (rect[0].l > h.bound[0].h || rect[1].l > h.bound[1].h || 
      rect[0].h < h.bound[0].l || rect[1].h < h.bound[1].l)
    return 0;
  return query (kdt, rect, 0, h.len);
}

typedef struct {
  long np, sp, lp;
} FilePointers;

static long query_sum (const Kdt * kdt,
		       KdtCheck includes, KdtCheck intersects, void * data, 
		       KdtRect bound, int index, long len,
		       FilePointers * f,
		       const KdtRect rect, KdtSum * sum)
{
  if (len > kdt->h.np) {
    Node node;
    long len1 = len/2;
    if (fseek (kdt->nodes, f->np, SEEK_SET))
      return -1;
    if (fread (&node, sizeof (Node), 1, kdt->nodes) != 1)
      return -1;
    f->np += sizeof (Node);
#if DEBUG
    fprintf (stderr, "read 1 node %ld\n", sizeof (Node));
#endif

    if ((* includes) (bound, data)) {
      KdtSum s;
      if (fseek (kdt->sums, f->sp, SEEK_SET))
	return -1;
      if (fread (&s, sizeof (KdtSum), 1, kdt->sums) != 1)
	return -1;
#if DEBUG
      fprintf (stderr, "read 1 sum %ld index %d\n", sizeof (KdtSum), index);
#endif
      f->sp += sizeof (KdtSum);
      sum_add_sum (rect, sum, bound, &s);
      return len;
    }
    f->sp += sizeof (KdtSum);

    long pos = f->np, lpos = f->lp, spos = f->sp;
    long n = 0;

    KdtRect bound1;
    bound1[0].l = bound[0].l; bound1[0].h = bound[0].h;
    bound1[1].l = bound[1].l; bound1[1].h = bound[1].h;
    bound1[index].h = node.split;
    if ((* intersects) (bound1, data)) {
      long n1 = query_sum (kdt, includes, intersects, data, 
			   bound1, (index + 1) % 2, len1, f, rect, sum);
      if (n1 < 0)
	return -1;
      n += n1;
    }

    bound1[0].l = bound[0].l; bound1[0].h = bound[0].h;
    bound1[1].l = bound[1].l; bound1[1].h = bound[1].h;
    bound1[index].l = node.split;
    if ((* intersects) (bound1, data)) {
      long snodes, ssums, sleaves;
      kdt_sizes (kdt, len1, &snodes, &ssums, &sleaves);
      f->np = pos + snodes;
      f->sp = spos + ssums;
      f->lp = lpos + sleaves;
      long n1 = query_sum (kdt, includes, intersects, data, 
			   bound1, (index + 1) % 2, len - len1, f, rect, sum);
      if (n1 < 0)
	return -1;
      n += n1;
    }
    return n;
  }
  else if (len > 0) {
    if (fseek (kdt->leaves, f->lp, SEEK_SET))
      return -1;
    if (fread (kdt->buffer, sizeof (KdtPoint), len, kdt->leaves) != len)
      return -1;
#if DEBUG
    fprintf (stderr, "read %ld leaves %ld\n", len, sizeof (KdtPoint));
#endif
    KdtRect boundp;
    KdtPoint * a = kdt->buffer;
    int i, n = 0;
    for (i = 0; i < len; i++, a++) {
      boundp[0].l = boundp[0].h = a->x;
      boundp[1].l = boundp[1].h = a->y;
      if ((* includes) (boundp, data)) {
	sum_add_point (rect, sum, a);
	n++;
      }
    }
    return n;
  }
  return 0;
}

long kdt_query_sum (const Kdt * kdt, 
		    KdtCheck includes, KdtCheck intersects, void * data, 
		    const KdtRect rect, KdtSum * sum)
{
  rewind (kdt->nodes);
  rewind (kdt->leaves);
  Header h;
  if (fread (&h, sizeof (Header), 1, kdt->nodes) != 1)
    return -1;
  FilePointers f;
  f.np = sizeof (Header);
  f.sp = f.lp = 0;
  return query_sum (kdt, includes, intersects, data, h.bound, 0, h.len, &f, rect, sum);
}

void kdt_sum_init (KdtSum * s)
{
  memset (s, 0, sizeof (KdtSum));
  s->Hmax = - 1e30;
  s->Hmin =   1e30;
}
