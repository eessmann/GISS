#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <glib.h>
#include <gts.h>

#ifndef PI
# define PI 3.14159265359
#endif

#define K 3.0
#define L 3.0

static void ** malloc2D (guint nx, guint ny, gulong size)
{
  void ** m = g_malloc (nx*sizeof (void *));
  guint i;

  for (i = 0; i < nx; i++)
    m[i] = g_malloc0 (ny*size);

  return m;
}

static void free2D (void ** m, guint nx)
{
  guint i;

  g_return_if_fail (m != NULL);

  for (i = 0; i < nx; i++)
    g_free (m[i]);
  g_free (m);
}

static void fill0 (gdouble **u, gint n)
{
  gint i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      u[i][j] = 0.;
}

static void symbc (gdouble **u, gint n)
{
  gint i;

  for (i=1;i<=n;i++) {
    u[i][1] = u[i][2];
    u[i][n] = u[i][n-1];
    u[1][i] = u[2][i];
    u[n][i] = u[n-1][i];
  }
}

static void relax (gdouble **u, gdouble **rhs, gint n)
{
  gint i,ipass,isw,j,jsw=1; 
  gdouble h,h2;

  h=1.0/(n-2); 
  h2=h*h; 

  for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
    isw=jsw; 
    for (j=2;j<n;j++,isw=3-isw)
      for (i=isw+1;i<n;i+=2)
	u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]		      
		      +u[i][j-1]-h2*rhs[i][j]);
  }
  symbc (u, n);
}

static void relax2 (gdouble **u, gdouble **rhs, gint n)
{
  gint i,j; 
  gdouble h,h2;

  h=1.0/(n-2); 
  h2=h*h; 

  for (i = 2; i < n; i++)
    for (j = 2; j < n; j++)
      u[i][j] = 0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]		      
		      +u[i][j-1]-h2*rhs[i][j]);
  symbc (u, n);
}

static void residual (double **res, double **u, double **rhs, int n)
{
  gint i,j; 
  gdouble h,h2i;

  h=1.0/(n-2);
  h2i=1.0/(h*h); 
  for (j=2;j<n;j++)
    for (i=2;i<n;i++)
      res[i][j] = rhs[i][j]-h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
				 4.0*u[i][j]); 
  for (i=1;i<=n;i++)
    res[i][1]=res[i][n]=res[1][i]=res[n][i]=0.0;
}

static void init_div (gdouble ** rhs, gint n)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++) {
      gdouble x = ((gdouble) i - (gdouble) (n + 1)/2.)/(gdouble) (n - 2);
      gdouble y = ((gdouble) j - (gdouble) (n + 1)/2.)/(gdouble) (n - 2);

      rhs[i][j] = - PI*PI*(K*K + L*L)*sin (PI*K*x)*sin (PI*L*y);
    }
}

static void write (gdouble ** val, gint n, FILE * fp)
{
  gint i, j;

  for (i = 2; i < n; i++) {
    for (j = 2; j < n; j++) {
      gdouble x = ((gdouble) i - (gdouble) (n + 1)/2.)/(gdouble) (n - 2);
      gdouble y = ((gdouble) j - (gdouble) (n + 1)/2.)/(gdouble) (n - 2);

      fprintf (fp, "%g %g %g\n", x, y, val[i][j]);      
    }      
    fprintf (fp, "\n");
  }
}

static GtsRange stats (gdouble ** val, gint n)
{
  GtsRange s;
  gint i, j;

  gts_range_init (&s);
  for (i = 2; i < n; i++)
    for (j = 2; j < n; j++)
      gts_range_add_value (&s, val[i][j]);
  gts_range_update (&s);
  
  return s;
}

static void restrictor (gdouble ** w, gdouble ** w1, gint n)
{
  gint i, j;

  for (i = 2; i < n; i++)
    for (j = 2; j < n; j++)
      w1[i][j] = (w[2*i - 2][2*j - 2] + w[2*i - 1][2*j - 2] +
		  w[2*i - 2][2*j - 1] + w[2*i - 1][2*j - 1])/4.;
}

static void interpolate (gdouble ** w, gdouble ** w1, gint n)
{
  gint i, j;

  for (i = 2; i < n; i++)
    for (j = 2; j < n; j++) {
      w[2*i - 2][2*j - 2] = w1[i][j];
      w[2*i - 1][2*j - 2] = w1[i][j];
      w[2*i - 2][2*j - 1] = w1[i][j];
      w[2*i - 1][2*j - 1] = w1[i][j];
    }
}

static void addp (gdouble ** p, gdouble ** p1, gint n)
{
  gint i, j;

  for (i = 2; i < n; i++)
    for (j = 2; j < n; j++) {
      p[2*i - 2][2*j - 2] += p1[i][j];
      p[2*i - 1][2*j - 2] += p1[i][j];
      p[2*i - 2][2*j - 1] += p1[i][j];
      p[2*i - 1][2*j - 1] += p1[i][j];
    }
}

static void addp1 (gdouble ** p, gdouble ** p1, gint n)
{
  gint i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      p[i][j] += p1[i][j];
}

static void multigrid (gdouble ** u, gdouble ** rhs, gdouble ** res, guint n,
		       guint vcycle, guint npre, guint npost)
{
  gdouble ** rhs1 = NULL, ** u1 = NULL, ** res1 = NULL;
  guint nc, i, j;

  g_return_if_fail (n % 2 == 0);

  nc = n/2 + 1;
  if (nc > 3) {
    rhs1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));
    u1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));
    res1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));
    restrictor (rhs, rhs1, nc);
    multigrid (u1, rhs1, res1, nc, vcycle, npre, npost);
    interpolate (u, u1, nc);
    symbc (u, n);
  }
  else
    fill0 (u, n);

  for (i = 0; i < npre; i++)
    relax (u, rhs, n);

  if (nc > 3) {
    for (j = 0; j < vcycle; j++) {
      residual (res, u, rhs, n);
      restrictor (res, rhs1, nc);
      multigrid (u1, rhs1, res1, nc, vcycle, npre, npost);
      addp (u, u1, nc);
      symbc (u, n);
      for (i = 0; i < npost; i++)
	relax (u, rhs, n);
    }
    free2D ((void **) rhs1, nc);
    free2D ((void **) u1, nc);
    free2D ((void **) res1, nc);
  }
}

static void multigrid1 (gdouble ** u, gdouble ** rhs, gdouble ** res, guint n,
			guint npre, guint npost)
{
  gdouble ** rhs1, ** u1, ** res1;
  guint nc, i;

  g_return_if_fail (n % 2 == 0);

  for (i = 0; i < npre; i++)
    relax (u, rhs, n);

  nc = n/2 + 1;
  if (nc > 3) {
    rhs1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));
    u1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));
    res1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));

    residual (res, u, rhs, n);
    restrictor (res, rhs1, nc);
    multigrid1 (u1, rhs1, res1, nc, npre, npost);
    addp (u, u1, nc);
    symbc (u, n);

    for (i = 0; i < npost; i++)
      relax (u, rhs, n);

    free2D ((void **) rhs1, nc);
    free2D ((void **) u1, nc);
    free2D ((void **) res1, nc);
  }
}

static void multigrid2 (gdouble ** u, gdouble ** rhs, 
			guint n, guint npre, guint npost)
{
  guint nc, i;

  g_return_if_fail (n % 2 == 0);

  nc = n/2 + 1;
  if (nc > 3) {
    gdouble ** u1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));
    gdouble ** rhs1 = (gdouble **) malloc2D (nc + 1, nc + 1, sizeof (gdouble));

    restrictor (rhs, rhs1, nc);
    multigrid2 (u1, rhs1, nc, npre, npost);
    interpolate (u, u1, nc);
    symbc (u, n);

    for (i = 0; i < npost; i++)
      relax (u, rhs, n);

    free2D ((void **) rhs1, nc);
    free2D ((void **) u1, nc);
  }
  else
    for (i = 0; i < 10*npost; i++)
      relax (u, rhs, n);
}

int main (int argc, char * argv[])
{
  gdouble ** p, ** rhs, ** res, ** dp;
  guint n = 1, vcycle, npre, npost, level, i;
  GtsRange s;
  GTimer * timer;

  if (argc != 5) {
    fprintf (stderr, "usage: simple LEVELS VCYCLE NPRE NPOST\n");
    return 1;
  }

  i = level = atoi (argv[1]);
  while (i) {
    n *= 2;
    i--;
  }
  n += 2;
  vcycle = atoi (argv[2]);
  npre = atoi (argv[3]);
  npost = atoi (argv[4]);
  
  p = (gdouble **) malloc2D (n + 1, n + 1, sizeof (gdouble));
  rhs = (gdouble **) malloc2D (n + 1, n + 1, sizeof (gdouble));
  res = (gdouble **) malloc2D (n + 1, n + 1, sizeof (gdouble));
  dp = (gdouble **) malloc2D (n + 1, n + 1, sizeof (gdouble));

  init_div (rhs, n);

  timer = g_timer_new ();
  g_timer_start (timer);
#if 0
  multigrid (p, rhs, res, n, vcycle, npre, npost);
#else
#if 0
  for (i = 0; i < vcycle; i++)
    multigrid1 (p, rhs, res, n, npre, npost);
#else /* as implemented in gerris (speed is identical to the previous one...
         (02/07/02) */
  for (i = 0; i < vcycle; i++) {
    residual (res, p, rhs, n);
    multigrid2 (dp, res, n, npre, npost);
    addp1 (p, dp, n);
    symbc (p, n);
  }
#endif
#endif
  g_timer_stop (timer);

  residual (res, p, rhs, n);
  s = stats (res, n);
  fprintf (stderr, "div min: %10.4e avg: %10.4e|%10.4e max: %10.4e n: %7d\n",
	   s.min, s.mean, s.stddev, s.max, s.n);
  s = stats (p, n);
  fprintf (stderr, "p   min: %10.4f avg: %10.4f|%10.4f max: %10.4f n: %7d\n",
	   s.min, s.mean, s.stddev, s.max, s.n);

  fprintf (stderr, "Time: %g s Speed: %.0f sites/s\n", 
	   g_timer_elapsed (timer, NULL),
	   s.n/g_timer_elapsed (timer, NULL));

  write (res, n, stdout);

  return 0;
}
