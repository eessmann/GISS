/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
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
#include <sys/wait.h>
#include <unistd.h>
#include <math.h>
#include "utils.h"
#include "solid.h"
#include "simulation.h"

/* GfsFunction: Object */

static gboolean load_module (GfsFunction * f, GtsFile * fp, gchar * mname)
{
  gchar * path;
  
  path = g_module_build_path (GFS_MODULES_DIR, mname);
  f->module = g_module_open (path, 0);
  g_free (path);
  if (f->module == NULL)
    f->module = g_module_open (mname, 0);
  if (f->module == NULL) {
    gts_file_error (fp, "cannot load module: %s", g_module_error ());
    return FALSE;
  }
  if (!g_module_symbol (f->module, "f", (gpointer) &f->f)) {
    gts_file_error (fp, "module `%s' does not export function `f'", mname);
    g_module_close (f->module);
    return FALSE;
  }
  return TRUE;
}

static void function_read (GtsObject ** o, GtsFile * fp)
{
  GfsFunction * f = GFS_FUNCTION (*o);

  if (GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  switch (fp->type) {
    /* constant value */
  case GTS_INT: case GTS_FLOAT:
    f->val = atof (fp->token->str);
    break;

    /* load module */
  case GTS_STRING:
    if (!g_module_supported ()) {
      gts_file_error (fp, "expecting a number (val)");
      return;
    }
    else
      load_module (f, fp, fp->token->str);
    break;

    /* compile C expression */
  case '{':
    if (!g_module_supported ()) {
      gts_file_error (fp, "expecting a number (val)");
      return;
    }
    else {
#if FTT_2D
      gchar cccommand[] = "gcc `pkg-config gerris2D --cflags --libs` -O -fPIC -shared -x c";
#elif FTT_2D3
      gchar cccommand[] = "gcc `pkg-config gerris2D3 --cflags --libs` -O -fPIC -shared -x c";
#else /* 3D */
      gchar cccommand[] = "gcc `pkg-config gerris3D --cflags --libs` -O -fPIC -shared -x c";
#endif 
     gchar finname[] = "/tmp/gfsXXXXXX";
      gchar foutname[] = "/tmp/gfsXXXXXX";
      gchar ferrname[] = "/tmp/gfsXXXXXX";
      gchar ftmpname[] = "/tmp/gfsXXXXXX";
      gint find, foutd, ferrd, ftmpd;
      FILE * fin;
      guint scope;
      gchar * cc;
      gint status, c;
      GfsVariable * v;

      find = mkstemp (finname);
      if (find < 0) {
	gts_file_error (fp, "cannot create temporary file");
	return;
      }
      fin = fdopen (find, "w");
      fputs ("#include <stdlib.h>\n"
	     "#include <stdio.h>\n"
	     "#include <math.h>\n"
	     "#include <gfs.h>\n"
	     "static double Dirichlet = 1.;\n"
	     "static double Neumann = 0.;\n"
	     "double f (FttCell * cell, double x, double y, double z, double t) {\n"
	     "  double ",
	     fin);
      v = GFS_DOMAIN (gfs_object_simulation (*o))->variables;
      fprintf (fin, "%s", v->name);
      while ((v = v->next)) {
	if (v->name)
	  fprintf (fin, ", %s", v->name);
      }
      fputs (";\n  if (cell) {\n", fin);
      v = GFS_DOMAIN (gfs_object_simulation (*o))->variables;
      while (v) {
	if (v->name)
	  fprintf (fin, "    %s = GFS_VARIABLE (cell, %d);\n", v->name, v->i);
	v = v->next;
      }
      fprintf (fin, "  }\n#line %d \"GfsFunction\"\n", fp->line);
      f->expr = g_string_new ("{");
      scope = fp->scope_max;
      c = gts_file_getc (fp);
      while (c != EOF && fp->scope > scope) {
	fputc (c, fin);
	g_string_append_c (f->expr, c);
	c = gts_file_getc (fp);
      }
      fputs ("}\n", fin);
      g_string_append_c (f->expr, '}');
      fclose (fin);
      
      if (fp->scope != scope) {
	gts_file_error (fp, "parse error");
	close (find);
	remove (finname);
	return;
      }

      foutd = mkstemp (foutname);
      ferrd = mkstemp (ferrname);
      ftmpd = mkstemp (ftmpname);
      if (foutd < 0 || ferrd < 0 || ftmpd < 0) {
	gts_file_error (fp, "cannot create temporary file");
	return;
      }
      cc = g_strjoin (" ",
		      cccommand, ftmpname, 
		      "-o", foutname,
		      "`awk '{"
                      "   if ($1 == \"#\" && $2 == \"link\") {"
		      "     for (i = 3; i <= NF; i++) printf (\"%s \", $i);"
		      "     print \"\" > \"/dev/stderr\";"
                      "   }"
                      "   else if ($1 == \"#link\") {"
		      "     for (i = 2; i <= NF; i++) printf (\"%s \", $i);"
		      "     print \"\" > \"/dev/stderr\";"
		      "   } else print $0 > \"/dev/stderr\";"
		      "}' <", finname, "2>", ftmpname, "` 2>",
		      ferrname, NULL);
      status = system (cc);
      g_free (cc);
      close (find);
      remove (finname);
      close (ftmpd);
      remove (ftmpname);
      if (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || WTERMSIG (status) == SIGQUIT)) {
	close (foutd);
	remove (foutname);
	close (ferrd);
	remove (ferrname);
	exit (0);
      }
      if (status == -1 || WEXITSTATUS (status) != 0) {
	GString * msg = g_string_new ("");
	FILE * ferr = fdopen (ferrd, "r");

	while ((c = fgetc (ferr)) != EOF)
	  g_string_append_c (msg, c);
	fclose (ferr);
	gts_file_error (fp, "error compiling expression\n%s", msg->str);
	g_string_free (msg, TRUE);
	close (foutd);
	remove (foutname);
	remove (ferrname);
	return;
      }
      load_module (f, fp, foutname);
      close (foutd);
      remove (foutname);
      close (ferrd);
      remove (ferrname);
    }
    break;

  default:
    gts_file_error (fp, "expecting an expression (val)");
    return;
  }

  gts_file_next_token (fp);
}

static void function_write (GtsObject * o, FILE * fp)
{
  GfsFunction * f = GFS_FUNCTION (o);

  if (GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->write) (o, fp);

  if (f->expr)
    fprintf (fp, " %s", f->expr->str);
  else if (f->module)
    fprintf (fp, " %s", g_module_name (f->module));
  else
    fprintf (fp, " %g", f->val);
}

static void function_destroy (GtsObject * object)
{
  GfsFunction * f = GFS_FUNCTION (object);

  if (f->module) g_module_close (f->module);
  if (f->expr) g_string_free (f->expr, TRUE);

  (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->destroy) 
    (object);
}

static void gfs_function_class_init (GfsFunctionClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = function_read;
  GTS_OBJECT_CLASS (klass)->write = function_write;
  GTS_OBJECT_CLASS (klass)->destroy = function_destroy;
}

GfsFunctionClass * gfs_function_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunction",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) gfs_function_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_function_info);
  }

  return klass;
}

GfsFunction * gfs_function_new (GfsFunctionClass * klass, 
				gdouble val)
{
  GfsFunction * object;

  object = GFS_FUNCTION (gts_object_new (GTS_OBJECT_CLASS (klass)));
  object->val = val;

  return object;
}

/**
 * gfs_function_value:
 * @f: a #GfsFunction.
 * @cell: a #FttCell or %NULL.
 * @p: a #FttVector.
 * @t: the time.
 *
 * Returns: the value of function @f at location @p of @cell.
 */
gdouble gfs_function_value (GfsFunction * f, FttCell * cell, FttVector * p, gdouble t)
{
  g_return_val_if_fail (f != NULL, 0.);

  if (f->f) {
    g_return_val_if_fail (p != NULL, 0.);
    return (* f->f) (cell, p->x, p->y, p->z, t);
  }
  else
    return f->val;
}

/**
 * gfs_function_face_value:
 * @f: a #GfsFunction.
 * @fa: a #FttCellFace.
 * @t: the time.
 *
 * Returns: the value of function @f at the center of face @fa.
 */
gdouble gfs_function_face_value (GfsFunction * f, FttCellFace * fa,
				 gdouble t)
{
  g_return_val_if_fail (f != NULL, 0.);

  if (f->f) {
    FttVector p;

    g_return_val_if_fail (fa != NULL, 0.);
    
    ftt_face_pos (fa, &p);
    return (* f->f) (NULL, p.x, p.y, p.z, t);
  }
  else
    return f->val;
}

/**
 * gfs_function_read:
 * @f: a #GfsFunction.
 * @fp: a #GtsFile.
 *
 * Calls the read() method of @f.
 */
void gfs_function_read (GfsFunction * f, GtsFile * fp)
{
  GtsObject * o = (GtsObject *) f;

  g_return_if_fail (f != NULL);
  g_return_if_fail (fp != NULL);

  (* GTS_OBJECT (f)->klass->read) (&o, fp);
}

/**
 * gfs_function_write:
 * @f: a #GfsFunction.
 * @fp: a file pointer.
 *
 * Calls the write() method of @f.
 */
void gfs_function_write (GfsFunction * f, FILE * fp)
{
  g_return_if_fail (f != NULL);
  g_return_if_fail (fp != NULL);

  (* GTS_OBJECT (f)->klass->write) (GTS_OBJECT (f), fp);
}

/**
 * gfs_object_class_from_name:
 * @name: the name of the class.
 *
 * Looks for a class called @name. If none is found append the "Gfs"
 * prefix and look again.
 *
 * Returns: the class or %NULL if none was found.
 */
GtsObjectClass * gfs_object_class_from_name (const gchar * name)
{
  GtsObjectClass * klass;

  g_return_val_if_fail (name != NULL, NULL);

  klass = gts_object_class_from_name (name);
  if (klass == NULL) {
    gchar * ename = g_strconcat ("Gfs", name, NULL);
    klass = gts_object_class_from_name (ename);
    g_free (ename);
  }
  return klass;
}

static void eigsrt (gdouble d[FTT_DIMENSION], gdouble v[FTT_DIMENSION][FTT_DIMENSION])
{
  gint k, j, i;
  gdouble p;

  for (i = 0; i < FTT_DIMENSION - 1; i++) {
    p = d[k = i];

    for (j = i + 1; j < FTT_DIMENSION; j++)
      if (d[j] >= p) 
	p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < FTT_DIMENSION; j++) {
	p = v[j][i];
	v[j][i] = v[j][k];
	v[j][k] = p;
      }
    }
  }
}

#define ROTATE(a,i,j,k,l) {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}

/**
 * gfs_eigenvalues:
 * @a: a symmetric matrix.
 * @d: a vector.
 * @v: another matrix.
 *
 * Fills @d (resp. @v) with the eigenvalues (resp. eigenvectors) of
 * matrix @a.
 */
void gfs_eigenvalues (gdouble a[FTT_DIMENSION][FTT_DIMENSION],
		      gdouble d[FTT_DIMENSION],
		      gdouble v[FTT_DIMENSION][FTT_DIMENSION])
{
  gint j, iq, ip, i;
  gdouble tresh, theta, tau, t, sm, s, h, g, c, b[FTT_DIMENSION], z[FTT_DIMENSION];

  for (ip = 0; ip < FTT_DIMENSION; ip++) {
    for (iq = 0; iq < FTT_DIMENSION; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < FTT_DIMENSION; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < FTT_DIMENSION - 1; ip++) {
      for (iq = ip + 1; iq < FTT_DIMENSION; iq++)
	sm += fabs (a[ip][iq]);
    }
    if (sm == 0.0) {
      eigsrt (d, v);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(FTT_DIMENSION*FTT_DIMENSION);
    else
      tresh = 0.0;
    for (ip = 0; ip < FTT_DIMENSION - 1; ip++) {
      for (iq = ip + 1; iq < FTT_DIMENSION; iq++) {
	g = 100.0*fabs (a[ip][iq]);
	if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) && fabs(d[iq]) + g == fabs(d[iq]))
	  a[ip][iq] = 0.0;
	else if (fabs (a[ip][iq]) > tresh) {
	  h = d[iq] - d[ip];
	  if (fabs(h) + g == fabs(h))
	    t = a[ip][iq]/h;
	  else {
	    theta = 0.5*h/a[ip][iq];
	    t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt (1 + t*t);
	  s = t*c;
	  tau = s/(1.0 + c);
	  h = t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.0;
	  for (j = 0; j <= ip - 1; j++)
	    ROTATE (a, j, ip, j, iq);
	  for (j = ip + 1; j <= iq - 1; j++)
	    ROTATE (a, ip, j, j, iq);
	  for (j = iq + 1; j < FTT_DIMENSION; j++)
	    ROTATE(a, ip, j, iq, j);
	  for (j = 0; j < FTT_DIMENSION; j++)
	    ROTATE(v, j, ip, j, iq);
	}
      }
    }
    for (ip = 0; ip < FTT_DIMENSION; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  /* Too many iterations */
  for (i = 0; i < FTT_DIMENSION; i++) {
    for (j = 0; j < FTT_DIMENSION; j++)
      fprintf (stderr, "%10.3g ", a[i][j]);
    fprintf (stderr, "\n");
  }
  g_assert_not_reached ();
}

/**
 * gfs_matrix_inverse:
 * @m: a square matrix.
 * @n: size of the matrix.
 * @pivmin: the minimum value of the pivoting coefficient.
 *
 * Replaces @m with its inverse.
 *
 * Returns: %FALSE if the inversion encounters a pivot coefficient
 * smaller than or equal to @pivmin (i.e. @m is non-invertible), %TRUE
 * otherwise.
 */
gboolean gfs_matrix_inverse (gdouble ** m, guint n, gdouble pivmin)
{
  gint * indxc, * indxr, * ipiv;
  gint i, icol = 0, irow = 0, j, k, l, ll;
  gdouble big, dum, pivinv, temp;

  g_return_val_if_fail (m != NULL, FALSE);
  g_return_val_if_fail (pivmin >= 0., FALSE);

  indxc = g_malloc (sizeof (gint)*n);
  indxr = g_malloc (sizeof (gint)*n);
  ipiv = g_malloc (sizeof (gint)*n);
  
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
	for (k = 0; k < n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs (m[j][k]) >= big) {
	      big = fabs (m[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++) 
	SWAP (m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin) {
      g_free (indxc);
      g_free (indxr);
      g_free (ipiv);
      return FALSE;
    }
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
	dum = m[ll][icol];
	m[ll][icol] = 0.0;
	for (l = 0; l < n; l++)
	  m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
	SWAP (m[k][indxr[l]], m[k][indxc[l]]);
  }
  g_free (indxc);
  g_free (indxr);
  g_free (ipiv);
  return TRUE;
}

/**
 * gfs_matrix_new:
 * @n: the size of the matrix.
 * @size: the size of the matrix elements.
 *
 * The matrix elements are initialised to zero.
 *
 * Returns: a newly allocated matrix.
 */
gpointer gfs_matrix_new (guint n, guint size)
{
  guint i;
  gpointer * m, a;
  
  g_return_val_if_fail (n > 0, NULL);
  g_return_val_if_fail (size > 0, NULL);

  m = g_malloc (n*sizeof (gpointer));
  a = g_malloc0 (n*n*size);
  for (i = 0; i < n; i++)
    m[i] = GUINT_TO_POINTER (GPOINTER_TO_UINT (a) + i*n*size);
  return m;
}

/**
 * gfs_matrix_free:
 * @m: a matrix allocated with gfs_matrix_new().
 *
 * Frees the memory occupied by @m.
 */
void gfs_matrix_free (gpointer m)
{
  g_return_if_fail (m != NULL);

  g_free (((gpointer *) m)[0]);
  g_free (m);
}
