/* Gerris - The GNU Flow Solver
 * Copyright (C) 2007 National Institute of Water and Atmospheric Research
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
#include "cartesian.h"

/* GfsCartesianGrid: Object */

static void gfs_cartesian_grid_read (GtsObject ** o, GtsFile * fp)
{
  GfsCartesianGrid * cgd = GFS_CARTESIAN_GRID (*o);
  guint i, j, taille = 1;

  /* call read method of parent */
  if (GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  /* do object-specific read here */

  while (fp->type == '\n') 
    gts_file_next_token (fp);
  if (fp->type != GTS_INT) {
     gts_file_error (fp, "expecting an integer (N)");
     return;
  }
  cgd->N = atoi (fp->token->str);
  gts_file_next_token (fp);

  cgd->name = g_malloc0 ((cgd->N + 1)*sizeof(char*));

  for (i = 0; i < cgd->N + 1; i++) {
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (name[%d])", i);
      return;
    } 
    cgd->name[i] = g_strdup (fp->token->str);
    gts_file_next_token (fp);
  }

  cgd->n = g_malloc(cgd->N*sizeof(guint));
  
  for (i = 0; i < cgd->N; i++) {
    while (fp->type == '\n') 
      gts_file_next_token (fp);
    if (fp->type != GTS_INT) {
      gts_file_error (fp, "expecting an integer (n[%d])", i);
      return;
    }
    
    cgd->n[i] = atoi (fp->token->str);
    gts_file_next_token (fp);
    taille *= cgd->n[i];
  }

  cgd->x = g_malloc0 (cgd->N*sizeof (gint));
  
  for (i = 0; i < cgd->N; i++) {
    cgd->x[i] = g_malloc (cgd->n[i]*sizeof (gdouble));
    for (j = 0; j < cgd->n[i]; j++) {
      if (fp->type == '\n')
	gts_file_next_token (fp);
      if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
        gts_file_error (fp, "expecting a number (x[%d][%d])", i, j);
        return;
      }
      cgd->x[i][j] = atof (fp->token->str);
      gts_file_next_token (fp);
    }
  }

  cgd->v = g_malloc (taille*sizeof (gdouble));
  
  for (i = 0; i < taille; i++) {
    if (fp->type == '\n')
      gts_file_next_token (fp);
    if (fp->type != GTS_FLOAT && fp->type != GTS_INT) {
      gts_file_error (fp, "expecting a double");
      return;
    }
    cgd->v[i] = atof (fp->token->str);
    gts_file_next_token (fp);
  }
}

static void gfs_cartesian_grid_write (GtsObject * o, FILE * fp)
{
  GfsCartesianGrid * cgd = GFS_CARTESIAN_GRID (o);
  guint i, j, taille = 1;

  /* call write method of parent */
  if (GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->write) 
      (o, fp);

  /* do object specific write here */
  for (i = 0; i < cgd->N; i++)
    taille *= cgd->n[i];

  fprintf (fp, "%d ", cgd->N);
  for (i = 0; i < cgd->N+1; i++)
    fprintf (fp, "%s ", cgd->name[i]);
  fputc ('\n', fp);
  for (i=0;i<cgd->N;i++)
    fprintf (fp,"%d\n",cgd->n[i]);

  for (i = 0; i < cgd->N; i++)
    for (j = 0; j < cgd->n[i]; j++)
      fprintf (fp,"%f\n",cgd->x[i][j]);

  for (i = 0; i < taille; i++)
    fprintf (fp,"%f\n", cgd->v[i]);  
}

static void gfs_cartesian_grid_destroy (GtsObject * object)
{
  /* do object-specific cleanup here */
  GfsCartesianGrid * cgd = GFS_CARTESIAN_GRID (object);  

  guint i;
  if (cgd->name) {
    for (i = 0; i < cgd->N+1; i++)
      g_free (cgd->name[i]);
    g_free (cgd->name);
  }
  g_free (cgd->n);
  if (cgd->x) {
    for (i = 0; i < cgd->N; i++)
      g_free (cgd->x[i]);
    g_free (cgd->x);
  }
  g_free (cgd->v);
 
  /* do not forget to call destroy method of the parent */
  (* GTS_OBJECT_CLASS (gfs_cartesian_grid_class ())->parent_class->destroy) (object);
}

static void gfs_cartesian_grid_class_init (GtsObjectClass * klass)
{
  /* define new methods and overload inherited methods here */
  GTS_OBJECT_CLASS (klass)->read = gfs_cartesian_grid_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_cartesian_grid_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_cartesian_grid_destroy;
}

GtsObjectClass * gfs_cartesian_grid_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_cartesian_grid_info = {
      "GfsCartesianGrid",
      sizeof (GfsCartesianGrid),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_cartesian_grid_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
 			  &gfs_cartesian_grid_info);
  }

  return klass;
}

GfsCartesianGrid * gfs_cartesian_grid_new (GtsObjectClass * klass)
{
  GfsCartesianGrid * object;

  object = GFS_CARTESIAN_GRID (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

static void slice (GfsCartesianGrid * g, guint p, GfsCartesianGrid * s)
{
  s->N = g->N - 1;
  s->n = &g->n[1];
  s->x = &g->x[1];
  guint i;
  gulong size = 1;
  for (i = 1; i < g->N; i++)
    size *= g->n[i];
  s->v = &g->v[size*p];
}

static gint lookup (GfsCartesianGrid * g, gdouble x)
{
  guint min = 0, max = g->n[0] - 1;
  if (x < g->x[0][min] || x > g->x[0][max])
    return -1;
  while (max > min + 1) {
    guint n = (min + max)/2;
    if (x > g->x[0][n])
      min = n;
    else
      max = n;
  }
  return min;
}

/**
 * gfs_cartesian_grid_interpolate:
 * @g: a Cartesian grid.
 * @p: a position vector of dimension @g->N.
 * @val: the interpolated value at position @p.
 *
 * Returns: %TRUE if @val has been computed, %FALSE if @p is not
 * contained within @g.
 */
gboolean gfs_cartesian_grid_interpolate (GfsCartesianGrid * g, gdouble * p, gdouble * val)
{
  g_return_val_if_fail (g != NULL, FALSE);
  g_return_val_if_fail (g->N > 0, FALSE);
  g_return_val_if_fail (p != NULL, FALSE);
  g_return_val_if_fail (val != NULL, FALSE);

  gint i = lookup (g, p[0]);
  if (i < 0)
    return FALSE;
  gdouble v1, v2;
  if (g->N > 1) {
    GfsCartesianGrid g1;
    slice (g, i, &g1);
    if (!gfs_cartesian_grid_interpolate (&g1, &p[1], &v1))
      return FALSE;
    slice (g, i + 1, &g1);
    if (!gfs_cartesian_grid_interpolate (&g1, &p[1], &v2))
      return FALSE;
  }
  else {
    v1 = g->v[i];
    v2 = g->v[i + 1];
  }

  g_assert (g->x[0][i + 1] -  g->x[0][i] != 0.);
  *val = v1 + (v2 - v1)*(p[0] - g->x[0][i])/(g->x[0][i + 1] -  g->x[0][i]);
  return TRUE;
}
