/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2004 National Institute of Water and Atmospheric
 * Research
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

#include "gfsgl.h"
#include <gerris/isocube.h>

/* GfsGlIsosurface: Object */

typedef struct {
  guint m;
  FttVector * v;
  FttVector * n;
  gdouble * sv;
  gboolean coarse_fine;
} polygon;

typedef struct {
  polygon * p;
  guint n;
} polygons;

#define param(v) (gls->max > gls->min ? ((v) - gls->min)/(gls->max - gls->min) : 0.5)
#define color(v) (gfs_colormap_color (gls->cmap, param (v)))

static void polygon_draw (polygon * p, GfsGlIsosurface * gl)
{
  guint i;

  glBegin (GL_POLYGON);
  for (i = 0; i < p->m; i++) {
    if (gl->use_scalar) {
      GfsGlScalar * gls = GFS_GL_SCALAR (gl);

      if (gfs_gl_vector_format (GFS_GL (gl))) {
	GtsColor c = color (p->sv[i]);
	glColor3f (c.r, c.g, c.b);
      }
      else
	glTexCoord1d (param (p->sv[i]));
    }
    if (gl->reversed)
      glNormal3d (-p->n[i].x, -p->n[i].y, -p->n[i].z);
    else
      glNormal3d (p->n[i].x, p->n[i].y, p->n[i].z);
    glVertex3d (p->v[i].x, p->v[i].y, p->v[i].z);
  }
  glEnd ();
}

static void polygons_draw (polygons * p, GfsGlIsosurface * gl)
{
  guint i;

  for (i = 0; i < p->n; i++)
    polygon_draw (&p->p[i], gl);
}

static void polygons_destroy (polygons * p)
{
  guint i;
  for (i = 0; i < p->n; i++) {
    g_free (p->p[i].v);
    g_free (p->p[i].n);
    g_free (p->p[i].sv);
  }
  g_free (p->p);
  g_free (p);
}

static polygons * polygons_add (polygons * p, guint n)
{
  if (p == NULL) {
    p = g_malloc (sizeof (polygons));
    p->p = g_malloc (sizeof (polygon));
    p->n = 0;
  }
  else
    p->p = g_realloc (p->p, sizeof (polygon)*(p->n + 1));
  p->p[p->n].v = g_malloc (n*sizeof (FttVector)); 
  p->p[p->n].n = g_malloc (n*sizeof (FttVector));
  p->p[p->n].sv = g_malloc (n*sizeof (gdouble));
  p->p[p->n].m = n;
  p->n++;
  return p;
}

static void free_polygons (FttCell * cell, GfsVariable * p)
{
  gpointer po = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, p));

  if (po) {
    polygons_destroy (po);
    GFS_VALUE (cell, p) = 0.;
  }
}

void gfs_gl_isosurface_reset (GfsGlIsosurface * gl)
{
  g_return_if_fail (gl != NULL);
  
  if (gl->p && GFS_GL (gl)->sim)
    gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) free_polygons, gl->p);
}

static void gl_isosurface_destroy (GtsObject * o)
{
  GfsGlIsosurface * gl = GFS_GL_ISOSURFACE (o);

  gfs_gl_var_func_destroy (gl->vf);
  g_string_free (gl->expr, TRUE);

  if (gl->min)
    gts_object_destroy (GTS_OBJECT (gl->min));
  if (gl->max)
    gts_object_destroy (GTS_OBJECT (gl->max));

  if (gl->p) {
    gfs_gl_isosurface_reset (gl);
    gts_object_destroy (GTS_OBJECT (gl->p));
  }

  (* GTS_OBJECT_CLASS (gfs_gl_isosurface_class ())->parent_class->destroy) (o);
}

static void gl_isosurface_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlIsosurface * gl = GFS_GL_ISOSURFACE (*o);
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "level",      TRUE},
    {GTS_INT,    "reversed",   TRUE},
    {GTS_INT,    "use_scalar", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_isosurface_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  g_string_free (gl->expr, TRUE);
  if (!(gl->expr = gfs_function_expression (fp, NULL)))
    return;
  gts_file_next_token (fp);

  var[0].data = &gl->level;
  var[1].data = &gl->reversed;
  var[2].data = &gl->use_scalar;
  gts_file_assign_variables (fp, var);
}

static void gl_isosurface_write (GtsObject * o, FILE * fp)
{
  GfsGlIsosurface * gl = GFS_GL_ISOSURFACE (o);

  (* GTS_OBJECT_CLASS (gfs_gl_isosurface_class ())->parent_class->write) (o, fp);

  fprintf (fp, " %s {\n"
	   "  level = %g\n"
	   "  reversed = %d\n"
	   "  use_scalar = %d\n"
	   "}",
	   gl->expr->str,
	   gl->level,
	   gl->reversed, 
	   (gl->use_scalar != NULL));
}

static void interpolated_normal (FttVector * n, GfsInterpolator * inter, GfsVariable * v)
{
  if (n->x == G_MAXDOUBLE) {
    guint i;
    
    n->x = n->y = n->z = 0.;
    for (i = 0; i < inter->n; i++) {
      gdouble w = inter->w[i];
      FttCell * c = inter->c[i];
      
      n->x += w*gfs_center_gradient (c, FTT_X, v->i);
      n->y += w*gfs_center_gradient (c, FTT_Y, v->i);
      n->z += w*gfs_center_gradient (c, FTT_Z, v->i);
    }
  }
}

static void cube_iso_intersection (FttCell * cell,
				   FttVector p[12],
				   FttVector n[12],
				   gint orient[12],
				   gdouble sv[12],
				   gdouble vc[8],
				   gdouble svc[8],
				   FttVector nc[8],
				   GfsInterpolator inter[8],
				   GfsVariable * var,
				   gdouble val,
				   GfsVariable * svar,
				   gint max_level)
{
  FttVector o;
  gdouble h = ftt_cell_size (cell)*SLIGHTLY_LARGER;
  guint i;

  for (i = 0; i < 8; i++) {
    guint j;

    gfs_cell_corner_interpolator (cell, corner[i], var->centered, max_level, &inter[i]);
    vc[i] = svc[i] = 0.;
    for (j = 0; j < inter[i].n; j++) {
      vc[i] += inter[i].w[j]*GFS_VALUE (inter[i].c[j], var);
      if (svar)
	svc[i] += inter[i].w[j]*GFS_VALUE (inter[i].c[j], svar);
    }
    nc[i].x = G_MAXDOUBLE;
  }

  ftt_cell_pos (cell, &o);
  o.x -= h/2.; o.y -= h/2.; o.z -= h/2.;
  for (i = 0; i < 12; i++) {
    guint j = edge1[i][0], k = edge1[i][1];
    if ((vc[j] >= val && vc[k] < val) || (vc[j] < val && vc[k] >= val)) {
      gdouble a = (val - vc[j])/(vc[k] - vc[j]);
      FttVector e, d;

      d.x = o.x + h*edge[i][0].x; d.y = o.y + h*edge[i][0].y; d.z = o.z + h*edge[i][0].z;
      e.x = o.x + h*edge[i][1].x; e.y = o.y + h*edge[i][1].y; e.z = o.z + h*edge[i][1].z;
      p[i].x = d.x + a*(e.x - d.x); p[i].y = d.y + a*(e.y - d.y); p[i].z = d.z + a*(e.z - d.z);
      sv[i] = svc[j] + a*(svc[k] - svc[j]);

      interpolated_normal (&nc[j], &inter[j], var);
      interpolated_normal (&nc[k], &inter[k], var);
      n[i].x = nc[j].x + a*(nc[k].x - nc[j].x);
      n[i].y = nc[j].y + a*(nc[k].y - nc[j].y); 
      n[i].z = nc[j].z + a*(nc[k].z - nc[j].z);
      orient[i] = (vc[j] >= val);
    }
    else
      orient[i] = -1;
  }
}

static gboolean face_is_coarse_fine (FttCell * cell, guint e, gint orient, gint max_depth)
{
  FttCell * n = ftt_cell_neighbor (cell, connect[e][orient][3]);

  if (n && !FTT_CELL_IS_LEAF (n)) {
    guint l = ftt_cell_level (n);
    return (l != max_depth && l == ftt_cell_level (cell));
  }
  return FALSE;
}

static gint next_vertex (guint s, 
			 gdouble v0, gdouble v1, gdouble v2, gdouble v3, 
			 FttVector p[4], FttVector n[4], gdouble sv[4],
			 gdouble val,
			 FttVector * p1, FttVector * n1, gdouble * sv1)
{
  gdouble v[4];
  guint i;
  v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
  for (i = 1; i < 4; i++) {
    guint k1 = (s + i) % 4;
    guint k2 = (s + i + 1) % 4;
    if ((v[k1] >= val && v[k2] < val) || (v[k1] < val && v[k2] >= val)) {
      gdouble a = (val - v[k1])/(v[k2] - v[k1]);
      p1->x = p[k1].x + a*(p[k2].x - p[k1].x);
      p1->y = p[k1].y + a*(p[k2].y - p[k1].y);
      p1->z = p[k1].z + a*(p[k2].z - p[k1].z);
      n1->x = n[k1].x + a*(n[k2].x - n[k1].x);
      n1->y = n[k1].y + a*(n[k2].y - n[k1].y);
      n1->z = n[k1].z + a*(n[k2].z - n[k1].z);
      *sv1 = sv[k1] + a*(sv[k2] - sv[k1]);
      return k1;
    }
  }
  return -1;
}

/* see doc/coarse_fine.fig */
static gint next_square (guint square, guint edge, 
			 gdouble v[9], gdouble sv[9], gdouble val,
			 FttVector p[9], FttVector n[9],
			 FttVector v3[24], FttVector n3[24], gdouble sv3[24],
			 guint * nv)
{
  static guint sqv[4][4] = {
    {3, 4, 8, 7}, {4, 0, 5, 8}, {7, 8, 6, 2}, {8, 5, 1, 6}
  };
  static gint sqn[4][4] = {
    {-1, 1, 2, -4}, {-1, -2, 3, 0}, {0, 3, -3, -4}, {1, -2, -3, 2}
  };
  static gint edn[4][4] = {
    {0, 3, 0, 0}, {0, 0, 0, 1}, {2, 3, 0, 0}, {2, 0, 0, 1}
  };
  FttVector p2[4], n2[4], p1, n1;  
  gdouble sv2[4], sv1;
  gint u, nsquare, nedge;

  p2[0] = p[sqv[square][0]];
  p2[1] = p[sqv[square][1]];
  p2[2] = p[sqv[square][2]];
  p2[3] = p[sqv[square][3]];

  n2[0] = n[sqv[square][0]];
  n2[1] = n[sqv[square][1]];
  n2[2] = n[sqv[square][2]];
  n2[3] = n[sqv[square][3]];

  sv2[0] = sv[sqv[square][0]];
  sv2[1] = sv[sqv[square][1]];
  sv2[2] = sv[sqv[square][2]];
  sv2[3] = sv[sqv[square][3]];

  u = next_vertex (edge,
		   v[sqv[square][0]], v[sqv[square][1]],
		   v[sqv[square][2]], v[sqv[square][3]], p2, n2, sv2, val, &p1, &n1, &sv1);
  if (u < 0) {
    *nv = 0;
    return -1;
  }
  nsquare = sqn[square][u];
  if (nsquare < 0)
    return - nsquare - 1;
  n3[*nv] = n1;
  sv3[*nv] = sv1;
  v3[(*nv)++] = p1;
  g_assert (*nv <= 24); 
  nedge = edn[square][u];
  return next_square (nsquare, nedge, v, sv, val, p, n, v3, n3, sv3, nv);
}

static void gl_isosurface (FttCell * cell, GfsGl * gl)
{
  GfsGlIsosurface * gls = GFS_GL_ISOSURFACE (gl);
  polygons * p = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gls->p));

  if (p) {
    polygons_draw (p, gls);
    gl->size++;
  }
  else {
    gboolean coarse_fine = FALSE;
    
    {
      FttVector _p[12], _n[12], v[24], n[24];
      gint _orient[12];
      gdouble sv[24], _sv[12], _svc[8], _vc[8], _h = ftt_cell_size (cell)*SLIGHTLY_LARGER;
      FttVector _nc[8];
      GfsInterpolator _inter[8];
      guint _i;
      cube_iso_intersection (cell, _p, _n, _orient, _sv, _vc, _svc, _nc, _inter, 
			     gls->v, gls->level,
			     gls->use_scalar ? GFS_GL_SCALAR (gls)->v : NULL,
			     GFS_GL (gl)->maxlevel);
      for (_i = 0; _i < 12; _i++) {
	guint nv = 0, _e = _i;
	while (_orient[_e] >= 0) {
	  guint _m = 0, * _ne = connect[_e][_orient[_e]];
	  n[nv] = _n[_e];
	  sv[nv] = _sv[_e];
	  v[nv++] = _p[_e];
	  if (face_is_coarse_fine (cell, _e, _orient[_e], gl->maxlevel)) {
	    guint _v0 = connectv[_e][_orient[_e]][0];
	    guint _v1 = connectv[_e][_orient[_e]][1];
	    guint _v2 = connectv[_e][_orient[_e]][2];
	    guint _v3 = connectv[_e][_orient[_e]][3];
	    FttVector _o, _p2[9], _n2[9];
	    gdouble _v[9], _vs[9];
	    guint _q, j;

	    _v[0] = _vc[_v0]; _v[1] = _vc[_v1]; _v[2] = _vc[_v2]; _v[3] = _vc[_v3];
	    _v[4] = (_v[3] + _v[0])/2.;
	    _v[5] = (_v[1] + _v[0])/2.;
	    _v[6] = (_v[1] + _v[2])/2.;
	    _v[7] = (_v[3] + _v[2])/2.;
	    _v[8] = (_v[0] + _v[1] + _v[2] + _v[3])/4.;

	    _vs[0] = _svc[_v0]; _vs[1] = _svc[_v1]; _vs[2] = _svc[_v2]; _vs[3] = _svc[_v3];
	    _vs[4] = (_vs[3] + _vs[0])/2.;
	    _vs[5] = (_vs[1] + _vs[0])/2.;
	    _vs[6] = (_vs[1] + _vs[2])/2.;
	    _vs[7] = (_vs[3] + _vs[2])/2.;
	    _vs[8] = (_vs[0] + _vs[1] + _vs[2] + _vs[3])/4.;

	    interpolated_normal (&_nc[_v0], &_inter[_v0], gls->v);
	    _n2[0] = _nc[_v0];
	    interpolated_normal (&_nc[_v1], &_inter[_v1], gls->v);
	    _n2[1] = _nc[_v1];
	    interpolated_normal (&_nc[_v2], &_inter[_v2], gls->v);
	    _n2[2] = _nc[_v2];
	    interpolated_normal (&_nc[_v3], &_inter[_v3], gls->v);
	    _n2[3] = _nc[_v3];
	    ftt_cell_pos (cell, &_o);
	    _o.x -= _h/2.; _o.y -= _h/2.; _o.z -= _h/2;
	    for (j = 0; j < 3; j++) {
	      (&_p2[0].x)[j] = (&_o.x)[j] + _h*(&cvertex[_v0].x)[j];
	      (&_p2[1].x)[j] = (&_o.x)[j] + _h*(&cvertex[_v1].x)[j];
	      (&_p2[2].x)[j] = (&_o.x)[j] + _h*(&cvertex[_v2].x)[j];
	      (&_p2[3].x)[j] = (&_o.x)[j] + _h*(&cvertex[_v3].x)[j];
	      (&_p2[4].x)[j] = ((&_p2[3].x)[j] + (&_p2[0].x)[j])/2.;
	      (&_p2[5].x)[j] = ((&_p2[1].x)[j] + (&_p2[0].x)[j])/2.;
	      (&_p2[6].x)[j] = ((&_p2[1].x)[j] + (&_p2[2].x)[j])/2.;
	      (&_p2[7].x)[j] = ((&_p2[3].x)[j] + (&_p2[2].x)[j])/2.;
	      (&_p2[8].x)[j] = ((&_p2[0].x)[j] + (&_p2[1].x)[j] + 
				(&_p2[2].x)[j] + (&_p2[3].x)[j])/4.;
	      (&_n2[4].x)[j] = ((&_n2[3].x)[j] + (&_n2[0].x)[j])/2.;
	      (&_n2[5].x)[j] = ((&_n2[1].x)[j] + (&_n2[0].x)[j])/2.;
	      (&_n2[6].x)[j] = ((&_n2[1].x)[j] + (&_n2[2].x)[j])/2.;
	      (&_n2[7].x)[j] = ((&_n2[3].x)[j] + (&_n2[2].x)[j])/2.;
	      (&_n2[8].x)[j] = ((&_n2[0].x)[j] + (&_n2[1].x)[j] + 
				(&_n2[2].x)[j] + (&_n2[3].x)[j])/4.;
	    }
	    _q = next_square ((gls->level - _v[3])/(_v[0] - _v[3]) > 0.5, 0, _v, _vs, gls->level,
			      _p2, _n2, v, n, sv, &nv);
	    if (_q <= 0)
	      break;
	    _orient[_e] = -1;
	    _e = _ne[_q - 1];
	    coarse_fine = FALSE;
	  }
	  else {
	    _orient[_e] = -1;
	    while (_m < 3 && _orient[_e] < 0)
	      _e = _ne[_m++];
	  }
	}
	if (nv > 2) {
	  guint i;
	  polygon * q;

	  p = polygons_add (p, nv);
	  q = &p->p[p->n - 1];

	  for (i = 0; i < nv; i++) {
	    q->v[i] = v[i];
	    q->n[i] = n[i];
	    q->sv[i] = sv[i];
	  }
	  q->coarse_fine = coarse_fine;
	}
	coarse_fine = FALSE;
      }
    }
    if (p) {
      polygons_draw (p, gls);
      GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gls->p)) = p;
      gl->size++;
    }
  }
}

static void gl_isosurface_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlIsosurface * gls = GFS_GL_ISOSURFACE (gl);
  
  if (gls->use_scalar && gls->use_scalar != GFS_GL_SCALAR (gl)->v) {
    gls->use_scalar = GFS_GL_SCALAR (gl)->v;
    gfs_gl_isosurface_reset (gls);
  }

  gl->size = 0;
  glShadeModel (GL_SMOOTH);
  if (gls->use_scalar && !gfs_gl_vector_format (gl)) {
    glEnable (GL_TEXTURE_1D);
    gfs_colormap_texture (GFS_GL_SCALAR (gl)->cmap);
    glColor3f (1., 1., 1.);
    gfs_gl_cell_traverse_visible_iso (gl, f, gls->min, gls->max, gls->level,
				      (FttCellTraverseFunc) gl_isosurface, gl);
    glDisable (GL_TEXTURE_1D);
  }
  else
    gfs_gl_cell_traverse_visible_iso (gl, f, gls->min, gls->max, gls->level,
				      (FttCellTraverseFunc) gl_isosurface, gl);

  if (gls->use_scalar)
    (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass->parent_class)->draw) (gl, f);  
}

static void min_max_iso (FttCell * cell, GfsGlIsosurface * gl)
{
  gdouble min = G_MAXDOUBLE, max = -G_MAXDOUBLE;
  guint i;
  gpointer p;

  if (FTT_CELL_IS_LEAF (cell))
    for (i = 0; i < 8; i++) {
      gdouble v = gfs_cell_corner_value (cell, corner[i], gl->v, -1);
      if (v < min) min = v;
      if (v > max) max = v;
    }
  else {
    FttCellChildren child;

    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i]) {
	gdouble vmin = GFS_VALUE (child.c[i], gl->min);
	gdouble vmax = GFS_VALUE (child.c[i], gl->max);
	if (vmin < min) min = vmin;
	if (vmax > max) max = vmax;
      }
  }
  GFS_VALUE (cell, gl->min) = min;
  GFS_VALUE (cell, gl->max) = max;
  if ((p = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gl->p)))) {
    polygons_destroy (p);
    GFS_VALUE (cell, gl->p) = 0.;
  }
}

static void min_max_iso1 (GfsBox * b, GfsGlIsosurface * gl)
{
  gdouble min = GFS_VALUE (b->root, gl->min);
  gdouble max = GFS_VALUE (b->root, gl->max);
  if (min < gl->minv) gl->minv = min;
  if (max > gl->maxv) gl->maxv = max;
}

GtsFile * gfs_gl_isosurface_set (GfsGlIsosurface * gl, const gchar * func)
{
  GtsFile * fp;

  g_return_val_if_fail (gl != NULL, NULL);
  g_return_val_if_fail (func != NULL, NULL);

  if ((fp = gfs_gl_var_func_set (gl->vf, GFS_GL (gl)->sim, func, gl->expr, NULL)))
    return fp;

  gl->v = gl->vf->v;
  gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim),
			    FTT_POST_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) min_max_iso, gl);
  gl->minv = G_MAXDOUBLE;
  gl->maxv = -G_MAXDOUBLE;
  gts_container_foreach (GTS_CONTAINER (GFS_GL (gl)->sim), (GtsFunc) min_max_iso1, gl);

  return NULL;
}

static void gl_isosurface_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsGlIsosurface * gls = GFS_GL_ISOSURFACE (object);
  GtsFile * fp = NULL;
  
  gfs_gl_isosurface_reset (gls);

  if (GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_isosurface_class ())->parent_class)->set_simulation)
    (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_isosurface_class ())->parent_class)->set_simulation)
      (object, sim);

  if (gls->min)
    gts_object_destroy (GTS_OBJECT (gls->min));
  gls->min = gfs_temporary_variable (domain);
  if (gls->max)
    gts_object_destroy (GTS_OBJECT (gls->max));
  gls->max = gfs_temporary_variable (domain);
  if (gls->p)
    gts_object_destroy (GTS_OBJECT (gls->p));
  gls->p = gfs_temporary_variable (domain);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, gls->p);

  if (gls->expr->str[0] == '\0' || (fp = gfs_gl_isosurface_set (gls, gls->expr->str))) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    
    if (domain->variables) {
      GfsVariable * v = domain->variables->data;
      gfs_gl_isosurface_set (gls, v->name);
    }
    else
      gfs_gl_isosurface_set (gls, "0");
  }
  if (fp)
    gts_file_destroy (fp);
}

static void gl_isosurface_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_isosurface_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl_isosurface_read;
  GTS_OBJECT_CLASS (klass)->write = gl_isosurface_write;
  klass->set_simulation = gl_isosurface_set_simulation;
  klass->draw = gl_isosurface_draw;
  klass->pick = NULL;
}

static void gl_isosurface_init (GfsGlIsosurface * gl)
{
  GtsColor c = { 1., 1., 1. };

  gl->expr = g_string_new ("");
  gl->vf = gfs_gl_var_func_new ();
  GFS_GL (gl)->lc = c;
}

GfsGlClass * gfs_gl_isosurface_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_isosurface_info = {
      "GfsGlIsosurface",
      sizeof (GfsGlIsosurface),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_isosurface_class_init,
      (GtsObjectInitFunc) gl_isosurface_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_scalar_class ()),
				  &gfs_gl_isosurface_info);
  }

  return klass;
}

