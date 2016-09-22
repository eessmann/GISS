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

static gboolean plane_cuts_cell (FttCell * cell, gpointer data)
{
  GfsGl2D * gl = data;
  FttVector p;
  ftt_cell_pos (cell, &p);
  gdouble r = ftt_cell_size (cell)*SLIGHTLY_LARGER*GFS_DIAGONAL;
  return (fabs ((p.x - gl->p[0].x)*gl->n.x + 
		(p.y - gl->p[0].y)*gl->n.y + 
		(p.z - gl->p[0].z)*gl->n.z)
	  <= r);
}

/**
 * gfs_gl_cell_traverse_visible_plane:
 * @gl: a #GfsGl2D.
 * @f: a view frustum.
 * @func: a used-defined function.
 * @data: user data to pass to @func.
 *
 * Traverse the cells of @gl which are visible and whose bounding
 * sphere is intersected by the plane defined by @gl.
 */
static void gfs_gl_cell_traverse_visible_plane (GfsGl * gl,
						GfsFrustum * f,
						FttCellTraverseFunc func,
						gpointer data)
{
  g_return_if_fail (gl != NULL);
  g_return_if_fail (f != NULL);
  g_return_if_fail (func != NULL);

  gfs_gl_cell_traverse_visible_condition (gl, f, plane_cuts_cell, gl, func, data);
}

/* GfsGlCells: Object */

static void gl_cell (FttCell * cell, GfsGl * gl)
{
  FttVector v[12];
  FttDirection d[12];
  guint nv = gfs_cut_cube_vertices (cell, gl->maxlevel,
				    &GFS_GL2D (gl)->p[0], &GFS_GL2D (gl)->n,
				    v, d,
				    NULL, NULL);

  if (nv > 2) {
    guint i;
    glBegin (GL_LINE_LOOP);
    for (i = 0; i < nv; i++)
      glVertex3d (v[i].x, v[i].y, v[i].z);
    glEnd ();
    gl->size++;
  }
}

/* GfsGlSquares: Object */

static void gl_square (FttCell * cell, GfsGl * gl)
{
  if (!GFS_HAS_DATA (cell, GFS_GL_SCALAR (gl)->v))
    return;
  FttVector v[12];
  FttDirection d[12];
  guint nv = gfs_cut_cube_vertices (cell, gl->maxlevel,
				    &GFS_GL2D (gl)->p[0], &GFS_GL2D (gl)->n,
				    v, d,
				    NULL, NULL);
  if (nv > 2) {
    GfsGlScalar * gls = GFS_GL_SCALAR (gl);
    GtsColor c;
    c = gfs_colormap_color (gls->cmap, gls->max > gls->min ?
			    (GFS_VALUE (cell, gls->v) - gls->min)/(gls->max - gls->min) :
			    0.5);
    glColor3f (c.r, c.g, c.b);
    guint i;
    glBegin (GL_POLYGON);
    for (i = 0; i < nv; i++)
      glVertex3d (v[i].x, v[i].y, v[i].z);
    glEnd ();
    gl->size++;
  }
}

static void gl_squares_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_square, gl);

  (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass->parent_class)->draw) (gl, f);
}

/* GfsGlLinear: Object */

#define param(v) (gls->max > gls->min ? ((v) - gls->min)/(gls->max - gls->min) : 0.5)
#define color(v) (gfs_colormap_color (gls->cmap, param (v)))

static void gl_linear_color (FttCell * cell, GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  if (!GFS_HAS_DATA (cell, gls->v))
    return;
  FttVector v[12];
  FttDirection d[12];
  gdouble val[12];
  guint nv = gfs_cut_cube_vertices (cell, gl->maxlevel,
				    &GFS_GL2D (gl)->p[0], &GFS_GL2D (gl)->n,
				    v, d,
				    gls->v, val);

  if (nv > 2) {
    guint i;

    glBegin (GL_POLYGON);
    for (i = 0; i < nv; i++) {
      GtsColor c = color (val[i]);
      glColor3f (c.r, c.g, c.b);
      glVertex3d (v[i].x, v[i].y, v[i].z);
    }
    glEnd ();
    gl->size++;
  }
}

static void gl_linear_texture (FttCell * cell, GfsGl * gl)
{
  FttVector v[12];
  FttDirection d[12];
  gdouble val[12];
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  guint nv = gfs_cut_cube_vertices (cell, gl->maxlevel,
				    &GFS_GL2D (gl)->p[0], &GFS_GL2D (gl)->n,
				    v, d,
				    gls->v, val);

  if (nv > 2) {
    guint i;
    FttVector c;
    gdouble vc = 0.;

    c.x = c.y = c.z = 0.;
    for (i = 0; i < nv; i++) {
      c.x += v[i].x; c.y += v[i].y; c.z += v[i].z;
      vc += val[i];
    }

    glBegin (GL_TRIANGLE_FAN);
    glTexCoord1d (param (vc/nv));
    glVertex3d (c.x/nv, c.y/nv, c.z/nv);
    for (i = 0; i < nv; i++) {
      glTexCoord1d (param (val[i]));
      glVertex3d (v[i].x, v[i].y, v[i].z);
    }
    glTexCoord1d (param (val[0]));
    glVertex3d (v[0].x, v[0].y, v[0].z);
    glEnd ();
    gl->size++;
  }
}

static void gl_linear_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glShadeModel (GL_SMOOTH);
  gfs_gl_normal (gl);
  if (gfs_gl_vector_format (gl))
    gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_linear_color, gl);
  else {
    glEnable (GL_TEXTURE_1D);
    gfs_colormap_texture (GFS_GL_SCALAR (gl)->cmap);
    glColor3f (1., 1., 1.);
    gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_linear_texture, gl);
    glDisable (GL_TEXTURE_1D);
  }

  (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass->parent_class)->draw) (gl, f);
}

/* GfsGlIsoline: Object */

static void gl_isoline (FttCell * cell, GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  if (!GFS_HAS_DATA (cell, gls->v))
    return;
  FttVector v[12];
  FttDirection d[12];
  gdouble val[12];
  guint nv = gfs_cut_cube_vertices (cell, gl->maxlevel,
				    &GFS_GL2D (gl)->p[0], &GFS_GL2D (gl)->n,
				    v, d,
				    gls->v, val);

  if (nv > 2) {
    GArray * levels = GFS_GL_ISOLINE (gl)->levels;
    guint i;

    for (i = 0; i < levels->len; i++) {
      gdouble z = g_array_index (levels, gdouble, i);
      guint j;
      guint n = 0;

      val[nv] = val[0];
      v[nv] = v[0];
      for (j = 0; j < nv; j++)
	if ((val[j] > z && val[j + 1] <= z) || (val[j] <= z && val[j + 1] > z)) {
	  gdouble a = (z - val[j])/(val[j + 1] - val[j]);
	  glVertex3d (v[j].x + a*(v[j + 1].x - v[j].x),
		      v[j].y + a*(v[j + 1].y - v[j].y),
		      v[j].z + a*(v[j + 1].z - v[j].z));
	  n++;
	}
      g_assert (n % 2 == 0);
    }
    gl->size++;
  }
}

static void gl_isoline_draw (GfsGl * gl, GfsFrustum * f)
{
  gl_isoline_update_levels (gl);
  
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glBegin (GL_LINES);
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_isoline, gl);
  glEnd ();
  glPopMatrix ();
}

/* GfsGlVOF: Object */

static void gl_vof (FttCell * cell, GfsGl * gl)
{
  FttVector v[12], m;
  guint nv = gfs_vof_facet (cell, GFS_VARIABLE_TRACER_VOF (GFS_GL_VOF (gl)->vf->v), v, &m);

  if (nv > 2) {
    guint i;
    glBegin (GL_POLYGON);
    if (GFS_GL_VOF (gl)->use_scalar) {
      GfsGlScalar * gls = GFS_GL_SCALAR (gl);
      gdouble value;
      GtsColor c;

      if (GFS_GL_VOF (gl)->interpolate) {
	FttVector c;
	gfs_vof_center (cell, GFS_VARIABLE_TRACER_VOF (GFS_GL_VOF (gl)->vf->v), &c);
	value = gfs_interpolate (cell, c, gls->v);
      }
      else
	value = GFS_VALUE (cell, gls->v);
      c = gfs_colormap_color (gls->cmap, gls->max > gls->min ?
			      (value - gls->min)/(gls->max - gls->min) : 0.5);
      glColor3f (c.r, c.g, c.b);
    }
    if (GFS_GL_VOF (gl)->reversed)
      glNormal3d (-m.x, -m.y, -m.z);
    else
      glNormal3d (m.x, m.y, m.z);
    for (i = 0; i < nv; i++)
      glVertex3d (v[i].x, v[i].y, v[i].z);
    glEnd ();
    gl->size++;
  }
}

static void gl_vof_edges (FttCell * cell, GfsGl * gl)
{
  FttVector v[12], m;
  guint nv = gfs_vof_facet (cell, GFS_VARIABLE_TRACER_VOF (GFS_GL_VOF (gl)->vf->v), v, &m);

  if (nv > 2) {
    guint i;
    glBegin (GL_LINE_LOOP);
    for (i = 0; i < nv; i++)
      glVertex3d (v[i].x, v[i].y, v[i].z);
    glEnd ();
  }
}

static gboolean is_vof (FttCell * cell, gpointer data)
{
  gdouble f = GFS_VALUE (cell, GFS_GL_VOF (data)->vf->v);

  return !GFS_IS_FULL (f);
}

static void gl_vof_draw (GfsGl * gl, GfsFrustum * f)
{
  if (GFS_IS_VARIABLE_TRACER_VOF (GFS_GL_VOF (gl)->vf->v)) {
    gl->size = 0;
    gfs_gl_cell_traverse_visible_condition (gl, f, is_vof, gl, (FttCellTraverseFunc) gl_vof, gl);
    if (GFS_GL_VOF (gl)->draw_edges) {
      glMatrixMode (GL_PROJECTION);
      glPushMatrix ();
      glTranslatef (0., 0., gl->p->lc);
      glColor3f (0., 0., 0.);
      gfs_gl_cell_traverse_visible_condition (gl, f, is_vof, gl,
					      (FttCellTraverseFunc) gl_vof_edges, gl);
      glPopMatrix ();
    }
  }
}

static gboolean plane_segment_intersection (GfsGl2D * plane, FttVector * A, FttVector * B,
					    FttVector * I)
{
  GtsVector AB, AP0;
  gts_vector_init (AB, A, B);
  gts_vector_init (AP0, A, &plane->p[0]);
  gdouble lambda = gts_vector_scalar (AB, &plane->n.x);
  if (lambda != 0.) {
    lambda = gts_vector_scalar (AP0, &plane->n.x)/lambda;
    if (lambda >= 0. && lambda < 1.) {
      I->x = A->x + lambda*AB[0];
      I->y = A->y + lambda*AB[1];
      I->z = A->z + lambda*AB[2];
      return TRUE;
    }
  }
  return FALSE;
}

static void gl_vof_cut (GfsGl * gl, FttCell * cell, GfsGl2D * plane)
{
  FttVector v[12], m;
  guint nv = gfs_vof_facet (cell, GFS_VARIABLE_TRACER_VOF (GFS_GL_VOF (gl)->vf->v), v, &m);

  if (nv > 2) {
    FttVector I[3];
    guint i, ni = 0;
    for (i = 0; i < nv - 1 && ni < 3; i++)
      if (plane_segment_intersection (plane, &v[i], &v[i+1], &I[ni]))
	ni++;
    if (plane_segment_intersection (plane, &v[nv - 1], &v[0], &I[ni]))
      ni++;
    if (ni == 2) {
      glVertex3d (I[0].x, I[0].y, I[0].z);
      glVertex3d (I[1].x, I[1].y, I[1].z);
    }
    GFS_GL (plane)->size++;
  }
}

/* GfsGlSolid: Object */

typedef struct {
  guint m;
  FttVector * v;
  FttVector * n;
  gdouble * val;
} polygon;

typedef struct {
  polygon * p;
  guint n;
} polygons;

#define param(v) (gls->max > gls->min ? ((v) - gls->min)/(gls->max - gls->min) : 0.5)
#define color(v) (gfs_colormap_color (gls->cmap, param (v)))

static void polygon_draw (polygon * p, GfsGlSolid * gl)
{
  guint i;

  glBegin (GL_POLYGON);
  for (i = 0; i < p->m; i++) {
    if (gl->use_scalar) {
      GfsGlScalar * gls = GFS_GL_SCALAR (gl);

      if (gfs_gl_vector_format (GFS_GL (gl))) {
	GtsColor c = color (p->val[i]);
	glColor3f (c.r, c.g, c.b);
      }
      else
	glTexCoord1d (param (p->val[i]));
    }
    if (gl->reversed)
      glNormal3d (-p->n[i].x, -p->n[i].y, -p->n[i].z);
    else
      glNormal3d (p->n[i].x, p->n[i].y, p->n[i].z);
    glVertex3d (p->v[i].x, p->v[i].y, p->v[i].z);
  }
  glEnd ();
}

static void polygons_draw (polygons * p, GfsGlSolid * gl)
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
    g_free (p->p[i].val);
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
  p->p[p->n].val = g_malloc (n*sizeof (gdouble));
  p->p[p->n].m = n;
  p->n++;
  return p;
}

static void free_polygons (FttCell * cell, GfsGlSolid * gl)
{
  gpointer po = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gl->p));
  gpointer s = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gl->s));

  if (po) {
    polygons_destroy (po);
    GFS_VALUE (cell, gl->p) = 0.;
  }
  if (s) {
    gts_object_destroy (s);
    GFS_VALUE (cell, gl->s) = 0.;
  }
}

void gfs_gl_solid_reset (GfsGlSolid * gl)
{
  g_return_if_fail (gl != NULL);
  
  if (gl->p && gl->s && GFS_GL (gl)->sim)
    gfs_domain_cell_traverse (GFS_DOMAIN (GFS_GL (gl)->sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			      (FttCellTraverseFunc) free_polygons, gl);
  gl->needs_updating = TRUE;
}

static void gl_solid_destroy (GtsObject * o)
{
  GfsGlSolid * gl = GFS_GL_SOLID (o);

  gfs_gl_solid_reset (gl);
  if (gl->p)
    gts_object_destroy (GTS_OBJECT (gl->p));
  if (gl->s)
    gts_object_destroy (GTS_OBJECT (gl->s));
  g_slist_free (gl->solids);

  (* GTS_OBJECT_CLASS (gfs_gl_solid_class ())->parent_class->destroy) (o);
}

static void gl_solid_read (GtsObject ** o, GtsFile * fp)
{
  GfsGlSolid * gl = GFS_GL_SOLID (*o);
  GtsFileVariable var[] = {
    {GTS_INT, "reversed",   TRUE},
    {GTS_INT, "use_scalar", TRUE},
    {GTS_NONE}
  };

  (* GTS_OBJECT_CLASS (gfs_gl_solid_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  var[0].data = &gl->reversed;
  var[1].data = &gl->use_scalar;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;
}

static void gl_solid_write (GtsObject * o, FILE * fp)
{
  GfsGlSolid * gl = GFS_GL_SOLID (o);

  (* GTS_OBJECT_CLASS (gfs_gl_solid_class ())->parent_class->write) (o, fp);

  fprintf (fp, " {\n"
	   "  reversed = %d\n"
	   "  use_scalar = %d\n"
	   "}",
	   gl->reversed, (gl->use_scalar != NULL));
}

typedef struct {
  GtsPoint p[8];
  GfsSegment s[12];
  gdouble v[12];
  GtsVector n1[12];
} CellCube;

static void cell_size (FttCell * cell, FttVector * h)
{
  h->x = h->y = ftt_cell_size (cell);
#if !FTT_2D
  h->z = h->x;
#endif /* 3D */
}

static gboolean cube_intersections (FttCell * cell,
				    GfsGenericSurface * s,
				    GfsVariable * v,
				    FttVector p[12],
				    FttVector n[12],
				    gint orient[12],
				    gdouble val[12],
				    gint max_level)
{
  CellCube cube;
  FttVector o, h;
  guint i;
  gint inside[8] = {0,0,0,0,0,0,0,0};
  gboolean cut = FALSE;

  ftt_cell_pos (cell, &o);
  cell_size (cell, &h);
  for (i = 0; i < FTT_DIMENSION; i++)
    (&o.x)[i] -= (&h.x)[i]/2.;
  for (i = 0; i < 8; i++) { /* for each vertex of the cube */
    cube.p[i].x = o.x + h.x*vertex[i].x;
    cube.p[i].y = o.y + h.y*vertex[i].y;
    cube.p[i].z = o.z + h.z*vertex[i].z;
  }

  for (i = 0; i < 12; i++) {
    GfsSegment * e = &cube.s[i];
    e->E = &cube.p[edge1[i][0]];
    e->D = &cube.p[edge1[i][1]];
    gfs_surface_segment_intersection (s, cell, e);
    if (e->n % 2 != 0) {
      gfs_surface_segment_normal (s, cell, e, cube.n1[i]);
      cut = TRUE;
    }
  }
  
  if (!cut)
    return FALSE;

  if (v)
    for (i = 0; i < 8; i++) /* for each vertex of the cube */
      cube.v[i] = gfs_cell_corner_value (cell, corner[i], v, max_level);
  
  for (i = 0; i < 12; i++) { /* for each edge of the cube */
    GfsSegment * e = &cube.s[i];
    if (e->n % 2 != 0) { /* only for odd number of intersections */
      guint j = edge1[i][0], k = edge1[i][1];

      /* intersection vertex position is the average of all the n[i] intersections */
      e->x /= e->n;

      p[i].x = (1. - e->x)*cube.p[j].x + e->x*cube.p[k].x;
      p[i].y = (1. - e->x)*cube.p[j].y + e->x*cube.p[k].y;
      p[i].z = (1. - e->x)*cube.p[j].z + e->x*cube.p[k].z;

      n[i].x = cube.n1[i][0];
      n[i].y = cube.n1[i][1];
      n[i].z = cube.n1[i][2];

      if (v)
	val[i] = (1. - e->x)*cube.v[j] + e->x*cube.v[k];

      g_assert (inside[j] == 0 || inside[j] == e->inside);
      g_assert (inside[k] == 0 || inside[k] == - e->inside);
      inside[j] = e->inside;
      inside[k] = - e->inside;
      orient[i] = (inside[j] > 0);
    }
    else
      orient[i] = -1;
  }

  return TRUE;
}

#define solid_cube_vertices(cut,s,var,max_level,block) {\
  FttVector _p[12], _n[12], * v[12], * n[12];\
  gdouble _val[12], val[12];\
  gint _orient[12];\
  guint _i;\
  if (cube_intersections (cell, s, var, _p, _n, _orient, _val, max_level)) \
    for (_i = 0; _i < 12; _i++) {					\
      guint nv = 0, _e = _i;						\
      while (_orient[_e] >= 0) {					\
	guint _m = 0, * _ne = connect[_e][_orient[_e]];			\
	n[nv] = &(_n[_e]);						\
	val[nv] = _val[_e];						\
	v[nv++] = &(_p[_e]);						\
	_orient[_e] = -1;						\
	while (_m < 3 && _orient[_e] < 0)				\
	  _e = _ne[_m++];						\
      }									\
      if (nv > 2) {							\
	block								\
      }								        \
    }\
}

static void gl_solid (FttCell * cell, GfsGl * gl)
{
  GfsGlSolid * gls = GFS_GL_SOLID (gl);
  polygons * p = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gls->p));

  if (p) {
    polygons_draw (p, gls);
    gl->size++;
  }
}

static void gl_solid_update (FttCell * cell, GfsGenericSurface * s, GfsGl * gl)
{
  GfsGlSolid * gls = GFS_GL_SOLID (gl);
  polygons * p = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gls->p));

  g_assert (!p);
  solid_cube_vertices (cut, s, gls->use_scalar ? GFS_GL_SCALAR (gl)->v : NULL, gl->maxlevel, {
      guint i;
      polygon * q;

      p = polygons_add (p, nv);
      q = &p->p[p->n - 1];

      for (i = 0; i < nv; i++) {
	q->v[i] = *(v[i]);
	q->n[i] = *(n[i]);
	q->val[i] = val[i];
      }
    });
  if (p) {
    polygons_draw (p, gls);
    GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, gls->p)) = p;
    gl->size++;
  }
}

static gboolean is_solid (FttCell * cell, gpointer data)
{
  return GFS_IS_MIXED (cell);
}

static void gl_solid_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlSolid * gls = GFS_GL_SOLID (gl);

  glShadeModel (GL_SMOOTH);
  if (gls->use_scalar && gls->use_scalar != GFS_GL_SCALAR (gl)->v) {
    gls->use_scalar = GFS_GL_SCALAR (gl)->v;
    gfs_gl_solid_reset (gls);
  }
  if (gls->use_scalar && !gfs_gl_vector_format (gl)) {
    glEnable (GL_TEXTURE_1D);
    gfs_colormap_texture (GFS_GL_SCALAR (gl)->cmap);
    glColor3f (1., 1., 1.);
  }

  gl->size = 0;
  if (gls->needs_updating) {
    GSList * i = gls->solids;
    while (i) {
      gfs_domain_traverse_cut (GFS_DOMAIN (gl->sim), GFS_SOLID (i->data)->s, 
			       FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			       (FttCellTraverseCutFunc) gl_solid_update, gl);
      i = i->next;
    }
    gls->needs_updating = FALSE;
  }
  else {
    gdouble res = f->res;
    f->res = 0.;
    gfs_gl_cell_traverse_visible_condition (gl, f, is_solid, NULL,
					    (FttCellTraverseFunc) gl_solid, gl);
    f->res = res;
  }

  if (gls->use_scalar && !gfs_gl_vector_format (gl))
    glDisable (GL_TEXTURE_1D);
}

static void reset_p_s (FttCell * cell, GfsGlSolid * gl)
{
  GFS_VALUE (cell, gl->p) = 0.;
  GFS_VALUE (cell, gl->s) = 0.;
}

static void gl_solid_set_simulation (GfsGl * object, GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  GfsGlSolid * gls = GFS_GL_SOLID (object);

  gfs_gl_solid_reset (gls);

  (*GFS_GL_CLASS (GTS_OBJECT_CLASS (gfs_gl_solid_class ())->parent_class)->set_simulation)
    (object, sim);

  if (gls->p)
    gts_object_destroy (GTS_OBJECT (gls->p));
  gls->p = gfs_temporary_variable (domain);
  if (gls->s)
    gts_object_destroy (GTS_OBJECT (gls->s));
  gls->s = gfs_temporary_variable (domain);
  g_slist_free (gls->solids);
  gls->solids = gfs_simulation_get_solids (sim);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) reset_p_s, gls);
}

static gboolean gl_solid_relevant (GfsSimulation * sim)
{
  return (sim->solids->items != NULL);
}

static void gl_solid_class_init (GfsGlClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gl_solid_destroy;
  GTS_OBJECT_CLASS (klass)->read = gl_solid_read;
  GTS_OBJECT_CLASS (klass)->write = gl_solid_write;
  klass->set_simulation = gl_solid_set_simulation;
  klass->draw = gl_solid_draw;
  klass->pick = NULL;
  klass->relevant = gl_solid_relevant;
}

static void gl_solid_init (GfsGl * gl)
{
  GtsColor c = { 1., 1., 1. };

  gl->lc = c;
  GFS_GL_SOLID (gl)->needs_updating = TRUE;
}

GfsGlClass * gfs_gl_solid_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_solid_info = {
      "GfsGlSolid",
      sizeof (GfsGlSolid),
      sizeof (GfsGlScalarClass),
      (GtsObjectClassInitFunc) gl_solid_class_init,
      (GtsObjectInitFunc) gl_solid_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_scalar_class ()),
				  &gfs_gl_solid_info);
  }

  return klass;
}

/* GfsGlFractions: Object */

static void gl_fractions (FttCell * cell, GfsGl * gl)
{
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  FttCellFace f;
  gdouble h = ftt_cell_size (cell)/2.;

  f.cell = cell;
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) 
    if (s->s[f.d] > 0. && s->s[f.d] < 1.) {
      FttComponent c = f.d/2;
      static FttVector o[FTT_DIMENSION][4] = {
	{{0.,1.,-1.}, {0.,1.,1.}, {0.,-1.,1.}, {0.,-1.,-1.}},
	{{1.,0.,-1.}, {1.,0.,1.}, {-1.,0.,1.}, {-1.,0.,-1.}},
	{{1.,-1.,0.}, {1.,1.,0.}, {-1.,1.,0.}, {-1.,-1.,0.}}
      };
      static FttVector n[FTT_NEIGHBORS] = {
	{ 1., 0., 0. }, { -1., 0., 0.},
	{ 0., 1., 0. }, { 0., -1., 0.},
	{ 0., 0., 1. }, { 0., 0., -1.},
      };
      FttVector p;
      gdouble l = h*sqrt (s->s[f.d]);

      glNormal3d (n[f.d].x, n[f.d].y, n[f.d].z);
      ftt_face_pos (&f, &p);
      glVertex3d (p.x + l*o[c][0].x, p.y + l*o[c][0].y, p.z + l*o[c][0].z);
      glVertex3d (p.x + l*o[c][1].x, p.y + l*o[c][1].y, p.z + l*o[c][1].z);
      glVertex3d (p.x + l*o[c][2].x, p.y + l*o[c][2].y, p.z + l*o[c][2].z);
      glVertex3d (p.x + l*o[c][3].x, p.y + l*o[c][3].y, p.z + l*o[c][3].z);
    }
  gl->size++;
}

static void gl_fractions_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glShadeModel (GL_CONSTANT);
  glBegin (GL_QUADS);
  gfs_gl_cell_traverse_visible_mixed (gl, f, (FttCellTraverseFunc) gl_fractions, gl);
  glEnd ();
}

/* GfsGlBoundaries: Object */

static void gl_boundaries (FttCell * cell, GfsGl * gl)
{
  if (!GFS_IS_MIXED (cell)) {
    FttCellNeighbors n;
    gdouble h = ftt_cell_size (cell);
    FttVector p;
    guint i;

    ftt_cell_neighbors (cell, &n);
    ftt_cell_pos (cell, &p);
    p.x -= h/2.; p.y -= h/2.; p.z -= h/2.;
    for (i = 0; i < 12; i++) {
      static FttDirection edge2[12][2] = {
	{FTT_BOTTOM,FTT_BACK},{FTT_BOTTOM,FTT_FRONT},
	{FTT_TOP,FTT_FRONT},  {FTT_TOP,FTT_BACK},
	{FTT_LEFT,FTT_BACK},  {FTT_LEFT,FTT_FRONT},
	{FTT_RIGHT,FTT_FRONT},{FTT_RIGHT,FTT_BACK},
	{FTT_BOTTOM,FTT_LEFT},{FTT_BOTTOM,FTT_RIGHT},
	{FTT_TOP,FTT_RIGHT},   {FTT_TOP,FTT_LEFT}
      };
      guint j = edge2[i][0], k = edge2[i][1];

      if ((!n.c[j] || GFS_CELL_IS_BOUNDARY (n.c[j])) &&
	  (!n.c[k] || GFS_CELL_IS_BOUNDARY (n.c[k]))) {
	gl->size++;
	glVertex3d (p.x + edge[i][0].x*h, 
		    p.y + edge[i][0].y*h,
		    p.z + edge[i][0].z*h);
	glVertex3d (p.x + edge[i][1].x*h, 
		    p.y + edge[i][1].y*h,
		    p.z + edge[i][1].z*h);
      }
    }
  }
}

/* GfsGlLevels: Object */

static void gl_face (FttCell * cell, GfsGlLevels * gl)
{
  FttVector v[12];
  FttDirection d[12];
  guint nv = gfs_cut_cube_vertices (cell, GFS_GL (gl)->maxlevel,
				    &GFS_GL2D (gl)->p[0], &GFS_GL2D (gl)->n,
				    v, d,
				    NULL, NULL);
  gboolean drawn = FALSE;

  if (nv > 2) {
    guint i;
    for (i = 0; i < nv - 1; i++) {
      FttCell * n = ftt_cell_neighbor (cell, d[i]);
      if (n && (floor (GFS_VALUE (cell, gl->v)) != 
		floor (GFS_VALUE (n, gl->v)))) {
	glVertex3d (v[i].x, v[i].y, v[i].z);
	glVertex3d (v[i+1].x, v[i+1].y, v[i+1].z);
	drawn = TRUE;
      }
    }
  }
  if (drawn)
    GFS_GL (gl)->size++;
}

/* GfsGlVectors: Object */

static void gl_vector (FttCell * cell, GfsGl * gl)
{
  GfsGlVectors * gls = GFS_GL_VECTORS (gl);
  if (gls->use_scalar && !GFS_HAS_DATA (cell, GFS_GL_SCALAR (gl)->v))
    return;
  if (gfs_plane_cuts_cell (GFS_GL2D (gl)->p, cell)) {
    GfsGl2D * gl2D = GFS_GL2D (gl);
    FttComponent c;
    FttVector pos, f;
    gdouble a;
    
    gl->size++;
    if (gls->use_scalar) {
      GfsGlScalar * gla = GFS_GL_SCALAR (gl);
      GtsColor c = gfs_colormap_color (gla->cmap, gla->max > gla->min ?
				       (GFS_VALUE (cell, gla->v) - gla->min)/
				       (gla->max - gla->min) :
				       0.5);
      glColor3f (c.r, c.g, c.b);
    }
    gfs_cell_cm (cell, &pos);
    a = ((gl2D->p[0].x - pos.x)*gl2D->n.x + 
	 (gl2D->p[0].y - pos.y)*gl2D->n.y + 
	 (gl2D->p[0].z - pos.z)*gl2D->n.z);
    pos.x += a*gl2D->n.x; pos.y += a*gl2D->n.y; pos.z += a*gl2D->n.z;
    for (c = 0; c < FTT_DIMENSION; c++)
      (&f.x)[c] = GFS_VALUE (cell, gls->v[c]);
    f.x *= gls->scale;
    f.y *= gls->scale;
    f.z *= gls->scale;
    glVertex3d (pos.x + f.x - (f.x - f.y/2.)/5., pos.y + f.y - (f.x/2. + f.y)/5., pos.z + f.z);
    glVertex3d (pos.x + f.x, pos.y + f.y, pos.z + f.z);
    glVertex3d (pos.x + f.x, pos.y + f.y, pos.z + f.z);
    glVertex3d (pos.x + f.x - (f.x + f.y/2.)/5., pos.y + f.y + (f.x/2. - f.y)/5., pos.z + f.z);
    glVertex3d (pos.x, pos.y, pos.z);
    glVertex3d (pos.x + f.x, pos.y + f.y, pos.z + f.z);
  }
}

/* GfsGlEllipses: Object */

static void gl_ellipse (FttCell * cell, GfsGl * gl)
{
  GfsGlEllipses * gls = GFS_GL_ELLIPSES (gl);
  if (gls->use_scalar && !GFS_HAS_DATA (cell, GFS_GL_SCALAR (gl)->v))
    return;
  if (gfs_plane_cuts_cell (GFS_GL2D (gl)->p, cell)) {
    GfsGl2D * gl2D = GFS_GL2D (gl);
    FttVector pos;
    gdouble a, t;
    
    gl->size++;
    if (gls->use_scalar) {
      GfsGlScalar * gla = GFS_GL_SCALAR (gl);
      GtsColor c = gfs_colormap_color (gla->cmap, gla->max > gla->min ?
				       (GFS_VALUE (cell, gla->v) - gla->min)/
				       (gla->max - gla->min) :
				       0.5);
      glColor3f (c.r, c.g, c.b);
    }
    gfs_cell_cm (cell, &pos);
    a = ((gl2D->p[0].x - pos.x)*gl2D->n.x + 
	 (gl2D->p[0].y - pos.y)*gl2D->n.y + 
	 (gl2D->p[0].z - pos.z)*gl2D->n.z);
    pos.x += a*gl2D->n.x; pos.y += a*gl2D->n.y; pos.z += a*gl2D->n.z;
    glBegin (GL_LINE_LOOP);
    for (t = 0.; t < 2.*M_PI; t += 2.*M_PI/20.) {
      gdouble cost = cos (t), sint = sin (t);
      
      glVertex3d (pos.x + gls->scale*(GFS_VALUE (cell, gls->v[0])*cost + 
				      GFS_VALUE (cell, gls->v[2])*sint),
		  pos.y + gls->scale*(GFS_VALUE (cell, gls->v[1])*cost + 
				      GFS_VALUE (cell, gls->v[3])*sint),
		  pos.z);
    }
    glEnd ();
  }
}

/* GfsGlLocate: Object */

static void gl_locate (FttCell * cell, GfsGl * gl)
{
  FttVector p;
  gdouble size = ftt_cell_size (cell)/2.;
  
  ftt_cell_pos (cell, &p);
  glBegin (GL_LINE_LOOP);
  glVertex3d (p.x - size, p.y - size, p.z - size);
  glVertex3d (p.x + size, p.y - size, p.z - size);
  glVertex3d (p.x + size, p.y + size, p.z - size);
  glVertex3d (p.x - size, p.y + size, p.z - size);
  glEnd ();
  glBegin (GL_LINE_LOOP);
  glVertex3d (p.x - size, p.y - size, p.z + size);
  glVertex3d (p.x + size, p.y - size, p.z + size);
  glVertex3d (p.x + size, p.y + size, p.z + size);
  glVertex3d (p.x - size, p.y + size, p.z + size);
  glEnd ();
  glBegin (GL_LINES);
  glVertex3d (p.x - size, p.y - size, p.z - size);
  glVertex3d (p.x - size, p.y - size, p.z + size);
  glVertex3d (p.x + size, p.y - size, p.z - size);
  glVertex3d (p.x + size, p.y - size, p.z + size);
  glVertex3d (p.x + size, p.y + size, p.z - size);
  glVertex3d (p.x + size, p.y + size, p.z + size);
  glVertex3d (p.x - size, p.y + size, p.z - size);
  glVertex3d (p.x - size, p.y + size, p.z + size);
  glEnd ();  
  gl->size++;
}

/* GfsGlStreamline: Object */

static GSList * circle_profile (GtsVertexClass * klass, 
				gdouble radius, guint np)
{
  GSList * lp = NULL;
  guint i;

  for (i = 0; i <= np; i++) {
    gdouble a = 2.*M_PI*i/(gdouble) np;
    gdouble cosa = cos (a), sina = sin (a);
    GtsPoint * p = gts_point_new (GTS_POINT_CLASS (klass), radius*cosa, radius*sina, 0.);
    gdouble * n = GTS_VERTEX_NORMAL (p)->n;

    n[0] = cosa; n[1] = sina; n[2] = 0.;
    lp = g_slist_prepend (lp, p);
  }
  return lp;
}

static void matrix_transpose (GtsMatrix * m)
{
  guint i, j;

  for (i = 1; i < 3; i++)
    for (j = 0; j < i; j++) {
      gdouble t = m[i][j];

      m[i][j] = m[j][i];
      m[j][i] = t;
    }
}

static void base (GtsMatrix * b, GtsPoint * p1, GtsPoint * p2)
{
  GtsVector x, y;

  x[0] = b[0][0];
  x[1] = b[1][0];
  x[2] = b[2][0];
  gts_vector_init (b[2], p2, p1);
  gts_vector_normalize (b[2]);
  gts_vector_cross (y, b[2], x);
  if (gts_vector_norm (y) > 1e-2) {
    b[1][0] = y[0];
    b[1][1] = y[1];
    b[1][2] = y[2];
    gts_vector_normalize (b[1]);
  }
  gts_vector_cross (b[0], b[1], b[2]);
  gts_vector_normalize (b[0]);
  matrix_transpose (b);
}

static void point_list (GtsMatrix * b, 
			GtsMatrix * c,
			GtsPoint * o,
			GSList * profile,
			GtsVertexNormal * pl, 
			guint np)
{
  guint i;

#if 0
  gboolean colored = FALSE;
  if (GTS_IS_COLORED_VERTEX (o) && 
      gts_object_class_is_from_class (GTS_OBJECT_CLASS (s->vertex_class),
				      GTS_OBJECT_CLASS (gts_colored_vertex_class ())))
    colored = TRUE;
#endif
  if (FALSE/*GFS_IS_TWISTED_VERTEX (o)*/) {
#if 0
    gdouble t = GFS_TWISTED_VERTEX (o)->theta;
    gdouble sint = sin (t), cost = cos (t);
    GtsMatrix * r = gts_matrix_new (cost, -sint, 0., 0.,
				    sint,  cost, 0., 0.,
				    0.,      0., 1., 0.,
				    0.,      0., 0., 0.);
    
    c = gts_matrix_product (b, r);
    gts_matrix_destroy (r);
#endif
  }
  else
    gts_matrix_assign (c,
		       b[0][0], b[0][1], b[0][2], 0.,
		       b[1][0], b[1][1], b[1][2], 0.,
		       b[2][0], b[2][1], b[2][2], 0.,
		       0., 0., 0., 0.);

  for (i = 0; i < np; i++, profile = profile->next, pl++) {
    GtsPoint * p = profile->data;
    gdouble * n = GTS_VERTEX_NORMAL (p)->n;
    GtsPoint n1;
    
    *GTS_POINT (pl) = *p;
#if 0
    if (colored)
      GTS_COLORED_VERTEX (v)->c = GTS_COLORED_VERTEX (o)->c;
#endif
    gts_point_transform (GTS_POINT (pl), c);
    GTS_POINT (pl)->x += o->x;
    GTS_POINT (pl)->y += o->y;
    GTS_POINT (pl)->z += o->z;

    n1.x = n[0]; n1.y = n[1]; n1.z = n[2];
    gts_point_transform (&n1, c);
    n = pl->n;
    n[0] = n1.x; n[1] = n1.y; n[2] = n1.z;
  }
}

static void draw_point (GtsPoint * p)
{
  glVertex3d (p->x, p->y, p->z);
}

static void draw_normal (GtsVertexNormal * v)
{
  glNormal3d (v->n[0], v->n[1], v->n[2]);
}

static void draw_faces (GtsVertexNormal * v1, GtsVertexNormal * v2, guint np)
{
  guint i;

  glBegin (GL_QUAD_STRIP);
  for (i = 0; i < np; i++, v1++, v2++) {
    draw_normal (v1);
    draw_point (GTS_POINT (v1));
    draw_normal (v2);
    draw_point (GTS_POINT (v2));    
  }
  glEnd ();
}

static GList * next_far_enough (GList * p, gdouble size)
{
  GtsPoint * ps;
  GList * pf = NULL;

  if (p == NULL)
    return NULL;
  ps = p->data;
  p = p->next;
  size *= size;
  while (p && !pf) {
    if (gts_point_distance2 (ps, p->data) > size)
      pf = p;
    p = p->next;
  }
  return pf;
}

static void extrude_profile (GSList * profile, GList * path)
{
  GtsMatrix * r, * c;
  GtsPoint * p0, * p1, * p2;
  GtsVertexNormal * pl1, * pl2, * tmp;
  GtsBBox * bb;
  gdouble size;
  guint np;

  g_return_if_fail (profile != NULL);
  g_return_if_fail (path != NULL);

  bb = gts_bbox_points (gts_bbox_class (), profile);
  size = bb->x2 - bb->x1;
  if (bb->y2 - bb->y1 > size)
    size = bb->y2 - bb->y1;
  gts_object_destroy (GTS_OBJECT (bb));

  size /= 4.;

  p0 = path->data;
  path = next_far_enough (path, size);
  if (path == NULL)
    return;
  p1 = path->data;
  r = gts_matrix_identity (NULL);
  c = gts_matrix_identity (NULL);
  np = g_slist_length (profile);
  pl1 = g_malloc (sizeof (GtsVertexNormal)*np);
  pl2 = g_malloc (sizeof (GtsVertexNormal)*np);

  base (r, p0, p1);
  point_list (r, c, p0, profile, pl1, np);
  do {
    path = next_far_enough (path, size);
    p2 = path ? path->data : NULL;
    if (p2)
      base (r, p0, p2);
    else
      base (r, p0, p1);
    point_list (r, c, p1, profile, pl2, np);
    draw_faces (pl1, pl2, np);
    tmp = pl1;
    pl1 = pl2;
    pl2 = tmp;
    p0 = p1;
    p1 = p2;
  } while (p1);

  g_free (pl1);
  g_free (pl2);
  gts_matrix_destroy (c);
  gts_matrix_destroy (r);
}

static void gl_streamline_draw (GfsGlStreamline * s,
				GfsGlStreamlines * gls)
{
  if (gls->radius > 0.) {
    GSList * profile = circle_profile (gts_vertex_normal_class (), gls->radius, 10);

    glShadeModel (GL_SMOOTH);
    extrude_profile (profile, s->l);

    g_slist_foreach (profile, (GFunc) gts_object_destroy, NULL);
    g_slist_free (profile);
  }
  else {
    glBegin (GL_LINE_STRIP);
    g_list_foreach (s->l, (GFunc) draw_point, NULL);
    glEnd ();
  }
}

