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

/* GfsGlCells: Object */

static void gl_cell (FttCell * cell, GfsGl * gl)
{
  FttVector p;
  gdouble size = ftt_cell_size (cell)/2.;
  
  ftt_cell_pos (cell, &p);
  glBegin (GL_LINE_LOOP);
  glVertex2d (p.x - size, p.y - size);
  glVertex2d (p.x + size, p.y - size);
  glVertex2d (p.x + size, p.y + size);
  glVertex2d (p.x - size, p.y + size);
  glEnd ();
  gl->size++;
}

/* GfsGlSquares: Object */

static void gl_square (FttCell * cell, GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  if (!GFS_HAS_DATA (cell, gls->v))
    return;
  FttVector p;
  gdouble size = ftt_cell_size (cell)/1.999;
  GtsColor c;

  gl->size++;
  ftt_cell_pos (cell, &p);
  c = gfs_colormap_color (gls->cmap, gls->max > gls->min ?
			  (GFS_VALUE (cell, gls->v) - gls->min)/(gls->max - gls->min) :
			  0.5);
  glColor3f (c.r, c.g, c.b);
  glVertex2d (p.x - size, p.y - size);
  glVertex2d (p.x + size, p.y - size);
  glVertex2d (p.x + size, p.y + size);
  glVertex2d (p.x - size, p.y + size);
}

static void gl_squares_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glBegin (GL_QUADS);
  gfs_gl_normal (gl);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_square, gl);
  glEnd ();

  (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass->parent_class)->draw) (gl, f);
}

/* GfsGlSolid: Object */

#define CUT(d) (s->s[d] > 0. && s->s[d] < 1.)

static void gl_solid (FttCell * cell, GfsGl * gl)
{
  if (GFS_IS_MIXED (cell)) {
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    guint n = 0;
    FttDirection d;
    FttVector c, p[4];
    gdouble h = ftt_cell_size (cell)/1.999;
    static guint etod[] = { 3, 0, 2, 1 };
    
    gl->size++;
    ftt_cell_pos (cell, &c);  
    p[0].x = c.x - h; p[0].y = c.y - h;
    p[1].x = c.x + h; p[1].y = c.y - h;
    p[2].x = c.x + h; p[2].y = c.y + h;
    p[3].x = c.x - h; p[3].y = c.y + h;
    
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (s->s[d] > 0. && s->s[d] < 1.)
	n++;
    switch (n) {
    case 2: {
      guint k, m, n1 = 0;
      FttVector r[2];
      gboolean ins;
      static guint ctoc[4][2] = {{0,2},{2,1},{1,3},{3,0}};
      
      for (m = 0; m < 4 && !CUT (etod[m]); m++);
      ins = CUT (ctoc[m][0]) ? s->s[ctoc[m][1]] : (s->s[ctoc[m][0]] != 1.);
      for (k = m; k < m + 4; k++) {
	guint i = k % 4, i1 = (i + 1) % 4;
	if (CUT (etod[i])) {
	  gdouble a = ins ? s->s[etod[i]] : 1. - s->s[etod[i]];
	  r[n1].x = (1. - a)*p[i].x + a*p[i1].x;
	  r[n1++].y = (1. - a)*p[i].y + a*p[i1].y;
	  ins = !ins;
	  if (n1 == 2) {
	    glVertex2d (r[0].x, r[0].y);
	    glVertex2d (r[1].x, r[1].y);
	    n1 = 0;
	  }
	}
      }
      break;
    }
    case 4:
      break;
    default:
#if 0
      g_assert_not_reached ();
#else
      {
	FttVector p;
	ftt_cell_pos (cell, &p);
	g_warning ("(%g,%g): %d", p.x, p.y, n);
      }
#endif
    }
  }
}

static void gl_solid_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glBegin (GL_LINES);
  glNormal3d (0., 0., 1.);
  gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc) gl_solid, gl);
  glEnd ();
  glPopMatrix ();
}

static void gl_solid_init (GfsGl * gl)
{
  GtsColor c = { 0., 1., 0. };

  gl->lc = c;
}

static gboolean gl_solid_relevant (GfsSimulation * sim)
{
  return (sim->solids->items != NULL);
}

static void gl_solid_class_init (GfsGlClass * klass)
{
  klass->draw = gl_solid_draw;
  klass->relevant = gl_solid_relevant;
}

GfsGlClass * gfs_gl_solid_class (void)
{
  static GfsGlClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_gl_solid_info = {
      "GfsGlSolid",
      sizeof (GfsGl),
      sizeof (GfsGlClass),
      (GtsObjectClassInitFunc) gl_solid_class_init,
      (GtsObjectInitFunc) gl_solid_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_gl_class ()), &gfs_gl_solid_info);
  }

  return klass;
}

/* GfsGlFractions: Object */

static void gl_fractions (FttCell * cell, GfsGl * gl)
{
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  FttVector c;
  gdouble h = ftt_cell_size (cell)*0.45;
  
  gl->size++;
  ftt_cell_pos (cell, &c);
  glVertex2d (c.x + h, c.y - h*s->s[0]);
  glVertex2d (c.x + h, c.y + h*s->s[0]);
  glVertex2d (c.x - h, c.y - h*s->s[1]);
  glVertex2d (c.x - h, c.y + h*s->s[1]);
  glVertex2d (c.x - h*s->s[2], c.y + h);
  glVertex2d (c.x + h*s->s[2], c.y + h);
  glVertex2d (c.x - h*s->s[3], c.y - h);
  glVertex2d (c.x + h*s->s[3], c.y - h);
}

static void gl_fractions_draw (GfsGl * gl, GfsFrustum * f)
{
  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  glBegin (GL_LINES);
  gfs_gl_cell_traverse_visible_mixed (gl, f, (FttCellTraverseFunc) gl_fractions, gl);
  glEnd ();
  glPopMatrix ();
}

/* GfsGlBoundaries: Object */

static void gl_boundaries (FttCell * cell, GfsGl * gl)
{
  if (!GFS_IS_MIXED (cell)) {
    FttDirection d;
    FttCellNeighbors n;
    gdouble h = ftt_cell_size (cell)/2.;
    FttVector p;

    ftt_cell_neighbors (cell, &n);
    ftt_cell_pos (cell, &p);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (!n.c[d] || GFS_CELL_IS_BOUNDARY (n.c[d])) {
	static FttVector dp[FTT_NEIGHBORS][2] = {
	  {{1.,-1.,0.},{1.,1.,0.}},
	  {{-1.,1.,0.},{-1.,-1,0.}},
	  {{1.,1.,0.},{-1.,1.,0.}},
	  {{-1.,-1.,0.},{1.,-1.,0.}}
	};
	gl->size++;
	glVertex2d (p.x + dp[d][0].x*h, p.y + dp[d][0].y*h);
	glVertex2d (p.x + dp[d][1].x*h, p.y + dp[d][1].y*h);
      }
  }
}

/* GfsGlLevels: Object */

/**
 * ftt_face_gl:
 * @face: a #FttCellFace.
 *
 * Creates an OpenGL representation of @face.  
 */
static void ftt_face_gl (const FttCellFace * face)
{
  gdouble size;
  FttVector p;
#if FTT_2D
  static FttVector dp[FTT_NEIGHBORS][2] = {
    {{1.,-1.,0.},{1.,1.,0.}},
    {{-1.,1.,0.},{-1.,-1,0.}},
    {{1.,1.,0.},{-1.,1.,0.}},
    {{-1.,-1.,0.},{1.,-1.,0.}}
  };
#else  /* FTT_3D */
  static FttVector dp[FTT_NEIGHBORS][4] = {
    {{1.,-1.,1.},{1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.}},
    {{-1.,-1.,1.},{-1.,-1.,-1.},{-1.,1.,-1.},{-1.,1.,1.}},
    {{1.,1.,1.},{1.,1.,-1.},{-1,1.,-1.},{-1.,1.,1.}},
    {{1.,-1.,1.},{1.,-1.,-1.},{-1,-1.,-1.},{-1.,-1.,1.}},
    {{1.,-1.,1.},{1.,1.,1.},{-1.,1.,1.},{-1.,-1.,1.}},
    {{1.,-1.,-1.},{1.,1.,-1.},{-1.,1.,-1.},{-1.,-1.,-1.}},
  };
#endif /* FTT_3D */

  g_return_if_fail (face != NULL);

  size = ftt_cell_size (face->cell)/2.;
  ftt_cell_pos (face->cell, &p);
#if FTT_2D
  glVertex2d (p.x + dp[face->d][0].x*size, p.y + dp[face->d][0].y*size);
  glVertex2d (p.x + dp[face->d][1].x*size, p.y + dp[face->d][1].y*size);
#else  /* FTT_3D */
#if 0
  fprintf (fp, 
	   "OFF 4 1 4 "
	   "%g %g %g "
	   "%g %g %g "
	   "%g %g %g "
	   "%g %g %g "
	   "4 0 1 2 3\n",
	   p.x + dp[face->d][0].x*size,
	   p.y + dp[face->d][0].y*size,
	   p.z + dp[face->d][0].z*size,
	   p.x + dp[face->d][1].x*size,
	   p.y + dp[face->d][1].y*size,
	   p.z + dp[face->d][1].z*size,
	   p.x + dp[face->d][2].x*size,
	   p.y + dp[face->d][2].y*size,
	   p.z + dp[face->d][2].z*size,
	   p.x + dp[face->d][3].x*size,
	   p.y + dp[face->d][3].y*size,
	   p.z + dp[face->d][3].z*size);
#endif
#endif /* FTT_3D */
}

static void gl_face (FttCell * cell, GfsGlLevels * gl)
{
  FttCellNeighbors n;
  FttCellFace f;
  gboolean drawn = FALSE;

  ftt_cell_neighbors (cell, &n);
  f.cell = cell;
  for (f.d = 0; f.d < FTT_NEIGHBORS; f.d++) {
    f.neighbor = n.c[f.d];
    if (f.neighbor &&
	floor (GFS_VALUE (cell, gl->v)) != 
	floor (GFS_VALUE (f.neighbor, gl->v))) {
      ftt_face_gl (&f);
      drawn = TRUE;
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

  FttComponent c;
  FttVector pos, f;

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
  for (c = 0; c < FTT_DIMENSION; c++)
    (&f.x)[c] = GFS_VALUE (cell, gls->v[c]);
  f.x *= gls->scale;
  f.y *= gls->scale;
  glVertex2d (pos.x + f.x - (f.x - f.y/2.)/5., pos.y + f.y - (f.x/2. + f.y)/5.);
  glVertex2d (pos.x + f.x, pos.y + f.y);
  glVertex2d (pos.x + f.x, pos.y + f.y);
  glVertex2d (pos.x + f.x - (f.x + f.y/2.)/5., pos.y + f.y + (f.x/2. - f.y)/5.);
  glVertex2d (pos.x, pos.y);
  glVertex2d (pos.x + f.x, pos.y + f.y);
}

/* GfsGlEllipses: Object */

static void gl_ellipse (FttCell * cell, GfsGl * gl)
{
  GfsGlEllipses * gls = GFS_GL_ELLIPSES (gl);
  if (gls->use_scalar && !GFS_HAS_DATA (cell, GFS_GL_SCALAR (gl)->v))
    return;
  FttVector pos;
  gdouble t;

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
  glBegin (GL_LINE_LOOP);
  for (t = 0.; t < 2.*M_PI; t += 2.*M_PI/20.) {
    gdouble cost = cos (t), sint = sin (t);

    glVertex2d (pos.x + gls->scale*(GFS_VALUE (cell, gls->v[0])*cost + 
				    GFS_VALUE (cell, gls->v[2])*sint),
		pos.y + gls->scale*(GFS_VALUE (cell, gls->v[1])*cost + 
				    GFS_VALUE (cell, gls->v[3])*sint));
  }
  glEnd ();
}

/* GfsGlLocate: Object */

static void gl_locate (FttCell * cell, GfsGl * gl)
{
  gl_cell (cell, gl);
}

/* GfsGlLinear: Object */

#define param(v) (gls->max > gls->min ? ((v) - gls->min)/(gls->max - gls->min) : 0.5)
#define color(v) (gfs_colormap_color (gls->cmap, param (v)))

static void gl_linear_color_constant (FttCell * cell, GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  if (!GFS_HAS_DATA (cell, gls->v))
    return;
  FttVector p;
  gdouble v, size = ftt_cell_size (cell)/1.999;
  GtsColor c;
  FttDirection d[2];

  gl->size++;
  ftt_cell_pos (cell, &p);
  d[0] = FTT_LEFT; d[1] = FTT_BOTTOM;
  v = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  c = color (v);
  glColor3f (c.r, c.g, c.b);
  glVertex2d (p.x - size, p.y - size);
  d[0] = FTT_RIGHT; d[1] = FTT_BOTTOM;
  v = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  c = color (v);
  glColor3f (c.r, c.g, c.b);
  glVertex2d (p.x + size, p.y - size);
  d[0] = FTT_RIGHT; d[1] = FTT_TOP;
  v = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  c = color (v);
  glColor3f (c.r, c.g, c.b);
  glVertex2d (p.x + size, p.y + size);
  d[0] = FTT_LEFT; d[1] = FTT_TOP;
  v = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  c = color (v);
  glColor3f (c.r, c.g, c.b);
  glVertex2d (p.x - size, p.y + size);
}

static void gl_linear_texture_constant (FttCell * cell, GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  if (!GFS_HAS_DATA (cell, gls->v))
    return;
  FttVector p;
  gdouble v1, v2, v3, v4, p1, p2, p3, p4;
  gdouble size = ftt_cell_size (cell)/1.999;
  FttDirection d[2];

  gl->size++;
  ftt_cell_pos (cell, &p);

  d[0] = FTT_LEFT; d[1] = FTT_BOTTOM;
  v1 = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  p1 = param (v1);
  d[0] = FTT_RIGHT; d[1] = FTT_BOTTOM;
  v2 = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  p2 = param (v2);
  d[0] = FTT_RIGHT; d[1] = FTT_TOP;
  v3 = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  p3 = param (v3);
  d[0] = FTT_LEFT; d[1] = FTT_TOP;
  v4 = gfs_cell_corner_value (cell, d, gls->v, gl->maxlevel);
  p4 = param (v4);

  glBegin (GL_TRIANGLE_FAN);
  glTexCoord1d ((p1 + p2 + p3 + p4)/4.);
  glVertex2d (p.x, p.y);
  glTexCoord1d (p1);
  glVertex2d (p.x - size, p.y - size);
  glTexCoord1d (p2);
  glVertex2d (p.x + size, p.y - size);
  glTexCoord1d (p3);
  glVertex2d (p.x + size, p.y + size);
  glTexCoord1d (p4);
  glVertex2d (p.x - size, p.y + size);
  glTexCoord1d (p1);
  glVertex2d (p.x - size, p.y - size);
  glEnd ();
}

typedef struct {
  gdouble p, z;
  FttVector n;
} Vertex;

static void new_vertex (FttCell * cell, FttDirection * d, GfsGl * gl, Vertex * v)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  GfsGlLinear * gll = GFS_GL_LINEAR (gl);
  GfsInterpolator inter;
  gfs_cell_corner_interpolator (cell, d, gl->maxlevel, FALSE, &inter);
  gdouble val = 0.;
  guint i;
  v->z = 0.;
  v->n.x = v->n.y = 0.;
  v->n.z = 1.;
  for (i = 0; i < inter.n; i++) {
    val += inter.w[i]*GFS_VALUE (inter.c[i], gls->v);
    v->z += inter.w[i]*GFS_VALUE (inter.c[i], gll->vf->v);
    v->n.x -= inter.w[i]*GFS_VALUE (inter.c[i], gll->nx);
    v->n.y -= inter.w[i]*GFS_VALUE (inter.c[i], gll->ny);
  }
  v->p = val;
  if (gll->reversed) {
    v->n.x = - v->n.x;
    v->n.y = - v->n.y;
    v->n.z = - 1.;
  }
}

static void gl_linear_color (FttCell * cell, GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  if (!GFS_HAS_DATA (cell, gls->v))
    return;
  gboolean use_scalar = GFS_GL_LINEAR (gl)->use_scalar != NULL;
  FttVector p;
  Vertex v;
  gdouble size = ftt_cell_size (cell)/1.999;
  GtsColor c;
  FttDirection d[2];

  gl->size++;
  ftt_cell_pos (cell, &p);
  d[0] = FTT_LEFT; d[1] = FTT_BOTTOM;
  new_vertex (cell, d, gl, &v);
  if (use_scalar) {
    c = color (v.p);
    glColor3f (c.r, c.g, c.b);
  }
  glNormal3d (v.n.x, v.n.y, v.n.z);
  glVertex3d (p.x - size, p.y - size, v.z);
  d[0] = FTT_RIGHT; d[1] = FTT_BOTTOM;
  new_vertex (cell, d, gl, &v);
  if (use_scalar) {
    c = color (v.p);
    glColor3f (c.r, c.g, c.b);
  }
  glNormal3d (v.n.x, v.n.y, v.n.z);
  glVertex3d (p.x + size, p.y - size, v.z);
  d[0] = FTT_RIGHT; d[1] = FTT_TOP;
  new_vertex (cell, d, gl, &v);
  if (use_scalar) {
    c = color (v.p);
    glColor3f (c.r, c.g, c.b);
  }
  glNormal3d (v.n.x, v.n.y, v.n.z);
  glVertex3d (p.x + size, p.y + size, v.z);
  d[0] = FTT_LEFT; d[1] = FTT_TOP;
  new_vertex (cell, d, gl, &v);
  if (use_scalar) {
    c = color (v.p);
    glColor3f (c.r, c.g, c.b);
  }
  glNormal3d (v.n.x, v.n.y, v.n.z);
  glVertex3d (p.x - size, p.y + size, v.z);
}

static void gl_linear_texture (FttCell * cell, GfsGl * gl)
{
  GfsGlScalar * gls = GFS_GL_SCALAR (gl);
  if (!GFS_HAS_DATA (cell, gls->v))
    return;
  FttVector p;
  Vertex v[4];
  gboolean use_scalar = GFS_GL_LINEAR (gl)->use_scalar != NULL;
  gdouble size = ftt_cell_size (cell)/1.999;
  FttDirection d[2];

  gl->size++;
  ftt_cell_pos (cell, &p);

  d[0] = FTT_LEFT; d[1] = FTT_BOTTOM;
  new_vertex (cell, d, gl, &v[0]);
  d[0] = FTT_RIGHT; d[1] = FTT_BOTTOM;
  new_vertex (cell, d, gl, &v[1]);
  d[0] = FTT_RIGHT; d[1] = FTT_TOP;
  new_vertex (cell, d, gl, &v[2]);
  d[0] = FTT_LEFT; d[1] = FTT_TOP;
  new_vertex (cell, d, gl, &v[3]);

  glBegin (GL_TRIANGLE_FAN);
  if (use_scalar)
    glTexCoord1d (param ((v[0].p + v[1].p + v[2].p + v[3].p)/4.));
  glNormal3d ((v[0].n.x + v[1].n.x + v[2].n.x + v[3].n.x)/4.,
	      (v[0].n.y + v[1].n.y + v[2].n.y + v[3].n.y)/4.,
	      (v[0].n.z + v[1].n.z + v[2].n.z + v[3].n.z)/4.);
  glVertex3d (p.x, p.y, (v[0].z + v[1].z + v[2].z + v[3].z)/4.);
  if (use_scalar)
    glTexCoord1d (param (v[0].p));
  glNormal3d (v[0].n.x, v[0].n.y, v[0].n.z);
  glVertex3d (p.x - size, p.y - size, v[0].z);
  if (use_scalar)
    glTexCoord1d (param (v[1].p));
  glNormal3d (v[1].n.x, v[1].n.y, v[1].n.z);
  glVertex3d (p.x + size, p.y - size, v[1].z);
  if (use_scalar)
    glTexCoord1d (param (v[2].p));
  glNormal3d (v[2].n.x, v[2].n.y, v[2].n.z);
  glVertex3d (p.x + size, p.y + size, v[2].z);
  if (use_scalar)
    glTexCoord1d (param (v[3].p));
  glNormal3d (v[3].n.x, v[3].n.y, v[3].n.z);
  glVertex3d (p.x - size, p.y + size, v[3].z);
  if (use_scalar)
    glTexCoord1d (param (v[0].p));
  glNormal3d (v[0].n.x, v[0].n.y, v[0].n.z);
  glVertex3d (p.x - size, p.y - size, v[0].z);
  glEnd ();
}

static void gl_linear_draw (GfsGl * gl, GfsFrustum * f)
{
  gboolean constant = (gfs_function_get_constant_value (GFS_GL_LINEAR (gl)->vf->f) == 0.);
  gl->size = 0;
  glShadeModel (GL_SMOOTH);
  if (constant)
    gfs_gl_normal (gl);
  if (gfs_gl_vector_format (gl)) {
    glBegin (GL_QUADS);
    gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc)
					(constant ? 
					 gl_linear_color_constant : 
					 gl_linear_color), gl);
    glEnd ();
  }
  else {
    if (GFS_GL_LINEAR (gl)->use_scalar) {
      glEnable (GL_TEXTURE_1D);
      gfs_colormap_texture (GFS_GL_SCALAR (gl)->cmap);
      glColor3f (1., 1., 1.);
      gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc)
					  (constant ?
					   gl_linear_texture_constant :
					   gl_linear_texture), gl);
      glDisable (GL_TEXTURE_1D);
    }
    else
      gfs_gl_cell_traverse_visible_plane (gl, f, (FttCellTraverseFunc)
					  (constant ?
					   gl_linear_texture_constant :
					   gl_linear_texture), gl);
  }

  if (GFS_GL_LINEAR (gl)->use_scalar)
    (* GFS_GL_CLASS (GTS_OBJECT (gl)->klass->parent_class)->draw) (gl, f);
}

/* GfsGlIsoline: Object */

typedef struct {
  FttVector p;
  gdouble v;
} VertexIsoline;

static gdouble inter (VertexIsoline * v1, VertexIsoline * v2, gdouble v)
{
  if ((v1->v > v && v2->v <= v) || (v1->v <= v && v2->v > v)) {
    gdouble a = (v - v1->v)/(v2->v - v1->v);
    glVertex3d (v1->p.x + a*(v2->p.x - v1->p.x),
		v1->p.y + a*(v2->p.y - v1->p.y),
		v1->p.z + a*(v2->p.z - v1->p.z));
#if DEBUG
    fprintf (stderr, "  %g %g %g\n",
	     v1->p.x + a*(v2->p.x - v1->p.x),
	     v1->p.y + a*(v2->p.y - v1->p.y),
	     v1->p.z + a*(v2->p.z - v1->p.z));
#endif
    return a;
  }
  return -1;
}

static void vertex_isoline (FttCell * cell, FttDirection * d, GfsGl * gl, VertexIsoline * v)
{
  ftt_corner_pos (cell, d, &v->p);
  if ((v->p.z = gfs_function_get_constant_value (GFS_GL_LINEAR (gl)->vf->f)) != 0.) {
    GfsInterpolator inter;
    guint i;
    gfs_cell_corner_interpolator (cell, d, gl->maxlevel, FALSE, &inter);
    v->v = v->p.z = 0.;
    for (i = 0; i < inter.n; i++) {
      v->v += inter.w[i]*GFS_VALUE (inter.c[i], GFS_GL_SCALAR (gl)->v);
      v->p.z += inter.w[i]*GFS_VALUE (inter.c[i], GFS_GL_LINEAR (gl)->vf->v);
    }
  }
  else
    v->v = gfs_cell_corner_value (cell, d, GFS_GL_SCALAR (gl)->v, gl->maxlevel);
}

/* the directions corresponding to the corners */
static FttDirection neighbor[4] = { FTT_BOTTOM, FTT_RIGHT, FTT_TOP, FTT_LEFT };
/* the corners corresponding to the directions */
static FttDirection icorner[4] = { 1, 3, 2, 0 };

static void set_used (FttCell * cell, gint i, GfsGl * gl)
{
  int v = GFS_VALUE (cell, GFS_GL_ISOLINE (gl)->used);
  v |= (1 << (i + 1));
  GFS_VALUE (cell, GFS_GL_ISOLINE (gl)->used) = v;
}

static gboolean get_used (FttCell * cell, gint i, GfsGl * gl)
{
  return (((int) GFS_VALUE (cell, GFS_GL_ISOLINE (gl)->used)) & (1 << (i + 1))) != 0;
}

static gboolean get_visible (FttCell * cell, GfsGl * gl)
{
  return GFS_VALUE (cell, GFS_GL_ISOLINE (gl)->used) != 0.;
}

static void follow_isoline (FttCell * cell, gint start, GfsGl * gl, gdouble val)
{
  while (cell &&
	 !GFS_CELL_IS_BOUNDARY (cell) &&
	 get_visible (cell, gl) && 
	 GFS_HAS_DATA (cell, GFS_GL_SCALAR (gl)->v)) {
    set_used (cell, start, gl);
    
    gint i1 = (start + 1) % 4;
    VertexIsoline v1, v2;
    vertex_isoline (cell, corner[i1], gl, &v1);
    FttCell * next = NULL;
    while (i1 != start && !get_used (cell, i1, gl)) {
      gint i2 = (i1 + 1) % 4;
      vertex_isoline (cell, corner[i2], gl, &v2);
      gdouble a;
      if ((a = inter (&v1, &v2, val)) >= 0.) {
	gl->size++;
	set_used (cell, i1, gl);
	FttDirection d = neighbor[i1];
	FttCell * n = ftt_cell_neighbor (cell, d);
	if (n) {
	  if (!FTT_CELL_IS_LEAF (n) && !get_visible (n, gl)) {
	    FttCellChildren child;
	    ftt_cell_children_direction (n, FTT_OPPOSITE_DIRECTION (d), &child);
	    gboolean visible0 = (child.c[0] && get_visible (child.c[0], gl));
	    gboolean visible1 = (child.c[1] && get_visible (child.c[1], gl));
	    if (visible0) {
	      if (visible1) {
		if (d == FTT_LEFT || d == FTT_BOTTOM)
		  n = child.c[a > 0.5];
		else
		  n = child.c[a < 0.5];
	      }
	      else
		n = child.c[0];
	    }
	    else if (visible1)
	      n = child.c[1];
	    else /* invisible children */
	      n = NULL;
	  }
	}
	next = n;
	start = (i1 + 2) % 4;
	break;
      }
      i1 = i2;
      v1 = v2;
    }
    cell = next;
  }
}

static void gl_isoline (FttCell * cell, GfsGl * gl)
{
  if (!GFS_HAS_DATA (cell, GFS_GL_SCALAR (gl)->v))
    return;
  gdouble val = GFS_GL_ISOLINE (gl)->val;

  gint i1 = 0;
  VertexIsoline v1, v2;
  vertex_isoline (cell, corner[i1], gl, &v1);
  do {
    gint i2 = (i1 + 1) % 4;
    vertex_isoline (cell, corner[i2], gl, &v2);
    if (!get_used (cell, i1, gl) && v1.v > val && v2.v <= val) {
#if DEBUG
      fprintf (stderr, "# glBegin (GL_LINE_LOOP) %g\n", val);
      glColor3f (rand()/(float)RAND_MAX, rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);
#endif
#if 0 /* fixme: this should work but the openGL drivers seem to be buggy... */
      FttCell * n = ftt_cell_neighbor (cell, neighbor[i1]);
      g_assert (n);
      GFS_VALUE (n, GFS_GL_ISOLINE (gl)->used[(i1 + 2) % 4]) = TRUE;
      glBegin (GL_LINE_LOOP);
#else
      glBegin (GL_LINE_STRIP); /* ... and this seems to work OK */
#endif
      inter (&v1, &v2, val);
      follow_isoline (cell, i1, gl, val);
      glEnd ();
#if DEBUG
      fprintf (stderr, "# glEnd ()\n\n");
#endif
    }
    i1 = i2;
    v1 = v2;
  } while (i1 != 0);
}

static gboolean neighbor_is_boundary (FttCell * cell, GfsGl * gl, FttDirection d)
{
  if (!cell || GFS_CELL_IS_BOUNDARY (cell))
    return TRUE;
  gdouble r = ftt_cell_size (cell)*GFS_DIAGONAL;
  FttVector p;
  ftt_cell_pos (cell, &p);
  return (gfs_sphere_in_frustum (&p, r, GFS_GL_ISOLINE (gl)->f) == GTS_OUT);
}

static void gl_boundary_isoline (FttCell * cell, GfsGl * gl)
{
  if (!GFS_HAS_DATA (cell, GFS_GL_SCALAR (gl)->v))
    return;
  gdouble val = GFS_GL_ISOLINE (gl)->val;

  FttCellNeighbors n;
  ftt_cell_neighbors (cell, &n);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    gint i1 = icorner[d];
    if (!get_used (cell, i1, gl) && neighbor_is_boundary (n.c[d], gl, d)) {
      VertexIsoline v1, v2;
      vertex_isoline (cell, corner[i1], gl, &v1);
      gint i2 = (i1 + 1) % 4;
      vertex_isoline (cell, corner[i2], gl, &v2);
      if (v1.v > val && v2.v <= val) {
#if DEBUG
	fprintf (stderr, "# glBegin (GL_LINE_STRIP) %g\n", val);
	glColor3f (rand()/(float)RAND_MAX, rand()/(float)RAND_MAX, rand()/(float)RAND_MAX);
#endif
	glBegin (GL_LINE_STRIP);
	inter (&v1, &v2, val);
	follow_isoline (cell, i1, gl, val);
	glEnd ();
#if DEBUG
	fprintf (stderr, "# glEnd ()\n\n");
#endif
      }
    }
  }
}

static void set_visible (FttCell * cell, GfsVariable * v)
{
  GFS_VALUE (cell, v) = 1.;
}

static void gl_isoline_draw (GfsGl * gl, GfsFrustum * f)
{
  GfsGlIsoline * gli = GFS_GL_ISOLINE (gl);
  gl_isoline_update_levels (gl);

  gli->used = gfs_temporary_variable (GFS_DOMAIN (gl->sim));
  gfs_domain_cell_traverse (GFS_DOMAIN (gl->sim), FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, gli->used);
  gli->f = f;

  gl->size = 0;
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glTranslatef (0., 0., gl->p->lc);
  gfs_gl_normal (gl);
  int i;
  for (i = 0; i < gli->levels->len; i++) {
    gli->val = g_array_index (gli->levels, gdouble, i);
    gfs_gl_cell_traverse_visible_iso (gl, f, gli->min, gli->max, gli->val,
				      (FttCellTraverseFunc) set_visible, gli->used);
    gfs_gl_cell_traverse_visible_iso (gl, f, gli->min, gli->max, gli->val,
				      (FttCellTraverseFunc) gl_boundary_isoline, gl);
    gfs_gl_cell_traverse_visible_iso (gl, f, gli->min, gli->max, gli->val,
				      (FttCellTraverseFunc) gl_isoline, gl);
    gfs_gl_cell_traverse_visible_iso (gl, f, gli->min, gli->max, gli->val,
				      (FttCellTraverseFunc) gfs_cell_reset, gli->used);
  }
  glPopMatrix ();
				      
  gts_object_destroy (GTS_OBJECT (gli->used));
}

/* GfsGlVOF: Object */

static void gl_vof (FttCell * cell, GfsGl * gl)
{
  FttVector p[2], m;
  if (gfs_vof_facet (cell, GFS_VARIABLE_TRACER_VOF (GFS_GL_VOF (gl)->vf->v), p, &m) == 2) {
    if (GFS_GL_VOF (gl)->use_scalar) {
      GfsGlScalar * gls = GFS_GL_SCALAR (gl);
      GtsColor c = gfs_colormap_color (gls->cmap, gls->max > gls->min ?
				       (GFS_VALUE (cell, gls->v) - gls->min)/
				       (gls->max - gls->min) :
				       0.5);
      glColor3f (c.r, c.g, c.b);
    }
    glVertex2d (p[0].x, p[0].y);
    glVertex2d (p[1].x, p[1].y);
    gl->size++;
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
    glMatrixMode (GL_PROJECTION);
    glPushMatrix ();
    glTranslatef (0., 0., gl->p->lc);
    glBegin (GL_LINES);
    gfs_gl_cell_traverse_visible_condition (gl, f, is_vof, gl, (FttCellTraverseFunc) gl_vof, gl);
    glEnd ();
    glPopMatrix ();
  }
}

static void gl_vof_cut (GfsGl * gl, FttCell * cell, GfsGl2D * plane)
{
}

/* GfsGlStreamline: Object */

static void draw_point (GtsPoint * p)
{
  glVertex3d (p->x, p->y, p->z);
}

static void gl_streamline_draw (GfsGlStreamline * s,
				GfsGlStreamlines * gls)
{
  glBegin (GL_LINE_STRIP);
  g_list_foreach (s->l, (GFunc) draw_point, NULL);
  glEnd ();
}
