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
#include "refine.h"
#include "solid.h"

/* GfsRefine: Object */

static gboolean refine_maxlevel (FttCell * cell, GfsFunction * maxlevel)
{
  FttVector p;

  ftt_cell_pos (cell, &p);
  return (ftt_cell_level (cell) < gfs_function_value (maxlevel, &p, 0.));
}

static void refine_box (GfsBox * box, GfsFunction * maxlevel)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_maxlevel, maxlevel,
		   (FttCellInitFunc) gfs_cell_init, gfs_box_domain (box));
}

static void gfs_refine_refine (GfsRefine * refine, GfsSimulation * sim)
{
  gts_container_foreach (GTS_CONTAINER (sim),
			 (GtsFunc) refine_box, refine->maxlevel);
}

static void gfs_refine_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_REFINE (o)->maxlevel));
  (* GTS_OBJECT_CLASS (gfs_refine_class ())->parent_class->destroy) (o);
}

static void gfs_refine_write (GtsObject * object, FILE * fp)
{
  fprintf (fp, "%s", object->klass->info.name);
  gfs_function_write (GFS_REFINE (object)->maxlevel, fp);
}

static void gfs_refine_read (GtsObject ** o, GtsFile * fp)
{
  GfsRefine * refine = GFS_REFINE (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsRefineClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_refine_class ())) {
    gts_file_error (fp, "`%s' is not a GfsRefine", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (refine));
    refine = GFS_REFINE (*o);
    class_changed = TRUE;
  }
  gts_file_next_token (fp);

  gfs_function_read (refine->maxlevel, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (class_changed && fp->type != '\n' && klass->read)
    (* klass->read) (o, fp);
}

static void gfs_refine_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_refine;

  GTS_OBJECT_CLASS (klass)->destroy = gfs_refine_destroy;
  GTS_OBJECT_CLASS (klass)->write = gfs_refine_write;
  GTS_OBJECT_CLASS (klass)->read =  gfs_refine_read;
}

static void gfs_refine_init (GfsRefine * object)
{
  object->maxlevel = gfs_function_new (gfs_function_class (), 1.);
}

GfsRefineClass * gfs_refine_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_info = {
      "GfsRefine",
      sizeof (GfsRefine),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_class_init,
      (GtsObjectInitFunc) gfs_refine_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = 
      gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
			    &gfs_refine_info);
  }

  return klass;
}

GfsRefine * gfs_refine_new (GfsRefineClass * klass)
{
  GfsRefine * object;

  object = GFS_REFINE (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/* GfsRefineSolid: Object */

static void refine_cut_cell (FttCell * cell, GtsSurface * s, gpointer * data)
{
  GfsRefine * refine = data[0];
  GfsDomain * domain = data[1];
  FttVector p;

  ftt_cell_pos (cell, &p);
  if (ftt_cell_level (cell) < gfs_function_value (refine->maxlevel, &p, 0.))
    ftt_cell_refine_single (cell, (FttCellInitFunc) gfs_cell_init, domain);
}

static void gfs_refine_solid_refine (GfsRefine * refine, GfsSimulation * sim)
{
  if (sim->surface) {
    gpointer data[2];

    data[0] = refine;
    data[1] = sim;
    gfs_domain_traverse_cut (GFS_DOMAIN (sim), sim->surface, 
			     FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseCutFunc) refine_cut_cell, data);
  }
}

static void gfs_refine_solid_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_solid_refine;
}

GfsRefineClass * gfs_refine_solid_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_solid_info = {
      "GfsRefineSolid",
      sizeof (GfsRefine),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_solid_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_solid_info);
  }

  return klass;
}

/* GfsRefineSurface: Object */

static void refine_surface_destroy (GtsObject * object)
{
  GfsRefineSurface * d = GFS_REFINE_SURFACE (object);

  if (d->surface)
    gts_object_destroy (GTS_OBJECT (d->surface));

  (* GTS_OBJECT_CLASS (gfs_refine_surface_class ())->parent_class->destroy) (object);
}

static void refine_surface_write (GtsObject * o, FILE * fp)
{
  GfsRefineSurface * d = GFS_REFINE_SURFACE (o);
  
  (* GTS_OBJECT_CLASS (gfs_refine_surface_class ())->parent_class->write) (o, fp);
  fprintf (fp, " { ");
  gts_surface_write (d->surface, fp);
  fputs ("}\n", fp);
}

static void refine_surface_read (GtsObject ** o, GtsFile * fp)
{
  GfsRefineSurface * refine;

  (* GTS_OBJECT_CLASS (gfs_refine_surface_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  refine = GFS_REFINE_SURFACE (*o);
  if (fp->type != '{') {
    FILE * f;
    GtsFile * gf;

    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a string (filename)\n");
      return;
    }
    f = fopen (fp->token->str, "rt");
    if (f == NULL) {
      gts_file_error (fp, "cannot open file `%s'\n", fp->token->str);
      return;
    }
    gf = gts_file_new (f);
    if (gts_surface_read (refine->surface, gf)) {
      gts_file_error (fp, 
		      "file `%s' is not a valid GTS file\n"
		      "%s:%d:%d: %s",
		      fp->token->str, fp->token->str,
		      gf->line, gf->pos, gf->error);
      gts_file_destroy (gf);
      fclose (f);
      return;
    }
    gts_file_destroy (gf);
    fclose (f);
  }
  else { /* embedded GTS file */
    fp->scope_max++;
    gts_file_next_token (fp);
    if (gts_surface_read (refine->surface, fp))
      return;
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
  }

  gts_file_next_token (fp);
}

static void gfs_refine_surface_refine (GfsRefine * refine, GfsSimulation * sim)
{
  gpointer data[2];
  
  data[0] = refine;
  data[1] = sim;
  gfs_domain_traverse_cut (GFS_DOMAIN (sim), GFS_REFINE_SURFACE (refine)->surface,
			   FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			   (FttCellTraverseCutFunc) refine_cut_cell, data);
}

static void gfs_refine_surface_class_init (GfsRefineClass * klass)
{
  klass->refine = gfs_refine_surface_refine;

  GTS_OBJECT_CLASS (klass)->destroy = refine_surface_destroy;
  GTS_OBJECT_CLASS (klass)->write = refine_surface_write;
  GTS_OBJECT_CLASS (klass)->read = refine_surface_read;
}

static void refine_surface_init (GfsRefineSurface * r)
{
  r->surface = gts_surface_new (gts_surface_class (), 
				gts_face_class (), 
				gts_edge_class (), 
				gts_vertex_class ());
}

GfsRefineClass * gfs_refine_surface_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_surface_info = {
      "GfsRefineSurface",
      sizeof (GfsRefineSurface),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_surface_class_init,
      (GtsObjectInitFunc) refine_surface_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_surface_info);
  }

  return klass;
}

/* GfsRefineDistance: Object */

static void refine_distance_destroy (GtsObject * object)
{
  GfsRefineDistance * d = GFS_REFINE_DISTANCE (object);

  if (d->stree)
    gts_bb_tree_destroy (d->stree, TRUE);  

  (* GTS_OBJECT_CLASS (gfs_refine_distance_class ())->parent_class->destroy) (object);
}

static void refine_distance_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_refine_distance_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GFS_REFINE_DISTANCE (*o)->stree = gts_bb_tree_surface (GFS_REFINE_SURFACE (*o)->surface);
}

static gboolean refine_distance_maxlevel (FttCell * cell, GfsRefine * refine)
{
  FttVector pos;
  GtsPoint p;
  gdouble d;

  ftt_cell_pos (cell, &pos);
  p.x = pos.x; p.y = pos.y; p.z = pos.z;
  d = gts_bb_tree_point_distance (GFS_REFINE_DISTANCE (refine)->stree, &p,
				  (GtsBBoxDistFunc) gts_point_triangle_distance, NULL);
  return (ftt_cell_level (cell) < gfs_function_value (refine->maxlevel, &pos, d));
}

static void refine_distance (GfsBox * box, gpointer data)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_distance_maxlevel, data,
		   (FttCellInitFunc) gfs_cell_init, gfs_box_domain (box));
}

static void refine_distance_refine (GfsRefine * refine, GfsSimulation * sim)
{
  gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) refine_distance, refine);
}

static void gfs_refine_distance_class_init (GfsRefineClass * klass)
{
  klass->refine = refine_distance_refine;
  
  GTS_OBJECT_CLASS (klass)->destroy = refine_distance_destroy;
  GTS_OBJECT_CLASS (klass)->read = refine_distance_read;
}

GfsRefineClass * gfs_refine_distance_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_distance_info = {
      "GfsRefineDistance",
      sizeof (GfsRefineDistance),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_distance_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_surface_class ()),
				  &gfs_refine_distance_info);
  }

  return klass;
}

/* GfsRefineHeight: Object */

static void refine_height_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_refine_distance_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  
  if (gts_surface_volume (GFS_REFINE_SURFACE (*o)->surface) < 0.)
    gts_surface_foreach_face (GFS_REFINE_SURFACE (*o)->surface, 
			      (GtsFunc) gts_triangle_revert, NULL);
}

static gboolean refine_height_maxlevel (FttCell * cell, gpointer * data)
{
  GfsRefine * refine = data[0];
  GtsFace * f, ** guess = data[1];
  FttVector pos;
  GtsPoint p;

  ftt_cell_pos (cell, &pos);
  p.x = pos.x; p.y = pos.y;
  if ((f = gts_point_locate (&p, GFS_REFINE_SURFACE (refine)->surface, NULL))) {
    *guess = f;
    gts_triangle_interpolate_height (GTS_TRIANGLE (f), &p);
    return (ftt_cell_level (cell) < gfs_function_value (refine->maxlevel, &pos, p.z));
  }
  return FALSE;
}

static void refine_height (GfsBox * box, gpointer data)
{
  gpointer datum[2];
  GtsFace * guess = NULL;

  datum[0] = data;
  datum[1] = &guess;
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_height_maxlevel, datum,
		   (FttCellInitFunc) gfs_cell_init, gfs_box_domain (box));
}

static void refine_height_refine (GfsRefine * refine, GfsSimulation * sim)
{
  gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) refine_height, refine);
}

static void gfs_refine_height_class_init (GfsRefineClass * klass)
{
  klass->refine = refine_height_refine;
  GTS_OBJECT_CLASS (klass)->read = refine_height_read;
}

GfsRefineClass * gfs_refine_height_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_height_info = {
      "GfsRefineHeight",
      sizeof (GfsRefineSurface),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_height_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_surface_class ()),
				  &gfs_refine_height_info);
  }

  return klass;
}
