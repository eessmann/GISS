/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric Research
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

#include "map.h"

/* GfsMap: Object */

static void gfs_map_read (GtsObject ** o, GtsFile * fp)
{
  GfsMap * map = GFS_MAP (*o);
  GtsObjectClass * klass;
  gboolean class_changed = FALSE;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsMapClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_map_class ())) {
    gts_file_error (fp, "`%s' is not a GfsMap", fp->token->str);
    return;
  }
  if (klass != (*o)->klass) {
    *o = gts_object_new (klass);
    gts_object_destroy (GTS_OBJECT (map));
    map = GFS_MAP (*o);
    class_changed = TRUE;
  }
  gts_file_next_token (fp);
}

static void gfs_map_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, "%s", o->klass->info.name);
}

static void gfs_map_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_write;
}

static void identity (GfsMap * map, const FttVector * src, FttVector * dest)
{
  *dest = *src;
}

static void gfs_map_init (GfsMap * map)
{
  map->transform = map->inverse = identity;
}

GfsMapClass * gfs_map_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_info = {
      "GfsMap",
      sizeof (GfsMap),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_class_init,
      (GtsObjectInitFunc) gfs_map_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
				  &gfs_map_info);
  }

  return klass;
}

/* GfsMapFunction: Object */

static void gfs_map_function_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);
  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return;
    }
    else {
      static gchar name[3][3] = { "x:", "y:", "z:" };
      FttComponent c;
      for (c = 0; c < FTT_DIMENSION; c++)
	if (!strcmp (fp->token->str, name[c]))
	  break;
      if (c == FTT_DIMENSION) {
	gts_file_error (fp, "unknown keyword '%s'", fp->token->str);
	return;
      }
      gts_file_next_token (fp);

      GfsMapFunction * map = GFS_MAP_FUNCTION (*o);
      if (!map->inverse[c])
	map->inverse[c] = gfs_function_new (gfs_function_map_class (), 0.);
      gfs_function_read (map->inverse[c], gfs_object_simulation (*o), fp);
      if (fp->type == GTS_ERROR)
	return;

      if (fp->type != '\n') {
	gts_file_error (fp, "expecting new line");
	return;
      }
      while (fp->type == '\n')
	gts_file_next_token (fp);
	
      if (!map->transform[c])
	map->transform[c] = gfs_function_new (gfs_function_map_class (), 0.);
      gfs_function_read (map->transform[c], gfs_object_simulation (*o), fp);
      if (fp->type == GTS_ERROR)
	return;
    }
  }
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void gfs_map_function_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_map_function_class ())->parent_class->write) (o, fp);
  fputs (" {\n", fp);
  GfsMapFunction * map = GFS_MAP_FUNCTION (o);
  static gchar name[3][5] = { "  x:", "  y:", "  z:" };
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (map->transform[c]) {
      fputs (name[c], fp);
      gfs_function_write (map->inverse[c], fp);
      fputs ("\n    ", fp);
      gfs_function_write (map->transform[c], fp);
      fputc ('\n', fp);
    }
  fputc ('}', fp);
}

static void gfs_map_function_destroy (GtsObject * o)
{
  GfsMapFunction * map = GFS_MAP_FUNCTION (o);
  FttComponent c;

  for (c = 0; c < FTT_DIMENSION; c++) {
    if (map->transform[c])
      gts_object_destroy (GTS_OBJECT (map->transform[c]));
    if (map->inverse[c])
      gts_object_destroy (GTS_OBJECT (map->inverse[c]));
  }

  (* GTS_OBJECT_CLASS (gfs_map_function_class ())->parent_class->destroy) (o);
}

static void gfs_map_function_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_function_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_function_write;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_map_function_destroy;
}

static void map_function_transform (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsMapFunction * mf = GFS_MAP_FUNCTION (map);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (mf->transform[c])
      (&dest->x)[c] = gfs_function_spatial_value (mf->transform[c], src);
    else
      (&dest->x)[c] = (&src->x)[c];
}

static void map_function_inverse (GfsMap * map, const FttVector * src, FttVector * dest)
{
  GfsMapFunction * mf = GFS_MAP_FUNCTION (map);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (mf->inverse[c])
      (&dest->x)[c] = gfs_function_spatial_value (mf->inverse[c], src);
    else
      (&dest->x)[c] = (&src->x)[c];
}

static void gfs_map_function_init (GfsMap * map)
{
  map->transform = map_function_transform;
  map->inverse =   map_function_inverse;
}

GfsMapClass * gfs_map_function_class (void)
{
  static GfsMapClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_map_function_info = {
      "GfsMapFunction",
      sizeof (GfsMapFunction),
      sizeof (GfsMapClass),
      (GtsObjectClassInitFunc) gfs_map_function_class_init,
      (GtsObjectInitFunc) gfs_map_function_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_map_class ()), &gfs_map_function_info);
  }

  return klass;
}
