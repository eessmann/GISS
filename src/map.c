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
#include "variable.h"
#include "utils.h"

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

static void identity (GfsMap * map, const FttVector * src, FttVector * dest)
{
  *dest = *src;
}

static void gfs_map_class_init (GfsMapClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_map_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_map_write;
  GFS_MAP_CLASS (klass)->transform = GFS_MAP_CLASS (klass)->inverse = identity;
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
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
				  &gfs_map_info);
  }

  return klass;
}
