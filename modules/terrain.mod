/* Gerris - The GNU Flow Solver                       (-*-C-*-)
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
#include <stdlib.h>
#include <postgresql/libpq-fe.h>
#include "refine.h"

#define DBMAX "1000"

#ifdef G_LITTLE_ENDIAN
# define DBYTE "NDR"
#else
# define DBYTE "XDR"
#endif
#define DOUBLE_OID 701
#define INT32_OID  23

static gboolean inverse (gdouble m[3][3], gdouble mi[3][3])
{
  gdouble det;
  
  det = (m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) -
	 m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) +
	 m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]));
  if (det == 0.)
    return FALSE;
  
  mi[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1])/det; 
  mi[0][1] = (m[2][1]*m[0][2] - m[0][1]*m[2][2])/det;
  mi[0][2] = (m[0][1]*m[1][2] - m[1][1]*m[0][2])/det; 
  mi[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2])/det; 
  mi[1][1] = (m[0][0]*m[2][2] - m[2][0]*m[0][2])/det; 
  mi[1][2] = (m[1][0]*m[0][2] - m[0][0]*m[1][2])/det; 
  mi[2][0] = (m[1][0]*m[2][1] - m[2][0]*m[1][1])/det; 
  mi[2][1] = (m[2][0]*m[0][1] - m[0][0]*m[2][1])/det; 
  mi[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0])/det; 

  return TRUE;
}

/* GfsRefineTerrain: Header */

typedef struct _GfsRefineTerrain         GfsRefineTerrain;

struct _GfsRefineTerrain {
  GfsRefine parent;
  PGconn * conn;
  int wkb_oid;

  gchar * name;
  gchar * host, * dbname, * user, * passwd, * tables;
  GfsVariable * h, * hx, * hy, * he, * hn;
  GfsFunction * criterion;
};

#define GFS_REFINE_TERRAIN(obj)            GTS_OBJECT_CAST (obj,\
					           GfsRefineTerrain,\
					           gfs_refine_terrain_class ())
#define GFS_IS_REFINE_TERRAIN(obj)         (gts_object_is_from_class (obj,\
						   gfs_refine_terrain_class ()))
     
GfsRefineClass * gfs_refine_terrain_class  (void);

/* GfsRefineTerrain: Object */

static void update_terrain (FttCell * cell, GfsRefineTerrain * t)
{
  if (GFS_VALUE (cell, t->h) == G_MAXDOUBLE) {
    FttVector p, q[4];
    ftt_cell_pos (cell, &p);
    gdouble h = ftt_cell_size (cell)/2.;
    GfsSimulation * sim = gfs_object_simulation (t);

    /* build the cell boundaries in the non-projected coordinates system */
    q[0].x = p.x + h; q[0].y = p.y + h;
    gfs_simulation_map_inverse (sim, &q[0]);
    q[1].x = p.x - h; q[1].y = p.y + h;
    gfs_simulation_map_inverse (sim, &q[1]);
    q[2].x = p.x - h; q[2].y = p.y - h;
    gfs_simulation_map_inverse (sim, &q[2]);
    q[3].x = p.x + h; q[3].y = p.y - h;
    gfs_simulation_map_inverse (sim, &q[3]);

    /* get all the terrain data within this polygon */
    gchar * polygon = 
      g_strdup_printf ("SRID=4326;POLYGON((%.8g %.8g,%.8g %.8g,%.8g %.8g,%.8g %.8g,%.8g %.8g))",
		       q[0].x, q[0].y,
		       q[1].x, q[1].y,
		       q[2].x, q[2].y,
		       q[3].x, q[3].y,
		       q[0].x, q[0].y);
    //    fprintf (stderr, "%s\n", polygon);
    PGresult * res = PQexecPrepared (t->conn, "get_points_in_polygon", 1, 
				     (const char *const*) &polygon, NULL, NULL, 1);
    g_free (polygon);
    if (!res || PQresultStatus (res) != PGRES_TUPLES_OK) {
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_CRITICAL,
	     "PostgreSQL command failed\n%s", PQerrorMessage (t->conn));
      PQclear (res);
      return;
    }

    /* least-squares fit */
    guint i, n = PQntuples (res);
    gdouble H[4] = {0.,0.,0.,0.}, m[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}, mi[3][3];
    m[0][0] = n;
    for (i = 0; i < n; i++) {
      gint32 height = g_ntohl (*((gint32 *) PQgetvalue (res, i, 0)));
      gchar * wkb = PQgetvalue (res, i, 1);
      FttVector r;
      memcpy (&r.x, &wkb[5], 8);
      memcpy (&r.y, &wkb[5+8], 8);
      gfs_simulation_map (sim, &r);
      r.x -= p.x; r.y -= p.y;
      m[0][1] += r.x; m[0][2] += r.y;
      m[1][1] += r.x*r.x; m[1][2] += r.x*r.y;
      m[2][2] += r.y*r.y;
      H[0] += height;
      H[1] += r.x*height;
      H[2] += r.y*height;
      H[3] += height*height;
    }
    PQclear (res);

    m[1][0] = m[0][1]; m[2][0] = m[0][2]; m[2][1] = m[1][2];
    if (n > 2 && inverse (m, mi)) {
      gdouble h0, hx, hy;
      h0 = GFS_VALUE (cell, t->h) = mi[0][0]*H[0] + mi[0][1]*H[1] + mi[0][2]*H[2];
      hx = GFS_VALUE (cell, t->hx) = mi[1][0]*H[0] + mi[1][1]*H[1] + mi[1][2]*H[2];
      hy = GFS_VALUE (cell, t->hy) = mi[2][0]*H[0] + mi[2][1]*H[1] + mi[2][2]*H[2];
      GFS_VALUE (cell, t->he) = 
	sqrt (fabs (h0*(h0*m[0][0] + 2.*(hx*m[0][1] + hy*m[0][2] - H[0])) +
		    hx*(hx*m[1][1] + 2.*(hy*m[1][2] - H[1])) +
		    hy*(hy*m[2][2] - 2.*H[2]) +
		    H[3]));
    }
    else { /* not enough points for fit, use parent info */
      FttCell * parent = ftt_cell_parent (cell);
      if (parent) {
	gdouble h0 = GFS_VALUE (parent, t->h);
	gdouble hx = GFS_VALUE (parent, t->hx);
	gdouble hy = GFS_VALUE (parent, t->hy);
	FttVector pp;
	ftt_cell_pos (parent, &pp);
	GFS_VALUE (cell, t->h) = h0 + (p.x - pp.x)*hx + (p.y - pp.y)*hy;
	GFS_VALUE (cell, t->hx) = hx;
	GFS_VALUE (cell, t->hy) = hy;
	GFS_VALUE (cell, t->he) = GFS_VALUE (parent, t->he);
      }
      else {
	g_warning ("cannot initialise terrain for cell at level %d", ftt_cell_level (cell));
	GFS_VALUE (cell, t->h) = G_MAXDOUBLE;
	GFS_VALUE (cell, t->hx) = G_MAXDOUBLE;
	GFS_VALUE (cell, t->hy) = G_MAXDOUBLE;
	GFS_VALUE (cell, t->he) = G_MAXDOUBLE;
      }
    }
    GFS_VALUE (cell, t->hn) = n;
  }
}

static gboolean refine_terrain (FttCell * cell, GfsRefine * refine)
{
  update_terrain (cell, GFS_REFINE_TERRAIN (refine));
  if (ftt_cell_level (cell) >= gfs_function_value (refine->maxlevel, cell))
    return FALSE;
  return gfs_function_value (GFS_REFINE_TERRAIN (refine)->criterion, cell);
}

static void refine_box (GfsBox * box, GfsRefine * refine)
{
  ftt_cell_refine (box->root, 
		   (FttCellRefineFunc) refine_terrain, refine,
		   (FttCellInitFunc) gfs_cell_fine_init, gfs_box_domain (box));
}

static void set_undefined (FttCell * cell, GfsVariable * h)
{
  GFS_VALUE (cell, h) = G_MAXDOUBLE;
}

static void set_undefined_coarse_fine (FttCell * parent, GfsVariable * h)
{
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n])
      GFS_VALUE (child.c[n], h) = G_MAXDOUBLE;
}

static void terrain_refine (GfsRefine * refine, GfsSimulation * sim)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (refine);
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) set_undefined, t->h);
  t->h->coarse_fine = set_undefined_coarse_fine;
  gts_container_foreach (GTS_CONTAINER (sim), (GtsFunc) refine_box, refine);
  t->h->coarse_fine = gfs_cell_coarse_fine;
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) update_terrain, refine);
}

static void refine_terrain_destroy (GtsObject * object)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (object);
  g_free (t->name);
  if (t->conn)
    PQfinish (t->conn);
  gts_object_destroy (GTS_OBJECT (t->criterion));
  g_free (t->host);
  g_free (t->dbname);
  g_free (t->user);
  g_free (t->passwd);
  g_free (t->tables);
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->destroy) (object);
}

static PGresult * query (PGconn * conn, const gchar * s, int status, GtsFile * fp)
{
  PGresult * res = PQexec (conn, s);
  if (!res || PQresultStatus (res) != status) {
    gts_file_error (fp, "PostgreSQL command failed\n%s", PQerrorMessage (conn));
    PQclear (res);
    PQfinish (conn);
    return NULL;
  }
  if (status != PGRES_TUPLES_OK)
    PQclear (res);
  return res;
}

static gboolean postgres_connection (GfsRefineTerrain * t, GtsFile * fp)
{
  gchar * q = g_strjoin (" ", 
			 "host =", t->host,
			 "dbname =", t->dbname,
			 "user =", t->user,
			 "password =", t->passwd,
			 NULL);
  PGconn * conn = PQconnectdb (q);
  g_free (q);
  if (PQstatus (conn) == CONNECTION_BAD) {
    gts_file_error (fp, "could not connect to database\n%s", PQerrorMessage (conn));
    PQfinish (conn);
    return FALSE;
  }
  PGresult * res = query (conn, "SELECT OID FROM pg_type WHERE typname = 'bytea'",
			  PGRES_TUPLES_OK, fp);
  if (!res)
    return FALSE;
  if (PQntuples (res) != 1) {
    gts_file_error (fp, "query to find OID of geometry did not return 1 row");
    PQclear (res);
    PQfinish (conn);
    return FALSE;
  }
  t->wkb_oid = atoi (PQgetvalue (res, 0, 0));
  PQclear (res);

  q = g_strjoin (NULL,
		 "BEGIN; DECLARE mycursor BINARY CURSOR FOR SELECT height, ASBINARY(geom,'"
		 DBYTE "') FROM ",
		 t->tables,
		 " LIMIT 1",
		 NULL);
  if (!query (conn, q, PGRES_COMMAND_OK, fp)) {
    g_free (q);
    return FALSE;  
  }
  g_free (q);
  if (!(res = query (conn, "FETCH ALL IN mycursor", PGRES_TUPLES_OK, fp)))
    return FALSE;
  if (PQftype (res, 0) != INT32_OID) {
    gts_file_error (fp, "invalid type (%d) for field 'height'", PQftype (res, 0));
    PQclear (res);
    PQfinish (conn);
    return FALSE;
  }
  if (PQftype (res, 1) != t->wkb_oid) {
    gts_file_error (fp, "invalid type (%d) for field 'geom'", PQftype (res, 1));
    PQclear (res);
    PQfinish (conn);
    return FALSE;
  }
  if (PQntuples (res) == 0) {
    gts_file_error (fp, "database seems to be empty");
    PQclear (res);
    PQfinish (conn);
    return FALSE;
  }
  gchar * wkb = PQgetvalue (res, 0, 1);
  guint32 type;
  memcpy (&type, wkb + 1, 4);
  PQclear (res);
#ifdef G_LITTLE_ENDIAN
  if (wkb[0] != 1) {
#else
  if (wkb[0] != 0) {
#endif
    gts_file_error (fp, "endianness do not match");
    PQfinish (conn);
    return FALSE;
  }
  if (type != 1) {
    gts_file_error (fp, "geometry type (%d) is not POINT", type);
    PQfinish (conn);
    return FALSE;
  }
  if ((type & 0x80000000) != 0) {
    gts_file_error (fp, "geometry type should be 2D");
    PQfinish (conn);
    return FALSE;
  }
  if (!query (conn, "CLOSE mycursor; COMMIT", PGRES_COMMAND_OK, fp))
    return FALSE;
  q = g_strjoin (NULL,
		 "PREPARE get_points_in_polygon (text) AS"
		 " SELECT height, ASBINARY(geom,'" DBYTE "') FROM ", t->tables,
		 " WHERE geom && $1 AND ST_Contains($1,geom)"
		 //		 " ORDER BY randid LIMIT " DBMAX,
		 " LIMIT " DBMAX,
		 NULL);
  if (!query (conn, q, PGRES_COMMAND_OK, fp)) {
    g_free (q);
    return FALSE;
  }
  g_free (q);
  t->conn = conn;
  return TRUE;
}

static void refine_terrain_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (name)");
    return;
  }
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (*o);
  t->name = g_strdup (fp->token->str);
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_STRING, "host",     TRUE},
      {GTS_STRING, "dbname",   TRUE},
      {GTS_STRING, "user",     TRUE},
      {GTS_STRING, "password", TRUE},
      {GTS_STRING, "tables",   TRUE},
      {GTS_NONE}
    };
    gchar * host = NULL, * dbname = NULL, * user = NULL, * passwd = NULL, * tables = NULL;
    var[0].data = &host;
    var[1].data = &dbname;
    var[2].data = &user;
    var[3].data = &passwd;
    var[4].data = &tables;
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    if (var[0].set) { g_free (t->host);   t->host = host; }
    if (var[1].set) { g_free (t->dbname); t->dbname = dbname; }
    if (var[2].set) { g_free (t->user);   t->user = user; }
    if (var[3].set) { g_free (t->passwd); t->passwd = passwd; }
    if (var[4].set) { g_free (t->tables); t->tables = tables; }
  }

  if (!postgres_connection (t, fp))
    return;

  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  gchar * name = g_strjoin (NULL, t->name, "0", NULL);
  t->h = gfs_domain_get_or_add_variable (domain, name, "Terrain height");
  g_free (name);
  name = g_strjoin (NULL, t->name, "x", NULL);
  t->hx = gfs_domain_get_or_add_variable (domain, name, "Terrain slope");
  g_free (name);
  name = g_strjoin (NULL, t->name, "y", NULL);
  t->hy = gfs_domain_get_or_add_variable (domain, name, "Terrain slope");
  g_free (name);
  name = g_strjoin (NULL, t->name, "e", NULL);
  t->he = gfs_domain_get_or_add_variable (domain, name, "Terrain RMS error");
  g_free (name);
  name = g_strjoin (NULL, t->name, "n", NULL);
  t->hn = gfs_domain_get_or_add_variable (domain, name, "Terrain samples #");
  g_free (name);

  gfs_function_read (t->criterion, domain, fp);
  if (fp->type == GTS_ERROR) {
    PQfinish (t->conn);
    t->conn = NULL;
  }    
}

static void refine_terrain_write (GtsObject * o, FILE * fp)
{
  GfsRefineTerrain * t = GFS_REFINE_TERRAIN (o);
  (* GTS_OBJECT_CLASS (gfs_refine_terrain_class ())->parent_class->write) (o, fp);
  fprintf (fp, " %s ", t->name);
  gfs_function_write (t->criterion, fp);
  fprintf (fp, " { host = %s dbname = %s user = %s tables = %s }",
	   t->host, t->dbname, t->user, t->tables);
}

static void gfs_refine_terrain_class_init (GfsRefineClass * klass)
{
  klass->refine = terrain_refine;

  GTS_OBJECT_CLASS (klass)->destroy = refine_terrain_destroy;
  GTS_OBJECT_CLASS (klass)->read = refine_terrain_read;
  GTS_OBJECT_CLASS (klass)->write = refine_terrain_write;
}

static void gfs_refine_terrain_init (GfsRefineTerrain * t)
{
  t->criterion = gfs_function_new (gfs_function_class (), 0.);
  t->host =   g_strdup ("localhost");
  t->dbname = g_strdup ("earth");
  t->user =   g_strdup ("popinet");
  t->passwd = g_strdup ("CloClo");
  t->tables = g_strdup ("etopo2");
}

GfsRefineClass * gfs_refine_terrain_class (void)
{
  static GfsRefineClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_refine_terrain_info = {
      "GfsRefineTerrain",
      sizeof (GfsRefineTerrain),
      sizeof (GfsRefineClass),
      (GtsObjectClassInitFunc) gfs_refine_terrain_class_init,
      (GtsObjectInitFunc) gfs_refine_terrain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_refine_class ()),
				  &gfs_refine_terrain_info);
  }

  return klass;
}

/* Initialize module */

const gchar * g_module_check_init (void);
void          gfs_module_read     (GtsFile * fp);
void          gfs_module_write    (FILE * fp);

const gchar * g_module_check_init (void)
{
  gfs_refine_terrain_class ();
  return NULL;
}

void gfs_module_read (GtsFile * fp)
{
  g_return_if_fail (fp != NULL);
}

void gfs_module_write (FILE * fp)
{
  g_return_if_fail (fp != NULL);
}
