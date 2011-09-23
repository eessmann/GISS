#include "output.h"
#include "cartesian.h"

#include <fftw3.h>
#include <stdlib.h>

/* GfsOutputSpectra: header */

typedef struct _GfsOutputSpectra                     GfsOutputSpectra;

struct _GfsOutputSpectra {
  /*< private >*/
  GfsOutput parent;

  /*< public >*/
  GfsVariable * v;
  FttVector L;
  GfsCartesianGrid cgd;
};


#define GFS_OUTPUT_SPECTRA(obj)            GTS_OBJECT_CAST (obj,\
    GfsOutputSpectra,\
    gfs_output_spectra_class ())
#define GFS_IS_OUTPUT_SPECTRA(obj)         (gts_object_is_from_class (obj,\
      gfs_output_spectra_class ()))

GfsOutputClass * gfs_output_spectra_class  (void);

/** \beginobject{GfsOutputSpectra} */

static void fill_cartesian_matrix_2D ( GfsCartesianGrid * cgd, GfsVariable * v, GfsDomain * domain )
{
  guint i,j;
  FttVector pos;
  FttCell * cell;

  for (i = 0; i < cgd->n[0]; i++) {
    for (j = 0; j < cgd->n[1]; j++) {
      pos.x = cgd->x[0][i];
      pos.y = cgd->x[1][j];
      pos.z = 0;

      cell = gfs_domain_locate (domain, pos, -1, NULL);
      cgd->v[i*cgd->n[1]+j] = gfs_interpolate (cell, pos, v);
    }
  }
}

static void fill_cartesian_matrix_3D ( GfsCartesianGrid * cgd, GfsVariable * v, GfsDomain * domain )
{
  guint i,j,k;
  FttVector pos;
  FttCell * cell;

  for (i = 0; i < cgd->n[0]; i++) {
    for (j = 0; j < cgd->n[1]; j++) {
      for (k = 0; k < cgd->n[2]; k++) {
        pos.x = cgd->x[0][i];
        pos.y = cgd->x[1][j];
        pos.z = cgd->x[2][j];

        cell = gfs_domain_locate (domain, pos, -1, NULL);
        cgd->v[k+cgd->n[2]*(i*cgd->n[1]+j)] = gfs_interpolate (cell, pos, v);
      }
    }
  }
}

static void write_spectra_2D ( GfsSimulation * sim, FILE * fp, fftw_complex *out, FttVector L, guint nx, guint ny){

  guint i,j; 
  gint aux;
  gdouble kx,ky;
  fputs ("# 1:kx 2:ky 3:kz 4:real 5:img \n", fp);
  gdouble kxmax = 2.*M_PI/L.x; 
  gdouble kymax = 2.*M_PI/L.y;

  for ( i = 0; i < nx; i++ )
  {
    if ( i < nx/2. +1 ) kx = kxmax*i;
    else {
      aux = i-nx;
      kx = kxmax*aux;
    }
    for ( j = 0; j < ny; j++ )
    {
      ky = kymax*j;
      fprintf (fp, "  %g  %g %g  %g  %g \n", 
          kx, ky, 0.0 , out[i*ny+j][0], out[i*ny+j][1] );
    }
  }
}

static void write_spectra_3D ( GfsSimulation * sim, FILE * fp, fftw_complex *out, FttVector L, guint nx, guint ny , guint nz){

  guint i,j,k; 
  gint aux;
  gdouble kx,ky,kz;
  fputs ("# 1:kx 2:ky 3:kz 4:real 5:img \n", fp);
  gdouble kxmax = 2.*M_PI/L.x; 
  gdouble kymax = 2.*M_PI/L.y;
  gdouble kzmax = 2.*M_PI/L.z;

  for ( i = 0; i < nx; i++ )
  {
    if ( i < nx/2. +1 ) kx = kxmax*i;
    else { 
      aux = i-nx;
      kx = kxmax*aux;
    }
    for ( j = 0; j < ny; j++ )
    {
      if ( j < ny/2. +1 ) ky = kymax*j;
      else {
        aux = j-ny;
        ky = kymax*aux;
      }
      for ( k = 0; k < nz; k++ )
        ky = kzmax*(gdouble) k;
      fprintf (fp, "  %g  %g %g  %g  %g \n", 
          kx, ky, kz , out[k+nz*(i*ny+j)][0], out[k+nz*(i*ny+j)][1] );
    }
  }
}

static gboolean output_spectra_event (GfsEvent * event, 
    GfsSimulation * sim) 
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (event);
    fftw_plan p;

#if FTT_2D
    fill_cartesian_matrix_2D( &(v->cgd), v->v, domain);
    guint nyh = ( v->cgd.n[1] / 2 ) + 1;
    fftw_complex *out =  fftw_malloc( sizeof(fftw_complex)*v->cgd.n[0]*nyh );
    p = fftw_plan_dft_r2c_2d( v->cgd.n[0], v->cgd.n[1], v->cgd.v, out,  FFTW_ESTIMATE);
    fftw_execute(p); 
    write_spectra_2D ( sim, GFS_OUTPUT (event)->file->fp, out, v->L, v->cgd.n[0], nyh);

#else
    fill_cartesian_matrix_3D( &(v->cgd), v->v, domain);
    guint nzh = ( v->cgd.n[2] / 2 ) + 1;
    fftw_complex *out =  fftw_malloc( sizeof(fftw_complex)*v->cgd.n[0]*v->cgd.n[1]*nzh );
    p = fftw_plan_dft_r2c_3d( v->cgd.n[0], v->cgd.n[1], v->cgd.n[2], v->cgd.v, out,  FFTW_ESTIMATE);
    fftw_execute(p); 
    write_spectra_3D ( sim, GFS_OUTPUT (event)->file->fp, out, v->L, v->cgd.n[0], v->cgd.n[1], nzh );

#endif /* FTT_3D */

    fftw_destroy_plan(p);
    fftw_free ( out );

    return TRUE;
  }
  return FALSE;
}

static void output_spectra_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;
  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (*o);
  FttVector pos;
  guint level;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }

  gts_file_next_token (fp);

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
    else if (!strcmp (fp->token->str, "x")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      pos.x = atof(fp->token->str);
    }
    else if (!strcmp (fp->token->str, "y")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      pos.y = atof(fp->token->str);
    }
    else if (!strcmp (fp->token->str, "z")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      pos.z = atof(fp->token->str);
    }
    else if (!strcmp (fp->token->str, "Lx")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      v->L.x = atof(fp->token->str);
    }
    else if (!strcmp (fp->token->str, "Ly")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      v->L.y = atof(fp->token->str);
    }
    else if (!strcmp (fp->token->str, "Lz")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      v->L.z = atof(fp->token->str);
    }
    else {
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
      return;
    }
    gts_file_next_token (fp);
  }
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integel (level)");
    return FALSE;
  }
  level = atoi(fp->token->str);
  gts_file_next_token (fp);

  guint i,j,size =1;
  v->cgd.n = g_malloc (v->cgd.N*sizeof (guint));  
  for (i = 0; i < v->cgd.N; i++) {
    v->cgd.n[i] = pow(2,level);
    size *= v->cgd.n[i];
  }

  //fixme: Generalize for tilted planes
  v->cgd.x = g_malloc0 (v->cgd.N*sizeof (gdouble *));
  for (i = 0; i < v->cgd.N; i++) {
    v->cgd.x[i] = g_malloc (v->cgd.n[i]*sizeof (gdouble));
    for (j = 0; j < v->cgd.n[i]; j++){ 
      v->cgd.x[i][j] = (&(pos.x))[i] + (&(v->L.x))[i]*(gdouble)j/((gdouble)(v->cgd.n[i]-1))-0.5;
    }
  }

  v->cgd.v = g_malloc0( sizeof ( gdouble ) * 2*(size/2+1)  );

}

static void output_spectra_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class->write) (o, fp); 

  fprintf (fp, " %s ", GFS_OUTPUT_SPECTRA (o)->v->name );
}

static void output_spectra_destroy ( GtsObject * o ) {

  guint i;
  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (o);
  //fixme: call gfs_cartesian_grid_destroy if possible
}

static void output_spectra_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = output_spectra_event;
  klass->read =  output_spectra_read;
  klass->write = output_spectra_write;
  klass->destroy = output_spectra_destroy;
}

GfsOutputClass * gfs_output_spectra_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_spectra_info = {
      "GfsOutputSpectra",
      sizeof (GfsOutputSpectra),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) output_spectra_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
        &gfs_output_spectra_info);
  }

  return klass;
}

/** \endobject{GfsOutputSpectra} */

/* GfsOutputSpectra2D: header */

#define GFS_OUTPUT_SPECTRA_2D(obj)            GTS_OBJECT_CAST (obj,\
    GfsOutputSpectra2D,\
    gfs_output_spectra_2D_class ())
#define GFS_IS_OUTPUT_SPECTRA_2D(obj)         (gts_object_is_from_class (obj,\
      gfs_output_spectra_2D_class ()))

GfsOutputClass * gfs_output_spectra_2D_class  (void);

/** \beginobject{GfsOutputSpectra2D} */

static void output_spectra_2D_init (GfsOutputSpectra *o)
{
  o->cgd.N = 2;
}

GfsOutputClass * gfs_output_spectra_2D_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_spectra_2D_info = {
      "GfsOutputSpectra2D",
      sizeof (GfsOutputSpectra),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) output_spectra_2D_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_spectra_class ()),
        &gfs_output_spectra_2D_info);
  }

  return klass;
}

/** \endobject{GfsOutputSpectra2D} */

/* GfsOutputSpectra3D: header */

#define GFS_OUTPUT_SPECTRA_3D(obj)            GTS_OBJECT_CAST (obj,\
    GfsOutputSpectra3D,\
    gfs_output_spectra_3D_class ())
#define GFS_IS_OUTPUT_SPECTRA_3D(obj)         (gts_object_is_from_class (obj,\
      gfs_output_spectra_3D_class ()))

GfsOutputClass * gfs_output_spectra_3D_class  (void);

/** \beginobject{GfsOutputSpectra3D} */

static void output_spectra_3D_init (GfsOutputSpectra *o)
{
  o->cgd.N = 3;
}

GfsOutputClass * gfs_output_spectra_3D_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_spectra_3D_info = {
      "GfsOutputSpectra3D",
      sizeof (GfsOutputSpectra),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) output_spectra_3D_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_spectra_class ()),
        &gfs_output_spectra_3D_info);
  }

  return klass;
}

/** \endobject{GfsOutputSpectra3D} */


/* Initialize module */

const gchar gfs_module_name[] = "fourier";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_output_spectra_class ();
  gfs_output_spectra_2D_class ();
  gfs_output_spectra_3D_class ();
  return NULL; 
}
