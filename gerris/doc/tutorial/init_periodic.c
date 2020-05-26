#include <math.h>
#include <stdlib.h>
#include <gfs.h>

/* InitPeriodic: Header */

typedef struct _InitPeriodic         InitPeriodic;

struct _InitPeriodic {
  /*< private >*/
  GfsInit parent;

  /*< public >*/
  gdouble m;
};

typedef struct _InitPeriodicClass    InitPeriodicClass;

struct _InitPeriodicClass {
  /*< private >*/
  GfsInitClass parent_class;

  /*< public >*/
  /* add extra methods here */
};

#define INIT_PERIODIC(obj)            GTS_OBJECT_CAST (obj,\
					         InitPeriodic,\
					         init_periodic_class ())
#define INIT_PERIODIC_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 InitPeriodicClass,\
						 init_periodic_class())
#define IS_INIT_PERIODIC(obj)         (gts_object_is_from_class (obj,\
						 init_periodic_class ()))

InitPeriodicClass * init_periodic_class  (void);
InitPeriodic * init_periodic_new    (InitPeriodicClass * klass);

/* InitPeriodic: Object */

static void init_periodic_read (GtsObject ** o, GtsFile * fp)
{
  /* call read method of parent */
  if (GTS_OBJECT_CLASS (init_periodic_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (init_periodic_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (m)");
    return;
  }
  INIT_PERIODIC (*o)->m = atof (fp->token->str);

  /* do not forget to prepare for next read */
  gts_file_next_token (fp);
}

static void init_periodic_write (GtsObject * o, FILE * fp)
{
  /* call write method of parent */
  if (GTS_OBJECT_CLASS (init_periodic_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (init_periodic_class ())->parent_class->write) 
      (o, fp);
  
  fprintf (fp, " %g", INIT_PERIODIC (o)->m);
}

static void init_velocity (FttCell * cell,
                           InitPeriodic * init)
{
  FttVector pos;

  ftt_cell_pos (cell, &pos);
  GFS_STATE (cell)->u = 
    - cos (2.*init->m*M_PI*pos.x)*sin (2.*init->m*M_PI*pos.y);
  GFS_STATE (cell)->v =   
    sin (2.*init->m*M_PI*pos.x)*cos (2.*init->m*M_PI*pos.y);
}

static gboolean init_periodic_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (init_periodic_class ())->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
                              FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                              (FttCellTraverseFunc) init_velocity,
                              event);
    return TRUE;
  }
  return FALSE;
}

static void init_periodic_class_init (InitPeriodicClass * klass)
{
  /* define new methods and overload inherited methods here */

  GFS_EVENT_CLASS (klass)->event = init_periodic_event;
  GTS_OBJECT_CLASS (klass)->read = init_periodic_read;
  GTS_OBJECT_CLASS (klass)->write = init_periodic_write;
}

static void init_periodic_init (InitPeriodic * object)
{
  object->m = 1.;
}

InitPeriodicClass * init_periodic_class (void)
{
  static InitPeriodicClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo init_periodic_info = {
      "InitPeriodic",
      sizeof (InitPeriodic),
      sizeof (InitPeriodicClass),
      (GtsObjectClassInitFunc) init_periodic_class_init,
      (GtsObjectInitFunc) init_periodic_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_init_class ()),
				  &init_periodic_info);
  }

  return klass;
}

InitPeriodic * init_periodic_new (InitPeriodicClass * klass)
{
  InitPeriodic * object;

  object = INIT_PERIODIC (gts_object_new (GTS_OBJECT_CLASS (klass)));

  return object;
}

/* Initialize module */

const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  init_periodic_class ();
  return NULL;
}
