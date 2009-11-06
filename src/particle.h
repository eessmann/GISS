#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "ftt.h"
#include "domain.h"
#include "simulation.h"

/* Particle */
typedef struct _GfsParticle GfsParticle;
struct _GfsParticle{
  //  GfsEvent parent;
  GfsEvent parent;
  FttVector pos;
  guint id;
};

typedef struct _GfsParticleClass    GfsParticleClass;

struct _GfsParticleClass {
  GfsEventClass parent_class;
  /*< private >*/
  /*< public >*/
  /* add extra methods here */
};

#define GFS_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsParticle,\
					         gfs_particle_class ())
#define GFS_PARTICLE_CLASS(klass)    GTS_OBJECT_CLASS_CAST (klass,\
						 GfsParticleClass,\
						 gfs_particle_class())
#define IS_GFS_PARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_particle_class ()))

GfsParticleClass * gfs_particle_class  (void);
