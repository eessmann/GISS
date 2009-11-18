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

#define GFS_PARTICLE(obj)            GTS_OBJECT_CAST (obj,\
					         GfsParticle,\
					         gfs_particle_class ())
#define GFS_IS_PARTICLE(obj)         (gts_object_is_from_class (obj,\
						 gfs_particle_class ()))

GfsEventClass * gfs_particle_class  (void);
