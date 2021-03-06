/* GfsOutputSolidMovingODE: Object */

GfsOutputClass * gfs_output_solid_moving_ode_class (void);

static gboolean gfs_output_solid_moving_ode_event (GfsEvent * event, 
					      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (gfs_output_class())->event) (event, sim)
	&& sim->advection_params.dt > 0.) {
  FILE * fp = GFS_OUTPUT (event)->file->fp;
  gint n;
  static gint IniKill=0;
  if (GFS_OUTPUT (event)->first_call) {
    fputs ("#", fp);
    for (n=1;n<BDNUM;n++) 
      fprintf (fp, "%d_T1 cx2 cy3 cz4 vx5 vy6 vz7 ex8 ey9 ez10 ox11 oy12 oz13 pfx14 pfy15 pfz16 vfx17 vfy18 vfz19 pmx20 pmy21 pmz22 vmx23 vmy24 vmz25 ", n);
      if(MKSWITCH==1) fprintf (fp, "mkx26 mky27 mkz28 ");
    fputs ("\n", fp);
  }
    if (K>0 && Step_SP==0.)return TRUE;
    if (IniKill==0){IniKill++; return TRUE;}
    for (n=1;n<BDNUM;n++) {
	if(BDInfo[n].OPInfo.buff_count==0) 
	{fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ",
        sim->time.t, 
	BDInfo[n].OPInfo.c.x, BDInfo[n].OPInfo.c.y, BDInfo[n].OPInfo.c.z,
	BDInfo[n].OPInfo.v.x, BDInfo[n].OPInfo.v.y, BDInfo[n].OPInfo.v.z,
	BDInfo[n].OPInfo.e.x, BDInfo[n].OPInfo.e.y, BDInfo[n].OPInfo.e.z,
	BDInfo[n].OPInfo.o.x, BDInfo[n].OPInfo.o.y, BDInfo[n].OPInfo.o.z,
	BDInfo[n].OPInfo.pf.x * Ln, BDInfo[n].OPInfo.pf.y * Ln, BDInfo[n].OPInfo.pf.z * Ln, 
	BDInfo[n].OPInfo.vf.x * Ln, BDInfo[n].OPInfo.vf.y * Ln, BDInfo[n].OPInfo.vf.z * Ln, 
	BDInfo[n].OPInfo.pm.x * Ln, BDInfo[n].OPInfo.pm.y * Ln, BDInfo[n].OPInfo.pm.z * Ln,
	BDInfo[n].OPInfo.vm.x * Ln, BDInfo[n].OPInfo.vm.y * Ln, BDInfo[n].OPInfo.vm.z * Ln);//25
	if(MKSWITCH==1) fprintf (fp, "%g %g %g ",BDInfo[n].OPInfo.Mark.x, BDInfo[n].OPInfo.Mark.y, BDInfo[n].OPInfo.Mark.z);
        if (n == BDNUM-1) fputs ("\n", fp);}
     }
    return TRUE;
  }
  return FALSE;
}

static void gfs_output_solid_moving_ode_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_solid_moving_ode_event;
}

GfsOutputClass * gfs_output_solid_moving_ode_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_solid_moving_ode_info = {
      "GfsOutputSolidMovingODE",
      sizeof (GfsOutput),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_solid_moving_ode_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_solid_moving_ode_info);
  }

  return klass;
}

