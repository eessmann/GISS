typedef struct {
    gdouble * f, * m;
    GfsVariable * v;
    GfsFunction * weight;
    GfsSourceDiffusion * d;
    GfsDomain * domain;
  } ForceODE;

static void getboundary_ODE (FttCell * cell, gdouble * bdy)
{
  FttVector p;
  ftt_cell_pos (cell, &p);
 /* gdouble d = 1./pow(2,cell->parent->level+2);
  if (p.x + d > bdy[0]) bdy[0] = p.x + d;
    else if (p.x - d < bdy[1]) bdy[1] = p.x - d;
  if (p.y + d > bdy[2]) bdy[2] = p.y + d;
    else if (p.y - d < bdy[3]) bdy[3] = p.y - d;
  if (p.z + d > bdy[4]) bdy[4] = p.z + d;
    else if (p.z - d < bdy[5]) bdy[5] = p.z - d;*/
  if(cell->parent->level+1. > bdy[6]) bdy[6] = cell->parent->level+1.;
}

static void getforce_ODE (gint n, gint m, gdouble * bdy, FttVector * force)
{
  gdouble dis,dis2, dr, Force;
  FttVector r;
  dr = 2* L / pow(2,bdy[6]);
  force[m].x = force[m].y = force[m].z = 0.;
  if (n==m) return;
  if (BDInfo[n].type == 1 || BDInfo[n].type == 2) return;
  if (m == BDNUM) {
    dis = bdy[0] - BDInfo[n].c.x;
    if (dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else force[m].x = -1 * AF_Stf * pow((BDInfo[n].cradius + dr -dis)/dr, 2);
    return;
  }
  if (m == BDNUM +1) {
    dis = BDInfo[n].c.x- bdy[1];
    if (dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else force[m].x = AF_Stf * pow((BDInfo[n].cradius + dr -dis)/dr, 2);
    return;
  }
  if (m == BDNUM +2) {
    dis = bdy[2] - BDInfo[n].c.y;
    if (dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else force[m].y = -1 * AF_Stf * pow((BDInfo[n].cradius + dr -dis)/dr, 2);
    return;
  }
  if (m == BDNUM +3) {
    dis = BDInfo[n].c.y - bdy[3];
    if (dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else force[m].y =  AF_Stf * pow((BDInfo[n].cradius + dr -dis)/dr, 2);
    return;
  }
  if (m == BDNUM +4) {
    dis = bdy[4] - BDInfo[n].c.z;
    if (dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else force[m].z = -1 * AF_Stf * pow((BDInfo[n].cradius + dr -dis)/dr, 2);
    return;
  }
  if (m == BDNUM +5) {
    dis = BDInfo[n].c.z - bdy[5];
    if ( dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else force[m].z =  AF_Stf * pow((BDInfo[n].cradius + dr -dis)/dr, 2);
    return;
  }
  r.x = BDInfo[n].c.x - BDInfo[m].c.x;
  r.y = BDInfo[n].c.y - BDInfo[m].c.y;
  r.z = BDInfo[n].c.z - BDInfo[m].c.z;
  if (BDInfo[m].type == 0) {
    dis = sqrt(pow(r.x, 2.)+pow(r.y, 2.)+pow(r.z, 2.));
    if (dis > (BDInfo[n].cradius + BDInfo[m].cradius +dr)) Force = 0.;
      else Force = AF_Stf * pow((BDInfo[n].cradius + BDInfo[m].cradius + dr - dis)/dr,2) / dis ;
    force[m].x = Force * r.x, force[m].y = Force * r.y, force[m].z = Force * r.z;
    printf("fx=%g, fy=%g, fz=%g\n",force[m].x,force[m].y,force[m].z);
    return;
  }
  if (BDInfo[m].type == 1) {//SP_ODE: OutsideWall    
    r.x = (BDInfo[n].c.x - BDInfo[m].c.x) * (BDInfo[m].axis==0? 0.:1.);
    r.y = (BDInfo[n].c.y - BDInfo[m].c.y) * (BDInfo[m].axis==1? 0.:1.);
    r.z = (BDInfo[n].c.z - BDInfo[m].c.z) * (BDInfo[m].axis==2? 0.:1.);
    dis2 = sqrt(pow(r.x, 2.)+pow(r.y, 2.)+pow(r.z, 2.));
    dis = BDInfo[m].cradius-dis2;
    if (dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else Force = AF_Stf * pow((BDInfo[n].cradius + dr - dis)/dr,2) / dis2;
    force[m].x = Force * r.x, force[m].y = Force * r.y, force[m].z = Force * r.z;
    return;
  }
  if (BDInfo[m].type == 2) {//SP_ODE: InsideWall    
    r.x = (BDInfo[n].c.x - BDInfo[m].c.x) * (BDInfo[m].axis==0? 0.:1.);
    r.y = (BDInfo[n].c.y - BDInfo[m].c.y) * (BDInfo[m].axis==1? 0.:1.);
    r.z = (BDInfo[n].c.z - BDInfo[m].c.z) * (BDInfo[m].axis==2? 0.:1.);
    dis2 = sqrt(pow(r.x, 2.)+pow(r.y, 2.)+pow(r.z, 2.));
    dis = dis2 - BDInfo[m].cradius;
    if (dis > (BDInfo[n].cradius + dr)) Force = 0.;
      else Force = AF_Stf * pow((BDInfo[n].cradius + dr - dis)/dr,2) / dis2;
    force[m].x = Force * r.x, force[m].y = Force * r.y, force[m].z = Force * r.z;
    return;
  }
}

static void collision_ODE (gdouble * bdy)
{
  int m,n;
  gdouble dis, dr;
  const dReal *tmp;  
  FttVector v1, v2, v1t, v2t, c1,c2;
  dr = 2* L / pow(2,bdy[6]);
  for (m = 1; m < BDNUM; m++) {
    c1.x = BDInfo[m].c.x, c1.y = BDInfo[m].c.y, c1.z = BDInfo[m].c.z;
    tmp = dBodyGetLinearVel (BDInfo[m].body);
    v1.x = tmp[0], v1.y = tmp[1], v1.z = tmp[2];
    for(n = 1; n < BDNUM; n++) { 
      if(n!=m){	
        c2.x = BDInfo[n].c.x, c2.y = BDInfo[n].c.y, c2.z = BDInfo[n].c.z;
        tmp = dBodyGetLinearVel (BDInfo[n].body);
        v2.x = tmp[0], v2.y = tmp[1], v2.z = tmp[2];
        dis = sqrt(pow(c2.x-c1.x, 2.)+pow(c2.y-c1.y, 2.)+pow(c2.z-c1.z, 2.));
        if (dis <= (BDInfo[m].cradius + BDInfo[n].cradius +dr)){
	  v1t.x = (v1.x*pow(c2.y-c1.y,2) + v2.x*pow(c2.x-c1.x,2) + (v2.y-v1.y)*(c2.x-c1.x)*(c2.y-c1.y)) / (pow(c2.x-c1.x,2) + pow(c2.y-c1.y,2));
          v1t.y = (v1.y*pow(c2.x-c1.x,2) + v2.y*pow(c2.y-c1.y,2) + (v2.x-v1.x)*(c2.x-c1.x)*(c2.y-c1.y)) / (pow(c2.x-c1.x,2) + pow(c2.y-c1.y,2));
          v1t.z = v1.z;
	  v2t.x = (v1.x*pow(c2.x-c1.x,2) + v2.x*pow(c2.y-c1.y,2) - (v2.y-v1.y)*(c2.x-c1.x)*(c2.y-c1.y)) / (pow(c2.x-c1.x,2) + pow(c2.y-c1.y,2));
          v2t.y = (v1.y*pow(c2.y-c1.y,2) + v2.y*pow(c2.x-c1.x,2) - (v2.x-v1.x)*(c2.x-c1.x)*(c2.y-c1.y)) / (pow(c2.x-c1.x,2) + pow(c2.y-c1.y,2));
          v2t.z = v2.z;
          dBodySetLinearVel (BDInfo[m].body, v1t.x, v1t.y, v1t.z);
          dBodySetLinearVel (BDInfo[n].body, v2t.x, v2t.y, v2t.z);
        }
      }
    }
  }
}

static void add_force_ODE (FttCell * cell, ForceODE * f)
{
  gint num = cell->bdnum;
  if(BDInfo[num].type != 0) return;
  gdouble tmp[3];
  gdouble D;
  GfsSolidVector * s = GFS_STATE (cell)->solid;
  gdouble * r = &s->ca.x;
  FttVector ff, mm, n, g;	
  tmp[0]=r[0]*L-BDInfo[num].c.x, tmp[1]=r[1]*L-BDInfo[num].c.y, tmp[2]=r[2]*L-BDInfo[num].c.z;
  g_assert (((cell)->flags & GFS_FLAG_DIRICHLET) != 0);
  FttComponent c;
  GfsVariable ** v;
  v = gfs_domain_velocity (f->domain);
  f->v = gfs_variable_from_name (f->domain->variables, "P");
  gfs_pressure_force (cell, f->v, &ff);
  gts_vector_cross (&mm.x, tmp, &ff.x);
  BDInfo[num].pf.x += ff.x,  BDInfo[num].pf.y += ff.y,  BDInfo[num].pf.z += ff.z;
  BDInfo[num].pm.x += mm.x,  BDInfo[num].pm.y += mm.y,  BDInfo[num].pm.z += mm.z;
  for (c = 0; c < FTT_DIMENSION; c++) {
    ff.x = ff.y = ff.z = mm.x =mm.y = mm.z =0.;
    f->v = v[c];
    GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v[c]->surface_bc)->klass)->bc(cell,v[c]->surface_bc);
    f->d = source_diffusion (v[c]);
    gfs_cell_dirichlet_gradient (cell, f->v->i, -1, s->fv, &g);
    D = - gfs_source_diffusion_cell (f->d, cell);
    n.x = s->s[1] - s->s[0];
    n.y = s->s[3] - s->s[2];
  #if FTT_2D
    ff.z = 0.;
    switch (f->v->component) {
      case FTT_X:
        ff.x = D*(2.*g.x*n.x + g.y*n.y);
        ff.y = D*g.y*n.x;
        break;
      case FTT_Y:
        ff.x = D*g.x*n.y;
        ff.y = D*(2.*g.y*n.y + g.x*n.x);
        break;
      default:
        g_assert_not_reached ();
    }
  #else /* 3D */
    n.z = s->s[5] - s->s[4];
    D *= ftt_cell_size (cell);
    switch (f->v->component) {
      case FTT_X:
        ff.x = D*(2.*g.x*n.x + g.y*n.y + g.z*n.z);
        ff.y = D*g.y*n.x;
        ff.z = D*g.z*n.x;
        break;
      case FTT_Y:
        ff.y = D*(2.*g.y*n.y + g.x*n.x + g.z*n.z);
        ff.x = D*g.x*n.y;
        ff.z = D*g.z*n.y;
        break;
      case FTT_Z:
        ff.z = D*(2.*g.z*n.z + g.x*n.x + g.y*n.y);
        ff.x = D*g.x*n.z;
        ff.y = D*g.y*n.z;
        break;
      default:
        g_assert_not_reached (); 
    }
#endif /* 3D */
    gts_vector_cross (&mm.x, tmp, &ff.x);
    BDInfo[num].vf.x += ff.x,  BDInfo[num].vf.y += ff.y,  BDInfo[num].vf.z += ff.z;
    BDInfo[num].vm.x += mm.x,  BDInfo[num].vm.y += mm.y,  BDInfo[num].vm.z += mm.z;
  }
    //locate marked cell
    if(MKSWITCH==1){
      gdouble MKDis_n2c = sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);//disctance between new cell and solid center
      gdouble MKDis_n2m = sqrt((r[0]*L-Mark.x)*(r[0]*L-Mark.x)+(r[1]*L-Mark.y)*(r[1]*L-Mark.y)+(r[2]*L-Mark.z)*(r[2]*L-Mark.z));//disctance between new cell and marked cell
      if((MKDis_n2c >= MKDisOld)&&(MKDis_n2m<BDInfo[num].cradius)) {//new cell distance is larger than the old cell distance and new cell is closed to the marked cell
         MKDisOld=MKDis_n2c; //new cell becomes old cell	 
         MarkOld.x=r[0]*L, MarkOld.y=r[1]*L, MarkOld.z=r[2]*L;//new cell becomes old cell
      }//find the farthest close mark cell
    }//locate marked cell
}

static void gfs_domain_solid_force_ODE (GfsDomain * domain, GfsFunction * weight)
{
  gint n;
  const dReal *tmp;
  g_return_if_fail (domain != NULL);
  if (GFS_IS_AXI (domain))	g_assert_not_implemented ();
  for(n=1;n<BDNUM;n++)  //SPODE: set up the body info index for the followint force calculation
    {	tmp = dBodyGetPosition (BDInfo[n].body);
	BDInfo[n].c.x = tmp[0], BDInfo[n].c.y = tmp[1], BDInfo[n].c.z = tmp[2];
       	BDInfo[n].pf.x = BDInfo[n].pf.y = BDInfo[n].pf.z  = BDInfo[n].vf.x = BDInfo[n].vf.y = BDInfo[n].vf.z = 0.;
        BDInfo[n].pm.x = BDInfo[n].pm.y = BDInfo[n].pm.z  = BDInfo[n].vm.x = BDInfo[n].vm.y = BDInfo[n].vm.z = 0.;
    }
  ForceODE f;
  f.weight = weight;
  f.domain = domain;
  if (weight)	gfs_catch_floating_point_exceptions ();
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, (FttCellTraverseFunc) add_force_ODE, &f);
  if(MKSWITCH==1){
    Mark.x=MarkOld.x, Mark.y=MarkOld.y, Mark.z=MarkOld.z;
    MKDisOld = 0.95*MKDisOld;
      	BDInfo[1].Mark.x = Mark.x;
      	BDInfo[1].Mark.y = Mark.y;
      	BDInfo[1].Mark.z = Mark.z;    
  }
  if (weight)	gfs_restore_fpe_for_function (weight);
}
