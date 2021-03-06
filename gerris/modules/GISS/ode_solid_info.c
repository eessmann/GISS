  if (fp->type == '{') {    
    /* variables definition */
    gdouble ex = 0., ey = 0., ez = 0.; //SP_ODE: initiating rotation quaternion 
    gdouble themass = 1.; //SP_ODE: mass of solid  
    gdouble I11 = 1., I22 = 1., I33 = 1., I12 = 0., I13 = 0., I23 = 0.;  //SP_ODE: moment of inertia tensor of solid  
    BDInfo[BDNUM].c.x = BDInfo[BDNUM].c.y = BDInfo[BDNUM].c.z = 0.;
    BDInfo[BDNUM].a_ex.x = BDInfo[BDNUM].a_ex.y = BDInfo[BDNUM].a_ex.z = 0.;
    BDInfo[BDNUM].m_ex.x = BDInfo[BDNUM].m_ex.y = BDInfo[BDNUM].m_ex.z = 0.;
    BDInfo[BDNUM].f_tx = BDInfo[BDNUM].f_ty = BDInfo[BDNUM].f_tz = 1;
    BDInfo[BDNUM].f_rx = BDInfo[BDNUM].f_ry = BDInfo[BDNUM].f_rz = 1;
    BDInfo[BDNUM].cradius = 0.;
    BDInfo[BDNUM].count = 0;
    BDInfo[BDNUM].type = 0;
    BDInfo[BDNUM].axis = 0;
    gdouble vx = 0., vy = 0., vz = 0.; //SP_ODE: initial linear velocity of solid 
    gdouble ox = 0., oy = 0., oz = 0.; //SP_ODE: initial angular velocity of solid 
    BDInfo[BDNUM].OPInfo.buff_count = 0;
    BDInfo[BDNUM].OPInfo.pf.x = BDInfo[BDNUM].OPInfo.pf.y = BDInfo[BDNUM].OPInfo.pf.z = 0.;
    BDInfo[BDNUM].OPInfo.vf.x = BDInfo[BDNUM].OPInfo.vf.y = BDInfo[BDNUM].OPInfo.vf.z = 0.;
    BDInfo[BDNUM].OPInfo.f.x = BDInfo[BDNUM].OPInfo.f.y = BDInfo[BDNUM].OPInfo.f.z = 0.;
    BDInfo[BDNUM].OPInfo.pm.x = BDInfo[BDNUM].OPInfo.pm.y = BDInfo[BDNUM].OPInfo.pm.z = 0.;
    BDInfo[BDNUM].OPInfo.vm.x = BDInfo[BDNUM].OPInfo.vm.y = BDInfo[BDNUM].OPInfo.vm.z = 0.;
    BDInfo[BDNUM].OPInfo.m.x = BDInfo[BDNUM].OPInfo.m.y = BDInfo[BDNUM].OPInfo.m.z = 0.;
    BDInfo[BDNUM].OPInfo.c.x = BDInfo[BDNUM].OPInfo.c.y = BDInfo[BDNUM].OPInfo.c.z = 0.;
    BDInfo[BDNUM].OPInfo.v.x = BDInfo[BDNUM].OPInfo.v.y = BDInfo[BDNUM].OPInfo.v.z = 0.;
    BDInfo[BDNUM].OPInfo.e.x = BDInfo[BDNUM].OPInfo.e.y = BDInfo[BDNUM].OPInfo.e.z = 0.;
    BDInfo[BDNUM].OPInfo.o.x = BDInfo[BDNUM].OPInfo.o.y = BDInfo[BDNUM].OPInfo.e.z = 0.;
    BDInfo[BDNUM].OPInfo.Mark.x = BDInfo[BDNUM].OPInfo.Mark.y = BDInfo[BDNUM].OPInfo.Mark.z = 0.;
    dQuaternion q; //SP_ODE: rotation matrix 
    Mark.x=Mark.y=Mark.z=0.;
    /* read from script file */
    GtsFileVariable var[] = {	//SP_ODE: reading solid parameters
      {GTS_DOUBLE, "cx", TRUE, &(BDInfo[BDNUM].c.x)}, {GTS_DOUBLE, "cy", TRUE, &(BDInfo[BDNUM].c.y)}, {GTS_DOUBLE, "cz", TRUE, &(BDInfo[BDNUM].c.z)},
      {GTS_DOUBLE, "ex", TRUE, &ex}, {GTS_DOUBLE, "ey", TRUE, &ey}, {GTS_DOUBLE, "ez", TRUE, &ez},
      {GTS_DOUBLE, "vx", TRUE, &vx}, {GTS_DOUBLE, "vy", TRUE, &vy}, {GTS_DOUBLE, "vz", TRUE, &vz},
      {GTS_DOUBLE, "ox", TRUE, &ox}, {GTS_DOUBLE, "oy", TRUE, &oy}, {GTS_DOUBLE, "oz", TRUE, &oz},
      {GTS_DOUBLE, "ax", TRUE, &(BDInfo[BDNUM].a_ex.x)}, {GTS_DOUBLE, "ay", TRUE, &(BDInfo[BDNUM].a_ex.y)}, {GTS_DOUBLE, "az", TRUE, &(BDInfo[BDNUM].a_ex.z)},
      {GTS_DOUBLE, "mx", TRUE, &(BDInfo[BDNUM].m_ex.x)}, {GTS_DOUBLE, "my", TRUE, &(BDInfo[BDNUM].m_ex.y)}, {GTS_DOUBLE, "mz", TRUE, &(BDInfo[BDNUM].m_ex.z)},
      {GTS_DOUBLE, "mass", TRUE, &themass},      {GTS_DOUBLE, "Step_SP", TRUE, &Step_SP},
      {GTS_DOUBLE, "I11", TRUE, &I11}, {GTS_DOUBLE, "I22", TRUE, &I22}, {GTS_DOUBLE, "I33", TRUE, &I33},
      {GTS_DOUBLE, "I12", TRUE, &I12}, {GTS_DOUBLE, "I13", TRUE, &I13}, {GTS_DOUBLE, "I23", TRUE, &I23},
      {GTS_INT, "f_tx", TRUE, &(BDInfo[BDNUM].f_tx)}, {GTS_INT, "f_ty", TRUE, &(BDInfo[BDNUM].f_ty)}, {GTS_INT, "f_tz", TRUE, &(BDInfo[BDNUM].f_tz)},
      {GTS_INT, "f_rx", TRUE, &(BDInfo[BDNUM].f_rx)}, {GTS_INT, "f_ry", TRUE, &(BDInfo[BDNUM].f_ry)}, {GTS_INT, "f_rz", TRUE, &(BDInfo[BDNUM].f_rz)},
      {GTS_DOUBLE, "cradius", TRUE, &(BDInfo[BDNUM].cradius)}, {GTS_INT, "type", TRUE, &(BDInfo[BDNUM].type)}, {GTS_INT, "axis", TRUE, &(BDInfo[BDNUM].axis)},
      {GTS_INT, "K", TRUE, &K},{GTS_INT, "BUFF", TRUE, &BUFF},{GTS_INT, "AF_Flag", TRUE, &AF_Flag},{GTS_DOUBLE, "AF_Stf", TRUE, &AF_Stf},
      {GTS_DOUBLE, "mk_x", TRUE, &(Mark.x)}, {GTS_DOUBLE, "mk_y", TRUE, &(Mark.y)}, {GTS_DOUBLE, "mk_z", TRUE, &(Mark.z)},{GTS_INT, "MKSWITCH", TRUE, &MKSWITCH},
      {GTS_NONE}
    };
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR) return;
    ex = ex /180.*PI, ey = ey /180.*PI, ez = ez /180.*PI;
    q[0] = cos(ex/2.) * cos(ey/2.) * cos(ez/2.) + sin(ex/2.) * sin(ey/2.) * sin(ez/2.),
    q[1] = sin(ex/2.) * cos(ey/2.) * cos(ez/2.) - cos(ex/2.) * sin(ey/2.) * sin(ez/2.),
    q[2] = cos(ex/2.) * sin(ey/2.) * cos(ez/2.) + sin(ex/2.) * cos(ey/2.) * sin(ez/2.); 
    q[3] = cos(ex/2.) * cos(ey/2.) * sin(ez/2.) - sin(ex/2.) * sin(ey/2.) * cos(ez/2.);   
    //locate marked cell
    if(MKSWITCH==1){
      MarkOld.x = Mark.x, MarkOld.y = Mark.y, MarkOld.z = MarkOld.z;
      MKDisOld = 0.95*sqrt((Mark.x-BDInfo[BDNUM].c.x)*(Mark.x-BDInfo[BDNUM].c.x)+(Mark.y-BDInfo[BDNUM].c.y)*(Mark.y-BDInfo[BDNUM].c.y)+(Mark.z-BDInfo[BDNUM].c.z)*(Mark.z-BDInfo[BDNUM].c.z));
    }//locate marked cell
    dBodySetPosition (BDInfo[BDNUM].body, BDInfo[BDNUM].c.x, BDInfo[BDNUM].c.y, BDInfo[BDNUM].c.z); //SP_ODE: set initial position     

    dBodySetQuaternion (BDInfo[BDNUM].body, q); //SP_ODE: set initial angular quaternion      
    
    dBodySetLinearVel (BDInfo[BDNUM].body, vx, vy, vz);  //SP_ODE: set initial linear velocity    
     
    dBodySetAngularVel (BDInfo[BDNUM].body, ox, oy, oz);  //SP_ODE: set initial angular velocity 

    dMassSetParameters (&BDInfo[BDNUM].mass, themass, 0., 0., 0., I11, I22, I33, I12, I13, I23); //SP_ODE: the center of gravity in body frame is (0,0,0)     
    
    dBodySetMass (BDInfo[BDNUM].body, &BDInfo[BDNUM].mass);  //SP_ODE: set mass and moment of inertia 
  }
