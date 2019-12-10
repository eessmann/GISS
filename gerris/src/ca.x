advection.c:    gdouble s1 = s->solid ? s->solid->s[dleft] : 1., s2;
advection.c:      s2 = GFS_IS_MIXED (n) ? GFS_STATE (n)->solid->s[dleft] : 1.;
advection.c:	s2 = GFS_IS_MIXED (n) ? GFS_STATE (n)->solid->s[dleft]/2. : 0.5;
advection.c:      if (n.c[d] && !GFS_CELL_IS_BOUNDARY (n.c[d]) && solid->s[d] > 0. && 
advection.c:	  solid->a/solid->s[d] < GFS_SMALL)
advection.c:    solid->merged = NULL;
advection.c:      if (neighbor.c[i] && !GFS_CELL_IS_BOUNDARY (neighbor.c[i]) && solid->s[i] > 0.) {
advection.c:	    gdouble a = GFS_STATE (neighbor.c[i])->solid->a;
advection.c:	      solid->merged = neighbor.c[i];
advection.c:	    solid->merged = neighbor.c[i];
advection.c:		gdouble a = GFS_STATE (child.c[j])->solid->a;
advection.c:		  solid->merged = child.c[j];
advection.c:		solid->merged = child.c[j];
advection.c:		 solid->a);
advection.c:    if (solid && solid->merged)
advection.c:      add_merged (merged, solid->merged);
advection.c:		GFS_STATE (child.c[j])->solid->merged == cell)
advection.c:		 GFS_STATE (neighbor.c[i])->solid->merged == cell)
advection.c:	if (solid->s[d] > 0. && 1./solid->s[d] < mins)
advection.c:	  mins = 1./solid->s[d];
advection.c:	 solid->a, mins, 
advection.c:	 GFS_VALUE (cell, par->fv)/(mins*solid->a));
advection.c:      if (mins*solid->a > 0.01)
advection.c:	GFS_VALUE (cell, par->v) += GFS_VALUE (cell, par->fv)/(mins*solid->a);
advection.c:      gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
domain.c:  GFS_STATE (cell)->solid->fv = 0.;
domain.c:  GFS_STATE (cell)->solid->fv = 0.;
domain.c:  gts_range_add_value (s, GFS_STATE (cell)->solid->a);
domain.c:    a += GFS_IS_MIXED (c) ? GFS_STATE (c)->solid->a : 1.;
domain.c:      fprintf (fp, " %g", s->solid->s[i]);
domain.c:    fprintf (fp, " %g", s->solid->a);
domain.c:      fprintf (fp, " %g", (&s->solid->cm.x)[i]);
domain.c:    gts_file_error (fp, "expecting a number (solid->s[0])");
domain.c:    gts_file_error (fp, "solid->s[0] must be positive");
domain.c:    s->solid->s[0] = s0;
domain.c:	gts_file_error (fp, "expecting a number (solid->s[%d])", i);
domain.c:      s->solid->s[i] = atof (fp->token->str);
domain.c:      gts_file_error (fp, "expecting a number (solid->a)");
domain.c:    s->solid->a = atof (fp->token->str);
domain.c:	gts_file_error (fp, "expecting a number (solid->cm[%d])", i);
domain.c:      (&s->solid->cm.x)[i] = atof (fp->token->str);
domain.c:    fwrite (s->solid->s, sizeof (gdouble), FTT_NEIGHBORS, fp);
domain.c:    fwrite (&s->solid->a, sizeof (gdouble), 1, fp);
domain.c:    fwrite (&s->solid->cm.x, sizeof (gdouble), FTT_DIMENSION, fp);
domain.c:    fwrite (&s->solid->ca.x, sizeof (gdouble), FTT_DIMENSION, fp);
domain.c:    gts_file_error (fp, "expecting a number (solid->s[0])");
domain.c:    gts_file_error (fp, "solid->s[0] must be positive");
domain.c:    s->solid->s[0] = s0;
domain.c:    if (gts_file_read (fp, &s->solid->s[1], sizeof (gdouble), FTT_NEIGHBORS - 1) 
domain.c:      gts_file_error (fp, "expecting numbers (solid->s[1..%d])", FTT_NEIGHBORS - 1);
domain.c:    if (gts_file_read (fp, &s->solid->a, sizeof (gdouble), 1) != 1) {
domain.c:      gts_file_error (fp, "expecting a number (solid->a)");
domain.c:    if (gts_file_read (fp, &s->solid->cm.x, sizeof (gdouble), FTT_DIMENSION) != FTT_DIMENSION) {
domain.c:      gts_file_error (fp, "expecting numbers (solid->cm[0..%d])", FTT_DIMENSION - 1);
domain.c:	gts_file_read (fp, &s->solid->ca.x, sizeof (gdouble), FTT_DIMENSION) != FTT_DIMENSION) {
domain.c:      gts_file_error (fp, "expecting numbers (solid->ca[0..%d])", FTT_DIMENSION - 1);
domain.c:    gdouble * r = &GFS_STATE (cell)->solid->ca.x;
domain.c:	(!solid || solid->s[d] > 0.)) {
domain.c:	      (!GFS_IS_MIXED (child.c[i]) || GFS_STATE (child.c[i])->solid->s[od] > 0.)) {
domain.c:	  solid->s[d] > 0. && solid->s[d] < 1.) {
domain.c:  return ftt_cell_volume (cell)*(solid ? solid->a : 1.)*gfs_function_value (f, cell);
domain.c:  if (!n || GFS_CELL_IS_BOUNDARY (n) || (solid && solid->s[data->d] == 0.)) {
domain.c:	     (!GFS_IS_MIXED (n) || GFS_STATE (n)->solid->s[data->d] > 0.));
domain.h:  gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
domain.h:  double v = ftt_cell_volume (cell)*(GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.);
匹配到二进制文件 domain.o
event.c:    gts_range_add_value (vol, size*size*GFS_STATE (cell)->solid->a);
event.c:    GFS_VALUE (cell, div) += size*size*GFS_STATE (cell)->solid->a*(*ddiv);
fluid.c:  if (!GFS_IS_MIXED (cell) || GFS_STATE (cell)->solid->s[d] > 0.)
fluid.c:	gdouble w = GFS_IS_MIXED (children.c[i]) ? GFS_STATE (children.c[i])->solid->s[od] : 1.;
fluid.c:    if (!cell_bilinear (cell, n, &GFS_STATE (cell)->solid->ca, 
fluid.c:  o = &GFS_STATE (cell)->solid->cm;
fluid.c:      o = &GFS_STATE (cell)->solid->ca;
fluid.c:      v0 = GFS_STATE (cell)->solid->fv;
fluid.c:  o = &GFS_STATE (cell)->solid->cm;
fluid.c:      o = &GFS_STATE (cell)->solid->ca;
fluid.c:      v0 = GFS_STATE (cell)->solid->fv;
fluid.c:    if (!cell_bilinear (cell, n, &GFS_STATE (cell)->solid->ca, 
fluid.c:	gdouble w = GFS_IS_MIXED (children.c[i]) ? GFS_STATE (children.c[i])->solid->s[od] : 1.;
fluid.c:	  gdouble w = GFS_IS_MIXED (children.c[i]) ? GFS_STATE (children.c[i])->solid->s[od] : 1.;
fluid.h:                               GFS_STATE ((fa)->cell)->solid->s[(fa)->d] : 1.)
fluid.h:                 GFS_STATE ((fa)->neighbor)->solid->s[FTT_OPPOSITE_DIRECTION ((fa)->d)] : 1.)
function.h:  GFS_STATE (cell)->solid->fv = 0.;
function.h:  gfs_cell_dirichlet_gradient (_cell, v->i, -1, GFS_STATE (_cell)->solid->fv, &g);
function.h:    flux = GFS_STATE (_cell)->solid->fv;
function.h:    flux = gfs_cell_dirichlet_gradient_flux (_cell, v->i, -1, GFS_STATE (_cell)->solid->fv);
匹配到二进制文件 libgfs2D_la-domain.o
匹配到二进制文件 libgfs2D_la-solid.o
log:moving.c:  old_solid->coarse_fine = (GfsVariableFineCoarseFunc) none;
log:moving.c:  old_solid->fine_coarse = (GfsVariableFineCoarseFunc) none;
moving2.c:    if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0. && d1 == -1 && d2 == -1)
moving2.c:    else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0 && d2 == -1)
moving2.c:    else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0.)
moving2.c:	    if (GFS_STATE(neighbors.c[d1])->solid->s[d1] > 0 && GFS_STATE(neighbors.c[d1])->solid->s[d1] < 1)
moving2.c:	    if (GFS_STATE(neighbors.c[d2])->solid->s[d2] > 0 && GFS_STATE(neighbors.c[d2])->solid->s[d2] < 1)
moving2.c:	    if (GFS_STATE(neighbors.c[d])->solid->s[d1] != 1.) {
moving2.c:		  s2 = GFS_STATE(neighbors.c[d])->solid->s[d1];
moving2.c:  s1 = 1.-GFS_STATE (cell)->solid->s[d1];
moving2.c:	    if (GFS_STATE(neighbors.c[d])->solid->s[d1] == 0. && SOLD2 (neighbors.c[d], d1) != 1.) {
moving2.c:  s1 = GFS_STATE (cell)->solid->s[d1];
moving2.c:	if ((GFS_IS_MIXED(neighbors.c[d]) && GFS_STATE(neighbors.c[d])->solid->s[d1] == 1.) ||
moving2.c:	else if ((GFS_STATE(cell)->solid->s[d1] == 0. && GFS_IS_MIXED(neighbors.c[d])) ) {
moving2.c:	  s2 = 1.-GFS_STATE(neighbors.c[d])->solid->s[d1];
moving2.c:      if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0. &&
moving2.c:      else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0 && d2 == -1)
moving2.c:      else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0.)
moving2.c:    if (GFS_IS_MIXED(cell) && GFS_STATE(cell)->solid->s[do1] == 0.)
moving2.c:    if (GFS_IS_MIXED(cell) && GFS_STATE(cell)->solid->s[do2] == 0.)
moving2.c:      OLD_SOLID (cell)->s[d1] = GFS_STATE(cell)->solid->s[d1]*(dt1-1.); 
moving2.c:      OLD_SOLID (cell)->s[d1] = (dt1-1.)*GFS_STATE(cell)->solid->s[d1]+2.-dt1;
moving2.c:      OLD_SOLID (cell)->s[d2] = GFS_STATE(cell)->solid->s[d2]*(dt2-1.); 
moving2.c:      OLD_SOLID (cell)->s[d2] = (dt2-1.)*GFS_STATE(cell)->solid->s[d2]+2.-dt2;
moving2.c:	  old_solid->a = 1.;
moving2.c:	    old_solid->s[c] = 1.;
moving2.c:	old_solid->s[ftt_opposite_direction[2*i]] += OLD_SOLID (cell)->s[2*i];
moving2.c:	  old_solid->a = 1.;
moving2.c:	    old_solid->s[c] = 1.;
moving2.c:	old_solid->s[ftt_opposite_direction[2*i+1]] += OLD_SOLID (cell)->s[2*i+1];	
moving2.c:	OLD_SOLID (cell)->merged = GFS_STATE (cell)->solid->merged;
moving2.c:	  if (solid->s[c] == 0.)
moving2.c:	    solid_old->s[c] = (solid_old->s[c]+solid->s[c])/2. ;	
moving2.c:      OLD_SOLID (cell)->merged = GFS_STATE (cell)->solid->merged;
moving2.c:	if (solid->s[c] == 0.)
moving2.c:	  solid_old->s[c] = (solid_old->s[c]+solid->s[c])/2. ;
moving2.c:      OLD_SOLID (cell)->merged = GFS_STATE (cell)->solid->merged;
moving2.c:      if (GFS_STATE(cell)->solid->s[c] < 0.) {
moving2.c:	    GFS_STATE(cell)->solid->s[c] = OLD_SOLID (cell)->s[c];
moving2.c:	    GFS_STATE(cell)->solid->s[c] = 1.;
moving2.c:	  GFS_STATE(cell)->solid->s[c] = 0.;
moving2.c:	  if (GFS_STATE(cell)->solid->s[c] >= 0.)
moving2.c:	    OLD_SOLID (cell)->s[c] = GFS_STATE(cell)->solid->s[c];
moving.c:      oldca = GFS_STATE (cell)->solid->ca;
moving.c:      GFS_STATE (cell)->solid->ca = *ca;
moving.c:      GFS_STATE (cell)->solid->ca = oldca;
moving.c:    val = GFS_STATE (cell)->solid->fv;
moving.c:    solid->a = 0.;
moving.c:	SOLD2 (child.c[n], k) = solid->s[k] = 0.;
moving.c:  gint maxlevel = gfs_function_value (solid->level, cell);
moving.c:  solid->nvertex = gts_surface_vertex_number (GFS_SURFACE (GFS_SOLID (solid)->s)->s);
moving.c:      gfs_function_read (solid->level, gfs_object_simulation (*o), fp);
moving.c:  gfs_function_write (solid->level, fp);
moving.c:    gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
moving.c:      gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
moving.c:      gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
moving.c:  solid->level = gfs_function_new (gfs_function_class (), 0.);
moving.c:      p.stmp = g_array_set_size ( g_array_new (FALSE, FALSE, sizeof (double)) , FTT_DIMENSION*solid->nvertex);
moving.c:      p.sall = g_array_set_size ( g_array_new (FALSE, FALSE, sizeof (double)) , FTT_DIMENSION*solid->nvertex);
moving.c:    GFS_STATE (cell)->solid->fv*(GFS_STATE (cell)->solid->s[2*p->c + 1] -
moving.c:				 GFS_STATE (cell)->solid->s[2*p->c])*ftt_cell_size (cell);
moving.c:      gdouble a = GFS_STATE (cell)->solid ? GFS_STATE (cell)->solid->a : 1.;
moving.c:      gdouble a = GFS_STATE (cell)->solid ? GFS_STATE (cell)->solid->a : 1.;
moving.c:  gdouble a = GFS_STATE (cell)->solid ? GFS_STATE (cell)->solid->a : 1.;
moving.c:  old_solid->coarse_fine = (GfsVariableFineCoarseFunc) none;
moving.c:  old_solid->fine_coarse = (GfsVariableFineCoarseFunc) none;
moving.c:  old_solid->cleanup = (FttCellCleanupFunc) old_solid_cleanup;
ocean.c:    c *= GFS_STATE (cell)->solid->s[FTT_FRONT];
ocean.c:      gdouble area = h*h*solid->s[FTT_FRONT];
ocean.c:      s->f[FTT_FRONT].un = w = GFS_STATE (c)->solid->s[FTT_FRONT] > 0. ? 
ocean.c:	wf/GFS_STATE (c)->solid->s[FTT_FRONT] : 0.;
ocean.c:      wf += (solid->s[FTT_RIGHT]*s->f[FTT_RIGHT].un - solid->s[FTT_LEFT]*s->f[FTT_LEFT].un +
ocean.c:	     solid->s[FTT_TOP]*s->f[FTT_TOP].un - solid->s[FTT_BOTTOM]*s->f[FTT_BOTTOM].un);
ocean.c:  gdouble f = GFS_STATE (cell)->solid->s[FTT_FRONT];
ocean.c:  return GFS_STATE (cell)->solid->a/f*(1 << (MAXLEVEL - level));
poisson.c:    g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, p->maxlevel, GFS_STATE (cell)->solid->fv);
poisson.c:    GFS_STATE (cell)->solid->v = v;
poisson.c:      GFS_STATE (cell)->solid->v.x += alpha;
poisson.c:      GFS_STATE (cell)->solid->v.x = diffusion;
poisson.c:      f = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, solid->fv);
poisson.c:      f = solid->fv*solid->v.x;
poisson.c:      g.b = gfs_cell_dirichlet_gradient_flux (cell, p->u, -1, solid->fv);
poisson.c:      g.b = solid->fv*solid->v.x;
simulation.c:    FttVector p = GFS_STATE (cell)->solid->ca;
simulation.c:    FttVector p = GFS_STATE (cell)->solid->ca;
simulation.c:    FttVector p = GFS_STATE (cell)->solid->ca;
simulation.c:  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
simulation.c:  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_RIGHT] : 1.;
simulation.c:  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_LEFT] : 1.;
simulation.c:  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_TOP] : 1.;
simulation.c:  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_BOTTOM] : 1.;
simulation.c:  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_FRONT] : 1.;
simulation.c:  return GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->s[FTT_BACK] : 1.;
simulation.c:  m->x = m->y = GFS_STATE (cell)->solid->ca.y;
solid.c:  solid->a = 0.;
solid.c:  solid->cm.x = solid->cm.y = solid->cm.z = 0.;
solid.c:  solid->ca.z = 0.;
solid.c:      solid->s[etod[i]] = ins ? f->s[i].x : 1. - f->s[i].x;
solid.c:      solid->a += (x1 + x2)*(y2 - y1);
solid.c:      solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
solid.c:      solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
solid.c:	solid->a += (x1 + x2)*(y2 - y1);
solid.c:	solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
solid.c:	solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
solid.c:	solid->ca.x = (x1 + x2)/2.;
solid.c:	solid->ca.y = (y1 + y2)/2.;
solid.c:      solid->s[etod[i]] = 1.;
solid.c:      solid->a += (x1 + x2)*(y2 - y1);
solid.c:      solid->cm.x += (x1 - x2)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
solid.c:      solid->cm.y += (y2 - y1)*(2.*(x1*y1 + x2*y2) + x1*y2 + x2*y1);
solid.c:      solid->s[etod[i]] = 0.;
solid.c:  a = solid->a < 0. ? 0. : solid->a/(2.*h->x*h->y);
solid.c:  solid->ca.x += x0;
solid.c:  solid->ca.y += y0;
solid.c:    solid->cm.x = x0 + solid->cm.x/(3.*solid->a);
solid.c:    solid->cm.y = y0 + solid->cm.y/(3.*solid->a);
solid.c:    solid->cm.x = solid->cm.y = 0.;
solid.c:	solid->cm.x += x;
solid.c:	solid->cm.y += y;
solid.c:	  solid->cm.x += x1;
solid.c:	  solid->cm.y += y1;
solid.c:	solid->cm.x += x1;
solid.c:	solid->cm.y += y1;
solid.c:    solid->cm.x = x0 + solid->cm.x/n;
solid.c:    solid->cm.y = y0 + solid->cm.y/n;
solid.c:  solid->a = a;
solid.c:    if (solid->a == 1.) {
solid.c:    sum += GFS_STATE (cell)->solid->s[d];
solid.c:      solid->s[d] = (solid->s[d] > 0.5);
solid.c:    solid->a = 1.;
solid.c:    ftt_cell_pos (cell, &solid->cm);
solid.c:    solid->ca = solid->cm;
solid.c:  else if (GFS_STATE (cell)->solid && GFS_STATE (cell)->solid->a == 0.)
solid.c:      solid->s[i] = ins > 0 ? 0. : 1.;
solid.c:      solid->s[i] = sol.a;
solid.c:  solid->ca = ca;
solid.c:      (&m.x)[c] = solid->s[2*c + 1] - solid->s[2*c];
solid.c:	g_assert (solid->s[c] == solid->s[0]);
solid.c:      if (solid->s[0] == 1.) { /* fluid */
solid.c:	solid->a = 0.;
solid.c:	solid->cm.x = solid->cm.y = solid->cm.z = 0.;
solid.c:      solid->a = gfs_plane_volume (&m, alpha);
solid.c:      gfs_plane_center (&m, alpha, solid->a, &solid->cm);
solid.c:      (&solid->cm.x)[c] = (&o.x)[c] + 
solid.c:	(sym[c] ? 1. - (&solid->cm.x)[c] : (&solid->cm.x)[c])*(&h.x)[c];
solid.c:  if (solid->a == 0.)
solid.c:	w += solid->a; wa += sa;
solid.c:	cm.x += solid->cm.x*solid->a;
solid.c:	cm.y += solid->cm.y*solid->a;
solid.c:	cm.z += solid->cm.z*solid->a;
solid.c:	ca.x += solid->ca.x*sa;
solid.c:	ca.y += solid->ca.y*sa;
solid.c:	ca.z += solid->ca.z*sa;
solid.c:    solid->a = w/FTT_CELLS;
solid.c:    solid->ca.x = ca.x/wa;
solid.c:    solid->ca.y = ca.y/wa;
solid.c:    solid->ca.z = ca.z/wa;
solid.c:      solid->cm.x = cm.x/w;
solid.c:      solid->cm.y = cm.y/w;
solid.c:      solid->cm.z = cm.z/w;
solid.c:      ftt_cell_pos (cell, &solid->cm);
solid.c:	  w += GFS_IS_FLUID (child.c[j]) ? 1. : GFS_STATE (child.c[j])->solid->s[i];
solid.c:      solid->s[i] = w/n;
solid.c:	if (solid->s[i] == 0. || solid->s[i] == 1.) {
solid.c:	  push_leaf (fifo, n, i, solid->s[i] + 1., status);
solid.c:	  paint_leaf (fifo, solid->s[i] + 1., status);
solid.c:		GFS_STATE (child.c[j])->solid->s[FTT_OPPOSITE_DIRECTION (i)];
solid.c:	    g_warning ("file %s: line %d (%s): w/k=%g solid->s[%d]=%g",
solid.c:		       w/k, i, solid->s[i]);
solid.c:	  solid->s[i] = w/k;
solid.c:	      s += GFS_IS_MIXED (child.c[i]) ? GFS_STATE (child.c[i])->solid->s[od] : 1.;
solid.c:	  solid->s[d] = s/n;
solid.c:	    solid->s[d] = 1.;
solid.c:	    solid->s[d] = GFS_STATE (neighbor.c[d])->solid->s[FTT_OPPOSITE_DIRECTION (d)];
solid.c:	  solid->s[d] = 0.;
solid.c:    if (!gts_bbox_point_is_inside (&bb, &solid->cm)) {
solid.c:		 solid->cm.x, solid->cm.y, solid->cm.z, ftt_cell_level (root),
solid.c:    if (!gts_bbox_point_is_inside (&bb, &solid->ca)) {
solid.c:		 solid->ca.x, solid->ca.y, solid->ca.z, ftt_cell_level (root),
solid.c:	    if (1. - nsolid->s[oi] >= 1e-10) {
solid.c:			 oi, nsolid->s[oi]);
solid.c:	    nsolid->s[oi] = 1.;
solid.c:	  if (fabs (solid->s[i] - nsolid->s[oi]) >= 1e-10) {
solid.c:		       i, solid->s[i],
solid.c:		       oi, nsolid->s[oi]);
solid.c:	  nsolid->s[oi] = solid->s[i];
solid.c:	  if (1. - solid->s[i] >= 1e-10) {
solid.c:		       i, solid->s[i]);
solid.c:	  solid->s[i] = 1.;
solid.c:	    if (1. - solid->s[i] >= 1e-10) {
solid.c:			 i, solid->s[i]);
solid.c:	    solid->s[i] = 1.;
solid.c:	else if (nsolid->s[oi] == 0.) {
solid.c:	  if (solid->s[i] >= 1e-10) {
solid.c:		       i, solid->s[i]);
solid.c:	  solid->s[i] = 0.;
solid.c:		   n, GFS_STATE (children.c[n])->solid->a);
solid.c:	  a += GFS_STATE (children.c[n])->solid->a;
solid.c:    if (fabs (GFS_STATE (cell)->solid->a - a) >= 1e-10) {
solid.c:		 a, GFS_STATE (cell)->solid->a);
solid.c:    GFS_VALUE (cell, c) = solid->a;
solid.c:    *cm = GFS_STATE (cell)->solid->cm;
solid.c:  m.x = solid->s[2*c1 + 1] - solid->s[2*c1]; nm = fabs (m.x);
solid.c:  m.y = solid->s[2*c2 + 1] - solid->s[2*c2]; nm += fabs (m.y);
solid.c:  gdouble alpha = gfs_line_alpha (&m, solid->s[d]);
solid.c:  if (fabs (solid->s[d] - ss/n) > 1e-5)
solid.c:    g_warning ("inconsistent surface fractions %d %f %f %f\n", d, solid->s[d], ss/n,
solid.c:	       fabs (solid->s[d] - ss/n));
solid.c:    (&m.x)[c] = solid->s[2*c + 1] - solid->s[2*c];
solid.c:  gdouble alpha = gfs_plane_alpha (&m, solid->a);
solid.c:	  s->s[d] = GFS_STATE (n.c[d])->solid->s[FTT_OPPOSITE_DIRECTION (d)];
匹配到二进制文件 solid.o
source.c:      f = gfs_cell_dirichlet_gradient_flux (cell, p->v->i, -1, solid->fv);
source.c:      f = solid->fv*solid->v.x;
source.c:  gdouble a = GFS_IS_MIXED (cell) ? GFS_STATE (cell)->solid->a : 1.;
timestep.c:    GFS_STATE (cell)->solid->fv = gfs_function_value (bc->val, cell);
timestep.c:    GFS_STATE (cell)->solid->fv = gfs_function_value (bc->val, cell)*ftt_vector_norm (&n)
