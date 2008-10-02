/* ----- RSTQuery.c ----- */


#include "RStarTree.h"
#include "RSTQuery.h"
#include "RSTUtil.h"
#include "RSTInOut.h"


/* declarations */

/************************************************************************/

boolean FoundRect(RSTREE R,
                  int depth,
                  typrect rectangle,
                  boolean isinsert,
                  refinfo *infoadr)

{
  refparameters par;
  refcount c;
  refDATAnode n;
  int instind;
  boolean found;
  int i;
  
  i= -1; instind= -1; found= FALSE;
  
  par= &(*R).parameters._;
  
  if (depth != (*par).height) {
    if (isinsert) {
      ChooseSubtree(R,rectangle,depth,&(*(*R).N[depth]).DIR,&instind);
      (*R).EInst[depth]= instind;
    }
    do {
      i++;
      if (Covers(R,(*(*R).N[depth]).DIR.entries[i].rect,rectangle)) {
        (*R).E[depth]= i;
        depth++;
        if ((*R).N[depth] == (*R).NInst[depth]) {
          if ((*R).Nmodified[depth]) {
            PutNode(R,(*R).N[depth],(*R).P[depth],depth);
            (*R).Nmodified[depth]= FALSE;
          }
          if (depth == (*par).height) {
            (*R).N[depth]= (refnode)malloc((*R).datanodelen);
          }
          else {
            (*R).N[depth]= (refnode)malloc((*R).dirnodelen);
          }
          (*R).P[depth]= (*(*R).N[depth-1]).DIR.entries[i].ptrtosub;
          GetNode(R,(*R).N[depth],(*R).P[depth],depth);
        }
        else if ((*(*R).N[depth-1]).DIR.entries[i].ptrtosub != (*R).P[depth]) {
          NewNode(R,depth);
        }
        if ( i == instind) {
          (*R).NInst[depth]= (*R).N[depth];
        }
        found= FoundRect(R,depth,rectangle,i==instind,infoadr);
        depth--;
      }
    } while (! found && i != (*(*R).N[depth]).DIR.nofentries - 1);

    c = &(*R).count;
      
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
  }
  else {
  
    n= &(*(*R).N[depth]).DATA;
    
    while (! found && i != (*n).nofentries - 1) {
      i++;
      found= RSTEqual(R,(*n).entries[i].rect,rectangle);
      if (found) {
        (*R).E[depth]= i;
        *infoadr= &(*n).entries[i].info;
      }
    }
    
    c = &(*R).count;
    if ((*c).countflag) {
      (*c).datavisitcount++;
    }
  }
  if (found) {
    (*R).EInst[depth]= -1;
    depth++;
    if ((*R).NInst[depth] != NULL) {
      if ((*R).NInst[depth] != (*R).N[depth]) {
        free((*R).NInst[depth]);
      }
      (*R).NInst[depth]= NULL;
    }
    depth--;
  }
  return found;
}

/************************************************************************/

void XstsRgn(RSTREE R,
             int depth,
             typrect rectangle1,
             typrect rectangle2,
             DirQueryProc DirQuery,
             DataQueryProc DataQuery,
             boolean *found)

{
  refcount c;
  refDIRnode DIN;
  refDATAnode DAN;
  boolean istoread;
  int i;
  
  if (depth != (*R).parameters._.height) {
    
    DIN= &(*(*R).N[depth]).DIR;
    
    i= -1;
    do {
      i++;
      if (DirQuery(R,(*DIN).entries[i].rect,rectangle1,rectangle2)) {
        (*R).E[depth]= i;
        istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
        if (istoread) {
          NewNode(R,depth+1);
        }
        XstsRgn(R,depth+1,rectangle1,rectangle2,DirQuery,DataQuery,found);
      }
    } while (! *found && i != (*DIN).nofentries - 1);

    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
  }
  else {
  
    DAN= &(*(*R).N[depth]).DATA;
    
    i= -1;
    while (! *found && i != (*DAN).nofentries - 1) {
      i++;
      if (DataQuery(R,(*DAN).entries[i].rect,rectangle1,rectangle2)) {
        (*R).E[depth]= i;
        *found= TRUE;
      }
    }
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datavisitcount++;
    }
  }
}

/************************************************************************/

void RgnCnt(RSTREE R,
            int depth,
            typrect rectangle1,
            typrect rectangle2,
            DirQueryProc DirQuery,
            DataQueryProc DataQuery,
            int *keysqualifying)

{
  refcount c;
  refDIRnode DIN;
  refDATAnode DAN;
  boolean istoread;
  int i;
  
  if (depth != (*R).parameters._.height) {
  
    DIN= &(*(*R).N[depth]).DIR;
    
    for (i= 0; i < (*DIN).nofentries; i++) {
      if (DirQuery(R,(*DIN).entries[i].rect,rectangle1,rectangle2)) {
        (*R).E[depth]= i;
        istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
        if (istoread) {
          NewNode(R,depth+1);
        }
        RgnCnt(R,depth+1,rectangle1,rectangle2,
               DirQuery, DataQuery,keysqualifying);
      }
    }
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
  }
  else {
  
    DAN= &(*(*R).N[depth]).DATA;
    
    for (i= 0; i < (*DAN).nofentries; i++) {
      if (DataQuery(R,(*DAN).entries[i].rect,rectangle1,rectangle2)) {
        (*R).E[depth]= i;
        (*keysqualifying)++;
      }
    }
    
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datavisitcount++;
    }
      
  }
}

/************************************************************************/

void RgnQuery(RSTREE R,
              int depth,
              typrect rectangle1,
              typrect rectangle2,
              DirQueryProc DirQuery,
              DataQueryProc DataQuery,
              QueryManageProc Manage,
              void *buf,
              boolean *finish)

{
  refcount c;
  refDIRnode DIN;
  refDATAnode DAN;
  boolean istoread;
  typrect rectfound;
  int i;
  
  if (depth != (*R).parameters._.height) {
  
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
    
    DIN= &(*(*R).N[depth]).DIR;
    
    for (i= 0; i < (*DIN).nofentries; i++) {
      if (*finish) {return;}
      
      if (DirQuery(R,(*DIN).entries[i].rect,rectangle1,rectangle2)) {
        (*R).E[depth]= i;
        istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
        if (istoread) {
          NewNode(R,depth+1);
        }
        RgnQuery(R,depth+1,rectangle1,rectangle2,
                 DirQuery,DataQuery,
                 Manage,buf,finish);
      }
    }
  }
  else {
  
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datavisitcount++;
    }
    
    DAN= &(*(*R).N[depth]).DATA;
    
    for (i= 0; i < (*DAN).nofentries; i++) {
      if (*finish) {return;}

      if (DataQuery(R,(*DAN).entries[i].rect,rectangle1,rectangle2)) {
        (*R).E[depth]= i;
        
        CopyRect(R,(*DAN).entries[i].rect,rectfound); /* avoid modification */
        Manage(R,rectfound,
               &(*DAN).entries[i].info,
               buf,
               &(*R).Nmodified[depth],
               finish);
        
      }
    }
  }
}

/************************************************************************/

void All(RSTREE R,
         int depth,
         QueryManageProc Manage,
         void *buf,
         boolean *finish)

{
  refDIRnode DIN;
  refDATAnode DAN;
  typrect rectfound;
  int i;
  
  if (depth != (*R).parameters._.height) {
  
    DIN= &(*(*R).N[depth]).DIR;
    
    for (i= 0; i < (*DIN).nofentries; i++) {
      if (*finish) {return;}
      
      (*R).E[depth]= i;
      if ((*DIN).entries[i].ptrtosub != (*R).P[depth+1]) {
        NewNode(R,depth+1);
      }
      All(R,depth+1,Manage,buf,finish);
    }
  }
  else {
  
    DAN= &(*(*R).N[depth]).DATA;
    
    for (i= 0; i < (*DAN).nofentries; i++) {
      if (*finish) {return;}

      (*R).E[depth]= i;
        
      CopyRect(R,(*DAN).entries[i].rect,rectfound); /* avoid modification */
      Manage(R,rectfound,
             &(*DAN).entries[i].info,
             buf,
             &(*R).Nmodified[depth],
             finish);
        
    }
  }
}

/************************************************************************/

static void dirinfo_init (typdirinfo * info)
{
  memset (info, 0, sizeof (typdirinfo));
  info->Hmin = 1e30;
  info->Hmax = - 1e30;
}

static void dirinfo_add_dir (typdirinfo * info, typdirinfo * a)
{
  info->m01 += a->m01;
  info->m02 += a->m02;
  info->m03 += a->m03;
  info->m11 += a->m11;
  info->m13 += a->m13;
  info->m22 += a->m22;
  info->m23 += a->m23;  
  info->m33 += a->m33;
  info->H0 += a->H0;
  info->H1 += a->H1;
  info->H2 += a->H2;
  info->H3 += a->H3;
  info->H4 += a->H4;
  info->n += a->n;
  if (a->Hmin < info->Hmin) info->Hmin = a->Hmin;
  if (a->Hmax > info->Hmax) info->Hmax = a->Hmax;
}

static void dirinfo_add_data (typdirinfo * info, typrect rect, refinfo data)
{
  double p[3];
  p[0] = rect[0].l; p[1] = rect[1].l; p[2] = data->height;
  info->m01 += p[0];
  info->m02 += p[1];
  info->m03 += p[0]*p[1];
  info->m11 += p[0]*p[0];
  info->m13 += p[0]*p[0]*p[1];
  info->m22 += p[1]*p[1];
  info->m23 += p[0]*p[1]*p[1];
  info->m33 += p[0]*p[0]*p[1]*p[1];
  info->H0 += p[2];
  info->H1 += p[0]*p[2];
  info->H2 += p[1]*p[2];
  info->H3 += p[0]*p[1]*p[2];
  info->H4 += p[2]*p[2];
  info->n++;
  if (p[2] < info->Hmin) info->Hmin = p[2];
  if (p[2] > info->Hmax) info->Hmax = p[2];
}

void UpdateAll(RSTREE R,
	       int depth,
	       typdirinfo * info)

{
  refDIRnode DIN;
  refDATAnode DAN;
  typrect rectfound;
  int i;
  
  if (depth != (*R).parameters._.height) {
  
    DIN= &(*(*R).N[depth]).DIR;
    
    for (i= 0; i < (*DIN).nofentries; i++) {
      (*R).E[depth]= i;
      if ((*DIN).entries[i].ptrtosub != (*R).P[depth+1]) {
        NewNode(R,depth+1);
      }
      dirinfo_init (&(*DIN).entries[i].info);
      (*R).Nmodified[depth] = 1;
      UpdateAll(R,depth+1,&(*DIN).entries[i].info);
      dirinfo_add_dir (info, &(*DIN).entries[i].info);
    }
  }
  else {
  
    DAN= &(*(*R).N[depth]).DATA;
    
    for (i= 0; i < (*DAN).nofentries; i++) {
      (*R).E[depth]= i;

      CopyRect(R,(*DAN).entries[i].rect,rectfound); /* avoid modification */
      dirinfo_add_data (info, rectfound, &(*DAN).entries[i].info);
    }
  }
}

/************************************************************************/

void RgnQueryInfo(RSTREE R,
		  int depth,
		  Check includes,
		  Check intersects,
		  void * data,
		  typdirinfo * info)

{
  refcount c;
  refDIRnode DIN;
  refDATAnode DAN;
  boolean istoread;
  typrect rectfound;
  int i;
  
  if (depth != (*R).parameters._.height) {
  
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).dirvisitcount++;
    }
    
    DIN= &(*(*R).N[depth]).DIR;
    
    for (i= 0; i < (*DIN).nofentries; i++) {
      if ((* includes) ((*DIN).entries[i].rect,data))
	dirinfo_add_dir (info, &(*DIN).entries[i].info);
	   else if ((* intersects) ((*DIN).entries[i].rect,data)) {
        (*R).E[depth]= i;
        istoread= (*DIN).entries[i].ptrtosub != (*R).P[depth+1];
        if (istoread) {
          NewNode(R,depth+1);
        }
        RgnQueryInfo(R,depth+1,includes,intersects,data,info);
      }
    }
  }
  else {
  
    c= &(*R).count;
    if ((*c).countflag) {
      (*c).datavisitcount++;
    }
    
    DAN= &(*(*R).N[depth]).DATA;
    
    for (i= 0; i < (*DAN).nofentries; i++) {

      if ((* includes) ((*DAN).entries[i].rect,data)) {
        (*R).E[depth]= i;
	        
        CopyRect(R,(*DAN).entries[i].rect,rectfound); /* avoid modification */
	dirinfo_add_data (info, rectfound, &(*DAN).entries[i].info);
      }
    }
  }
}

/************************************************************************/
