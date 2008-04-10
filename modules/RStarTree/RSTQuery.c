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
