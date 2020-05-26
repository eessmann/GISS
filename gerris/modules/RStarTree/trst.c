/* -----  trst.c  ----- */


#include "RStarTree.h"
#include "Set.h"

/* constants */

#define STDMODE 0644
#define RSTsuffix ".RSF"


/* types */

typedef
  struct {
         typatomkey  center[NumbOfDim], ext[NumbOfDim];
         } data;


/* declarations */

    int GetLOF(File f);
   void getname(char *name);
   void DataAndTreeName(void);
   void QueryName(void);
   void DisplErr(char *name);
boolean opendata(void);
   void closedata(void);
boolean openquery(void);
   void closequery(void);
   void exeCreate(void);
   void exeRemove(void);
   void exeOpen(RSTREE *RST);
   void exeClose(void);
   void GetTestParameters(File file,
                          char *name);
   void exeInsert(void);
   void exeDelete(void);
   void exeSearch(void);
   void exeExistsRegion(void);
   void exeRegionCount(void);
   void exeRegionQueries(void);
   void exeJoinCountNv(void);
   void exeJoinNv(void);
   void exeJoinCountX(void);
   void exeJoinX(void);
   void exeInquire(void);
   void exeComputeNodeLengths(void);
   void ShortInstruction(void);
boolean Equal(RSTREE rst,
              typrect RSTrect,
              typrect qrect1,
              typrect qrect2);
boolean Intersects(RSTREE rst,
                   typrect RSTrect,
                   typrect qrect1,
                   typrect qrect2);
boolean Encloses(RSTREE rst,
                 typrect RSTrect,
                 typrect qrect1,
                 typrect qrect2);
boolean IsContained(RSTREE rst,
                    typrect RSTrect,
                    typrect qrect1,
                    typrect qrect2);
boolean AlwaysTrue(RSTREE rst,
                   typrect RSTrect,
                   typrect qrect1,
                   typrect qrect2);
void    CompCount0();
void ManageQuery(RSTREE R,
                 typrect rectangle,
                 refinfo infoptr,
                 void *buf,
                 boolean *modify,
                 boolean *finish);
void ManageJoin(RSTREE R, RSTREE R2,
                typrect rectangle, typrect rectangle2,
                refinfo infoptr, refinfo infoptr2,
                void *buf, void *buf2,
                boolean *finish);


/* global variables */

char distname[80];
char quername[80];
char rfilename[80];
char ch, dummy;

File datafile, querfile;

int ferr, nbytes, pos;
int begin, end;
int constante;
int height;
typrect rectangle;
data datarect;
boolean notdmsg, donemsg;
RSTREE R, R2;

int compcount;
int globnumfound;
int globpaircount;

int i;


/***********************************************************************/

void main(int argc, char *argv[])

{
  char mainchoose, choose;
  boolean WorkMenu;
  NoRSTree(&R);
  NoRSTree(&R2);
  ShortInstruction();
  do {
    printf("\n");
    printf("%20s,   %20s,\n","(forced)Create = \"c\"","Remove = \"R\"");
    printf("%20s,   %20s,\n","Open = \"o\"","Close = \".\"");
    printf("%30s,\n","Inquire Description = \"I\"");
    printf("%30s,\n","Node lenghts computation = \"`\"");
    printf("%30s,\n","leave this menu = \"-\"");
    printf("%30s.\n","quit = \"q\"");
    printf("Input: ");
    do {
      mainchoose= getchar();
    } while ((mainchoose != 'c') && (mainchoose != 'R') &&
             (mainchoose != 'o') && (mainchoose != '.') &&
             (mainchoose != 'I') &&
             (mainchoose != '`') &&
             (mainchoose != '-') &&
             (mainchoose != 'q'));
    dummy= getchar();
    WorkMenu= mainchoose != '`';
    if (mainchoose != 'q') {
      if (mainchoose == '-') {
      }
      else if (mainchoose == 'c') {
        if (R != NULL) {
          exeClose();
        }
        DataAndTreeName();
        exeCreate();
        exeOpen(&R);
      }
      else if (mainchoose == 'R') {
        if (R != NULL) {
          exeClose();
        }
        exeRemove();
      }
      else if (mainchoose == 'o') {
        if (R != NULL) {
          exeClose();
        }
        DataAndTreeName();
        exeOpen(&R);
        printf("\n     ===========================================\n\n");
        printf("Prepare for Join? (y/n) ");
        do {
          ch= getchar();
        } while ((ch != 'y') && (ch != 'n'));
        dummy= getchar();
        if (ch == 'y') {
          printf("R2=R? (y/n) ");
          do {
            ch= getchar();
          } while ((ch != 'y') && (ch != 'n'));
          dummy= getchar();
          if (ch == 'y') {
            R2= R;
          }
          else {
            DataAndTreeName();
            exeOpen(&R2);
          }
        }
      }
      else if (mainchoose == '.') {
        exeClose();
      }
      else if (mainchoose == 'I') {
        exeInquire();
      }
      else if (mainchoose == '`') {
        exeComputeNodeLengths();
      }
      if (WorkMenu) {
        printf("\n");
        printf("%20s,   %20s,\n","Insert = \"i\"","Delete = \"d\"");
        printf("%30s,\n","Search = \"s\"");
        printf("%30s,\n","ExistsRegion = \"e\"");
        printf("%20s,   %20s,\n","RegionCount = \"r\"","RegionQueries = \"[\"");
        printf("%20s,   %20s,\n","JoinCountNv = \"j\"","JoinNv = \"=\"");
        printf("%20s,   %20s,\n","JoinCountX = \"J\"","JoinX = \"#\"");
        printf("%30s.\n","leave this menu = \"-\"");
        printf("Input: ");
        do {
          choose= getchar();
        } while ((choose != 'i') && (choose != 'd') &&
                 (choose != 's') &&
                 (choose != 'e') &&
                 (choose != 'r') && (choose != '[') &&
                 (choose != 'j') && (choose != '=') &&
                 (choose != 'J') && (choose != '#') &&
                 (choose != '-'));
        dummy= getchar();
        if (choose == '-') {
        }
        else if (choose == 'i') {
          if (opendata()) {
            exeInsert();
            closedata();
          }
        }
        else if (choose == 'd') {
          if (opendata()) {
            exeDelete();
            closedata();
          }
        }
        else if (choose == 's') {
          if (opendata()) {
            exeSearch();
            closedata();
          }
        }
        else if (choose == 'e') {
          QueryName();
          if (openquery()) {
            exeExistsRegion();
            closequery();
          }
        }
        else if (choose == 'r') {
          QueryName();
          if (openquery()) {
            exeRegionCount();
            closequery();
          }
        }
        else if (choose == '[') {
          QueryName();
          if (openquery()) {
            exeRegionQueries();
            closequery();
          }
        }
        else if (choose == 'j') {
          exeJoinCountNv();
        }
        else if (choose == '=') {
          exeJoinNv();
        }
        else if (choose == 'J') {
          exeJoinCountX();
        }
        else if (choose == '#') {
          exeJoinX();
        }
      }
    }
    else {
      exeClose();
    }
  } while (mainchoose != 'q');
}

/***********************************************************************/

void ShortInstruction()

{
  printf("\n  --- TestRStarTree ---\n");
  printf("This testprogram principally works on one \"main\" tree only.\n");
  printf("Exception: JoinCount.\n");
  printf("Joins work on the main tree and a secondary tree.\n");
  printf("Joins are normally initialized by \"Open\".\n");
  printf("If an error occurs the corresponding procedure does nothing\n");
  printf("but issuing a message, though actions may be retried.\n");
  printf("\"Create\", \"Open\" and \"Remove\" automatically perform\n");
  printf("\"Close\" on the main tree (if it is not NULL)\n");
  printf("before performing the concerning action.\n");
  printf("\"quit\" will also try to close the main tree.\n");
  printf("The secondary tree is automatically closed after the join.\n");
}

/***********************************************************************/

int GetLOF(File f) /* GetLengthOfFile */

{
  struct stat *refstatus;
  
  refstatus= (struct stat *)malloc(sizeof(struct stat));
  ferr= fstat(f,refstatus);
  if (ferr == -1) {
    return 0;
  }
  else {
    return (*refstatus).st_size;
  }
}

/***********************************************************************/

void getname(char *name)

{
  printf("\nData File (without suffix): ");
  scanf("%s",name);
}

/***********************************************************************/

void DataAndTreeName()

{
  getname(distname);
  strcpy(rfilename,distname);
  strcat(rfilename,RSTsuffix);
}

/***********************************************************************/

void QueryName()

{
  getname(quername);
}

/***********************************************************************/

void DisplErr(char *name)

{
  char message[80];
  
  printf("\n--------------------\n");
  printf("%s %d\n","FileError:",errno);
  strcpy(message,"Error handling file:\n");
  strcat(message,name);
  perror(message);
  printf("--------------------\n");
  errno= 0;
}
  
/***********************************************************************/

boolean opendata()

{
  printf("Open Data File:\n");
  printf("%s\n",distname);
  datafile= open(distname,O_RDONLY,STDMODE);
  if (datafile == -1) {
    DisplErr(distname);
    return FALSE;
  }
  else {
    GetTestParameters(datafile,distname);
    return TRUE;
  }
}

/***********************************************************************/

void closedata()

{
  printf("Close Data File:\n");
  printf("%s\n",distname);
  ferr= close(datafile);
  if (ferr == -1) {
    DisplErr(distname);
  }
}

/***********************************************************************/

boolean openquery()

{
  printf("Open Query File:\n");
  printf("%s\n",quername);
  querfile= open(quername,O_RDONLY,STDMODE);
  if (querfile == -1) {
    DisplErr(quername);
    return FALSE;
  }
  else {
    GetTestParameters(querfile,quername);
    return TRUE;
  }
}

/***********************************************************************/

void closequery()

{
  printf("Close Query File:\n");
  printf("%s\n",quername);
  ferr= close(querfile);
  if (ferr == -1) {
    DisplErr(quername);
  }
}

/***********************************************************************/

void exeCreate()

{
  int pagelen;
  boolean unique;
  
  printf("Length of a page: ");
  scanf("%d",&pagelen);
  printf("Unique?, (y/n): ");
  do {
    ch= getchar();
  } while ((ch != 'y') && (ch != 'n'));
  dummy= getchar();
  unique= ch == 'y';
  printf("%s%s%s","\nRemoveRST(",rfilename,"):\n");
  if (RemoveRST(rfilename)) {
    printf("Done\n");
  }
  printf("%s%s,%d,%d%s","CreateRST(",rfilename,pagelen,unique,"):\n");
  if (CreateRST(rfilename,pagelen,unique)) {
    printf("Done\n");
  }
  else {
    printf("FAILURE\n");
  }
}

/***********************************************************************/

void exeRemove()

{
  printf("%s%s%s","\nRemoveRST(",rfilename,"):\n");
  if (RemoveRST(rfilename)) {
    printf("Done\n");
  }
  else {
    printf("FAILURE\n");
  }
}

/***********************************************************************/

void exeOpen(RSTREE *r)

{
  printf("%s%s%s","OpenRST(",rfilename,"):\n");
  if (OpenRST(r,rfilename)) {
    printf("Done\n");
  }
  else {
    printf("FAILURE\n");
  }
}

/***********************************************************************/

void exeClose()

{
  printf("%s%s%s","CloseRST(",rfilename,"):\n");
  if (CloseRST(&R)) {
    printf("Done\n");
  }
  else {
    printf("FAILURE\n");
  }
}

/***********************************************************************/

void GetTestParameters(File file, char *name)

{
  printf("\nExecution from/up to entrynumber\n");
  printf("(\"0 0\": Execution on all entries of the file.)\n");
  printf("Input: ");
  scanf("%d",&begin); scanf("%d",&end);
  if (begin != 0) { begin--; }
  if (end == 0) {
    end= GetLOF(file) / sizeof(data);
  }
  pos= lseek(file,begin*sizeof(data),SEEK_SET);
  if (pos == -1) {
    DisplErr(name);
  }
  printf("Message if done =\"d\",");
  printf(" not done =\"n\", message off =\"o\"\n");
  printf("Input: ");
  do {
    ch= getchar();
  } while ((ch != 'd') && (ch != 'n') && (ch != 'o'));
  dummy= getchar();
  notdmsg= ch == 'n';
  donemsg= ch == 'd';
  printf("Echo resp. performance-measurement(Insert) at ");
  printf("entrynumber MOD how much?\n");
  printf("(\"0\": none, not \"0\": mandatory at the end)\n");
  printf("Input: ");
  scanf("%d",&constante);
  dummy= getchar(); /* scanf leaves trailing white characters */
  if (constante == 0) {constante= MAXINT;}
}

/***********************************************************************/

void exeInsert()

{
  int xdirv, xdatav, xdirr, xdatar,
      xdirm, xdatam, xdirw, xdataw;
  boolean inserted, success;
  int height;
  int d;
  typinfo info;
  
  CountsOn0(R);
  i= begin;
  do {
    nbytes= read(datafile,&datarect,sizeof(data));
    if (nbytes <= 0) {
      DisplErr(distname);
    }
    i++;
    for (d= 0; d < NumbOfDim; d++) {
      rectangle[d].l= datarect.center[d] - datarect.ext[d];
      rectangle[d].h= datarect.center[d] + datarect.ext[d];
    }
    /* Begin set info.contents to record number */
    /**/
    info.contents= i;
    /**/
    /* End   set info.contents to record number */
    success= InsertRecord(R,rectangle,&info,&inserted);
    if (notdmsg) {
      if (! inserted) {
        printf("%s%d\n"," NOT INSERTED: ",i);
      }
    }
    else if (donemsg) {
      if (inserted) {
        printf("%s%d\n"," INSERTED: ",i);
      }
    }
    if ((i % constante == 0) || (i == end)) {
      printf("%5d%s",i,"::  ");
      GetHeight(R,&height);
      printf("%s%3d\n","height:",height);
    }
    if (! success) {printf("FAILURE\n"); return;}
  } while (i != end);
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
  GetCountWrite(R,&xdirm,&xdatam,&xdirw,&xdataw);
  printf("%s%d\n","dirm: ",xdirm);
  printf("%s%d\n","datam: ",xdatam);
  printf("%s%d\n","dirw: ",xdirw);
  printf("%s%d\n","dataw: ",xdataw);
}

/***********************************************************************/

void exeDelete()

{
  int xdirv, xdatav, xdirr, xdatar,
      xdirm, xdatam, xdirw, xdataw;
  boolean found, success;
  int d;
  
  CountsOn0(R);
  i= begin;
  do {
    nbytes= read(datafile,&datarect,sizeof(data));
    if (nbytes <= 0) {
      DisplErr(distname);
    }
    i++;
    for (d= 0; d < NumbOfDim; d++) {
      rectangle[d].l= datarect.center[d] - datarect.ext[d];
      rectangle[d].h= datarect.center[d] + datarect.ext[d];
    }
    success= DeleteRecord(R,rectangle,&found);    
    if (notdmsg) {
      if (! found) {
        printf("%s%d\n"," NOT FOUND: ",i);
      }
    }
    else if (donemsg) {
      if (found) {
        printf("%s%d\n"," FOUND: ",i);
      }
    }
    if ((i % constante == 0) || (i == end)) {
      printf("%5d%s",i,"::  ");
      GetHeight(R,&height);
      printf("%s%3d\n","height:",height);
    }
    if (! success) {printf("FAILURE\n"); return;}
  } while (i != end);
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
  GetCountWrite(R,&xdirm,&xdatam,&xdirw,&xdataw);
  printf("%s%d\n","dirm: ",xdirm);
  printf("%s%d\n","datam: ",xdatam);
  printf("%s%d\n","dirw: ",xdirw);
  printf("%s%d\n","dataw: ",xdataw);
}

/***********************************************************************/

void exeSearch()

{
  int xdirv, xdatav, xdirr, xdatar;
  int allfound;
  boolean found, success;
  typrect unused;
  typinfo check;
  int d;
  
  allfound= 0;
  CountsOn0(R);
  
  i= begin;
  do {
    nbytes= read(datafile,&datarect,sizeof(data));
    if (nbytes <= 0) {
      DisplErr(distname);
    }
    i++;
    for (d= 0; d < NumbOfDim; d++) {
      rectangle[d].l= datarect.center[d] - datarect.ext[d];
      rectangle[d].h= datarect.center[d] + datarect.ext[d];
    }
    success= Find(R,rectangle,&found,&check,sizeof(check));
    if (found) {
      allfound++;
    }
    if (notdmsg) {
      if (! found) {
        printf("%s%d\n"," NOT FOUND: ",i);
      }
    }
    else if (donemsg) {
      if (found) {
        printf("%s%d\n"," FOUND: ",i);
      }
    }
    if (found) {
      /* Begin check correctness of an integer infopart */
      /**/
      if (check.contents != i) {
        printf("%d%s%d\n",i,
                        ": ODD INFOPART: ",
                        check.contents);
      }
      /**/
      /* End   check correctness of an integer infopart */
    }
    if ((i % constante == 0) || (i == end)) {
      printf("%6d%s",i,"::  ");
      GetHeight(R,&height);
      printf("%s%3d","height:",height);
      printf("%s%6d\n","  found:",allfound);
    }
    if (! success) {printf("FAILURE\n"); return;}
  } while (i != end);
  printf("\n");
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
}

/***********************************************************************/

void exeExistsRegion()

{
  int xdirv, xdatav, xdirr, xdatar;
  int allfound;
  boolean found, success;
  typrect unused;
  int d;
  
  allfound= 0;
  CountsOn0(R);
  CompCount0();
  
  i= begin;
  do {
    nbytes= read(querfile,&datarect,sizeof(data));
    if (nbytes <= 0) {
      DisplErr(quername);
    }
    i++;
    for (d= 0; d < NumbOfDim; d++) {
      rectangle[d].l= datarect.center[d] - datarect.ext[d];
      rectangle[d].h= datarect.center[d] + datarect.ext[d];
    }
    success= ExistsRegion(R,rectangle,unused,Encloses,Equal,&found);
    if (found) {
      allfound++;
    }
    if (notdmsg) {
      if (! found) {
        printf("%s%d\n"," NOT FOUND: ",i);
      }
    }
    else if (donemsg) {
      if (found) {
        printf("%s%d\n"," FOUND: ",i);
      }
    }
    if ((i % constante == 0) || (i == end)) {
      printf("%5d%s",i,"::  ");
      GetHeight(R,&height);
      printf("%s%3d","height:",height);
      printf("%s%6d\n","  found:",allfound);
    }
    if (! success) {printf("FAILURE\n"); return;}
  } while (i != end);
  printf("\n");
  printf("%s%d\n","comp: ",compcount);
  printf("\n");
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
}

/***********************************************************************/

void exeRegionCount()

{
  int xdirv, xdatav, xdirr, xdatar;
  int numfound, allfound, emptyquery;
  typrect unused;
  boolean success;
  int d;
  
  allfound= 0;
  emptyquery= 0;
  CountsOn0(R);
  CompCount0();

  i= begin;
  do {
    nbytes= read(querfile,&datarect,sizeof(data));
    if (nbytes <= 0) {
      DisplErr(quername);
    }
    i++;
    for (d= 0; d < NumbOfDim; d++) {
      rectangle[d].l= datarect.center[d] - datarect.ext[d];
      rectangle[d].h= datarect.center[d] + datarect.ext[d];
    }
    success= RegionCount(R,rectangle,unused,Intersects,Intersects,&numfound);
    if (numfound == 0) {
      emptyquery++;
    }
    else {
      allfound= allfound+numfound;
    }
    if (notdmsg) {
      if (numfound == 0) {
        printf("%s%d\n"," NOT FOUND: ",i);
      }
    }
    else if (donemsg) {
      if (numfound != 0) {
        printf("%s%d\n"," FOUND: ",i);
      }
    }
    if ((i % constante == 0) || (i == end)) {
      printf("%5d%s",i,"::  ");
      GetHeight(R,&height);
      printf("%s%3d","height:",height);
      printf("%s%6d\n","  found:",allfound);
    }
    if (! success) {printf("FAILURE\n"); return;}
  } while (i != end);
  printf("\n");
  printf("%s%d\n","comp: ",compcount);
  printf("\n");
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","empty queries: ",emptyquery);
  printf("%s%d\n","rectangles found: ",allfound);
  printf("\n");
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
}

/***********************************************************************/

void exeRegionQueries()

{
  int xdirv, xdatav, xdirr, xdatar;
  int allfound, emptyquery;
  typrect unused;
  boolean success;
  Set buffer;
  int d;
  
  allfound= 0;
  emptyquery= 0;
  CountsOn0(R);
  CompCount0();
  
  buffer= EmptySet();
  
  i= begin;
  do {
    nbytes= read(querfile,&datarect,sizeof(data));
    if (nbytes <= 0) {
      DisplErr(quername);
    }
    i++;
    for (d= 0; d < NumbOfDim; d++) {
      rectangle[d].l= datarect.center[d] - datarect.ext[d];
      rectangle[d].h= datarect.center[d] + datarect.ext[d];
    }
    globnumfound= 0;
    success= RegionQuery(R,
                         rectangle,
                         unused,
                         Intersects,
                         Intersects,
                         ManageQuery,
                         &buffer);
    if (globnumfound == 0) {
      emptyquery++;
    }
    else {
      allfound= allfound+globnumfound;
    }
    if (notdmsg) {
      if (globnumfound == 0) {
        printf("%s%d\n"," NOT FOUND: ",i);
      }
    }
    else if (donemsg) {
      if (globnumfound != 0) {
        printf("%s%d\n"," FOUND: ",i);
      }
    }
    if ((i % constante == 0) || (i == end)) {
      printf("%5d%s",i,"::  ");
      GetHeight(R,&height);
      printf("%s%3d","height:",height);
      printf("%s%6d\n","  found:",allfound);
    }
    if (! success) {printf("FAILURE\n"); return;}
  } while (i != end);
  printf("\n");
  printf("%s%d\n","comp: ",compcount);
  printf("\n");
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","empty queries: ",emptyquery);
  printf("%s%d\n","rectangles found: ",allfound);
  printf("\n");
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
    
  WriteSet(buffer);
  printf("\n");
  WriteTree(buffer);
  printf("\n");
}

/***********************************************************************/

void ManageQuery(RSTREE R,
                 typrect rectangle,
                 refinfo infoptr,
                 void *buf,
                 boolean *modify,
                 boolean *finish)

{
  Set *set;
  char dummy;
  int d;
  
  /***** ----- count ----- *****/
  
  globnumfound++;
  
  /***** ----- prompt the user ----- *****/
  
  printf(">"); dummy= getchar();
  
  /***** ----- print rectangle ----- *****/
  
  for (d= 0; d < NumbOfDim; d++) {
    printf("%25.17f",rectangle[d].l);
    printf("%25.17f",rectangle[d].h);
    printf("\n");
  }
  
  /***** ----- print info part ----- *****/
  
  printf("%15d\n",(*infoptr).contents);
  
  /***** ----- get info part ----- *****/
  
  set= buf;
  *set= IncludeSet(*set,(*infoptr).contents);
  
  /***** ----- modify record ----- *****/
  
  (*infoptr).contents= 42;
  *modify= TRUE;
  
  /***** ----- finish after finding one ----- *****/
  
  *finish= TRUE;
  
}

/***********************************************************************/

void exeJoinCountNv()

{
  int xdirv, xdatav, xdirr, xdatar;
  int paircount;
  typrect unused;
  boolean success;
  
  CountsOn0(R);
  CountsOn0(R2);
  CompCount0();
  
  globpaircount= 0;
  success= JoinCountNv(R,R2,
                       unused,unused,unused,unused,
                       AlwaysTrue,AlwaysTrue,AlwaysTrue,AlwaysTrue,
                       Intersects,Intersects,
                       &paircount);
  if (! success) {printf("FAILURE\n"); return;}
  printf("\n%s%d\n\n","Number of Pairs found: ",paircount);
  printf("\n");
  printf("%s%d\n","comp: ",compcount);
  printf("\n");
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
  GetCountRead(R2,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);

  if (R2 != R) {
    printf("Close R2:\n");
    if (CloseRST(&R2)) {
      printf("Done\n");
    }
    else {
      printf("FAILURE\n");
    }
  }
}

/***********************************************************************/

void exeJoinNv()

{
  int xdirv, xdatav, xdirr, xdatar;
  typrect unused;
  boolean success;
  void *buffer, *buffer2;

  CountsOn0(R);
  CountsOn0(R2);
  CompCount0();
  
  buffer= malloc(sizeof(typinfo));
  buffer2= malloc(sizeof(typinfo));

  success= JoinNv(R,R2,
                  unused,unused,unused,unused,
                  AlwaysTrue,AlwaysTrue,AlwaysTrue,AlwaysTrue,
                  Intersects,Intersects,
                  ManageJoin,buffer,buffer2);
  if (! success) {printf("FAILURE\n"); return;}
  printf("\n%s%d\n\n","Number of Pairs found: ",globpaircount);
  printf("\n");
  printf("%s%d\n","comp: ",compcount);
  printf("\n");
  GetCountRead(R,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);
  GetCountRead(R2,&xdirv,&xdatav,&xdirr,&xdatar);
  printf("%s%d\n","dirv: ",xdirv);
  printf("%s%d\n","datav: ",xdatav);
  printf("%s%d\n","dirr: ",xdirr);
  printf("%s%d\n","datar: ",xdatar);

  if (R2 != R) {
    printf("Close R2:\n");
    if (CloseRST(&R2)) {
      printf("Done\n");
    }
    else {
      printf("FAILURE\n");
    }
  }
}

/***********************************************************************/

void ManageJoin(RSTREE R, RSTREE R2,
                typrect rectangle, typrect rectangle2,
                refinfo infoptr, refinfo infoptr2,
                void *buf, void *buf2,
                boolean *finish)

{
# define MAXPAIRS 42
  
  char dummy;
  int d;
  
  /***** ----- count pairs ----- *****/
  
  globpaircount++;
  
  /***** ----- prompt the user ----- *****/
  
  printf(">"); dummy= getchar();
  
  /***** ----- print rectangle of R ----- *****/
  
  for (d= 0; d < NumbOfDim; d++) {
    printf("%25.17f",rectangle[d].l);
    printf("%25.17f",rectangle[d].h);
    printf("\n");
  }
  
  /***** ----- print info part of R ----- *****/
  
  printf("%15d\n",(*infoptr).contents);
  
  /***** ----- print rectangle of R2 ----- *****/
  
  for (d= 0; d < NumbOfDim; d++) {
    printf("%25.17f",rectangle2[d].l);
    printf("%25.17f",rectangle2[d].h);
    printf("\n");
  }
  
  /***** ----- print info part of R2 ----- *****/
  
  printf("%15d\n",(*infoptr2).contents);
  
  /***** ----- finish after finding MAXPAIRS pairs ----- *****/
  /*
  if (globpaircount == MAXPAIRS) {
    *finish= TRUE;
  }
  */
  
# undef MAXPAIRS
}

/***********************************************************************/

void exeJoinCountX()

{
  printf("\nJoinCountX NOT IMPLEMENTED\n");
}

/***********************************************************************/

void exeJoinX()

{
  printf("\nJoinX NOT IMPLEMENTED\n");
}

/***********************************************************************/

void exeInquire()

{
  char name[80];
  int NumberOfDimensions,
      SIZEdirentry, SIZEdataentry, SIZEinfo,
      dirM, dataM,
      PageLength,
      NumbDirPages, NumbDataPages;
  int PagesPerLevel[50];
  int NumbRecords,
      height;
  boolean unique;
  double spaceutil[50];
  double sumspaceutil;
  boolean success;
  int i;
  
  success= InquireRSTDesc(R,
                          name,
                          &NumberOfDimensions,
                          &PageLength,
                          &SIZEdirentry,
                          &SIZEdataentry,
                          &SIZEinfo,
                          &dirM,
                          &dataM,
                          &NumbDirPages,
                          &NumbDataPages,
                          PagesPerLevel,
                          &NumbRecords,
                          &height,
                          &unique);
  if (! success) {printf("FAILURE\n"); return;}
  printf("%20s%s\n","name: ",name);
  printf("%20s%d\n","NumberOfDimensions: ",NumberOfDimensions);
  printf("%20s%d\n","PageLength: ",PageLength);
  printf("%20s%d\n","SIZEdirentry: ",SIZEdirentry);
  printf("%20s%d\n","SIZEdataentry: ",SIZEdataentry);
  printf("%20s%d\n","SIZEinfo: ",SIZEinfo);
  printf("%20s%d\n","dirM: ",dirM);
  printf("%20s%d\n","dataM: ",dataM);
  printf("%20s%d\n","NumbDirPages: ",NumbDirPages);
  printf("%20s%d\n","NumbDataPages: ",NumbDataPages);
  printf("%20s%d\n","NumbPages: ",NumbDirPages+NumbDataPages);
  printf("%20s%d\n","NumbRecords: ",NumbRecords);
  printf("%20s%d\n","height: ",height);
  printf("%20s%d\n","unique: ",unique);
  printf("pages per level:\n");
  for(i= 0; i < height; i++) {
    printf("%7d",PagesPerLevel[i]);
  }
  printf("\n");
  for (i= 0; i < height-1; i++) {
    spaceutil[i]= (double)PagesPerLevel[i+1] /
                  (double)(PagesPerLevel[i]*dirM);
  }
  spaceutil[height-1]= (double)(NumbRecords) /
                       (double)(PagesPerLevel[height-1]*dataM);
  printf("space utilization:\n");
  for(i= 0; i < height; i++) {
    printf("%.2e ",spaceutil[i]);
  }
  printf("\n");
  sumspaceutil= 0.0;
  for (i= 1; i < height-1; i++) {
    sumspaceutil= sumspaceutil+spaceutil[i];
  }
  printf("%s%.2e\n","    avg spc util dir (without root): ",
  sumspaceutil / (double)(height-2));
  sumspaceutil= sumspaceutil+spaceutil[height-1];
  printf("%s%.2e\n","avg spc util overall (without root): ",
  sumspaceutil / (double)(height-1));
}

/***********************************************************************/

void exeComputeNodeLengths()

{
  int infolen;
  int entryqty, bytesqty;
  int keylen, direntrylen, dataentrylen;
  int dirlen, datalen;
  int dirnumb, datanumb;
  char again;
  do {
    printf("Length of the info part: ");
    scanf("%d",&infolen);
    printf("Number of entries: ");
    scanf("%d",&entryqty);
    printf("Computing from sizeof(data),\n");
    keylen= sizeof(data);
    direntrylen= keylen+sizeof(int);
    dataentrylen= keylen+infolen;
    dirlen= sizeof(int)+entryqty*direntrylen;
    datalen= sizeof(int)+entryqty*dataentrylen;
    printf("directory page length is ");
    printf("%5d\n",dirlen);
    printf("     data page length is ");
    printf("%5d\n",datalen);
    printf("  Number of bytes: ");
    scanf("%d",&bytesqty);
    printf("Computing from sizeof(data),\n");
    dirnumb= (bytesqty-sizeof(int)) / direntrylen;
    datanumb= (bytesqty-sizeof(int)) / dataentrylen;
    printf("number of directory page entries is ");
    printf("%5d\n",dirnumb);
    printf("     number of data page entries is ");
    printf("%5d\n",datanumb);
    printf("\nagain?, (y/n): ");
    do {
      again= getchar();
    } while ((again != 'y') && (again != 'n'));
    dummy= getchar();
  } while (again == 'y');
}

/***********************************************************************/

boolean Equal(RSTREE R,
              typrect RSTrect,
              typrect queryrect,
              typrect unused)

{
  int maxdim= NumbOfDim -1;
  boolean eql;
  int d;
  
  compcount++;
  d= -1;
  do {
    d++;
    eql= RSTrect[d].l == queryrect[d].l &&
         RSTrect[d].h == queryrect[d].h;
  } while (eql && d != maxdim);
  return eql;
}

/***********************************************************************/

boolean Intersects(RSTREE R,
                   typrect RSTrect,
                   typrect queryrect,
                   typrect unused)

{
  int maxdim= NumbOfDim -1;
  boolean inter;
  int d;
  
  compcount++;
  d= -1;
  do {
    d++;
    inter= RSTrect[d].l <= queryrect[d].h &&
           RSTrect[d].h >= queryrect[d].l;
  } while (inter && d != maxdim);
  return inter;
}

/***********************************************************************/

boolean Encloses(RSTREE R,
                 typrect RSTrect,
                 typrect queryrect,
                 typrect unused)

{
  int maxdim= NumbOfDim -1;
  boolean encl;
  int d;
  
  compcount++;
  d= -1;
  do {
    d++;
    encl= RSTrect[d].l <= queryrect[d].l &&
          RSTrect[d].h >= queryrect[d].h;
  } while (encl && d != maxdim);
  return encl;
}

/***********************************************************************/

boolean IsContained(RSTREE R,
                    typrect RSTrect,
                    typrect queryrect,
                    typrect unused)

{
  int maxdim= NumbOfDim -1;
  boolean iscont;
  int d;
  
  compcount++;
  d= -1;
  do {
    d++;
    iscont= RSTrect[d].l >= queryrect[d].l &&
            RSTrect[d].h <= queryrect[d].h;
  } while (iscont && d != maxdim);
  return iscont;
}

/***********************************************************************/

boolean AlwaysTrue(RSTREE R,
                   typrect unused1,
                   typrect unused2,
                   typrect unused3)

{
  return TRUE;
}

/***********************************************************************/

void CompCount0()

{
  compcount= 0;
}

/***********************************************************************/
