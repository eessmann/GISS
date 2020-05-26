/* ----- Set.c ----- */


#include "Set.h"


/* declarations */

Set AddLeft(Set set, Set branch);
void OutputSet(Set set, int *cardcount);


/*********************************************************************/

Set EmptySet()

{
  return NULL;
}

/*********************************************************************/

Set IncludeSet(Set set, int value)

{
  Set newson;
  
  if (set == NULL) {
    newson= (Set)malloc(sizeof(SetRec));
    (*newson).value= value;
    (*newson).lson= NULL;
    (*newson).rson= NULL;
    return newson;
  }
  else {
    if (value < (*set).value) {
      (*set).lson= IncludeSet((*set).lson,value);
    }
    else if (value > (*set).value) {
      (*set).rson= IncludeSet((*set).rson,value);
    }
  }
  return set;
}

/*********************************************************************/

Set ExcludeSet(Set set, int value)

{
  Set help;
  
  if (set != NULL) {
    if (value == (*set).value) {
      if ((*set).rson == NULL) {
        help= (*set).lson;
        free(set);
        set= help;
      }
      else {
        (*set).rson= AddLeft((*set).rson,(*set).lson);
        help= (*set).rson;
        free(set);
        set= help;
      }
    }
    else {
      if (value < (*set).value) {
        (*set).lson= ExcludeSet((*set).lson,value);
      }
      else {
        (*set).rson= ExcludeSet((*set).rson,value);
      }
    }
  }
  return set;
}

/*********************************************************************/

Set AddLeft(Set set, Set branch)

{
  if ((*set).lson == NULL) {
    (*set).lson= branch;
  }
  else {
    (*set).lson= AddLeft((*set).lson,branch);
  }
  return set;
}

/*********************************************************************/

int Cardinality(Set set)

{
  if (set == NULL) {
    return 0;
  }
  else {
    return 1+Cardinality((*set).lson)+Cardinality((*set).rson);
  }
}

/*********************************************************************/

void WriteSet(Set set)

{
  int cardcount;
  
  cardcount= 0;
  OutputSet(set,&cardcount);
}

/*********************************************************************/

void OutputSet(Set set, int *cardcount)

{
  if (set != NULL) {
    OutputSet((*set).lson,cardcount);
    printf("%5d",(*set).value);
    (*cardcount)++;
    if (*cardcount % 10 == 0)
      printf("\n");
    OutputSet((*set).rson,cardcount);
  }
}

/*********************************************************************/

void WriteTree(Set set)

{
  if (set == 0) {
    printf(".");
  }
  else {
    printf("[");
    WriteTree((*set).lson);
    printf("%d",(*set).value);
    WriteTree((*set).rson);
    printf("]");
  }
}

/*********************************************************************/
