/* -----  Set.h  ----- */
#ifndef __Set_h
#define __Set_h


#include <stdio.h>
#include <math.h>
#include <malloc.h>


/* private includes */

#include "SetDef.h"


/* ----- types ----- */

typedef SetRec *Set;


/* ----- procedure declarations ----- */

Set EmptySet(void);

Set IncludeSet(Set set, int value);

Set ExcludeSet(Set set, int value);

int Cardinality(Set set);

void WriteSet(Set set);

void WriteTree(Set set);


#endif /* !__Set_h */
