#include "init.h"

int main (int argc, char * argv[])
{
  GtsObjectClass ** klass = gfs_classes ();
  printf ("klass = {\\\n");
  while (*klass) {
    printf ("'%s' : 'http://gfs.sf.net/wiki/index.php/%s',\\\n", (*klass)->info.name, (*klass)->info.name);
    klass++;
  }
  printf ("}\n");
  return 0;
}
