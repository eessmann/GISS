#include "init.h"
#define WIKI "http://gfs.sf.net/wiki/index.php/"

static void key_value_pair (const char * key)
{
  printf ("'%s' : '" WIKI "%s',\\\n", key, key);
}

int main (int argc, char * argv[])
{
  GtsObjectClass ** klass;

  gfs_init (&argc, &argv);
  klass = gfs_classes ();
  printf ("klass = {\\\n");
  key_value_pair ("Define");
  while (*klass) {
    key_value_pair ((*klass)->info.name);
    klass++;
  }
  printf ("}\n");
  return 0;
}
