#include <string.h>
#include "init.h"
#define WIKI "http\\://gfs.sf.net/wiki/index.php/"

static void key_value_pair (const char * key, FILE * lang)
{
  fprintf (lang, "gfs_keyword = \"%s\"\n", key);
  /* keywords must start with Gfs */
  g_assert (strstr (key, "Gfs") == key);
  fprintf (lang, "gfs_keyword = \"%s\"\n", &(key[3]));
}

int main (int argc, char * argv[])
{
  GtsObjectClass ** klass;

  klass = gfs_classes ();

  printf ("# Language file for source-highlight\n"
	  "# Generated automatically by classes.c\n"
	  "\n");

  key_value_pair ("GfsDefine", stdout);
  key_value_pair ("GfsProjectionParams", stdout);
  key_value_pair ("GfsApproxProjectionParams", stdout);
  key_value_pair ("GfsPhysicalParams", stdout);
  key_value_pair ("GfsAdvectionParams", stdout);

  /* Map module  */
  key_value_pair ("GfsMapProjection", stdout);

  while (*klass) {
    key_value_pair ((*klass)->info.name, stdout);
    klass++;
  }
  
  printf ("\n"
	  "include \"cpp.lang\"\n"
	  "\n"
	  "comment start \"#\"\n"
	  "\n"
	  "redef preproc = \"C preprocessor command is not compatible with"
	  " the use of # as comment character in GTS\"\n");
  return 0;
}
