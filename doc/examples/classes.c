#include <string.h>
#include "init.h"
#define WIKI "http\\://gfs.sf.net/wiki/index.php/"

static void key_value_pair (const char * key, FILE * lang)
{
  static int first = 1;
  if (first) {
    fprintf (lang, "gfs_keyword = ");
    first = 0;
  }
  else
    fprintf (lang, ",\n");

  /* keywords must start with Gfs */
  g_assert (strstr (key, "Gfs") == key);
  fprintf (lang, "    '(Gfs){0,1}%s'", &(key[3]));
}

int main (int argc, char * argv[])
{
  GtsObjectClass ** klass;

  gfs_init (&argc, &argv);
  klass = gfs_classes ();

  printf ("# Language file for source-highlight\n"
	  "# Generated automatically by classes.c\n"
	  "\n"
	  "include \"cpp.lang\"\n"
	  "\n"
	  "comment start \"#\"\n"
	  "\n"
	  "redef preproc = \"C preprocessor command is not compatible with"
	  " the use of # as comment character in GTS\"\n"
	  "\n");
  key_value_pair ("GfsDefine", stdout);
  key_value_pair ("GfsProjectionParams", stdout);
  key_value_pair ("GfsApproxProjectionParams", stdout);
  while (*klass) {
    key_value_pair ((*klass)->info.name, stdout);
    klass++;
  }
  return 0;
}
