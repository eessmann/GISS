/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <stdlib.h>
#include <ctype.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <math.h>
#include <sys/times.h>
#include "config.h"
#include "solid.h"
#include "simulation.h"
#include "cartesian.h"

/**
 * @c: a character.
 * @s: a string.
 *
 * Returns: %TRUE if @c belongs to @s, %FALSE otherwise.
 */
gboolean gfs_char_in_string (char c, const char * s)
{
  if (s == NULL)
    return FALSE;
  while (*s != '\0')
    if (*(s++) == c)
      return TRUE;
  return FALSE;
}

/**
 * gfs_file_statement:
 * @fp: a #GtsFile.
 *
 * Reads the next brackets-delimited ({...}) statemement in @fp,
 * including all comments.
 *
 * Returns: a newly allocated string containing the text of the next
 * statement in @fp, or %NULL if an error occured in which case
 * @fp->error is set.
 */
gchar * gfs_file_statement (GtsFile * fp)
{
  g_return_val_if_fail (fp != NULL, NULL);

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return NULL;
  }
  GString * s = g_string_new ("");
  gchar empty[] = "", * comments = fp->comments;
  fp->comments = empty;
  guint scope = fp->scope_max;
  gint c = gts_file_getc (fp);
  while (c != EOF && fp->scope > scope) {
    g_string_append_c (s, c);
    c = gts_file_getc (fp);
  }
  fp->comments = comments;
  if (fp->scope != scope) {
    gts_file_error (fp, "parse error");
    g_string_free (s, TRUE);
    return NULL;
  }
  gchar * ret = s->str;
  g_string_free (s, FALSE);
  return ret;
}

typedef gdouble (* GfsFunctionFunc) (const FttCell * cell,
				     const FttCellFace * face,
				     GfsSimulation * sim);
typedef gdouble (* GfsFunctionDerivedFunc) (const FttCell * cell,
					    const FttCellFace * face,
					    GfsSimulation * sim,
					    gpointer data);

static GfsDerivedVariable * lookup_derived_variable (const gchar * name,
						     GSList * i)
{
  while (i) {
    GfsDerivedVariable * v = i->data;
    if (!strcmp (v->name, name))
      return v;
    i = i->next;
  }
  return NULL;
}

/* GfsGlobal: Object */

struct _GfsGlobal {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  gchar * s;
  guint line;
};

static void global_destroy (GtsObject * object)
{
  g_free (GFS_GLOBAL (object)->s);
  (* gfs_global_class ()->parent_class->destroy) (object);
}

static void global_write (GtsObject * object, FILE * fp)
{
  fprintf (fp, "%s {", object->klass->info.name);
  fputs (GFS_GLOBAL (object)->s, fp);
  fputs ("}\n", fp);
}

static void global_read (GtsObject ** object, GtsFile * fp)
{
  GfsGlobal * global = GFS_GLOBAL (*object);
  GtsObjectClass * klass;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsGlobalClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_global_class ())) {
    gts_file_error (fp, "`%s' is not a GfsGlobal", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  global->line = fp->line;
  g_free (global->s);
  if ((global->s = gfs_file_statement (fp)))
    gts_file_next_token (fp);
}

static void gfs_global_class_init (GtsObjectClass * klass)
{
  klass->destroy = global_destroy;
  klass->read =    global_read;
  klass->write =   global_write;
}

GtsObjectClass * gfs_global_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_global_info = {
      "GfsGlobal",
      sizeof (GfsGlobal),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_global_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_global_info);
  }

  return klass;
}

/* GfsFunction: Object */

struct _GfsFunction {
  GtsObject parent;
  GString * expr;
  gboolean isexpr;
  GModule * module;
  GfsFunctionFunc f;
  gchar * sname;
  GtsSurface * s;
  GfsCartesianGrid * g;
  guint index[4];
  GfsVariable * v;
  GfsDerivedVariable * dv;
  gdouble val;
  gboolean spatial, constant;
  GtsFile fpd;
  gdouble units;
};

static GtsSurface * read_surface (gchar * name, GtsFile * fp)
{
  FILE * fptr = fopen (name, "r");
  GtsSurface * s;
  GtsFile * fp1;

  if (fptr == NULL) {
    gts_file_error (fp, "cannot open file `%s'", name);
    return NULL;
  }
  fp1 = gts_file_new (fptr);
  s = gts_surface_new (gts_surface_class (), gts_face_class (), 
		       gts_edge_class (), gts_vertex_class ());
  if (gts_surface_read (s, fp1)) {
    gts_file_error (fp, "%s:%d:%d: %s", name, fp1->line, fp1->pos, fp1->error);
    gts_object_destroy (GTS_OBJECT (s));
    s = NULL;
  }
  gts_file_destroy (fp1);
  fclose (fptr);
  return s;
}

static GfsCartesianGrid * read_cartesian_grid (gchar * name, GtsFile * fp)
{
  FILE * fptr = fopen (name, "r");
  GtsFile * fp1;
  GfsCartesianGrid * grid;
  GtsObjectClass * klass;

  if (fptr == NULL) {
    gts_file_error (fp, "cannot open file `%s'", name);
    return NULL;
  }

  fp1 = gts_file_new (fptr);

  klass = gfs_cartesian_grid_class ();

  grid = gfs_cartesian_grid_new (klass);
  GtsObject * o = GTS_OBJECT (grid);
  (* klass->read) (&o, fp1);

  if (fp1->type == GTS_ERROR) {
    gts_file_error (fp, "%s:%d:%d: %s", name, fp1->line, fp1->pos, fp1->error);
    gts_object_destroy (GTS_OBJECT(grid));
    grid = NULL;
  }
  gts_file_destroy (fp1);
  fclose (fptr);
  return grid;
}

static gboolean fit_index_dimension (GfsCartesianGrid * grid, guint * val, GtsFile * fp)
{
  guint i, j;
  gchar liste[] = {'x','y','z','t'};

  if (grid->N > 4) {
    gts_file_error (fp, "Cartesian grids can only use four dimensions or less");
    return FALSE;
  }

  for(i = 0; i < grid->N; i++) {
    for (j = 0; j < 4 && *grid->name[i] != liste[j]; j++);
    if (j == 4) {
      gts_file_error (fp, "unknown Cartesian grid index `%s'", grid->name[i]);
      return FALSE;
    }
    val[i] = j;
  }
  return TRUE;
}

static gdouble interpolated_cgd (GfsFunction * f, FttVector * p)
{
  gdouble vecteur[4];
  gdouble val;
  guint i;

  gfs_simulation_map_inverse (gfs_object_simulation (f), p);
  for (i = 0; i < f->g->N; i++)
    switch (f->index[i]) {
    case 0: vecteur[i] = p->x; break;
    case 1: vecteur[i] = p->y; break;
    case 2: vecteur[i] = p->z; break;
    case 3: vecteur[i] = gfs_object_simulation (f)->time.t; break;
    default: g_assert_not_reached ();
    }

  if (!gfs_cartesian_grid_interpolate (f->g, vecteur, &val))
    return 0.;
  return val;
}

static gboolean load_module (GfsFunction * f, GtsFile * fp, gchar * mname)
{
  gchar * path;
  
  path = g_module_build_path (GFS_MODULES_DIR, mname);
  f->module = g_module_open (path, 0);
  g_free (path);
  if (f->module == NULL)
    f->module = g_module_open (mname, 0);
  if (f->module == NULL) {
    gts_file_error (fp, "cannot load module: %s", g_module_error ());
    return FALSE;
  }
  if (!g_module_symbol (f->module, "f", (gpointer) &f->f)) {
    gts_file_error (fp, "module `%s' does not export function `f'", mname);
    g_module_close (f->module);
    return FALSE;
  }
  if (f->constant) {
    f->val = (* f->f) (NULL, NULL, NULL);
    f->f = NULL;
    g_module_close (f->module);
    f->module = NULL;
    if (f->expr) g_string_free (f->expr, TRUE);
    f->expr = NULL;
  }
  return TRUE;
}

/**
 * gfs_function_expression:
 * @fp: a #GtsFile.
 * @is_expression: a pointer to a boolean or %NULL.
 *
 * Reads the expression (in which case @is_expression is set to %TRUE)
 * or function from @fp.
 *
 * Returns: a newly allocated GString containing the result or %NULL
 * in case of error.
 */
GString * gfs_function_expression (GtsFile * fp, gboolean * is_expression)
{
  GString * expr = NULL;

  g_return_val_if_fail (fp != NULL, NULL);

  if (is_expression)
    *is_expression = TRUE;
  if (fp->type == '{') {
    gchar * s = gfs_file_statement (fp);
    if (fp->type == GTS_ERROR)
      return NULL;
    expr = g_string_new ("{");
    g_string_append (expr, s);
    g_free (s);
    g_string_append_c (expr, '}');
    if (is_expression)
      *is_expression = FALSE;
    return expr;
  }
  else {
    static gchar spaces[] = " \t\f\r";
    static gchar operators[] = "+-*/%<>=&^|?:!";
    gint c, scope = 0;
    gchar * s;
    gchar empty[] = "", * comments = fp->comments;

    fp->comments = empty;
    expr = g_string_new (fp->token->str);
    s = expr->str;
    while (*s != '\0') {
      if (*s == '(') scope++;
      else if (*s == ')') scope--;
      s++;
    }
    if (fp->next_token != '\0')
      c = fp->next_token;
    else {
      if (fp->type != '(')
	c = ' ';
      else
	c = gts_file_getc (fp);
    }
    if (strlen (expr->str) == 1 && gfs_char_in_string (expr->str[0], operators))
      while (c != EOF && gfs_char_in_string (c, spaces)) {
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
      }
    while (c != EOF) {
      if (gfs_char_in_string (c, "{}\n")) {
	fp->next_token = c;
	break;
      }
      else if (scope > 0) {
	while (c != EOF && scope > 0) {
	  if (c == '(') scope++;
	  else if (c == ')') scope--;
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
      }
      else if (gfs_char_in_string (c, spaces)) {
	while (c != EOF && gfs_char_in_string (c, spaces)) {
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
	if (c == '(') {
	  scope++;
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
	else {
	  if (!gfs_char_in_string (c, operators)) {
	    fp->next_token = c;
	    break;
	  }
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	  while (c != EOF && gfs_char_in_string (c, spaces)) {
	    g_string_append_c (expr, c);
	    c = gts_file_getc (fp);
	  }
	}
      }
      else if (gfs_char_in_string (c, operators)) {
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
	while (c != EOF && gfs_char_in_string (c, spaces)) {
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
      }
      else {
	if (c == '(') scope++;
	else if (c == ')') scope--;
	if (scope < 0) {
	  fp->next_token = c;
	  break;
	}
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
      }
    }
    g_strchomp (expr->str);
    fp->comments = comments;
    return expr;
  }
}

static gint compile (GtsFile * fp, GfsFunction * f, const gchar * finname)
{
  gchar foutname[] = "/tmp/gfsXXXXXX";
  gchar ferrname[] = "/tmp/gfsXXXXXX";
  gchar ftmpname[] = "/tmp/gfsXXXXXX";
  gint foutd, ferrd, ftmpd;
  gchar * cc;
  gint status;
  gchar cccommand[] = "gcc `pkg-config "
#if FTT_2D
    "gerris2D"
#elif FTT_2D3
    "gerris2D3"
#else /* 3D */
    "gerris3D"
#endif
    " --cflags --libs` -O -Wall -Wno-unused -Werror "
    MODULES_FLAGS;
  
  foutd = mkstemp (foutname);
  ferrd = mkstemp (ferrname);
  ftmpd = mkstemp (ftmpname);
  if (foutd < 0 || ferrd < 0 || ftmpd < 0) {
    gts_file_error (fp, "cannot create temporary file");
    return SIGABRT;
  }
  cc = g_strjoin (" ",
		  cccommand, ftmpname, 
		  "-o", foutname,
                  "`sed 's/@/#/g' <", finname,
		  "| awk '{"
		  "   if ($1 == \"#\" && $2 == \"link\") {"
		  "     for (i = 3; i <= NF; i++) printf (\"%s \", $i);"
		  "     print \"\" > \"/dev/stderr\";"
		  "   }"
		  "   else if ($1 == \"#link\") {"
		  "     for (i = 2; i <= NF; i++) printf (\"%s \", $i);"
		  "     print \"\" > \"/dev/stderr\";"
		  "   } else print $0 > \"/dev/stderr\";"
		  "}' 2>", ftmpname, "` 2>",
		  ferrname, NULL);
  status = system (cc);
  g_free (cc);
  close (ftmpd);
  remove (ftmpname);
  if (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || WTERMSIG (status) == SIGQUIT))
    status = SIGQUIT;
  else if (status == -1 || WEXITSTATUS (status) != 0) {
    GString * msg = g_string_new ("");
    FILE * ferr = fdopen (ferrd, "r");
    gchar * needle;
    gint c;

    while ((c = fgetc (ferr)) != EOF)
      g_string_append_c (msg, c);
    fclose (ferr);
    while ((needle = strstr (msg->str, "GfsFunction:")))
      g_string_erase (msg, needle - msg->str, strlen ("GfsFunction:"));
    gts_file_error (fp, "error compiling expression\n%s", msg->str);
    g_string_free (msg, TRUE);
    status = SIGABRT;
  }
  else {
    if (load_module (f, fp, foutname))
      status = SIGCONT;
    else
      status = SIGABRT;
  }
  close (foutd);
  remove (foutname);
  close (ferrd);
  remove (ferrname);
  return status;
}

static gchar * find_identifier (const gchar * s, const gchar * i)
{
  gchar * f = strstr (s, i);
  static gchar allowed[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_1234567890";

  while (f) {
    if (gfs_char_in_string (f[strlen(i)], allowed) ||
	(f > s && gfs_char_in_string (f[-1], allowed)))
      f = strstr (++f, i);
    else
      return f;
  }
  return NULL;
}

static void function_compile (GfsFunction * f, GtsFile * fp)
{
  if (!HAVE_PKG_CONFIG) {
    gts_file_error (fp, "expecting a number, variable or GTS surface (val)\n"
		    "(functions are not supported on this system)");
    return;
  }
  else {
    GfsSimulation * sim = gfs_object_simulation (f);
    GfsDomain * domain = GFS_DOMAIN (sim);
    gchar finname[] = "/tmp/gfsXXXXXX";
    gint find, status;
    FILE * fin;
    GSList * lv = NULL, * ldv = NULL, * i;

    find = mkstemp (finname);
    if (find < 0) {
      gts_file_error (fp, "cannot create temporary file");
      return;
    }
    fin = fdopen (find, "w");
    fputs ("#include <stdlib.h>\n"
	   "#include <stdio.h>\n"
	   "#include <math.h>\n"
	   "#include <gfs.h>\n",
	   fin);
    if (f->spatial)
      fputs ("#include <gerris/spatial.h>\n", fin);
    else if (!f->constant)
      fputs ("#include <gerris/function.h>\n", fin);
    i = sim->globals;
    while (i) {
      fprintf (fin, "#line %d \"GfsGlobal\"\n", GFS_GLOBAL (i->data)->line);
      fputs (GFS_GLOBAL (i->data)->s, fin);
      fputc ('\n', fin);
      i = i->next;
    }
    if (f->spatial)
      fputs ("double f (double x, double y, double z, double t) {\n"
	     "  _x = x; _y = y; _z = z;\n", 
	     fin);
    else if (f->constant)
      fputs ("double f (void) {\n", fin);
    else {
      fputs ("typedef double (* Func) (const FttCell * cell,\n"
	     "                         const FttCellFace * face,\n"
	     "                         GfsSimulation * sim,\n"
	     "                         gpointer data);\n"
	     "double f (FttCell * cell, FttCellFace * face, GfsSimulation * sim) {\n"
	     "  _sim = sim; _cell = cell;\n",
	     fin);
      i = domain->variables;
      while (i) {
	if (GFS_VARIABLE1 (i->data)->name && 
	    find_identifier (f->expr->str, GFS_VARIABLE1 (i->data)->name))
	  lv = g_slist_prepend (lv, i->data);
	i = i->next;
      }
      i = domain->derived_variables;
      while (i) {
	GfsDerivedVariable * v = i->data;
	if (find_identifier (f->expr->str, v->name))
	  ldv = g_slist_prepend (ldv, v);
	i = i->next;
      }
      if (lv || ldv) {
	GSList * i = lv;

	while (i) {
	  GfsVariable * v = i->data;
	  fprintf (fin, "  double %s;\n", v->name);
	  i = i->next;
	}
	i = ldv;
	while (i) {
	  GfsDerivedVariable * v = i->data;
	  fprintf (fin, "  double %s;\n", v->name);
	  i = i->next;
	}
	if (lv) {
	  fputs ("  if (cell) {\n", fin);
	  i = lv;
	  while (i) {
	    GfsVariable * v = i->data;
	    fprintf (fin, 
		     "    %s = gfs_dimensional_value (GFS_VARIABLE1 (%p),\n"
		     "           GFS_VALUE (cell, GFS_VARIABLE1 (%p)));\n", 
		     v->name, v, v);
	    i = i->next;
	  }
	  fputs ("  } else {\n", fin);
	  i = lv;
	  while (i) {
	    GfsVariable * v = i->data;
	    fprintf (fin, 
		     "    %s = gfs_dimensional_value (GFS_VARIABLE1 (%p),\n"
		     "           gfs_face_interpolated_value (face, GFS_VARIABLE1 (%p)->i));\n", 
		     v->name, v, v);
	    i = i->next;
	  }
	  fputs ("  }\n", fin);
	  g_slist_free (lv);
	}
	if (ldv) {
	  i = ldv;
	  while (i) {
	    GfsDerivedVariable * v = i->data;
	    fprintf (fin, "  %s = (* (Func) %p) (cell, face, sim, ((GfsDerivedVariable *) %p)->data);\n", 
		     v->name, v->func, v);
	    i = i->next;
	  }
	  g_slist_free (ldv);
	}
      }
    }
    fprintf (fin, "#line %d \"GfsFunction\"\n", fp->line);

    if (f->isexpr)
      fprintf (fin, "return %s;\n}\n", f->expr->str);
    else {
      gchar * s = f->expr->str;
      guint len = strlen (s);
      g_assert (s[0] == '{' && s[len-1] == '}');
      s[len-1] = '\0';
      fprintf (fin, "%s\n}\n", &s[1]);
      s[len-1] = '}';
    }
    fclose (fin);
    close (find);

    status = compile (fp, f, finname);
    remove (finname);
    switch (status) {
    case SIGQUIT: exit (0);
    case SIGABRT: return;
    }
  }  
}

#define DEFERRED_COMPILATION ((GfsFunctionFunc) 0x1)

static void check_for_deferred_compilation (GfsFunction * f)
{
  if (f->f == DEFERRED_COMPILATION) {
    function_compile (f, &f->fpd);
    if (f->fpd.type == GTS_ERROR) {
      g_log (G_LOG_DOMAIN, G_LOG_LEVEL_CRITICAL, 
	     "error in deferred compilation\n%s", 
	     f->fpd.error);
      exit (1);
    }
  }
}

static void function_read (GtsObject ** o, GtsFile * fp)
{
  GfsFunction * f = GFS_FUNCTION (*o);
  GfsSimulation * sim;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  sim = gfs_object_simulation (*o);
  domain = GFS_DOMAIN (sim);
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT && fp->type != GTS_STRING &&
      fp->type != '(' && fp->type != '{') {
    gts_file_error (fp, "expecting an expression (val)");
    return;
  }

  if ((f->expr = gfs_function_expression (fp, &f->isexpr)) == NULL)
    return;

  if (f->isexpr) {
    if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
      if (!strcmp (fp->token->str, f->expr->str)) {
	f->val = atof (fp->token->str);
	gts_file_next_token (fp);
	return;
      }
    }
    else if (fp->type == GTS_STRING && !f->spatial && !f->constant) {
      if (strlen (f->expr->str) > 3 &&
	  !strcmp (&(f->expr->str[strlen (f->expr->str) - 4]), ".gts")) {
	if ((f->s = read_surface (f->expr->str, fp)) == NULL)
	  return;
	f->sname = g_strdup (f->expr->str);
	gts_file_next_token (fp);
	return;
      }
      else if (strlen (f->expr->str) > 3 &&
	       !strcmp (&(f->expr->str[strlen (f->expr->str) - 4]), ".cgd")) {
	if ((f->g = read_cartesian_grid (f->expr->str, fp)) == NULL)
	  return;
	if (!fit_index_dimension (f->g, f->index, fp))
	  return;
	f->sname = g_strdup (f->expr->str);
	gts_file_next_token (fp);
	return;
      }
      else if ((f->v = gfs_variable_from_name (domain->variables, f->expr->str))) {
	gts_file_next_token (fp);
	return;
      }
      else if ((f->dv = lookup_derived_variable (f->expr->str, domain->derived_variables))) {
	gts_file_next_token (fp);
	return;
      }
    }
  }

  if (sim->deferred_compilation) {
    f->f = DEFERRED_COMPILATION;
    f->fpd = *fp;
  }
  else
    function_compile (f, fp);

  if (fp->type == GTS_ERROR)
    return;
  gts_file_next_token (fp);
}

static void function_write (GtsObject * o, FILE * fp)
{
  GfsFunction * f = GFS_FUNCTION (o);

  if (GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->write) (o, fp);

  if (f->expr)
    fprintf (fp, " %s", f->expr->str);
  else if (f->module)
    fprintf (fp, " %s", g_module_name (f->module));
  else if (f->v)
    fprintf (fp, " %s", f->v->name);
  else if (f->s || f->g)
    fprintf (fp, " %s", f->sname);
  else
    fprintf (fp, " %g", f->val);
}

static void function_destroy (GtsObject * object)
{
  GfsFunction * f = GFS_FUNCTION (object);

  if (f->module) g_module_close (f->module);
  if (f->expr) g_string_free (f->expr, TRUE);
  if (f->s) {
    gts_object_destroy (GTS_OBJECT (f->s));
    g_free (f->sname);
  }
  if (f->g) {
    gts_object_destroy (GTS_OBJECT (f->g));
    g_free (f->sname);
  }

  (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->destroy) 
    (object);
}

static void gfs_function_class_init (GfsFunctionClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = function_read;
  GTS_OBJECT_CLASS (klass)->write = function_write;
  GTS_OBJECT_CLASS (klass)->destroy = function_destroy;
}

GfsFunctionClass * gfs_function_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunction",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) gfs_function_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_function_info);
  }

  return klass;
}

/**
 * gfs_function_new:
 * @klass: a #GfsFunctionClass.
 * @val: a value.
 *
 * Returns: a new #GfsFunction with constant value @val.
 */
GfsFunction * gfs_function_new (GfsFunctionClass * klass, 
				gdouble val)
{
  GfsFunction * object;

  object = GFS_FUNCTION (gts_object_new (GTS_OBJECT_CLASS (klass)));
  object->val = val;

  return object;
}

GfsFunction * gfs_function_new_from_variable (GfsFunctionClass * klass, 
					      GfsVariable * v)
{
  GfsFunction * object;

  g_return_val_if_fail (v != NULL, NULL);

  object = GFS_FUNCTION (gts_object_new (GTS_OBJECT_CLASS (klass)));
  object->v = v;

  return object;
}

/**
 * gfs_function_set_units:
 * @f: a #GfsFunction.
 * @units: the units of @f.
 *
 * Sets the units of @f.
 */
void gfs_function_set_units (GfsFunction * f, 
			     gdouble units)
{
  g_return_if_fail (f != NULL);
  f->units = units;
}

static gdouble interpolated_value (GfsFunction * f, FttVector * p)
{
  GtsPoint q;
  GtsFace * t;

  gfs_simulation_map_inverse (gfs_object_simulation (f), p);
  q.x = p->x; q.y = p->y;
  t = gts_point_locate (&q, f->s, NULL);
  if (t == NULL)
    return 0.;
  gts_triangle_interpolate_height (GTS_TRIANGLE (t), &q);
  return q.z;
}

/**
 * gfs_function_description:
 * @f: a #GfsFunction.
 * @truncate: whether to truncate long descriptions.
 *
 * Returns: a newly allocated string describing @f.
 */
gchar * gfs_function_description (GfsFunction * f,
				  gboolean truncate)
{
  gchar * s;

  g_return_val_if_fail (f != NULL, NULL);

  if (f->s)
    s = g_strdup (f->sname);
  else if (f->v)
    s = g_strdup (f->v->name);
  else if (f->expr) {
    s = g_strdup (f->expr->str);    
    if (truncate) {
      gchar * c = s;
      guint n = 0;
      
      while (*c != '\0' && !isspace (*c))
	c++;
      while (*c != '\0' && n < 3) {
	*c = '.';
	c++; n++;
      }
      *c = '\0';
    }
  }
  else
    s = g_strdup_printf ("%g", f->val);
  return s;
}

static gdouble adimensional_value (GfsFunction * f, gdouble v)
{
  gdouble L;
  if (v == G_MAXDOUBLE || f->units == 0. || 
      (L = gfs_object_simulation (f)->physical_params.L) == 1.)
    return v;
  return v*pow (L, - f->units);
}

/**
 * gfs_function_value:
 * @f: a #GfsFunction.
 * @cell: a #FttCell or %NULL.
 *
 * Returns: the value of function @f in @cell.
 */
gdouble gfs_function_value (GfsFunction * f, FttCell * cell)
{
  g_return_val_if_fail (f != NULL, 0.);

  gdouble dimensional;
  if (f->s) {
    FttVector p;
    gfs_cell_cm (cell, &p);
    dimensional = interpolated_value (f, &p);
  }
  else if (f->g) {
    FttVector p;
    gfs_cell_cm (cell, &p);
    dimensional = interpolated_cgd (f, &p);
  }
  else if (f->v)
    dimensional = gfs_dimensional_value (f->v, GFS_VALUE (cell, f->v));
  else if (f->dv)
    dimensional = (* (GfsFunctionDerivedFunc) f->dv->func) (cell, NULL,
							    gfs_object_simulation (f),
							    f->dv->data);
  else if (f->f) {
    check_for_deferred_compilation (f);
    dimensional = (* f->f) (cell, NULL, gfs_object_simulation (f));
  }
  else
    dimensional = f->val;
  return adimensional_value (f, dimensional);
}

/**
 * gfs_function_face_value:
 * @f: a #GfsFunction.
 * @fa: a #FttCellFace.
 *
 * Returns: the value of function @f at the center of face @fa.
 */
gdouble gfs_function_face_value (GfsFunction * f, FttCellFace * fa)
{
  g_return_val_if_fail (f != NULL, 0.);
  g_return_val_if_fail (fa != NULL, 0.);

  gdouble dimensional;
  if (f->s) {
    FttVector p;
    ftt_face_pos (fa, &p);
    dimensional = interpolated_value (f, &p);
  }
  else if (f->g) {
    FttVector p;
    ftt_face_pos (fa, &p);
    dimensional = interpolated_cgd (f, &p);
  }
  else if (f->v)
    dimensional = gfs_dimensional_value (f->v, gfs_face_interpolated_value (fa, f->v->i));
  else if (f->dv)
    dimensional = (* (GfsFunctionDerivedFunc) f->dv->func) (NULL, fa,
							    gfs_object_simulation (f),
							    f->dv->data);
  else if (f->f) {
    check_for_deferred_compilation (f);
    dimensional = (* f->f) (NULL, fa, gfs_object_simulation (f));
  }
  else
    dimensional = f->val;
  return adimensional_value (f, dimensional);
}

/**
 * gfs_function_set_constant_value:
 * @f: a #GfsFunction.
 * @val: the value.
 *
 * Sets the value of the constant function @f to @val.
 */
void gfs_function_set_constant_value (GfsFunction * f, gdouble val)
{
  g_return_if_fail (f != NULL);
  g_return_if_fail (!f->f && !f->s && !f->v && !f->dv);

  f->val = val;
}

/**
 * gfs_function_get_constant_value:
 * @f: a #GfsFunction.
 *
 * Returns: the value of function @f if @f is constant, G_MAXDOUBLE
 * otherwise.
 */
gdouble gfs_function_get_constant_value (GfsFunction * f)
{
  g_return_val_if_fail (f != NULL, G_MAXDOUBLE);

  check_for_deferred_compilation (f);
  if (f->f || f->s || f->v || f->dv)
    return G_MAXDOUBLE;
  else
    return adimensional_value (f, f->val);
}

/**
 * gfs_function_get_variable:
 * @f: a #GfsFunction.
 *
 * Returns: the variable containing the value of @f if @f is a simple
 * variable, NULL otherwise.
 */
GfsVariable * gfs_function_get_variable (GfsFunction * f)
{
  g_return_val_if_fail (f != NULL, NULL);

  return f->v;
}

/**
 * gfs_function_read:
 * @f: a #GfsFunction.
 * @domain: a #GfsDomain.
 * @fp: a #GtsFile.
 *
 * Calls the read() method of @f.
 */
void gfs_function_read (GfsFunction * f, gpointer domain, GtsFile * fp)
{
  GtsObject * o = (GtsObject *) f;

  g_return_if_fail (f != NULL);
  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  GTS_OBJECT (f)->reserved = domain;
  (* GTS_OBJECT (f)->klass->read) (&o, fp);
}

/**
 * gfs_function_write:
 * @f: a #GfsFunction.
 * @fp: a file pointer.
 *
 * Calls the write() method of @f.
 */
void gfs_function_write (GfsFunction * f, FILE * fp)
{
  g_return_if_fail (f != NULL);
  g_return_if_fail (fp != NULL);

  (* GTS_OBJECT (f)->klass->write) (GTS_OBJECT (f), fp);
}

/* GfsFunctionSpatial: object */

static void gfs_function_spatial_init (GfsFunction * f)
{
  f->spatial = TRUE;
}

GfsFunctionClass * gfs_function_spatial_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunctionSpatial",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_function_spatial_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_function_class ()),
				  &gfs_function_info);
  }

  return klass;
}

typedef gdouble (* GfsFunctionSpatialFunc) (double x, double y, double z, double t);

/**
 * gfs_function_spatial_value:
 * @f: a #GfsFunction.
 * @p: a #FttVector.
 *
 * Returns: the value of function @f at location @p.
 */
gdouble gfs_function_spatial_value (GfsFunction * f, FttVector * p)
{
  g_return_val_if_fail (f != NULL, 0.);
  g_return_val_if_fail (GFS_IS_FUNCTION_SPATIAL (f), 0.);
  g_return_val_if_fail (p != NULL, 0.);

  if (f->f) {
    GfsSimulation * sim = gfs_object_simulation (f);
    FttVector q = *p;
    check_for_deferred_compilation (f);
    gfs_simulation_map_inverse (sim, &q);
    return (* (GfsFunctionSpatialFunc) f->f) (q.x, q.y, q.z, sim->time.t);
  }
  else
    return f->val;
}

/* GfsFunctionConstant: object */

static void gfs_function_constant_init (GfsFunction * f)
{
  f->constant = TRUE;
}

GfsFunctionClass * gfs_function_constant_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunctionConstant",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_function_constant_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_function_class ()),
				  &gfs_function_info);
  }

  return klass;
}

/**
 * gfs_read_constant:
 * @fp: a #GtsFile.
 * @domain: a #GfsDomain.
 *
 * Reads a constant value from @fp.
 *
 * Returns: the value of the constant or G_MAXDOUBLE if an error
 * occured.
 */
gdouble gfs_read_constant (GtsFile * fp, gpointer domain)
{
  g_return_val_if_fail (fp != NULL, G_MAXDOUBLE);
  g_return_val_if_fail (domain != NULL, G_MAXDOUBLE);

  GfsFunction * f = gfs_function_new (gfs_function_constant_class (), 0.);
  gfs_function_read (f, domain, fp);
  if (fp->type == GTS_ERROR)
    return G_MAXDOUBLE;
  gdouble val = gfs_function_get_constant_value (f);
  gts_object_destroy (GTS_OBJECT (f));
  if (val == G_MAXDOUBLE)
    gts_file_error (fp, "expecting a constant");
  return val;
}

/**
 * gfs_object_class_from_name:
 * @name: the name of the class.
 *
 * Looks for a class called @name. If none is found append the "Gfs"
 * prefix and look again.
 *
 * Returns: the class or %NULL if none was found.
 */
GtsObjectClass * gfs_object_class_from_name (const gchar * name)
{
  GtsObjectClass * klass;

  g_return_val_if_fail (name != NULL, NULL);

  if ((klass = gts_object_class_from_name (name)))
    return klass;
  /* for backward parameter file compatibility */
  if (!strcmp (name, "GtsSurfaceFile"))
    return GTS_OBJECT_CLASS (gfs_solid_class ());
  gchar * ename = g_strconcat ("Gfs", name, NULL);
  klass = gts_object_class_from_name (ename);
  g_free (ename);
  return klass;
}

static void eigsrt (gdouble d[FTT_DIMENSION], gdouble v[FTT_DIMENSION][FTT_DIMENSION])
{
  gint k, j, i;
  gdouble p;

  for (i = 0; i < FTT_DIMENSION - 1; i++) {
    p = d[k = i];

    for (j = i + 1; j < FTT_DIMENSION; j++)
      if (d[j] >= p) 
	p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < FTT_DIMENSION; j++) {
	p = v[j][i];
	v[j][i] = v[j][k];
	v[j][k] = p;
      }
    }
  }
}

#define ROTATE(a,i,j,k,l) {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}

/**
 * gfs_eigenvalues:
 * @a: a symmetric matrix.
 * @d: a vector.
 * @v: another matrix.
 *
 * Fills @d (resp. @v) with the eigenvalues (resp. eigenvectors) of
 * matrix @a.
 */
void gfs_eigenvalues (gdouble a[FTT_DIMENSION][FTT_DIMENSION],
		      gdouble d[FTT_DIMENSION],
		      gdouble v[FTT_DIMENSION][FTT_DIMENSION])
{
  gint j, iq, ip, i;
  gdouble tresh, theta, tau, t, sm, s, h, g, c, b[FTT_DIMENSION], z[FTT_DIMENSION];

  for (ip = 0; ip < FTT_DIMENSION; ip++) {
    for (iq = 0; iq < FTT_DIMENSION; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < FTT_DIMENSION; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < FTT_DIMENSION - 1; ip++) {
      for (iq = ip + 1; iq < FTT_DIMENSION; iq++)
	sm += fabs (a[ip][iq]);
    }
    if (sm == 0.0) {
      eigsrt (d, v);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(FTT_DIMENSION*FTT_DIMENSION);
    else
      tresh = 0.0;
    for (ip = 0; ip < FTT_DIMENSION - 1; ip++) {
      for (iq = ip + 1; iq < FTT_DIMENSION; iq++) {
	g = 100.0*fabs (a[ip][iq]);
	if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) && fabs(d[iq]) + g == fabs(d[iq]))
	  a[ip][iq] = 0.0;
	else if (fabs (a[ip][iq]) > tresh) {
	  h = d[iq] - d[ip];
	  if (fabs(h) + g == fabs(h))
	    t = a[ip][iq]/h;
	  else {
	    theta = 0.5*h/a[ip][iq];
	    t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt (1 + t*t);
	  s = t*c;
	  tau = s/(1.0 + c);
	  h = t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.0;
	  for (j = 0; j <= ip - 1; j++)
	    ROTATE (a, j, ip, j, iq);
	  for (j = ip + 1; j <= iq - 1; j++)
	    ROTATE (a, ip, j, j, iq);
	  for (j = iq + 1; j < FTT_DIMENSION; j++)
	    ROTATE(a, ip, j, iq, j);
	  for (j = 0; j < FTT_DIMENSION; j++)
	    ROTATE(v, j, ip, j, iq);
	}
      }
    }
    for (ip = 0; ip < FTT_DIMENSION; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  /* Too many iterations */
  for (i = 0; i < FTT_DIMENSION; i++) {
    for (j = 0; j < FTT_DIMENSION; j++)
      fprintf (stderr, "%10.3g ", a[i][j]);
    fprintf (stderr, "\n");
  }
  g_assert_not_reached ();
}

/**
 * gfs_matrix_inverse:
 * @m: a square matrix.
 * @n: size of the matrix.
 * @pivmin: the minimum value of the pivoting coefficient.
 *
 * Replaces @m with its inverse.
 *
 * Returns: 0. if the inversion encounters a pivot coefficient smaller
 * than or equal to @pivmin (i.e. @m is non-invertible), the minimum
 * absolute value of the pivoting coefficient otherwise.
 */
gdouble gfs_matrix_inverse (gdouble ** m, guint n, gdouble pivmin)
{
  gint * indxc, * indxr, * ipiv;
  gint i, icol = 0, irow = 0, j, k, l, ll;
  gdouble big, dum, pivinv, temp, minpiv = G_MAXDOUBLE;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (pivmin >= 0., 0.);

  indxc = g_malloc (sizeof (gint)*n);
  indxr = g_malloc (sizeof (gint)*n);
  ipiv = g_malloc (sizeof (gint)*n);
  
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
	for (k = 0; k < n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs (m[j][k]) >= big) {
	      big = fabs (m[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++) 
	SWAP (m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin) {
      g_free (indxc);
      g_free (indxr);
      g_free (ipiv);
      return 0.;
    }
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
	dum = m[ll][icol];
	m[ll][icol] = 0.0;
	for (l = 0; l < n; l++)
	  m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
	SWAP (m[k][indxr[l]], m[k][indxc[l]]);
  }
  g_free (indxc);
  g_free (indxr);
  g_free (ipiv);
  return minpiv;
}

/**
 * gfs_matrix_new:
 * @n: the size of the matrix.
 * @p: the size of the matrix.
 * @size: the size of the matrix elements.
 *
 * The matrix elements are initialised to zero.
 *
 * Returns: a newly allocated matrix.
 */
gpointer gfs_matrix_new (guint n, guint p, guint size)
{
  guint i;
  gpointer * m, a;
  
  g_return_val_if_fail (n > 0, NULL);
  g_return_val_if_fail (p > 0, NULL);
  g_return_val_if_fail (size > 0, NULL);

  m = g_malloc (n*sizeof (gpointer));
  a = g_malloc0 (n*p*size);
  for (i = 0; i < n; i++)
    m[i] = GUINT_TO_POINTER (GPOINTER_TO_UINT (a) + i*p*size);
  return m;
}

/**
 * gfs_matrix_free:
 * @m: a matrix allocated with gfs_matrix_new().
 *
 * Frees the memory occupied by @m.
 */
void gfs_matrix_free (gpointer m)
{
  g_return_if_fail (m != NULL);

  g_free (((gpointer *) m)[0]);
  g_free (m);
}

/**
 * gfs_clock_new:
 *
 * Returns: a new #GfsClock.
 */
GfsClock * gfs_clock_new (void)
{
  GfsClock * t = g_malloc (sizeof (GfsClock));

  t->start = -1;
  t->started = FALSE;
  return t;
}

/**
 * gfs_clock_start:
 * @t: a #GfsClock.
 *
 * Starts clock @t.
 */
void gfs_clock_start (GfsClock * t)
{
  struct tms tm;

  g_return_if_fail (t != NULL);
  g_return_if_fail (!t->started);

  if (times (&tm) < 0)
    g_warning ("cannot read clock");
  t->start = tm.tms_utime;
  t->started = TRUE;
}

/**
 * gfs_clock_stop:
 * @t: a #GfsClock.
 *
 * Stops clock @t.
 */
void gfs_clock_stop (GfsClock * t)
{
  struct tms tm;

  g_return_if_fail (t != NULL);
  g_return_if_fail (t->started);

  if (times (&tm) < 0)
    g_warning ("cannot read clock");
  t->stop = tm.tms_utime;
  t->started = FALSE;
}

/**
 * gfs_clock_elapsed:
 * @t: a #GfsClock.
 *
 * Returns: the time elapsed in seconds since @t was started.
 */
gdouble gfs_clock_elapsed (GfsClock * t)
{
  g_return_val_if_fail (t != NULL, 0.);
  g_return_val_if_fail (t->start >= 0, 0.);

  if (t->started == FALSE)
    return (t->stop - t->start)/(gdouble) sysconf (_SC_CLK_TCK);
  else {
    struct tms tm;
    if (times (&tm) < 0)
      g_warning ("cannot read clock");
    return (tm.tms_utime - t->start)/(gdouble) sysconf (_SC_CLK_TCK);
  }
}

/**
 * gfs_clock_destroy:
 * @t: a #GfsClock.
 *
 * Destroys the clock, freeing the memory allocated for it.
 */
void gfs_clock_destroy (GfsClock * t)
{
  g_return_if_fail (t != NULL);

  g_free (t);
}
