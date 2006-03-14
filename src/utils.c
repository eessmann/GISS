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

/* GfsFunction: Object */

struct _GfsFunction {
  GtsObject parent;
  GString * expr;
  GModule * module;
  GfsFunctionFunc f;
  gchar * sname;
  GtsSurface * s;
  GfsVariable * v;
  GfsDerivedVariable * dv;
  gdouble val;
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
    gint c, scope;

    expr = g_string_new ("{");
    scope = fp->scope_max;
    c = gts_file_getc (fp);
    while (c != EOF && fp->scope > scope) {
      g_string_append_c (expr, c);
      c = gts_file_getc (fp);
    }
    g_string_append_c (expr, '}');
    if (fp->scope != scope) {
      gts_file_error (fp, "parse error");
      g_string_free (expr, TRUE);
      return NULL;
    }
    if (is_expression)
      *is_expression = FALSE;
    return expr;
  }
  else {
    static gchar spaces[] = " \t\f\r";
    static gchar operators[] = "+-*/%<>=&^|?:";
    gint c, scope = 0;
    gchar * s;

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
    while (c != EOF) {
      if (gfs_char_in_string (c, "{}\n")) {
	fp->next_token = c;
	g_strchomp (expr->str);
	return expr;
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
	if (!gfs_char_in_string (c, operators)) {
	  fp->next_token = c;
	  g_strchomp (expr->str);
	  return expr;
	}
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
	while (c != EOF && gfs_char_in_string (c, spaces)) {
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
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
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
      }
    }
    g_strchomp (expr->str);
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
#if defined(__APPLE__) && defined(__MACH__)
    " --cflags --libs` -O -fPIC -bundle -x c";
#else  /* not MACOSX */
    " --cflags --libs` -O -fPIC -shared -x c";
#endif /* not MACOSX*/
  
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

static void function_read (GtsObject ** o, GtsFile * fp)
{
  GfsFunction * f = GFS_FUNCTION (*o);
  GfsSimulation * sim;
  GfsDomain * domain;
  gboolean isexpr;

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

  if ((f->expr = gfs_function_expression (fp, &isexpr)) == NULL)
    return;

  if (isexpr) {
    if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
      if (!strcmp (fp->token->str, f->expr->str)) {
	f->val = atof (fp->token->str);
	gts_file_next_token (fp);
	return;
      }
    }
    else if (fp->type == GTS_STRING) {
      if (strlen (f->expr->str) > 3 &&
	  !strcmp (&(f->expr->str[strlen (f->expr->str) - 4]), ".gts")) {
	if ((f->s = read_surface (f->expr->str, fp)) == NULL)
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

  if (!HAVE_PKG_CONFIG) {
    gts_file_error (fp, "expecting a number, variable or GTS surface (val)");
    return;
  }
  else {
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
	   "#include <gfs.h>\n"
	   "typedef double (* Func) (const FttCell * cell,\n"
	   "                         const FttCellFace * face,\n"
	   "                         GfsSimulation * sim,\n"
	   "                         gpointer data);\n"
	   "static double Dirichlet = 1.;\n"
	   "static double Neumann = 0.;\n"
	   "double f (FttCell * cell, FttCellFace * face, GfsSimulation * sim) {\n",
	   fin);
    i = domain->variables;
    while (i) {
      if (find_identifier (f->expr->str, GFS_VARIABLE1 (i->data)->name))
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
	fputs ("  g_return_val_if_fail (cell != NULL, 0.);\n", fin);
	i = lv;
	while (i) {
	  GfsVariable * v = i->data;
	  fprintf (fin, "  %s = GFS_VARIABLE (cell, %d);\n", v->name, v->i);
	  i = i->next;
	}
	g_slist_free (lv);
      }
      i = ldv;
      while (i) {
	GfsDerivedVariable * v = i->data;
	fprintf (fin, "  %s = (* (Func) %p) (cell, face, sim, ((GfsDerivedVariable *) %p)->data);\n", 
		 v->name, v->func, v);
	i = i->next;
      }
      g_slist_free (ldv);
    }
    fprintf (fin, "#line %d \"GfsFunction\"\n", fp->line);

    if (isexpr)
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
  else if (f->s)
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

static gdouble interpolated_value (GfsFunction * f, FttVector * p)
{
  GtsPoint q;
  GtsFace * t;

  q.x = p->x; q.y = p->y;
  t = gts_point_locate (&q, f->s, NULL);
  if (t == NULL) {
    g_warning ("%s: cannot locate point (%g,%g)", f->sname, p->x, p->y);
    return 0.;
  }
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

  if (f->s) {
    FttVector p;

    gfs_cell_cm (cell, &p);
    return interpolated_value (f, &p);
  }
  else if (f->v)
    return GFS_VARIABLE (cell, f->v->i);
  else if (f->dv)
    return (* (GfsFunctionDerivedFunc) f->dv->func) (cell, NULL, 
						     gfs_object_simulation (f), 
						     f->dv->data);
  else if (f->f)
    return (* f->f) (cell, NULL, gfs_object_simulation (f));
  else
    return f->val;
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

  if (f->s) {
    FttVector p;

    ftt_face_pos (fa, &p);
    return interpolated_value (f, &p);
  }
  else if (f->v)
    return gfs_face_interpolated_value (fa, f->v->i);
  else if (f->dv)
    return (* (GfsFunctionDerivedFunc) f->dv->func) (fa->cell, fa,
						     gfs_object_simulation (f), 
						     f->dv->data);
  else if (f->f)
    return (* f->f) (fa->cell, fa, gfs_object_simulation (f));
  else
    return f->val;
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

  if (f->f || f->s || f->v || f->dv)
    return G_MAXDOUBLE;
  else
    return f->val;
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

/**
 * gfs_derived_variable_from_name:
 * @i: a list of #GfsDerivedVariable.
 * @name: a name.
 *
 * Returns: the #GfsDerivedVariable @name of @list or %NULL.
 */
GfsDerivedVariable * gfs_derived_variable_from_name (GSList * i, const gchar * name)
{
  g_return_val_if_fail (name != NULL, NULL);

  while (i) {
    GfsDerivedVariable * v = i->data;
    if (!strcmp (v->name, name))
      return v;
    i = i->next;
  }
  return NULL;
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

  klass = gts_object_class_from_name (name);
  if (klass == NULL) {
    gchar * ename = g_strconcat ("Gfs", name, NULL);
    klass = gts_object_class_from_name (ename);
    g_free (ename);
  }
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
 * Returns: %FALSE if the inversion encounters a pivot coefficient
 * smaller than or equal to @pivmin (i.e. @m is non-invertible), %TRUE
 * otherwise.
 */
gboolean gfs_matrix_inverse (gdouble ** m, guint n, gdouble pivmin)
{
  gint * indxc, * indxr, * ipiv;
  gint i, icol = 0, irow = 0, j, k, l, ll;
  gdouble big, dum, pivinv, temp;

  g_return_val_if_fail (m != NULL, FALSE);
  g_return_val_if_fail (pivmin >= 0., FALSE);

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
      return FALSE;
    }
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
  return TRUE;
}

/**
 * gfs_matrix_new:
 * @n: the size of the matrix.
 * @size: the size of the matrix elements.
 *
 * The matrix elements are initialised to zero.
 *
 * Returns: a newly allocated matrix.
 */
gpointer gfs_matrix_new (guint n, guint size)
{
  guint i;
  gpointer * m, a;
  
  g_return_val_if_fail (n > 0, NULL);
  g_return_val_if_fail (size > 0, NULL);

  m = g_malloc (n*sizeof (gpointer));
  a = g_malloc0 (n*n*size);
  for (i = 0; i < n; i++)
    m[i] = GUINT_TO_POINTER (GPOINTER_TO_UINT (a) + i*n*size);
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
