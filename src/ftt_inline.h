/**
 * ftt_level_size:
 * @level: a guint.
 *
 * Returns: the size of a cell of level @level.
 */
G_INLINE_FUNC
gdouble ftt_level_size (guint level)
{
  gdouble size = 1.;

  while (level) {
    size /= 2.;
    level--;
  }

  return size;
}

/**
 * ftt_cell_size:
 * @cell: a #FttCell.
 *
 * Returns: the size of @cell.
 */
G_INLINE_FUNC
gdouble ftt_cell_size (const FttCell * cell)
{
  g_return_val_if_fail (cell != NULL, 0.);

  return ftt_level_size (ftt_cell_level (cell));
}

/**
 * ftt_cell_volume:
 * @cell: a #FttCell.
 *
 * Returns: the volume (area in 2D) of @cell.
 */
G_INLINE_FUNC
gdouble ftt_cell_volume (const FttCell * cell)
{
  gdouble size;

  g_return_val_if_fail (cell != NULL, 0.);

  size = ftt_level_size (ftt_cell_level (cell));
#if (FTT_2D || FTT_2D3)
  return size*size;
#else  /* FTT_3D */
  return size*size*size;
#endif /* FTT_3D */
}

/**
 * ftt_cell_children:
 * @cell: a #FttCell.
 * @children: a #FttCellChildren.
 *
 * Fills @children with the children of @cell.
 * 
 * This function fails if @cell is a leaf.
 */
G_INLINE_FUNC
void ftt_cell_children (const FttCell * cell,
			FttCellChildren * children)
{
  struct _FttOct * oct;
  guint i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (cell));
  g_return_if_fail (children != NULL);

  oct = cell->children;
  for (i = 0; i < FTT_CELLS; i++)
    children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[i])) ? 
      NULL : &(oct->cell[i]);
}

/**
 * ftt_cell_children_direction:
 * @cell: a #FttCell.
 * @d: a direction.
 * @children: a #FttCellChildren.
 *
 * Fills @children with the children (2 in 2D, 4 in 3D, 2 or 4 in 2D3)
 * of @cell in direction @d.
 * 
 * This function fails if @cell is a leaf.
 *
 * Returns: the number of children in direction @d.
 */
G_INLINE_FUNC
guint ftt_cell_children_direction (const FttCell * cell,
				   FttDirection d,
				   FttCellChildren * children)
{
  struct _FttOct * oct;
  guint i;
#if (FTT_2D || FTT_2D3)
  static gint index[FTT_NEIGHBORS_2D][FTT_CELLS/2] =
  {{1, 3},
   {0, 2},
   {0, 1},
   {2, 3}};
#else  /* FTT_3D */
  static gint index[FTT_NEIGHBORS][FTT_CELLS/2] =
  {{1, 3, 5, 7},
   {0, 2, 4, 6},
   {0, 1, 4, 5},
   {2, 3, 6, 7},
   {0, 1, 2, 3},
   {4, 5, 6, 7}};
#endif /* FTT_3D */

  g_return_val_if_fail (cell != NULL, 0);
  g_return_val_if_fail (!FTT_CELL_IS_LEAF (cell), 0);
  g_return_val_if_fail (d < FTT_NEIGHBORS, 0);
  g_return_val_if_fail (children != NULL, 0);

  oct = cell->children;

#if FTT_2D3
  if (d >= FTT_NEIGHBORS_2D) {
    for (i = 0; i < FTT_CELLS; i++)
      children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[i])) ? NULL : &(oct->cell[i]);
    return FTT_CELLS;
  }
#endif /* 2D3 */

  for (i = 0; i < FTT_CELLS/2; i++)
    children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[index[d][i]])) ? 
      NULL : &(oct->cell[index[d][i]]);
  return FTT_CELLS/2;
}

/**
 * ftt_cell_child_corner:
 * @cell: a #FttCell.
 * @d: a set of perpendicular directions.
 *
 * This function fails if @cell is a leaf.  
 *
 * Returns: the children of @cell in the corner defined by directions @d.
 */
G_INLINE_FUNC
FttCell * ftt_cell_child_corner (const FttCell * cell,
				 FttDirection d[FTT_DIMENSION])
{
#if (FTT_2D || FTT_2D3)
  static gint index[FTT_NEIGHBORS_2D][FTT_NEIGHBORS_2D] = {
    {-1,-1,1,3},
    {-1,-1,0,2},
    {1,0,-1,-1},
    {3,2,-1,-1}
  };
  gint i;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (!FTT_CELL_IS_LEAF (cell), NULL);

  g_return_val_if_fail (d[0] < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (d[1] < FTT_NEIGHBORS, NULL);

#  if FTT_2D3
  if (d[0] >= FTT_NEIGHBORS_2D)
    i = index[d[1]][d[2]];
  else if (d[1] >= FTT_NEIGHBORS_2D)
    i = index[d[0]][d[2]];
  else
#  endif
    i = index[d[0]][d[1]];
#else  /* FTT_3D */
  static gint index[FTT_NEIGHBORS][FTT_NEIGHBORS][FTT_NEIGHBORS] = {
    {{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {-1,-1,-1,-1,1,5},{-1,-1,-1,-1,3,7},
     {-1,-1,1,3,-1,-1},{-1,-1,5,7,-1,-1}},
    {{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {-1,-1,-1,-1,0,4},{-1,-1,-1,-1,2,6},
     {-1,-1,0,2,-1,-1},{-1,-1,4,6,-1,-1}},
    {{-1,-1,-1,-1,1,5},{-1,-1,-1,-1,0,4},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {1,0,-1,-1,-1,-1},{5,4,-1,-1,-1,-1}},
    {{-1,-1,-1,-1,3,7},{-1,-1,-1,-1,2,6},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1},
     {3,2,-1,-1,-1,-1},{7,6,-1,-1,-1,-1}},
    {{-1,-1,1,3,-1,-1},{-1,-1,0,2,-1,-1},
     {1,0,-1,-1,-1,-1},{3,2,-1,-1,-1,-1},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1}},
    {{-1,-1,5,7,-1,-1},{-1,-1,4,6,-1,-1},
     {5,4,-1,-1,-1,-1},{7,6,-1,-1,-1,-1},
     {-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,-1}},
  };
  gint i;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (!FTT_CELL_IS_LEAF (cell), NULL);
  g_return_val_if_fail (d[0] < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (d[1] < FTT_NEIGHBORS, NULL);
  g_return_val_if_fail (d[2] < FTT_NEIGHBORS, NULL);

  i = index[d[0]][d[1]][d[2]];
#endif /* FTT_3D */

  g_return_val_if_fail (i >= 0, NULL);

  return FTT_CELL_IS_DESTROYED (&(cell->children->cell[i])) ? NULL:
    &(cell->children->cell[i]);
}

/**
 * ftt_cell_neighbors_not_cached:
 * @cell: a #FttCell.
 * @neighbors: a #FttCellNeighbors.
 *
 * Fills @neighbors with the neighbors of @cell (does not use saved
 * values even if available).  
 */
G_INLINE_FUNC
void ftt_cell_neighbors_not_cached (const FttCell * cell,
				    FttCellNeighbors * neighbors)
{
  static gint neighbor_index[FTT_NEIGHBORS][FTT_CELLS]
#if FTT_2D
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2}};
#elif FTT_2D3
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2},
       {-1,-2,-3,-4},
       {-1,-2,-3,-4}};
#else  /* FTT_3D */
    = {{1,-1,3,-3,5,-5,7,-7},
       {-2,0,-4,2,-6,4,-8,6},
       {-3,-4,0,1,-7,-8,4,5},
       {2,3,-1,-2,6,7,-5,-6},
       {-5,-6,-7,-8,0,1,2,3},
       {4,5,6,7,-1,-2,-3,-4}};
#endif /* FTT_3D */
  guint n, d;
  struct _FttOct * parent;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (neighbors != NULL);

  if (FTT_CELL_IS_ROOT (cell)) {
    memcpy (neighbors, &((struct _FttRootCell *) cell)->neighbors,
	    sizeof (FttCellNeighbors));
    return;
  }

  parent = cell->parent;
  n = FTT_CELL_ID (cell);
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    gint nn = neighbor_index[d][n];
    FttCell * c;

    if (nn >= 0) /* neighbor belongs to same Oct */
      c = &(parent->cell[nn]);
    else {       /* neighbor belongs to neighboring Cell or Oct */
      c = parent->neighbors.c[d];
      if (c != NULL && c->children != NULL)
	c = &(c->children->cell[- nn - 1]);
    }
    if (c == NULL || FTT_CELL_IS_DESTROYED (c))
      neighbors->c[d] = NULL;
    else
      neighbors->c[d] = c;
  }
}

/**
 * ftt_cell_neighbor_not_cached:
 * @cell: a #FttCell.
 * @d: a direction.
 *
 * Returns: the neighbor of @cell in direction @d or %NULL if @cell
 * has no neighbor in this direction (does not use saved values even
 * if available).  
 */
G_INLINE_FUNC
FttCell * ftt_cell_neighbor_not_cached (const FttCell * cell,
					FttDirection d)
{
  static gint neighbor_index[FTT_NEIGHBORS][FTT_CELLS]
#if FTT_2D
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2}};
#elif FTT_2D3
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2},
       {-1,-2,-3,-4},
       {-1,-2,-3,-4}};
#else  /* FTT_3D */
    = {{1,-1,3,-3,5,-5,7,-7},
       {-2,0,-4,2,-6,4,-8,6},
       {-3,-4,0,1,-7,-8,4,5},
       {2,3,-1,-2,6,7,-5,-6},
       {-5,-6,-7,-8,0,1,2,3},
       {4,5,6,7,-1,-2,-3,-4}};
#endif /* FTT_3D */
  gint n;
  FttCell * c;

  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (d < FTT_NEIGHBORS, NULL);

  if (FTT_CELL_IS_ROOT (cell))
    return ((struct _FttRootCell *) cell)->neighbors.c[d];

  n = neighbor_index[d][FTT_CELL_ID (cell)];
  if (n >= 0) /* neighbor belongs to same Oct */
    c = &(cell->parent->cell[n]);
  else {      /* neighbor belongs to neighboring Cell or Oct */
    c = cell->parent->neighbors.c[d];
    if (c != NULL && c->children != NULL)
      c = &(c->children->cell[- n - 1]);
  }
  if (c == NULL || FTT_CELL_IS_DESTROYED (c))
    return NULL;
  else
    return c;
}

/**
 * ftt_cell_neighbors:
 * @cell: a #FttCell.
 * @neighbors: a #FttCellNeighbors.
 *
 * Fills @neighbors with the neighbors of @cell.
 */
G_INLINE_FUNC
void ftt_cell_neighbors (const FttCell * cell,
			 FttCellNeighbors * neighbors)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (neighbors != NULL);

  if (!FTT_CELL_IS_LEAF (cell) && neighbors != &cell->children->neighbors) {
    memcpy (neighbors, &cell->children->neighbors, sizeof (FttCellNeighbors));
    return;
  }

  ftt_cell_neighbors_not_cached (cell, neighbors);
}

/**
 * ftt_cell_neighbor:
 * @cell: a #FttCell.
 * @d: a direction.
 *
 * Returns: the neighbor of @cell in direction @d or %NULL if @cell
 * has no neighbor in this direction.  
 */
G_INLINE_FUNC
FttCell * ftt_cell_neighbor (const FttCell * cell,
			     FttDirection d)
{
  g_return_val_if_fail (cell != NULL, NULL);
  g_return_val_if_fail (d < FTT_NEIGHBORS, NULL);

  if (!FTT_CELL_IS_LEAF (cell))
    return cell->children->neighbors.c[d];

  return ftt_cell_neighbor_not_cached (cell, d);
}

/**
 * ftt_cell_face:
 * @cell: a #FttCell.
 * @d: a direction.
 *
 * Returns: the face of @cell in direction @d.
 */
G_INLINE_FUNC
FttCellFace ftt_cell_face (FttCell * cell,
			   FttDirection d)
{
  FttCellFace f;

  g_return_val_if_fail (cell != NULL, f);

  f.cell = cell;
  f.neighbor = ftt_cell_neighbor (cell, d);
  f.d = d;

  return f;
}

/**
 * ftt_face_type:
 * @face: a #FttCellFace.
 *
 * Returns: the type of @face.
 */
G_INLINE_FUNC
FttFaceType ftt_face_type (const FttCellFace * face)
{
  g_return_val_if_fail (face != NULL, 0);

  if (face->neighbor == NULL)
    return FTT_BOUNDARY;
  if (ftt_cell_level (face->cell) > ftt_cell_level (face->neighbor))
    return FTT_FINE_COARSE;
  g_assert (ftt_cell_level (face->cell) == ftt_cell_level (face->neighbor));
  return FTT_FINE_FINE;
}

/**
 * ftt_cell_neighbor_is_brother:
 * @cell: a #FttCell.
 * @d: a #FttDirection.
 *
 * Returns: %TRUE if a (potential) neighbor of @cell in direction @d
 * and @cell would have the same parent, %FALSE otherwise.
 */
G_INLINE_FUNC
gboolean ftt_cell_neighbor_is_brother (FttCell * cell, 
				       FttDirection d)
{
  static gboolean b[FTT_CELLS][FTT_NEIGHBORS] = {
#if FTT_2D
    {1,0,0,1}, {0,1,0,1}, {1,0,1,0}, {0,1,1,0}
#elif FTT_2D3
    {1,0,0,1,0,0}, {0,1,0,1,0,0}, {1,0,1,0,0,0}, {0,1,1,0,0,0}
#else  /* 3D */
    {1,0,0,1,0,1}, {0,1,0,1,0,1}, {1,0,1,0,0,1}, {0,1,1,0,0,1},
    {1,0,0,1,1,0}, {0,1,0,1,1,0}, {1,0,1,0,1,0}, {0,1,1,0,1,0}
#endif /* 3D */
  };

  g_return_val_if_fail (cell != NULL, FALSE);
  
  if (FTT_CELL_IS_ROOT (cell))
    return FALSE;
  return b[FTT_CELL_ID (cell)][d];
}
