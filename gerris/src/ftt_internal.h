//
// Created by erich on 18/06/2020.
//

#ifndef GISS_FTT_INTERNAL_H
#define GISS_FTT_INTERNAL_H

#include "ftt.h"

static void traverse_face(FttCell *cell, gpointer *datum);

static void traverse_all_faces(FttCell *cell, gpointer *datum);

static void traverse_all_direct_faces(FttCell *cell, gpointer *datum);

static void traverse_face_direction(FttCell *cell, gpointer *datum);

static void traverse_face_component(FttCell *cell, gpointer *datum);

static void reset_flag(FttCell *cell);

#endif //GISS_FTT_INTERNAL_H
