/*
===========================================================================

Project:   Generic Polygon Clipper

           A new algorithm for calculating the difference, intersection,
           exclusive-or or union of arbitrary polygon sets.

File:      gpc.h
Author:    Alan Murta (email: gpc@cs.man.ac.uk)
Version:   2.32
Date:      17th December 2004

Copyright: (C) Advanced Interfaces Group,
           University of Manchester.

           This software is free for non-commercial use. It may be copied,
           modified, and redistributed provided that this copyright notice
           is preserved on all copies. The intellectual property rights of
           the algorithms used reside with the University of Manchester
           Advanced Interfaces Group.

           You may not use this software, in whole or in part, in support
           of any commercial product without the express consent of the
           author.

           There is no warranty or other guarantee of fitness of this
           software for any purpose. It is provided solely "as is".

===========================================================================
*/

#ifndef __gpc_h
#define __gpc_h

#include <stdio.h>

//for geos
#include <fstream>
#include <iostream>
#include "geos/geom/Coordinate.h"

using namespace std;
using namespace geos::geom;


/*
===========================================================================
                               Constants
===========================================================================
*/

/* Increase GPC_EPSILON to encourage merging of near coincident edges    */

#define GPC_EPSILON (DBL_EPSILON)

#define GPC_VERSION "2.32"

/*
===========================================================================
                           Public Data Types
===========================================================================
*/

typedef enum                        /* Set operation type                */
{
  GPC_DIFF,                         /* Difference                        */
  GPC_INT,                          /* Intersection                      */
  GPC_XOR,                          /* Exclusive or                      */
  GPC_UNION                         /* Union                             */
} gpc_op;

typedef struct                      /* Polygon vertex structure          */
{
  double              x;            /* Vertex x component                */
  double              y;            /* vertex y component                */
} gpc_vertex;

/*
		Some types are my data types
*/

typedef struct                      
{
  long long              x;            
  long long              y;            
} gpc_vertex2;

typedef struct                      
{
  long long              x;          
  long long              y;           
  int z;
} gpc_vertex3;

typedef struct                     
{
	gpc_vertex p1;
	gpc_vertex p2;
} Line;

typedef struct                     
{
	Coordinate p1;
	Coordinate p2;
} Line2;

typedef struct                     
{
	gpc_vertex2 p1;
	gpc_vertex2 p2;
} Line3;

typedef struct                      /* Vertex list structure             */
{
  int                 maxAllocated;
  int                 num_vertices; /* Number of vertices in list        */
  gpc_vertex         *vertex;       /* Vertex array pointer              */
} gpc_vertex_list;

typedef struct                      /* Polygon set structure             */
{
  int                num_contours; /* Number of contours in polygon     */
  int                *hole;         /* Hole / external contour flags     */
  gpc_vertex_list    *contour;      /* Contour array pointer             */
  gpc_vertex lowerBBox;             /* Lower bounding box                */
  gpc_vertex upperBBox;             /* Upper bounding box                */
  char *fID;                        /* Face id                           */ 
} gpc_polygon;

typedef struct
{
 int num_elementsBasePoly;
 int num_allocatedBasePoly;
 gpc_polygon basePolygon;
 int overlayPolyCount;             /* Number of overlay polygon count      */
 gpc_polygon *overlayPolygons;     /* Intersecting polygons                */
 int *overlayPolyIndices;          /* Index of the intersecting clip polys */
}gpc_overlayPoly;

typedef struct                      /* Tristrip set structure            */
{
  int                 num_strips;   /* Number of tristrips               */
  gpc_vertex_list    *strip;        /* Tristrip array pointer            */
} gpc_tristrip;


typedef struct                      /* Vertex list structure             */
{
  int                 maxAllocated;
  int                 num_vertices; /* Number of vertices in list        */
  Coordinate         *vertex;       /* Vertex array pointer              */
} gpc_vertex_list2;

typedef struct                      /* Polygon set structure             */
{
  int                num_contours; /* Number of contours in polygon     */
  int                *hole;         /* Hole / external contour flags     */
  gpc_vertex_list2    *contour;      /* Contour array pointer             */
  Coordinate lowerBBox;             /* Lower bounding box                */
  Coordinate upperBBox;             /* Upper bounding box                */
  char *fID;                        /* Face id                           */ 
} gpc_polygon2;

typedef struct
{
 int num_elementsBasePoly;
 int num_allocatedBasePoly;
 gpc_polygon2 basePolygon;
 int overlayPolyCount;             /* Number of overlay polygon count      */
 gpc_polygon2 *overlayPolygons;     /* Intersecting polygons                */
 int *overlayPolyIndices;          /* Index of the intersecting clip polys */
}gpc_overlayPoly2;

/*
===========================================================================
                       Public Function Prototypes
===========================================================================
*/

void gpc_read_polygon        (FILE            *infile_ptr, 
                              int              read_hole_flags,
                              gpc_polygon     *polygon);

void gpc_read_overlaypolygon(FILE *fp, int read_hole_flags, gpc_overlayPoly *p);
                              
void gpc_write_polygon       (FILE            *outfile_ptr,
                              int              write_hole_flags,
                              gpc_polygon     *polygon);

void gpc_add_contour         (gpc_polygon     *polygon,
                              gpc_vertex_list *contour,
                              int              hole);

void gpc_polygon_clip        (gpc_op           set_operation,
                              gpc_polygon     *subject_polygon,
                              gpc_polygon     *clip_polygon,
                              gpc_polygon     *result_polygon);

void gpc_tristrip_clip       (gpc_op           set_operation,
                              gpc_polygon     *subject_polygon,
                              gpc_polygon     *clip_polygon,
                              gpc_tristrip    *result_tristrip);

void gpc_polygon_to_tristrip (gpc_polygon     *polygon,
                              gpc_tristrip    *tristrip);

void gpc_free_polygon        (gpc_polygon     *polygon);

void gpc_free_tristrip       (gpc_tristrip    *tristrip);

#endif

/*
===========================================================================
                           End of file: gpc.h
===========================================================================
*/
