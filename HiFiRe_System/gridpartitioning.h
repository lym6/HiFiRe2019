#ifndef GRIDPARTITION_H
#define GRIDPARTITION_H

#include <stdlib.h>
#include <stdio.h>

#include "rTreeIndex.h"

#define NUM_COORDINATES 4
#define PACKING_Box_Size 6

typedef struct Rect_int{
	int boundary[4];
} Rect_int;

/* MBR of a polygon */
typedef struct PolyRect
{
    int poly_id;
    int cellid;       
                      
    struct Rect mbr;   /* Rect is defined as a structure with 4
                          doubles: xmin,ymin,xmax,ymax */
} PolyRect;

/* bBox holds MBRs of base/overlay polygons for a given cell */
typedef struct bBox
{
 int count;     
 int allocated; 
 int processorIncharge;  
 struct PolyRect *rects;   
} bBox;

typedef struct searchCallBack{
	int id;
	bBox * boxes;
	PolyRect *rect;
}searchCallBack;

int gridPartition(int id,int ncells,int nprocs,PolyRect *baseRect,int basePcount,PolyRect *overlayRect,int overlayPcount,bBox **baseBoxes,bBox **overlayBoxes,double lowX,double lowY,double upX,double upY,double **grid);

#endif
