#ifndef __CPU_JOIN_H_INCLUDED__
#define __CPU_JOIN_H_INCLUDED__

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "shapefil.h"
#include "gridpartitioning.h"
#include "rTreeIndex.h"
#include "gpc.h"
#include "cpu_gpu_interface.h"

typedef struct
{
  gpc_vertex v;
  float alphas,alphac;
  short int bIndex;
  short int cIndex;
}gpc_intersection;

typedef struct 
{
  gpc_vertex start;
  gpc_vertex end;
}Edge;

typedef struct
{
	int taskId;
	Edge* edgeList;
	int edgeCount;
	int isValid;
}TASK_INFO, *LP_TASK_INFO;

typedef struct Rect_double{
        double boundary[4];
} Rect_double;

int readShapefile(int argc, char *argv[]);

double my_difftime();

void rtreeBuildingAndSearch(int num_contours, double *rect1, double *rect2, double *rect3,double *rect4, 
 int *id_base, int numOfQuerys,double *rect1_query,double *rect2_query, double *rect3_query, double *rect4_query);

int MySearchCallback(int id, void* arg);

int SHPReadMBR1( SHPHandle psSHP, int startIndex, int endIndex, PolyRect ** mbrs);

void test_gpu_polygon_kernel();

void distRead(int index,int size,SHPHandle	hSHP,SHPObject	***psShape,int *num);

void destoryObjects(SHPObject **psShape,int num);

void printObjects(SHPObject	**psShape,int num);

void convert(double *rect1, double *rect2, double *rect3, double *rect4,SHPObject **psShape_base, 
int *num_base, int * prefix_base,bBox *baseBoxes,int cellsPerProcess,int sum_mbrs_overlay,
SHPHandle hSHP_base,double * minX,double *minY);

void parse(char * filename);
int SHPReadMBR1( SHPHandle hSHP, int startIndex, int endIndex, PolyRect ** mbrs);
void convertCPU(Rect_double * rects_data,SHPObject **psShape_base, int *num_base, int * prefix_base,bBox *baseBoxes,
int cellsPerProcess,int sum_mbrs_overlay,SHPHandle hSHP_base,double * minX,double * minY);
                        
#endif		
 