#ifndef __CPU_GPU_INTERFACE_H_INCLUDED__
#define __CPU_GPU_INTERFACE_H_INCLUDED__

#include "gpc.h"
#include "cpu_join.h"

//for geos
#include <fstream>
#include <iostream>
#include "geos/geom/Coordinate.h"

using namespace std;
using namespace geos::geom;


typedef struct
{
  double *xminArr;
  double *yminArr;
  double *xmaxArr;
  double *ymaxArr;
}MBR;

int ST_Intersect(int L1PolNum, int L2PolNum, int* L1VNum, int* L2VNum, int *L1VPrefixSum, int *L2VPrefixSum, 
					gpc_vertex *L1Coords, gpc_vertex *L2Coords,  MBR *L1MBR, MBR* L2MBR, 
					int numTasks, int *L1TaskId, int *L2TaskId, int *taskResult);
 			  
int* processPNPTasks16(int pnpTaskPairs, Line3 *hA, Line3 *hB, int* L1VNum, int* L2VNum, int *L1VPrefixSum, int *L2VPrefixSum,
					gpc_vertex2 *newHa, gpc_vertex2 *newHb, int *naSum, int *nbSum, int *naPrefix, int *nbPrefix);
						
int* processPNPTasks19(int pnpTaskPairs, int *pnpL1TaskId, int *pnpL2TaskId, int L1PolNum, int L2PolNum, int* L1VNum, int* L2VNum,
                      int *L1VPrefixSum, int *L2VPrefixSum, gpc_vertex2 *hA, gpc_vertex2 *hB,int *newNum, int nNum,
					  int* refineOne, int polygonInsideT, int* refineTwo, int polygonInsideT2);
			   
int wakeGPUup(int tasks);					  			
						
gpc_vertex2* transfer_pointsvalues(int L1PolNum, int* L1VNum, int *L1VPrefixSum, gpc_vertex *hA);
gpc_vertex2* transfer_pointsvalues2(int L1PolNum, int* L1VNum, long *L1VPrefixSum, Coordinate *hA);	

long long* transfer_boundingbox(int tasks, double *xy);					

long long* calculateCBBys(int tasks, int *htaskSubId, int *htaskClipId, long long *rectn, long long* rect_queryn,int L1PolNum, int L2PolNum);
long long* calculateCBByb(int tasks, int *htaskSubId, int *htaskClipId, long long *rectn, long long* rect_queryn,int L1PolNum, int L2PolNum);
					  
int* calculateMbrPQ3(int tasks, int *pnpL1TaskId, int *pnpL2TaskId, int L1PolNum, int L2PolNum, int* L1VNum, int* L2VNum,
                      int *L1VPrefixSum, int *L2VPrefixSum, gpc_vertex2 *ha1, gpc_vertex2 *hb1, 
					  int *numOfPartL1,int *numOfPartL2,int *lastNumL1,int *lastNumL2, long long *prefixPQ1,long long *prefixPQ2,
					  int *cellsizeL1, int *cellsizeL2);
	  
int* calculateMbrPQ6(int tasks, int *pnpL1TaskId, int *pnpL2TaskId, int L1PolNum, int L2PolNum, int* L1VNum, int* L2VNum,
                      long *L1VPrefixSum, long *L2VPrefixSum, gpc_vertex2 *ha1, gpc_vertex2 *hb1, 
					  int *numOfPartL1,int *numOfPartL2,int *lastNumL1,int *lastNumL2, long long *prefixPQ1,long long *prefixPQ2,
					  int *cellsizeL1, int *cellsizeL2);
								  
int verifyLayer(int PolNum, int* VNum, int *VPrefixSum, gpc_vertex *vertices);

#endif

