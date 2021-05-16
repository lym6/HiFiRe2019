#include "cpu_join.h"
#include "cpu_gpu_interface.h"
#include <openacc.h>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <math.h>

#define minn(X,Y) ((X) < (Y) ? (X) : (Y))
#define maxx(X,Y) ((X) > (Y) ? (X) : (Y))

using namespace std;

// Returns true if two rectangles (l1, r1) and (l2, r2) overlap 
bool doOverlap2(long long l1x,long long l1y, long long r1x, long long r1y, 
				long long l2x, long long l2y, long long r2x, long long r2y) 
{ 
    if ((l1x > r2x) || (l2x > r1x)) {
		return false;} 
  
    if ((r1y < l2y) || (r2y < l1y)) {
		return false;} 
	
	return true;
} 

int wakeGPUup(int tasks)
{
	int *A = (int *)malloc(sizeof(int) * tasks); 
	for(int i=0; i<tasks; i++){
		A[i]=0;
	}
	#pragma acc data copy(A[0:tasks])
	{	
		#pragma acc kernels
		{
			for(int i=0; i<tasks; i++)
			{
				A[i] = 4+ i*2;
			}
		}
	}
	return 0;
}

//transfer from doube to long long for vertices
gpc_vertex2* transfer_pointsvalues(int L1PolNum, int* L1VNum, int *L1VPrefixSum, gpc_vertex *hA)
{
    int lastL1PolVCount = L1VNum[L1PolNum - 1];
	int L1VCount = L1VPrefixSum[L1PolNum - 1] + lastL1PolVCount;
	
	gpc_vertex2 *a = (gpc_vertex2 *)malloc(sizeof(gpc_vertex2) * L1VCount);	 //inside for #1
	
	#pragma acc data copy(a[0:L1VCount]) copyin(hA[0:L1VCount])
	{	
		#pragma acc parallel loop
		for(int i = 0; i<L1VCount;i++){
			a[i].x = hA[i].x *10000000;
			a[i].y = hA[i].y *10000000; 
		}
	}
	return a;
}

//for wkt data
gpc_vertex2* transfer_pointsvalues2(int L1PolNum, int* L1VNum, long *L1VPrefixSum, Coordinate *hA)
{
    int lastL1PolVCount = L1VNum[L1PolNum - 1];
	long L1VCount = L1VPrefixSum[L1PolNum - 1] + lastL1PolVCount;
	
	gpc_vertex2 *a = (gpc_vertex2 *)malloc(sizeof(gpc_vertex2) * L1VCount);	 //inside for #1
	
	#pragma acc data copy(a[0:L1VCount]) copyin(hA[0:L1VCount])
	{	
		#pragma acc parallel loop
		for(int i = 0; i<L1VCount;i++){
			a[i].x = hA[i].x *10000000;
			a[i].y = hA[i].y *10000000; 
		}
	}
	return a;
}

//transfer from double to long long for bounding box
long long* transfer_boundingbox(int tasks, double *xy)
{	
	long long *a = (long long *)malloc(sizeof(long long) * tasks);
	#pragma acc data copy(a[0:tasks]) copyin(xy[0:tasks])
	{	
		#pragma acc parallel loop
		for(int i = 0; i<tasks;i++){
			a[i] = xy[i] *10000000;
		}
	}
	return a;
}

// calcaulte CBB
long long* calculateCBBys(int tasks, int *htaskSubId, int *htaskClipId, long long *rectn, long long* rect_queryn,
			int L1PolNum, int L2PolNum)
{	
	long long *b = (long long *)malloc(sizeof(long long) * tasks);
	#pragma acc data copy(b[0:tasks]) copyin(htaskSubId[0:tasks],htaskClipId[0:tasks],rectn[0:L1PolNum],rect_queryn[0:L2PolNum])
	{	
		#pragma acc parallel loop
		for(int i = 0; i<tasks;i++){
			int l1PolyId = htaskSubId[i];		
			int l2PolyId = htaskClipId[i];
			b[i] = max(rectn[l1PolyId], rect_queryn[l2PolyId]);
		}
	}
	return b;
}

long long* calculateCBByb(int tasks, int *htaskSubId, int *htaskClipId, long long *rectn, long long* rect_queryn,
			int L1PolNum, int L2PolNum)
{	
	long long *b = (long long *)malloc(sizeof(long long) * tasks);
	#pragma acc data copy(b[0:tasks]) copyin(htaskSubId[0:tasks],htaskClipId[0:tasks],rectn[0:L1PolNum],rect_queryn[0:L2PolNum])
	{	
		#pragma acc parallel loop
		for(int i = 0; i<tasks;i++){
			int l1PolyId = htaskSubId[i];		
			int l2PolyId = htaskClipId[i];
			b[i] = min(rectn[l1PolyId], rect_queryn[l2PolyId]);
		}
	}
	return b;
}

//PNP Function for Task that includes intersection points
int* processPNPTasks16(int pnpTaskPairs, Line3 *hA, Line3 *hB, int* L1VNum, int* L2VNum, int *L1VPrefixSum, 
			int *L2VPrefixSum, gpc_vertex2 *newHa, gpc_vertex2 *newHb, int *naSum, int *nbSum, int *naPrefix,
			int *nbPrefix)
{
    int lastL1PolVCount = L1VNum[pnpTaskPairs - 1];
	int L1VCount = L1VPrefixSum[pnpTaskPairs - 1] + lastL1PolVCount;
	int lastL2PolVCount = L2VNum[pnpTaskPairs - 1];
	int L2VCount = L2VPrefixSum[pnpTaskPairs - 1] + lastL2PolVCount;
	
    int lastL1PolVCount2 = naSum[pnpTaskPairs - 1];
	int L1VCount2 = naPrefix[pnpTaskPairs - 1] + lastL1PolVCount2;
	int lastL2PolVCount2 = nbSum[pnpTaskPairs - 1];
	int L2VCount2 = nbPrefix[pnpTaskPairs - 1] + lastL2PolVCount2;
	
	int *a = (int *)malloc(sizeof(int) * pnpTaskPairs);	 //inside for #1
	int *b = (int *)malloc(sizeof(int) * pnpTaskPairs);	 //inside for #2
	int *c = (int *)malloc(sizeof(int) * pnpTaskPairs);  //a+b	
	
	#pragma acc data copy(c[0:pnpTaskPairs],a[0:pnpTaskPairs],b[0:pnpTaskPairs]) copyin(hA[0:L1VCount], hB[0:L2VCount], L1VNum[0:pnpTaskPairs], L2VNum[0:pnpTaskPairs], L1VPrefixSum[0:pnpTaskPairs], L2VPrefixSum[0:pnpTaskPairs],newHa[0:L1VCount2],newHb[0:L2VCount2],naSum[0:pnpTaskPairs],nbSum[0:pnpTaskPairs],naPrefix[0:pnpTaskPairs],nbPrefix[0:pnpTaskPairs])
	{	
		#pragma acc parallel loop
		for (int i = 0; i < pnpTaskPairs; i++)      
		{	
			int L2PolVCount = L2VNum[i];	
			
			//char d;  // for the result of function InPoly
			
			//r means inside
			int r = 0;
			int L1PolVCount=naSum[i];
			
			#pragma acc loop reduction(+:r)
			for (int j = 0; j < L1PolVCount; j++)
			{   
				gpc_vertex2 test = newHa[naPrefix[i] + j];			 
				
				//d='z';

				long long x = 0;
				int Rcross = 0; 
				int Lcross = 0; 
				int Zcross = 0; 
				   				   
				#pragma acc loop reduction(+:Rcross) reduction(+:Lcross) reduction(+:Zcross) 
				for(int k = 0; k < L2PolVCount; k++ ) 
				{	
					//choose two points from second polygon
				    gpc_vertex2 vi = hB[L2VPrefixSum[i] + k].p1;
				    gpc_vertex2 vj = hB[L2VPrefixSum[i] + k].p2;
					
					if ( (vi.x==test.x) && (vi.y==test.y) )
					{
						Zcross++;  
					}

					if( ( (vi.y-test.y) > 0 ) != ( (vj.y-test.y) > 0 ) ) 
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) * (vi.y-test.y))
								/ ((vj.y-test.y) - (vi.y-test.y));
								
						if (x > 0) Rcross++;
					}

					if ( ( (vi.y-test.y) < 0 ) != ( (vj.y-test.y) < 0 ) )
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) * (vi.y-test.y))
								/ ((vj.y-test.y) - (vi.y-test.y)); 
								
						if (x < 0) Lcross++;
					}	
				}
				if(Zcross != 0 ) {
					//d = 'v';
				}
				else 
				{
					if(( Rcross % 2 )!=(Lcross % 2 )){
						//d='e';
					}
					else if( (Rcross % 2) == 1){
						//d='i';
						r++; 
					}
					else {
						//d='o';
					}
				}
			}
			//a[i] = r;	//inside
			c[i] = r;
		}
		
		#pragma acc parallel loop
		for (int i = 0; i < pnpTaskPairs; i++)      
		{	
			int L1PolVCount = L1VNum[i];	
			
			//char d;  // for the result of function InPoly
			
			//s means inside
			int s = 0;
			int L2PolVCount = nbSum[i];
			
			#pragma acc loop reduction(+:s) 
			for (int j = 0; j < L2PolVCount; j++)
			{   
				gpc_vertex2 test = newHb[nbPrefix[i] + j];			 
				
				//d='z';
				
				long long x = 0;
				int R2cross = 0; 
				int L2cross = 0; 
				int Z2cross = 0; 
				   				   
				#pragma acc loop reduction(+:R2cross) reduction(+:L2cross) reduction(+:Z2cross) 
				for(int k = 0; k < L1PolVCount; k++ ) 
				{	
					//choose two points from second polygon
				    gpc_vertex2 vi = hA[L1VPrefixSum[i] + k].p1;
				    gpc_vertex2 vj = hA[L1VPrefixSum[i] + k].p2;
				  
					if ( (vi.x==test.x) && (vi.y==test.y) )
					{
						Z2cross++; 
					}

					if( ( (vi.y-test.y) > 0 ) != ( (vj.y-test.y) > 0 ) ) 
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) * (vi.y-test.y))
								/ ((vj.y-test.y) - (vi.y-test.y));
								
						if (x > 0) R2cross++;
					}

					if ( ( (vi.y-test.y) < 0 ) != ( (vj.y-test.y) < 0 ) )
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) * (vi.y-test.y))
								/ ((vj.y-test.y) - (vi.y-test.y)); 
								
						if (x < 0) L2cross++;
					}	
				}
				if(Z2cross != 0 ) {
					//d = 'v';
				}
				else 
				{
					if(( R2cross % 2 )!=(L2cross % 2 )){
						//d='e';
					}
					else if( (R2cross % 2) == 1){
						//d='i';
						s++;
					}
					else {
						//d='o';
					}
				}
			}
			//b[i] = s;	//inside
			c[i] = c[i]+s;
		} 	
	}
	free(a);
	free(b);
	return c;
}

//PNP Function for Task that does not include intersection point
int* processPNPTasks19(int pnpTaskPairs, int *pnpL1TaskId, int *pnpL2TaskId, int L1PolNum, int L2PolNum, int* L1VNum, 
			int* L2VNum, int *L1VPrefixSum, int *L2VPrefixSum, gpc_vertex2 *hA, gpc_vertex2 *hB,int *newNum, int nNum,
			int* refineOne, int polygonInsideT, int* refineTwo, int polygonInsideT2)
{
    int lastL1PolVCount = L1VNum[L1PolNum - 1];
	int L1VCount = L1VPrefixSum[L1PolNum - 1] + lastL1PolVCount;
	int lastL2PolVCount = L2VNum[L2PolNum - 1];
	int L2VCount = L2VPrefixSum[L2PolNum - 1] + lastL2PolVCount;
	
	int realNumofInside = polygonInsideT+polygonInsideT2;
	int *c = (int *)malloc(sizeof(int) * realNumofInside);
	
	#pragma acc data copy(c[0:realNumofInside]) copyin(refineTwo[0:polygonInsideT2],refineOne[0:polygonInsideT], newNum[0:nNum],pnpL1TaskId[0:pnpTaskPairs], pnpL2TaskId[0:pnpTaskPairs], hA[0:L1VCount], hB[0:L2VCount], L1VNum[0:L1PolNum], L2VNum[0:L2PolNum], L1VPrefixSum[0:L1PolNum], L2VPrefixSum[0:L2PolNum])
	{	
		#pragma acc parallel loop
		for (int i = 0; i < polygonInsideT; i++)      
		{	
			int l2PolyId = pnpL2TaskId[refineOne[i]];
			int L2PolVCount = L2VNum[l2PolyId]-1;	
			
			int l1PolyId = pnpL1TaskId[refineOne[i]];
			int L1PolVCount = L1VNum[l1PolyId]-1;
			
			//char d;  // for the result of function InPoly
			
			//r means inside
			int r = 0;
			
			int newNNNum =  L1PolVCount;
			if(newNNNum >  10) 
			{
				newNNNum = 10;
			}
			
			#pragma acc loop reduction(+:r)
			for (int j = 0; j < newNNNum ; j++)    
			{   
				gpc_vertex2 test = hA[L1VPrefixSum[l1PolyId] + j];			 
				
				//d='z';

				long long x = 0;
				int Rcross = 0; // number of right edge/ray crossings 
				int Lcross = 0; // number of left edge/ray crossings 
				int Zcross = 0; // For the point is a vertex of a polygon
				 				   
				#pragma acc loop reduction(+:Rcross) reduction(+:Lcross) reduction(+:Zcross) 
				for(int k = 0; k < L2PolVCount; k++ ) 
				{
					//choose two points from second polygon
				    gpc_vertex2 vi = hB[L2VPrefixSum[l2PolyId] + k + 1];
				    gpc_vertex2 vj = hB[L2VPrefixSum[l2PolyId] + k ];
					
					if ( (vi.x==test.x) && (vi.y==test.y) )
					{
						Zcross++;  
					}
					
					if( ( (vi.y-test.y) > 0 ) != ( (vj.y-test.y) > 0 ) ) 
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) * (vi.y-test.y))
								/ ((vj.y-test.y) -(vi.y-test.y));
								
						if (x > 0) Rcross++;
					}

					if ( ( (vi.y-test.y) < 0 ) != ( (vj.y-test.y) < 0 ) )
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) *(vi.y-test.y))
								/ ((vj.y-test.y) - (vi.y-test.y)); 
								
						if (x < 0) Lcross++;
					}	
				}
				if(Zcross != 0 ) {
					//d = 'v';
				}
				else 
				{
					if(( Rcross % 2 )!=(Lcross % 2 )){
						//d='e';
					}
					else if( (Rcross % 2) == 1){
						//d='i';
						r++;
					}
					else {
						//d='o';
					}
				}
			}  
			if(r>3){
				c[i] = L1PolVCount;	//inside
			}
			else{c[i] =0;}
		}
		
		#pragma acc parallel loop
		for (int i = 0; i < polygonInsideT2; i++)      
		{	
			int l1PolyId = pnpL1TaskId[refineTwo[i]];
			int L1PolVCount = L1VNum[l1PolyId]-1;	
			
			int l2PolyId = pnpL2TaskId[refineTwo[i]];
			int L2PolVCount = L2VNum[l2PolyId]-1;	
			
			//char d;  // for the result of function InPoly
			
			//s means inside
			int s = 0;
			
			int newNNNum2 =  L2PolVCount;
			if(newNNNum2 >  10)
			{
				newNNNum2 = 10;
			}
			
			#pragma acc loop reduction(+:s) 
			for (int j = 0; j < newNNNum2; j++)    
			{   
				gpc_vertex2 test = hB[L2VPrefixSum[l2PolyId] + j];			 
				
				//d='z';
				
				long long x = 0;
				int R2cross = 0; // number of right edge/ray crossings 
				int L2cross = 0; // number of left edge/ray crossings 
				int Z2cross = 0; // For the point is a vertex of a polygon
				   				   
				#pragma acc loop reduction(+:R2cross) reduction(+:L2cross) reduction(+:Z2cross) 
				for(int k = 0; k < L1PolVCount; k++ ) 
				{	
					//choose two points from second polygon
				    gpc_vertex2 vi = hA[L1VPrefixSum[l1PolyId] + k + 1];
				    gpc_vertex2 vj = hA[L1VPrefixSum[l1PolyId] + k];
				  
					if ( (vi.x==test.x) && (vi.y==test.y) )
					{
						Z2cross++;  
					}

					if( ( (vi.y-test.y) > 0 ) != ( (vj.y-test.y) > 0 ) ) 
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) * (vi.y-test.y))
								/ ((vj.y-test.y) - (vi.y-test.y));
								
						if (x > 0) R2cross++;
					}

					if ( ( (vi.y-test.y) < 0 ) != ( (vj.y-test.y) < 0 ) )
					{
						x = ((vi.x-test.x) * (vj.y-test.y) - (vj.x-test.x) * (vi.y-test.y))
								/ ((vj.y-test.y) - (vi.y-test.y)); 
								
						if (x < 0) L2cross++;
					}	
				}
				if(Z2cross != 0 ) {
					//d = 'v';
				}
				else 
				{
					if(( R2cross % 2 )!=(L2cross % 2 )){
						//d='e';
					}
					else if( (R2cross % 2) == 1){
						//d='i';
						s++;
					}
					else {
						//d='o';
					}
				}
			}
			if(s>3){
				c[polygonInsideT+i] =L2PolVCount;	//inside
			}
			else{c[polygonInsideT+i]=0;}
		} 	
	}
	return c;
}

//LSI function with PolySketch
int* calculateMbrPQ3(int tasks, int *pnpL1TaskId, int *pnpL2TaskId, int L1PolNum, int L2PolNum, int* L1VNum, 
			int* L2VNum, int *L1VPrefixSum, int *L2VPrefixSum, gpc_vertex2 *ha1, gpc_vertex2 *hb1, 
			int *numOfPartL1,int *numOfPartL2,int *lastNumL1,int *lastNumL2, long long *prefixPQ1,
			long long *prefixPQ2, int *cellsizeL1, int *cellsizeL2)
{	
	int lastL1PolVCount = L1VNum[L1PolNum - 1];
	int L1VCount = L1VPrefixSum[L1PolNum - 1] + lastL1PolVCount;
	int lastL2PolVCount = L2VNum[L2PolNum - 1];
	int L2VCount = L2VPrefixSum[L2PolNum - 1] + lastL2PolVCount;
	
	int numfpoints = 0;
	int *f = (int *)malloc(sizeof(int) * tasks);
	//assign memory for intersection points, may need be changed 
	gpc_vertex2 *fpoints = (gpc_vertex2 *)malloc(sizeof(gpc_vertex2) *80*tasks);
	
	long long tasknum1 = prefixPQ1[tasks-1]+ numOfPartL1[tasks-1];
	long long tasknum2 = prefixPQ2[tasks-1]+ numOfPartL2[tasks-1];
	
	long long *xmax1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *xmin1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *ymax1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *ymin1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *xmax2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	long long *xmin2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	long long *ymax2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	long long *ymin2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	
	#pragma acc data copy(f[0:tasks],fpoints[0:tasks*80]) create(xmax1[0:tasknum1],xmin1[0:tasknum1],ymax1[0:tasknum1],ymin1[0:tasknum1],xmax2[0:tasknum2],xmin2[0:tasknum2],ymax2[0:tasknum2],ymin2[0:tasknum2]) copyin(pnpL1TaskId[0:tasks], pnpL2TaskId[0:tasks], ha1[0:L1VCount], hb1[0:L2VCount], L1VNum[0:L1PolNum], L2VNum[0:L2PolNum], L1VPrefixSum[0:L1PolNum], L2VPrefixSum[0:L2PolNum],numOfPartL1[0:tasks],numOfPartL2[0:tasks],lastNumL1[0:tasks],lastNumL2[0:tasks],prefixPQ1[0:tasks],prefixPQ2[0:tasks],cellsizeL1[0:tasks],cellsizeL2[0:tasks])
	{	
		#pragma acc parallel
		#pragma acc loop
		for(int i =0; i<tasks;i++)
		{
			int l1PolyId = pnpL1TaskId[i];		 
			int l2PolyId = pnpL2TaskId[i];		
			
			int partL1 = numOfPartL1[i];
			int L1Last = lastNumL1[i];
			int partL2 = numOfPartL2[i];
			int L2Last = lastNumL2[i];
			int tempprefix1 = prefixPQ1[i];
			int tempprefix2 = prefixPQ2[i];

			int cellsize1 =  cellsizeL1[i];
			int cellsize2 =  cellsizeL2[i];
			
			#pragma acc loop
			for(int j = 0; j<partL1; j++)
			{
				long long xmaxtemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].x;
				long long xmintemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].x;
				long long ymaxtemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].y;
				long long ymintemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].y;
				
				int maxtemp = cellsize1;
				if (j == (partL1-1)){maxtemp = L1Last;}
				
				#pragma acc loop seq
				for(int k = 0; k<maxtemp;k++)
				{
					gpc_vertex2 aPoint = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j + k ]; 
					if(xmaxtemp<aPoint.x){xmaxtemp=aPoint.x;}
					if(xmintemp>aPoint.x){xmintemp=aPoint.x;}
					if(ymaxtemp<aPoint.y){ymaxtemp=aPoint.y;}
					if(ymintemp>aPoint.y){ymintemp=aPoint.y;}
				}
				xmax1[tempprefix1+j] = xmaxtemp;
				xmin1[tempprefix1+j] = xmintemp;
				ymax1[tempprefix1+j] = ymaxtemp;
				ymin1[tempprefix1+j] = ymintemp;
			}
			
			#pragma acc loop			
			for(int j = 0; j<partL2; j++)
			{
				long long xmaxtemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].x;
				long long xmintemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].x;
				long long ymaxtemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].y;
				long long ymintemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].y;
				
				int maxtemp = cellsize2;
				if (j == (partL2-1)){maxtemp = L2Last;}
				
				#pragma acc loop seq
				for(int k = 0; k<maxtemp;k++)
				{
					gpc_vertex2 aPoint = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j + k ];
					if(xmaxtemp<aPoint.x){xmaxtemp=aPoint.x;}
					if(xmintemp>aPoint.x){xmintemp=aPoint.x;}
					if(ymaxtemp<aPoint.y){ymaxtemp=aPoint.y;}
					if(ymintemp>aPoint.y){ymintemp=aPoint.y;}
				}
				xmax2[tempprefix2+j] = xmaxtemp;
				xmin2[tempprefix2+j] = xmintemp;
				ymax2[tempprefix2+j] = ymaxtemp;
				ymin2[tempprefix2+j] = ymintemp;
			}	
			
			int q = 0;
			#pragma acc loop reduction(+:q)
			for(int j = 0; j<partL1; j++)
			{
				int maxtemp1 = cellsize1;
				if (j == (partL1-1)){maxtemp1 = L1Last;}
				
				#pragma acc loop 
				for(int k = 0; k <partL2; k++)
				{		
					int maxtemp2 = cellsize2;
					if (k == (partL2-1)){maxtemp2 = L2Last;}
					
					if(doOverlap2(xmin1[tempprefix1+j],ymin1[tempprefix1+j],xmax1[tempprefix1+j],ymax1[tempprefix1+j],xmin2[tempprefix2+k], ymin2[tempprefix2+k], xmax2[tempprefix2+k], ymax2[tempprefix2+k]))
					{
						#pragma acc loop seq
						for(int jjj=0;jjj<(maxtemp1-1);jjj++)
						{
							gpc_vertex2 test = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j + jjj ]; 
							gpc_vertex2 test2 = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j + jjj+1 ]; 
							#pragma acc loop seq
							for(int kkk=0;kkk<(maxtemp2-1);kkk++)
							{
								gpc_vertex2 vi = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*k + kkk ];
								gpc_vertex2 vj = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*k + kkk+1 ];
								
								long long o1 = (vi.x-test.x)*(test2.y-test.y) - (vi.y-test.y)*(test2.x-test.x); 
								long long o2 = (vj.x-test.x)*(test2.y-test.y) - (vj.y-test.y)*(test2.x-test.x); 
								long long o3 = (test.x-vi.x)*(vj.y-vi.y) - (test.y-vi.y)*(vj.x-vi.x);			  
								long long o4 = (test2.x-vi.x)*(vj.y-vi.y) - (test2.y-vi.y)*(vj.x-vi.x);		 
				 
								//check intersections
								if(((o1 < 0)!= (o2 < 0))&&((o3 < 0)!= (o4 < 0)))
								{
									q++;
							
									gpc_vertex2 P_intersection;
							
									//if(test2.x == test.x){test2.x = test2.x+1;}
									//if(vj.x == vi.x){vj.x = vj.x+1;}
										
									long long l1m = ((test2.y - test.y) / (test2.x - test.x));
									long long l1c = (test.y) - l1m*(test.x);
									long long l2m = ((vj.y - vi.y) / (vj.x - vi.x));
									long long l2c = (vi.y) - l2m*(vi.x);
										
									//if(l1m == l2m){l1m = l1m+1;}
							
									P_intersection.x = (l2c - l1c)/(l1m - l2m);
									P_intersection.y = l1m*P_intersection.x + l1c;
							
									#pragma acc atomic
									{
										fpoints[numfpoints]=P_intersection;
										numfpoints++;
									}
								}
								else if((o1==0)&&(vi.x<=max(test.x,test2.x))&&(vi.x>=min(test.x,test2.x))&&
										(vi.y<=max(test.y,test2.y))&&(vi.y>=min(test.y,test2.y))){
									q++;
								}
								else if((o3==0)&&(test.x<=max(vi.x,vj.x))&&(test.x>=min(vi.x,vj.x))&&
											(test.y<=max(vi.y,vj.y))&&(test.y>=min(vi.y,vj.y))){
									q++;
								}
							}
						}
					}
				}
			}
			f[i] = q;
		}	 
	}	
	free(xmax1);
	free(xmin1);
	free(ymax1);
	free(ymin1);
	free(xmax2);
	free(xmin2);
	free(ymax2);
	free(ymin2);
	return f;
}

//LSI function with PolySketch for wkt data
int* calculateMbrPQ6(int tasks, int *pnpL1TaskId, int *pnpL2TaskId, int L1PolNum, int L2PolNum, int* L1VNum, 
			int* L2VNum, long *L1VPrefixSum, long *L2VPrefixSum, gpc_vertex2 *ha1, gpc_vertex2 *hb1, 
			int *numOfPartL1,int *numOfPartL2,int *lastNumL1,int *lastNumL2, long long *prefixPQ1,
			long long *prefixPQ2, int *cellsizeL1, int *cellsizeL2)
{	
	int lastL1PolVCount = L1VNum[L1PolNum - 1];
	long L1VCount = L1VPrefixSum[L1PolNum - 1] + lastL1PolVCount;
	int lastL2PolVCount = L2VNum[L2PolNum - 1];
	long L2VCount = L2VPrefixSum[L2PolNum - 1] + lastL2PolVCount;
	
	int numfpoints = 0;
	int *f = (int *)malloc(sizeof(int) * tasks);
	//assign memory for intersection points, may need be changed 
	gpc_vertex2 *fpoints = (gpc_vertex2 *)malloc(sizeof(gpc_vertex2) *50*tasks);
	
	long long tasknum1 = prefixPQ1[tasks-1]+ numOfPartL1[tasks-1];
	long long tasknum2 = prefixPQ2[tasks-1]+ numOfPartL2[tasks-1];
	
	long long *xmax1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *xmin1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *ymax1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *ymin1 = (long long *)malloc(tasknum1 * sizeof(long long)); 
	long long *xmax2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	long long *xmin2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	long long *ymax2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	long long *ymin2 = (long long *)malloc(tasknum2 * sizeof(long long)); 
	
	#pragma acc data copy(f[0:tasks],fpoints[0:tasks*50]) create(xmax1[0:tasknum1],xmin1[0:tasknum1],ymax1[0:tasknum1],ymin1[0:tasknum1],xmax2[0:tasknum2],xmin2[0:tasknum2],ymax2[0:tasknum2],ymin2[0:tasknum2]) copyin(pnpL1TaskId[0:tasks], pnpL2TaskId[0:tasks], ha1[0:L1VCount], hb1[0:L2VCount], L1VNum[0:L1PolNum], L2VNum[0:L2PolNum], L1VPrefixSum[0:L1PolNum], L2VPrefixSum[0:L2PolNum],numOfPartL1[0:tasks],numOfPartL2[0:tasks],lastNumL1[0:tasks],lastNumL2[0:tasks],prefixPQ1[0:tasks],prefixPQ2[0:tasks],cellsizeL1[0:tasks],cellsizeL2[0:tasks])
	{	
		#pragma acc parallel
		#pragma acc loop
		for(int i =0; i<tasks;i++)
		{
			int l1PolyId = pnpL1TaskId[i];		 
			int l2PolyId = pnpL2TaskId[i];		 
			
			int partL1 = numOfPartL1[i];
			int L1Last = lastNumL1[i];
			int partL2 = numOfPartL2[i];
			int L2Last = lastNumL2[i];
			int tempprefix1 = prefixPQ1[i];
			int tempprefix2 = prefixPQ2[i];

			int cellsize1 =  cellsizeL1[i];
			int cellsize2 =  cellsizeL2[i];
			
			#pragma acc loop
			for(int j = 0; j<partL1; j++)
			{
				long long xmaxtemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].x;
				long long xmintemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].x;
				long long ymaxtemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].y;
				long long ymintemp = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j ].y;
				
				int maxtemp = cellsize1;
				if (j == (partL1-1)){maxtemp = L1Last;}
				
				#pragma acc loop seq
				for(int k = 0; k<maxtemp;k++)
				{
					gpc_vertex2 aPoint = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j + k ]; 
					if(xmaxtemp<aPoint.x){xmaxtemp=aPoint.x;}
					if(xmintemp>aPoint.x){xmintemp=aPoint.x;}
					if(ymaxtemp<aPoint.y){ymaxtemp=aPoint.y;}
					if(ymintemp>aPoint.y){ymintemp=aPoint.y;}
				}
				xmax1[tempprefix1+j] = xmaxtemp;
				xmin1[tempprefix1+j] = xmintemp;
				ymax1[tempprefix1+j] = ymaxtemp;
				ymin1[tempprefix1+j] = ymintemp;
			}
			
			#pragma acc loop			
			for(int j = 0; j<partL2; j++)
			{
				long long xmaxtemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].x;
				long long xmintemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].x;
				long long ymaxtemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].y;
				long long ymintemp = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j ].y;
				
				int maxtemp = cellsize2;
				if (j == (partL2-1)){maxtemp = L2Last;}
				
				#pragma acc loop seq
				for(int k = 0; k<maxtemp;k++)
				{
					gpc_vertex2 aPoint = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*j + k ];
					if(xmaxtemp<aPoint.x){xmaxtemp=aPoint.x;}
					if(xmintemp>aPoint.x){xmintemp=aPoint.x;}
					if(ymaxtemp<aPoint.y){ymaxtemp=aPoint.y;}
					if(ymintemp>aPoint.y){ymintemp=aPoint.y;}
				}
				xmax2[tempprefix2+j] = xmaxtemp;
				xmin2[tempprefix2+j] = xmintemp;
				ymax2[tempprefix2+j] = ymaxtemp;
				ymin2[tempprefix2+j] = ymintemp;
			}	
			
			int q = 0;
			#pragma acc loop reduction(+:q)
			for(int j = 0; j<partL1; j++)
			{
				int maxtemp1 = cellsize1;
				if (j == (partL1-1)){maxtemp1 = L1Last;}
				
				#pragma acc loop 
				for(int k = 0; k <partL2; k++)
				{		
					int maxtemp2 = cellsize2;
					if (k == (partL2-1)){maxtemp2 = L2Last;}
					
					if(doOverlap2(xmin1[tempprefix1+j],ymin1[tempprefix1+j],xmax1[tempprefix1+j],ymax1[tempprefix1+j],xmin2[tempprefix2+k], ymin2[tempprefix2+k], xmax2[tempprefix2+k], ymax2[tempprefix2+k]))
					{	
						#pragma acc loop seq
						for(int jjj=0;jjj<(maxtemp1-1);jjj++)
						{
							gpc_vertex2 test = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j + jjj ]; 
							gpc_vertex2 test2 = ha1[L1VPrefixSum[l1PolyId] + (cellsize1-1)*j + jjj+1 ]; 
							#pragma acc loop seq
							for(int kkk=0;kkk<(maxtemp2-1);kkk++)
							{
								gpc_vertex2 vi = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*k + kkk ];
								gpc_vertex2 vj = hb1[L2VPrefixSum[l2PolyId] + (cellsize2-1)*k + kkk+1 ];
									
								long long o1 = (vi.x-test.x)*(test2.y-test.y) - (vi.y-test.y)*(test2.x-test.x); 
								long long o2 = (vj.x-test.x)*(test2.y-test.y) - (vj.y-test.y)*(test2.x-test.x); 
								long long o3 = (test.x-vi.x)*(vj.y-vi.y) - (test.y-vi.y)*(vj.x-vi.x);			 
								long long o4 = (test2.x-vi.x)*(vj.y-vi.y) - (test2.y-vi.y)*(vj.x-vi.x);		  
				  
								//check intersections
								if(((o1 < 0)!= (o2 < 0))&&((o3 < 0)!= (o4 < 0)))
								{
									q++;
							
									gpc_vertex2 P_intersection;
							
									long long l1m = ((test2.y - test.y) / (test2.x - test.x));
									long long l1c = (test.y) - l1m*(test.x);
									long long l2m = ((vj.y - vi.y) / (vj.x - vi.x));
									long long l2c = (vi.y) - l2m*(vi.x);
						
									P_intersection.x = (l2c - l1c)/(l1m - l2m);
									P_intersection.y = l1m*P_intersection.x + l1c;
							
									#pragma acc atomic
									{
										fpoints[numfpoints]=P_intersection;
										numfpoints++;
									}
								}
								else if((o1==0)&&(vi.x<=max(test.x,test2.x))&&(vi.x>=min(test.x,test2.x))&&
										(vi.y<=max(test.y,test2.y))&&(vi.y>=min(test.y,test2.y))){
									q++;
								}
								else if((o3==0)&&(test.x<=max(vi.x,vj.x))&&(test.x>=min(vi.x,vj.x))&&
											(test.y<=max(vi.y,vj.y))&&(test.y>=min(vi.y,vj.y))){
									q++;
								}
							}
						}
					}
				}
			}
			f[i] = q;
		}	 
	}	
	
	free(xmax1);
	free(xmin1);
	free(ymax1);
	free(ymin1);
	free(xmax2);
	free(xmin2);
	free(ymax2);
	free(ymin2);
	free(fpoints);
	return f;
}
