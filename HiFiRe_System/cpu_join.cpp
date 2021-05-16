#include "cpu_join.h"
#include <fstream>
#include <list>
#include <iostream>
#include <vector> 
#include <omp.h>
#include <cstring> 	
#include <string> 	
#include <math.h>
#include <openacc.h>
#include <mutex>
#include <thread>

#include "geos/geom/Geometry.h"
#include "geos/io/WKTReader.h"
#include "geos/geom/Coordinate.h"
#include "geos/geom/Envelope.h"
#include "geos/geom/LineString.h"
#include "geos/geom/Point.h"
#include "geos/geom/Polygon.h"

using namespace std;
using namespace geos::geom;


bool doOverlapyyy(long long l1y, long long r1y, long long l2y, long long r2y) 
{ 
    if ((r1y < l2y) || (r2y < l1y)) {
		return false;} 
	
	return true;
} 

//for geos data sets
geos::io::WKTReader wktreader;
gpc_overlayPoly2 *subjectPoly2 = NULL; 
gpc_polygon2 *clipPoly2 = NULL;

std::mutex push_mutex;
//the class is for layer(for geos)
typedef struct
{
int polNum;
vector<int> * vNum;  
vector<long> * prefixSum; 
vector<Coordinate> * coord; 
vector<Coordinate> * mbr; 
vector<Envelope> * mbr2; 
}polygonLayerTmp;

typedef struct
{
int polNum;
int* vNum;  
long* prefixSum; 
Coordinate *coord; 
Coordinate* mbr; 
Envelope* mbr2; 
}polygonLayer;

/* base layer polygons */
gpc_overlayPoly *subjectPoly = NULL; 

int    num_elementsSubPoly1;
int    num_allocatedSubPoly = 0; 

/* overlay layer polygons */
gpc_polygon *clipPoly = NULL;
 
int    num_elementsClipPoly = 0;
int    num_allocatedClipPoly = 0; 
int    totalMyHit = 0;


int SHPReadMBR1( SHPHandle psSHP, int startIndex, int endIndex, PolyRect ** mbrs){
	PolyRect * bounding_boxes;
	int num=(endIndex-startIndex+1);
	bounding_boxes=(PolyRect *) malloc(num*sizeof(PolyRect));
	int i,j;

/* -------------------------------------------------------------------- */
/*      Read the record.                                                */
/* -------------------------------------------------------------------- */
    for(i=startIndex;i<=endIndex;i++){
		if( psSHP->sHooks.FSeek( psSHP->fpSHP, psSHP->panRecOffset[i]+12, 0 ) != 0 )
		{
			char str[128];
			sprintf( str,
					 "Error in fseek() reading object from .shp file at offset %u",
					 psSHP->panRecOffset[i]+12);

			psSHP->sHooks.Error( str );
			return -1;
		}
		for(j=0;j<4;j++){
			if( psSHP->sHooks.FRead( &(bounding_boxes[i-startIndex].mbr.boundary[j]), sizeof(double), 1, psSHP->fpSHP ) != 1 )
			{
				char str[128];
				sprintf( str,
						 "Error in fread() reading object of size %u at offset %u from .shp file",
						 4*sizeof(double), psSHP->panRecOffset[i]+12 );

				psSHP->sHooks.Error( str );
				return -1;
			}			
		}
		bounding_boxes[i-startIndex].poly_id=i;
    }
    (*mbrs)=bounding_boxes;
    return 0;
}

void convert(double * rect1, double * rect2, double * rect3, double * rect4,SHPObject **psShape_base, 
int *num_base, int * prefix_base,bBox *baseBoxes,int cellsPerProcess,int sum_mbrs_overlay,
SHPHandle hSHP_base,double * minX, double *minY)
{
	int i;
	int prefix=0;
	
	for(i =0 ; i< baseBoxes[0].count;i++)
	{             
		rect1[i+prefix]=baseBoxes[0].rects[i].mbr.boundary[0];
		rect2[i+prefix]=baseBoxes[0].rects[i].mbr.boundary[1];
		rect3[i+prefix]=baseBoxes[0].rects[i].mbr.boundary[2];
		rect4[i+prefix]=baseBoxes[0].rects[i].mbr.boundary[3];
		if(rect1[i+prefix]<(*minX))
		{
			(*minX)=rect1[i+prefix];
		}
		if(rect2[i+prefix]<(*minY))
		{
			(*minY)=rect2[i+prefix];
		}
		psShape_base[i+prefix] = SHPReadObject( hSHP_base, baseBoxes[0].rects[i].poly_id);
		if( psShape_base[i+prefix] == NULL )
		{
			fprintf( stderr,"Unable to read shape %d, terminating object reading.\n",baseBoxes[0].rects[i].poly_id );
			exit(1);
		}
		num_base[i+prefix]=psShape_base[i+prefix]->nVertices;
		if((i+prefix)==0)
		{
			prefix_base[i+prefix]=0;
		}	
		else
		{
			prefix_base[i+prefix]=prefix_base[i+prefix-1]+num_base[i+prefix-1];
		}
	}
	prefix+=baseBoxes[0].count;
}

void destoryObjects(SHPObject **psShape,int num){
	int i;
	for(i=0;i<num;i++){
	    SHPDestroyObject(psShape[i]);
	}
}

double my_difftime ()
{
    struct timeval tp;
    struct timezone tzp;
    int i;

    i = gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

int MySearchCallback(int id, void* arg) 
 {
       int* idOfbase = (int *)arg;
	
       if(subjectPoly[*idOfbase].num_elementsBasePoly == subjectPoly[*idOfbase].num_allocatedBasePoly)
       {
			if (subjectPoly[*idOfbase].num_allocatedBasePoly == 0)
            subjectPoly[*idOfbase].num_allocatedBasePoly = 100; 
			else
			{
				subjectPoly[*idOfbase].num_allocatedBasePoly += 100; 
			}
			void *_tmp = realloc(subjectPoly[*idOfbase].overlayPolyIndices, (subjectPoly[*idOfbase].num_allocatedBasePoly 
			* sizeof(int)));
			if (!_tmp)
			{
				printf("************ ERROR: Couldn't realloc memory!***********\n");
				fflush(stdout);
			}
			subjectPoly[*idOfbase].overlayPolyIndices = (int *)_tmp;   
       }
	   subjectPoly[*idOfbase].overlayPolyIndices[subjectPoly[*idOfbase].num_elementsBasePoly] = id-1;
       subjectPoly[*idOfbase].num_elementsBasePoly++;
       return 1;
}

void rtreeBuildingAndSearch(int num_contours, double *rect1, double *rect2, double *rect3,double *rect4, 
 int *id_base, int numOfQuerys,double *rect1_query,double *rect2_query, double *rect3_query, double *rect4_query)
{
   subjectPoly = (gpc_overlayPoly *)malloc(num_contours * 
			  sizeof(gpc_overlayPoly));
  
    struct Node* root = RTreeNewIndex();
    int nhits,i;
			  
	for(i=0; i<numOfQuerys; i++)
	{
	    struct Rect rect = {rect1_query[i], rect2_query[i], rect3_query[i], rect4_query[i]};
        RTreeInsertRect(&rect, i+1, &root, 0);
    }
	int c = 0;
	for (c= 0; c < num_contours; c++)
    {
      gpc_overlayPoly p;
      p.num_elementsBasePoly = 0;
      p.num_allocatedBasePoly = 0;
      p.overlayPolyIndices = NULL;
      subjectPoly[c] = p;
	}
					
 	int a;
	for(a = 0; a <num_contours; a++)
    {
	  nhits = 0;
	  struct Rect search_rect = {rect1[a], rect2[a], rect3[a], rect4[a]};
	  nhits = RTreeSearch(root, &search_rect, MySearchCallback, &a);
	}
}

//function for geos
polygonLayerTmp* populateLayerData(list<Geometry*>* geoms);
void readGeomsFromStr(vector<string> *vstr, list<Geometry*> *geoms);

void combineLayers(polygonLayerTmp* send, polygonLayerTmp* recv){
	recv->polNum += send->polNum;

	for(vector<int >::iterator itr = send->vNum->begin(); itr != send->vNum->end(); ++itr){
		recv->vNum->push_back(*itr);
	}

    int tmpPrefix = 0;
    if(! recv->prefixSum->empty()){
		tmpPrefix = recv->prefixSum->back();
	}
	
	for(vector<long >::iterator itr = send->prefixSum->begin(); itr != send->prefixSum->end(); ++itr){
		if(!((tmpPrefix!=0) && ((*itr)==0)))
		recv->prefixSum->push_back(*itr+tmpPrefix);
	}
	
	for(vector<Coordinate >::iterator itr = send->coord->begin(); itr != send->coord->end(); ++itr){
		recv->coord->push_back(*itr);
	}
	
	for(vector<Coordinate>::iterator itr = send->mbr->begin(); itr != send->mbr->end(); ++itr){
		recv->mbr->push_back(*itr);
	}

	for(vector<Envelope>::iterator itr = send->mbr2->begin(); itr != send->mbr2->end(); ++itr){
		recv->mbr2->push_back(*itr);
	}
}

void parsing(vector<string> *vstr1, int start1, int end1, vector<string> *vstr2, int start2, int end2, polygonLayerTmp* layer1, polygonLayerTmp* layer2){

    vector<string> *localStr1 = new vector<string>(vstr1->begin()+start1, vstr1->begin()+end1-1);

    vector<string> *localStr2 = new vector<string>(vstr2->begin()+start2, vstr2->begin()+end2-1);

    list<Geometry*>* lGeos1 = new list<Geometry*>;
    list<Geometry*>* lGeos2 = new list<Geometry*>;

    readGeomsFromStr(localStr1, lGeos1);
    readGeomsFromStr(localStr2, lGeos2);

    polygonLayerTmp* tmpLayer1 = populateLayerData(lGeos1);
	polygonLayerTmp* tmpLayer2 = populateLayerData(lGeos2);
	
	push_mutex.lock();
	combineLayers(tmpLayer1,layer1);
	combineLayers(tmpLayer2,layer2);
	push_mutex.unlock();
}

void convertMBRToFloats(const Envelope* v, vector<Coordinate> *vertVect)
{	
	Coordinate P1,P2;
	P1.x =(double)v->getMinX();
	P1.y =(double)v->getMinY();
	P2.x =(double)v->getMaxX();
	P2.y =(double)v->getMaxY();
    vertVect->push_back(P1);
    vertVect->push_back(P2);
}		

void convertToFloats(const LineString* vertices, vector<Coordinate> *vertVect)
{
	size_t numPoints = vertices->getNumPoints();
  
	for(size_t i = 0; i<numPoints; i++) 
	{
		Point* pt = vertices->getPointN(i);
		Coordinate P3;
		P3.x = pt->getX();
		P3.y = pt->getY();
		vertVect->push_back(P3);
	} 
}

vector<long> * prefixSum(polygonLayerTmp *layer)
{
	long numPoly = layer->polNum;
	vector<int>* vNum = layer->vNum;
   	vector<long>* prefixsum = new vector<long>;
	prefixsum->push_back(0);
   
	for(long i = 1; i < numPoly+1; i++) 
	{
		prefixsum->push_back(prefixsum->at(i-1)+vNum->at(i-1));
	} 
	return prefixsum;
}

void gpuHelperForMultiPolygon(Geometry *geom, vector<Coordinate> *verticesVec, 
		vector<Coordinate> *envVec, vector<Envelope> *gpuEnvInLongVector, vector<int> *vNumVector)
{
	size_t numGeoms = geom->getNumGeometries();
         
	for(size_t i = 0; i < numGeoms; i++) 
	{
		const Geometry* inner = geom->getGeometryN(i);
		
        const Polygon* poly = dynamic_cast<const Polygon*>(inner);      
		const LineString *innerLinestring = poly->getExteriorRing();
        
		vNumVector->push_back(innerLinestring->getNumPoints());		
		convertToFloats(innerLinestring, verticesVec);   	
		convertMBRToFloats(poly->getEnvelopeInternal(), envVec);		
	}
}

void gpuHelperForPolygon(Geometry *geom, vector<Coordinate> *verticesVec, 
		vector<Coordinate> *envVec, vector<Envelope> *gpuEnvInLongVector, vector<int> *vNumVector)
{
	Polygon* poly = dynamic_cast<Polygon*>(geom); 
	const LineString *linestring = poly->getExteriorRing();
   
	vNumVector->push_back(linestring->getNumPoints());
	convertToFloats(linestring, verticesVec);     
	convertMBRToFloats(poly->getEnvelopeInternal(), envVec);		
}

polygonLayerTmp* populateLayerData(list<Geometry*>* geoms)
{	
	polygonLayerTmp* layer = (polygonLayerTmp*)malloc(1* sizeof(polygonLayerTmp)); 
   
	vector<Coordinate> *runningVector = new vector<Coordinate>();	
	vector<Coordinate> *runningEnvVector = new vector<Coordinate>();  
	vector<Envelope> *gpuEnvInLongVector = new vector<Envelope>();
	vector<int> *vNumVector = new vector<int>();  

	int numPolygons = 0;  
      
	for(list<Geometry*>::iterator it = geoms->begin() ; it != geoms->end(); ++it) 
	{
		Geometry *geom = *it;
		GeometryTypeId typeId = geom->getGeometryTypeId();
    
		switch(typeId)
		{
            case GEOS_POLYGON:
			{
				gpuHelperForPolygon(geom, runningVector, runningEnvVector, gpuEnvInLongVector, vNumVector);
				numPolygons++;
			} 
			break;
       
			case GEOS_MULTIPOLYGON:
			{
				gpuHelperForMultiPolygon(geom, runningVector, runningEnvVector, gpuEnvInLongVector, vNumVector);
			}
			break;
       
			case GEOS_GEOMETRYCOLLECTION:
			{
				size_t numGeoms = geom->getNumGeometries();
				numPolygons += numGeoms;       
       
				for(size_t i = 0; i < numGeoms; i++)
				{
					GeometryTypeId typeId = geom->getGeometryTypeId();
         
					switch(typeId)
					{
						case GEOS_POLYGON:
						{
							gpuHelperForPolygon(geom, runningVector, runningEnvVector, gpuEnvInLongVector, 
								vNumVector);
							numPolygons++;
						} 
						break;
       
						case GEOS_MULTIPOLYGON:
						{
							gpuHelperForMultiPolygon(geom, runningVector, runningEnvVector, gpuEnvInLongVector, 
								vNumVector);
						}
						break;
					} 
				} 
			}
			break;  
		}
	}
	layer->polNum = vNumVector->size();     
	layer->vNum = vNumVector;     
	layer->prefixSum = prefixSum(layer);	
	layer->coord = runningVector;   
	layer->mbr = runningEnvVector;  
	layer->mbr2 = gpuEnvInLongVector;
	return layer;
}

void readGeomsFromStr(vector<string> *vstr, list<Geometry*> *geoms){
    geos::io::WKTReader wktreader;

    for(vector<string >::iterator itr = vstr->begin(); itr != vstr->end(); ++itr){
        string tmpStr = *itr;
        try{
		    Geometry* tmpGeo = NULL;
		    tmpGeo = wktreader.read(tmpStr);
		    if(tmpGeo != NULL && tmpGeo->isValid())
			    geoms->push_back(tmpGeo);
		}catch(exception &e){
		}
    }
    
}
//above functions get the basic information for one data set 

int MySearchCallback2(int id, void* arg) 
 {
       int* idOfbase = (int *)arg;
	
       if(subjectPoly2[*idOfbase].num_elementsBasePoly == subjectPoly2[*idOfbase].num_allocatedBasePoly)
       {
			if (subjectPoly2[*idOfbase].num_allocatedBasePoly == 0)
            subjectPoly2[*idOfbase].num_allocatedBasePoly = 100; 
			else
			{
				subjectPoly2[*idOfbase].num_allocatedBasePoly += 100; 
			}
			void *_tmp = realloc(subjectPoly2[*idOfbase].overlayPolyIndices, (subjectPoly2[*idOfbase].num_allocatedBasePoly 
			* sizeof(int)));
			if (!_tmp)
			{
				printf("************ ERROR: Couldn't realloc memory!***********\n");
				fflush(stdout);
			}
			subjectPoly2[*idOfbase].overlayPolyIndices = (int *)_tmp;   
       }
	   subjectPoly2[*idOfbase].overlayPolyIndices[subjectPoly2[*idOfbase].num_elementsBasePoly] = id-1;
       subjectPoly2[*idOfbase].num_elementsBasePoly++;
       return 1;
}

void rtreeBuildingAndSearch2(int num_contours, double *rect1, double *rect2, double *rect3,double *rect4, 
 int numOfQuerys,double *rect1_query,double *rect2_query, double *rect3_query, double *rect4_query)
{
   subjectPoly2 = (gpc_overlayPoly2 *)malloc(num_contours * 
			  sizeof(gpc_overlayPoly2));
  
    struct Node* root = RTreeNewIndex();
    int nhits,i;
			  
	for(i=0; i<numOfQuerys; i++)
	{
	    struct Rect rect = {rect1_query[i], rect2_query[i], rect3_query[i], rect4_query[i]};
        RTreeInsertRect(&rect, i+1, &root, 0);
    }
	int c = 0;
	
	for (c= 0; c < num_contours; c++)
    {
      gpc_overlayPoly2 p;
      p.num_elementsBasePoly = 0;
      p.num_allocatedBasePoly = 0;
      p.overlayPolyIndices = NULL;
      subjectPoly2[c] = p;
	}
					
 	int a;
	
	for(a = 0; a <num_contours; a++)
    {
	  nhits = 0;
	  struct Rect search_rect = {rect1[a], rect2[a], rect3[a], rect4[a]};
	  nhits = RTreeSearch(root, &search_rect, MySearchCallback2, &a);
	}
}


int localProcessing35(int num_contours, double *rect1, double *rect2, double *rect3,double *rect4, int *id_base,
 int numOfQuerys,double *rect1_query,double *rect2_query, double *rect3_query,
 double *rect4_query,int *id_query,int *num_base,int *prefix_base,gpc_vertex *vertex_base,
 int *num_overlay,int *prefix_overlay,gpc_vertex *vertex_overlay,int *hole_base_cpu,
 int *hole_overlay_cpu,int m,int M,int nprocs,int minX,int minY,char * name)
{
	printf("%d,%d\n",num_contours,numOfQuerys); 
	
	double starttime0, endtime0;
    double difference0;	
    starttime0 = my_difftime();
	
    rtreeBuildingAndSearch(num_contours, rect1, rect2, rect3, rect4,
     id_base, numOfQuerys, rect1_query, rect2_query, rect3_query, rect4_query);

	int tasks = 0; 
    int clipIndex;
    int a;
             
	for(a = 0; a < num_contours; a++)
	{
        tasks = tasks + subjectPoly[a].num_elementsBasePoly;
    }
       
    int *htaskSubId = (int *)malloc(tasks * sizeof(int));
    int *htaskClipId = (int *)malloc(tasks * sizeof(int)); 
    int counter = 0;
	
    for(a = 0; a < num_contours; a++)
	{
        for(clipIndex = 0; clipIndex < subjectPoly[a].num_elementsBasePoly; clipIndex++)
        {
          htaskSubId[counter] = a;
          htaskClipId[counter] = subjectPoly[a].overlayPolyIndices[clipIndex];
          counter = counter + 1;
        }
    }
	
	endtime0 = my_difftime();
	difference0 = endtime0 - starttime0;
	printf("Rtree time taken =  %f\t \n",difference0);
	printf("total tasks = %d \n", tasks);
	
	//transfer from double to long long 
	gpc_vertex2* ha1 = transfer_pointsvalues(num_contours, num_base, prefix_base, vertex_base);
	gpc_vertex2* hb1 = transfer_pointsvalues(numOfQuerys, num_overlay, prefix_overlay, vertex_overlay);
	long long* rect1n = transfer_boundingbox(num_contours, rect1);
	long long* rect2n = transfer_boundingbox(num_contours, rect2);
	long long* rect3n = transfer_boundingbox(num_contours, rect3);
	long long* rect4n = transfer_boundingbox(num_contours, rect4);
	long long* rect1_queryn = transfer_boundingbox(numOfQuerys, rect1_query);
	long long* rect2_queryn = transfer_boundingbox(numOfQuerys, rect2_query);
	long long* rect3_queryn = transfer_boundingbox(numOfQuerys, rect3_query);
	long long* rect4_queryn = transfer_boundingbox(numOfQuerys, rect4_query);
	
	double starttime1, endtime1;
    double difference1;	
    starttime1 = my_difftime();
	
	//get basic information for sketch
	int *numOfPartL1 = (int *)malloc(tasks * sizeof(int));   //how many parts of 1st polygon
	int *numOfPartL2 = (int *)malloc(tasks * sizeof(int)); 
	int *lastNumL1 = (int *)malloc(tasks * sizeof(int)); 	//for the final part, how many points are inside
	int *lastNumL2 = (int *)malloc(tasks * sizeof(int)); 	
	int *cellsizeL1 = (int *)malloc(tasks * sizeof(int));   
	int *cellsizeL2 = (int *)malloc(tasks * sizeof(int)); 
	
	//get sketch basic information
	#pragma omp parallel num_threads(32)
	{
		#pragma omp for schedule(static)
		for(int i =0; i<tasks;i++)
		{
			int l1PolyId = htaskSubId[i];		 
			int L1PolVCount = num_base[l1PolyId];  
			int l2PolyId = htaskClipId[i];		
			int L2PolVCount = num_overlay[l2PolyId];  
			
			int cellsize1 = 15;
			int cellsize2 = 15;
			if(L1PolVCount<400){cellsize1=5;}
			if(L2PolVCount<400){cellsize2=5;}
			
			cellsizeL1[i] = cellsize1;
			cellsizeL2[i] = cellsize2;
			
			int partL1 = ((L1PolVCount-1)/(cellsize1-1))+1;
			int lastNumL1sub = L1PolVCount - (cellsize1-1)*(partL1-1);
			if((partL1==0)||(partL1==1)){partL1=1;lastNumL1sub=L1PolVCount;}
			if(lastNumL1sub<=3)
			{
				partL1 = partL1-1;
				lastNumL1sub = lastNumL1sub+(cellsize1-1);
			}
			numOfPartL1[i] = partL1;
			lastNumL1[i] = lastNumL1sub;
		
			int partL2 = ((L2PolVCount-1)/(cellsize2-1))+1;
			int lastNumL2sub = L2PolVCount - (cellsize2-1)*(partL2-1);
			if((partL2==0)||(partL2==1)){partL2=1;lastNumL2sub=L2PolVCount;}
			if(lastNumL2sub<=3)
			{
				partL2 = partL2-1;
				lastNumL2sub = lastNumL2sub+(cellsize2-1);
			}
			numOfPartL2[i] = partL2;
			lastNumL2[i] = lastNumL2sub;
		}
	}

	long long *prefixPQ1 = (long long *)malloc(tasks * sizeof(long long)); 	//the number of numOfPartL1 
	long long *prefixPQ2 = (long long *)malloc(tasks * sizeof(long long)); 
	
	prefixPQ1[0] = 0;
	prefixPQ2[0] = 0;
	
	for(int i =1;i<tasks;i++)
	{	
		prefixPQ1[i] = prefixPQ1[i-1]+ numOfPartL1[i-1];
		prefixPQ2[i] = prefixPQ2[i-1]+ numOfPartL2[i-1];
	}
	
	endtime1 = my_difftime();
	difference1 = endtime1 - starttime1;
	printf("Get sketch basic information =  %f\t \n",difference1);
	
	wakeGPUup(tasks);
	
	double starttime2, endtime2;
    double difference2;	
    starttime2 = my_difftime();
	
	//LSI function with PolySketch 
	int *SIresult = calculateMbrPQ3(tasks, htaskSubId, htaskClipId, num_contours, numOfQuerys, num_base, num_overlay, 
						prefix_base, prefix_overlay, ha1, hb1, numOfPartL1, numOfPartL2, lastNumL1, lastNumL2, 
						prefixPQ1, prefixPQ2,cellsizeL1, cellsizeL2);
						
	endtime2 = my_difftime();
	difference2 = endtime2 - starttime2;
	printf("Segment Intersection Function time taken =  %f\t \n",difference2);
	
	//for department gpu
	wakeGPUup(tasks);
	
	double starttime3, endtime3;
    double difference3;	
    starttime3 = my_difftime();
	
	int lastL1PolVCount = num_base[num_contours - 1];
	int L1VCount = prefix_base[num_contours - 1] + lastL1PolVCount;
	int lastL2PolVCount = num_overlay[numOfQuerys - 1];
	int L2VCount = prefix_overlay[numOfQuerys - 1] + lastL2PolVCount;
	
	long long* ys = calculateCBBys(tasks, htaskSubId, htaskClipId, rect2n, rect2_queryn,num_contours,numOfQuerys);
	long long* yb = calculateCBByb(tasks, htaskSubId, htaskClipId, rect4n, rect4_queryn,num_contours,numOfQuerys);
	
	//for checking points inside or outside
	int nNum = 0;  //the number of tasks some points of one ploygon maybe inside another polygon
	int *newNum = (int *)malloc(tasks * sizeof(int));  
	int polygonInsideT = 0;
	int *refineOne = (int *)malloc(tasks * sizeof(int));
	int polygonInsideT2 = 0;
	int *refineTwo = (int *)malloc(tasks * sizeof(int));
	int dDiscard = 0 ;
	
	for(int i =0; i < tasks; i++) 
	{	
		int l1PolyId = htaskSubId[i];		
		int l2PolyId = htaskClipId[i];
	
		if(SIresult[i] != 0){
			newNum[nNum] = i;
			nNum++;
		}	
		else if((rect1n[l1PolyId]>=rect1_queryn[l2PolyId])&&(rect2n[l1PolyId]>=rect2_queryn[l2PolyId])&&(rect3n[l1PolyId]<=rect3_queryn[l2PolyId])&&(rect4n[l1PolyId]<=rect4_queryn[l2PolyId]))
		{		
			refineOne[polygonInsideT] = i;
			polygonInsideT++;  //layer 1 inside layer 2 
		}
		else if((rect1n[l1PolyId]<=rect1_queryn[l2PolyId])&&(rect2n[l1PolyId]<=rect2_queryn[l2PolyId])&&(rect3n[l1PolyId]>=rect3_queryn[l2PolyId])&&(rect4n[l1PolyId]>=rect4_queryn[l2PolyId]))
		{
			refineTwo[polygonInsideT2] = i;
			polygonInsideT2++;  //layer 2 inside layer 1 
		}
		else{dDiscard++;}
	}
	//printf("%d,%d,%d\n",nNum,polygonInsideT+polygonInsideT2,dDiscard);
	
// cell = 8........................
	vector<Line3> vecSegmentA11;   
	vector<Line3> vecSegmentA21;   
	vector<Line3> vecSegmentA12;   
	vector<Line3> vecSegmentA22;   
	vector<Line3> vecSegmentA13;   
	vector<Line3> vecSegmentA23;   
	vector<Line3> vecSegmentA14;   
	vector<Line3> vecSegmentA24;
	
	vector<Line3> vecSegmentA15;   
	vector<Line3> vecSegmentA25;   
	vector<Line3> vecSegmentA16;   
	vector<Line3> vecSegmentA26;   
	vector<Line3> vecSegmentA17;   
	vector<Line3> vecSegmentA27;   
	vector<Line3> vecSegmentA18;   
	vector<Line3> vecSegmentA28; 

	vecSegmentA11.reserve(L1VCount);   
	vecSegmentA21.reserve(L2VCount);
	vecSegmentA12.reserve(L1VCount);   
	vecSegmentA22.reserve(L2VCount);
	vecSegmentA13.reserve(L1VCount);   
	vecSegmentA23.reserve(L2VCount); 
	vecSegmentA14.reserve(L1VCount);   
	vecSegmentA24.reserve(L2VCount); 
	vecSegmentA15.reserve(L1VCount);   
	vecSegmentA25.reserve(L2VCount); 
	vecSegmentA16.reserve(L1VCount);   
	vecSegmentA26.reserve(L2VCount); 
	vecSegmentA17.reserve(L1VCount);   
	vecSegmentA27.reserve(L2VCount); 
	vecSegmentA18.reserve(L1VCount);   
	vecSegmentA28.reserve(L2VCount); 
 	
	vector<gpc_vertex2> vecPointP11;  
	vector<gpc_vertex2> vecPointP21;
	vector<gpc_vertex2> vecPointP12;  
	vector<gpc_vertex2> vecPointP22; 
	vector<gpc_vertex2> vecPointP13;  
	vector<gpc_vertex2> vecPointP23;   
	vector<gpc_vertex2> vecPointP14;  
	vector<gpc_vertex2> vecPointP24;   
	
	vector<gpc_vertex2> vecPointP15;  
	vector<gpc_vertex2> vecPointP25;   
	vector<gpc_vertex2> vecPointP16;  
	vector<gpc_vertex2> vecPointP26;   
	vector<gpc_vertex2> vecPointP17;  
	vector<gpc_vertex2> vecPointP27;   
	vector<gpc_vertex2> vecPointP18;  
	vector<gpc_vertex2> vecPointP28;  

	vecPointP11.reserve(L1VCount);  
	vecPointP21.reserve(L2VCount);
	vecPointP12.reserve(L1VCount);  
	vecPointP22.reserve(L2VCount);	
	vecPointP13.reserve(L1VCount);  
	vecPointP23.reserve(L2VCount);	
	vecPointP14.reserve(L1VCount);  
	vecPointP24.reserve(L2VCount);	
	vecPointP15.reserve(L1VCount);  
	vecPointP25.reserve(L2VCount);	
	vecPointP16.reserve(L1VCount);  
	vecPointP26.reserve(L2VCount);	
	vecPointP17.reserve(L1VCount);  
	vecPointP27.reserve(L2VCount);	
	vecPointP18.reserve(L1VCount);  
	vecPointP28.reserve(L2VCount);	
	
	vector<int> ngSumAsub1;  
	vector<int> nhSumAsub1;  
	vector<int> ngSumAsub2;  
	vector<int> nhSumAsub2;  
	vector<int> ngSumAsub3;  
	vector<int> nhSumAsub3;  
	vector<int> ngSumAsub4;  
	vector<int> nhSumAsub4;  
	
	vector<int> ngSumAsub5;  
	vector<int> nhSumAsub5;  
	vector<int> ngSumAsub6;  
	vector<int> nhSumAsub6;  
	vector<int> ngSumAsub7;  
	vector<int> nhSumAsub7;  
	vector<int> ngSumAsub8;  
	vector<int> nhSumAsub8;  
	
	vector<int> naSumAsub1;  
	vector<int> nbSumAsub1; 
	vector<int> naSumAsub2;  
	vector<int> nbSumAsub2; 
	vector<int> naSumAsub3;  
	vector<int> nbSumAsub3; 
	vector<int> naSumAsub4;  
	vector<int> nbSumAsub4; 

	vector<int> naSumAsub5;  
	vector<int> nbSumAsub5; 
	vector<int> naSumAsub6;  
	vector<int> nbSumAsub6; 
	vector<int> naSumAsub7;  
	vector<int> nbSumAsub7; 
	vector<int> naSumAsub8;  
	vector<int> nbSumAsub8; 
	
	vector<int> realTaskNumsub2;
	
	#pragma omp parallel num_threads(32)
	{	
		#pragma omp for schedule(static)
		for(int idx =0; idx<nNum; idx++){
		
			int currentid = newNum[idx]; //  this means current id
			int originalid = currentid;
			
			vector<int> realTaskNumsubsub2;  //every new task number has current id
			vector<Line3> ngA1;   
			vector<Line3> nhA1;   
			vector<Line3> ngA2;   
			vector<Line3> nhA2;   
			vector<Line3> ngA3;   
			vector<Line3> nhA3;   
			vector<Line3> ngA4;   
			vector<Line3> nhA4;   
			
			vector<Line3> ngA5;   
			vector<Line3> nhA5;   
			vector<Line3> ngA6;   
			vector<Line3> nhA6;   
			vector<Line3> ngA7;   
			vector<Line3> nhA7;   
			vector<Line3> ngA8;   
			vector<Line3> nhA8; 

			vector<gpc_vertex2> naA1; 
			vector<gpc_vertex2> nbA1; 
			vector<gpc_vertex2> naA2; 
			vector<gpc_vertex2> nbA2;  
			vector<gpc_vertex2> naA3; 
			vector<gpc_vertex2> nbA3;  
			vector<gpc_vertex2> naA4; 
			vector<gpc_vertex2> nbA4;  
			
			vector<gpc_vertex2> naA5; 
			vector<gpc_vertex2> nbA5;  
			vector<gpc_vertex2> naA6; 
			vector<gpc_vertex2> nbA6;  
			vector<gpc_vertex2> naA7; 
			vector<gpc_vertex2> nbA7;  
			vector<gpc_vertex2> naA8; 
			vector<gpc_vertex2> nbA8;  
		
			vector<int> ngSumAsubsub1;  
			vector<int> nhSumAsubsub1;  
			vector<int> ngSumAsubsub2;  
			vector<int> nhSumAsubsub2;  
			vector<int> ngSumAsubsub3; 
			vector<int> nhSumAsubsub3; 
			vector<int> ngSumAsubsub4;  
			vector<int> nhSumAsubsub4; 
			
			vector<int> ngSumAsubsub5;  
			vector<int> nhSumAsubsub5;  
			vector<int> ngSumAsubsub6;  
			vector<int> nhSumAsubsub6;  
			vector<int> ngSumAsubsub7; 
			vector<int> nhSumAsubsub7; 
			vector<int> ngSumAsubsub8;  
			vector<int> nhSumAsubsub8; 
			
			vector<int> naSumAsubsub1;  
			vector<int> nbSumAsubsub1;
			vector<int> naSumAsubsub2;  
			vector<int> nbSumAsubsub2; 
			vector<int> naSumAsubsub3;  
			vector<int> nbSumAsubsub3; 
			vector<int> naSumAsubsub4;  
			vector<int> nbSumAsubsub4; 
			
			vector<int> naSumAsubsub5;  
			vector<int> nbSumAsubsub5; 
			vector<int> naSumAsubsub6;  
			vector<int> nbSumAsubsub6; 
			vector<int> naSumAsubsub7;  
			vector<int> nbSumAsubsub7; 
			vector<int> naSumAsubsub8;  
			vector<int> nbSumAsubsub8; 
			
			int l1PolyId = htaskSubId[originalid];		 
			int L1PolVCount = num_base[l1PolyId]-1;  
			int l2PolyId = htaskClipId[originalid];		
			int L2PolVCount = num_overlay[l2PolyId]-1;   
			
			long long y1 = ys[originalid];
			long long y9 = yb[originalid];  //y9>y1
			long long intervaly = (y9-y1)/8;
			
			long long y2 = y1+intervaly;
			long long y3 = y1+2*intervaly;
			long long y4 = y1+3*intervaly;
			long long y5 = y1+4*intervaly;
			long long y6 = y1+5*intervaly;
			long long y7 = y1+6*intervaly;
			long long y8 = y1+7*intervaly;
			
			if(intervaly==0)
			continue;
			
			for (int j = 0; j < L1PolVCount; j++)   //check how many points of first polygon inside common box
			{
				gpc_vertex2 test = ha1[prefix_base[l1PolyId] + j + 1];    //choose the first vertex
				gpc_vertex2 testTwo = ha1[prefix_base[l1PolyId] + j ];  //choose the second vertex
				
				if(((test.y<y1)&&(testTwo.y<y1))||((test.y>y9)&&(testTwo.y>y9)))
				continue;
			
				if((testTwo.y>=y1)&&(testTwo.y<=y9))
				{
					if((testTwo.y>=y1)&&(testTwo.y<=y2)){naA1.push_back(testTwo);}
					else if((testTwo.y>y2)&&(testTwo.y<=y3)){naA2.push_back(testTwo);}
					else if((testTwo.y>y3)&&(testTwo.y<=y4)){naA3.push_back(testTwo);}
					else if((testTwo.y>y4)&&(testTwo.y<=y5)){naA4.push_back(testTwo);}
					else if((testTwo.y>y5)&&(testTwo.y<=y6)){naA5.push_back(testTwo);}
					else if((testTwo.y>y6)&&(testTwo.y<=y7)){naA6.push_back(testTwo);}
					else if((testTwo.y>y7)&&(testTwo.y<=y8)){naA7.push_back(testTwo);}
					else if((testTwo.y>y8)&&(testTwo.y<=y9)){naA8.push_back(testTwo);}
				}
			
				long long testymax = test.y;
				long long testymin = testTwo.y;
				if(testymax<testymin)
				{
					long long temp11 = testymax;
					testymax = testymin;
					testymin = temp11;
				}
				
				Line3 lineNew;
				lineNew.p1 = test;
				lineNew.p2 = testTwo;
				if(doOverlapyyy(testymin,testymax,y1,y2)){ngA1.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y2,y3)){ngA2.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y3,y4)){ngA3.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y4,y5)){ngA4.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y5,y6)){ngA5.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y6,y7)){ngA6.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y7,y8)){ngA7.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y8,y9)){ngA8.push_back(lineNew);}
			} 
			
			for (int j = 0; j < L2PolVCount; j++)   
			{
				gpc_vertex2 test = hb1[prefix_overlay[l2PolyId] + j + 1];    
				gpc_vertex2 testTwo = hb1[prefix_overlay[l2PolyId] + j]; 
								
				if(((test.y<y1)&&(testTwo.y<y1))||((test.y>y9)&&(testTwo.y>y9)))
				continue;
			
				if((testTwo.y>=y1)&&(testTwo.y<=y9))
				{
					if((testTwo.y>=y1)&&(testTwo.y<=y2)){nbA1.push_back(testTwo);}
					else if((testTwo.y>y2)&&(testTwo.y<=y3)){nbA2.push_back(testTwo);}
					else if((testTwo.y>y3)&&(testTwo.y<=y4)){nbA3.push_back(testTwo);}
					else if((testTwo.y>y4)&&(testTwo.y<=y5)){nbA4.push_back(testTwo);}
					else if((testTwo.y>y5)&&(testTwo.y<=y6)){nbA5.push_back(testTwo);}
					else if((testTwo.y>y6)&&(testTwo.y<=y7)){nbA6.push_back(testTwo);}
					else if((testTwo.y>y7)&&(testTwo.y<=y8)){nbA7.push_back(testTwo);}
					else if((testTwo.y>y8)&&(testTwo.y<=y9)){nbA8.push_back(testTwo);}
				}
				
				long long testymax = test.y;
				long long testymin = testTwo.y;
				if(testymax<testymin)
				{
					long long temp22 = testymax;
					testymax = testymin;
					testymin = temp22;
				}
				
				Line3 lineNew;
				lineNew.p1 = test;
				lineNew.p2 = testTwo;
				if(doOverlapyyy(testymin,testymax,y1,y2)){nhA1.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y2,y3)){nhA2.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y3,y4)){nhA3.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y4,y5)){nhA4.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y5,y6)){nhA5.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y6,y7)){nhA6.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y7,y8)){nhA7.push_back(lineNew);}
				if(doOverlapyyy(testymin,testymax,y8,y9)){nhA8.push_back(lineNew);}
			
			} 
			ngSumAsubsub1.push_back(ngA1.size());
			nhSumAsubsub1.push_back(nhA1.size());
			ngSumAsubsub2.push_back(ngA2.size());
			nhSumAsubsub2.push_back(nhA2.size());
			ngSumAsubsub3.push_back(ngA3.size());
			nhSumAsubsub3.push_back(nhA3.size());
			ngSumAsubsub4.push_back(ngA4.size());
			nhSumAsubsub4.push_back(nhA4.size());
			
			ngSumAsubsub5.push_back(ngA5.size());
			nhSumAsubsub5.push_back(nhA5.size());
			ngSumAsubsub6.push_back(ngA6.size());
			nhSumAsubsub6.push_back(nhA6.size());
			ngSumAsubsub7.push_back(ngA7.size());
			nhSumAsubsub7.push_back(nhA7.size());
			ngSumAsubsub8.push_back(ngA8.size());
			nhSumAsubsub8.push_back(nhA8.size());
			
			naSumAsubsub1.push_back(naA1.size());
			nbSumAsubsub1.push_back(nbA1.size());
			naSumAsubsub2.push_back(naA2.size());
			nbSumAsubsub2.push_back(nbA2.size());
			naSumAsubsub3.push_back(naA3.size());
			nbSumAsubsub3.push_back(nbA3.size());
			naSumAsubsub4.push_back(naA4.size());
			nbSumAsubsub4.push_back(nbA4.size());
		
			naSumAsubsub5.push_back(naA5.size());
			nbSumAsubsub5.push_back(nbA5.size());
			naSumAsubsub6.push_back(naA6.size());
			nbSumAsubsub6.push_back(nbA6.size());
			naSumAsubsub7.push_back(naA7.size());
			nbSumAsubsub7.push_back(nbA7.size());
			naSumAsubsub8.push_back(naA8.size());
			nbSumAsubsub8.push_back(nbA8.size());
			
			realTaskNumsubsub2.push_back(originalid);
			
			#pragma omp critical
			{
				vecSegmentA11.insert(vecSegmentA11.end(),ngA1.begin(),ngA1.end());
				vecSegmentA21.insert(vecSegmentA21.end(),nhA1.begin(),nhA1.end());
				vecSegmentA12.insert(vecSegmentA12.end(),ngA2.begin(),ngA2.end());
				vecSegmentA22.insert(vecSegmentA22.end(),nhA2.begin(),nhA2.end());
				vecSegmentA13.insert(vecSegmentA13.end(),ngA3.begin(),ngA3.end());
				vecSegmentA23.insert(vecSegmentA23.end(),nhA3.begin(),nhA3.end());
				vecSegmentA14.insert(vecSegmentA14.end(),ngA4.begin(),ngA4.end());
				vecSegmentA24.insert(vecSegmentA24.end(),nhA4.begin(),nhA4.end());
				
				vecSegmentA15.insert(vecSegmentA15.end(),ngA5.begin(),ngA5.end());
				vecSegmentA25.insert(vecSegmentA25.end(),nhA5.begin(),nhA5.end());
				vecSegmentA16.insert(vecSegmentA16.end(),ngA6.begin(),ngA6.end());
				vecSegmentA26.insert(vecSegmentA26.end(),nhA6.begin(),nhA6.end());
				vecSegmentA17.insert(vecSegmentA17.end(),ngA7.begin(),ngA7.end());
				vecSegmentA27.insert(vecSegmentA27.end(),nhA7.begin(),nhA7.end());
				vecSegmentA18.insert(vecSegmentA18.end(),ngA8.begin(),ngA8.end());
				vecSegmentA28.insert(vecSegmentA28.end(),nhA8.begin(),nhA8.end());
				
				vecPointP11.insert(vecPointP11.end(),naA1.begin(),naA1.end());
				vecPointP21.insert(vecPointP21.end(),nbA1.begin(),nbA1.end());
				vecPointP12.insert(vecPointP12.end(),naA2.begin(),naA2.end());
				vecPointP22.insert(vecPointP22.end(),nbA2.begin(),nbA2.end());
				vecPointP13.insert(vecPointP13.end(),naA3.begin(),naA3.end());
				vecPointP23.insert(vecPointP23.end(),nbA3.begin(),nbA3.end());
				vecPointP14.insert(vecPointP14.end(),naA4.begin(),naA4.end());
				vecPointP24.insert(vecPointP24.end(),nbA4.begin(),nbA4.end());
				
				vecPointP15.insert(vecPointP15.end(),naA5.begin(),naA5.end());
				vecPointP25.insert(vecPointP25.end(),nbA5.begin(),nbA5.end());
				vecPointP16.insert(vecPointP16.end(),naA6.begin(),naA6.end());
				vecPointP26.insert(vecPointP26.end(),nbA6.begin(),nbA6.end());
				vecPointP17.insert(vecPointP17.end(),naA7.begin(),naA7.end());
				vecPointP27.insert(vecPointP27.end(),nbA7.begin(),nbA7.end());
				vecPointP18.insert(vecPointP18.end(),naA8.begin(),naA8.end());
				vecPointP28.insert(vecPointP28.end(),nbA8.begin(),nbA8.end());
				
				ngSumAsub1.insert(ngSumAsub1.end(),ngSumAsubsub1.begin(),ngSumAsubsub1.end());
				nhSumAsub1.insert(nhSumAsub1.end(),nhSumAsubsub1.begin(),nhSumAsubsub1.end());
				ngSumAsub2.insert(ngSumAsub2.end(),ngSumAsubsub2.begin(),ngSumAsubsub2.end());
				nhSumAsub2.insert(nhSumAsub2.end(),nhSumAsubsub2.begin(),nhSumAsubsub2.end());
				ngSumAsub3.insert(ngSumAsub3.end(),ngSumAsubsub3.begin(),ngSumAsubsub3.end());
				nhSumAsub3.insert(nhSumAsub3.end(),nhSumAsubsub3.begin(),nhSumAsubsub3.end());
				ngSumAsub4.insert(ngSumAsub4.end(),ngSumAsubsub4.begin(),ngSumAsubsub4.end());
				nhSumAsub4.insert(nhSumAsub4.end(),nhSumAsubsub4.begin(),nhSumAsubsub4.end());
				
				ngSumAsub5.insert(ngSumAsub5.end(),ngSumAsubsub5.begin(),ngSumAsubsub5.end());
				nhSumAsub5.insert(nhSumAsub5.end(),nhSumAsubsub5.begin(),nhSumAsubsub5.end());
				ngSumAsub6.insert(ngSumAsub6.end(),ngSumAsubsub6.begin(),ngSumAsubsub6.end());
				nhSumAsub6.insert(nhSumAsub6.end(),nhSumAsubsub6.begin(),nhSumAsubsub6.end());
				ngSumAsub7.insert(ngSumAsub7.end(),ngSumAsubsub7.begin(),ngSumAsubsub7.end());
				nhSumAsub7.insert(nhSumAsub7.end(),nhSumAsubsub7.begin(),nhSumAsubsub7.end());
				ngSumAsub8.insert(ngSumAsub8.end(),ngSumAsubsub8.begin(),ngSumAsubsub8.end());
				nhSumAsub8.insert(nhSumAsub8.end(),nhSumAsubsub8.begin(),nhSumAsubsub8.end());
				
				naSumAsub1.insert(naSumAsub1.end(),naSumAsubsub1.begin(),naSumAsubsub1.end());
				nbSumAsub1.insert(nbSumAsub1.end(),nbSumAsubsub1.begin(),nbSumAsubsub1.end());
				naSumAsub2.insert(naSumAsub2.end(),naSumAsubsub2.begin(),naSumAsubsub2.end());
				nbSumAsub2.insert(nbSumAsub2.end(),nbSumAsubsub2.begin(),nbSumAsubsub2.end());
				naSumAsub3.insert(naSumAsub3.end(),naSumAsubsub3.begin(),naSumAsubsub3.end());
				nbSumAsub3.insert(nbSumAsub3.end(),nbSumAsubsub3.begin(),nbSumAsubsub3.end());
				naSumAsub4.insert(naSumAsub4.end(),naSumAsubsub4.begin(),naSumAsubsub4.end());
				nbSumAsub4.insert(nbSumAsub4.end(),nbSumAsubsub4.begin(),nbSumAsubsub4.end());
				
				naSumAsub5.insert(naSumAsub5.end(),naSumAsubsub5.begin(),naSumAsubsub5.end());
				nbSumAsub5.insert(nbSumAsub5.end(),nbSumAsubsub5.begin(),nbSumAsubsub5.end());
				naSumAsub6.insert(naSumAsub6.end(),naSumAsubsub6.begin(),naSumAsubsub6.end());
				nbSumAsub6.insert(nbSumAsub6.end(),nbSumAsubsub6.begin(),nbSumAsubsub6.end());
				naSumAsub7.insert(naSumAsub7.end(),naSumAsubsub7.begin(),naSumAsubsub7.end());
				nbSumAsub7.insert(nbSumAsub7.end(),nbSumAsubsub7.begin(),nbSumAsubsub7.end());
				naSumAsub8.insert(naSumAsub8.end(),naSumAsubsub8.begin(),naSumAsubsub8.end());
				nbSumAsub8.insert(nbSumAsub8.end(),nbSumAsubsub8.begin(),nbSumAsubsub8.end());
				
				realTaskNumsub2.insert(realTaskNumsub2.end(),realTaskNumsubsub2.begin(),realTaskNumsubsub2.end());
			}
		}
	}
	Line3 *newAHc1 = vecSegmentA11.data();    
	Line3 *newAHd1 = vecSegmentA21.data();
	Line3 *newAHc2 = vecSegmentA12.data();    
	Line3 *newAHd2 = vecSegmentA22.data();
	Line3 *newAHc3 = vecSegmentA13.data();    
	Line3 *newAHd3 = vecSegmentA23.data();
	Line3 *newAHc4 = vecSegmentA14.data();    
	Line3 *newAHd4 = vecSegmentA24.data();
	
	Line3 *newAHc5 = vecSegmentA15.data();    
	Line3 *newAHd5 = vecSegmentA25.data();
	Line3 *newAHc6 = vecSegmentA16.data();    
	Line3 *newAHd6 = vecSegmentA26.data();
	Line3 *newAHc7 = vecSegmentA17.data();    
	Line3 *newAHd7 = vecSegmentA27.data();
	Line3 *newAHc8 = vecSegmentA18.data();    
	Line3 *newAHd8 = vecSegmentA28.data();
	
	gpc_vertex2 *newAHa1 = vecPointP11.data();    
	gpc_vertex2 *newAHb1 = vecPointP21.data();
	gpc_vertex2 *newAHa2 = vecPointP12.data();    
	gpc_vertex2 *newAHb2 = vecPointP22.data();
	gpc_vertex2 *newAHa3 = vecPointP13.data();    
	gpc_vertex2 *newAHb3 = vecPointP23.data();
	gpc_vertex2 *newAHa4 = vecPointP14.data();    
	gpc_vertex2 *newAHb4 = vecPointP24.data();
	
	gpc_vertex2 *newAHa5 = vecPointP15.data();    
	gpc_vertex2 *newAHb5 = vecPointP25.data();
	gpc_vertex2 *newAHa6 = vecPointP16.data();    
	gpc_vertex2 *newAHb6 = vecPointP26.data();
	gpc_vertex2 *newAHa7 = vecPointP17.data();    
	gpc_vertex2 *newAHb7 = vecPointP27.data();
	gpc_vertex2 *newAHa8 = vecPointP18.data();    
	gpc_vertex2 *newAHb8 = vecPointP28.data();
	
	int *ngSumA1 = ngSumAsub1.data();
	int *nhSumA1 = nhSumAsub1.data();
	int *ngSumA2 = ngSumAsub2.data();
	int *nhSumA2 = nhSumAsub2.data();
	int *ngSumA3 = ngSumAsub3.data();
	int *nhSumA3 = nhSumAsub3.data();
	int *ngSumA4 = ngSumAsub4.data();
	int *nhSumA4 = nhSumAsub4.data();
	
	int *ngSumA5 = ngSumAsub5.data();
	int *nhSumA5 = nhSumAsub5.data();
	int *ngSumA6 = ngSumAsub6.data();
	int *nhSumA6 = nhSumAsub6.data();
	int *ngSumA7 = ngSumAsub7.data();
	int *nhSumA7 = nhSumAsub7.data();
	int *ngSumA8 = ngSumAsub8.data();
	int *nhSumA8 = nhSumAsub8.data();
	
	int *naSumA1 = naSumAsub1.data();
	int *nbSumA1 = nbSumAsub1.data();
	int *naSumA2 = naSumAsub2.data();
	int *nbSumA2 = nbSumAsub2.data();
	int *naSumA3 = naSumAsub3.data();
	int *nbSumA3 = nbSumAsub3.data();
	int *naSumA4 = naSumAsub4.data();
	int *nbSumA4 = nbSumAsub4.data();
	
	int *naSumA5 = naSumAsub5.data();
	int *nbSumA5 = nbSumAsub5.data();
	int *naSumA6 = naSumAsub6.data();
	int *nbSumA6 = nbSumAsub6.data();
	int *naSumA7 = naSumAsub7.data();
	int *nbSumA7 = nbSumAsub7.data();
	int *naSumA8 = naSumAsub8.data();
	int *nbSumA8 = nbSumAsub8.data();
	
	int *realTaskNum2 = realTaskNumsub2.data();
	
	int *ngPrefixA1 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA1 = (int *)malloc(tasks * sizeof(int));
	int *ngPrefixA2 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA2 = (int *)malloc(tasks * sizeof(int));
	int *ngPrefixA3 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA3 = (int *)malloc(tasks * sizeof(int));
	int *ngPrefixA4 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA4 = (int *)malloc(tasks * sizeof(int));
	
	int *ngPrefixA5 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA5 = (int *)malloc(tasks * sizeof(int));
	int *ngPrefixA6 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA6 = (int *)malloc(tasks * sizeof(int));
	int *ngPrefixA7 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA7 = (int *)malloc(tasks * sizeof(int));
	int *ngPrefixA8 = (int *)malloc(tasks * sizeof(int));
	int *nhPrefixA8 = (int *)malloc(tasks * sizeof(int));
	
	int *naPrefixA1 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA1 = (int *)malloc(tasks * sizeof(int));
	int *naPrefixA2 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA2 = (int *)malloc(tasks * sizeof(int));
	int *naPrefixA3 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA3 = (int *)malloc(tasks * sizeof(int));
	int *naPrefixA4 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA4 = (int *)malloc(tasks * sizeof(int));
	
	int *naPrefixA5 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA5 = (int *)malloc(tasks * sizeof(int));
	int *naPrefixA6 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA6 = (int *)malloc(tasks * sizeof(int));
	int *naPrefixA7 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA7 = (int *)malloc(tasks * sizeof(int));
	int *naPrefixA8 = (int *)malloc(tasks * sizeof(int));
	int *nbPrefixA8 = (int *)malloc(tasks * sizeof(int));
	
	ngPrefixA1[0] = 0;
	nhPrefixA1[0] = 0;
	ngPrefixA2[0] = 0;
	nhPrefixA2[0] = 0;
	ngPrefixA3[0] = 0;
	nhPrefixA3[0] = 0;
	ngPrefixA4[0] = 0;
	nhPrefixA4[0] = 0;
	
	ngPrefixA5[0] = 0;
	nhPrefixA5[0] = 0;
	ngPrefixA6[0] = 0;
	nhPrefixA6[0] = 0;
	ngPrefixA7[0] = 0;
	nhPrefixA7[0] = 0;
	ngPrefixA8[0] = 0;
	nhPrefixA8[0] = 0;
	
	naPrefixA1[0] = 0;
	nbPrefixA1[0] = 0;
	naPrefixA2[0] = 0;
	nbPrefixA2[0] = 0;
	naPrefixA3[0] = 0;
	nbPrefixA3[0] = 0;
	naPrefixA4[0] = 0;
	nbPrefixA4[0] = 0;
	
	naPrefixA5[0] = 0;
	nbPrefixA5[0] = 0;
	naPrefixA6[0] = 0;
	nbPrefixA6[0] = 0;
	naPrefixA7[0] = 0;
	nbPrefixA7[0] = 0;
	naPrefixA8[0] = 0;
	nbPrefixA8[0] = 0;
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA1[idx]= ngPrefixA1[idx-1]+ngSumA1[idx-1];	
		nhPrefixA1[idx]= nhPrefixA1[idx-1]+nhSumA1[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA2[idx]= ngPrefixA2[idx-1]+ngSumA2[idx-1];	
		nhPrefixA2[idx]= nhPrefixA2[idx-1]+nhSumA2[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA3[idx]= ngPrefixA3[idx-1]+ngSumA3[idx-1];	
		nhPrefixA3[idx]= nhPrefixA3[idx-1]+nhSumA3[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA4[idx]= ngPrefixA4[idx-1]+ngSumA4[idx-1];	
		nhPrefixA4[idx]= nhPrefixA4[idx-1]+nhSumA4[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA5[idx]= ngPrefixA5[idx-1]+ngSumA5[idx-1];	
		nhPrefixA5[idx]= nhPrefixA5[idx-1]+nhSumA5[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA6[idx]= ngPrefixA6[idx-1]+ngSumA6[idx-1];	
		nhPrefixA6[idx]= nhPrefixA6[idx-1]+nhSumA6[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA7[idx]= ngPrefixA7[idx-1]+ngSumA7[idx-1];	
		nhPrefixA7[idx]= nhPrefixA7[idx-1]+nhSumA7[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		ngPrefixA8[idx]= ngPrefixA8[idx-1]+ngSumA8[idx-1];	
		nhPrefixA8[idx]= nhPrefixA8[idx-1]+nhSumA8[idx-1];
	}
	
	//for points
	for(int idx=1; idx<nNum;idx++){
		naPrefixA1[idx]= naPrefixA1[idx-1]+naSumA1[idx-1];	
		nbPrefixA1[idx]= nbPrefixA1[idx-1]+nbSumA1[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		naPrefixA2[idx]= naPrefixA2[idx-1]+naSumA2[idx-1];	
		nbPrefixA2[idx]= nbPrefixA2[idx-1]+nbSumA2[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		naPrefixA3[idx]= naPrefixA3[idx-1]+naSumA3[idx-1];	
		nbPrefixA3[idx]= nbPrefixA3[idx-1]+nbSumA3[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		naPrefixA4[idx]= naPrefixA4[idx-1]+naSumA4[idx-1];	
		nbPrefixA4[idx]= nbPrefixA4[idx-1]+nbSumA4[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		naPrefixA5[idx]= naPrefixA5[idx-1]+naSumA5[idx-1];	
		nbPrefixA5[idx]= nbPrefixA5[idx-1]+nbSumA5[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		naPrefixA6[idx]= naPrefixA6[idx-1]+naSumA6[idx-1];	
		nbPrefixA6[idx]= nbPrefixA6[idx-1]+nbSumA6[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		naPrefixA7[idx]= naPrefixA7[idx-1]+naSumA7[idx-1];	
		nbPrefixA7[idx]= nbPrefixA7[idx-1]+nbSumA7[idx-1];
	}
	
	for(int idx=1; idx<nNum;idx++){
		naPrefixA8[idx]= naPrefixA8[idx-1]+naSumA8[idx-1];	
		nbPrefixA8[idx]= nbPrefixA8[idx-1]+nbSumA8[idx-1];
	}
	
	endtime3 = my_difftime();
	difference3 = endtime3 - starttime3;
	printf("pre-process for PNP function =  %f\t \n",difference3);

	acc_init(acc_device_nvidia);
	int num_gpus = acc_get_num_devices(acc_device_nvidia); 
	//printf("This system has %d GPUs\n", num_gpus);
	int tid;
	#pragma omp parallel num_threads(num_gpus)
	{
		#pragma omp for private(tid)
		for (tid = 0; tid < num_gpus; tid++)
		{
			int threadNum = omp_get_thread_num();
			acc_set_device_num(threadNum, acc_device_nvidia);

			int gpu_num = acc_get_device_num(acc_device_nvidia); 
			//printf("Thread # %d is going to use GPU # %d \n", threadNum, gpu_num);
			//printf("tid:%d\n",tid);

			if(gpu_num==0)
			{
				wakeGPUup(tasks);
				//printf("0 done\n");
			}
			else if(gpu_num==1)
			{
				wakeGPUup(tasks);
				//printf("1 done\n");
			}
			else if(gpu_num==2)
			{
				wakeGPUup(tasks);
				//printf("2 done\n");
			}
			else if(gpu_num==3)
			{
				wakeGPUup(tasks);
				//printf("3 done\n");
			}
		}
	}
	
	double starttime4, endtime4;
    double difference4;	
    starttime4 = my_difftime();
	
	int realNumofInside = polygonInsideT+polygonInsideT2;
	int *fa1=(int *)malloc(nNum*sizeof(int));
	int *fa2=(int *)malloc(nNum*sizeof(int));
	int *fa3=(int *)malloc(nNum*sizeof(int));
	int *fa4=(int *)malloc(nNum*sizeof(int));
	int *fa5=(int *)malloc(nNum*sizeof(int));
	int *fa6=(int *)malloc(nNum*sizeof(int));
	int *fa7=(int *)malloc(nNum*sizeof(int));
	int *fa8=(int *)malloc(nNum*sizeof(int));
	int *fa9=(int *)malloc(realNumofInside*sizeof(int));
	#pragma omp parallel num_threads(num_gpus)
	{
		#pragma omp for private(tid)
		for (tid = 0; tid < num_gpus; tid++)
		{
			int threadNum = omp_get_thread_num();
			acc_set_device_num(threadNum, acc_device_nvidia);

			int gpu_num = acc_get_device_num(acc_device_nvidia); 
			//printf("Thread # %d is going to use GPU # %d \n", threadNum, gpu_num);
			//printf("tid:%d\n",tid);

			if(gpu_num==0)
			{
				fa2= processPNPTasks16(nNum, newAHc2, newAHd2, ngSumA2, nhSumA2, ngPrefixA2, nhPrefixA2,
						newAHa2, newAHb2, naSumA2, nbSumA2, naPrefixA2, nbPrefixA2);
				fa7= processPNPTasks16(nNum, newAHc7, newAHd7, ngSumA7, nhSumA7, ngPrefixA7, nhPrefixA7,
						newAHa7, newAHb7, naSumA7, nbSumA7, naPrefixA7, nbPrefixA7);
				fa9= processPNPTasks19(tasks, htaskSubId, htaskClipId, num_contours, numOfQuerys, num_base, num_overlay, prefix_base, prefix_overlay, ha1, hb1,
						newNum, nNum, refineOne, polygonInsideT,refineTwo,polygonInsideT2);
				
				//printf("0 done\n");
			}
			else if(gpu_num==1)
			{
				fa1= processPNPTasks16(nNum, newAHc1, newAHd1, ngSumA1, nhSumA1, ngPrefixA1, nhPrefixA1,
						newAHa1, newAHb1, naSumA1, nbSumA1, naPrefixA1, nbPrefixA1);
				fa6= processPNPTasks16(nNum, newAHc6, newAHd6, ngSumA6, nhSumA6, ngPrefixA6, nhPrefixA6,
						newAHa6, newAHb6, naSumA6, nbSumA6, naPrefixA6, nbPrefixA6);
				fa8= processPNPTasks16(nNum, newAHc8, newAHd8, ngSumA8, nhSumA8, ngPrefixA8, nhPrefixA8,
						newAHa8, newAHb8, naSumA8, nbSumA8, naPrefixA8, nbPrefixA8);
				//printf("1 done\n");
			}
			else if(gpu_num==2)
			{
				fa4= processPNPTasks16(nNum, newAHc4, newAHd4, ngSumA4, nhSumA4, ngPrefixA4, nhPrefixA4,
						newAHa4, newAHb4, naSumA4, nbSumA4, naPrefixA4, nbPrefixA4);
				fa5= processPNPTasks16(nNum, newAHc5, newAHd5, ngSumA5, nhSumA5, ngPrefixA5, nhPrefixA5,
						newAHa5, newAHb5, naSumA5, nbSumA5, naPrefixA5, nbPrefixA5);
				
				//printf("2 done\n");
			}
			else if(gpu_num==3)
			{
				fa3= processPNPTasks16(nNum, newAHc3, newAHd3, ngSumA3, nhSumA3, ngPrefixA3, nhPrefixA3,
						newAHa3, newAHb3, naSumA3, nbSumA3, naPrefixA3, nbPrefixA3);
				//printf("3 done\n");
			}	
		}
	}
	
	endtime4 = my_difftime();
	difference4 = endtime4 - starttime4;
	printf("Total PNP Function time taken =  %f\t \n",difference4);		
	printf("Total time taken =  %f\t \n",difference1+difference2+difference3+difference4);	 
	
	return 0;
}

//............................................................................................................................//
//for WKT
int localProcessing38(int num_contours, double *rect1, double *rect2, double *rect3,double *rect4, 
		int numOfQuerys,double *rect1_query,double *rect2_query, double *rect3_query,
		double *rect4_query,int *num_base,long *prefix_base,Coordinate *vertex_base,
		int *num_overlay,long *prefix_overlay,Coordinate *vertex_overlay)
{
    rtreeBuildingAndSearch2(num_contours, rect1, rect2, rect3, rect4,
     numOfQuerys, rect1_query, rect2_query, rect3_query, rect4_query);

	int tasks = 0; 
    int clipIndex;
    int a;
             
	for(a = 0; a < num_contours; a++)
	{
        tasks = tasks + subjectPoly2[a].num_elementsBasePoly;
    }
    
    int *htaskSubId = (int *)malloc(tasks * sizeof(int));
    int *htaskClipId = (int *)malloc(tasks * sizeof(int)); 
    int counter = 0;
	
    for(a = 0; a < num_contours; a++)
	{
        for(clipIndex = 0; clipIndex < subjectPoly2[a].num_elementsBasePoly; clipIndex++)
        {
          htaskSubId[counter] = a;
          htaskClipId[counter] = subjectPoly2[a].overlayPolyIndices[clipIndex];
          counter = counter + 1;
        }
    }
	printf("total tasks = %d \n", tasks);
	
	//transfer from double to long long 
	gpc_vertex2* ha1 = transfer_pointsvalues2(num_contours, num_base, prefix_base, vertex_base);
	gpc_vertex2* hb1 = transfer_pointsvalues2(numOfQuerys, num_overlay, prefix_overlay, vertex_overlay);
	long long* rect1n = transfer_boundingbox(num_contours, rect1);
	long long* rect2n = transfer_boundingbox(num_contours, rect2);
	long long* rect3n = transfer_boundingbox(num_contours, rect3);
	long long* rect4n = transfer_boundingbox(num_contours, rect4);
	long long* rect1_queryn = transfer_boundingbox(numOfQuerys, rect1_query);
	long long* rect2_queryn = transfer_boundingbox(numOfQuerys, rect2_query);
	long long* rect3_queryn = transfer_boundingbox(numOfQuerys, rect3_query);
	long long* rect4_queryn = transfer_boundingbox(numOfQuerys, rect4_query);
	
	double starttime1, endtime1;
    double difference1;	
    starttime1 = my_difftime();
	
	//get basic information for sketch
	int *numOfPartL1 = (int *)malloc(tasks * sizeof(int));   //how many parts of 1st polygon
	int *numOfPartL2 = (int *)malloc(tasks * sizeof(int)); 
	int *lastNumL1 = (int *)malloc(tasks * sizeof(int)); 	//for the final part, how many points are inside
	int *lastNumL2 = (int *)malloc(tasks * sizeof(int)); 	
	int *cellsizeL1 = (int *)malloc(tasks * sizeof(int));   
	int *cellsizeL2 = (int *)malloc(tasks * sizeof(int)); 
	
	//get sketch basic information
	#pragma omp parallel num_threads(32)
	{
		#pragma omp for schedule(static)
		for(int i =0; i<tasks;i++)
		{
			int l1PolyId = htaskSubId[i];		 
			int L1PolVCount = num_base[l1PolyId];  
			int l2PolyId = htaskClipId[i];		
			int L2PolVCount = num_overlay[l2PolyId];  
			
			int cellsize1 = 20;
			int cellsize2 = 20;
			if(L1PolVCount<400){cellsize1=10;}
			if(L2PolVCount<400){cellsize2=10;}
			
			cellsizeL1[i] = cellsize1;
			cellsizeL2[i] = cellsize2;
			
			int partL1 = ((L1PolVCount-1)/(cellsize1-1))+1;
			int lastNumL1sub = L1PolVCount - (cellsize1-1)*(partL1-1);
			if((partL1==0)||(partL1==1)){partL1=1;lastNumL1sub=L1PolVCount;}
			if(lastNumL1sub<=3)
			{
				partL1 = partL1-1;
				lastNumL1sub = lastNumL1sub+(cellsize1-1);
			}
			numOfPartL1[i] = partL1;
			lastNumL1[i] = lastNumL1sub;
		
			int partL2 = ((L2PolVCount-1)/(cellsize2-1))+1;
			int lastNumL2sub = L2PolVCount - (cellsize2-1)*(partL2-1);
			if((partL2==0)||(partL2==1)){partL2=1;lastNumL2sub=L2PolVCount;}
			if(lastNumL2sub<=3)
			{
				partL2 = partL2-1;
				lastNumL2sub = lastNumL2sub+(cellsize2-1);
			}
			numOfPartL2[i] = partL2;
			lastNumL2[i] = lastNumL2sub;
		}
	}

	long long *prefixPQ1 = (long long *)malloc(tasks * sizeof(long long)); 	
	long long *prefixPQ2 = (long long *)malloc(tasks * sizeof(long long)); 
	
	prefixPQ1[0] = 0;
	prefixPQ2[0] = 0;
	
	for(int i =1;i<tasks;i++)
	{	
		prefixPQ1[i] = prefixPQ1[i-1]+ numOfPartL1[i-1];
		prefixPQ2[i] = prefixPQ2[i-1]+ numOfPartL2[i-1];
	}
	
	endtime1 = my_difftime();
	difference1 = endtime1 - starttime1;
	printf("Get sketch basic information =  %f\t \n",difference1);
	
	acc_set_device_num(0, acc_device_nvidia);
	int gpu_num2 = acc_get_device_num(acc_device_nvidia); 
	//printf("IT is going to use GPU # %d \n",  gpu_num2);
	wakeGPUup(10204580);
	
	double starttime2, endtime2;
    double difference2;	
    starttime2 = my_difftime();
	
	//LSI function with PolySketch
	int *SIresult = calculateMbrPQ6(tasks, htaskSubId, htaskClipId, num_contours, numOfQuerys, num_base, num_overlay, 
						prefix_base, prefix_overlay, ha1, hb1, numOfPartL1, numOfPartL2, lastNumL1, lastNumL2, 
						prefixPQ1, prefixPQ2,cellsizeL1, cellsizeL2);
						
	endtime2 = my_difftime();
	difference2 = endtime2 - starttime2;
	printf("Segment Intersection Function time taken =  %f\t \n",difference2);
	
	return 0;
}

int readShapefile(int argc, char *argv[])
{    
	int i,j;

	SHPHandle	hSHP_base;
	hSHP_base = SHPOpen(argv[1],"rb");
	PolyRect *mbrs_base;
	
	int startIndex_base=0;
    int endIndex_base=hSHP_base->nRecords-1;
    int count_base=endIndex_base-startIndex_base+1;
    
	SHPReadMBR1(hSHP_base, startIndex_base, endIndex_base, &mbrs_base);
	double xmin_base=hSHP_base->adBoundsMin[0];
	double ymin_base=hSHP_base->adBoundsMin[1];
	double xmax_base= hSHP_base->adBoundsMax[0];
	double ymax_base=hSHP_base->adBoundsMax[1];
	
	SHPHandle	hSHP_overlay;
	hSHP_overlay = SHPOpen(argv[2],"rb");
	PolyRect *mbrs_overlay;
	int startIndex_overlay=0;
    int endIndex_overlay=hSHP_overlay->nRecords-1;
    int count_overlay=endIndex_overlay-startIndex_overlay+1;
	SHPReadMBR1(hSHP_overlay, startIndex_overlay, endIndex_overlay, &mbrs_overlay);

	double xmin_overlay=hSHP_overlay->adBoundsMin[0];
	double ymin_overlay=hSHP_overlay->adBoundsMin[1];
	double xmax_overlay= hSHP_overlay->adBoundsMax[0];
	double ymax_overlay=hSHP_overlay->adBoundsMax[1];
	
	double xmin=(xmin_base>xmin_overlay)?xmin_overlay:xmin_base;
	double ymin=(ymin_base>ymin_overlay)?ymin_overlay:ymin_base;
	double xmax=(xmax_base>xmax_overlay)?xmax_base:xmax_overlay;
	double ymax=(ymax_base>ymax_overlay)?ymax_base:ymax_overlay;
	
	bBox *baseBoxes;
	bBox *overlayBoxes;
	
	baseBoxes=(bBox *)malloc(sizeof(bBox));
	overlayBoxes=(bBox *)malloc(sizeof(bBox));
	baseBoxes->count=count_base;
	baseBoxes->rects=mbrs_base;
	overlayBoxes->count=count_overlay;
	overlayBoxes->rects=mbrs_overlay;	
	
	int sum_mbrs_base=baseBoxes[0].count ;
	int sum_mbrs_overlay=overlayBoxes[0].count;
	
	double minX_base,minY_base;
	minX_base=10000000;
	minY_base=10000000;
	SHPObject	**psShape_base;
	psShape_base=(SHPObject **)malloc((sum_mbrs_base)*sizeof(SHPObject*));
	
	int *num_base;
	int * prefix_base;
	
	int * hole_base_cpu=(int *)malloc((sum_mbrs_base)*sizeof(int));
	for(i=0;i<sum_mbrs_base;i++){
		hole_base_cpu[i]=0;
	}
	
	num_base=(int *)malloc((sum_mbrs_base)*sizeof(int));
	prefix_base=(int *)malloc((sum_mbrs_base)*sizeof(int));
	
	double 	*rect1=(double *)malloc(sum_mbrs_base*sizeof(double));
 	double	*rect2=(double *)malloc(sum_mbrs_base*sizeof(double));
 	double	*rect3=(double *)malloc(sum_mbrs_base*sizeof(double));
 	double	*rect4=(double *)malloc(sum_mbrs_base*sizeof(double));	
	
	convert(rect1,rect2,rect3,rect4,psShape_base,num_base,prefix_base,baseBoxes,0,
	sum_mbrs_base,hSHP_base,&minX_base,&minY_base);

	gpc_vertex *vertex_base_cpu=(gpc_vertex *)malloc((prefix_base[sum_mbrs_base-1]+num_base[sum_mbrs_base-1])
	*sizeof(gpc_vertex));
	for(i =0 ; i< sum_mbrs_base;i++)
	{
		for(j=0;j<num_base[i];j++)
		{
			vertex_base_cpu[prefix_base[i]+j].x=psShape_base[i]->padfX[j];
			vertex_base_cpu[prefix_base[i]+j].y=psShape_base[i]->padfY[j];
		}
	}
	
	destoryObjects(psShape_base,sum_mbrs_base);

	double minX_overlay,minY_overlay;
	minX_overlay=10000000;
	minY_overlay=10000000;   
	SHPObject	**psShape_overlay;
	psShape_overlay=(SHPObject **)malloc((sum_mbrs_overlay)*sizeof(SHPObject*));
	
	int *num_overlay;
	int * prefix_overlay;
	
	int * hole_overlay_cpu=(int *)malloc((sum_mbrs_overlay)*sizeof(int));
	for(i=0;i<overlayBoxes[0].count;i++){
		hole_overlay_cpu[i]=0;
	}
	
	num_overlay=(int *)malloc((sum_mbrs_overlay)*sizeof(int));
	prefix_overlay=(int *)malloc((sum_mbrs_overlay)*sizeof(int));
	double * rect1_query_cpu;
	double * rect2_query_cpu;
	double * rect3_query_cpu;
	double * rect4_query_cpu;
	rect1_query_cpu=(double *)malloc(sum_mbrs_overlay*sizeof(double));
 	rect2_query_cpu=(double *)malloc(sum_mbrs_overlay*sizeof(double));
 	rect3_query_cpu=(double *)malloc(sum_mbrs_overlay*sizeof(double));
 	rect4_query_cpu=(double *)malloc(sum_mbrs_overlay*sizeof(double));	
	convert(rect1_query_cpu,rect2_query_cpu,rect3_query_cpu,rect4_query_cpu,
        psShape_overlay,num_overlay,prefix_overlay,overlayBoxes,0,sum_mbrs_overlay,
        hSHP_overlay,&minX_overlay,&minY_overlay);
		
	gpc_vertex *vertex_overlay_cpu=(gpc_vertex *)malloc((prefix_overlay[sum_mbrs_overlay-1]+
	num_overlay[sum_mbrs_overlay-1])*sizeof(gpc_vertex));
	for(i =0 ; i< sum_mbrs_overlay;i++)
	{
		for(j=0;j<num_overlay[i];j++)
		{   
			vertex_overlay_cpu[prefix_overlay[i]+j].x=psShape_overlay[i]->padfX[j];
			vertex_overlay_cpu[prefix_overlay[i]+j].y=psShape_overlay[i]->padfY[j];
		}
	}
	destoryObjects(psShape_overlay,sum_mbrs_overlay);
			
	
	localProcessing35(sum_mbrs_base, rect1,rect2,rect3,rect4,NULL,sum_mbrs_overlay,
		rect1_query_cpu,rect2_query_cpu,rect3_query_cpu,rect4_query_cpu,NULL,num_base,
		prefix_base,vertex_base_cpu,num_overlay,prefix_overlay,vertex_overlay_cpu,
		hole_base_cpu,hole_overlay_cpu,6,12,0,minX_base,minY_base,NULL);
	//printf("final code with one and multi-GPUs end\n");
	//printf("  \n");
	
	SHPClose(hSHP_base);	
	SHPClose(hSHP_overlay);	
	
	return 0;
} 

int readNewFile(string filePath1, string filePath2){
	vector<string> *vstr1 = new vector<string>;
    vector<string> *vstr2 = new vector<string>;
    string tmpStr;

    ifstream file1(filePath1.c_str());
    while(std::getline(file1, tmpStr)){
        vstr1->push_back(tmpStr);
    }
    file1.close();

    ifstream file2(filePath2.c_str());
    while(std::getline(file2, tmpStr)){
        vstr2->push_back(tmpStr);
    }
    file2.close();

	polygonLayerTmp* tmpLayer1 = (polygonLayerTmp*)malloc(1* sizeof(polygonLayerTmp));
	tmpLayer1->polNum = 0;
    tmpLayer1->vNum = new vector<int>;
    tmpLayer1->prefixSum = new vector<long>;
    tmpLayer1->coord = new vector<Coordinate>;
    tmpLayer1->mbr = new vector<Coordinate>;
    tmpLayer1->mbr2 = new vector<Envelope>;
    
    polygonLayerTmp* tmpLayer2 = (polygonLayerTmp*)malloc(1* sizeof(polygonLayerTmp));
	tmpLayer2->polNum = 0;
    tmpLayer2->vNum = new vector<int>;
    tmpLayer2->prefixSum = new vector<long>;
    tmpLayer2->coord = new vector<Coordinate>;
    tmpLayer2->mbr = new vector<Coordinate>;
    tmpLayer2->mbr2 = new vector<Envelope>;
    
    int numThreads = 36;
    thread * parseThread = new thread[numThreads];

    int numOfStrPerThread1 = vstr1->size() / numThreads;
	int numOfStrPerThread2 = vstr2->size() / numThreads;
    
    for (int i = 0; i < numThreads; i++){
        int start1 = i * numOfStrPerThread1;
        int end1 = (i+1) * numOfStrPerThread1;
        int start2 = i * numOfStrPerThread2;
        int end2 = (i+1) * numOfStrPerThread2;
        if(numThreads - 1 == i){ 
            end1 = vstr1->size();
            end2 = vstr2->size();
        }
        parseThread[i] = thread(parsing, vstr1, start1, end1, vstr2, start2, end2, tmpLayer1, tmpLayer2);
    }

    for (int i = 0; i < numThreads; i++){
         parseThread[i].join();
     }
     delete [] parseThread;

	int polNum1 = tmpLayer1->polNum;
	int* vNum1 = tmpLayer1->vNum->data();
	long* prefixSum1 = tmpLayer1->prefixSum->data();		
	Coordinate* coord1 = tmpLayer1->coord->data();		
	vector<Coordinate>* mbr21 = tmpLayer1->mbr;		
	
	int polNum2 = tmpLayer2->polNum;
	int* vNum2 = tmpLayer2->vNum->data();
	long* prefixSum2 = tmpLayer2->prefixSum->data();		
	Coordinate* coord2 = tmpLayer2->coord->data();		
	vector<Coordinate> *mbr22 = tmpLayer2->mbr;			

	double *rect1 = (double*)malloc(polNum1*sizeof(double));
    double *rect2 = (double*)malloc(polNum1*sizeof(double));
    double *rect3 = (double*)malloc(polNum1*sizeof(double));
    double *rect4 = (double*)malloc(polNum1*sizeof(double));
    double *rect1_query = (double*)malloc(polNum2*sizeof(double));
    double *rect2_query = (double*)malloc(polNum2*sizeof(double));
    double *rect3_query = (double*)malloc(polNum2*sizeof(double));
    double *rect4_query = (double*)malloc(polNum2*sizeof(double));

	for(int idx=0;idx<(polNum1);idx++){
		rect1[idx]= mbr21->at(idx*2).x;
		rect2[idx]= mbr21->at(idx*2).y;
		rect3[idx]= mbr21->at(idx*2+1).x;
		rect4[idx]= mbr21->at(idx*2+1).y;
	}
	for(int idx=0;idx<(polNum2);idx++){
		rect1_query[idx]= mbr22->at(idx*2).x;
		rect2_query[idx]= mbr22->at(idx*2).y;
		rect3_query[idx]= mbr22->at(idx*2+1).x;
		rect4_query[idx]= mbr22->at(idx*2+1).y;
	}
	
	localProcessing38(polNum1, rect1,rect2,rect3,rect4,polNum2,
		rect1_query,rect2_query,rect3_query,rect4_query,vNum1,
		prefixSum1,coord1,vNum2,prefixSum2,coord2);
	printf("final code with one GPU end\n");
	printf("  \n");
	
	return 0;
}

//.....................................................................................................................................//
int main(int argc, char *argv[])
{ 		
	//for normal data
	//readShapefile(argc, argv);
	
	//for wkt data  
	string A1 = "/home/yliu0204/lakes_data";
	string A2 = "/home/yliu0204/sports_data";
	
	readNewFile(A1,A2);  
	
	return 0;
}

