/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkIntersectionPolyDataFilterMine.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkIntersectionPolyDataFilterMine.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDelaunay2D.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkObjectFactory.h"
#include "vtkOBBTree.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPointLocator.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkTransform.h"
#include "vtkTriangle.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkLongArray.h"
#include "vtkConnectivityFilter.h"
#include "vtkStripper.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkThreshold.h"
#include "vtkIdList.h"
#include "vtkCleanPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkMergeCells.h"
#include "vtkPolygon.h"
#include "vtkTransformPolyDataFilter.h"

#include <map>
#include <queue>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <list>
#include "math.h"

struct simPoint
{
  vtkIdType id;
  double pt[3];
  int concave;
};

struct simPolygon
{
  std::list<simPoint> points;
  int orientation;
};

//Function to turn an integer into a string
std::string int2String(int i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

double round_to_digits(double value, int digits)
{
  double factor = pow(10.0, digits - ceil(log10(fabs(value))));
  return round(value * factor) / factor;   
}

//Function to write the polydata to a vtp
void WriteVTPFile(vtkPolyData *writePolyData,std::string attachName)
{
  std::string outputFilename;

  vtkSmartPointer<vtkXMLPolyDataWriter> writer  = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  
  outputFilename = "/Users/adamupdegrove/Documents/Intersect/Tests/"+attachName+".vtp";

  writer->SetFileName(outputFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(writePolyData);
#else
  writer->SetInputData(writePolyData);
#endif
  //writer->SetDataModeToAscii();

  writer->Write();
}

//----------------------------------------------------------------------------
// Helper typedefs and data structure.
typedef std::multimap< vtkIdType, vtkIdType >    IntersectionMapType;
typedef IntersectionMapType::iterator            IntersectionMapIteratorType;

//typedef std::pair< vtkIdType, vtkIdType >            CellEdgePairType;
typedef struct _CellEdgeLine {
  vtkIdType CellId;
  vtkIdType EdgeId;
  vtkIdType LineId;
} CellEdgeLineType;

typedef std::multimap< vtkIdType, CellEdgeLineType > PointEdgeMapType;
typedef PointEdgeMapType::iterator                   PointEdgeMapIteratorType;


//----------------------------------------------------------------------------
// Private implementation to hide STL.
//----------------------------------------------------------------------------
class vtkIntersectionPolyDataFilterMine::Impl
{
public:
  Impl();
  virtual ~Impl();

  static int FindTriangleIntersections(vtkOBBNode *node0, vtkOBBNode *node1,
                                       vtkMatrix4x4 *transform, void *arg);

  int SplitMesh(int inputIndex, vtkPolyData *output,
                vtkPolyData *intersectionLines);

  void CleanAndCheckSurface(vtkPolyData *pd);
  
protected:

  vtkCellArray* SplitCell(vtkPolyData *input, vtkIdType cellId,
		  		     vtkIdType *cellPts, IntersectionMapType *map,
				     vtkPolyData *interLines, int inputIndex,
				     int numCurrCells);

  void GetLoops(vtkPolyData *pd,std::vector<simPolygon> *loops);

  void Triangulate(vtkPolyData *pd,simPolygon *loop,vtkCellArray *polys);

  int Clip(double newp1[2], double newp2[2], simPolygon *loop,
      std::list<simPoint>::iterator it1, std::list<simPoint>::iterator it2);

  int TestConcave(double line1[3],double line2[3],int orientation);

  void GetSingleLoop(vtkPolyData *pd,simPolygon *loop,vtkIdType nextCell,
      bool *interPtBool,bool *lineBool);

  void SetLoopOrientation(vtkPolyData *pd,simPolygon *loop,vtkIdType *nextCell,
    vtkIdType nextPt, vtkIdType prevPt, vtkIdList *pointCells);

  void FollowLoopOrientation(vtkPolyData *pd,simPolygon *loop,vtkIdType *nextCell,
    vtkIdType nextPt, vtkIdType prevPt, vtkIdList *pointCells);

  void Orient(vtkPolyData *pd,vtkTransform *transform, vtkPolyData *boundary,
		  vtkPolygon *boundarypoly);

  int GetLoopOrientation(vtkPolyData *pd,vtkIdType cell,vtkIdType ptId1,
      vtkIdType ptId2);

  int CheckLine(vtkPolyData *pd, vtkIdType ptId1,vtkIdType ptId2);

  int AddToPointEdgeMap(int index, vtkIdType ptId, double x[3],
                         vtkPolyData *mesh, vtkIdType cellId,
                         vtkIdType edgeId, vtkIdType lineId,
                         vtkIdType triPts[3]);

  void AddToNewCellMap(int inputIndex, int interPtCount, int interPts[3],
      vtkPolyData *interLines,int numCurrCells);

public:
  vtkPolyData         *Mesh[2];
  vtkOBBTree          *OBBTree1;

  // Stores the intersection lines.
  vtkCellArray        *IntersectionLines;

  vtkIdTypeArray      *SurfaceId;
  vtkIdTypeArray      *CaseId;
  vtkIdTypeArray      *NewCellIds[2];

  // Cell data that indicates in which cell each intersection
  // lies. One array for each output surface.
  vtkIdTypeArray      *CellIds[2];

  // Cell data that indicates on which surface the intersection point lies.

  // Map from points to the cells that contain them. Used for point
  // data interpolation. For points on the edge between two cells, it
  // does not matter which cell is recorded bcause the interpolation
  // will be the same.  One array for each output surface.
  vtkIdTypeArray      *PointCellIds[2];
  vtkIntArray      *BoundaryPoints[2];

  // Merging filter used to convert intersection lines from "line
  // soup" to connected polylines.
  vtkPointLocator     *PointMerger;

  // Map from cell ID to intersection line.
  IntersectionMapType *IntersectionMap[2];
  IntersectionMapType *IntersectionPtsMap[2];
  IntersectionMapType *PointMapper;

  // Map from point to an edge on which it resides, the ID of the
  // cell, and the ID of the line.
  PointEdgeMapType    *PointEdgeMap[2];

protected:
  Impl(const Impl&); // purposely not implemented
  void operator=(const Impl&); // purposely not implemented

};

//----------------------------------------------------------------------------
vtkIntersectionPolyDataFilterMine::Impl::Impl() :
  OBBTree1(0), IntersectionLines(0), PointMerger(0), SurfaceId(0), CaseId(0)
{
  for (int i = 0; i < 2; i++)
    {
    this->Mesh[i]                 = NULL;
    this->CellIds[i]              = NULL;
    this->IntersectionMap[i]      = new IntersectionMapType();
    this->IntersectionPtsMap[i]   = new IntersectionMapType();
    this->PointEdgeMap[i]         = new PointEdgeMapType();
    }
    this->PointMapper             = new IntersectionMapType();
}

//----------------------------------------------------------------------------
vtkIntersectionPolyDataFilterMine::Impl::~Impl()
{
  for (int i = 0; i < 2; i++)
    {
    delete this->IntersectionMap[i];
    delete this->IntersectionPtsMap[i];
    delete this->PointEdgeMap[i];
    }
    delete this->PointMapper;
}


//----------------------------------------------------------------------------
int vtkIntersectionPolyDataFilterMine::Impl
::FindTriangleIntersections(vtkOBBNode *node0, vtkOBBNode *node1,
                            vtkMatrix4x4 *transform, void *arg)
{
  vtkIntersectionPolyDataFilterMine::Impl *info =
    reinterpret_cast<vtkIntersectionPolyDataFilterMine::Impl*>(arg);

  vtkPolyData     *mesh0                 = info->Mesh[0];
  vtkPolyData     *mesh1                 = info->Mesh[1];
  vtkOBBTree      *obbTree1              = info->OBBTree1;
  vtkCellArray    *intersectionLines     = info->IntersectionLines;
  vtkIdTypeArray  *intersectionSurfaceId = info->SurfaceId;
  vtkIdTypeArray  *intersectionCaseId    = info->CaseId;
  vtkIdTypeArray  *intersectionCellIds0  = info->CellIds[0];
  vtkIdTypeArray  *intersectionCellIds1  = info->CellIds[1];
  vtkPointLocator *pointMerger           = info->PointMerger;

  static int total0 = 0;
  static int total1 = 0;
  static int cpt = 0;
  static int mypt = 0;
  static int upt = 0;

  int numCells0 = node0->Cells->GetNumberOfIds();
  total0 += numCells0;
  //std::cout<<"Number o cells: "<<numCells0<<endl;
  //std::cout<<"Total0: "<<total0<<endl;
  int retval = 0;

  for (vtkIdType id0 = 0; id0 < numCells0; id0++)
    {
    vtkIdType cellId0 = node0->Cells->GetId(id0);
    int type0 = mesh0->GetCellType(cellId0);

    if (type0 == VTK_TRIANGLE)
      {
      vtkIdType npts0, *triPtIds0;
      mesh0->GetCellPoints(cellId0, npts0, triPtIds0);
      double triPts0[3][3];
      for (vtkIdType id = 0; id < npts0; id++)
        {
        mesh0->GetPoint(triPtIds0[id], triPts0[id]);
        }

      if (obbTree1->TriangleIntersectsNode
          (node1, triPts0[0], triPts0[1], triPts0[2], transform))
        {
        int numCells1 = node1->Cells->GetNumberOfIds();
	total1 += numCells1;
        //std::cout<<"Total1: "<<total1<<endl;
        for (vtkIdType id1 = 0; id1 < numCells1; id1++)
          {
          vtkIdType cellId1 = node1->Cells->GetId(id1);
          int type1 = mesh1->GetCellType(cellId1);
          if (type1 == VTK_TRIANGLE)
            {
            // See if the two cells actually intersect. If they do,
            // add an entry into the intersection maps and add an
            // intersection line.
            vtkIdType npts1, *triPtIds1;
            mesh1->GetCellPoints(cellId1, npts1, triPtIds1);

            double triPts1[3][3];
            for (vtkIdType id = 0; id < npts1; id++)
              {
              mesh1->GetPoint(triPtIds1[id], triPts1[id]);
              }

            int coplanar = 0;
	    int pointCase = 0;
            double outpt0[3], outpt1[3];
	    double surfaceid[2];
	    if (cellId0 == 41 && cellId1 == 1)
	      std::cout<<"Tell me whats happenin!"<<endl;
            int intersects =
              vtkIntersectionPolyDataFilterMine::TriangleTriangleIntersection
              (triPts0[0], triPts0[1], triPts0[2],
               triPts1[0], triPts1[1], triPts1[2],
               coplanar, outpt0, outpt1, surfaceid, pointCase);
	    std::cout<<"Intersects "<<intersects<<endl;
	    std::cout<<"The points are 0: " <<outpt0[0]<<" "<<outpt0[1]<<" "<<outpt0[2]<<endl;
	    std::cout<<"1: " <<outpt1[0]<<" "<<outpt1[1]<<" "<<outpt1[2]<<endl;

            if ( coplanar )
              {
              // Coplanar triangle intersection is not handled.
              // This intersection will not be included in the output. TODO
		std::cout<<"Coplanar"<<endl;
              intersects = 0;
              continue;
              }

            if ( intersects)// &&
                // ( outpt0[0] != outpt1[0] ||
                //   outpt0[1] != outpt1[1] ||
                //   outpt0[2] != outpt1[2] ) )
              {
              vtkIdType lineId = intersectionLines->GetNumberOfCells(); 

	      if (lineId == 300 || lineId == 301)
	      {
		std::cout<<"TELL MEEEEEEEEEE"<<endl;
	      }
              vtkIdType ptId0, ptId1;
	      int unique1,unique2;
	      vtkIdType check1,check2;
	      check1 = pointMerger->IsInsertedPoint(outpt0);
	      check2 = pointMerger->IsInsertedPoint(outpt1);
              unique1 = pointMerger->InsertUniquePoint(outpt0, ptId0);
              unique2 = pointMerger->InsertUniquePoint(outpt1, ptId1);

	      if (check1 == -1) 
	      {
		cpt++;
	        std::cout<<"Point 1 inserted is yes"<<endl;
	      }
	      if (check2 == -1)
	      {
		cpt++;
	        std::cout<<"Point 2 inserted is yes"<<endl;
	      }
	      if (unique1 == 1 && unique2 == 1) 
	      {
		upt++;
		std::cout<<"Point 1 inserted is yes"<<endl;
	      }
	      else if (unique1 == 1 || unique2 == 1)
	      {
		upt++;
	      }

	      int addline = 1;
	      if (ptId0 == ptId1)
	        addline = 0;

	      if (ptId0 == ptId1 && surfaceid[0] != surfaceid[1])
	      {
		std::cout<<"Both Points are the same! "<<ptId0<<endl;
		std::cout<<"TELL ME IMPORTANT "<<surfaceid[0]<<" "<<surfaceid[1]<<endl;
		intersectionSurfaceId->InsertValue(ptId0,3);
	      }
	      else
	      {
		if (unique1)
		{
		  intersectionSurfaceId->InsertValue(ptId0,surfaceid[0]);
		  std::cout<<"Surface Id for pt "<<ptId0<<" is "<<surfaceid[0]<<endl;
		}
		else
		{
	          if (intersectionSurfaceId->GetValue(ptId0) != 3)
		  {
		    intersectionSurfaceId->InsertValue(ptId0,surfaceid[0]);
		    std::cout<<"Surface Id for pt "<<ptId0<<" is "<<surfaceid[0]<<endl;
		  }
		}
		if (unique2)
		{
		  intersectionSurfaceId->InsertValue(ptId1,surfaceid[1]);
		  std::cout<<"Surface Id for pt "<<ptId1<<" is "<<surfaceid[1]<<endl;
		}
		else
		{
		  if (intersectionSurfaceId->GetValue(ptId1) != 3)
		  {
		    intersectionSurfaceId->InsertValue(ptId1,surfaceid[1]);
		    std::cout<<"Surface Id for pt "<<ptId1<<" is "<<surfaceid[1]<<endl;
		  }
		}
	      }
	      info->IntersectionPtsMap[0]->insert(std::make_pair(ptId0, cellId0));
	      info->IntersectionPtsMap[1]->insert(std::make_pair(ptId0, cellId1));
	      info->IntersectionPtsMap[0]->insert(std::make_pair(ptId1, cellId0));
	      info->IntersectionPtsMap[1]->insert(std::make_pair(ptId1, cellId1));
              
	      //Check to see if duplicate line. Line can only be a duplicate
	      //line if both points are not unique and they don't equal eachother
	      if (!unique1 && !unique2 && ptId0 != ptId1)
	      {
		vtkSmartPointer<vtkPolyData> test = vtkSmartPointer<vtkPolyData>::New();
		test->SetPoints(pointMerger->GetPoints());
		test->SetLines(intersectionLines);
		test->BuildLinks();
		std::cout<<"Both Points are not unique"<<endl;
		vtkSmartPointer<vtkIdList> iCells1 = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> iCells2 = vtkSmartPointer<vtkIdList>::New();
		test->GetPointCells(ptId0,iCells1);
		test->GetPointCells(ptId1,iCells2);
		iCells1->IntersectWith(iCells2);
		if (iCells1->GetNumberOfIds() > 0)
		{
		  if (iCells1->GetNumberOfIds() > 1)
		    std::cout<<"More than one cell herere!!"<<endl;
	          addline = 0;
		}
	      }

	      if (addline)
	      {
		intersectionLines->InsertNextCell(2);
		intersectionLines->InsertCellPoint(ptId0);
		intersectionLines->InsertCellPoint(ptId1);
		intersectionCaseId->InsertNextValue(pointCase);

		intersectionCellIds0->InsertNextValue(cellId0);
		intersectionCellIds1->InsertNextValue(cellId1);

		info->PointCellIds[0]->InsertValue( ptId0, cellId0 );
		info->PointCellIds[0]->InsertValue( ptId1, cellId0 );
		info->PointCellIds[1]->InsertValue( ptId0, cellId1 );
		info->PointCellIds[1]->InsertValue( ptId1, cellId1 );

		info->IntersectionMap[0]->insert(std::make_pair(cellId0, lineId));
		info->IntersectionMap[1]->insert(std::make_pair(cellId1, lineId));

		// Check which edges of cellId0 and cellId1 outpt0 and
		// outpt1 are on, if any.
		int isOnEdge=0;
		int m0p0=0,m0p1=0,m1p0=0,m1p1=0;
		for (vtkIdType edgeId = 0; edgeId < 3; edgeId++)
		  {
		    std::cout<<"Edge "<<edgeId<<endl;
		  isOnEdge = info->AddToPointEdgeMap(0, ptId0, outpt0, mesh0, cellId0,
					  edgeId, lineId, triPtIds0);
		  if (isOnEdge != -1)
                    m0p0++;
		  std::cout<<"pt 0 Is on edge of mesh0? "<<isOnEdge<<endl;
		  isOnEdge = info->AddToPointEdgeMap(0, ptId1, outpt1, mesh0, cellId0,
					  edgeId, lineId, triPtIds0);
		  if (isOnEdge != -1)
                    m0p1++;
		  std::cout<<"pt 1 Is on edge of mesh0? "<<isOnEdge<<endl;
		  isOnEdge = info->AddToPointEdgeMap(1, ptId0, outpt0, mesh1, cellId1,
					  edgeId, lineId, triPtIds1);
		  if (isOnEdge != -1)
                    m1p0++;
		  std::cout<<"pt 0 Is on edge of mesh1? "<<isOnEdge<<endl;
		  isOnEdge = info->AddToPointEdgeMap(1, ptId1, outpt1, mesh1, cellId1,
					  edgeId, lineId, triPtIds1);
		  if (isOnEdge != -1)
                    m1p1++;
		  std::cout<<"pt 1 Is on edge of mesh1? "<<isOnEdge<<endl;
		  }

		if (m0p0 > 0 && m1p0 > 0)
		{
		  std::cout<<"Special CASE! check!"<<endl;
		  intersectionSurfaceId->InsertValue(ptId0,3);
		}
		if (m0p1 > 0 && m1p1 > 0)
		{
		  std::cout<<"Special CASE! check!"<<endl;
		  intersectionSurfaceId->InsertValue(ptId1,3);
		}
                  
	      }
	      else
	      {
		std::cout<<"I am not adding line with points "<<ptId0<<" and "<<ptId1<<endl;
	      }
	      retval++;
              }
	      else 
	      {
		vtkIdType check1,check2;
		check1 = pointMerger->IsInsertedPoint(outpt0);
		check2 = pointMerger->IsInsertedPoint(outpt1);
		std::cout<<"Points would have been! "<<check1<<" and "<<check2<<endl;
	      }
            }
          }
        }
      }
    }
  std::cout<<"Number unique pts "<<upt<<endl;
  std::cout<<"Number check pts "<<cpt<<endl;
  std::cout<<"Number my pts "<<mypt<<endl;

  return retval;
}

void vtkIntersectionPolyDataFilterMine::Impl
::CleanAndCheckSurface(vtkPolyData *pd)
{
  vtkIdType npts,p0,p1;
  vtkIdType *pts;
  int badcell = 0;
  int freeedgecell = 0;
  vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  vtkSmartPointer<vtkIntArray> bad = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> freeedge = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIdList> edgeneigh = vtkSmartPointer<vtkIdList>::New();

  //Clean the input surface
  cleaner->SetInputData(pd);
  cleaner->Update();
  pd->DeepCopy(cleaner->GetOutput());
  pd->BuildLinks();

  //Loop through the surface and find edges with cells that have either more 
  //than one neighbor or no neighbors. No neighbors can be okay,as this can 
  //indicate a free edge. However, for a polydata surface, multiple neighbors 
  //indicates a bad cell with possible intersecting facets!
  for (int i = 0;i<pd->GetNumberOfCells();i++)
  {
    pd->GetCellPoints(i,npts,pts);
    badcell = 0;
    freeedgecell = 0;
    for (int j=0;j<npts;j++)
    {
      p0 = pts[j];
      p1 = pts[(j+1)%npts];

      pd->GetCellEdgeNeighbors(i,p0,p1,edgeneigh);
      if (edgeneigh->GetNumberOfIds() > 1)
      {
	badcell = badcell + 1;
	std::cout<<"Cell on first surface has more than one!!: Numba "<<edgeneigh->GetNumberOfIds()<<" for cell "<<i<<endl;
      }
      else if (edgeneigh->GetNumberOfIds() < 1)
      {
	freeedgecell = freeedgecell + 1;
	std::cout<<"Cell on first surface has a free edge for cell "<<i<<endl;
      }

    }
    bad->InsertValue(i,badcell);
    freeedge->InsertValue(i,freeedgecell);
  }

  bad->SetName("BadTri");
  pd->GetCellData()->AddArray(bad);
  pd->GetCellData()->SetActiveScalars("BadTri");

  freeedge->SetName("FreeEdge");
  pd->GetCellData()->AddArray(freeedge);
  pd->GetCellData()->SetActiveScalars("FreeEdge");
}

//----------------------------------------------------------------------------
int vtkIntersectionPolyDataFilterMine::Impl
::SplitMesh(int inputIndex, vtkPolyData *output, vtkPolyData *intersectionLines)
{
  vtkSmartPointer<vtkIntArray> TouchedCellArray = 
    vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> TouchedPointArray = 
    vtkSmartPointer<vtkIntArray>::New();
  TouchedCellArray->SetNumberOfComponents(1);
  TouchedPointArray->SetNumberOfComponents(1);

  vtkPolyData *input = this->Mesh[inputIndex];
  IntersectionMapType *intersectionMap = this->IntersectionMap[inputIndex];
  vtkCellData *inCD  = input->GetCellData();
  vtkCellData *outCD = output->GetCellData();
  vtkIdType numCells = input->GetNumberOfCells();
  vtkIdType cellIdX = 0;

  //
  // Process points
  //
  vtkIdType inputNumPoints = input->GetPoints()->GetNumberOfPoints();
  vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
  points->Allocate(100);
  output->SetPoints(points);

  //
  // Split intersection lines. The lines structure is constructed
  // using a vtkPointLocator. However, some lines may have an endpoint
  // on a cell edge that has no neighbor. We need to duplicate a line
  // point in such a case and update the point ID in the line cell.
  //
  vtkSmartPointer< vtkPolyData > splitLines =
    vtkSmartPointer <vtkPolyData >::New();
  splitLines->DeepCopy(intersectionLines);

  vtkPointData *inPD  = input->GetPointData();
  vtkPointData *outPD = output->GetPointData();
  outPD->CopyAllocate( inPD, input->GetNumberOfPoints() );

  // Copy over the point data from the input
  for (vtkIdType ptId = 0; ptId < inputNumPoints; ptId++)
    {
    double pt[3];
    input->GetPoints()->GetPoint(ptId, pt);
    output->GetPoints()->InsertNextPoint(pt);
    outPD->CopyData(inPD, ptId, ptId);
    TouchedPointArray->InsertTuple1(ptId,-2);
    BoundaryPoints[inputIndex]->InsertValue(ptId,0);
    }

  // Copy the points from splitLines to the output, interpolating the
  // data as we go.
  std::cout<<"NumberOfPoints "<<splitLines->GetNumberOfPoints()<<endl;
  std::cout<<"NumberOfLines "<<splitLines->GetNumberOfLines()<<endl;
  for (vtkIdType id = 0; id < splitLines->GetNumberOfPoints(); id++)
    {
    double pt[3];
    splitLines->GetPoint(id, pt);
    vtkIdType newPtId = output->GetPoints()->InsertNextPoint(pt);

    // Retrieve the cell ID from splitLines
    vtkIdType cellId = this->PointCellIds[inputIndex]->GetValue( id );

    double closestPt[3], pcoords[3], dist2, weights[3];
    int subId;
    vtkCell *cell = input->GetCell(cellId);
    cell->EvaluatePosition(pt, closestPt, subId, pcoords, dist2, weights);
    outPD->InterpolatePoint(input->GetPointData(), newPtId, cell->PointIds,
                            weights);
    TouchedPointArray->InsertTuple1(newPtId,-2);
    BoundaryPoints[inputIndex]->InsertValue(newPtId,0);
    }

  //
  // Process cells
  //
  outCD->CopyAllocate(inCD, numCells);

  if ( input->GetPolys()->GetNumberOfCells() > 0 )
    {
    vtkCellArray *cells = input->GetPolys();
    vtkIdType newId = output->GetNumberOfCells();

    vtkSmartPointer< vtkCellArray > newPolys =
      vtkSmartPointer< vtkCellArray >::New();

    newPolys->EstimateSize(cells->GetNumberOfCells(),3);
    output->SetPolys(newPolys);

    vtkSmartPointer< vtkIdList > edgeNeighbors =
      vtkSmartPointer< vtkIdList >::New();
    vtkIdType nptsX = 0;
    vtkIdType *pts = 0;
    vtkSmartPointer< vtkIdList > cellsToCheck =
      vtkSmartPointer< vtkIdList >::New();
    for (cells->InitTraversal(); cells->GetNextCell(nptsX, pts); cellIdX++)
      {
      if ( nptsX != 3 )
        {
        vtkGenericWarningMacro( << "vtkIntersectionPolyDataFilterMine only works with "
                                << "triangle meshes." );
        continue;
        }

      cellsToCheck->Reset();
      cellsToCheck->Allocate(nptsX+1);
      cellsToCheck->InsertNextId(cellIdX);

      // Collect the cells relevant for splitting this cell.  If the
      // cell is in the intersection map, split. If not, one of its
      // edges may be split by an intersection line that splits a
      // neighbor cell. Mark the cell as needing a split if this is
      // the case.
      bool needsSplit = intersectionMap->find( cellIdX ) != intersectionMap->end();
      for (vtkIdType ptId = 0; ptId < nptsX; ptId++)
        {
        vtkIdType pt0Id = pts[ptId];
        vtkIdType pt1Id = pts[(ptId+1) % nptsX];
        edgeNeighbors->Reset();
        input->GetCellEdgeNeighbors(cellIdX, pt0Id, pt1Id, edgeNeighbors);
        for (vtkIdType nbr = 0; nbr < edgeNeighbors->GetNumberOfIds(); nbr++)
          {
          vtkIdType nbrId = edgeNeighbors->GetId(nbr);
          cellsToCheck->InsertNextId(nbrId);

          if ( intersectionMap->find( nbrId ) != intersectionMap->end() )
            {
            needsSplit = true;
            }
          } // for (vtkIdType nbr = 0; ...
        } // for (vtkIdType pt = 0; ...

      // Splitting occurs here
      if ( !needsSplit )
        {
        // Just insert the cell and copy the cell data
        newId = newPolys->InsertNextCell(3, pts);
        outCD->CopyData(inCD, cellIdX, newId);
	TouchedCellArray->InsertTuple1(newId,-2);
        //BoundaryPoints[inputIndex]->InsertValue(newId,0);
	
        }
      else
        {
  	int numCurrCells = newPolys->GetNumberOfCells();
        vtkCellArray *splitCells = this->SplitCell
          (input, cellIdX, pts, intersectionMap, splitLines, inputIndex,numCurrCells);
        

        double pt0[3], pt1[3], pt2[3], normal[3];
        points->GetPoint(pts[0], pt0);
        points->GetPoint(pts[1], pt1);
        points->GetPoint(pts[2], pt2);
        vtkTriangle::ComputeNormal(pt0, pt1, pt2, normal);
        vtkMath::Normalize(normal);

        vtkIdType npts, *ptIds, subCellId = 0;
        for (splitCells->InitTraversal(); splitCells->GetNextCell(npts, ptIds); subCellId++)
          {
          // Check for reversed cells. I'm not sure why, but in some
          // cases, cells are reversed.
          double subCellNormal[3];
          points->GetPoint(ptIds[0], pt0);
          points->GetPoint(ptIds[1], pt1);
          points->GetPoint(ptIds[2], pt2);
          vtkTriangle::ComputeNormal(pt0, pt1, pt2, subCellNormal);
          vtkMath::Normalize(subCellNormal);

          if ( vtkMath::Dot(normal, subCellNormal) > 0 )
            {
            newId = newPolys->InsertNextCell(npts, ptIds);
            }
          else
            {
            newId = newPolys->InsertNextCell(npts);
            for ( int i = 0; i < npts; i++)
              {
              newPolys->InsertCellPoint( ptIds[ npts-i-1 ] );
              }
            }

          outCD->CopyData(inCD, cellIdX, newId); // Duplicate cell data
	  TouchedCellArray->InsertTuple1(newId,1);
	  for (int i=0;i<npts;i++)
	  {
	    TouchedPointArray->InsertTuple1(ptIds[i],1);
	    std::cout<<"Touched Point Adding Point!: "<<ptIds[i]<<endl;
	  }

          }

        splitCells->Delete();
        }
      } // for (cells->InitTraversal(); ...
    }
  TouchedCellArray->SetName("VisitedCells");
  outCD->AddArray(TouchedCellArray);
  outCD->SetActiveScalars("VisitedCells");

  TouchedPointArray->SetName("VisitedPoints");
  outPD->AddArray(TouchedPointArray);
  outPD->SetActiveScalars("VisitedPoints");

  return 1;
}

vtkCellArray* vtkIntersectionPolyDataFilterMine::Impl
::SplitCell(vtkPolyData *input, vtkIdType cellId, vtkIdType *cellPts,
            IntersectionMapType *map,
	    vtkPolyData *interLines, int inputIndex,
	    int numCurrCells)
{
  // Copy down the SurfaceID array that tells which surface the point belongs
  // to
  vtkSmartPointer<vtkLongArray> surfaceMapper = 
    vtkSmartPointer<vtkLongArray>::New();
  surfaceMapper = static_cast<vtkLongArray*>(interLines->GetPointData()->GetArray("SurfaceID"));

  //Array to keep track of which points are on the boundary of the cell
  vtkSmartPointer<vtkIdTypeArray> cellBoundaryPt = 
    vtkSmartPointer<vtkIdTypeArray>::New();
  //Array to tell wheter the original cell points lie on the intersecting 
  //line
  int CellPointOnInterLine[3] = {0,0,0};

  // Gather points from the cell
  vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
  vtkSmartPointer< vtkPointLocator > merger =
    vtkSmartPointer< vtkPointLocator >::New();
  merger->SetTolerance(1e-6);
  merger->InitPointInsertion( points, input->GetBounds() );

  double xyz[3];
  for ( int i = 0; i < 3; i++)
    {
    if ( cellPts[i] >= input->GetNumberOfPoints() )
      {
      vtkGenericWarningMacro(<< "invalid point read 1");
      }
    input->GetPoint(cellPts[i], xyz);
    merger->InsertNextPoint( xyz );
    cellBoundaryPt->InsertNextValue(1);
    }

  // Set up line cells and array to track the just the intersecting lines 
  // on the cell.
  vtkSmartPointer< vtkCellArray > lines =
    vtkSmartPointer< vtkCellArray >::New();
  vtkSmartPointer< vtkCellArray > interceptlines =
    vtkSmartPointer< vtkCellArray >::New();

  double p0[3], p1[3], p2[3];
  input->GetPoint(cellPts[0], p0);
  input->GetPoint(cellPts[1], p1);
  input->GetPoint(cellPts[2], p2);

  // This maps the point IDs for the vtkPolyData passed to
  // vtkDelaunay2D back to the original IDs in interLines. NOTE: The
  // point IDs from the cell are not stored here.
  std::map< vtkIdType, vtkIdType > ptIdMap;

  IntersectionMapIteratorType iterLower = map->lower_bound( cellId );
  IntersectionMapIteratorType iterUpper = map->upper_bound( cellId );
  //Get all the lines associated with the original cell
  while ( iterLower != iterUpper )
  {
    vtkIdType lineId = iterLower->second;
    vtkIdType nLinePts, *linePtIds;
    interLines->GetLines()->GetCell( 3*lineId, nLinePts, linePtIds );

    interceptlines->InsertNextCell(2);
    lines->InsertNextCell(2);
    //Loop through the points of each line
    for (vtkIdType i = 0; i < nLinePts; i++)
    {
      std::map< vtkIdType, vtkIdType >::iterator location =
	ptIdMap.find( linePtIds[i] );
      //If point already isn't in list
      if ( location == ptIdMap.end() )
      {
	interLines->GetPoint( linePtIds[i], xyz );
	if ( linePtIds[i] >= interLines->GetNumberOfPoints() )
	{
	  vtkGenericWarningMacro(<< "invalid point read 2");
	}
	//Check to see if point is unique
	int unique = merger->InsertUniquePoint( xyz, ptIdMap[ linePtIds[i] ]);
	if (unique)
	{
          //If point is unique, check to see if it is actually a
	  //point originating from this input surface or on both surfaces
	  //Don't mark as boundary point if it originates from other surface
	  if (surfaceMapper->GetValue(linePtIds[i]) == inputIndex + 1 ||
	      surfaceMapper->GetValue(linePtIds[i]) == 3)
	    cellBoundaryPt->InsertValue(ptIdMap[linePtIds[i]],1);
	  else
	    cellBoundaryPt->InsertValue(ptIdMap[linePtIds[i]],0);
	}
	else 
	{
	  //Obviously if the pointid is less than three, it is one of the 
	  //original cell points and can be added to the inter cell point arr
	  if (ptIdMap[linePtIds[i] ] < 3)
	    CellPointOnInterLine[ptIdMap[linePtIds[i]]] = 1;
	}
	interceptlines->InsertCellPoint( ptIdMap[ linePtIds[i] ] );
	std::cout<<"First Adding Line "<<lines->GetNumberOfCells()<<" with point "<<ptIdMap[ linePtIds[i] ]<<endl;
	lines->InsertCellPoint( ptIdMap[ linePtIds[i] ] );
      }
      //Point is already in list, so run through checks with its value
      else
      {
	interceptlines->InsertCellPoint( location->second );
	std::cout<<"First Adding Line "<<lines->GetNumberOfCells()<<" with point "<<ptIdMap[ linePtIds[i] ]<<endl;
	lines->InsertCellPoint( location->second );
	if (location->second < 3)
	  CellPointOnInterLine[location->second] = 1;
       }
    }
    ++iterLower;
  }
  
  // Now check the neighbors of the cell
  IntersectionMapIteratorType ptIterLower;  
  IntersectionMapIteratorType ptIterUpper;
  IntersectionMapIteratorType cellIterLower;
  IntersectionMapIteratorType cellIterUpper;
  vtkSmartPointer< vtkIdList > nbrCellIds =
    vtkSmartPointer< vtkIdList >::New();
  for (vtkIdType i = 0; i < 3; i++)
  {
    //Get Points belonging to each edge of this cell
    vtkIdType edgePtId0 = cellPts[i];
    vtkIdType edgePtId1 = cellPts[(i+1) % 3];

    double edgePt0[3], edgePt1[3];
    if ( edgePtId0 >= input->GetNumberOfPoints() )
    {
      vtkGenericWarningMacro( << "invalid point read 3");
    }
    if ( edgePtId1 >= input->GetNumberOfPoints() )
    {
      vtkGenericWarningMacro( << "invalid point read 4");
    }
    input->GetPoint(edgePtId0, edgePt0);
    input->GetPoint(edgePtId1, edgePt1);

    nbrCellIds->Reset();
    input->GetCellEdgeNeighbors(cellId, edgePtId0, edgePtId1, nbrCellIds);
    //Loop through attached neighbor cells and check for split edges
    for (vtkIdType j = 0; j < nbrCellIds->GetNumberOfIds(); j++)
    {
      vtkIdType nbrCellId = nbrCellIds->GetId( j );
      iterLower = map->lower_bound( nbrCellId );
      iterUpper = map->upper_bound( nbrCellId );
      while ( iterLower != iterUpper )
      {
        vtkIdType lineId = iterLower->second;
        vtkIdType nLinePts, *linePtIds;
        interLines->GetLines()->GetCell( 3*lineId, nLinePts, linePtIds );
        for (vtkIdType k = 0; k < nLinePts; k++)
        {
          if ( linePtIds[k] >= interLines->GetNumberOfPoints() )
          {
            vtkGenericWarningMacro( << "invalid point read 5");
          }
          interLines->GetPoint( linePtIds[k], xyz );
	  ptIterLower = this->PointMapper->lower_bound(linePtIds[k]);
	  ptIterUpper = this->PointMapper->upper_bound(linePtIds[k]);
	  //Find all points withing this neighbor cell
	  while (ptIterLower != ptIterUpper)
	  {
	    vtkIdType mappedPtId = ptIterLower->second;
	    std::cout<<"Mapped Pt Id: "<<mappedPtId<<endl;
	    cellIterLower = this->IntersectionPtsMap[inputIndex]->lower_bound(mappedPtId);
	    cellIterUpper = this->IntersectionPtsMap[inputIndex]->upper_bound(mappedPtId);
	    //Check all cell values associated with this point
	    while (cellIterLower != cellIterUpper)
	    {
	      vtkIdType checkCellId = cellIterLower->second;
	      std::cout<<"Attached Cell Id: "<<checkCellId<<endl;

	      //If this cell id is the same as the current cell id, this 
	      //means the point is a split edge, need to add to list!!
	      if (checkCellId == cellId)
	      {
	        std::cout<<"Adding because cellId is: "<<cellId<<endl;
	        int unique=0;
	        if ( ptIdMap.find( linePtIds[k] ) == ptIdMap.end() )
	        {
	          unique = merger->InsertUniquePoint( xyz, ptIdMap[ linePtIds[k] ] );
	        }
	        else
	        {
		  //Point is less than 3, original cell point
		  if (ptIdMap[ linePtIds[k] ] < 3)
		    CellPointOnInterLine[ptIdMap[ linePtIds[k] ]] = 1;
	        }
	        if (unique)
	        {
		  //Check to see what surface point originates from. Don't
		  //mark if point is from other surface
		  if (surfaceMapper->GetValue(linePtIds[k]) == inputIndex + 1 ||
		      surfaceMapper->GetValue(linePtIds[k]) == 3)
		    cellBoundaryPt->InsertValue(ptIdMap[linePtIds[k]],1);
		  else
		  {
		  cellBoundaryPt->InsertValue(ptIdMap[linePtIds[k]],0);
		  }

	        }
		else
		{
		  if (ptIdMap[ linePtIds[k] ] < 3)
		    CellPointOnInterLine[ptIdMap[ linePtIds[k] ]] = 1;
		}
	      }
	      ++cellIterLower;
	    }
	    ++ptIterLower;
	  }
	}
        ++iterLower;
      }
    }
  }

  // Set up reverse ID map
  std::map< vtkIdType, vtkIdType > reverseIdMap;
  std::map< vtkIdType, vtkIdType > reverseLineIdMap;
  std::map< vtkIdType, vtkIdType >::iterator iter = ptIdMap.begin();
  while ( iter != ptIdMap.end() )
    {
    // If we have more than one point mapping back to the same point
    // in the input mesh, just use the first one. This will give a
    // preference for using cell points when an intersection line shares
    // a point with a a cell and prevent introducing accidental holes
    // in the mesh.
    if ( reverseIdMap.find( iter->second ) == reverseIdMap.end() )
      {
      reverseIdMap[ iter->second ] = iter->first + input->GetNumberOfPoints();
      }
    if ( reverseLineIdMap.find( iter->second ) == reverseLineIdMap.end() )
      {
      reverseLineIdMap[ iter->second ] = iter->first;
      }
    ++iter;
    }
  std::cout<<"Reverse Id Map:"<<endl;
  for (std::map<vtkIdType,vtkIdType>::iterator mapit=reverseIdMap.begin(); mapit!=reverseIdMap.end(); ++mapit)
        std::cout << mapit->first << " => " << mapit->second << '\n';

  std::cout<<"Reverse Line Id Map:"<<endl;
  for (std::map<vtkIdType,vtkIdType>::iterator mapit=reverseLineIdMap.begin(); mapit!=reverseLineIdMap.end(); ++mapit)
        std::cout << mapit->first << " => " << mapit->second << '\n';

  double v0[3], v1[3], n[3], c[3];
  vtkTriangle::TriangleCenter( p0, p1, p2, c );
  vtkTriangle::ComputeNormal( p0, p1, p2, n );
  vtkMath::Perpendiculars( n, v0, v1, 0.0 );

  // For each point on an edge, compute it's relative angle about n.
  vtkSmartPointer< vtkIdTypeArray > edgePtIdList =
    vtkSmartPointer< vtkIdTypeArray >::New();
  vtkSmartPointer< vtkIdTypeArray > interPtIdList =
    vtkSmartPointer< vtkIdTypeArray >::New();
  edgePtIdList->Allocate( points->GetNumberOfPoints() );
  vtkSmartPointer< vtkDoubleArray > angleList =
    vtkSmartPointer< vtkDoubleArray >::New();
  angleList->Allocate( points->GetNumberOfPoints() );
  bool *interPtBool;
  interPtBool = new bool[points->GetNumberOfPoints()];

  std::cout<<"Number Of POINTS!!! "<<points->GetNumberOfPoints()<<endl;
  for ( vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ptId++ )
    {
    double x[3], closestPt[3], t0, t1, t2;
    points->GetPoint( ptId, x );
    
    interPtBool[ptId] = false;
    if (cellBoundaryPt->GetValue(ptId)) 
      { 
//      std::cout<<"Point has been added: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
      // Point is on line. Add its id to id list and add its angle to
      // angle list.
      edgePtIdList->InsertNextValue( ptId );
      double d[3];
      vtkMath::Subtract( x, c, d );
      angleList->InsertNextValue( atan2( vtkMath::Dot(d, v0),
                                         vtkMath::Dot(d, v1) ) );
	if (ptId > 2)
	{
//	  std::cout<<"Intersecting Point!!"<<endl;
	  interPtIdList->InsertNextValue( ptId );
	  interPtBool[ptId] = true;
	}
      }
    if (ptId > 2)
    {
	    std::cout<<"Setting Boundary Point "<<reverseIdMap[ptId]<<endl;
          BoundaryPoints[inputIndex]->InsertValue(reverseIdMap[ptId],1);
    }
    else if (CellPointOnInterLine[ptId])
    {
	    std::cout<<"Setting Boundary Point "<<cellPts[ptId]<<endl;
          BoundaryPoints[inputIndex]->InsertValue(cellPts[ptId],1);
    }
    else
    {
	    std::cout<<"Setting Boundary Point Bad "<<cellPts[ptId]<<endl;
          BoundaryPoints[inputIndex]->InsertValue(cellPts[ptId],0);
    }

    }
  // Sort the edgePtIdList according to the angle list. The starting
  // point doesn't matter. We just need to generate boundary lines in
  // a consistent order.
  vtkSortDataArray::Sort( angleList, edgePtIdList );

  vtkSmartPointer<vtkPolyData> checkPD = 
    vtkSmartPointer<vtkPolyData>::New();
  checkPD->SetPoints(points);
  checkPD->SetLines(lines);
  checkPD->BuildLinks();
  vtkIdType id;
  for ( id = 0; id < edgePtIdList->GetNumberOfTuples()-1; id++ )
    {
      int unique = this->CheckLine(checkPD,edgePtIdList->GetValue(id),
	  edgePtIdList->GetValue(id+1));
      if (unique)
      {
	std::cout<<"Second Adding Line "<<lines->GetNumberOfCells()<<" with points "<<edgePtIdList->GetValue( id )<<" and "<<edgePtIdList->GetValue( id + 1)<<endl;
	lines->InsertNextCell(2);
	lines->InsertCellPoint( edgePtIdList->GetValue( id ) );
	lines->InsertCellPoint( edgePtIdList->GetValue( id + 1 ) );
      }
    }
  int unique = this->CheckLine(checkPD, 
      edgePtIdList->GetValue(edgePtIdList->GetNumberOfTuples()-1), 
      edgePtIdList->GetValue(0));
  if (unique)
  {
    std::cout<<"Second Adding Line "<<lines->GetNumberOfCells()<<" with points "<<edgePtIdList->GetValue( id )<<" and "<<edgePtIdList->GetValue( 0)<<endl;
    lines->InsertNextCell(2);
    lines->InsertCellPoint
      ( edgePtIdList->GetValue( edgePtIdList->GetNumberOfTuples()-1 ) );
    lines->InsertCellPoint( edgePtIdList->GetValue( 0 ) );
  }
    
  // Set up a transform that will rotate the points to the
  // XY-plane (normal aligned with z-axis).
  vtkSmartPointer< vtkTransform > transform =
    vtkSmartPointer< vtkTransform >::New();
  double zaxis[3] = {0, 0, 1};
  double rotationAxis[3], normal[3], center[3], rotationAngle;

  double pt0[3], pt1[3], pt2[3];
  points->GetPoint( 0, pt0 );
  points->GetPoint( 1, pt1 );
  points->GetPoint( 2, pt2 );
  vtkTriangle::ComputeNormal(pt0, pt1, pt2, normal);

  double dotZAxis = vtkMath::Dot( normal, zaxis );
  if ( fabs(1.0 - dotZAxis) < 1e-6 )
    {
    // Aligned with z-axis
    rotationAxis[0] = 1.0;
    rotationAxis[1] = 0.0;
    rotationAxis[2] = 0.0;
    rotationAngle = 0.0;
    }
  else if ( fabs( 1.0 + dotZAxis ) < 1e-6 )
    {
    // Co-linear with z-axis, but reversed sense.
    // Aligned with z-axis
    rotationAxis[0] = 1.0;
    rotationAxis[1] = 0.0;
    rotationAxis[2] = 0.0;
    rotationAngle = 180.0;
    }
  else
    {
    // The general case
    vtkMath::Cross(normal, zaxis, rotationAxis);
    vtkMath::Normalize(rotationAxis);
    rotationAngle =
      vtkMath::DegreesFromRadians(acos(vtkMath::Dot(zaxis, normal)));
    }

  transform->PreMultiply();
  transform->Identity();

  transform->RotateWXYZ(rotationAngle,
			rotationAxis[0],
			rotationAxis[1],
			rotationAxis[2]);

  vtkTriangle::TriangleCenter(pt0, pt1, pt2, center);
  transform->Translate(-center[0], -center[1], -center[2]);

    vtkCellArray *splitCells = vtkCellArray::New();
#ifdef NEW_SPLIT
    vtkSmartPointer<vtkPolyData> interpd = vtkSmartPointer<vtkPolyData>::New();
    interpd->SetPoints(points);
    interpd->SetLines(interceptlines);
    interpd->BuildLinks();

    std::cout<<"NUMBER OF INTERSECTING POINTS! "<<interPtIdList->GetNumberOfTuples()<<endl;
    std::cout<<"NUMBER OF INTERSECTING LINES! "<<interceptlines->GetNumberOfCells()<<endl;
    std::cout<<"FULLCELLID: "<<cellId<<endl;
    vtkSmartPointer<vtkPolyData> fullpd = vtkSmartPointer<vtkPolyData>::New();
    fullpd->SetPoints(points);
    fullpd->SetLines(lines);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformer = 
      vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    vtkSmartPointer<vtkPolyData> transformedpd = 
      vtkSmartPointer<vtkPolyData>::New();
    transformer->SetInputData(fullpd);
    transformer->SetTransform(transform);
    transformer->Update();
    transformedpd = transformer->GetOutput();
    transformedpd->BuildLinks();
    //WriteVTPFile(transformedpd,"TransformedCell_"+int2String(cellId));

    //WriteVTPFile(fullpd,"FullCell"+int2String(cellId));
    //WriteVTPFile(interpd,"InterceptCell"+int2String(cellId));

    if (interPtIdList->GetNumberOfTuples()%2 != 0)
    {
      std::cout<<"There is a big problem here because triangle is not completely intersected"<<endl;
    }

    if (interPtIdList->GetNumberOfTuples() > 0 && interceptlines->GetNumberOfCells() > 0)
    {
      std::vector<simPolygon> loops;
      std::cout<<"BOOL: ";
      for (int k=0;k<points->GetNumberOfPoints();k++)
      {
	std::cout<<interPtBool[k]<<" ";
      }
      std::cout<<" end"<<endl;
      this->GetLoops(transformedpd,&loops);
      for (int k=0;k<loops.size();k++)
      {
        vtkCellArray *polys;
#ifdef EAR_CLIP
	polys = vtkCellArray::New();
	loops[k].points.pop_back();
	this->Triangulate(transformedpd,&loops[k],polys);

        vtkIdType npts, *ptIds;
	interLines->BuildLinks();
        for (polys->InitTraversal(); polys->GetNextCell(npts, ptIds); )
        {
          if ( ptIds[0] >= points->GetNumberOfPoints() ||
               ptIds[1] >= points->GetNumberOfPoints() ||
               ptIds[2] >= points->GetNumberOfPoints() )
          {
            vtkGenericWarningMacro( << "Invalid point ID!!!");
          }

          splitCells->InsertNextCell( npts );
	  int interPtCount=0;
	  int interPts[3];
          for (int i = 0; i < npts; i++)
          {
            vtkIdType remappedPtId;
            if ( ptIds[i] < 3) // Point from the cell
            {
              remappedPtId = cellPts[ ptIds[i]] ;
	      std::cout<<"Cell has point "<<remappedPtId<<endl;
	      if (CellPointOnInterLine[ptIds[i]])
	      {
		interPts[interPtCount++] = reverseLineIdMap[ptIds[i]];
	      }
            }
            else
            {
              remappedPtId = reverseIdMap[ptIds[i]];
	      std::cout<<"Cell has point "<<remappedPtId<<endl;
	      interPts[interPtCount++] = reverseLineIdMap[ptIds[i]];
            }
            splitCells->InsertCellPoint( remappedPtId );
          }
          std::cout<<"Inter Point Count! "<<interPtCount<<endl;
	  if (interPtCount >= 2)
	  {
	    this->AddToNewCellMap(inputIndex, interPtCount,interPts,
		interLines,numCurrCells);
	  }
          std::cout<<"Num Current Cells "<<numCurrCells<<endl;
	  numCurrCells++;
        }
	vtkSmartPointer<vtkPolyData> justcheck = 
	  vtkSmartPointer<vtkPolyData>::New();
	justcheck->SetPoints(points);
	justcheck->SetPolys(polys);
          //WriteVTPFile(justcheck,"NewLoop"+int2String(cellId)+"_"+int2String(k));
	  polys->Delete();
#else      
        vtkSmartPointer<vtkPolyData> newpd = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();
	std::list<simPoint>::iterator it;
	vtkIdType ptiter=0;
	int *pointMapper;
	pointMapper = new int[loops[k].points.size()];
	for (it = loops[k].points.begin();it != loops[k].points.end();++it)
	{
	  if (ptiter < loops[k].points.size()-1)
	  {
	    newPoints->InsertNextPoint(points->GetPoint((it)->id));
	    pointMapper[ptiter] = (it)->id;
	  }
	  if (ptiter < loops[k].points.size()-2)
	  {
	    newLines->InsertNextCell(2);
	    newLines->InsertCellPoint(ptiter);
	    newLines->InsertCellPoint(ptiter+1);
	  }
	  ptiter++;
	}
	newLines->InsertNextCell(2);
	newLines->InsertCellPoint(ptiter-2);
	newLines->InsertCellPoint(0);

	newpd->SetPoints(newPoints);
	newpd->SetLines(newLines);
	vtkSmartPointer<vtkPolyData> boundary = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolygon> boundaryPoly = vtkSmartPointer<vtkPolygon>::New();
	this->Orient(newpd,transform,boundary,boundaryPoly);
	if (cellId == 371 && inputIndex == 1)
	{
          WriteVTPFile(newpd,"PreDelaunay"+int2String(cellId)+"_"+int2String(k));
        //  //WriteVTPFile(boundary,"PreDelaunayBoundary"+int2String(cellId)+"_"+int2String(k));
	}

	vtkSmartPointer< vtkDelaunay2D > del2D =
	  vtkSmartPointer< vtkDelaunay2D >::New();
	del2D->SetInputData(newpd);
	del2D->SetSourceData(boundary);
	del2D->SetTolerance(0.0);
	del2D->SetAlpha(0.0);
	del2D->SetOffset(0);
        del2D->SetProjectionPlaneMode(VTK_SET_TRANSFORM_PLANE);
        del2D->SetTransform(transform);
	del2D->BoundingTriangulationOff();
	del2D->Update();
	if (cellId == 371 && inputIndex == 1)
	{
	  WriteVTPFile(del2D->GetOutput(),"NewLoopDelaunay"+int2String(cellId)+"_"+int2String(k));
	}

	vtkSmartPointer<vtkTriangleFilter> triangulator = 
          vtkSmartPointer<vtkTriangleFilter>::New();
	vtkSmartPointer< vtkDelaunay2D > del2Dxy =
	  vtkSmartPointer< vtkDelaunay2D >::New();
	if (del2D->GetOutput()->GetPolys()->GetNumberOfCells() != newpd->GetNumberOfPoints() - 2)
	{
	  vtkSmartPointer<vtkTransformPolyDataFilter> deltransformer = 
	    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	  vtkSmartPointer<vtkPoints> transpoints = 
	    vtkSmartPointer<vtkPoints>::New();
	  vtkSmartPointer<vtkPolyData> transpd = 
	    vtkSmartPointer<vtkPolyData>::New();
	  std::cout<<"Del2D failed, transform to xy plane, try again"<<endl;
	  deltransformer->SetInputData(newpd);
	  deltransformer->SetTransform(transform);
	  deltransformer->Update();

	  std::cout<<"Rotation angle "<<rotationAngle<<endl;
	  double pt[3];
	  for (vtkIdType i=0;i<deltransformer->GetOutput()->GetNumberOfPoints();i++)
	  {
	    deltransformer->GetOutput()->GetPoint(i,pt);
	    for (int j=0;j<2;j++)
	    {
	      if (j==0)
	        pt[j] = round_to_digits(pt[j],2);
	      else
	      if (j!=0)
	        pt[j] = round_to_digits(pt[j],3);
	    }
	    pt[2] = 0.;
	    transpoints->InsertNextPoint(pt);
	  }
	  transpd->SetPoints(transpoints);
	  transpd->SetLines(newLines);

	  del2Dxy->SetInputData(transpd);
	  del2Dxy->SetSourceData(transpd);
	  del2Dxy->SetTolerance(0.0);
	  del2Dxy->SetAlpha(0.0);
	  del2Dxy->SetOffset(0);
	  del2Dxy->BoundingTriangulationOff();
	  del2Dxy->Update();

          polys = del2Dxy->GetOutput()->GetPolys();
	  if (polys->GetNumberOfCells() != newpd->GetNumberOfPoints() - 2)
	  {
	    std::cout<<"Transform also failed, attempting ear clipping"<<endl;
	    triangulator->SetInputData(boundary);
	    triangulator->Update();
            polys = triangulator->GetOutput()->GetPolys();
	    if (polys->GetNumberOfCells() != newpd->GetNumberOfPoints() - 2)
	    {
	      std::cout<<"Clipping also failed :("<<endl;
	    }
	  }
	  if (cellId == 371 && inputIndex == 1)
	  {
	    WriteVTPFile(triangulator->GetOutput(),"NewLoopEarClip"+int2String(cellId)+"_"+int2String(k));
	  }
	}
	else
	{
          polys = del2D->GetOutput()->GetPolys();
	}
	

       // // Renumber the point IDs.
        vtkIdType npts, *ptIds;
	interLines->BuildLinks();
        for (polys->InitTraversal(); polys->GetNextCell(npts, ptIds); )
        {
          if ( pointMapper[ptIds[0]] >= points->GetNumberOfPoints() ||
               pointMapper[ptIds[1]] >= points->GetNumberOfPoints() ||
               pointMapper[ptIds[2]] >= points->GetNumberOfPoints() )
          {
            vtkGenericWarningMacro( << "Invalid point ID!!!");
          }

          splitCells->InsertNextCell( npts );
	  int interPtCount=0;
	  int interPts[3];
          for (int i = 0; i < npts; i++)
          {
            vtkIdType remappedPtId;
            if ( pointMapper[ptIds[i]] < 3) // Point from the cell
            {
              remappedPtId = cellPts[ pointMapper[ptIds[i]] ];
	      std::cout<<"Cell has point "<<remappedPtId<<endl;
	      if (CellPointOnInterLine[pointMapper[ptIds[i]]])
	      {
		interPts[interPtCount++] = reverseLineIdMap[pointMapper[ptIds[i]] ];
	      }
            }
            else
            {
              remappedPtId = reverseIdMap[ pointMapper[ptIds[i]] ];
	      std::cout<<"Cell has point "<<remappedPtId<<endl;
	      interPts[interPtCount++] = reverseLineIdMap[ pointMapper[ptIds[i]] ];
            }
            splitCells->InsertCellPoint( remappedPtId );
          }
          std::cout<<"Inter Point Count! "<<interPtCount<<endl;
	  if (interPtCount >= 2)
	  {
	    this->AddToNewCellMap(inputIndex, interPtCount,interPts,
		interLines,numCurrCells);
	  }
          std::cout<<"Num Current Cells "<<numCurrCells<<endl;
	  numCurrCells++;
        }
	delete [] pointMapper;
#endif
      }
    }
    else
    {
      vtkSmartPointer< vtkDelaunay2D > del2D =
	vtkSmartPointer< vtkDelaunay2D >::New();
      del2D->SetInputData(fullpd);
      del2D->SetSourceData(fullpd);
      del2D->SetTolerance(0.0);
      del2D->SetAlpha(0.0);
      del2D->SetOffset(0);
      del2D->SetProjectionPlaneMode(VTK_SET_TRANSFORM_PLANE);
      del2D->SetTransform(transform);
      del2D->BoundingTriangulationOff();
      del2D->Update();

//      WriteVTPFile(del2D->GetOutput(),"NewLoop"+int2String(cellId));
      vtkCellArray *polys = del2D->GetOutput()->GetPolys();

      // Renumber the point IDs.
      vtkIdType npts, *ptIds;
      for (polys->InitTraversal(); polys->GetNextCell(npts, ptIds); )
      {
	if ( ptIds[0] >= points->GetNumberOfPoints() ||
	     ptIds[1] >= points->GetNumberOfPoints() ||
	     ptIds[2] >= points->GetNumberOfPoints() )
	{
	  vtkGenericWarningMacro( << "Invalid point ID!!!");
	}

	splitCells->InsertNextCell( npts );
	int interPtCount=0;
	int interPts[3];
	for (int i = 0; i < npts; i++)
	{
	  vtkIdType remappedPtId;
	  if ( ptIds[i] < 3) // Point from the cell
	  {
	    remappedPtId = cellPts[ ptIds[i] ];
	    std::cout<<"Cell has point "<<remappedPtId<<endl;
	    if (CellPointOnInterLine[ptIds[i]])
	    {
	      interPts[interPtCount++] = reverseLineIdMap[ptIds[i] ];
	    }
	  }
	  else
	  {
	    remappedPtId = reverseIdMap[ ptIds[i] ];
	    std::cout<<"Cell has point "<<remappedPtId<<endl;
	    interPts[interPtCount++] = reverseLineIdMap[ptIds[i] ];
	  }
	  splitCells->InsertCellPoint( remappedPtId );
	}
	std::cout<<"Inter Point Count! "<<interPtCount<<endl;
	if (interPtCount >= 2)
	{
	  this->AddToNewCellMap(inputIndex, interPtCount,interPts,
	      interLines,numCurrCells);
	}
        std::cout<<"Num Current Cells "<<numCurrCells<<endl;
        numCurrCells++;
      }
    }

#else
  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->SetPoints(points);
  pd->SetLines(lines);

  //
  // Set up vtkPolyData to feed to vtkDelaunay2D
  //
  vtkSmartPointer< vtkDelaunay2D > del2D =
    vtkSmartPointer< vtkDelaunay2D >::New();
  del2D->SetInputData(pd);
  del2D->SetSourceData(pd);
  del2D->SetTolerance(0.0);
  del2D->SetAlpha(0.0);
  del2D->SetOffset(0);
  del2D->SetProjectionPlaneMode(VTK_SET_TRANSFORM_PLANE);
  del2D->SetTransform(transform);
  del2D->BoundingTriangulationOff();
  del2D->Update();

  if (del2D->GetOutput()->GetNumberOfCells() == 0)
  {
    std::cout<<"That bad thing happened!!!"<<endl;
    std::cout<<"Rotation angle: "<<rotationAngle<<endl;
    std::cout<<"Rotation axis: "<<rotationAxis[0]<<" "<<rotationAxis[1]<<" "<<rotationAxis[2]<<endl;
    std::cout<<"Normal: "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<<endl;
    std::cout<<"Center: "<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;
  }

  vtkCellArray *polys = del2D->GetOutput()->GetPolys();
  splitCells->Allocate( 3*polys->GetNumberOfCells() );

  // Renumber the point IDs.
  vtkIdType npts, *ptIds;
  for (polys->InitTraversal(); polys->GetNextCell(npts, ptIds); )
  {
    if ( ptIds[0] >= points->GetNumberOfPoints() ||
         ptIds[1] >= points->GetNumberOfPoints() ||
         ptIds[2] >= points->GetNumberOfPoints() )
    {
      vtkGenericWarningMacro( << "Invalid point ID!!!");
    }

    splitCells->InsertNextCell( npts );
    int interPtCount=0;
    int interPts[3];
    for (int i = 0; i < npts; i++)
    {
      vtkIdType remappedPtId;
      if ( ptIds[i] < 3) // Point from the cell
      {
        remappedPtId = cellPts[ ptIds[i] ];
	std::cout<<"Cell has point "<<remappedPtId<<endl;
	if (CellPointOnInterLine[ptIds[i]])
	{
	  interPts[interPtCount++] = reverseLineIdMap[ptIds[i] ];
	}
      }
      else
      {
        remappedPtId = reverseIdMap[ ptIds[i] ];
	std::cout<<"Cell has point "<<remappedPtId<<endl;
        interPts[interPtCount++] = reverseLineIdMap[ptIds[i] ];
      }
      splitCells->InsertCellPoint( remappedPtId );
    }
    std::cout<<"Inter Point Count! "<<interPtCount<<endl;
    if (interPtCount >= 2)
    {
      this->AddToNewCellMap(inputIndex, interPtCount,interPts,
	  interLines,numCurrCells);
    }
    std::cout<<"Num Current Cells "<<numCurrCells<<endl;
    numCurrCells++;
  }
#endif
  delete [] interPtBool;
  return splitCells;
}

void vtkIntersectionPolyDataFilterMine::Impl
::Orient(vtkPolyData *pd,vtkTransform *transform,vtkPolyData *boundary,
		vtkPolygon *boundarypoly)
{
  vtkSmartPointer<vtkTransformPolyDataFilter> transformer = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  vtkSmartPointer<vtkPolyData> transformedpd = 
    vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> newcells = 
    vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> cellarray = 
    vtkSmartPointer<vtkCellArray>::New();

  transformer->SetInputData(pd);
  transformer->SetTransform(transform);
  transformer->Update();
  transformedpd = transformer->GetOutput();

  double area = 0;
  double tedgept1[3];
  double tedgept2[3];
  int iter=0;
  vtkIdType nextPt;
  std::cout<<"Number of Points in Transform: "<<pd->GetNumberOfPoints()<<endl;
  for(nextPt=0;nextPt<pd->GetNumberOfPoints()-1;nextPt++)
  {
      transformedpd->GetPoint(nextPt,tedgept1);
      transformedpd->GetPoint(nextPt+1,tedgept2);
      area = area + (tedgept1[0]*tedgept2[1])-(tedgept2[0]*tedgept1[1]);
      std::cout<<"Area is "<<area<<" for "<<nextPt<<" and "<<nextPt+1<<endl;
  }
  transformedpd->GetPoint(nextPt,tedgept1);
  transformedpd->GetPoint(0,tedgept2);
  area = area + (tedgept1[0]*tedgept2[1])-(tedgept2[0]*tedgept1[1]);
  std::cout<<"Area is "<<area<<" for "<<nextPt<<" and "<<"0"<<endl;

  std::cout<<"AREA "<<area<<endl;
  if (area < 0)
  {
    std::cout<<"Flipping!"<<endl;
    for (nextPt=pd->GetNumberOfPoints()-1;nextPt>-1;nextPt--)
      boundarypoly->GetPointIds()->InsertNextId(nextPt);
  }
  else 
  {
    for (nextPt=0;nextPt<pd->GetNumberOfPoints();nextPt++)
      boundarypoly->GetPointIds()->InsertNextId(nextPt);
  }
  cellarray->InsertNextCell(boundarypoly);
  boundary->SetPoints(pd->GetPoints());
  boundary->SetPolys(cellarray);
}

void vtkIntersectionPolyDataFilterMine::Impl
::GetLoops(vtkPolyData *pd,std::vector<simPolygon> *loops)
{
  vtkSmartPointer<vtkIdList> pointCells = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
  simPoint nextPt;
  vtkIdType nextCell;
  double maxangle = 0;
  int numLoops=0;
  bool *ptBool;
  ptBool = new bool[pd->GetNumberOfPoints()];
  bool *lineBool;
  lineBool = new bool[pd->GetNumberOfCells()];

  int iter = 0;
  int used = -1;
  vtkIdType npts;
  vtkIdType *pts;

  for (vtkIdType ptId = 0;ptId <pd->GetNumberOfPoints();ptId++)
  {
    ptBool[ptId] = false;
  }
  std::cout<<"Number Of Cells: "<<pd->GetNumberOfCells()<<endl;
  for (vtkIdType lineId = 0;lineId <pd->GetNumberOfCells();lineId++)
  {
    lineBool[lineId] = false;
  }

  for (vtkIdType ptId = 0;ptId <pd->GetNumberOfPoints();ptId++)
  {
    if (ptBool[ptId] == false)
    {
      ptBool[ptId] = true;
      nextPt.id = ptId;
      pd->GetPoint(nextPt.id,nextPt.pt);
      simPolygon interloop;
      interloop.points.push_back(nextPt);
      pd->GetPointCells(nextPt.id,pointCells);
      if (pointCells->GetNumberOfIds() > 2)
      {
	std::cout<<"Tell me that there are more than two cells!"<<endl;
      }
      nextCell = pointCells->GetId(0);
      lineBool[nextCell] = true;
      this->GetSingleLoop(pd,&interloop,nextCell,ptBool,lineBool);
      loops->push_back(interloop);
      std::cout<<"INCREMETING LOOP COUNT!"<<endl;
      numLoops++;
    }
  }
  for (vtkIdType lineId = 0;lineId <pd->GetNumberOfCells();lineId++)
  {
    if (lineBool[lineId] == false)
    {
      std::cout<<"LINE FALSE: Find extra loop/s"<<endl;
      lineBool[lineId] = true;
      nextCell = lineId;
      pd->GetCellPoints(lineId,cellPoints);
      nextPt.id = cellPoints->GetId(0);
      pd->GetPoint(nextPt.id,nextPt.pt);
      simPolygon interloop;
      interloop.points.push_back(nextPt);
      this->GetSingleLoop(pd,&interloop,nextCell,ptBool,lineBool);
      loops->push_back(interloop);
      std::cout<<"INCREMETING LOOP COUNT!"<<endl;
      numLoops++;
    }
  }

  delete [] ptBool;
  delete [] lineBool;
}

void vtkIntersectionPolyDataFilterMine::Impl
::Triangulate(vtkPolyData *pd,simPolygon *loop, vtkCellArray *polys)
{     
  int foundear = 1;
  int sz;
  std::list<simPoint>::iterator it;
  std::list<simPoint>::iterator nit;
  std::list<simPoint>::iterator pit;
  std::cout<<"ORDER: ";
  for (it = loop->points.begin();it != loop->points.end();++it)
  {
    std::cout<<(it)->id<<" ";
  }
  std::cout<<endl;
  int iter = 0;
  while ((sz = loop->points.size()) > 3)
  {
    std::cout<<"Number of points in loop: "<<sz<<endl;
    double p1[3],p2[3],p3[3];
    double l1[3],l2[3];

    if (foundear)
    {
      it = loop->points.begin();
      iter = 0;
    }
    if (iter == 0)
    {
      std::cout<<"ITER is 0"<<endl;
      nit = std::next(it,1);
      pit = std::next(it,sz-1);
      std::cout<<"ID is "<<it->id<<endl;
      std::cout<<"Next id is "<<nit->id<<endl;
      std::cout<<"Previous id is "<<pit->id<<endl;
    }
    else if (iter == sz-1)
    {
      nit = std::prev(it,sz-1);
      pit = std::prev(it,1);
    }
    else
    {
      nit = std::next(it,1);
      pit = std::prev(it,1);
    }
    for (int j=0;j<2;j++)
    {
      p1[j] = (pit)->pt[j];
      p2[j] = (it)->pt[j];
      p3[j] = (nit)->pt[j];
      l1[j] = p2[j]-p1[j];
      l2[j] = p2[j]-p3[j];
    }
    l1[2] = 0.;  l2[2] = 0.;
    p1[2] = 0.;  p2[2] = 0.;  p3[2] = 0.;
    int concave = this->TestConcave(l1,l2,loop->orientation);
    std::cout<<"NORMAL POINT "<<(it)->id<<" IS CONCAVE "<<concave<<endl;
    (it)->concave = concave;
    if (!((it)->concave))
    {
      std::cout<<"New line should be with points :"<<(pit)->id<<" and "<< (nit)->id<<endl;
      int ear = this->Clip(p1,p3,loop,pit,it); 
      if (ear)
      {
	std::cout<<"Found an ear son!!"<<endl;
	polys->InsertNextCell(3);
	polys->InsertCellPoint((pit)->id);
	polys->InsertCellPoint((it)->id);
	polys->InsertCellPoint((nit)->id);
	std::cout<<"Removing point"<<endl;
	loop->points.erase(it);
	foundear = 1;
      }
      else
	foundear = 0;

    }
    else
      foundear = 0;
    //std::cout<<" POINT IS "<<(it)->pt[0]<<" "<<(it)->pt[1]<<endl;
    iter++;
    ++it;
  }
  pit = loop->points.begin();
  it = std::next(pit,1);
  nit = std::next(pit,2);
  polys->InsertNextCell(3);
  polys->InsertCellPoint((pit)->id);
  polys->InsertCellPoint((it)->id);
  polys->InsertCellPoint((nit)->id);
}

int vtkIntersectionPolyDataFilterMine::Impl
::Clip(double newp1[3], double newp2[3],simPolygon *loop,
    std::list<simPoint>::iterator it1, std::list<simPoint>::iterator it2)
{
  int sz = loop->points.size();
  int iter = 0;
  int val = 0;
  int ear = 1;
  std::list<simPoint>::iterator cit;
  std::list<simPoint>::iterator cnit;
  double t1,t2;
  double ip1[3],ip2[3];
  for (cit = loop->points.begin();cit != loop->points.end();++cit)
  {
    if (cit != it1 && cit != it2)
    {
      if (iter == sz-1)
	cnit = std::prev(cit,sz-1);
      else
	cnit = std::next(cit,1);

      for (int j=0;j<2;j++)
      {
	ip1[j] = (cit)->pt[j];
	ip2[j] = (cnit)->pt[j];
      }
      ip1[2] = 0.;  ip2[2] = 0.;

      //val equals 0 if no intersection (parallel lines!)
      //val equals 2 if there is an intersection, check parametric ts
      //val equals 3 if the lines are colinear
      val = vtkLine::Intersection(newp1,newp2,ip1,ip2,t1,t2);
      std::cout<<"Points "<<cit->id<<" and "<<cnit->id<<" INTERSECTS: "<<val<<endl;
      std::cout<<"t vals "<<t1<<" "<<t2<<endl;
      if ((val == 2 && (t1 > 0 && t1 < 1) && (t2 > 0 && t2 < 1)) ||
	  (val == 3))
      { 
	std::cout<<"Intersects!!!"<<endl;
	ear = 0;
      }
    }
    iter++;
  }

  return ear;
}

int vtkIntersectionPolyDataFilterMine::Impl
::TestConcave(double l1[3],double l2[3],int orientation)
{
  int concave = 0;
  vtkMath::Normalize( l1 );
  vtkMath::Normalize( l2 );

  std::cout<<"Line 1: "<<l1[0]<<" "<<l1[1]<<endl;
  std::cout<<"Angle 1 "<<atan2(l1[1],l1[0])<<endl;
  std::cout<<"Line 2: "<<l2[0]<<" "<<l2[1]<<endl;
  std::cout<<"Angle 2 "<<atan2(l2[1],l2[0])<<endl;
  double radangle = atan2(l1[1],l1[0]) - atan2(l2[1],l2[0]);
  std::cout<<"In between ang: "<<vtkMath::DegreesFromRadians(radangle);
  std::cout<<"And orientation "<<orientation<<endl;
  double angle = orientation * vtkMath::DegreesFromRadians(radangle);
  if (angle < 0)
    angle += 360;
  std::cout<<"Angle is "<<angle<<endl;
  if (angle <= 0 || angle >= 180)
    concave = 1;

   return concave;
}

void vtkIntersectionPolyDataFilterMine::Impl
::GetSingleLoop(vtkPolyData *pd,simPolygon *loop,vtkIdType nextCell,bool *interPtBool, bool *lineBool)
{
  vtkSmartPointer<vtkIdList> pointCells = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
  double mincell = 0;
  double l0pt0[3],l0pt1[3],l1pt0[3],l1pt1[3];
  vtkIdType prevPt;
  int intertype = 0;

  vtkIdType nextPt = loop->points.front().id;
  vtkIdType startPt = nextPt;
  std::cout<<"Start Point is: "<<startPt<<endl;
  interPtBool[nextPt] = true;
  pd->GetCellPoints(nextCell,cellPoints);
  if (cellPoints->GetId(0) == nextPt)
  {
    simPoint newpoint; newpoint.id = cellPoints->GetId(1);
    pd->GetPoint(newpoint.id,newpoint.pt);
    loop->points.push_back(newpoint);
    prevPt = nextPt;
    nextPt = cellPoints->GetId(1);
  }
  else 
  {
    simPoint newpoint; newpoint.id = cellPoints->GetId(0);
    pd->GetPoint(newpoint.id,newpoint.pt);
    loop->points.push_back(newpoint);
    prevPt = nextPt;
    nextPt = cellPoints->GetId(0);
  }
  interPtBool[nextPt] = true;
  std::cout<<"Next Point is: "<<nextPt<<endl;

  while(nextPt != startPt)
  {
    pd->GetPointCells(nextPt,pointCells);
    if (pointCells->GetNumberOfIds() > 2)
    {
      std::cout<<"Number of cells are: "<<pointCells->GetNumberOfIds()<<endl;
      std::cout<<"Inter type is: "<<intertype<<endl;
      if (intertype == 0)
      {
	this->SetLoopOrientation(pd,loop,&nextCell,nextPt,prevPt,pointCells);
	intertype = 1;
      }
      else 
      {
	this->FollowLoopOrientation(pd,loop,&nextCell,nextPt,prevPt,pointCells);
      }
      std::cout<<"While looping, point has more than two cells"<<endl;
    }
    else
    {
      if (pointCells->GetId(0) == nextCell)
	nextCell = pointCells->GetId(1);
      else
	nextCell = pointCells->GetId(0);
    }
    lineBool[nextCell] = true;

    pd->GetCellPoints(nextCell,cellPoints);
    if (cellPoints->GetId(0) == nextPt)
    {
      simPoint newpoint; newpoint.id = cellPoints->GetId(1);
      pd->GetPoint(newpoint.id,newpoint.pt);
      loop->points.push_back(newpoint);
      prevPt = nextPt;
      nextPt = cellPoints->GetId(1);
    }
    else 
    {
      simPoint newpoint; newpoint.id = cellPoints->GetId(0);
      pd->GetPoint(newpoint.id,newpoint.pt);
      loop->points.push_back(newpoint);
      prevPt = nextPt;
      nextPt = cellPoints->GetId(0);
    }
    std::cout<<"Next Point is :"<<nextPt<<endl;
    interPtBool[nextPt] = true;
  }
  //Cell is boring; i.e. it only has boundary points. set the orientation
  if (intertype == 0)
  {
    nextPt = 0;
    pd->GetPointCells(nextPt,pointCells);
    nextCell = pointCells->GetId(0);
    pd->GetCellPoints(pointCells->GetId(1),cellPoints);
    if (cellPoints->GetId(0) == nextPt)
      prevPt = cellPoints->GetId(1);
    else
      prevPt = cellPoints->GetId(0);

    loop->orientation = this->GetLoopOrientation(pd,nextCell,prevPt,nextPt);
  }
}

//----------------------------------------------------------------------------
void vtkIntersectionPolyDataFilterMine::Impl
::SetLoopOrientation(vtkPolyData *pd,simPolygon *loop,vtkIdType *nextCell,
    vtkIdType nextPt, vtkIdType prevPt, vtkIdList *pointCells)
{
double l0pt0[3],l0pt1[3],l1pt0[3],l1pt1[3];
double mincell = 0;
double minangle = VTK_DOUBLE_MAX;
for (vtkIdType i=0;i<pointCells->GetNumberOfIds();i++)
{
  vtkIdType cellId = pointCells->GetId(i);
  if (*nextCell != cellId)
  {
    pd->GetPoint(prevPt,l0pt0);
    pd->GetPoint(nextPt,l0pt1);
    std::cout<<"Edge 1 pt 1: "<<prevPt<<endl;
    std::cout<<"Edge 1 pt 2: "<<nextPt<<endl;
    vtkSmartPointer<vtkIdList> specialCellPoints = 
      vtkSmartPointer<vtkIdList>::New();
    pd->GetCellPoints(cellId,specialCellPoints);
    if (specialCellPoints->GetId(0) == nextPt)
    {
      pd->GetPoint(specialCellPoints->GetId(1),l1pt0);
      pd->GetPoint(specialCellPoints->GetId(0),l1pt1);
      std::cout<<"Edge 2 pt 1: "<<specialCellPoints->GetId(1)<<endl;
      std::cout<<"Edge 2 pt 2: "<<specialCellPoints->GetId(0)<<endl;
    }
    else
    {
      pd->GetPoint(specialCellPoints->GetId(0),l1pt0);
      pd->GetPoint(specialCellPoints->GetId(1),l1pt1);
      std::cout<<"Edge 2 pt 1: "<<specialCellPoints->GetId(0)<<endl;
      std::cout<<"Edge 2 pt 2: "<<specialCellPoints->GetId(1)<<endl;
    }
    double edge1[3],edge2[3];
    for (int j=0;j<2;j++)
    {
      edge1[j] = l0pt1[j]-l0pt0[j];
      edge2[j] = l1pt1[j]-l1pt0[j];
    }
    edge1[2] = 0.;
    edge2[2] = 0.;
    vtkMath::Normalize( edge1 );
    vtkMath::Normalize( edge2 );
    double angle = 
      vtkMath::DegreesFromRadians(acos(vtkMath::Dot(edge1, edge2)));
    std::cout<<"Edge 1: "<<edge1[0]<<" "<<edge1[1]<<" "<<edge1[2]<<endl;
    std::cout<<"Edge 2: "<<edge2[0]<<" "<<edge2[1]<<" "<<edge2[2]<<endl;
    std::cout<<"Angle check "<<angle<<endl;

    if (angle < minangle)
    {
      minangle = angle;
      mincell = cellId;
    }
  }
}
*nextCell = mincell;
loop->orientation = this->GetLoopOrientation(pd,*nextCell,prevPt,nextPt);
std::cout<<"Orientation set as: "<<loop->orientation<<endl;
}

//----------------------------------------------------------------------------
void vtkIntersectionPolyDataFilterMine::Impl
::FollowLoopOrientation(vtkPolyData *pd,simPolygon *loop,vtkIdType *nextCell,
    vtkIdType nextPt, vtkIdType prevPt, vtkIdList *pointCells)
{
  double l0pt0[3],l0pt1[3],l1pt0[3],l1pt1[3];
  double newcell = 0;
  double minangle = VTK_DOUBLE_MAX;
  for (vtkIdType i=0;i<pointCells->GetNumberOfIds();i++)
  {
    vtkIdType cellId = pointCells->GetId(i);
    if (*nextCell != cellId)
    {
      int neworient = this->GetLoopOrientation(pd,cellId,prevPt,nextPt);

      std::cout<<"Cell orientation is: "<<neworient<<endl;
      if (neworient == loop->orientation)
      {
	pd->GetPoint(prevPt,l0pt0);
	pd->GetPoint(nextPt,l0pt1);
	vtkSmartPointer<vtkIdList> specialCellPoints = 
	  vtkSmartPointer<vtkIdList>::New();
	pd->GetCellPoints(cellId,specialCellPoints);
	if (specialCellPoints->GetId(0) == nextPt)
	{
	  pd->GetPoint(specialCellPoints->GetId(1),l1pt0);
	  pd->GetPoint(specialCellPoints->GetId(0),l1pt1);
	  std::cout<<"Edge 2 pt 1: "<<specialCellPoints->GetId(1)<<endl;
	  std::cout<<"Edge 2 pt 2: "<<specialCellPoints->GetId(0)<<endl;
	}
	else
	{
	  pd->GetPoint(specialCellPoints->GetId(0),l1pt0);
	  pd->GetPoint(specialCellPoints->GetId(1),l1pt1);
	  std::cout<<"Edge 2 pt 1: "<<specialCellPoints->GetId(0)<<endl;
	  std::cout<<"Edge 2 pt 2: "<<specialCellPoints->GetId(1)<<endl;
	}
	double edge1[3],edge2[3];
	for (int j=0;j<2;j++)
	{
	  edge1[j] = l0pt1[j]-l0pt0[j];
	  edge2[j] = l1pt1[j]-l1pt0[j];
	}
	edge1[2] = 0.;
	edge2[2] = 0.;
	vtkMath::Normalize( edge1 );
	vtkMath::Normalize( edge2 );
	double angle = 
	  vtkMath::DegreesFromRadians(acos(vtkMath::Dot(edge1, edge2)));
	std::cout<<"Found same orientation cell!!!!"<<endl;
	if (angle < minangle)
	{
	  minangle = angle;
	  newcell = cellId;
	}
      }
    }
  }
  *nextCell = newcell;
}
//----------------------------------------------------------------------------
int vtkIntersectionPolyDataFilterMine::Impl
::AddToPointEdgeMap(int index, vtkIdType ptId, double x[3], vtkPolyData *mesh,
                    vtkIdType cellId, vtkIdType edgeId, vtkIdType lineId,
                    vtkIdType triPtIds[3])
{
  int val = -1;
  vtkIdType edgePtId0 = triPtIds[edgeId];
  vtkIdType edgePtId1 = triPtIds[(edgeId+1) % 3];
  double pt0[3], pt1[3];

  mesh->GetPoint(edgePtId0, pt0);
  mesh->GetPoint(edgePtId1, pt1);

  // Check to see if this point-cell combo is already in the list
  PointEdgeMapIteratorType iterLower =
    this->PointEdgeMap[index]->lower_bound( ptId );
  PointEdgeMapIteratorType iterUpper =
    this->PointEdgeMap[index]->upper_bound( ptId );

  while (iterLower != iterUpper)
    {
    if ( iterLower->second.CellId == cellId )
      {
      return iterLower->second.EdgeId;
      }
    ++iterLower;
    }

  double t, dist, closestPt[3];
  dist = vtkLine::DistanceToLine(x, pt0, pt1, t, closestPt);
  if (fabs(dist) < 1e-15 && t >= 0.0 && t <= 1.0)
    {
    CellEdgeLineType cellEdgeLine;
    cellEdgeLine.CellId = cellId;
    cellEdgeLine.EdgeId = edgeId;
    cellEdgeLine.LineId = lineId;
    this->PointEdgeMap[index]->insert(std::make_pair(ptId, cellEdgeLine));
    val = edgeId;
    }
  return val;
}

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkIntersectionPolyDataFilterMine);

//----------------------------------------------------------------------------
vtkIntersectionPolyDataFilterMine::vtkIntersectionPolyDataFilterMine()
  : SplitFirstOutput(1), SplitSecondOutput(1)
{
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(3);
}

//----------------------------------------------------------------------------
vtkIntersectionPolyDataFilterMine::~vtkIntersectionPolyDataFilterMine()
{
}

//----------------------------------------------------------------------------
void vtkIntersectionPolyDataFilterMine::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "SplitFirstOutput: " << this->SplitFirstOutput << "\n";
  os << indent << "SplitSecondOutput: " << this->SplitSecondOutput << "\n";
}

//----------------------------------------------------------------------------
int vtkIntersectionPolyDataFilterMine::TriangleTriangleIntersection(double p1[3], 
    					double q1[3], double r1[3], double p2[3], 
					double q2[3], double r2[3], int &coplanar, 
					double pt1[3], double pt2[3], double surfaceid[2],
					int &pointCase)
{
  double n1[3], n2[3];

  // Compute supporting plane normals.
  vtkTriangle::ComputeNormal(p1, q1, r1, n1);
  vtkTriangle::ComputeNormal(p2, q2, r2, n2);
  double s1 = -vtkMath::Dot(n1, p1);
  double s2 = -vtkMath::Dot(n2, p2);

  // Compute signed distances of points p1, q1, r1 from supporting
  // plane of second triangle.
  double dist1[3];
  dist1[0] = vtkMath::Dot(n2, p1) + s2;
  dist1[1] = vtkMath::Dot(n2, q1) + s2;
  dist1[2] = vtkMath::Dot(n2, r1) + s2;

  // If signs of all points are the same, all the points lie on the
  // same side of the supporting plane, and we can exit early.
  if ((dist1[0]*dist1[1] > 0.0) && (dist1[0]*dist1[2] > 0.0))
    {
    //std::cout<<"Same side supporting plane 1!"<<endl;
    return 0;
    }
  // Do the same for p2, q2, r2 and supporting plane of first
  // triangle.
  double dist2[3];
  dist2[0] = vtkMath::Dot(n1, p2) + s1;
  dist2[1] = vtkMath::Dot(n1, q2) + s1;
  dist2[2] = vtkMath::Dot(n1, r2) + s1;

  // If signs of all points are the same, all the points lie on the
  // same side of the supporting plane, and we can exit early.
  if ((dist2[0]*dist2[1] > 0.0) && (dist2[0]*dist2[2] > 0.0))
    {
    //std::cout<<"Same side supporting plane 2!"<<endl;
    return 0;
    }
  // Check for coplanarity of the supporting planes.
  if ( fabs( n1[0] - n2[0] ) < 1e-9 &&
       fabs( n1[1] - n2[1] ) < 1e-9 &&
       fabs( n1[2] - n2[2] ) < 1e-9 &&
       fabs( s1 - s2 ) < 1e-9 )
    {
    coplanar = 1;
    //std::cout<<"Coplanar!"<<endl;
    return 0;
    }

  coplanar = 0;

  // There are more efficient ways to find the intersection line (if
  // it exists), but this is clear enough.
  double *pts1[3] = {p1, q1, r1}, *pts2[3] = {p2, q2, r2};

  // Find line of intersection (L = p + t*v) between two planes.
  double n1n2 = vtkMath::Dot(n1, n2);
  double a = (s1 - s2*n1n2) / (n1n2*n1n2 - 1.0);
  double b = (s2 - s1*n1n2) / (n1n2*n1n2 - 1.0);
  double p[3], v[3];
  p[0] = a*n1[0] + b*n2[0];
  p[1] = a*n1[1] + b*n2[1];
  p[2] = a*n1[2] + b*n2[2];
  vtkMath::Cross(n1, n2, v);
  vtkMath::Normalize( v );

  int index1 = 0, index2 = 0;
  double t1[3], t2[3];
  int ts1,ts2;
  for (int i = 0; i < 3; i++)
    {
    double t, x[3];
    int id1 = i, id2 = (i+1) % 3;

    // Find t coordinate on line of intersection between two planes.
    if (vtkPlane::IntersectWithLine( pts1[id1], pts1[id2], n2, p2, t, x ))
      {
	 if (t == 1)
	   ts1 = index1;
         t1[index1++] = vtkMath::Dot(x, v) - vtkMath::Dot(p, v);
      std::cout<<"T val 1 "<<t<<endl;
      }

    if (vtkPlane::IntersectWithLine( pts2[id1], pts2[id2], n1, p1, t, x ))
      {
	if (t == 1)
	  ts2 = index2;
        t2[index2++] = vtkMath::Dot(x, v) - vtkMath::Dot(p, v);
      std::cout<<"T val 2 "<<t<<endl;
      }
    }

  // Check if only one edge or all edges intersect the supporting
  // planes intersection.
  // TODO Check here first for weird intersection lines!
   std::cout<<"Index 1 "<<index1<<" and Index 2 "<<index2<<endl;
  if ( index1 > 2)
  {
    double tmp;
    index1--;
    tmp = t1[ts1];
    t1[ts1] = t1[2];
    t1[2] = tmp;
  } 
  if ( index2 > 2)
  {
    double tmp;
    index2--;
    tmp = t2[ts2];
    t2[ts2] = t2[2];
    t2[2] = tmp;
  }
  if ( index1 != 2 || index2 != 2 )
    {
    std::cout<<"Only one edge intersecting!"<<endl;
    return 0;
    }
  else 
  {

  }

  // Check for NaNs
  if ( vtkMath::IsNan( t1[0] ) || vtkMath::IsNan( t1[1] ) ||
       vtkMath::IsNan( t2[0] ) || vtkMath::IsNan( t2[1] ) )
    {
    std::cout<<"NaNs!"<<endl;
    return 0;
    }

  if ( t1[0] > t1[1] )
    {
    std::swap( t1[0], t1[1] );
    }
  if ( t2[0] > t2[1] )
    {
    std::swap( t2[0], t2[1] );
    }
  // Handle the different interval configuration cases.
  double tt1, tt2;
  if ( t1[1] < t2[0] || t2[1] < t1[0] )
    {
    std::cout<<"No Overlap!"<<endl;
    return 0; // No overlap
    }
  else if ( t1[0] < t2[0] )
    {
    if ( t1[1] < t2[1] )
      {
      //First point on surface 2, second point on surface 1
      pointCase = 1;
      surfaceid[0] = 2;
      surfaceid[1] = 1;
      tt1 = t2[0];
      tt2 = t1[1];
      }
    else
      { 
      std::cout<<"THE T!: "<<t1[0]<<" "<<t1[1]<<endl;
      std::cout<<"THE T!: "<<t2[0]<<" "<<t2[1]<<endl;
      //Both points belong to lines on surface 2
      pointCase = 2;
      surfaceid[0] = 2;
      surfaceid[1] = 2;
      tt1 = t2[0];
      tt2 = t2[1];
      }
    }
  else // t1[0] >= t2[0]
    {
    if ( t1[1] < t2[1] )
      {
      //Both points belong to lines on surface 1
      pointCase = 3;
      surfaceid[0] = 1;
      surfaceid[1] = 1;
      tt1 = t1[0];
      tt2 = t1[1];
      }
    else
      {
      std::cout<<"THE T!: "<<t1[0]<<" "<<t1[1]<<endl;
      std::cout<<"THE T!: "<<t2[0]<<" "<<t2[1]<<endl;
      //First point on surface 1, second point on surface 2
      pointCase = 4;
      surfaceid[0] = 1;
      surfaceid[1] = 2;
      tt1 = t1[0];
      tt2 = t2[1];
      }
    }

  // Create actual intersection points.
  pt1[0] = p[0] + tt1*v[0];
  pt1[1] = p[1] + tt1*v[1];
  pt1[2] = p[2] + tt1*v[2];

  pt2[0] = p[0] + tt2*v[0];
  pt2[1] = p[1] + tt2*v[1];
  pt2[2] = p[2] + tt2*v[2];

  return 1;
}

//----------------------------------------------------------------------------
int vtkIntersectionPolyDataFilterMine::RequestData(vtkInformation*        vtkNotUsed(request),
                                               vtkInformationVector** inputVector,
                                               vtkInformationVector*  outputVector)
{
  vtkInformation* inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkInformation* inInfo1 = inputVector[1]->GetInformationObject(0);
  vtkInformation* outIntersectionInfo = outputVector->GetInformationObject(0);
  vtkInformation* outPolyDataInfo0 = outputVector->GetInformationObject(1);
  vtkInformation* outPolyDataInfo1 = outputVector->GetInformationObject(2);

  vtkPolyData *input0 = vtkPolyData::SafeDownCast(
    inInfo0->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *input1 = vtkPolyData::SafeDownCast(
    inInfo1->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *outputIntersection = vtkPolyData::SafeDownCast(
    outIntersectionInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer< vtkPoints > outputIntersectionPoints =
    vtkSmartPointer< vtkPoints >::New();
  outputIntersection->SetPoints(outputIntersectionPoints);

  vtkPolyData *outputPolyData0 = vtkPolyData::SafeDownCast(
    outPolyDataInfo0->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *outputPolyData1 = vtkPolyData::SafeDownCast(
    outPolyDataInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // Set up new poly data for the inputs to build cells and links.
  vtkSmartPointer< vtkPolyData > mesh0 = vtkSmartPointer< vtkPolyData >::New();
  mesh0->DeepCopy(input0);

  vtkSmartPointer< vtkPolyData > mesh1 = vtkSmartPointer< vtkPolyData >::New();
  mesh1->DeepCopy(input1);

  // Find the triangle-triangle intersections between mesh0 and mesh1
  vtkSmartPointer< vtkOBBTree > obbTree0 = vtkSmartPointer< vtkOBBTree >::New();
  obbTree0->SetDataSet(mesh0);
  obbTree0->SetNumberOfCellsPerNode(10);
  obbTree0->SetMaxLevel(1000000);
  obbTree0->SetTolerance(1e-6);
  obbTree0->AutomaticOn();
  obbTree0->BuildLocator();

  vtkSmartPointer<vtkPolyData> obbTree0Rep = vtkSmartPointer<vtkPolyData>::New();
  obbTree0->GenerateRepresentation(-1,obbTree0Rep);
  //WriteVTPFile(obbTree0Rep,"Tree0Rep");

  vtkSmartPointer< vtkOBBTree > obbTree1 = vtkSmartPointer< vtkOBBTree >::New();
  obbTree1->SetDataSet(mesh1);
  obbTree1->SetNumberOfCellsPerNode(10);
  obbTree1->SetMaxLevel(1000000);
  obbTree1->SetTolerance(1e-6);
  obbTree1->AutomaticOn();
  obbTree1->BuildLocator();

  vtkSmartPointer<vtkPolyData> obbTree1Rep = vtkSmartPointer<vtkPolyData>::New();
  obbTree1->GenerateRepresentation(-1,obbTree1Rep);
  //WriteVTPFile(obbTree1Rep,"Tree1Rep");

  // Set up the structure for determining exact triangle-triangle
  // intersections.
  vtkIntersectionPolyDataFilterMine::Impl *impl = new vtkIntersectionPolyDataFilterMine::Impl();
  impl->Mesh[0]  = mesh0;
  impl->Mesh[1]  = mesh1;
  impl->OBBTree1 = obbTree1;

  vtkSmartPointer< vtkCellArray > lines = vtkSmartPointer< vtkCellArray >::New();
  outputIntersection->SetLines(lines);
  impl->IntersectionLines = lines;

  // Add cell data arrays that map the intersection line to the cells
  // it splits.
  impl->CellIds[0] = vtkIdTypeArray::New();
  impl->CellIds[0]->SetName("Input0CellID");
  outputIntersection->GetCellData()->AddArray(impl->CellIds[0]);
  impl->CellIds[0]->Delete();
  impl->CellIds[1] = vtkIdTypeArray::New();
  impl->CellIds[1]->SetName("Input1CellID");
  outputIntersection->GetCellData()->AddArray(impl->CellIds[1]);
  impl->CellIds[1]->Delete();

  impl->PointCellIds[0] = vtkIdTypeArray::New();
  impl->PointCellIds[0]->SetName("PointCellsIDs");
  impl->PointCellIds[1] = vtkIdTypeArray::New();
  impl->PointCellIds[1]->SetName("PointCellsIDs");

  impl->SurfaceId = vtkIdTypeArray::New();
  impl->SurfaceId->SetName("SurfaceID");
  outputIntersection->GetPointData()->AddArray(impl->SurfaceId);
  impl->SurfaceId->Delete();

  impl->CaseId = vtkIdTypeArray::New();
  impl->CaseId->SetName("CaseID");
  outputIntersection->GetCellData()->AddArray(impl->CaseId);
  impl->CaseId->Delete();

  impl->NewCellIds[0] = vtkIdTypeArray::New();
  impl->NewCellIds[0]->SetNumberOfComponents(2);
  impl->NewCellIds[1] = vtkIdTypeArray::New();
  impl->NewCellIds[1]->SetNumberOfComponents(2);

  double bounds0[6], bounds1[6];
  mesh0->GetBounds(bounds0);
  mesh1->GetBounds(bounds1);
  for (int i = 0; i < 3; i++)
    {
    int minIdx = 2*i;
    int maxIdx = 2*i+1;
    if (bounds1[minIdx] < bounds0[minIdx])
      {
      bounds0[minIdx] = bounds1[minIdx];
      }
    if (bounds1[maxIdx] > bounds0[maxIdx])
      {
      bounds0[maxIdx] = bounds1[maxIdx];
      }
    }

  vtkSmartPointer< vtkPointLocator > pointMerger =
    vtkSmartPointer< vtkPointLocator >::New();
  pointMerger->SetTolerance(1e-6);
  pointMerger->InitPointInsertion(outputIntersection->GetPoints(), bounds0);
  impl->PointMerger = pointMerger;
  
  // This performs the triangle intersection search
  obbTree0->IntersectWithOBBTree
    (obbTree1, 0, vtkIntersectionPolyDataFilterMine::Impl::FindTriangleIntersections,
     impl);

  for (int i=0;i<2;i++)
  {
    for (vtkIdType interCellId=0;interCellId<outputIntersection->GetNumberOfLines();interCellId++)
    {
      impl->NewCellIds[i]->InsertTuple2(interCellId,-1,-1);
    }
  }

  std::cout<<"LINEPTSBEFORE "<<outputIntersection->GetNumberOfPoints()<<endl;
  //WriteVTPFile(outputIntersection,"BeforeIntersectionLines");
  vtkSmartPointer<vtkPolyData> tmpLines = vtkSmartPointer<vtkPolyData>::New();
  tmpLines->DeepCopy(outputIntersection);

  tmpLines->BuildLinks();
  vtkSmartPointer<vtkCleanPolyData> lineCleaner = 
	  vtkSmartPointer<vtkCleanPolyData>::New();
  lineCleaner->SetInputData(outputIntersection);
  lineCleaner->Update();
  outputIntersection->DeepCopy(lineCleaner->GetOutput());
  vtkSmartPointer< vtkPointLocator > linePtMapper =
    vtkSmartPointer< vtkPointLocator >::New();
  linePtMapper->SetDataSet(outputIntersection);
  linePtMapper->BuildLocator();
  double newpt[3];
  vtkIdType mapPtId=0;
  for (vtkIdType ptId=0;ptId<tmpLines->GetNumberOfPoints();ptId++)
  {
    tmpLines->GetPoint(ptId,newpt);
    mapPtId = linePtMapper->FindClosestPoint(newpt);
    impl->PointMapper->insert(std::make_pair(mapPtId,ptId));
  }
  //WriteVTPFile(tmpLines,"FullLinesNotClean");
  std::cout<<"LINEPTSAFTER "<<outputIntersection->GetNumberOfPoints()<<endl;

  impl->BoundaryPoints[0] = vtkIntArray::New();
  impl->BoundaryPoints[1] = vtkIntArray::New();
  // Split the first output if so desired
  if ( this->SplitFirstOutput )
    {
    mesh0->BuildLinks();
    impl->SplitMesh(0, outputPolyData0, outputIntersection);
    impl->BoundaryPoints[0]->SetName("BoundaryPoint");
    //outputPolyData0->GetPointData()->AddArray(impl->BoundaryPoints[0]);
    //outputPolyData0->GetPointData()->SetActiveScalars("BoundaryPoint");
    impl->CleanAndCheckSurface(outputPolyData0);
    outputPolyData0->BuildLinks();
    }
  else
    {
    outputPolyData0->ShallowCopy( mesh0 );
    }

  // Split the second output if desired
  if ( this->SplitSecondOutput )
    {
    mesh1->BuildLinks();
    impl->SplitMesh(1, outputPolyData1, outputIntersection);
    impl->BoundaryPoints[1]->SetName("BoundaryPoint");
    //outputPolyData1->GetPointData()->AddArray(impl->BoundaryPoints[1]);
    //outputPolyData1->GetPointData()->SetActiveScalars("BoundaryPoint");
    impl->CleanAndCheckSurface(outputPolyData1);
    outputPolyData1->BuildLinks();
    }
  else
    {
    outputPolyData1->ShallowCopy( mesh1 );
    }

  impl->NewCellIds[0]->SetName("NewCell0ID");
  outputIntersection->GetCellData()->AddArray(impl->NewCellIds[0]);
  impl->NewCellIds[0]->Delete();
  impl->NewCellIds[1]->SetName("NewCell1ID");
  outputIntersection->GetCellData()->AddArray(impl->NewCellIds[1]);
  impl->NewCellIds[1]->Delete();

  impl->BoundaryPoints[0]->Delete();
  impl->BoundaryPoints[1]->Delete();
  impl->PointCellIds[0]->Delete();
  impl->PointCellIds[1]->Delete();
  delete impl;

  return 1;
}

//----------------------------------------------------------------------------
int vtkIntersectionPolyDataFilterMine::FillInputPortInformation(int port,
                                                            vtkInformation *info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
  else if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 0);
    }
  return 1;
}

int vtkIntersectionPolyDataFilterMine::Impl::GetLoopOrientation(
    vtkPolyData *pd, vtkIdType cell, vtkIdType ptId1, vtkIdType ptId2)
{
  double area = 0;
  double pt1[3],pt2[3],pt3[3];
  vtkIdType ptId3;
  vtkSmartPointer<vtkIdList> cellPoints = 
    vtkSmartPointer<vtkIdList>::New();

  pd->GetCellPoints(cell,cellPoints);
  if (cellPoints->GetId(0) == ptId2)
  {
    ptId3 = cellPoints->GetId(1);
  }
  else
  {                              
    ptId3 = cellPoints->GetId(0);
  }

  pd->GetPoint(ptId1,pt1); pd->GetPoint(ptId2,pt2); pd->GetPoint(ptId3,pt3);
  area = area + (pt1[0]*pt2[1])-(pt2[0]*pt1[1]);
  area = area + (pt2[0]*pt3[1])-(pt3[0]*pt2[1]);
  area = area + (pt3[0]*pt1[1])-(pt1[0]*pt3[1]);

  int orientation = 1;

  if (area < 0)
    orientation = -1;

  return orientation;
}

int vtkIntersectionPolyDataFilterMine::Impl::CheckLine(
    vtkPolyData *pd, vtkIdType ptId1,vtkIdType ptId2)
{
  vtkSmartPointer<vtkIdList> pointCells1 = 
    vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> pointCells2 = 
    vtkSmartPointer<vtkIdList>::New();

  pd->GetPointCells(ptId1,pointCells1);
  pd->GetPointCells(ptId2,pointCells2);

  pointCells1->IntersectWith(pointCells2);

  int unique = 1;
  if (pointCells1->GetNumberOfIds() > 0)
    unique = 0;

  return unique;
}

void vtkIntersectionPolyDataFilterMine::Impl::AddToNewCellMap(
    int inputIndex, int interPtCount, int interPts[3],
    vtkPolyData *interLines,int numCurrCells)
{
  std::cout<<"POINTS OF INTER LINES ARE "<<interPts[0]<<" AND "<<interPts[1]<<endl;
  //std::cout<<"INTER POINTS ARE "<<interPts[0]<<" AND "<<interPts[1]<<endl;
  vtkIdList *cellIds[interPtCount];
  for (int i=0;i<interPtCount;i++)
  {
    cellIds[i] = vtkIdList::New();
    vtkSmartPointer<vtkIdList> temp = vtkSmartPointer<vtkIdList>::New();
    interLines->GetPointCells(interPts[i],cellIds[i]);
//	      std::cout<<"Connected Cell Ids are: ";
//	      for(vtkIdType j=0;j<cellIds[i]->GetNumberOfIds();j++)
//	      {
//		std::cout<<cellIds[i]->GetId(j)<<" ";
//	      }	
//	      std::cout<<" for cell given number "<<numCurrCells<<endl;
    if (i > 0)
    {
      temp->DeepCopy(cellIds[i-1]);
      temp->IntersectWith(cellIds[i]);
    }
    std::cout<<"Number Of Intersected "<<temp->GetNumberOfIds()<<endl;
    if (temp->GetNumberOfIds() > 0)
    {
      for (int j=0;j<temp->GetNumberOfIds();j++)
      {
	if (NewCellIds[inputIndex]->GetComponent(temp->GetId(j),0) == -1)
	{
	    std::cout<<"SETTTINGINGING TO "<<numCurrCells<<endl;
	  NewCellIds[inputIndex]->InsertComponent(temp->GetId(j),0,numCurrCells);
	}
	else
	{
	    std::cout<<"SETTTINGINGING TO "<<numCurrCells<<endl;
	  NewCellIds[inputIndex]->InsertComponent(temp->GetId(j),1,numCurrCells);
	}
      }
    }
  }
  if (interPtCount > 2)
  {
    cellIds[0]->IntersectWith(cellIds[interPtCount-1]);
    if (cellIds[0]->GetNumberOfIds() > 0)
    {
      for (int j=0;j<cellIds[0]->GetNumberOfIds();j++)
      {
	if (NewCellIds[inputIndex]->GetComponent(cellIds[0]->GetId(j),0) == -1)
	{
	    std::cout<<"SETTTINGINGING TO "<<numCurrCells<<endl;
	  NewCellIds[inputIndex]->InsertComponent(cellIds[0]->GetId(j),0,numCurrCells);
	}
	else
	{
	    std::cout<<"SETTTINGINGING TO "<<numCurrCells<<endl;
	  NewCellIds[inputIndex]->InsertComponent(cellIds[0]->GetId(j),1,numCurrCells);
	}
      }
    }
    for (int i=0;i<interPtCount;i++)
    {
      cellIds[i]->Delete();
    }
  }
}
