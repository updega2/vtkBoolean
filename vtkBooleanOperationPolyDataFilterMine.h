/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBooleanOperationPolyDataFilterMine.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBooleanOperationPolyDataFilterMine
// .SECTION Description
//
// Computes the boundary of the union, intersection, or difference
// volume computed from the volumes defined by two input surfaces. The
// two surfaces do not need to be manifold, but if they are not,
// unexpected results may be obtained. The resulting surface is
// available in the first output of the filter. The second output
// contains a set of polylines that represent the intersection between
// the two input surfaces.
//
// This code was contributed in the VTK Journal paper:
// "Boolean Operations on Surfaces in VTK Without External Libraries"
// by Cory Quammen, Chris Weigle C., Russ Taylor
// http://hdl.handle.net/10380/3262
// http://www.midasjournal.org/browse/publication/797

#ifndef __vtkBooleanOperationPolyDataFilterMine_h
#define __vtkBooleanOperationPolyDataFilterMine_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

#include "vtkDataSetAttributes.h" // Needed for CopyCells() method

class vtkIdList;

class VTKFILTERSGENERAL_EXPORT vtkBooleanOperationPolyDataFilterMine : public vtkPolyDataAlgorithm
{
public:
  // Description:
  // Construct object that computes the boolean surface.
  static vtkBooleanOperationPolyDataFilterMine *New();

  vtkTypeMacro(vtkBooleanOperationPolyDataFilterMine,
               vtkPolyDataAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  enum OperationType
  {
    VTK_UNION=0,
    VTK_INTERSECTION,
    VTK_DIFFERENCE
  };

  // Description:
  // Set the boolean operation to perform. Defaults to union.
  vtkSetClampMacro( Operation, int, VTK_UNION, VTK_DIFFERENCE );
  vtkGetMacro( Operation, int );
  void SetOperationToUnion()
  { this->SetOperation( VTK_UNION ); }
  void SetOperationToIntersection()
  { this->SetOperation( VTK_INTERSECTION ); }
  void SetOperationToDifference()
  { this->SetOperation( VTK_DIFFERENCE ); }

protected:
  vtkBooleanOperationPolyDataFilterMine();
  ~vtkBooleanOperationPolyDataFilterMine();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int FillInputPortInformation(int, vtkInformation*);

private:
  vtkBooleanOperationPolyDataFilterMine(const vtkBooleanOperationPolyDataFilterMine&); // Not implemented
  void operator=(const vtkBooleanOperationPolyDataFilterMine&); // Not implemented

  // Description:
  // PolyDatas for each surface out of intersection and also the intersection
  // lines
  vtkPolyData *OutputSurface;
  // Description:
  // Which operation to perform.
  // Can be VTK_UNION, VTK_INTERSECTION, or VTK_DIFFERENCE.
  int Operation;

  class Impl;

};

#endif
