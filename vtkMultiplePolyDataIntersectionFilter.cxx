/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMultiplePolyDataIntersectionFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMultiplePolyDataIntersectionFilter.h"

#include "vtkAlgorithmOutput.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataSetAttributes.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkBooleanOperationPolyDataFilter2.h"
#include "vtkTrivialProducer.h"
#include "vtkSmartPointer.h"
#include "vtkBoundingBox.h"
#include "vtkIntArray.h"

vtkStandardNewMacro(vtkMultiplePolyDataIntersectionFilter);

//----------------------------------------------------------------------------
vtkMultiplePolyDataIntersectionFilter::vtkMultiplePolyDataIntersectionFilter()
{
  this->ParallelStreaming = 0;
  this->UserManagedInputs = 0;
  this->NoIntersectionOutput = 1;
  this->PassInfoAsGlobal = 0;
  this->AssignSurfaceIds = 0;

  this->BooleanObject = vtkPolyData::New();
  this->IntersectionTable = NULL;
}

//----------------------------------------------------------------------------
vtkMultiplePolyDataIntersectionFilter::~vtkMultiplePolyDataIntersectionFilter()
{
  if (this->BooleanObject)
    BooleanObject->Delete();
}

//----------------------------------------------------------------------------
// Add a dataset to the list of data to append.
void vtkMultiplePolyDataIntersectionFilter::AddInputData(vtkPolyData *ds)
{
  if (this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "AddInput is not supported if UserManagedInputs is true");
    return;
    }
  this->Superclass::AddInputData(ds);
}

//----------------------------------------------------------------------------
// Remove a dataset from the list of data to append.
void vtkMultiplePolyDataIntersectionFilter::RemoveInputData(vtkPolyData *ds)
{
  if (this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "RemoveInput is not supported if UserManagedInputs is true");
    return;
    }

  if (!ds)
    {
    return;
    }
  int numCons = this->GetNumberOfInputConnections(0);
  for(int i=0; i<numCons; i++)
    {
    if (this->GetInput(i) == ds)
      {
      this->RemoveInputConnection(0,
        this->GetInputConnection(0, i));
      }
    }
}

//----------------------------------------------------------------------------
// make ProcessObject function visible
// should only be used when UserManagedInputs is true.
void vtkMultiplePolyDataIntersectionFilter::SetNumberOfInputs(int num)
{
  if (!this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "SetNumberOfInputs is not supported if UserManagedInputs is false");
    return;
    }

  // Ask the superclass to set the number of connections.
  this->SetNumberOfInputConnections(0, num);
}

//----------------------------------------------------------------------------
void vtkMultiplePolyDataIntersectionFilter::
SetInputDataByNumber(int num, vtkPolyData* input)
{
  vtkTrivialProducer* tp = vtkTrivialProducer::New();
  tp->SetOutput(input);
  this->SetInputConnectionByNumber(num, tp->GetOutputPort());
  tp->Delete();
}

//----------------------------------------------------------------------------
// Set Nth input, should only be used when UserManagedInputs is true.
void vtkMultiplePolyDataIntersectionFilter::
SetInputConnectionByNumber(int num,vtkAlgorithmOutput *input)
{
  if (!this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "SetInputConnectionByNumber is not supported if UserManagedInputs "<<
      "is false");
    return;
    }

  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, num, input);
}

int vtkMultiplePolyDataIntersectionFilter::BuildIntersectionTable(
    vtkPolyData* inputs[], int numInputs)
{
  for (int i = 0;i < numInputs;i++)
    {
    inputs[i]->ComputeBounds();
    }

  int totalIntersections=0;
  for (int i = 0;i < numInputs;i++)
    {
    if (this->AssignSurfaceIds)
      this->SetSurfaceId(inputs[i],i+1);
    int objectIntersections=0;
    double bounds0[6];
    inputs[i]->GetBounds(bounds0);
    vtkBoundingBox boundingBox0; 
    boundingBox0.SetBounds(bounds0);
    for (int j = 0;j < numInputs; j++)
      {
      if (i != j)
	{
	double bounds1[6];
	inputs[j]->GetBounds(bounds1);
	vtkBoundingBox boundingBox1; 
	boundingBox1.SetBounds(bounds1);
	int intersects = boundingBox0.Intersects(boundingBox1);
	if (intersects)
	  {
	  this->IntersectionTable[i][j] = 1;
	  objectIntersections++;
	  totalIntersections++;
	  }
	}
      }
      if (objectIntersections == 0)
        {
        vtkGenericWarningMacro( << "Input object "<<i<<" doesn't intersect "
	                        << "with any other input object." );
        }
    }
  return totalIntersections;
}

int vtkMultiplePolyDataIntersectionFilter::ExecuteIntersection(
    vtkPolyData* inputs[], int numInputs)
{                                            
  int numChecks = 0;  
  int totalIntersections = 0;
  vtkSmartPointer<vtkIdList> checkInputArray = 
    vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> checkInputArray2 = 
    vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> tmp = 
    vtkSmartPointer<vtkIdList>::New();

  this->BooleanObject->DeepCopy(inputs[0]);
  checkInputArray->InsertNextId(0);
  while ((numChecks = checkInputArray->GetNumberOfIds()) > 0)
    {
    for(int c = 0;c < numChecks; c++)
      {
      int i = checkInputArray->GetId(c);
      for (int j = 0;j < numInputs; j++)
        {
	  std::cout<<"INTERSECTION of "<<i<<" and "<<j<<endl;
	//Bounding boxes intersect!
        if (this->IntersectionTable[i][j] == 1)
          {
	    this->IntersectionTable[i][j] = -1;
	    this->IntersectionTable[j][i] = -1;

	  vtkSmartPointer<vtkBooleanOperationPolyDataFilter2> boolean = 
	    vtkSmartPointer<vtkBooleanOperationPolyDataFilter2>::New();
	  if (this->PassInfoAsGlobal && totalIntersections != 0)
	    this->PreSetGlobalArrays(inputs[j]);

	  boolean->SetInputData(0,this->BooleanObject);
	  boolean->SetInputData(1,inputs[j]);
	  boolean->SetNoIntersectionOutput(0);
	  boolean->SetOperationToUnion();
	  boolean->Update();

	  int numPts = boolean->GetNumberOfIntersectionPoints();
	  int numLines = boolean->GetNumberOfIntersectionLines();

	  //Objects actually don't intersect
	  if (numPts == 0 || numLines == 0)
	    {
	      std::cout<<"NO INTERSECTION"<<endl;
	    }
	  else 
	    {
	      totalIntersections++;
	      this->BooleanObject->DeepCopy(boolean->GetOutput());
	      if (this->PassInfoAsGlobal)
		this->PostSetGlobalArrays(totalIntersections);

	      checkInputArray2->InsertNextId(j);
	      for (int k = 0;k < numInputs; k++)
		{
		this->IntersectionTable[k][j] = -1;
		}
	    }
          }
        }
        this->PrintTable(numInputs);
      }
      tmp = checkInputArray;
      checkInputArray = checkInputArray2;
      checkInputArray2 = tmp;
      tmp->Reset();
    }
  return 1;
}

void vtkMultiplePolyDataIntersectionFilter::PreSetGlobalArrays(
    vtkPolyData *input)
{
  vtkSmartPointer<vtkIntArray> newPointArray = 
    vtkSmartPointer<vtkIntArray>::New();
  int numPts = input->GetNumberOfPoints();
  for (int i = 0;i < numPts; i++)
  {
    newPointArray->InsertValue(i,0);
  }
  newPointArray->SetName("GlobalBoundaryPoints");
  input->GetPointData()->AddArray(newPointArray);
}

void vtkMultiplePolyDataIntersectionFilter::PostSetGlobalArrays(
    int numIntersections)
{
  std::cout<<"Passing Data"<<endl;
  if (numIntersections == 1)
  {
    vtkSmartPointer<vtkIntArray> currentPointArray = 
      vtkSmartPointer<vtkIntArray>::New();
    currentPointArray = vtkIntArray::SafeDownCast(
	this->BooleanObject->GetPointData()->GetArray("BoundaryPoints"));
    
    currentPointArray->SetName("GlobalBoundaryPoints");
    this->BooleanObject->GetPointData()->AddArray(currentPointArray);
  }
  else
  {
    vtkSmartPointer<vtkIntArray> currentPointArray = 
      vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> globalPointArray = 
      vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> newPointArray = 
      vtkSmartPointer<vtkIntArray>::New();
    currentPointArray = vtkIntArray::SafeDownCast(
	this->BooleanObject->GetPointData()->GetArray("BoundaryPoints"));
    globalPointArray = vtkIntArray::SafeDownCast(
	this->BooleanObject->GetPointData()->
	GetArray("GlobalBoundaryPoints"));

    int numPts = this->BooleanObject->GetNumberOfPoints();
    for (int i = 0; i< numPts; i++)
    {
      newPointArray->InsertValue(i,0);
      if (globalPointArray->GetValue(i) == 1 || 
	  currentPointArray->GetValue(i) == 1)
	newPointArray->InsertValue(i,1);
    } 
    this->BooleanObject->GetPointData()->RemoveArray("GlobalBoundaryPoints");
    newPointArray->SetName("GlobalBoundaryPoints");
    this->BooleanObject->GetPointData()->AddArray(newPointArray);
  }
}

void vtkMultiplePolyDataIntersectionFilter::SetSurfaceId(
    vtkPolyData *input,int surfaceid)
{
  vtkSmartPointer<vtkIntArray> surfaceIdArray = 
    vtkSmartPointer<vtkIntArray>::New();
  int numCells = input->GetNumberOfCells();
  for (int i = 0;i < numCells; i++)
  {
    surfaceIdArray->InsertValue(i,surfaceid);
  }
  surfaceIdArray->SetName("ModelFaceID");
  input->GetCellData()->AddArray(surfaceIdArray);
}

void vtkMultiplePolyDataIntersectionFilter::PrintTable(int numInputs)
{
  std::cout<<"INTERSECTION TABLE"<<endl;
  for (int i = 0; i < numInputs; i++)
    {
    std::cout<<" ";
    for (int j = 0; j < numInputs;j++)
      {
	std::cout<<this->IntersectionTable[i][j]<<" ";
      }
    std::cout<<" "<<endl;
    }
}

//----------------------------------------------------------------------------
// This method is much too long, and has to be broken up!
// Append data sets into single polygonal data set.
int vtkMultiplePolyDataIntersectionFilter::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // get the info object
  // get the ouptut
  vtkPolyData *output = vtkPolyData::GetData(outputVector, 0);

  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  if (numInputs == 1)
    {
    output->ShallowCopy(vtkPolyData::GetData(inputVector[0], 0));
    return 1;
    }

  this->IntersectionTable = new int*[numInputs];
  vtkPolyData** inputs = new vtkPolyData*[numInputs];
  for (int idx = 0; idx < numInputs; ++idx)
    {
    inputs[idx] = vtkPolyData::GetData(inputVector[0], idx);
    this->IntersectionTable[idx] = new int[numInputs];
    for (int idy = 0; idy < numInputs; ++idy)
      {
	this->IntersectionTable[idx][idy] = -1;
      }
    }

  int intersections = this->BuildIntersectionTable(inputs, numInputs);
  if (intersections == 0)
    vtkGenericWarningMacro( << "No intersections!");
  this->PrintTable(numInputs);

  int retVal = this->ExecuteIntersection(inputs,numInputs);

  output->DeepCopy(this->BooleanObject);

  for (int idx = 0; idx < numInputs; ++idx)
    {
      delete [] this->IntersectionTable[idx];
    }
  delete [] this->IntersectionTable;
  delete [] inputs;
  return retVal;
}

//----------------------------------------------------------------------------
int vtkMultiplePolyDataIntersectionFilter::RequestUpdateExtent(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // get the output info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int piece, numPieces, ghostLevel;
  int idx;

  piece = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevel = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

  // make sure piece is valid
  if (piece < 0 || piece >= numPieces)
    {
    return 0;
    }

  int numInputs = this->GetNumberOfInputConnections(0);
  if (this->ParallelStreaming)
    {
    piece = piece * numInputs;
    numPieces = numPieces * numInputs;
    }

  vtkInformation *inInfo;
  // just copy the Update extent as default behavior.
  for (idx = 0; idx < numInputs; ++idx)
    {
    inInfo = inputVector[0]->GetInformationObject(idx);
    if (inInfo)
      {
      if (this->ParallelStreaming)
        {
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
	    piece + idx);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                    numPieces);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
                    ghostLevel);
        }
      else
        {
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
                    piece);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                    numPieces);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
                    ghostLevel);
        }
      }
    }

  return 1;
}

//----------------------------------------------------------------------------
vtkPolyData *vtkMultiplePolyDataIntersectionFilter::GetInput(int idx)
{
  return vtkPolyData::SafeDownCast(
    this->GetExecutive()->GetInputData(0, idx));
}

//----------------------------------------------------------------------------
void vtkMultiplePolyDataIntersectionFilter::PrintSelf(ostream& os, 
    vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << "ParallelStreaming:" << (this->ParallelStreaming?"On":"Off") << endl;
  os << "UserManagedInputs:" << (this->UserManagedInputs?"On":"Off") << endl;
  os << "AssignSurfaceIds:" << (this->AssignSurfaceIds?"On":"Off") << endl;
  os << "PassInfoAsGlobal:" << (this->PassInfoAsGlobal?"On":"Off") << endl;
}

//----------------------------------------------------------------------------
int vtkMultiplePolyDataIntersectionFilter::FillInputPortInformation(
    int port, vtkInformation *info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}
