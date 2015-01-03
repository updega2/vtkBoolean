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
#include "vtkBooleanOperationPolyDataFilterMine.h"
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
  this->OutputPointsPrecision = vtkAlgorithm::DEFAULT_PRECISION;
  this->NoIntersectionOutput = 1;
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
      "SetInputConnectionByNumber is not supported if UserManagedInputs is false");
    return;
    }

  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, num, input);
}

int vtkMultiplePolyDataIntersectionFilter::BuildIntersectionTable(vtkPolyData* inputs[], int numInputs)
{
  for (int i = 0;i < numInputs;i++)
    {
    inputs[i]->ComputeBounds();
    }

  int totalIntersections=0;
  for (int i = 0;i < numInputs;i++)
    {
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

int vtkMultiplePolyDataIntersectionFilter::ExecuteIntersection(vtkPolyData* inputs[], int numInputs)
{                                            
  int numChecks = 0;  
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
        if (this->IntersectionTable[i][j] == 1)
          {
	  for (int k = 0;k < numInputs; k++)
	    {
	    this->IntersectionTable[k][i] = -1;
	    }
#ifdef USE_MINE
	  vtkSmartPointer<vtkBooleanOperationPolyDataFilterMine> boolean = 
	    vtkSmartPointer<vtkBooleanOperationPolyDataFilterMine>::New();
#else
	  vtkSmartPointer<vtkBooleanOperationPolyDataFilter> boolean = 
	    vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
#endif

	  boolean->SetInputData(0,this->BooleanObject);
	  boolean->SetInputData(1,inputs[j]);
	  boolean->SetOperationToUnion();
	  boolean->Update();

	  this->BooleanObject->DeepCopy(boolean->GetOutput());
	  checkInputArray2->InsertNextId(j);
          }
        }
      }
      tmp = checkInputArray;
      checkInputArray = checkInputArray2;
      checkInputArray2 = tmp;
      tmp->Reset();
    }
  return 1;
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
int vtkMultiplePolyDataIntersectionFilter::RequestData(vtkInformation *vtkNotUsed(request),
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
int vtkMultiplePolyDataIntersectionFilter::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector)
{
  // get the output info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int piece, numPieces, ghostLevel;
  int idx;

  piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevel = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

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
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
                    piece + idx);
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                    numPieces);
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
                    ghostLevel);
        }
      else
        {
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
                    piece);
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                    numPieces);
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
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
void vtkMultiplePolyDataIntersectionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << "ParallelStreaming:" << (this->ParallelStreaming?"On":"Off") << endl;
  os << "UserManagedInputs:" << (this->UserManagedInputs?"On":"Off") << endl;
  os << indent << "Output Points Precision: " << this->OutputPointsPrecision
     << endl;
}

//----------------------------------------------------------------------------
int vtkMultiplePolyDataIntersectionFilter::FillInputPortInformation(int port, vtkInformation *info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}
