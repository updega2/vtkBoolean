//
//  TetgenInterface.cxx
//  
//
//  Created by Adam Updegrove on 1/31/14.
//
//

/*=========================================================================
 
 Program:   SimVascular
 
 This is a program to read in a solid file, extract the boundaries, call 
 Tetgen, and output the necessary files for the SimVascular presolver

 =========================================================================*/

#include "vtkSTLReader.h"
#include "vtkOBBTree.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkDistancePolyDataFilter.h"
#include "vtkIntersectionPolyDataFilterMine.h"
#include "vtkBooleanOperationPolyDataFilterMine.h"
#include "vtkIntersectionPolyDataFilter.h"
#include "vtkBooleanOperationPolyDataFilter.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkIntArray.h"
#include "vtkGetBoundaryFaces.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkDataWriter.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkFillHolesFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"

#include <string>
#include <sstream>
#include <iostream>

//Function to turn an integer into a string
std::string intToString(int i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

//Function to get the directory from the input File Name
//For example, /User/Adam.stl returns /User
std::string getPath(std::string fullName)
{
  std::string pathName;
  unsigned split = fullName.find_last_of("/\\");
  pathName = fullName.substr(0,split);
  return pathName;
}

//Function to get the raw file name from the input File name
//For example, Adam.stl returns Adam
std::string getRawName(std::string fullName)
{
  std::string rawName;
  unsigned split = fullName.find_last_of("/\\");
  rawName = fullName.substr(split+1);
  rawName.erase(rawName.find_last_of("."),std::string::npos);
  return rawName;
}

//Function to read in the STL file, extract the boundaries and pass the input 
//Poly Data information 
void ReadSTLFile(std::string inputFilename, vtkPolyData *polydata)
{
  //Create an STL reader for reading the file
  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();

  //Save the output information from the boundary filter to a Poly Data 
  //structure
  polydata->DeepCopy(reader->GetOutput());
  polydata->BuildLinks();
}

//Function to write the polydata to a vtp
void WriteVTPFile(std::string inputFilename,vtkPolyData *writePolyData,std::string attachName)
{
  std::string rawName, pathName, outputFilename;

  vtkSmartPointer<vtkXMLPolyDataWriter> writer  = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  
  pathName = getPath(inputFilename);
  rawName = getRawName(inputFilename);

  outputFilename = pathName+"/"+rawName+attachName+".vtp";

  writer->SetFileName(outputFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(writePolyData);
#else
  writer->SetInputData(writePolyData);
#endif
  //writer->SetDataModeToAscii();

  writer->Write();
}

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
      std::cout << "Input Filenames Required, Need Boolean!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];
  int op = (int) (argv[3][0] - '0');
                                                                             
  //Create pointers for reading the STL, creating the full unstructured grid,
  //creating the full poly data, and create the region poly data sets
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> outBoolean = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> outNewBoolean = vtkSmartPointer<vtkPolyData>::New();
#ifdef USE_MINE
  vtkSmartPointer<vtkIntersectionPolyDataFilterMine> PolyDataIntersection = 
	  vtkSmartPointer<vtkIntersectionPolyDataFilterMine>::New();
  vtkSmartPointer<vtkDistancePolyDataFilter> PolyDataDistance =
    vtkSmartPointer<vtkDistancePolyDataFilter>::New();

  vtkSmartPointer<vtkBooleanOperationPolyDataFilterMine> myBoolean = 
	  vtkSmartPointer<vtkBooleanOperationPolyDataFilterMine>::New();
#else
  vtkSmartPointer<vtkIntersectionPolyDataFilter> PolyDataIntersection = 
	  vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
  vtkSmartPointer<vtkDistancePolyDataFilter> PolyDataDistance =
    vtkSmartPointer<vtkDistancePolyDataFilter>::New();

  vtkSmartPointer<vtkBooleanOperationPolyDataFilter> myBoolean = 
	  vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
#endif

  vtkSmartPointer<vtkCleanPolyData> cleaner =
	  vtkSmartPointer<vtkCleanPolyData>::New();
  vtkSmartPointer<vtkCleanPolyData> cleaner2 =
	  vtkSmartPointer<vtkCleanPolyData>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  ReadSTLFile(inputFilename1,pd1); 
  ReadSTLFile(inputFilename2,pd2); 

 // std::cout<<"Performing Intersection..."<<endl;
 // PolyDataIntersection->AddInputData(0,pd1);
 // PolyDataIntersection->AddInputData(1,pd2);
 // PolyDataIntersection->SplitFirstOutputOn();
 // PolyDataIntersection->SplitSecondOutputOn();
 // PolyDataIntersection->Update();

 // std::cout<<"Done...Writing Files..."<<endl;
 // WriteVTPFile(inputFilename1,PolyDataIntersection->GetOutput(0),"_IntersectionLines");
 // WriteVTPFile(inputFilename1,PolyDataIntersection->GetOutput(1),"_IntersectionObject1");
 // WriteVTPFile(inputFilename2,PolyDataIntersection->GetOutput(2),"_IntersectionObject2");
 // std::cout<<"Done"<<endl;

//  std::cout<<"Getting Regions 1..."<<endl;
//  getSurfaces->SetInputData(0,PolyDataIntersection->GetOutput(1));
//  getSurfaces->SetInputData(1,PolyDataIntersection->GetOutput(0));
//  getSurfaces->Update();
//  std::cout<<"Got Regions 1..."<<endl;
//  std::cout<<"Getting Regions 2..."<<endl;
//  getSurfaces2->SetInputData(0,PolyDataIntersection->GetOutput(2));
//  getSurfaces2->SetInputData(1,PolyDataIntersection->GetOutput(0));
//  getSurfaces2->Update();
//  std::cout<<"Got Regions 2..."<<endl;
//
//  std::cout<<"Done...Writing Files..."<<endl;
//  WriteVTPFile(inputFilename1,getSurfaces->GetOutput(),"_Surf1WithRegions");
//  WriteVTPFile(inputFilename1,getSurfaces2->GetOutput(),"_Surf2WithRegions");
//  std::cout<<"Done"<<endl;
//  std::cout<<"Getting Regions 1..."<<endl;
//  othergetSurfaces->SetInputData(0,PolyDataIntersection->GetOutput(1));
//  othergetSurfaces->SetInputData(1,PolyDataIntersection->GetOutput(0));
//  othergetSurfaces->Update();
//  std::cout<<"Got Regions 1..."<<endl;
//  std::cout<<"Getting Regions 2..."<<endl;
//  othergetSurfaces2->SetInputData(0,PolyDataIntersection->GetOutput(2));
//  othergetSurfaces2->SetInputData(1,PolyDataIntersection->GetOutput(0));
//  othergetSurfaces2->Update();
//  std::cout<<"Got Regions 2..."<<endl;
//
//  std::cout<<"Done...Writing Files..."<<endl;
//  WriteVTPFile(inputFilename1,othergetSurfaces->GetOutput(),"_Surf1WithRegions");
//  WriteVTPFile(inputFilename1,othergetSurfaces2->GetOutput(),"_Surf2WithRegions");
//  std::cout<<"Done"<<endl;

  myBoolean->SetInputData(0,pd1);
  myBoolean->SetInputData(1,pd2);

  if (op == 0)
    myBoolean->SetOperationToUnion();
  else if (op == 1)
    myBoolean->SetOperationToIntersection();
  else if (op == 2)
    myBoolean->SetOperationToDifference();

  myBoolean->Update();
  vtkSmartPointer<vtkPolyData> fullpd = 
	  vtkSmartPointer<vtkPolyData>::New();
  fullpd->DeepCopy(myBoolean->GetOutput());

  fullpd->GetCellData()->RemoveArray("BadTri");
  fullpd->GetCellData()->RemoveArray("FreeEdge");

  vtkIntersectionPolyDataFilterMine::CleanAndCheckSurface(fullpd);
  double fullbadtri[2], fullfreeedge[2];
  fullpd->GetCellData()->GetArray("BadTri")->GetRange(fullbadtri,0);
  fullpd->GetCellData()->GetArray("FreeEdge")->GetRange(fullfreeedge,0);

  std::cout<<"FULL SURFACE BAD TRI MIN: "<<fullbadtri[0]<<" MAX: "<<fullbadtri[1]<<endl;
  std::cout<<"FULL SURFACE FREE EDGE MIN: "<<fullfreeedge[0]<<" MAX: "<<fullfreeedge[1]<<endl;


  WriteVTPFile(inputFilename1,fullpd,"_FullBoolean");
  //Exit the program without errors
  return EXIT_SUCCESS;
}




