cmake_minimum_required(VERSION 2.8)
 
PROJECT(Intersect)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

option(WITH_MINE "" OFF)
option(NEW_SPLIT "" OFF)
option(EAR_CLIP "" OFF)
 
add_executable(Intersect Intersect.cxx 
			 vtkIntersectionPolyDataFilterMine.cxx)

add_executable(Boolean Boolean.cxx 
			 vtkIntersectionPolyDataFilterMine.cxx
			 vtkBooleanOperationPolyDataFilterMine.cxx)

add_executable(VesselIntersection VesselIntersection.cxx 
			 vtkIntersectionPolyDataFilterMine.cxx
			 vtkBooleanOperationPolyDataFilterMine.cxx
			 vtkMultiplePolyDataIntersectionFilter.cxx)
 
if(VTK_LIBRARIES)
	target_link_libraries(Intersect ${VTK_LIBRARIES})
	target_link_libraries(Boolean ${VTK_LIBRARIES})
	target_link_libraries(VesselIntersection ${VTK_LIBRARIES})
else()
	target_link_libraries(Intersect vtkHybrid)
	target_link_libraries(Boolean vtkHybrid)
	target_link_libraries(VesselIntersection vtkHybrid)
endif()

if(WITH_MINE)
	add_definitions(-DUSE_MINE)
endif()
if (NEW_SPLIT)
	add_definitions(-DNEW_SPLIT)
endif()
if (EAR_CLIP)
  	add_definitions(-DEAR_CLIP)
endif()
