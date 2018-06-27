#include <vtkOBJReader.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkSTLWriter.h>

using namespace std;

int main()
{
	auto objfile = string(R"(D:\projects\ZV_WORK\src\work\Production\msvs\Zed\subPoly.obj)");
	//auto objfile = string(R"(d:\sphere.obj)");

	auto objreader = vtkSmartPointer<vtkOBJReader>::New();
	objreader->SetFileName(objfile.c_str());
	objreader->Update();

	auto objpoly = objreader->GetOutput();

	auto linearextrusion = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
	linearextrusion->SetInputData(objpoly);
	linearextrusion->SetCapping(1);
	linearextrusion->SetExtrusionTypeToPointExtrusion();
	//linearextrusion->SetExtrusionPoint(0., 0., 0.);
	linearextrusion->SetExtrusionPoint(91.350445360311966, 165.46352360347873, -93.616526598563397);
	//linearextrusion->SetScaleFactor(-1. * .1);
	linearextrusion->SetScaleFactor(.1);
	linearextrusion->Update();

	auto outputlinearextrusion = linearextrusion->GetOutput();

	auto stlwriter = vtkSmartPointer<vtkSTLWriter>::New();
	stlwriter->SetFileName(R"(d:\linextsphere.stl)");
	stlwriter->SetInputData(outputlinearextrusion);
	stlwriter->Update();

	return 0;
}

