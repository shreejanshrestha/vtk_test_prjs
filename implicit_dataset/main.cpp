#include <vtkOBJReader.h>
#include <vtkImplicitPolyDataDistance.h>

using namespace std;

int main()
{
	auto objFile = string(R"(d:\sphere.obj)");

	auto objReader = vtkSmartPointer<vtkOBJReader>::New();
	objReader->SetFileName(objFile.c_str());
	objReader->Update();

	auto outObjReader = objReader->GetOutput();

	auto implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicitPolyDataDistance->SetInput(outObjReader);

	for(int i = 0; i < 291600; ++i)
	{
		auto v0 = implicitPolyDataDistance->FunctionValue(.0, .0, .0);
	}

	return 0;
}
