#include <vtkSmartPointer.h>
#include <vtkBox.h>
#include <vtkSphere.h>
#include <vtkImplicitBoolean.h>

using namespace std;

int main()
{
	auto sphere = vtkSmartPointer<vtkSphere>::New();
	sphere->SetCenter(0., 0., 0.);
	sphere->SetRadius(20.);

	auto sphere0 = vtkSmartPointer<vtkSphere>::New();
	sphere0->SetCenter(0., 0., 0.);
	sphere0->SetRadius(10.);

	auto implicitBoolean = vtkSmartPointer<vtkImplicitBoolean>::New();
	implicitBoolean->AddFunction(sphere);
	implicitBoolean->AddFunction(sphere0);
	implicitBoolean->SetOperationTypeToDifference();

	auto a = implicitBoolean->FunctionValue(5., 5., 5.); // r = 8.66
	auto b = implicitBoolean->FunctionValue(8., 8., 8.); // r = 13.85
	auto c = implicitBoolean->FunctionValue(13., 13., 13.); // r = 22.51

	auto box = vtkSmartPointer<vtkBox>::New();
	box->SetBounds(0., 20., 0., 20., 0., 20.);

	auto a0 = box->FunctionValue(5., 5., 5.);
	auto b0 = box->FunctionValue(8., 8., 8.);
	auto c0 = box->FunctionValue(13., 13., 13.);

	auto implicitIntersection = vtkSmartPointer<vtkImplicitBoolean>::New();
	implicitIntersection->AddFunction(implicitBoolean);
	implicitIntersection->AddFunction(box);
	implicitIntersection->SetOperationTypeToIntersection();

	auto a1 = implicitIntersection->FunctionValue(5., 5., 5.); // +, +, +
	auto b1 = implicitIntersection->FunctionValue(8., 8., 8.); // +, +, +
	auto c1 = implicitIntersection->FunctionValue(13., 13., 13.); // +, +, +

	auto a2 = implicitIntersection->FunctionValue(-5., -5., -5.); // -, -, -
	auto b2 = implicitIntersection->FunctionValue(-8., -8., -8.); // -, -, -
	auto c2 = implicitIntersection->FunctionValue(-13., -13., -13.); // -, -, -

	auto a3 = implicitIntersection->FunctionValue(-5., 5., 5.); // -, +, +
	auto b3 = implicitIntersection->FunctionValue(-8., 8., 8.); // -, +, +
	auto c3 = implicitIntersection->FunctionValue(-13., 13., 13.); // -, +, +

	auto a4 = implicitIntersection->FunctionValue(5., -5., 5.); // +, -, +
	auto b4 = implicitIntersection->FunctionValue(8., -8., 8.); // +, -, +
	auto c4 = implicitIntersection->FunctionValue(13., -13., 13.); // +, -, +

	auto a5 = implicitIntersection->FunctionValue(5., 5., -5.); // +, +, -
	auto b5 = implicitIntersection->FunctionValue(8., 8., -8.); // +, +, -
	auto c5 = implicitIntersection->FunctionValue(13., 13., -13.); // +, +, -

	auto a6 = implicitIntersection->FunctionValue(-5., -5., 5.); // -, -, +
	auto b6 = implicitIntersection->FunctionValue(-8., -8., 8.); // -, -, +
	auto c6 = implicitIntersection->FunctionValue(-13., -13., 13.); // -, -, +

	auto a7 =  implicitIntersection->FunctionValue(-5., 5., -5.); // -, +, -
	auto b7 =  implicitIntersection->FunctionValue(-8., 8., -8.); // -, +, -
	auto c7 =  implicitIntersection->FunctionValue(-13., 13., -13.); // -, +, -

	auto a8 = implicitIntersection->FunctionValue(5., -5., -5.); // +, -, -
	auto b8 = implicitIntersection->FunctionValue(8., -8., -8.); // +, -, -
	auto c8 = implicitIntersection->FunctionValue(13., -13., -13.); // +, -, -

	string fileName(R"(d:\abc.obj)");
	ofstream output(fileName);
	output.precision(7);

	auto spacing = .1;
	for(auto z = -22.; z < 44.; z += spacing)
	{
		for(auto y = -22.; y < 44.; y += spacing)
		{
			for(auto x = -22.; x < 44.; x += spacing)
			{
				if(implicitIntersection->FunctionValue(x, y, z) > 0.)
				{
					continue;
				}
				output << "v " << x << " " << y << " " << z << endl;
			}
		}
	}

	output.close();

	return 0;
}
