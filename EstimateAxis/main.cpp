#include <vtkOBJReader.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkMath.h>
#include <vtkVector.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/EigenValues>

#define ISEQUAL_TOL(a, b, tol) (((a)-(b)) > -(tol) && ((a)-(b)) < (tol))

using namespace std;

typedef vtkSmartPointer<vtkPolyData> PolyDataSPtr;
typedef vector<vtkVector3d> vtk3dVectors;

// https://gist.github.com/ialhashim/0a2554076a6cf32831ca
std::pair<vtkVector3d, vtkVector3d> FitBestLine(const vtk3dVectors& c)
{
	// copy coordinates to  matrix in Eigen format
	size_t num_atoms = c.size();
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic > centers(num_atoms, 3);
	for(size_t i = 0; i < num_atoms; ++i)
	{
		centers.row(i).x() = c[i].GetX();
		centers.row(i).y() = c[i].GetY();
		centers.row(i).z() = c[i].GetZ();
	}

	Eigen::Vector3d origin = centers.colwise().mean();
	Eigen::MatrixXd centered = centers.rowwise() - origin.transpose();
	Eigen::MatrixXd cov = centered.adjoint() * centered;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	Eigen::Vector3d axis = eig.eigenvectors().col(2).normalized();

	return make_pair(vtkVector3d(origin.data()), vtkVector3d(axis.data()));
}

void WriteToCSV(const vtk3dVectors& pts, const wstring& path)
{
	ofstream out(path);

	out.precision(5);
	out.setf(ios::fixed, ios::floatfield);
	for_each(pts.begin(), pts.end(), [&out](const auto& pt)
	{
		out << pt[0] << "; " << pt[1] << "; " << pt[2] << endl;
	});

	out.close();
}

pair<double, double> CalculateRange(PolyDataSPtr poly, const vtkVector3d& pos, const vtkVector3d& norm)
{
	double min = DBL_MAX, max = -1. * DBL_MAX;
	auto nPts = poly->GetNumberOfPoints();
	for(decltype(nPts) i = 0; i < nPts; ++i)
	{
		vtkVector3d pt(poly->GetPoint(i));

		pt[0] = pt[0] - pos[0];
		pt[1] = pt[1] - pos[1];
		pt[2] = pt[2] - pos[2];

		auto ip = norm.Dot(pt);

		if(ip < min)
		{
			min = ip;
		}

		if(ip > max)
		{
			max = ip;
		}
	}

	auto len =  max - min;

	min = min + 0.1 * len;
	max = max - 0.1 * len;

	return{min, max};
}

vtk3dVectors CalculateCOGs(PolyDataSPtr poly, const vtkVector3d& pos, const vtkVector3d& norm)
{
	auto range = CalculateRange(poly, pos, norm);

	double min = range.first, max = range.second;

	auto len = max - min;

	auto delta = .1;

	vtkVector3d start(pos[0] + min * norm[0], pos[1] + min * norm[1], pos[2] + min * norm[2]);

	int iter = len / delta + .5;

	vtk3dVectors cogs;
	cogs.reserve(iter);

	for(decltype(iter) i = 0; i < iter; ++i)
	{
		auto aPlane = vtkSmartPointer<vtkPlane>::New();
		aPlane->SetOrigin(start.GetData());
		aPlane->SetNormal(norm.GetData()[0], norm.GetData()[1], norm.GetData()[2]);

		auto aCutter = vtkSmartPointer<vtkCutter>::New();
		aCutter->SetCutFunction(aPlane);
		aCutter->SetInputData(poly);
		aCutter->Update();

		auto out = aCutter->GetOutput();

		vtkVector3d cog(0, 0, 0);
		auto nPts = out->GetNumberOfPoints();
		for(decltype(nPts) j = 0; j < nPts; ++j)
		{
			auto pt = out->GetPoint(j);
			cog[0] += pt[0];
			cog[1] += pt[1];
			cog[2] += pt[2];
		}

		cog[0] /= nPts;
		cog[1] /= nPts;
		cog[2] /= nPts;

		cogs.push_back(cog);

		start[0] = start[0] + delta * norm[0];
		start[1] = start[1] + delta * norm[1];
		start[2] = start[2] + delta * norm[2];
	}

	return cogs;
}

int main()
{
	//auto objFile = string(R"(F:\Vc++PRJs\VTKTests\test_prjs\x64\Debug\cylinder.obj)");
	//auto objFile = string(R"(F:\Vc++PRJs\VTKTests\test_prjs\x64\Debug\cone.obj)");
	auto objFile = string(R"(F:\Vc++PRJs\VTKTests\test_prjs\x64\Debug\vase2.obj)");

	auto objReader = vtkSmartPointer<vtkOBJReader>::New();
	objReader->SetFileName(objFile.c_str());
	objReader->Update();

	auto cylinder = objReader->GetOutput();

	vtkVector3d pos(0, 0, 0);
	vtkVector3d norm(-0.061628, 0.996195, 0.061628); // (0, 1, 0) rotated by 5° about (0.707106, 0, 0.707106)
	//vtkVector3d norm(0, 1, 0);

	int iter = 0;
	while(true)
	{
		auto cogs = CalculateCOGs(cylinder, pos, norm);

		auto line = FitBestLine(cogs);

		wstringstream ss;
		ss << L"pts" << iter++ << ".txt";
		WriteToCSV(cogs, ss.str());

		auto tol = 1e-5;
		auto ip = norm.Dot(line.second);

		if(ISEQUAL_TOL(ip, 1, tol) || ISEQUAL_TOL(ip, -1, tol))
		{
			break;
		}

		pos = line.first;
		norm = line.second;
	};

	return 0;
}
