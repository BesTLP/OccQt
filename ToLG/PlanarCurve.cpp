#include "PlanarCurve.h"
#include "MathTool.h"
PlanarCurve::PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance)
	: curveType(CurveType::NOTPLANAR), curve(theCurve), line(), plane()
{
	IsPlanarCurve(curve);
}

PlanarCurve::PlanarCurve()
	:curveType(CurveType::NOTPLANAR), curve(), line(), plane()
{
	curveType = CurveType::NOTPLANAR;
}

bool PlanarCurve::IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	// 检查曲线是否为空
	if (theCurve.IsNull())
	{
		std::cerr << "输入的曲线为空。" << std::endl;
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	bool isLinear = IsBSplineCurveLinear(theCurve);
	if (isLinear)
	{
		curveType = CurveType::LINEAR;
		line = gp_Lin(theCurve->StartPoint(), gp_Vec(theCurve->StartPoint(), theCurve->EndPoint()));
		return true;
	}

	bool isPoint = IsBSplineCurvePoint(theCurve);
	if (isPoint)
	{
		curveType = CurveType::POINT;
		return true;
	}

	// 采样曲线上的点
	std::vector<gp_Pnt> boundary_Sampling;
	int numSamples = 100; // 采样点数量，可根据需要调整
	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();
	Standard_Real step = (lastParam - firstParam) / (numSamples - 1);

	boundary_Sampling.reserve(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		Standard_Real param = firstParam + i * step;
		gp_Pnt pnt;
		theCurve->D0(param, pnt); // 获取曲线上的点
		boundary_Sampling.push_back(pnt);
	}

	// 检查是否成功采样
	if (boundary_Sampling.empty())
	{
		std::cerr << "曲线采样失败，没有采样到任何点。" << std::endl;
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	// 声明法向量和坐标轴
	gp_XYZ N, X, Y;
	// 计算质心 P	
	gp_Pnt P;
	P.SetCoord(0, 0, 0);
	int num = boundary_Sampling.size();
	for (const auto& point : boundary_Sampling)
	{
		P.SetCoord(P.X() + point.X(), P.Y() + point.Y(), P.Z() + point.Z());
	}
	P.SetCoord(P.X() / num, P.Y() / num, P.Z() / num);

	// 替换协方差矩阵构建和特征值分解部分
	int m = 3;
	Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(m, m);

	double c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
	gp_Pnt point; // 声明 'point'

	for (int i = 0; i < num; i++)
	{
		point = boundary_Sampling[i];

		c00 += (point.X() - P.X()) * (point.X() - P.X());
		c01 += (point.X() - P.X()) * (point.Y() - P.Y());
		c02 += (point.X() - P.X()) * (point.Z() - P.Z());

		c11 += (point.Y() - P.Y()) * (point.Y() - P.Y());
		c12 += (point.Y() - P.Y()) * (point.Z() - P.Z());

		c22 += (point.Z() - P.Z()) * (point.Z() - P.Z());
	}

	A1(0, 0) = c00;
	A1(0, 1) = c01;
	A1(0, 2) = c02;

	A1(1, 0) = c01;
	A1(1, 1) = c11;
	A1(1, 2) = c12;

	A1(2, 0) = c02;
	A1(2, 1) = c12;
	A1(2, 2) = c22;


	Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A1);
	Eigen::VectorXcd E1 = eigensolver.eigenvalues();
	auto  E2 = eigensolver.eigenvectors();

	auto eigen1 = E1.col(0)[0];
	auto eigen2 = E1.col(0)[1];
	auto eigen3 = E1.col(0)[2];

	double eigenvalue1 = E1.col(0)[0].real();
	double eigenvalue2 = E1.col(0)[1].real();
	double eigenvalue3 = E1.col(0)[2].real();

	// 根据最小特征值确定法向量和坐标轴
	if (eigenvalue1 < eigenvalue2 && eigenvalue1 < eigenvalue3)
	{
		N.SetCoord(E2.col(0)[0].real(), E2.col(0)[1].real(), E2.col(0)[2].real());
		X.SetCoord(E2.col(1)[0].real(), E2.col(1)[1].real(), E2.col(1)[2].real());
		Y.SetCoord(E2.col(2)[0].real(), E2.col(2)[1].real(), E2.col(2)[2].real());
	}
	else if (eigenvalue2 < eigenvalue1 && eigenvalue2 < eigenvalue3)
	{
		N.SetCoord(E2.col(1)[0].real(), E2.col(1)[1].real(), E2.col(1)[2].real());
		Y.SetCoord(E2.col(0)[0].real(), E2.col(0)[1].real(), E2.col(0)[2].real());
		X.SetCoord(E2.col(2)[0].real(), E2.col(2)[1].real(), E2.col(2)[2].real());
	}
	else
	{
		N.SetCoord(E2.col(2)[0].real(), E2.col(2)[1].real(), E2.col(2)[2].real());
		Y.SetCoord(E2.col(1)[0].real(), E2.col(1)[1].real(), E2.col(1)[2].real());
		X.SetCoord(E2.col(0)[0].real(), E2.col(0)[1].real(), E2.col(0)[2].real());
	}

	// 创建 gp_Dir 对象用于 gp_pln
	gp_Dir dirN(N);
	plane = gp_Pln(P, dirN);
	// 评估平面性
	double maxDistance = 0.0;
	double sumDistance = 0.0;
	for (const auto& point : boundary_Sampling)
	{
		double distance = MathTool::ComputeDistancePointToPlane(point, plane);
		sumDistance += distance;
		if (distance > maxDistance)
		{
			maxDistance = distance;
		}
	}

	double averageDistance = sumDistance / boundary_Sampling.size();

	// 设置平面性阈值（根据具体需求调整）
	if (maxDistance < tolerance)
	{
		curveType = CurveType::PLANAR;
		return true;
	}
	else
	{
		curveType = CurveType::NOTPLANAR;
		return false;
	}
}

bool PlanarCurve::IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	if (theCurve.IsNull())
	{
		std::cerr << "The B-Spline curve is null." << std::endl;
		return false;
	}

	// 获取控制点数量
	Standard_Integer numPoles = theCurve->NbPoles();
	if (numPoles < 2)
	{
		// 少于两个控制点，无法确定直线
		std::cerr << "The B-Spline curve has fewer than two control points." << std::endl;
		return false;
	}

	// 获取第一个和第二个控制点
	gp_Pnt P0 = theCurve->Pole(1); // 注意：OpenCASCADE的索引从1开始
	gp_Pnt P1 = theCurve->Pole(2);

	// 计算方向向量 v = P1 - P0
	gp_Vec v(P0, P1);
	if (v.Magnitude() < tolerance)
	{
		// 前两个点重合，无法定义方向
		std::cerr << "The first two control points are coincident; cannot define a direction." << std::endl;
		return false;
	}

	// 对于每一个后续的控制点，检查是否与方向向量共线
	for (Standard_Integer i = 3; i <= numPoles; ++i)
	{
		gp_Pnt Pi = theCurve->Pole(i);
		gp_Vec u(P0, Pi);

		// 计算叉积 v × u
		gp_Vec cross = v.Crossed(u);

		// 检查叉积的模长是否在容差范围内
		if (cross.Magnitude() > tolerance)
		{
			// 不共线
			return false;
		}
	}
	// 所有控制点共线
	return true;
}

bool PlanarCurve::IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	// 获取控制点数量
	Standard_Integer numPoles = theCurve->NbPoles();
	// 获取第一个和第二个控制点
	gp_Pnt P0 = theCurve->Pole(1); // 注意：OpenCASCADE的索引从1开始
	gp_Pnt P1 = theCurve->Pole(2);

	// 计算方向向量 v = P1 - P0
	gp_Vec v(P0, P1);
	if (numPoles == 2 && v.Magnitude() < tolerance)
	{
		return true;
	}
	return false;
}