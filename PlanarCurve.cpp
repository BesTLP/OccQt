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
	// ��������Ƿ�Ϊ��
	if (theCurve.IsNull())
	{
		std::cerr << "���������Ϊ�ա�" << std::endl;
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

	// ���������ϵĵ�
	std::vector<gp_Pnt> boundary_Sampling;
	int numSamples = 100; // �������������ɸ�����Ҫ����
	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();
	Standard_Real step = (lastParam - firstParam) / (numSamples - 1);

	boundary_Sampling.reserve(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		Standard_Real param = firstParam + i * step;
		gp_Pnt pnt;
		theCurve->D0(param, pnt); // ��ȡ�����ϵĵ�
		boundary_Sampling.push_back(pnt);
	}

	// ����Ƿ�ɹ�����
	if (boundary_Sampling.empty())
	{
		std::cerr << "���߲���ʧ�ܣ�û�в������κε㡣" << std::endl;
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	// ������������������
	gp_XYZ N, X, Y;
	// �������� P	
	gp_Pnt P;
	P.SetCoord(0, 0, 0);
	int num = boundary_Sampling.size();
	for (const auto& point : boundary_Sampling)
	{
		P.SetCoord(P.X() + point.X(), P.Y() + point.Y(), P.Z() + point.Z());
	}
	P.SetCoord(P.X() / num, P.Y() / num, P.Z() / num);

	// �滻Э������󹹽�������ֵ�ֽⲿ��
	int m = 3;
	Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(m, m);

	double c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
	gp_Pnt point; // ���� 'point'

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

	// ������С����ֵȷ����������������
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

	// ���� gp_Dir �������� gp_pln
	gp_Dir dirN(N);
	plane = gp_Pln(P, dirN);
	// ����ƽ����
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

	// ����ƽ������ֵ�����ݾ������������
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

	// ��ȡ���Ƶ�����
	Standard_Integer numPoles = theCurve->NbPoles();
	if (numPoles < 2)
	{
		// �����������Ƶ㣬�޷�ȷ��ֱ��
		std::cerr << "The B-Spline curve has fewer than two control points." << std::endl;
		return false;
	}

	// ��ȡ��һ���͵ڶ������Ƶ�
	gp_Pnt P0 = theCurve->Pole(1); // ע�⣺OpenCASCADE��������1��ʼ
	gp_Pnt P1 = theCurve->Pole(2);

	// ���㷽������ v = P1 - P0
	gp_Vec v(P0, P1);
	if (v.Magnitude() < tolerance)
	{
		// ǰ�������غϣ��޷����巽��
		std::cerr << "The first two control points are coincident; cannot define a direction." << std::endl;
		return false;
	}

	// ����ÿһ�������Ŀ��Ƶ㣬����Ƿ��뷽����������
	for (Standard_Integer i = 3; i <= numPoles; ++i)
	{
		gp_Pnt Pi = theCurve->Pole(i);
		gp_Vec u(P0, Pi);

		// ������ v �� u
		gp_Vec cross = v.Crossed(u);

		// �������ģ���Ƿ����ݲΧ��
		if (cross.Magnitude() > tolerance)
		{
			// ������
			return false;
		}
	}
	// ���п��Ƶ㹲��
	return true;
}

bool PlanarCurve::IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	// ��ȡ���Ƶ�����
	Standard_Integer numPoles = theCurve->NbPoles();
	// ��ȡ��һ���͵ڶ������Ƶ�
	gp_Pnt P0 = theCurve->Pole(1); // ע�⣺OpenCASCADE��������1��ʼ
	gp_Pnt P1 = theCurve->Pole(2);

	// ���㷽������ v = P1 - P0
	gp_Vec v(P0, P1);
	if (numPoles == 2 && v.Magnitude() < tolerance)
	{
		return true;
	}
	return false;
}