#include "PlanarCurve.h"

PlanarCurve::PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance)
	: curveType(CurveType::NOTPLANAR), curve(theCurve), line(), plane()
{
	IsPlanarCurve(curve, 50);
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
		double distance = ComputeDistancePointToPlane(point, plane);
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

double PlanarCurve::ComputeAngleBetweenPlanarCurves(const PlanarCurve& curve1, const PlanarCurve& curve2)
{
	double angle = INT_MAX;
	double distance = INT_MAX;
	if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::LINEAR)
	{
		// ����ֱ��ֱ����
		angle = ComputeAngleBetweenLines(curve1.GetLine(), curve2.GetLine());
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::PLANAR)
	{
		// ����ƽ��ƽ����
		angle = ComputeAngleBetweenPlanes(curve1.plane, curve2.plane);
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::LINEAR)
	{
		// ����ƽ���ֱ����
		angle = ComputeAngleBetweenLineAndPlane(curve2.line, curve1.plane);
	}
	else if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::PLANAR)
	{
		// ����ƽ���ֱ����
		angle = ComputeAngleBetweenLineAndPlane(curve1.line, curve2.plane);
	}
	else if (curve1.GetCurveType() == CurveType::POINT && curve2.GetCurveType() == CurveType::POINT)
	{
		distance = curve1.pnt.Distance(curve2.pnt);
	}
	else if (curve1.GetCurveType() == CurveType::POINT && curve2.GetCurveType() == CurveType::PLANAR)
	{
		distance = ComputeDistancePointToPlane(curve1.pnt, curve2.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::POINT && curve2.GetCurveType() == CurveType::LINEAR)
	{
		distance = ComputeDistancePointToLine(curve1.pnt, curve2.GetLine());
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::POINT)
	{
		distance = ComputeDistancePointToPlane(curve2.pnt, curve1.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::POINT)
	{
		distance = ComputeDistancePointToLine(curve2.pnt, curve1.GetLine());
	}
	
	if (distance < 1e-3) angle = 0.0;

	return angle;
}

double PlanarCurve::ComputeAngleBetweenLineAndPlane(const gp_Lin& line, const gp_Pln& plane)
{
	// ��ȡֱ�ߵķ�������
	gp_Vec lineDirection = line.Direction();  // ֱ�߷�������

	// ��ȡƽ��ķ�����
	gp_Vec planeNormal = plane.Axis().Direction();  // ƽ�淨����

	// ����ֱ�߷���������ƽ�淨�����ĵ��
	double dotProduct = lineDirection.Dot(planeNormal);

	// ����ֱ�߷��������ͷ�������ģ��
	double lineMagnitude = lineDirection.Magnitude();
	double planeNormalMagnitude = planeNormal.Magnitude();

	// ����нǣ����ȣ�
	double angle = std::acos(dotProduct / (lineMagnitude * planeNormalMagnitude));

	return angle;
}

double PlanarCurve::ComputeAngleBetweenLines(const gp_Lin& line1, const gp_Lin& line2)
{
	// ��ȡ����ֱ�ߵķ�������
	gp_Vec direction1 = line1.Direction();
	gp_Vec direction2 = line2.Direction();

	// ���㷽�������ĵ��
	double dotProduct = direction1.Dot(direction2);

	// ���㷽��������ģ��
	double magnitude1 = direction1.Magnitude();
	double magnitude2 = direction2.Magnitude();

	// ����нǣ����ȣ�
	double angle = std::acos(dotProduct / (magnitude1 * magnitude2));

	return angle;
}

double PlanarCurve::ComputeAngleBetweenPlanes(const gp_Pln& plane1, const gp_Pln& plane2)
{
	// ��ȡƽ��ķ�����
	gp_Dir normal1 = plane1.Axis().Direction();
	gp_Dir normal2 = plane2.Axis().Direction();

	// ���㷨����֮��ĵ��
	double dotProduct = normal1.Dot(normal2);

	// ���� dotProduct �ķ�Χ�� [-1, 1] ֮�䣬�Է�ֹ��ֵ���µ� acos �������
	dotProduct = std::max(-1.0, std::min(1.0, dotProduct));

	// ����нǵĻ���ֵ
	double angleRad = std::acos(dotProduct);

	// ������ת��Ϊ����
	double angleDeg = angleRad * 180.0 / M_PI;

	// ��ȡ��ǣ�0�㵽90�㣩
	if (angleDeg > 90.0)
	{
		angleDeg = 180.0 - angleDeg;
	}

	return angleDeg;
}

double PlanarCurve::ComputeDistancePointToLine(const gp_Pnt& point, const gp_Lin& line)
{
	// ֱ�ߵķ�������
	gp_Vec lineDirection = line.Direction();

	// �㵽ֱ���ϵ�����һ�㣨����ѡ��ֱ���ϵ�һ������Ϊ�ο��㣩
	gp_Pnt linePoint = line.Location(); // ֱ���ϵ�һ����

	// ����ֱ���ϵĵ������
	gp_Vec pointToLine = point.XYZ() - linePoint.XYZ();

	// ����㵽ֱ�ߵĴ�ֱ����
	gp_Vec crossProduct = pointToLine.Crossed(lineDirection);

	// ���ؾ��룬ʹ�ò�˵�ģ�����Է���������ģ��
	return crossProduct.Magnitude() / lineDirection.Magnitude();
}

double PlanarCurve::ComputeDistancePointToPlane(const gp_Pnt& p, const gp_Pln& plane)
{
	// ��ȡƽ��ķ������͵�
	gp_Dir normal = plane.Axis().Direction();
	gp_Pnt planeP = plane.Location();

	// ���� (p - planeP)
	gp_Vec vec(p.X() - planeP.X(), p.Y() - planeP.Y(), p.Z() - planeP.Z());

	// �㵽ƽ��ľ��� = |vec . normal| / |normal|
	// ���� normal �ǵ�λ��������ĸΪ1
	double distance = std::abs(vec.Dot(normal));
	return distance;
}
