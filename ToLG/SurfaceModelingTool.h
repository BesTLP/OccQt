#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Geom_BSplineCurve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Lin.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>
#include <gp_Pln.hxx>
#include <gp_Dir.hxx>
#include <iostream>

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>

#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <TopoDS_Face.hxx>
#include <Geom_BSplineCurve.hxx>
#include <gp_Pln.hxx>
#include <gp_Lin.hxx>

// ö�����ͣ���ʾ���ߵĲ�ͬ����
enum class CurveType
{
	LINEAR,    // ��������
	PLANAR,    // ƽ������
	NOTPLANAR, // ��ƽ������
	POINT      // ��
};

// PlanarCurve �࣬�����жϺʹ���ƽ������
class PlanarCurve
{
public:
	// Ĭ�Ϲ��캯��
	PlanarCurve();

	// �������Ĺ��캯��������һ�� B-Spline ���߾�����ݲ�ֵ
	PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);

	// ��ȡ��ǰ���ߵ�����
	CurveType GetCurveType() const { return curveType; }

	// ��ȡ�������ߣ�������������Ϊ LINEAR ʱ��Ч��
	gp_Lin GetLine() const { return line; }

	// ��ȡƽ�棨������������Ϊ PLANAR ʱ��Ч��
	gp_Pln GetPlane() const { return plane; }

	// ��ȡ�㣨������������Ϊ POINT ʱ��Ч��
	gp_Pnt GetPoint() const { return pnt; }

	// ��ȡ B-Spline ����
	Handle(Geom_BSplineCurve) GetCurve() const { return curve; }

	// ���� B-Spline ����
	void SetCurve(Handle(Geom_BSplineCurve) theCurve, double tolerance = 10)
	{
		curve = theCurve;
		IsPlanarCurve(theCurve, tolerance);
	}

private:
	CurveType curveType; // ��ǰ���ߵ�����
	gp_Lin line;         // �������ߵı�ʾ����������Ϊ����ʱ��Ч��
	gp_Pln plane;        // ƽ��ı�ʾ����������Ϊƽ������ʱ��Ч��
	gp_Pnt pnt;          // ��ı�ʾ����������Ϊ��ʱ��Ч��
	Handle(Geom_BSplineCurve) curve; // B-Spline ���ߵľ��

	// �жϸ����� B-Spline �����Ƿ�Ϊƽ������
	bool IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);

	// �жϸ����� B-Spline �����Ƿ�Ϊ��������
	bool IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);

	// �жϸ����� B-Spline �����Ƿ�Ϊһ����
	bool IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);
};

class MathTool
{
public:
	// �������� B-Spline ����֮�����С����
	static double ComputeCurveCurveDistance(const Handle(Geom_BSplineCurve)& curve, const Handle(Geom_BSplineCurve)& boundaryCurve);

	// ���������ϲ������ƽ������
	static gp_Pnt ComputeAverageSamplePoint(const Handle(Geom_BSplineCurve)& curve, int numSamples);

	// ������������֮��ļнǣ��Ի���Ϊ��λ��
	static double ComputeAngleWithAxis(const gp_Vec& vec, const gp_Vec& axis);

	// ���һ�� B-Spline �����Ƿ�����Խ�
	static void CheckSelfIntersect(std::vector<Handle(Geom_BSplineCurve)> theBSplineCurvesArray);

	// �������ߵ�ƽ������������
	static gp_Dir ComputeAverageTangent(const Handle(Geom_BSplineCurve)& curve, int numSamples);

	// ������������֮��ļнǣ��Ի���Ϊ��λ��
	static double ComputeAngleBetweenCurves(Handle(Geom_BSplineCurve)& curve1,
		Handle(Geom_BSplineCurve)& curve2,
		int numSamples = 10);

	// ���ݲο����߶� B-Spline ���߽�������
	static void SortBSplineCurves(std::vector<Handle(Geom_BSplineCurve)>& theCurves,
		Handle(Geom_BSplineCurve) referCurve);

	// ������Ҫ��ת���ߵķ�����ȷ��һ����
	static void ReverseIfNeeded(std::vector<Handle(Geom_BSplineCurve)>& curves);

	// ����㵽ƽ��ľ���
	static double ComputeDistancePointToPlane(const gp_Pnt& p, const gp_Pln& plane);

	// ����ֱ�ߺ�ƽ��֮��ļн�
	static double ComputeAngleBetweenLineAndPlane(const gp_Lin& line, const gp_Pln& plane);

	// ��������ֱ��֮��ļн�
	static double ComputeAngleBetweenLines(const gp_Lin& line1, const gp_Lin& line2);

	// ��������ƽ��֮��ļн�
	static double ComputeAngleBetweenPlanes(const gp_Pln& plane1, const gp_Pln& plane2);

	// ����㵽ֱ�ߵľ���
	static double ComputeDistancePointToLine(const gp_Pnt& point, const gp_Lin& line);

	// ��������ƽ������֮��ļн�
	static double ComputeAngleBetweenPlanarCurves(const PlanarCurve& curve1, const PlanarCurve& curve2);

	static std::vector<gp_Pnt> GetSamplePointsOnCurve(const Handle(Geom_Curve)& curve, int numPoints = 50);

	static void SortPoints(std::vector<gp_Pnt>& thePoints, const gp_Pnt& referPoint);
};

enum ReferSurfaceType
{
	GORDEN_ONE_DIRECTION_GORDEN,
	GORDEN_TWO_DIRECTION_GORDEN
};

class SurfaceModelingTool
{
public:

	static bool GetInternalCurves(
		std::vector<Handle(Geom_BSplineCurve)>& aBoundarycurveArray,
		std::vector<Handle(Geom_BSplineCurve)>& anInternalBSplineCurves,
		std::vector<Handle(Geom_BSplineCurve)>& uInternalCurve,
		std::vector<Handle(Geom_BSplineCurve)>& vInternalCurve,
		double& uAngleSum,
		double& vAngleSum,
		double AngleTolerance = 5);

	static Handle(Geom_BSplineSurface) GenerateReferSurface(
		std::vector<Handle(Geom_BSplineCurve)> aBoundarycurveArray,
		std::vector<Handle(Geom_BSplineCurve)>& uInternalCurve,
		std::vector<Handle(Geom_BSplineCurve)>& vInternalCurve,
		double uAngleSum,
		double vAngleSum,
		int isoCount,
		ReferSurfaceType referSurfaceType);
private:
	
};
