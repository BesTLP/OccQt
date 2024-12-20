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

// 枚举类型，表示曲线的不同类型
enum class CurveType
{
	LINEAR,    // 线性曲线
	PLANAR,    // 平面曲线
	NOTPLANAR, // 非平面曲线
	POINT      // 点
};

// PlanarCurve 类，用于判断和处理平面曲线
class PlanarCurve
{
public:
	// 默认构造函数
	PlanarCurve();

	// 带参数的构造函数，接受一个 B-Spline 曲线句柄和容差值
	PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);

	// 获取当前曲线的类型
	CurveType GetCurveType() const { return curveType; }

	// 获取线性曲线（仅当曲线类型为 LINEAR 时有效）
	gp_Lin GetLine() const { return line; }

	// 获取平面（仅当曲线类型为 PLANAR 时有效）
	gp_Pln GetPlane() const { return plane; }

	// 获取点（仅当曲线类型为 POINT 时有效）
	gp_Pnt GetPoint() const { return pnt; }

	// 获取 B-Spline 曲线
	Handle(Geom_BSplineCurve) GetCurve() const { return curve; }

	// 设置 B-Spline 曲线
	void SetCurve(Handle(Geom_BSplineCurve) theCurve, double tolerance = 10)
	{
		curve = theCurve;
		IsPlanarCurve(theCurve, tolerance);
	}

private:
	CurveType curveType; // 当前曲线的类型
	gp_Lin line;         // 线性曲线的表示（仅当曲线为线性时有效）
	gp_Pln plane;        // 平面的表示（仅当曲线为平面曲线时有效）
	gp_Pnt pnt;          // 点的表示（仅当曲线为点时有效）
	Handle(Geom_BSplineCurve) curve; // B-Spline 曲线的句柄

	// 判断给定的 B-Spline 曲线是否为平面曲线
	bool IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);

	// 判断给定的 B-Spline 曲线是否为线性曲线
	bool IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);

	// 判断给定的 B-Spline 曲线是否为一个点
	bool IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);
};

class MathTool
{
public:
	// 计算两条 B-Spline 曲线之间的最小距离
	static double ComputeCurveCurveDistance(const Handle(Geom_BSplineCurve)& curve, const Handle(Geom_BSplineCurve)& boundaryCurve);

	// 计算曲线上采样点的平均坐标
	static gp_Pnt ComputeAverageSamplePoint(const Handle(Geom_BSplineCurve)& curve, int numSamples);

	// 计算向量与轴之间的夹角（以弧度为单位）
	static double ComputeAngleWithAxis(const gp_Vec& vec, const gp_Vec& axis);

	// 检查一组 B-Spline 曲线是否存在自交
	static void CheckSelfIntersect(std::vector<Handle(Geom_BSplineCurve)> theBSplineCurvesArray);

	// 计算曲线的平均切向量方向
	static gp_Dir ComputeAverageTangent(const Handle(Geom_BSplineCurve)& curve, int numSamples);

	// 计算两条曲线之间的夹角（以弧度为单位）
	static double ComputeAngleBetweenCurves(Handle(Geom_BSplineCurve)& curve1,
		Handle(Geom_BSplineCurve)& curve2,
		int numSamples = 10);

	// 根据参考曲线对 B-Spline 曲线进行排序
	static void SortBSplineCurves(std::vector<Handle(Geom_BSplineCurve)>& theCurves,
		Handle(Geom_BSplineCurve) referCurve);

	// 根据需要反转曲线的方向，以确保一致性
	static void ReverseIfNeeded(std::vector<Handle(Geom_BSplineCurve)>& curves);

	// 计算点到平面的距离
	static double ComputeDistancePointToPlane(const gp_Pnt& p, const gp_Pln& plane);

	// 计算直线和平面之间的夹角
	static double ComputeAngleBetweenLineAndPlane(const gp_Lin& line, const gp_Pln& plane);

	// 计算两条直线之间的夹角
	static double ComputeAngleBetweenLines(const gp_Lin& line1, const gp_Lin& line2);

	// 计算两个平面之间的夹角
	static double ComputeAngleBetweenPlanes(const gp_Pln& plane1, const gp_Pln& plane2);

	// 计算点到直线的距离
	static double ComputeDistancePointToLine(const gp_Pnt& point, const gp_Lin& line);

	// 计算两个平面曲线之间的夹角
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
