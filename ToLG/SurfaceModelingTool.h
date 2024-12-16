#ifndef SURFACEMODELINGTOOL_H
#define SURFACEMODELINGTOOL_H
#pragma once
#include "Geom_BSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Edge.hxx"
#include "vector"


#include "gp_Pln.hxx"
#include "TopoDS_Face.hxx"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>

// OpenCASCADE相关头文件
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <TopoDS_Face.hxx>
#include <Geom_BSplineCurve.hxx>
#include <vector>
#include <Eigen/Dense>
#include <gp_Pnt.hxx>
#include <gp_Lin.hxx>
#include <gp_XYZ.hxx>
#include <gp_Pln.hxx>
#include <iostream>
enum ReferSurfaceType
{
	GORDEN_ONE_DIRECTION_COONS,
	GORDEN_ONE_DIRECTION_GORDEN,
	GORDEN_TWO_DIRECTION
};

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
	void SetCurve(Handle(Geom_BSplineCurve) theCurve)
	{
		curve = theCurve;
		IsPlanarCurve(theCurve);
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

class SurfaceModelingTool
{
public:
	//约定输入四条边满足如下的方向要求
	// 四条边连接和方向满足要求 G0连续的Coons曲面构造
	// 
	//                -------------curve1--->-------
	//                |                             |
	//              curve4                        curve 2
	//                v                             v
	//                |                             |
	//                |-----------curve3--->--------
	static void Coons_G0(Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4, Handle(Geom_BSplineSurface)& mySurface_coons);

	//              约定输入四条边和四条边的跨界导矢满足如下的方向要求
	//                           G1连续的Coons曲面构造
	// 
	//                ----------------------c1--->------------
	//                |                     |                |
	//               c4-->c4_d            c1_d              c2 -->-c2_d-->
	//                v                     v                v
	//                |                     |                |
	//                |--------------------c3--->-------------
	//                                      |
	//                                      v
	//                                    c3_d  
	//                                      | 
	static void Coons_G1(Handle(Geom_BSplineCurve)& c1, Handle(Geom_BSplineCurve)& c2, Handle(Geom_BSplineCurve)& c3, Handle(Geom_BSplineCurve)& c4, Handle(Geom_BSplineCurve)& c1_derivative, Handle(Geom_BSplineCurve)& curve2_derivative, Handle(Geom_BSplineCurve)& curve3_derivative, Handle(Geom_BSplineCurve)& curve4_derivative, Handle(Geom_BSplineSurface)& mySurface_coons);

	//make the curve same degree and knots
	static Standard_Integer SetSameDistribution(Handle(Geom_BSplineCurve)& C1, Handle(Geom_BSplineCurve)& C2);

	//make the curve arranged for the G0 construction
	// Input: curveArray, the curve array
	// Output: bslpineCurve1-bslpineCurve4 the arranged four boundary curves
	// IsModify: If the curve are not connected end to end, whether change the end point of the curve
	// Tol: the tolerance to check whether the four curves are connected or not

	static int Arrange_Coons_G0(std::vector<Handle(Geom_BSplineCurve)>& curveArray, Handle(Geom_BSplineCurve)& bslpineCurve1, Handle(Geom_BSplineCurve)& bslpineCurve2, Handle(Geom_BSplineCurve)& bslpineCurve3, Handle(Geom_BSplineCurve)& bslpineCurve4, double Tol, int IsModify);

	static void ClassifyAndSortISOcurves(const std::vector<Handle(Geom_BSplineCurve)>& anISOcurvesArray, std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray, std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray);

	static void CreateLoftingSurface(const std::vector<Handle(Geom_BSplineCurve)>& curvesArray,
		const std::vector<gp_Vec>& normals,
		std::vector<TopoDS_Shape>& loftingSurfaces);

	static void LoftSurfaceIntersectWithCurve(
		const std::vector<TopoDS_Shape>& LoftingSur,
		const std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_Initial,
		const std::vector<Handle(Geom_BSplineCurve)>& anInternalBSplineCurves,
		std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_New,
		Standard_Integer isoCount,
		std::vector<std::vector<gp_Pnt>>& InterpolatePoints,
		std::vector<TopoDS_Edge>& TangentArray1,
		std::vector<TopoDS_Edge>& TangentArray2,
		Handle(Geom_BSplineSurface) CoonsSurface);

	static void CreateFinalISOCurves(
		const std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_New,
		const std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_New,
		std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Final,
		std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Final,
		std::vector<std::vector<gp_Pnt>>& uInterpolatePoints,
		std::vector<std::vector<gp_Pnt>>& vInterpolatePoints,
		std::vector<std::vector<double> >& uKnots,
		std::vector<std::vector<double> >& vKnots,
		std::vector<gp_Pnt>& boundaryPoints,
		std::vector<gp_Pnt>& interPoints,
		Standard_Integer isoCount,
		std::vector<TopoDS_Edge>& TangentArray,
		std::vector<Handle(Geom_BSplineSurface)>& surfaceArr);



	static void UpdateFinalCurves(const std::vector<Handle(Geom_BSplineCurve)>& aBoundarycurveArray,
		std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Final,
		std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Final);

	static void LoadBSplineCurves(const std::string& filePath, std::vector<Handle(Geom_BSplineCurve)>& curveArray);
	static void LoadBSplineSurfaces(const std::string& filePath, std::vector<Handle(Geom_BSplineSurface)>& surfaceArray);
	static void GetISOCurveWithNormal(
		const Handle(Geom_BSplineSurface)& surfacecoons,
		std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Initial,
		std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Initial,
		std::vector<gp_Vec>& normalsOfUISOLines,
		std::vector<gp_Vec>& normalsOfVISOLines,
		int numIsoCurves = 10);


	void setKnotsOutputPath(std::string knotsOuputPath)
	{
		this->knotsOutputPath = knotsOuputPath;
	}
	std::string getKnotsOuputPath()
	{
		return this->knotsOutputPath;
	}


	// 将节点数据输出到txt文件
	Standard_Boolean KnotsToTxt(const std::vector<double>& knots) const 
	{
		// 打开文件
		std::ofstream outFile(knotsOutputPath, std::ios::app);

		if (!outFile.is_open())
		{
			return Standard_False;
		}

		// 设置输出到小数点后三位
		outFile << std::fixed << std::setprecision(3); 

		// 写入数据到文件
		outFile << "[";
		for (size_t i = 0; i < knots.size(); ++i) 
		{
			outFile << knots[i];
			if (i < knots.size() - 1)
			{
				outFile << ",";
			}
		}
		outFile << "]\n";

		outFile.close();
		return Standard_True;
	}
	Standard_Boolean ContextToTxt(const std::string context) const
	{
		std::ofstream outFile(knotsOutputPath, std::ios::app);

		if (!outFile.is_open())
		{
			return Standard_False;
		}

		outFile << context << "\n";

		outFile.close();

		return Standard_True;
	}

	// 函数声明
	static bool ExportBSplineCurves(
		const std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_Final,
		const std::string& Filename);

	static void ApproximateBoundaryCurves(std::vector<Handle(Geom_BSplineCurve)>& curves, int samplingNum = 50);

	/**
	 * @brief 获取并分类内部曲线，根据与边界曲线的角度和距离进行分类。
	 *
	 * @param aBoundarycurveArray 边界曲线数组。
	 * @param anInternalBSplineCurves 内部 B-Spline 曲线数组。
	 * @param uInternalCurve 分类后的 u 方向内部曲线数组。
	 * @param vInternalCurve 分类后的 v 方向内部曲线数组。
	 * @param uAngleSum 累计的 u 方向角度和。
	 * @param vAngleSum 累计的 v 方向角度和。
	 * @param AngleTolerance 角度容差值，默认值为 5。
	 * @return true 如果成功获取并分类内部曲线。
	 * @return false 如果无法满足条件。
	 */
	static bool GetInternalCurves(
		std::vector<Handle(Geom_BSplineCurve)>& aBoundarycurveArray,
		std::vector<Handle(Geom_BSplineCurve)>& anInternalBSplineCurves,
		std::vector<Handle(Geom_BSplineCurve)>& uInternalCurve,
		std::vector<Handle(Geom_BSplineCurve)>& vInternalCurve,
		double& uAngleSum,
		double& vAngleSum,
		double AngleTolerance = 5);

	// 根据 ReferSurfaceType 生成参考 B-Spline 曲面，根据边界曲线和内部曲线构造
	static Handle(Geom_BSplineSurface) GenerateReferSurface(
		std::vector<Handle(Geom_BSplineCurve)> aBoundarycurveArray,
		std::vector<Handle(Geom_BSplineCurve)>& uInternalCurve,
		std::vector<Handle(Geom_BSplineCurve)>& vInternalCurve,
		double uAngleSum,
		double vAngleSum,
		int isoCount,
		ReferSurfaceType referSurfaceType);


	bool CompatibleWithInterPoints(const std::vector<Handle(Geom_BSplineCurve)>& interCurves, std::vector<Handle(Geom_BSplineCurve)>& compatibleCurves, Standard_Real toler = 1.e-3);

	std::vector<Standard_Real> CalSameKnotFromCurves(std::vector< Handle(Geom_BSplineCurve) >& curves, Standard_Real toler = 0.1);



private:
	std::string knotsOutputPath;
	
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




#endif