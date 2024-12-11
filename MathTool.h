#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>

// OpenCASCADE相关头文件
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
#include <PlanarCurve.h>

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

#endif // MATHTOOL_H
