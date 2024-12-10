// PlanarCurve.h
#ifndef PLANARCURVE_H
#define PLANARCURVE_H

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

// 枚举类型，表示曲线类型
enum class CurveType
{
    LINEAR,    // 线性曲线
    PLANAR,    // 平面曲线
    NOTPLANAR, // 非平面曲线
    POINT
};

// PlanarCurve 类，用于判断和处理平面曲线
class PlanarCurve
{
public:
    // 构造函数
    PlanarCurve();

    PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);

    static double ComputeAngleBetweenPlanarCurves(const PlanarCurve& curve1, const PlanarCurve& curve2);
    static double ComputeDistancePointToPlane(const gp_Pnt& p, const gp_Pln& plane);
    static double ComputeAngleBetweenLineAndPlane(const gp_Lin& line, const gp_Pln& plane);
    static double ComputeAngleBetweenLines(const gp_Lin& line1, const gp_Lin& line2);
    static double ComputeAngleBetweenPlanes(const gp_Pln& plane1, const gp_Pln& plane2);
    // 计算点到直线的距离
    static double ComputeDistancePointToLine(const gp_Pnt& point, const gp_Lin& line);
    // 获取曲线类型
    CurveType GetCurveType() const { return curveType; }

    // 获取线性曲线（仅当曲线类型为 LINEAR 时有效）
    gp_Lin GetLine() const { return line; }

    // 获取平面（仅当曲线类型为 PLANAR 时有效）
    gp_Pln GetPlane() const { return plane; }

    gp_Pnt GetPoint() const { return pnt; }
    Handle(Geom_BSplineCurve) GetCurve() const { return curve; }

    void SetCurve(Handle(Geom_BSplineCurve) theCurve) { curve = theCurve; }

private:
    CurveType curveType; // 曲线类型
    gp_Lin line;         // 线性曲线表示
    gp_Pln plane;        // 平面表示
    gp_Pnt pnt;
    Handle(Geom_BSplineCurve) curve;

    bool IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);
    bool IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);
    bool IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);
};

#endif // PLANARCURVE_H
