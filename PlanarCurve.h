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

#endif // PLANARCURVE_H
