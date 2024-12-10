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

// ö�����ͣ���ʾ��������
enum class CurveType
{
    LINEAR,    // ��������
    PLANAR,    // ƽ������
    NOTPLANAR, // ��ƽ������
    POINT
};

// PlanarCurve �࣬�����жϺʹ���ƽ������
class PlanarCurve
{
public:
    // ���캯��
    PlanarCurve();

    PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);

    static double ComputeAngleBetweenPlanarCurves(const PlanarCurve& curve1, const PlanarCurve& curve2);
    static double ComputeDistancePointToPlane(const gp_Pnt& p, const gp_Pln& plane);
    static double ComputeAngleBetweenLineAndPlane(const gp_Lin& line, const gp_Pln& plane);
    static double ComputeAngleBetweenLines(const gp_Lin& line1, const gp_Lin& line2);
    static double ComputeAngleBetweenPlanes(const gp_Pln& plane1, const gp_Pln& plane2);
    // ����㵽ֱ�ߵľ���
    static double ComputeDistancePointToLine(const gp_Pnt& point, const gp_Lin& line);
    // ��ȡ��������
    CurveType GetCurveType() const { return curveType; }

    // ��ȡ�������ߣ�������������Ϊ LINEAR ʱ��Ч��
    gp_Lin GetLine() const { return line; }

    // ��ȡƽ�棨������������Ϊ PLANAR ʱ��Ч��
    gp_Pln GetPlane() const { return plane; }

    gp_Pnt GetPoint() const { return pnt; }
    Handle(Geom_BSplineCurve) GetCurve() const { return curve; }

    void SetCurve(Handle(Geom_BSplineCurve) theCurve) { curve = theCurve; }

private:
    CurveType curveType; // ��������
    gp_Lin line;         // �������߱�ʾ
    gp_Pln plane;        // ƽ���ʾ
    gp_Pnt pnt;
    Handle(Geom_BSplineCurve) curve;

    bool IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance = 0.1);
    bool IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);
    bool IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance = 1e-6);
};

#endif // PLANARCURVE_H
