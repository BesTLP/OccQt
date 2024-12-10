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

    // ��������ƽ������֮��ļн�
    static double ComputeAngleBetweenPlanarCurves(const PlanarCurve& curve1, const PlanarCurve& curve2);

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
    void SetCurve(Handle(Geom_BSplineCurve) theCurve) { curve = theCurve; }

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

#endif // PLANARCURVE_H
