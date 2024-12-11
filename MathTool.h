#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>

// OpenCASCADE���ͷ�ļ�
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

#endif // MATHTOOL_H
