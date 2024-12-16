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

// OpenCASCADE���ͷ�ļ�
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
	void SetCurve(Handle(Geom_BSplineCurve) theCurve)
	{
		curve = theCurve;
		IsPlanarCurve(theCurve);
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

class SurfaceModelingTool
{
public:
	//Լ�������������������µķ���Ҫ��
	// ���������Ӻͷ�������Ҫ�� G0������Coons���湹��
	// 
	//                -------------curve1--->-------
	//                |                             |
	//              curve4                        curve 2
	//                v                             v
	//                |                             |
	//                |-----------curve3--->--------
	static void Coons_G0(Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4, Handle(Geom_BSplineSurface)& mySurface_coons);

	//              Լ�����������ߺ������ߵĿ�絼ʸ�������µķ���Ҫ��
	//                           G1������Coons���湹��
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


	// ���ڵ����������txt�ļ�
	Standard_Boolean KnotsToTxt(const std::vector<double>& knots) const 
	{
		// ���ļ�
		std::ofstream outFile(knotsOutputPath, std::ios::app);

		if (!outFile.is_open())
		{
			return Standard_False;
		}

		// ���������С�������λ
		outFile << std::fixed << std::setprecision(3); 

		// д�����ݵ��ļ�
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

	// ��������
	static bool ExportBSplineCurves(
		const std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_Final,
		const std::string& Filename);

	static void ApproximateBoundaryCurves(std::vector<Handle(Geom_BSplineCurve)>& curves, int samplingNum = 50);

	/**
	 * @brief ��ȡ�������ڲ����ߣ�������߽����ߵĽǶȺ;�����з��ࡣ
	 *
	 * @param aBoundarycurveArray �߽��������顣
	 * @param anInternalBSplineCurves �ڲ� B-Spline �������顣
	 * @param uInternalCurve ������ u �����ڲ��������顣
	 * @param vInternalCurve ������ v �����ڲ��������顣
	 * @param uAngleSum �ۼƵ� u ����ǶȺ͡�
	 * @param vAngleSum �ۼƵ� v ����ǶȺ͡�
	 * @param AngleTolerance �Ƕ��ݲ�ֵ��Ĭ��ֵΪ 5��
	 * @return true ����ɹ���ȡ�������ڲ����ߡ�
	 * @return false ����޷�����������
	 */
	static bool GetInternalCurves(
		std::vector<Handle(Geom_BSplineCurve)>& aBoundarycurveArray,
		std::vector<Handle(Geom_BSplineCurve)>& anInternalBSplineCurves,
		std::vector<Handle(Geom_BSplineCurve)>& uInternalCurve,
		std::vector<Handle(Geom_BSplineCurve)>& vInternalCurve,
		double& uAngleSum,
		double& vAngleSum,
		double AngleTolerance = 5);

	// ���� ReferSurfaceType ���ɲο� B-Spline ���棬���ݱ߽����ߺ��ڲ����߹���
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




#endif