#pragma once
#include "Geom_BSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Edge.hxx"
#include "vector"


#include "gp_Pln.hxx"
#include "TopoDS_Face.hxx"
#include <Interpolate.h>

// OpenCASCADE���ͷ�ļ�
#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pln.hxx>
#include <gp_Lin.hxx>
#include <gp_XYZ.hxx>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
// ö�����ͣ���ʾ���ߵĲ�ͬ����
enum class CurveType
{
	LINEAR, // ��������
	PLANAR, // ƽ������
	NOTPLANAR, // ��ƽ������
	POINT // ��
};
// PlanarCurve �࣬�����жϺʹ���ƽ������
class PlanarCurve
{
public:
	// Ĭ�Ϲ��캯��
	PlanarCurve();
	// �������Ĺ��캯��������һ�� B-Spline ���߾�����ݲ�ֵ
	PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, Standard_Real theTolerance = 10);
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
		IsPlanarCurve(theCurve, 10);
	}
private:
	CurveType curveType; // ��ǰ���ߵ�����
	gp_Lin line; // �������ߵı�ʾ����������Ϊ����ʱ��Ч��
	gp_Pln plane; // ƽ��ı�ʾ����������Ϊƽ������ʱ��Ч��
	gp_Pnt pnt; // ��ı�ʾ����������Ϊ��ʱ��Ч��
	Handle(Geom_BSplineCurve) curve; // B-Spline ���ߵľ��
	// �жϸ����� B-Spline �����Ƿ�Ϊƽ������
	bool IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, Standard_Real
		theTolerance = 0.1);
	// �жϸ����� B-Spline �����Ƿ�Ϊ��������
	bool IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve,
		Standard_Real theTolerance = 1e-6);
	// �жϸ����� B-Spline �����Ƿ�Ϊһ����
	bool IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve,
		Standard_Real theTolerance = 1e-6);
};
class MathTool
{
public:
	// �������� B-Spline ����֮�����С����
	static Standard_Real ComputeCurveCurveDistance(const
		Handle(Geom_BSplineCurve)& theCurve, const Handle(Geom_BSplineCurve)&
		theBoundaryCurve);
	// ���������ϲ������ƽ������
	static gp_Pnt ComputeAverageSamplePoint(const Handle(Geom_BSplineCurve)&
		theCurve, Standard_Integer theNumSamples);
	// ������������֮��ļнǣ��Ի���Ϊ��λ��
	static Standard_Real ComputeAngleWithAxis(const gp_Vec& theVec, const gp_Vec&
		theAxis);
	// ���һ�� B-Spline �����Ƿ�����Խ�
	static void CheckSelfIntersect(std::vector<Handle(Geom_BSplineCurve)>
		theBSplineCurvesArray);
	// �������ߵ�ƽ������������
	static gp_Dir ComputeAverageTangent(const Handle(Geom_BSplineCurve)&
		theCurve, Standard_Integer theNumSamples);
	// ������������֮��ļнǣ��Ի���Ϊ��λ��
	static Standard_Real ComputeAngleBetweenCurves(Handle(Geom_BSplineCurve)&
		theCurve1,
		Handle(Geom_BSplineCurve)& theCurve2,
		Standard_Integer theNumSamples = 10);
		// ���ݲο����߶� B-Spline ���߽�������
		static void SortBSplineCurves(std::vector<Handle(Geom_BSplineCurve)>&
			theCurves,
			Handle(Geom_BSplineCurve) theReferCurve);
	// ������Ҫ��ת���ߵķ�����ȷ��һ����
	static void ReverseIfNeeded(std::vector<Handle(Geom_BSplineCurve)>&
		theCurves);
	// ����㵽ƽ��ľ���
	static Standard_Real ComputeDistancePointToPlane(const gp_Pnt& theP, const
		gp_Pln& thePlane);
	// ����ֱ�ߺ�ƽ��֮��ļн�
	static Standard_Real ComputeAngleBetweenLineAndPlane(const gp_Lin& theLine,
		const gp_Pln& thePlane);
	// ��������ֱ��֮��ļн�
	static Standard_Real ComputeAngleBetweenLines(const gp_Lin& theLine1, const
		gp_Lin& theLine2);
	// ��������ƽ��֮��ļн�
	static Standard_Real ComputeAngleBetweenPlanes(const gp_Pln& thePlane1, const
		gp_Pln& thePlane2);
	// ����㵽ֱ�ߵľ���
	static Standard_Real ComputeDistancePointToLine(const gp_Pnt& thePoint, const
		gp_Lin& theLine);
	// ��������ƽ������֮��ļн�
	static Standard_Real ComputeAngleBetweenPlanarCurves(const PlanarCurve&
		theCurve1, const PlanarCurve& theCurve2);
	static std::vector<gp_Pnt> GetSamplePointsOnCurve(const Handle(Geom_Curve)&
		theCurve, Standard_Integer theNumPoints = 50);
	static void SortPoints(std::vector<gp_Pnt>& thePoints, const gp_Pnt&
		theReferPoint);
	static void TrimInternalCurves(
		std::vector<Handle(Geom_BSplineCurve)>& theInternalBSplineCurves,
		const std::vector<Handle(Geom_BSplineCurve)>& theBoundaryCurveArray,
		Standard_Real theToleranceDistance = 10);
};

enum ReferSurfaceType
{
	GORDEN_ONE_DIRECTION_GORDEN,
	GORDEN_TWO_DIRECTION_GORDEN
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

	

private:
	std::string knotsOutputPath;
	
};

class CurveOperate
{
public:
	static Standard_Boolean CompatibleWithInterPoints(const std::vector<Handle(Geom_BSplineCurve)>& theInterCurves, std::vector<Handle(Geom_BSplineCurve)>& theCompatibleCurves, Standard_Real theTolerance = 0.01);
	static std::tuple<std::vector<gp_Pnt>, std::vector<Standard_Real>> CalCurvesInterPointsParamsToCurve(const std::vector<Handle(Geom_BSplineCurve)>& theCurves, const Handle(Geom_BSplineCurve)& theCurve, Standard_Real theTolerance = 0.1);
};