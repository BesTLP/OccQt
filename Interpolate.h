#pragma once

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <gp_Pln.hxx>
#include <gp_Cylinder.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <STEPControl_Writer.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <Math_matrix.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <TopoDS_Face.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Eigen/Dense>

namespace InterPolateTool
{
	/// <summary>
	/// 插值曲线
	/// </summary>
	/// <param name="DataPnts">数据点</param>
	/// <param name="PntParams">数据点相对应的参数化</param>
	/// <param name="knotSeqs">节点向量</param>
	/// <param name="degree">次数</param>
	Handle(Geom_BSplineCurve) Interpolate(const std::vector<gp_Pnt>& DataPnts, const std::vector<double>& PntParams, const std::vector<double>& knotSeqs, int degree);

	Handle(Geom_BSplineSurface) Interpolate(const std::vector<gp_Pnt>& Pnts, const std::vector<gp_Pnt2d>& PntParams,
		std::vector<double>& UKnots, std::vector<double>& VKnots, std::vector<int>& UMuti, std::vector<int>& VMuti, int Udegree, int Vdegree);

	Handle(Geom_BSplineSurface) Loft(const std::vector<Handle(Geom_BSplineCurve)>& isoCurves, int perpendDegree = 3);

	std::vector<double> CalKnots(const std::vector<double>& params, int degree);

	std::vector<double> ChordLenParam(const std::vector<std::vector<gp_Pnt>>& PointMatrix);

	std::vector<std::vector<gp_Pnt>> ExtractControlPoints(const std::vector<Handle(Geom_BSplineCurve)>& bsplineCurves);

	double OneBasicFun(double u, int i, int p, const std::vector<double>& Knots);

	void sequenceToKnots(const std::vector<double>& sequence, std::vector<double>& knots, std::vector<int>& multiplicities);
	std::vector<double> knotsToSequence(const std::vector<double>& knots, const std::vector<int>& multiplicities);

	Eigen::MatrixXd composeTensowMat(const Eigen::MatrixXd& matA, const Eigen::MatrixXd& matB);
	TColgp_Array2OfPnt vectorToOCCMatrix(const std::vector<gp_Pnt>& pnts, int m);
	Handle(Geom_BSplineSurface) LoftV(const std::vector<Handle(Geom_BSplineCurve)>& isoCurves, int perpendDegree);
	bool isEqual(double x, double y, double epsilon = 0.0000000000001);
	bool isGreaterThan(double x, double y, double epsilon = 0.0000000000001);
	bool isLessThan(double x, double y, double epsilon = 0.0000000000001);
	bool isGreaterThanOrEqual(double x, double y, double epsilon = 0.0000000000001);
	bool isLessThanOrEqual(double x, double y, double epsilon = 0.0000000000001);
}