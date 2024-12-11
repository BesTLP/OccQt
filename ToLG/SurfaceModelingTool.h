#pragma once
#include "gp_Pln.hxx"
#include "TopoDS_Face.hxx"
#include <Interpolate.h>
#include <PlanarCurve.h>
#include <GordenSurface.h>
#include <MathTool.h>
enum ReferSurfaceType
{
	GORDEN_ONE_DIRECTION_COONS,
	GORDEN_ONE_DIRECTION_GORDEN,
	GORDEN_TWO_DIRECTION
};

class SurfaceModelingTool
{
public:

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

private:
	
};
