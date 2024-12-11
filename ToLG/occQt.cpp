Handle(Geom_BSplineCurve) bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4;
SurfaceModelingTool::Arrange_Coons_G0(tempArray, bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4, 10, Standard_True);

std::vector<Handle(Geom_BSplineCurve)> aBoundarycurveArray = { bslpineCurve1 , bslpineCurve2, bslpineCurve3, bslpineCurve4 };

std::vector<Handle(Geom_BSplineCurve)> anInternalBSplineCurves;
std::string internalPath = filename + "internal.brep";
SurfaceModelingTool::LoadBSplineCurves(internalPath, anInternalBSplineCurves);

int isoCount = 20;
std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_Initial, vISOcurvesArray_Initial;
std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_Final, vISOcurvesArray_Final;
std::vector<gp_Vec> normalsOfUISOLines, normalsOfVISOLines;
std::vector<std::vector<double>> uKnots;
std::vector<std::vector<double>> vKnots;

std::vector<Handle(Geom_BSplineCurve)> uInternalCurve, vInternalCurve;
double uAngleSum = 0, vAngleSum = 0;

if (SurfaceModelingTool::GetInternalCurves(aBoundarycurveArray, anInternalBSplineCurves, uInternalCurve, vInternalCurve, uAngleSum, vAngleSum))
{
    Handle(Geom_BSplineSurface) referSurface;
    referSurface = SurfaceModelingTool::GenerateReferSurface(aBoundarycurveArray,
        uInternalCurve, vInternalCurve,
        uAngleSum, vAngleSum,
        isoCount, ReferSurfaceType::GORDEN_ONE_DIRECTION_GORDEN);
}