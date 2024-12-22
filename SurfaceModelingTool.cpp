#include "SurfaceModelingTool.h"
#include "GeomFill_BSplineCurves.hxx"
#include "BSplCLib.hxx"
#include "Geom_BSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <algorithm>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_IntCS.hxx>
#include <TopExp_Explorer.hxx>
#include <GeomAPI_Interpolate.hxx>
#include "GeomAPI_ExtremaCurveCurve.hxx"
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAdaptor_Curve.hxx>
#include "TColgp_HArray1OfPnt.hxx"
#include "Geom_TrimmedCurve.hxx"
#include <BRepTools.hxx>
#include "BRep_Builder.hxx"
#include <GeomConvert.hxx>
#include "Geom_Line.hxx"
#include <STEPControl_Reader.hxx>
//3rd part
#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include <cmath>
#include <memory>
#include <algorithm>
#include <functional>
#include <numeric> 
#include <IGESControl_Reader.hxx>
#include <KnotUpdate.h>
#include <RealCompare.h>
#include <GeomAPI_ProjectPointOnSurf.hxx>

#include <filesystem>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <GordenSurface.h>
#include <occQt.h>
//=======================================================================
//function : SetSameDistribution
//purpose  : Internal Use Only
//=======================================================================

Standard_Integer SurfaceModelingTool::SetSameDistribution(Handle(
	Geom_BSplineCurve)& C1,
	Handle(Geom_BSplineCurve)& C2) 
{

	Standard_Integer C1_Degree = C1->Degree();
	Standard_Integer C2_Degree = C2->Degree();

	// ���Degree����ͬ�Ļ����Եͽ׽�������
	if (C1_Degree < C2_Degree)
	{
		C1->IncreaseDegree(C2_Degree);
	}
	else
	{
		C2->IncreaseDegree(C1_Degree);
	}
	Standard_Integer nbp1 = C1->NbPoles();
	Standard_Integer nbk1 = C1->NbKnots();
	TColgp_Array1OfPnt      P1(1, nbp1);
	TColStd_Array1OfReal    W1(1, nbp1);
	W1.Init(1.);
	TColStd_Array1OfReal    K1(1, nbk1);
	TColStd_Array1OfInteger M1(1, nbk1);

	C1->Poles(P1);
	if (C1->IsRational())
		C1->Weights(W1);
	C1->Knots(K1);
	C1->Multiplicities(M1);

	Standard_Integer nbp2 = C2->NbPoles();
	Standard_Integer nbk2 = C2->NbKnots();
	TColgp_Array1OfPnt      P2(1, nbp2);
	TColStd_Array1OfReal    W2(1, nbp2);
	W2.Init(1.);
	TColStd_Array1OfReal    K2(1, nbk2);
	TColStd_Array1OfInteger M2(1, nbk2);

	C2->Poles(P2);
	if (C2->IsRational())
		C2->Weights(W2);
	C2->Knots(K2);
	C2->Multiplicities(M2);

	Standard_Real K11 = K1(1);
	Standard_Real K12 = K1(nbk1);
	Standard_Real K21 = K2(1);
	Standard_Real K22 = K2(nbk2);

	if ((K12 - K11) > (K22 - K21)) {
		BSplCLib::Reparametrize(K11, K12, K2);
		C2->SetKnots(K2);
	}
	else if ((K12 - K11) < (K22 - K21)) {
		BSplCLib::Reparametrize(K21, K22, K1);
		C1->SetKnots(K1);
	}
	else if (Abs(K12 - K11) > 1.e-15) {
		BSplCLib::Reparametrize(K11, K12, K2);
		C2->SetKnots(K2);
	}

	Standard_Integer NP, NK;
	if (BSplCLib::PrepareInsertKnots(C1->Degree(), Standard_False,
		K1, M1, K2, &M2, NP, NK, 1.e-15,
		Standard_False)) {
		TColgp_Array1OfPnt      NewP(1, NP);
		TColStd_Array1OfReal    NewW(1, NP);
		TColStd_Array1OfReal    NewK(1, NK);
		TColStd_Array1OfInteger NewM(1, NK);
		BSplCLib::InsertKnots(C1->Degree(), Standard_False,
			P1, &W1, K1, M1, K2, &M2,
			NewP, &NewW, NewK, NewM, 1.e-15,
			Standard_False);
		if (C1->IsRational()) {
			C1 = new Geom_BSplineCurve(NewP, NewW, NewK, NewM, C1->Degree());
		}
		else {
			C1 = new Geom_BSplineCurve(NewP, NewK, NewM, C1->Degree());
		}
		BSplCLib::InsertKnots(C2->Degree(), Standard_False,
			P2, &W2, K2, M2, K1, &M1,
			NewP, &NewW, NewK, NewM, 1.e-15,
			Standard_False);
		if (C2->IsRational()) {
			C2 = new Geom_BSplineCurve(NewP, NewW, NewK, NewM, C2->Degree());
		}
		else {
			C2 = new Geom_BSplineCurve(NewP, NewK, NewM, C2->Degree());
		}
	}
	else {
		throw Standard_ConstructionError(" ");
	}

	return C1->NbPoles();
}
void SurfaceModelingTool::Coons_G0(Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4, Handle(Geom_BSplineSurface)& mySurface_coons)
{
	// First make the opposite boundary curves compatible,
	// Second generate the ruled surface in the u/v direction and of the four corner points
	// Third generate the sum surface 
	//1: Make Curve Compatible
	//curve 1 and curve 3 are opposite curves, make them compatible
	SetSameDistribution(curve1, curve3);
	//curve 2 and curve 4 are opposite curves, make them compatible
	SetSameDistribution(curve2, curve4);

	//2: Generate the coons Bspline surface

	//Get the control point number in the u/v direction
	Standard_Integer NbUPoles = curve1->NbPoles();
	Standard_Integer NbVPoles = curve2->NbPoles();

	//control points of the coons surface
	TColgp_Array2OfPnt Poles_result(1, NbUPoles, 1, NbVPoles);

	//control points of ruled surface in the v direction
	TColgp_Array2OfPnt Poles_vruled(1, NbUPoles, 1, 2);

	//control points of ruled surface in the u direction
	TColgp_Array2OfPnt Poles_uruled(1, 2, 1, NbVPoles);

	//control points of ruled surface of the four corner points
	TColgp_Array2OfPnt Poles_uruled_vruled(1, 2, 1, 2);

	// Get the u knot vector
	Standard_Integer NbUKnot = curve1->NbKnots();
	TColStd_Array1OfReal    UKnots(1, NbUKnot);
	TColStd_Array1OfInteger UMults(1, NbUKnot);
	curve1->Knots(UKnots);
	curve1->Multiplicities(UMults);

	// Get the v knot vector
	Standard_Integer NbVKnot = curve2->NbKnots();
	TColStd_Array1OfReal    VKnots(1, NbVKnot);
	TColStd_Array1OfInteger VMults(1, NbVKnot);
	curve2->Knots(VKnots);
	curve2->Multiplicities(VMults);

	// Set v knots of the v ruled surface
	TColStd_Array1OfReal    VKnots_RuledInVdirection(1, 2);
	TColStd_Array1OfInteger VMults_RuledInVdirection(1, 2);
	VKnots_RuledInVdirection(1) = curve2->FirstParameter();
	VKnots_RuledInVdirection(2) = curve2->LastParameter();

	VMults_RuledInVdirection(1) = 2;
	VMults_RuledInVdirection(2) = 2;

	// Set u knots of the u ruled surface
	TColStd_Array1OfReal    UKnots_RuledInUdirection(1, 2);
	TColStd_Array1OfInteger UMults_RuledInUdirection(1, 2);
	UKnots_RuledInUdirection(1) = curve1->FirstParameter();
	UKnots_RuledInUdirection(2) = curve1->LastParameter();

	UMults_RuledInUdirection(1) = 2;
	UMults_RuledInUdirection(2) = 2;

	// 
	//2.1 generate the ruled surface in the v direction
	// 
	// 
	// Set control points of the v ruled surface
	for (Standard_Integer i = 1; i <= NbUPoles; i++)
	{
		Poles_vruled(i, 1) = curve1->Pole(i);
		Poles_vruled(i, 2) = curve3->Pole(i);
	}

	// Generate the bspline surface
	Handle(Geom_BSplineSurface) mySurface_VRuled = new Geom_BSplineSurface(Poles_vruled,
		UKnots, VKnots_RuledInVdirection,
		UMults, VMults_RuledInVdirection,
		curve1->Degree(), 1);

	// Make the surface compatible with the result coons surface
	mySurface_VRuled->IncreaseDegree(curve1->Degree(), curve2->Degree());

	for (Standard_Integer i = 1; i <= NbVKnot; i++)
		mySurface_VRuled->InsertVKnot(VKnots(i), VMults(i), Precision::Confusion());

	// 
	//2.2 generate the ruled surface in the u direction
	// 

	// Set control points of the u ruled surface
	for (Standard_Integer i = 1; i <= NbVPoles; i++)
	{
		Poles_uruled(1, i) = curve4->Pole(i);
		Poles_uruled(2, i) = curve2->Pole(i);
	}

	Handle(Geom_BSplineSurface) mySurface_URuled = new Geom_BSplineSurface(Poles_uruled,
		UKnots_RuledInUdirection, VKnots,
		UMults_RuledInUdirection, VMults,
		1, curve2->Degree());

	//increase degree
	mySurface_URuled->IncreaseDegree(curve1->Degree(), curve2->Degree());

	for (Standard_Integer i = 1; i <= NbUKnot; i++)
		mySurface_URuled->InsertUKnot(UKnots(i), UMults(i), Precision::Confusion());

	// 
	//2.3 generate the ruled surface of the four corner points
	// 
	// 
	Poles_uruled_vruled(1, 1) = curve1->Pole(1);
	Poles_uruled_vruled(1, 2) = curve3->Pole(1);
	Poles_uruled_vruled(2, 1) = curve1->Pole(NbUPoles);
	Poles_uruled_vruled(2, 2) = curve3->Pole(NbUPoles);

	Handle(Geom_BSplineSurface) mySurface_RuledSurfaceof4CornerPoints = new Geom_BSplineSurface(Poles_uruled_vruled,
		UKnots_RuledInUdirection, VKnots_RuledInVdirection,
		UMults_RuledInUdirection, VMults_RuledInVdirection,
		1, 1);

	mySurface_RuledSurfaceof4CornerPoints->IncreaseDegree(curve1->Degree(), curve2->Degree());

	for (Standard_Integer i = 1; i <= NbVKnot; i++)
		mySurface_RuledSurfaceof4CornerPoints->InsertVKnot(VKnots(i), VMults(i), Precision::Confusion());

	for (Standard_Integer i = 1; i <= NbUKnot; i++)
		mySurface_RuledSurfaceof4CornerPoints->InsertUKnot(UKnots(i), UMults(i), Precision::Confusion());

	//2.4 Generate the sum surface mySurface_VRuled + mySurface_URuled - mySurface_RuledSurfaceof4CornerPoints 
	//increase degree

	for (Standard_Integer i = 1; i <= NbUPoles; i++)
	{
		for (Standard_Integer j = 1; j <= NbVPoles; j++)
		{
			Poles_result(i, j).SetXYZ(mySurface_VRuled->Pole(i, j).XYZ() + mySurface_URuled->Pole(i, j).XYZ() - mySurface_RuledSurfaceof4CornerPoints->Pole(i, j).XYZ());
		}
	}

	mySurface_coons = new Geom_BSplineSurface(Poles_result,
		UKnots, VKnots,
		UMults, VMults,
		curve1->Degree(), curve2->Degree());
}

void SurfaceModelingTool::Coons_G1(Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4, Handle(Geom_BSplineCurve)& curve1_derivative, Handle(Geom_BSplineCurve)& curve2_derivative, Handle(Geom_BSplineCurve)& curve3_derivative, Handle(Geom_BSplineCurve)& curve4_derivative, Handle(Geom_BSplineSurface)& mySurface_coons)
{
	// First make the opposite boundary curves compatible,
	// Second generate the ruled surface in the u/v direction and of the four corner points
	// Third generate the sum surface 

	//
	//1.1 Make Curve and the Derivative Compatible
	SetSameDistribution(curve1, curve1_derivative);
	SetSameDistribution(curve2, curve2_derivative);
	SetSameDistribution(curve3, curve3_derivative);
	SetSameDistribution(curve4, curve4_derivative);

	//1.2: Make opposite Curves Compatible
	//curve 1 and curve 3 are opposite curves, make them compatible
	SetSameDistribution(curve1, curve3);
	//curve 2 and curve 4 are opposite curves, make them compatible
	SetSameDistribution(curve2, curve4);

	//1.3 distribute to the derivatives
	//     
	SetSameDistribution(curve1, curve1_derivative);
	SetSameDistribution(curve2, curve2_derivative);
	SetSameDistribution(curve3, curve3_derivative);
	SetSameDistribution(curve4, curve4_derivative);


	//2: Generate the coons Bspline surface

	//Get the control point number in the u/v direction
	Standard_Integer NbUPoles = curve1->NbPoles();
	Standard_Integer NbVPoles = curve2->NbPoles();

	//control points of the coons surface
	TColgp_Array2OfPnt Poles_result(1, NbUPoles, 1, NbVPoles);

	//control points of G1 surface in the v direction
	TColgp_Array2OfPnt Poles_v(1, NbUPoles, 1, 4);

	//control points of G1 surface in the u direction
	TColgp_Array2OfPnt Poles_u(1, 4, 1, NbVPoles);

	//control points of G1 surface of the four corner points
	TColgp_Array2OfPnt Poles_u_v(1, 4, 1, 4);

	// Get the u knot vector
	Standard_Integer NbUKnot = curve1->NbKnots();
	TColStd_Array1OfReal    UKnots(1, NbUKnot);
	TColStd_Array1OfInteger UMults(1, NbUKnot);
	curve1->Knots(UKnots);
	curve1->Multiplicities(UMults);

	// Get the v knot vector
	Standard_Integer NbVKnot = curve2->NbKnots();
	TColStd_Array1OfReal    VKnots(1, NbVKnot);
	TColStd_Array1OfInteger VMults(1, NbVKnot);
	curve2->Knots(VKnots);
	curve2->Multiplicities(VMults);

	// Set v knots of the v ruled surface
	TColStd_Array1OfReal    VKnots_G1InVdirection(1, 2);
	TColStd_Array1OfInteger VMults_G1InVdirection(1, 2);
	VKnots_G1InVdirection(1) = curve2->FirstParameter();
	VKnots_G1InVdirection(2) = curve2->LastParameter();

	VMults_G1InVdirection(1) = 4;
	VMults_G1InVdirection(2) = 4;

	// Set u knots of the u ruled surface
	TColStd_Array1OfReal    UKnots_G1InUdirection(1, 2);
	TColStd_Array1OfInteger UMults_G1InUdirection(1, 2);
	UKnots_G1InUdirection(1) = curve1->FirstParameter();
	UKnots_G1InUdirection(2) = curve1->LastParameter();

	UMults_G1InUdirection(1) = 4;
	UMults_G1InUdirection(2) = 4;

	// 
	//2.1 generate the ruled surface in the v direction
	// 
	// 
	gp_Pnt p1, p2, p1d, p2d;
	// Set control points of the v ruled surface
	for (Standard_Integer i = 1; i <= NbUPoles; i++)
	{
		p1 = curve1->Pole(i);
		p1d = curve1_derivative->Pole(i);
		p2 = curve3->Pole(i);
		p2d = curve3_derivative->Pole(i);

		Poles_v(i, 1) = p1;
		Poles_v(i, 2).SetXYZ(p1.XYZ() + p1d.XYZ() / 3);
		Poles_v(i, 3).SetXYZ(p2.XYZ() - p2d.XYZ() / 3);
		Poles_v(i, 4) = p2;
	}

	// Generate the bspline surface
	Handle(Geom_BSplineSurface) mySurface_V = new Geom_BSplineSurface(Poles_v,
		UKnots, VKnots_G1InVdirection,
		UMults, VMults_G1InVdirection,
		curve1->Degree(), 3);

	// Make the surface compatible with the result coons surface
	mySurface_V->IncreaseDegree(curve1->Degree(), curve2->Degree());

	for (Standard_Integer i = 1; i <= NbVKnot; i++)
		mySurface_V->InsertVKnot(VKnots(i), VMults(i), Precision::Confusion());

	// 
	//2.2 generate the ruled surface in the u direction
	// 

	// Set control points of the u ruled surface
	for (Standard_Integer i = 1; i <= NbVPoles; i++)
	{
		p1 = curve4->Pole(i);
		p1d = curve4_derivative->Pole(i);
		p2 = curve2->Pole(i);
		p2d = curve2_derivative->Pole(i);

		Poles_u(1, i) = p1;
		Poles_u(2, i).SetXYZ(p1.XYZ() + p1d.XYZ() / 3);
		Poles_u(3, i).SetXYZ(p2.XYZ() - p2d.XYZ() / 3);
		Poles_u(4, i) = p2;
	}

	Handle(Geom_BSplineSurface) mySurface_U = new Geom_BSplineSurface(Poles_u,
		UKnots_G1InUdirection, VKnots,
		UMults_G1InUdirection, VMults,
		3, curve2->Degree());

	//increase degree
	mySurface_U->IncreaseDegree(curve1->Degree(), curve2->Degree());

	for (Standard_Integer i = 1; i <= NbUKnot; i++)
		mySurface_U->InsertUKnot(UKnots(i), UMults(i), Precision::Confusion());

	// 
	//2.3 generate the ruled surface of the four corner points
	// 
	// 
	gp_Vec T00_1, T00_2, T10_1, T10_2, T01_1, T01_2, T11_1, T11_2;
	gp_Pnt point;
	curve1_derivative->D1(curve1_derivative->FirstParameter(), point, T00_1);
	curve4_derivative->D1(curve4_derivative->FirstParameter(), point, T00_2);
	//Twist Compatiblity: T00_1 should coicide with T00_2
	curve1_derivative->D1(curve1_derivative->LastParameter(), point, T10_1);
	curve2_derivative->D1(curve2_derivative->FirstParameter(), point, T10_2);

	curve3_derivative->D1(curve3_derivative->FirstParameter(), point, T01_1);
	curve4_derivative->D1(curve4_derivative->LastParameter(), point, T01_2);

	curve2_derivative->D1(curve2_derivative->LastParameter(), point, T11_1);
	curve3_derivative->D1(curve3_derivative->LastParameter(), point, T11_2);

	Poles_u_v(1, 1) = curve1->Pole(1);
	Poles_u_v(2, 1).SetXYZ(curve1->Pole(1).XYZ() + curve4_derivative->Pole(1).XYZ() / 3);
	Poles_u_v(3, 1).SetXYZ(curve1->Pole(NbUPoles).XYZ() - curve2_derivative->Pole(1).XYZ() / 3);
	Poles_u_v(4, 1) = curve1->Pole(NbUPoles);

	Poles_u_v(1, 2).SetXYZ(curve1->Pole(1).XYZ() + curve1_derivative->Pole(1).XYZ() / 3);
	//Poles_u_v(2, 2).SetXYZ(curve1->Pole(1).XYZ() + curve1_derivative->Pole(1).XYZ() / 3);;
	//Poles_u_v(3, 2) = curve3->Pole(1);
	Poles_u_v(4, 2).SetXYZ(curve1->Pole(NbUPoles).XYZ() + curve1_derivative->Pole(NbUPoles).XYZ() / 3);

	Poles_u_v(1, 3).SetXYZ(curve3->Pole(1).XYZ() - curve3_derivative->Pole(1).XYZ() / 3);
	//Poles_u_v(2, 3).SetXYZ(curve1->Pole(1).XYZ() + curve1_derivative->Pole(1).XYZ() / 3);;
	//Poles_u_v(3, 3) = curve3->Pole(1);
	Poles_u_v(4, 3).SetXYZ(curve3->Pole(NbUPoles).XYZ() - curve3_derivative->Pole(NbUPoles).XYZ() / 3);

	Poles_u_v(1, 4).SetXYZ(curve3->Pole(1).XYZ());
	Poles_u_v(2, 4).SetXYZ(curve3->Pole(1).XYZ() + curve4_derivative->Pole(NbVPoles).XYZ() / 3);
	Poles_u_v(3, 4).SetXYZ(curve3->Pole(NbUPoles).XYZ() - curve2_derivative->Pole(NbVPoles).XYZ() / 3);
	Poles_u_v(4, 4).SetXYZ(curve3->Pole(NbUPoles).XYZ());

	Poles_u_v(2, 2).SetXYZ(Poles_u_v(2, 1).XYZ() + Poles_u_v(1, 2).XYZ() - Poles_u_v(1, 1).XYZ() + (T00_1.XYZ() + T00_2.XYZ()) / 18);
	Poles_u_v(3, 2).SetXYZ(Poles_u_v(3, 1).XYZ() + Poles_u_v(4, 2).XYZ() - Poles_u_v(4, 1).XYZ() + (T10_1.XYZ() + T10_2.XYZ()) / 18);
	Poles_u_v(2, 3).SetXYZ(Poles_u_v(1, 3).XYZ() + Poles_u_v(2, 4).XYZ() - Poles_u_v(1, 4).XYZ() + (T01_1.XYZ() + T01_2.XYZ()) / 18);
	Poles_u_v(3, 3).SetXYZ(Poles_u_v(3, 4).XYZ() + Poles_u_v(4, 3).XYZ() - Poles_u_v(4, 4).XYZ() + (T11_1.XYZ() + T11_2.XYZ()) / 18);

	Handle(Geom_BSplineSurface) mySurface_G1Surfaceof4CornerPoints = new Geom_BSplineSurface(Poles_u_v,
		UKnots_G1InUdirection, VKnots_G1InVdirection,
		UMults_G1InUdirection, VMults_G1InVdirection,
		3, 3);

	mySurface_G1Surfaceof4CornerPoints->IncreaseDegree(curve1->Degree(), curve2->Degree());

	for (Standard_Integer i = 1; i <= NbVKnot; i++)
		mySurface_G1Surfaceof4CornerPoints->InsertVKnot(VKnots(i), VMults(i), Precision::Confusion());

	for (Standard_Integer i = 1; i <= NbUKnot; i++)
		mySurface_G1Surfaceof4CornerPoints->InsertUKnot(UKnots(i), UMults(i), Precision::Confusion());

	//2.4 Generate the sum surface mySurface_VRuled + mySurface_URuled - mySurface_RuledSurfaceof4CornerPoints 
	//increase degree

	for (Standard_Integer i = 1; i <= NbUPoles; i++)
	{
		for (Standard_Integer j = 1; j <= NbVPoles; j++)
		{
			Poles_result(i, j).SetXYZ(mySurface_V->Pole(i, j).XYZ() + mySurface_U->Pole(i, j).XYZ() - mySurface_G1Surfaceof4CornerPoints->Pole(i, j).XYZ());
		}
	}

	mySurface_coons = new Geom_BSplineSurface(Poles_result,
		UKnots, VKnots,
		UMults, VMults,
		curve1->Degree(), curve2->Degree());
}

Standard_Integer SurfaceModelingTool:: Arrange_Coons_G0(std::vector<Handle(Geom_BSplineCurve)>& curveArray, Handle(Geom_BSplineCurve)& bslpineCurve1, Handle(Geom_BSplineCurve)& bslpineCurve2, Handle(Geom_BSplineCurve)& bslpineCurve3, Handle(Geom_BSplineCurve)& bslpineCurve4, double Tol, Standard_Integer IsModify)
{
	//Standard_Real Tol = 2;
	//Currently only work for four curves
	std::vector<Handle(Geom_BSplineCurve)> curveArraybak(curveArray);

	if (curveArraybak.size() != 4)
		return 0;

	bslpineCurve1 = curveArraybak[0];
	bslpineCurve2 = curveArraybak[1];
	bslpineCurve3 = curveArraybak[2];
	bslpineCurve4 = curveArraybak[3];

	Standard_Real maxDis = -1;
	Standard_Real dis;
	gp_Pnt curve1startpoint = bslpineCurve1->StartPoint();
	gp_Pnt curve1endpoint = bslpineCurve1->EndPoint();

	curveArraybak.erase(curveArraybak.begin());
	gp_Pnt curve2startpoint = bslpineCurve2->StartPoint();
	gp_Pnt curve2endpoint = bslpineCurve2->EndPoint();

	gp_Pnt curve3startpoint = bslpineCurve3->StartPoint();
	gp_Pnt curve3endpoint = bslpineCurve3->EndPoint();

	gp_Pnt curve4startpoint = bslpineCurve4->StartPoint();
	gp_Pnt curve4endpoint = bslpineCurve4->EndPoint();

	gp_Pnt point;
	//fine curve2
	for (Standard_Integer i = 0; i < curveArraybak.size(); i++)
	{
		bslpineCurve2 = curveArraybak[i];

		curve2startpoint = bslpineCurve2->StartPoint();
		curve2endpoint = bslpineCurve2->EndPoint();

		if (curve1endpoint.IsEqual(curve2startpoint, Tol))
		{
			dis = curve1endpoint.Distance(curve2startpoint);
			if (dis > maxDis)
				maxDis = dis;

			if (IsModify && dis < Tol)
			{
				curve1endpoint.BaryCenter(0.5, curve2startpoint, 0.5);
				bslpineCurve2->SetPole(1, curve1endpoint);
				bslpineCurve1->SetPole(bslpineCurve1->NbPoles(), curve1endpoint);
			}

			curveArraybak.erase(curveArraybak.begin() + i);
			break;
		}
		if (curve1endpoint.IsEqual(curve2endpoint, Tol))
		{
			dis = curve1endpoint.Distance(curve2endpoint);
			if (dis > maxDis)
				maxDis = dis;

			if (IsModify && dis < Tol)
			{
				curve1endpoint.BaryCenter(0.5, curve2endpoint, 0.5);
				bslpineCurve2->SetPole(bslpineCurve2->NbPoles(), curve1endpoint);
				bslpineCurve1->SetPole(bslpineCurve1->NbPoles(), curve1endpoint);
			}

			curveArraybak.erase(curveArraybak.begin() + i);
			bslpineCurve2->Reverse();

			curve2startpoint = bslpineCurve2->StartPoint();
			curve2endpoint = bslpineCurve2->EndPoint();
			break;
		}
	}

	//fine curve3
	for (Standard_Integer i = 0; i < curveArraybak.size(); i++)
	{
		bslpineCurve3 = curveArraybak[i];

		curve3startpoint = bslpineCurve3->StartPoint();
		curve3endpoint = bslpineCurve3->EndPoint();

		if (curve2endpoint.IsEqual(curve3endpoint, Tol))
		{
			dis = curve2endpoint.Distance(curve3endpoint);
			if (dis > maxDis)
				maxDis = dis;

			if (IsModify && dis < Tol)
			{
				curve2endpoint.BaryCenter(0.5, curve3endpoint, 0.5);
				bslpineCurve2->SetPole(bslpineCurve2->NbPoles(), curve2endpoint);
				bslpineCurve3->SetPole(bslpineCurve3->NbPoles(), curve2endpoint);
			}

			curveArraybak.erase(curveArraybak.begin() + i);
			break;
		}
		if (curve2endpoint.IsEqual(curve3startpoint, Tol))
		{
			dis = curve2endpoint.Distance(curve3startpoint);
			if (dis > maxDis)
				maxDis = dis;

			if (IsModify && dis < Tol)
			{
				curve2endpoint.BaryCenter(0.5, curve3startpoint, 0.5);
				bslpineCurve2->SetPole(bslpineCurve2->NbPoles(), curve2endpoint);
				bslpineCurve3->SetPole(1, curve2endpoint);
			}

			curveArraybak.erase(curveArraybak.begin() + i);
			bslpineCurve3->Reverse();

			curve3startpoint = bslpineCurve3->StartPoint();
			curve3endpoint = bslpineCurve3->EndPoint();
			break;
		}
	}
	//

	if (curveArraybak.size() != 1)
		return 0;

	bslpineCurve4 = curveArraybak[0];
	curve4startpoint = bslpineCurve4->StartPoint();
	curve4endpoint = bslpineCurve4->EndPoint();
	if (curve3startpoint.IsEqual(curve4endpoint, Tol))
	{
		dis = curve3startpoint.Distance(curve4endpoint);
		if (dis > maxDis)
			maxDis = dis;

		if (IsModify && dis < Tol)
		{
			curve3startpoint.BaryCenter(0.5, curve4endpoint, 0.5);
			bslpineCurve3->SetPole(1, curve3startpoint);
			bslpineCurve4->SetPole(bslpineCurve4->NbPoles(), curve3startpoint);
		}
	}
	else
		if (curve3startpoint.IsEqual(curve4startpoint, Tol))
		{
			dis = curve3startpoint.Distance(curve4startpoint);
			if (dis > maxDis)
				maxDis = dis;

			if (IsModify && dis < Tol)
			{
				curve3startpoint.BaryCenter(0.5, curve4startpoint, 0.5);
				bslpineCurve3->SetPole(1, curve3startpoint);
				bslpineCurve4->SetPole(1, curve3startpoint);
			}

			bslpineCurve4->Reverse();
			curve4startpoint = bslpineCurve4->StartPoint();
			curve4endpoint = bslpineCurve4->EndPoint();
		}
		else
			return 0;
	return 1;
}

void SurfaceModelingTool::ClassifyAndSortISOcurves(const std::vector<Handle(Geom_BSplineCurve)>&
	anISOcurvesArray,
	std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray,
	std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray)
{
	if (anISOcurvesArray.empty())
	{
		return;
	}

	// �������еȲ��ߵĲ����е㣬�����еȲ��ߵ����ĵ�
	double x = 0, y = 0, z = 0;
	for (const auto& curve : anISOcurvesArray)
	{
		gp_Pnt avgPnt = MathTool::ComputeAverageSamplePoint(curve, 10);
		x += avgPnt.X() / anISOcurvesArray.size();
		y += avgPnt.Y() / anISOcurvesArray.size();
		z += avgPnt.Z() / anISOcurvesArray.size();
	}
	// ���еȲ��ߵ����ĵ㣬�������Ϊ���������
	gp_Pnt middlePoint(x, y, z);
	// ȡ����һ��������ȷ��U����
	const Handle(Geom_BSplineCurve)& firstCurve = anISOcurvesArray.front();
	uISOcurvesArray.push_back(firstCurve);
	for (Standard_Integer i = 1; i < anISOcurvesArray.size(); i++)
	{
		const auto& curve = anISOcurvesArray[i];
		// ����Ȳ���֮�����С������� 0����ô����Ϊͬ��Ȳ��ߣ�������Ϊ u ��
		if (MathTool::ComputeCurveCurveDistance(curve, firstCurve) > Precision::Confusion())
			uISOcurvesArray.push_back(curve); // ��ӵ�U�����еȲ�������
		else
			vISOcurvesArray.push_back(curve); // ��ӵ�V�����еȲ�������
	}
	// ���ݵȲ��߾�������е�ľ����������
	// �����ҵ��߽磬Ȼ���ٸ���ÿһ���ߵ��߽�ľ���������
	double maxDist = RealFirst();
	Handle(Geom_BSplineCurve) uBoundary; // u ����߽�
	for (const auto& uCurve : uISOcurvesArray)
	{
		// ����ͶӰ����
		GeomAPI_ProjectPointOnCurve projector(middlePoint, uCurve);
		// ����Ƿ���ͶӰ���
		if (projector.NbPoints() > 0)
		{
			// ��ȡ������̵�ͶӰ��
			double distance = projector.LowerDistance();
			if (distance > maxDist)
			{
				maxDist = distance;
				uBoundary = uCurve;
			}
		}
	}
	maxDist = RealFirst();
	Handle(Geom_BSplineCurve) vBoundary;
	for (const auto& vCurve : vISOcurvesArray)
	{
		// ����ͶӰ����
		GeomAPI_ProjectPointOnCurve projector(middlePoint, vCurve);
		// ����Ƿ���ͶӰ���
		if (projector.NbPoints() > 0)
		{
			// ��ȡ������̵�ͶӰ��
			double distance = projector.LowerDistance();
			if (distance > maxDist)
			{
				maxDist = distance;
				vBoundary = vCurve;
			}
		}
	}
	// ��U���V��ĵȲ��߷ֱ��������
	std::sort(uISOcurvesArray.begin(), uISOcurvesArray.end(),
		[&](const Handle(Geom_BSplineCurve)& c1, const Handle(Geom_BSplineCurve)&
			c2) {
				return MathTool::ComputeCurveCurveDistance(c1, uBoundary) <
					MathTool::ComputeCurveCurveDistance(c2, uBoundary);
		});
	// ��V��ĵȲ��߸��ݾ���߽�ľ����������
	std::sort(vISOcurvesArray.begin(), vISOcurvesArray.end(),
		[&](const Handle(Geom_BSplineCurve)& c1, const Handle(Geom_BSplineCurve)&
			c2) {
				return MathTool::ComputeCurveCurveDistance(c1, vBoundary) <
					MathTool::ComputeCurveCurveDistance(c2, vBoundary);
		});
}

double CompareCurve(const std::vector<gp_Pnt>& points, std::vector<double>& params, Handle(Geom_BSplineCurve)& BsplineCurve
	, std::vector<double>& DisCovParams, bool isSingle, double Toler)
{
	if (points.size() != params.size())
	{
		std::cout << "SIZE NOT EQUAL :IN CXFUNC!" << std::endl;
		return 0;
	}
	double maxError1 = 0;
	double maxParam = 0;
	bool NeedIter = false;
	for (Standard_Integer i = 0; i < params.size(); i++) {
		gp_Pnt pnt = BsplineCurve->Value(params[i]);
		Standard_Real x, y, z;
		pnt.Coord(x, y, z);
		double disTemp = pnt.Distance(points[i]);
		if (disTemp > maxError1)
		{
			maxError1 = disTemp;
			maxParam = params[i];
		}
		if (!isSingle && disTemp > Toler)
		{
			NeedIter = true;
			DisCovParams.push_back(params[i]);
		}
	}
	if (isSingle && maxError1 > Toler)
	{
		DisCovParams.push_back(maxParam);
	}
	return maxError1;
}

Standard_Integer FindSpan(Standard_Integer n, Standard_Integer p, double u, const std::vector<double>& Knots)
{
	if (u == Knots[n + 1])
		return n;
	Standard_Integer low, hign;
	low = p;
	hign = n + 1;
	Standard_Integer mid = (low + hign) / 2;
	while (u < Knots[mid] || u >= Knots[mid + 1])
	{
		if (u < Knots[mid])
		{
			hign = mid;
		}
		else
		{
			low = mid;
		}
		mid = (low + hign) / 2;
	}
	return mid;
}

std::vector<double> KnotGernerationByMergeKnots(const std::vector<double>& KnotsA, const std::vector<double>& KnotsB)
{
	std::vector<double> C;
	std::vector<double> TempB;
	for (size_t i = 0; i < KnotsB.size(); i++)
	{
		if (isEqual(KnotsB[i], 1))
		{
			break;
		}
		if (!isEqual(KnotsB[i], 0))
		{
			TempB.push_back(KnotsB[i]);
		}
	}
	Standard_Integer i = 0, j = 0;
	while (i < KnotsA.size() && j < TempB.size()) {
		if (islessequal(KnotsA[i], TempB[j])) {
			C.push_back(KnotsA[i]);
			i++;
		}
		else {
			C.push_back(TempB[j]);
			j++;
		}
	}

	// ���ʣ��Ԫ��
	while (i < KnotsA.size()) {
		C.push_back(KnotsA[i]);
		i++;
	}
	while (j < TempB.size()) {
		C.push_back(TempB[j]);
		j++;
	}

	return C;
}

//Uniform Curve Knot to [0,1]
void UniformCurve(Handle(Geom_BSplineCurve)& curve)
{
	TColStd_Array1OfReal curveKnots(1, curve->NbKnots());
	curve->Knots(curveKnots);
	if (!(curveKnots(curveKnots.Lower()) == 0 && curveKnots(curveKnots.Upper()) == 1))
	{
		BSplCLib::Reparametrize(0, 1, curveKnots);
		curve->SetKnots(curveKnots);
	}
}

//To compute the value of a b-spline basic function value 
double OneBasicFun(double u, Standard_Integer i, Standard_Integer p, std::vector<double>& Knots)
{
	double Nip, uleft, uright, saved, temp;
	Standard_Integer m = Knots.size() - 1;
	std::vector<double>N(p + 1);
	if ((i == 0 && isEqual(u, Knots[0])) || (i == m - p - 1 && isEqual(u, Knots[m])))
	{
		return 1.0;
	}

	if (isLessThan(u, Knots[i]) || isGreaterThanOrEqual(u, Knots[i + p + 1]))
	{
		return 0.0;
	}
	for (size_t j = 0; j <= p; j++)
	{
		if (isGreaterThanOrEqual(u, Knots[i + j]) && isLessThan(u, Knots[i + j + 1]))
		{
			N[j] = 1.0;
		}
		else
		{
			N[j] = 0.0;
		}
	}
	for (size_t k = 1; k <= p; k++)
	{
		if (N[0] == 0.0)
		{
			saved = 0.0;
		}
		else
		{
			saved = ((u - Knots[i]) * N[0]) / (Knots[i + k] - Knots[i]);
		}
		for (size_t j = 0; j < p - k + 1; j++)
		{
			uleft = Knots[i + j + 1];
			uright = Knots[i + j + k + 1];
			if (N[j + 1] == 0.0)
			{
				N[j] = saved;
				saved = 0.0;
			}
			else
			{
				temp = N[j + 1] / (uright - uleft);
				N[j] = saved + (uright - u) * temp;
				saved = (u - uleft) * temp;
			}
		}
	}
	Nip = N[0];
	return Nip;
}

//To compute the Res Point Value
gp_Vec CalResPnt(Standard_Integer k, const std::vector<gp_Pnt>& dataPoints, std::vector<double>& parameters, Standard_Integer p,
	std::vector<double>& Knots, Standard_Integer CtrlPntNum) {
	Standard_Real aCoeff1 = OneBasicFun(parameters[k], 0, p, Knots);
	Standard_Real aCoeff2 = OneBasicFun(parameters[k], CtrlPntNum, p, Knots);
	gp_Vec vecTemp0(dataPoints[0].Coord());
	gp_Vec vecTempm(dataPoints[dataPoints.size() - 1].Coord());
	gp_Vec vecTempk(dataPoints[k].Coord());
	gp_Vec vectemp = vecTempk - aCoeff1 * vecTemp0 - aCoeff2 * vecTempm;
	return vectemp;
}

//To compute the Chord length parameterization
void ComputeChordLength(const std::vector<gp_Pnt>& points, std::vector<double>& parameters)
{
	parameters.push_back(0.0);  // Start parameter

	double totalLength = 0.0;
	std::vector<double> segmentLengths;

	// Calculate distances between consecutive points
	for (size_t i = 1; i < points.size(); ++i) {
		double dist = points[i - 1].Distance(points[i]);
		segmentLengths.push_back(dist);
		totalLength += dist;
	}

	// Calculate parameter for each point based on the cumulative length
	double cumulativeLength = 0.0;
	for (size_t i = 0; i < segmentLengths.size(); ++i) {
		cumulativeLength += segmentLengths[i];
		parameters.push_back(cumulativeLength / totalLength);  // Normalize by total length
	}
}

//To trans Sequence to Knots
void sequenceToKnots(const std::vector<double>& sequence, std::vector<double>& knots, std::vector<Standard_Integer>& multiplicities)
{
	if (sequence.empty()) return;

	std::map<double, Standard_Integer> knotMap;

	// ʹ��map��ͳ��ÿ���ڵ���ظ�����
	for (double value : sequence) 
	{
		bool found = false;
		for (auto& knot : knotMap)
		{
			if (isEqual(value, knot.first, 1e-6))
			{
				knot.second++;
				found = true;
				break;
			}
		}
		if (!found) {
			knotMap[value] = 1;
		}
	}

	// ��map������ת�Ƶ�knots��multiplicities����
	for (const auto& knot : knotMap) 
	{
		knots.push_back(knot.first);
		multiplicities.push_back(knot.second);
	}
}

//To approximate a Bspline curve by a set of points
//in this case, the knot vector(the number of contral point), the degree is known before approximate.
//Assume there are m data points.We need to use n( CtrlPntNum)+1 to approximate the m data points .
//CtrlPntNum n
Handle(Geom_BSplineCurve) ApproximateC(const std::vector<gp_Pnt>& Pnts, std::vector<double>& PntParams, TColStd_Array1OfReal& Knots, TColStd_Array1OfInteger& Mutis, std::vector<double>& FKnots,
	Standard_Integer degree)
{
	Standard_Integer n = FKnots.size() - degree - 2;
	Standard_Integer m = Pnts.size() - 1;
	//1.Chord Parameterized
	if (PntParams.size() == 0)
	{
		ComputeChordLength(Pnts, PntParams);
	}
	//(N^T*N)P=R

	//2.Construct matrix N
	Eigen::MatrixXd matN(m - 1, n - 1);
	for (size_t i = 0; i < m - 1; i++)
	{
		for (size_t j = 0; j < n - 1; j++)
		{
			double value = OneBasicFun(PntParams[i + 1], j + 1, degree, FKnots);
			matN(i, j) = value;
		}
	}
	Eigen::MatrixXd matTransposeN = matN.transpose();
	Eigen::MatrixXd NTN = matTransposeN * matN;

	//3.Construct matrix R
	//3.1 Compute Ri(2-(m-1))
	Eigen::VectorXd VRx(m - 1);
	Eigen::VectorXd VRy(m - 1);
	Eigen::VectorXd VRz(m - 1);
	for (size_t i = 1; i <= m - 1; i++)
	{
		gp_Vec VecTemp = CalResPnt(i, Pnts, PntParams, degree, FKnots, n);
		double x, y, z;
		VecTemp.Coord(x, y, z);
		VRx(i - 1) = x;
		VRy(i - 1) = y;
		VRz(i - 1) = z;
	}
	//3.2 Compute the component of R
	Eigen::VectorXd Rx = matTransposeN * VRx;
	Eigen::VectorXd Ry = matTransposeN * VRy;
	Eigen::VectorXd Rz = matTransposeN * VRz;

	//4.solve NtN P = R
	Eigen::VectorXd Sx = NTN.lu().solve(Rx);
	Eigen::VectorXd Sy = NTN.lu().solve(Ry);
	Eigen::VectorXd Sz = NTN.lu().solve(Rz);


	//5.construct bspline curve
	TColgp_Array1OfPnt ctrlPnts(1, n + 1);
	ctrlPnts[1] = Pnts[0];
	ctrlPnts[n + 1] = Pnts[m];
	for (size_t i = 2; i <= n; i++)
	{
		gp_Pnt pntTemp(Sx(i - 2), Sy(i - 2), Sz(i - 2));
		ctrlPnts[i] = pntTemp;
	}
	Handle(Geom_BSplineCurve) bspline = new Geom_BSplineCurve(ctrlPnts, Knots, Mutis, degree);
	return bspline;
}

//To approximate a Bspline curve by a set of points
	//in this case, the knot vector(the number of contral point), the degree is known before approximate.
	//Assume there are m data points.We need to use n( CtrlPntNum)+1 to approximate the m data points .
	//CtrlPntNum n
Handle(Geom_BSplineCurve) ApproximateC(const std::vector<gp_Pnt>& Pnts, std::vector<double>& params, std::vector<double>& FKnots, Standard_Integer degree)
{
	std::vector<double> Knots;
	std::vector<Standard_Integer> Mutis;
	sequenceToKnots(FKnots, Knots, Mutis);
	TColStd_Array1OfReal Knots_OCC(1, Knots.size());
	TColStd_Array1OfInteger Mutis_OCC(1, Mutis.size());
	for (size_t i = 0; i < Knots.size(); i++)
	{
		Knots_OCC[i + 1] = Knots[i];
		Mutis_OCC[i + 1] = Mutis[i];
	}
	return ApproximateC(Pnts, params, Knots_OCC, Mutis_OCC, FKnots, degree);
}

void SurfaceModelingTool::CreateLoftingSurface(const std::vector<Handle(Geom_BSplineCurve)>& curvesArray, const std::vector<gp_Vec>& normals, std::vector<TopoDS_Shape>& loftingSurfaces)
{
	double offsetDistance = 400;
	for (Standard_Integer i = 0; i < curvesArray.size(); i++)
	{
		TColgp_Array1OfPnt aPnts1(1, curvesArray[i]->NbPoles());
		TColgp_Array1OfPnt aPnts2(1, curvesArray[i]->NbPoles());

		for (Standard_Integer j = 0; j < curvesArray[i]->NbPoles(); j++)
		{
			aPnts1.SetValue(j + 1, curvesArray[i]->Pole(j + 1).Translated(normals[i] * offsetDistance));
			aPnts2.SetValue(j + 1, curvesArray[i]->Pole(j + 1).Translated(-normals[i] * offsetDistance));
		}

		// ���� B �������߶���
		Handle(Geom_BSplineCurve) aBSplineCurve1 = new Geom_BSplineCurve(
			aPnts1,
			curvesArray[i]->Knots(),
			curvesArray[i]->Multiplicities(),
			curvesArray[i]->Degree()
		);
		Handle(Geom_BSplineCurve) aBSplineCurve2 = new Geom_BSplineCurve(
			aPnts2,
			curvesArray[i]->Knots(),
			curvesArray[i]->Multiplicities(),
			curvesArray[i]->Degree()
		);

		BRepBuilderAPI_MakeWire aWireMaker0, aWireMaker1, aWireMaker2;
		BRepBuilderAPI_MakeVertex aStartVertex0(curvesArray[i]->Pole(1));
		BRepBuilderAPI_MakeVertex aEndVertex0(curvesArray[i]->Pole(curvesArray[i]->NbPoles()));
		TopoDS_Edge aEdge0 = BRepBuilderAPI_MakeEdge(curvesArray[i], aStartVertex0, aEndVertex0);
		aWireMaker0.Add(aEdge0);

		BRepBuilderAPI_MakeVertex aStartVertex1(aBSplineCurve1->Pole(1));
		BRepBuilderAPI_MakeVertex aEndVertex1(aBSplineCurve1->Pole(aBSplineCurve1->NbPoles()));
		TopoDS_Edge aEdge1 = BRepBuilderAPI_MakeEdge(aBSplineCurve1, aStartVertex1, aEndVertex1);
		aWireMaker1.Add(aEdge1);

		BRepBuilderAPI_MakeVertex aStartVertex2(aBSplineCurve2->Pole(1));
		BRepBuilderAPI_MakeVertex aEndVertex2(aBSplineCurve2->Pole(aBSplineCurve2->NbPoles()));
		TopoDS_Edge aEdge2 = BRepBuilderAPI_MakeEdge(aBSplineCurve2, aStartVertex2, aEndVertex2);
		aWireMaker2.Add(aEdge2);

		BRepOffsetAPI_ThruSections aLoftSurface(Standard_True, Standard_False);
		aLoftSurface.AddWire(aWireMaker1.Wire());
		aLoftSurface.AddWire(aWireMaker0.Wire());
		aLoftSurface.AddWire(aWireMaker2.Wire());

		aLoftSurface.Build();
		if (!aLoftSurface.IsDone())
		{
			// ����쳣
		}
		else
		{
			loftingSurfaces.emplace_back(aLoftSurface.Shape());
		}
	}
}

Handle(Geom_BSplineCurve) IterateApproximate(std::vector<double>& InsertKnots, const std::vector<gp_Pnt>& Pnts, std::vector<double>& PntsParams, std::vector<double>& InitKnots, Standard_Integer degree, Standard_Integer MaxIterNum, double toler)
{
	Standard_Integer itNum = 1;
	double currentMaxError = 100;
	Handle(Geom_BSplineCurve) IterBspineCurve;
	std::vector<double> CurrentKnots = InitKnots;
	std::cout << std::endl;
	std::cout << "----------------Strat Iteration----------------" << std::endl;
	std::cout << "----------------Iteration Infor----------------" << std::endl;
	std::cout << "Max Iteration Num = " << MaxIterNum << " " << std::endl;
	std::cout << "Iteration degree = " << degree << " " << std::endl;
	std::cout << "Iteration toler = " << toler << " " << std::endl;
	std::cout << "----------------Strat----------------" << std::endl;
	while (currentMaxError > toler && itNum <= MaxIterNum)
	{
		auto start = std::chrono::high_resolution_clock::now();

		IterBspineCurve = ApproximateC(Pnts, PntsParams, CurrentKnots, degree);

		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> duration = end - start;

		KnotUpdate knotUpdate(IterBspineCurve, CurrentKnots, Pnts, PntsParams);
		auto newKnot = knotUpdate.SelfSingleUpdate(PARAM_BASED_BY_INTERVAL_ERROR);
		currentMaxError = knotUpdate.getMaxError();
		std::cout << "Iteration Num = " << itNum << " The Current Error = " << currentMaxError << " The Target Error = " << toler << " The Time = " << duration.count() << " ms " << std::endl;
		if (currentMaxError > toler)
		{
			//update knot vector
			InsertKnots.push_back(newKnot);
			CurrentKnots = knotUpdate.getSequences();
		}
		else
		{
			std::cout << "----------------End----------------" << std::endl;
			std::cout << "Find Curve,iteration ending--------------the current error is " << currentMaxError << std::endl;
			std::cout << "----------------End Iteration----------------" << std::endl << std::endl;
			return IterBspineCurve;
		}
		itNum++;
	}
	//set check params
	std::cout << "----------------End----------------" << std::endl;
	std::cout << "Arrive max iteration number, iteration ending--------------the current error is " << currentMaxError << std::endl;
	std::cout << "----------------End Iteration----------------" << std::endl << std::endl;

	//return bspline;
	return IterBspineCurve;
}
double CalPointsChordLen(const std::vector<gp_Pnt>& points)
{
	if (points.size() < 2) {
		return 0.0;
	}
	double totalLength = 0.0;
	for (size_t i = 1; i < points.size(); ++i) {
		totalLength += points[i - 1].Distance(points[i]);
	}
	return totalLength;
}
gp_Vec CalTangent(const std::vector<gp_Pnt>& points, double TanMagnitude)
{
	if (points.size() != 3) {
		throw std::invalid_argument("��������������Ϊ3����");
	}

	gp_Pnt P0 = points[0];
	gp_Pnt P1 = points[1];
	gp_Pnt P2 = points[2];

	// B1 = 2 * P1 - 0.5 * P0 - 0.5 * P2
	gp_Vec vec_P0_P1(P0, P1); // P1 - P0
	gp_Vec vec_P2_P1(P2, P1); // P1 - P2

	gp_Pnt B1(
		2.0 * P1.X() - 0.5 * P0.X() - 0.5 * P2.X(),
		2.0 * P1.Y() - 0.5 * P0.Y() - 0.5 * P2.Y(),
		2.0 * P1.Z() - 0.5 * P0.Z() - 0.5 * P2.Z()
	);

	// ���������� B'(0) = 2 * (B1 - B0)
	gp_Vec tangentVec(P0, B1); // B1 - B0
	tangentVec.Multiply(2.0);

	double magnitude = tangentVec.Magnitude();
	if (magnitude == 0.0) {
		throw std::runtime_error("�������ĳ���Ϊ�㣬�޷����š�");
	}

	tangentVec.Normalize();
	tangentVec.Multiply(TanMagnitude);

	return tangentVec;
}
void SurfaceModelingTool::LoftSurfaceIntersectWithCurve(const std::vector<TopoDS_Shape>& LoftingSur, 
	const std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_Initial, 
	const std::vector<Handle(Geom_BSplineCurve)>& anInternalBSplineCurves, 
	std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_New, 
	Standard_Integer isoCount,
	std::vector<std::vector<gp_Pnt>>& InterpolatePoints,
	std::vector<TopoDS_Edge>& TangentArray1,
	std::vector<TopoDS_Edge>& TangentArray2,
	Handle(Geom_BSplineSurface) CoonsSurface)
{
	std::vector<std::vector<gp_Pnt>> debugPoints;
	for (Standard_Integer i = 0; i < LoftingSur.size(); i++)
	{
		std::vector<gp_Pnt> aPntsVector; // ��̬�洢����Ͷ˵�
		std::vector<gp_Pnt> aInterPnts;  // ÿ��LoftingSur���ڲ��ߵĽ���

		// ����Lofting����
		TopExp_Explorer explorer2(LoftingSur[i], TopAbs_FACE);
		for (; explorer2.More(); explorer2.Next())
		{
			const TopoDS_Face& aFace = TopoDS::Face(explorer2.Current());
			Handle(Geom_Surface) aSur = BRep_Tool::Surface(aFace);
			Handle(Geom_BSplineSurface) aLoftingSur = Handle(Geom_BSplineSurface)::DownCast(aSur);

			if (!aLoftingSur.IsNull())
			{
				// �����ڲ���B�������ߣ����㽻��
				for (Standard_Integer j = 0; j < anInternalBSplineCurves.size(); j++)
				{
					GeomAPI_IntCS anInterCS(anInternalBSplineCurves[j], aLoftingSur);
					if (anInterCS.IsDone())
					{
						if (anInterCS.NbPoints())
						{
							for (Standard_Integer k = 1; k <= anInterCS.NbPoints(); k++)
							{
								aInterPnts.emplace_back(anInterCS.Point(k));
							}
						}
					}
				}
			}
		}

		// ������㡢������յ�
		if (aInterPnts.size() > 0)
		{
			gp_Pnt startPoint = ISOcurvesArray_Initial[i]->StartPoint();
			gp_Pnt endPoint = ISOcurvesArray_Initial[i]->EndPoint();
			aPntsVector.insert(aPntsVector.end(), aInterPnts.begin(), aInterPnts.end()); // ����

			// ������startPoint�ľ����С����������Ҫ�±�֤����
			std::sort(aPntsVector.begin(), aPntsVector.end(), [&startPoint](const gp_Pnt& p1, const gp_Pnt& p2)
			{
					return p1.Distance(startPoint) < p2.Distance(startPoint);
			});

			aPntsVector.insert(aPntsVector.begin(), startPoint);  // ���
			aPntsVector.push_back(endPoint);    // �յ�

			for (Standard_Integer j = 1; j < aPntsVector.size() - 1; j++)
			{
				gp_Pnt lastPnt = aPntsVector[j - 1];
				gp_Pnt Pnt = aPntsVector[j];
				gp_Pnt nextPnt = aPntsVector[j + 1];
				if (!(Pnt.Distance(lastPnt) > startPoint.Distance(endPoint) / (isoCount * 2) &&
					Pnt.Distance(nextPnt) > startPoint.Distance(endPoint) / (isoCount * 2)))
				{
					aPntsVector.erase(aPntsVector.begin() + j);
					j--;
				}
			}
			InterpolatePoints.push_back(aPntsVector);

			// startPoint���ߺ�endPoint������
			gp_Vec FirstD1, LastD1;
			gp_Vec startDirection = aPntsVector[1].XYZ() - aPntsVector[0].XYZ();
			startDirection.Normalize();
			gp_Vec endDirection = aPntsVector[aPntsVector.size() - 1].XYZ() - aPntsVector[aPntsVector.size() - 2].XYZ();
			endDirection.Normalize();

			Standard_Real uMin, uMax, vMin, vMax;
			CoonsSurface->Bounds(uMin, uMax, vMin, vMax); // ��ȡ������Χ
			if (!CoonsSurface.IsNull())
			{
				gp_Vec U_StartTangent, V_StartTangent;
				GeomAPI_ProjectPointOnSurf projector(aPntsVector.front(), CoonsSurface);

				if (projector.NbPoints() > 0)
				{
					gp_Pnt closestPoint = projector.NearestPoint();
					double uParam, vParam;
					projector.LowerDistanceParameters(uParam, vParam);

					// ��ȡU��V�����ƫ��������������
					CoonsSurface->D1(uParam, vParam, closestPoint, U_StartTangent, V_StartTangent);

					// ����Ƿ���Ҫ��ֵ�ƽ���
					if (!U_StartTangent.Magnitude()|| !V_StartTangent.Magnitude())
					{
						// ȷ�� uParam ����Ч������
						if (uParam - Precision::Confusion() < uMin) 
						{
							uParam += Precision::Confusion();  // ��������
						}
						if (uParam + Precision::Confusion() > uMax) 
						{
							uParam -= Precision::Confusion();  // ��С�������
						}
						else 
						{
							uParam += Precision::Confusion();  // ��������
						}

						// ȷ�� vParam ����Ч������
						if (vParam - Precision::Confusion() < vMin) 
						{
							vParam += Precision::Confusion();  // ��������
						}
						if (vParam + Precision::Confusion() > vMax)
						{
							vParam -= Precision::Confusion();  // ��С�������
						}
						else 
						{
							vParam += Precision::Confusion();  // ��������
						}
						CoonsSurface->D1(uParam, vParam, closestPoint, U_StartTangent, V_StartTangent);
					}
					// ���㷨������U��Vƫ�����Ĳ����
					gp_Vec normalVector = U_StartTangent.Crossed(V_StartTangent);
					normalVector = normalVector.Normalized();
					if (startDirection.Dot(normalVector) < 0)
					{
						normalVector.Reverse(); // ��������
					}

					// ����startDirection�ڷ������ϵ�ͶӰ����
					double projectionMagnitude = startDirection.Dot(normalVector) / normalVector.Magnitude();
					gp_Vec projection = normalVector * projectionMagnitude; // ����ͶӰ����

					// ����startDirection��ƽ���ϵ�ͶӰ��ȥ���������ķ�����
					gp_Vec projectedDirection = startDirection - projection;
					projectedDirection.Normalize();


					// ����ͶӰ�������
					FirstD1 = projectedDirection;
					FirstD1.Multiply(CalPointsChordLen(aPntsVector) / FirstD1.Magnitude());

					// �������߽߱�
					TangentArray1.push_back(BRepBuilderAPI_MakeEdge(closestPoint, closestPoint.Translated(FirstD1 * 0.1)).Edge());
				}
			}

			if (!CoonsSurface.IsNull())
			{
				gp_Vec U_EndTangent, V_EndTangent;
				GeomAPI_ProjectPointOnSurf projector(aPntsVector.back(), CoonsSurface);

				if (projector.NbPoints() > 0)
				{
					gp_Pnt closestPoint = projector.NearestPoint();
					double uParam, vParam;
					projector.LowerDistanceParameters(uParam, vParam);

					// ��ȡU��V�����ƫ��������������
					CoonsSurface->D1(uParam, vParam, closestPoint, U_EndTangent, V_EndTangent);
					// ����Ƿ���Ҫ��ֵ�ƽ���
					if (!U_EndTangent.Magnitude() || !V_EndTangent.Magnitude())
					{
						// ȷ�� uParam ����Ч������
						if (uParam - Precision::Confusion() < uMin) 
						{
							uParam += Precision::Confusion();  // ��������
						}
						if (uParam + Precision::Confusion() > uMax) 
						{
							uParam -= Precision::Confusion();  // ��С�������
						}
						else
						{
							uParam += Precision::Confusion();  // ��������
						}

						// ȷ�� vParam ����Ч������
						if (vParam - Precision::Confusion() < vMin)
						{
							vParam += Precision::Confusion();  // ��������
						}
						if (vParam + Precision::Confusion() > vMax) 
						{
							vParam -= Precision::Confusion();  // ��С�������
						}
						else
						{
							vParam += Precision::Confusion();  // ��������
						}
						CoonsSurface->D1(uParam, vParam, closestPoint, U_EndTangent, V_EndTangent);
					}
					// ���㷨������U��Vƫ�����Ĳ����
					gp_Vec normalVector = U_EndTangent.Crossed(V_EndTangent);
					normalVector = normalVector.Normalized();
					if (normalVector.Dot(endDirection) < 0)
					{
						normalVector.Reverse();
					}
					// ����endDirection�ڷ������ϵ�ͶӰ����
					double projectionMagnitude = endDirection.Dot(normalVector) / normalVector.Magnitude();
					gp_Vec projection = normalVector * projectionMagnitude; // ����ͶӰ����

					// ����endDirection��ƽ���ϵ�ͶӰ��ȥ���������ķ�����
					gp_Vec projectedDirection = endDirection - projection;
					projectedDirection.Normalize();
					// ����ͶӰ�������
					LastD1 = projectedDirection;
					LastD1.Multiply(CalPointsChordLen(aPntsVector) / LastD1.Magnitude());

					// �������߽߱�
					TangentArray1.push_back(BRepBuilderAPI_MakeEdge(closestPoint, closestPoint.Translated(LastD1 * 0.1)).Edge());
				}
			}

			// ����startDirection��endDirection�ĳ���
			startDirection.Multiply(CalPointsChordLen(aPntsVector) / startDirection.Magnitude());
			endDirection.Multiply(CalPointsChordLen(aPntsVector) / endDirection.Magnitude());

			Handle(TColgp_HArray1OfPnt) points = new TColgp_HArray1OfPnt(1, aPntsVector.size());
			for (Standard_Integer j = 0; j < aPntsVector.size(); j++)
			{
				points->SetValue(j + 1, aPntsVector[j]);
			}
			debugPoints.push_back(aPntsVector);
			GeomAPI_Interpolate interpolate(points, Standard_False, 0.1);
			//interpolate.Load(FirstD1, LastD1, Standard_True);
			// ִ�в�ֵ����
			interpolate.Perform();
			// ����Ƿ�ɹ���ɲ�ֵ
			if (interpolate.IsDone())
			{
				// ��ȡ��ֵ������߶���
				Handle(Geom_BSplineCurve) aBSplineCurve = interpolate.Curve();
				ISOcurvesArray_New.emplace_back(aBSplineCurve);
			}
		}
		else
		{
			ISOcurvesArray_New.emplace_back(ISOcurvesArray_Initial[i]);
			aPntsVector.emplace_back(ISOcurvesArray_Initial[i]->StartPoint());
			aPntsVector.emplace_back(ISOcurvesArray_Initial[i]->EndPoint());
			InterpolatePoints.emplace_back(aPntsVector);
		}

	}
}

std::vector<double> GetKnotsSequence(Handle(Geom_BSplineCurve) curve)
{
	std::vector<double> debugKnots;
	TColStd_Array1OfReal KnotsSequence = curve->KnotSequence();
	for (Standard_Integer i = KnotsSequence.Lower(); i <= KnotsSequence.Upper(); i++)
	{
		debugKnots.push_back(KnotsSequence.Value(i));
	}
	return debugKnots;
}

std::vector<double> ComputeUniformParam(Standard_Integer numSamples, double left, double right) {
	std::vector<double> parameters;
	if (numSamples == 0) {
		return parameters;
	}
	for (Standard_Integer i = 1; i <= numSamples; i++) {
		Standard_Real param = left + (right - left) * (i - 1) / (numSamples - 1);
		parameters.push_back(param);
	}
	return parameters;
}

std::vector<double> KnotGernerationByParams(const std::vector<double>& params, Standard_Integer n, Standard_Integer p)
{
	Standard_Integer m = params.size() - 1;
	double d = (m + 1) / (n - p + 1);
	std::vector<double> Knots(n + p + 2);
	Standard_Integer temp;
	double alpha;
	for (size_t i = 0; i <= p; i++)
	{
		Knots[i] = 0.0;
	}
	for (size_t j = 1; j <= n - p; j++)
	{
		temp = Standard_Integer(j * d);
		alpha = j * d - temp;
		Knots[p + j] = (1 - alpha) * params[temp - 1] + alpha * params[temp];
	}
	for (size_t i = n + 1; i <= n + p + 1; i++)
	{
		Knots[i] = 1;
	}
	return Knots;
}
//To compute the Res Point Value - first d1 case
gp_Vec CalResPnt(Standard_Integer k, const std::vector<gp_Pnt>& dataPoints, const gp_Pnt& SecondPoint, const gp_Pnt& LastSecondPoint, const std::vector<double>& parameters, Standard_Integer p,
	std::vector<double>& Knots, Standard_Integer CtrlPntNum) {
	Standard_Real aCoeff0 = OneBasicFun(parameters[k], 0, p, Knots);
	Standard_Real aCoeff1 = OneBasicFun(parameters[k], 1, p, Knots);
	Standard_Real aCoeffms1 = OneBasicFun(parameters[k], CtrlPntNum - 1, p, Knots);
	Standard_Real aCoeffm = OneBasicFun(parameters[k], CtrlPntNum, p, Knots);
	gp_Vec vecTemp0(dataPoints[0].Coord());
	gp_Vec vecTemp1(SecondPoint.Coord());
	gp_Vec vecTempms1(LastSecondPoint.Coord());
	gp_Vec vecTempm(dataPoints.back().Coord());
	gp_Vec vecTempk(dataPoints[k].Coord());
	gp_Vec vectemp = vecTempk - aCoeff0 * vecTemp0 - aCoeff1 * vecTemp1 - aCoeffms1 * vecTempms1 - aCoeffm * vecTempm;
	return vectemp;
}

Handle(Geom_BSplineCurve) EndDerivaConstraintBsplineCurveAppro(const std::vector<gp_Pnt>& Pnts, const gp_Vec& FirstD1, const gp_Vec& LastD1,
	std::vector<double>& Params, std::vector<double>& KnotSequences, Standard_Integer degree) {
	//check the number must large than 3
	if (Pnts.size() < 3) {
		std::cerr << "the number of data points must great than 3!" << std::endl;
		return nullptr;
	}
	//compute the second point
	gp_Vec first_point_vectype = gp_Vec(Pnts[0].XYZ());
	gp_Vec second_point_vectype = first_point_vectype + (KnotSequences[degree + 1] / (double)degree) * FirstD1;
	gp_Pnt second_point = gp_Pnt(second_point_vectype.XYZ());

	//compute the second point
	Standard_Integer index = KnotSequences.size() - degree - 2;
	gp_Vec last_point_vectype = gp_Vec(Pnts.back().XYZ());
	gp_Vec last_second_point_vectype = last_point_vectype - ((1 - KnotSequences[index]) / (double)degree) * LastD1;
	gp_Pnt last_second_point = gp_Pnt(last_second_point_vectype.XYZ());

	Standard_Integer n = KnotSequences.size() - degree - 2;
	Standard_Integer m = Pnts.size() - 1;

	// Construct matrix N
	Eigen::MatrixXd matN = Eigen::MatrixXd::Zero(m - 1, n - 3);
	for (Standard_Integer i = 0; i < m - 1; ++i) {
		for (Standard_Integer j = 0; j < n - 3; ++j) {
			matN(i, j) = OneBasicFun(Params[i + 1], j + 2, degree, KnotSequences);
		}
	}

	// Construct matrix R for x, y, z components
	Eigen::MatrixXd VR(3, m - 1);
	for (Standard_Integer i = 1; i <= m - 1; ++i) {
		gp_Vec VecTemp = CalResPnt(i, Pnts, second_point, last_second_point, Params, degree, KnotSequences, n);
		double x, y, z;
		VecTemp.Coord(x, y, z);
		VR(0, i - 1) = x;
		VR(1, i - 1) = y;
		VR(2, i - 1) = z;
	}

	// Solve NtN * P = R for x, y, z components
	Eigen::MatrixXd S = matN.householderQr().solve(VR.transpose()).transpose();

	// Construct BSpline control points
	TColgp_Array1OfPnt ctrlPnts(1, n + 1);
	ctrlPnts.SetValue(1, Pnts[0]);
	ctrlPnts.SetValue(2, second_point);
	ctrlPnts.SetValue(n, last_second_point);
	ctrlPnts.SetValue(n + 1, Pnts[m]);
	for (Standard_Integer i = 3; i <= n - 1; ++i) {
		gp_Pnt pntTemp(S(0, i - 3), S(1, i - 3), S(2, i - 3));
		ctrlPnts.SetValue(i, pntTemp);
	}

	std::vector<double> Knots;
	std::vector<Standard_Integer> Mutis;
	sequenceToKnots(KnotSequences, Knots, Mutis);
	TColStd_Array1OfReal Knots_OCC(1, Knots.size());
	TColStd_Array1OfInteger Mutis_OCC(1, Mutis.size());
	for (size_t i = 0; i < Knots.size(); ++i) {
		Knots_OCC.SetValue(static_cast<Standard_Integer>(i + 1), Knots[i]);
		Mutis_OCC.SetValue(static_cast<Standard_Integer>(i + 1), Mutis[i]);
	}
	Handle(Geom_BSplineCurve) bspline = new Geom_BSplineCurve(ctrlPnts, Knots_OCC, Mutis_OCC, degree);
	return bspline;
}

Handle(Geom_BSplineCurve) IterateApproximate(std::vector<double>& InsertKnots,
	const std::vector<gp_Pnt>& Pnts, 
	const gp_Vec& FirstD1, const gp_Vec& LastD1, 
	std::vector<double>& PntsParams, std::vector<double>& InitKnots, 
	Standard_Integer degree, Standard_Integer MaxIterNum, double toler) {
	Standard_Integer itNum = 1;
	double currentMaxError = 100;
	Handle(Geom_BSplineCurve) IterBspineCurve;
	std::vector<double> CurrentKnots = InitKnots;
	std::cout << std::endl;
	std::cout << "----------------Strat Iteration----------------" << std::endl;
	std::cout << "----------------Iteration Infor----------------" << std::endl;
	std::cout << "Max Iteration Num = " << MaxIterNum << " " << std::endl;
	std::cout << "Iteration degree = " << degree << " " << std::endl;
	std::cout << "Iteration toler = " << toler << " " << std::endl;
	std::cout << "----------------Strat----------------" << std::endl;
	while (currentMaxError > toler && itNum <= MaxIterNum) {
		auto start = std::chrono::high_resolution_clock::now();

		IterBspineCurve = EndDerivaConstraintBsplineCurveAppro(Pnts, FirstD1, LastD1, PntsParams, CurrentKnots, degree);

		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> duration = end - start;

		KnotUpdate knotUpdate(IterBspineCurve, CurrentKnots, Pnts, PntsParams);
		auto newKnot = knotUpdate.SelfSingleUpdate(PARAM_BASED_BY_INTERVAL_ERROR);
		currentMaxError = knotUpdate.getMaxError();
		std::cout << "Iteration Num = " << itNum << " The Current Error = " << currentMaxError << " The Target Error = " << toler << " The Time = " << duration.count() << " ms " << std::endl;
		if (currentMaxError > toler) {
			//update knot vector
			InsertKnots.push_back(newKnot);
			CurrentKnots = knotUpdate.getSequences();
		}
		else {
			std::cout << "----------------End----------------" << std::endl;
			std::cout << "Find Curve,iteration ending--------------the current error is " << currentMaxError << std::endl;
			std::cout << "----------------End Iteration----------------" << std::endl << std::endl;
			return IterBspineCurve;
		}
		itNum++;
	}
	//set check params
	std::cout << "----------------End----------------" << std::endl;
	std::cout << "Arrive max iteration number, iteration ending--------------the current error is " << currentMaxError << std::endl;
	std::cout << "----------------End Iteration----------------" << std::endl << std::endl;

	//return bspline;
	return IterBspineCurve;
}

//To compute the Res Point Value - first d1 case
gp_Vec CalResPnt(Standard_Integer k, const std::vector<gp_Pnt>& dataPoints, const gp_Pnt& SecondPoint, const std::vector<double>& parameters, Standard_Integer p,
	std::vector<double>& Knots, Standard_Integer CtrlPntNum) {
	Standard_Real aCoeff0 = OneBasicFun(parameters[k], 0, p, Knots);
	Standard_Real aCoeff1 = OneBasicFun(parameters[k], 1, p, Knots);
	Standard_Real aCoeffm = OneBasicFun(parameters[k], CtrlPntNum, p, Knots);
	gp_Vec vecTemp0(dataPoints[0].Coord());
	gp_Vec vecTemp1(SecondPoint.Coord());
	gp_Vec vecTempm(dataPoints.back().Coord());
	gp_Vec vecTempk(dataPoints[k].Coord());
	gp_Vec vectemp = vecTempk - aCoeff0 * vecTemp0 - aCoeff1 * vecTemp1 - aCoeffm * vecTempm;
	return vectemp;
}

std::pair<double, double> processPoints(const gp_Pnt& P1, const gp_Pnt& P2, const std::vector<gp_Pnt>& isoInterpolatePoints, const std::vector<gp_Pnt>& oppsiteInterpolatePoints)
{
	std::vector<std::pair<double, gp_Pnt>> distancesP1;
	for (const auto& point : isoInterpolatePoints)
	{
		distancesP1.push_back({ P1.Distance(point), point });
	}

	std::sort(distancesP1.begin(), distancesP1.end(), [](std::pair<double, gp_Pnt> p1, std::pair<double, gp_Pnt> p2) {return p1.first < p2.first; });

	gp_Pnt nearestPoint1_P1 = distancesP1[0].second;
	gp_Pnt nearestPoint2_P1 = distancesP1[1].second;

	// ������������
	gp_Vec lineVec(nearestPoint1_P1, nearestPoint2_P1);
	gp_Vec vecToP1(nearestPoint1_P1, P1);

	// ����ͶӰ
	double projectionLength = vecToP1.Dot(lineVec) / lineVec.Magnitude();

	// �ж�ͶӰ�Ƿ�������������
	bool isOppositeSide = false;
	if (projectionLength >= 0 && projectionLength <= lineVec.Magnitude())
	{
		isOppositeSide = true;  // ���ͶӰ�����������ڣ���Ϊ����������
	}

	// ���������ͬһ�࣬����������һ����
	Standard_Integer next = 2;
	while (!isOppositeSide && next < distancesP1.size())  // ��������Ƿ�Խ��
	{
		nearestPoint2_P1 = distancesP1[next++].second;
		lineVec = gp_Vec(nearestPoint1_P1, nearestPoint2_P1);
		vecToP1 = gp_Vec(nearestPoint1_P1, P1);
		projectionLength = vecToP1.Dot(lineVec) / lineVec.Magnitude();

		if (projectionLength >= 0 && projectionLength <= lineVec.Magnitude())
		{
			isOppositeSide = true;  // �ҵ����������ĵ�
		}
	}

	// ����P2�������
	std::vector<std::pair<double, gp_Pnt>> distancesP2;
	for (const auto& point : oppsiteInterpolatePoints)
	{
		distancesP2.push_back({ P2.Distance(point), point });
	}

	std::sort(distancesP2.begin(), distancesP2.end(), [](std::pair<double, gp_Pnt> p1, std::pair<double, gp_Pnt> p2) {return p1.first < p2.first; });

	gp_Pnt nearestPoint1_P2 = distancesP2[0].second;
	gp_Pnt nearestPoint2_P2 = distancesP2[1].second;

	// ������������
	gp_Vec lineVec_P2(nearestPoint1_P2, nearestPoint2_P2);
	gp_Vec vecToP2(nearestPoint1_P2, P2);

	// ����ͶӰ
	double projectionLength_P2 = vecToP2.Dot(lineVec_P2) / lineVec_P2.Magnitude();  // ��һ��ͶӰ���ȣ���ֹ�����߶εĳ���

	// �ж�ͶӰ�Ƿ������������ڣ�0 <= projectionLength_P2 <= lineVec_P2.Magnitude()��
	bool isOppositeSide_P2 = false;
	if (projectionLength_P2 >= 0 && projectionLength_P2 <= lineVec_P2.Magnitude())
	{
		isOppositeSide_P2 = true;  // ���ͶӰ�����������ڣ���Ϊ����������
	}

	// ���������ͬһ�࣬����������һ����
	Standard_Integer nextP2 = 2;
	while (!isOppositeSide_P2 && nextP2 < distancesP2.size())
	{
		nearestPoint2_P2 = distancesP2[nextP2++].second;
		lineVec_P2 = gp_Vec(nearestPoint1_P2, nearestPoint2_P2);
		vecToP2 = gp_Vec(nearestPoint1_P2, P2);
		projectionLength_P2 = vecToP2.Dot(lineVec_P2) / lineVec_P2.Magnitude();

		if (projectionLength_P2 >= 0 && projectionLength_P2 <= lineVec_P2.Magnitude())
		{
			isOppositeSide_P2 = true;  // �ҵ����������ĵ�
		}
	}

	// ����뾶
	double L1 = P1.Distance(nearestPoint1_P1) + P1.Distance(nearestPoint2_P1);
	double L2 = P2.Distance(nearestPoint1_P2) + P2.Distance(nearestPoint2_P2);
	double searchRadius = (std::max(L1,L2) - std::min(L1, L2)) / 2;

	Standard_Integer M = 0;
	double w1, w2;
	if (L1 > L2)
	{
		for (const auto& point : oppsiteInterpolatePoints)
		{
			if (nearestPoint1_P2.Distance(point) <= searchRadius 
				|| nearestPoint2_P2.Distance(point) <= searchRadius)
				M++;
		}

		w2 = (double)M / (M + 2);
		w1 = 1 - w2;
		return std::make_pair(w1, w2);
	}
	else
	{
		for (const auto& point : isoInterpolatePoints)
		{
			if (nearestPoint1_P1.Distance(point) <= searchRadius 
				|| nearestPoint2_P1.Distance(point) <= searchRadius)
				M++;
		}

		w1 = (double)M / (M + 2);
		w2 = 1 - w1;
		return std::make_pair(w1, w2);
	}
}

void ProcessISOCurvesWithTangent(
	const std::vector<Handle(Geom_BSplineCurve)>& isoCurvesArray_New,
	const std::vector<Handle(Geom_BSplineCurve)>& oppsiteISOcurvesArray_New,
	std::vector<std::vector<gp_Pnt>>& isoInterpolatePoints,
	std::vector<std::vector<gp_Pnt>>& oppsiteInterpolatePoints,
	std::vector<Handle(Geom_BSplineCurve)>& isoCurvesArray_Final,
	std::vector<std::vector<double>>& knotsArray,
	std::vector<gp_Pnt>& boundaryPoints,
	std::vector<gp_Pnt>& interPoints,
	Standard_Integer isoCount,
	std::vector<TopoDS_Edge>& TangentArray,
	const std::vector<Handle(Geom_BSplineSurface)>& surfaceArr)
{
	Standard_Integer degree = isoCurvesArray_New[0]->Degree();
	for (Standard_Integer i = 0; i < isoCurvesArray_New.size(); i++)
	{
		auto curve = isoCurvesArray_New[i];

		gp_Pnt startPoint = curve->StartPoint();
		gp_Pnt endPoint = curve->EndPoint();
		boundaryPoints.push_back(startPoint);
		boundaryPoints.push_back(endPoint);
		std::vector<gp_Pnt> intersectionPoints;

		// ������Է���ĵȲ��ߣ����㽻��
		for (Standard_Integer j = 0; j < oppsiteISOcurvesArray_New.size(); j++)
		{
			Handle(Geom_BSplineCurve) oppositeCurve = oppsiteISOcurvesArray_New[j];
			GeomAPI_ExtremaCurveCurve extrema(curve, oppositeCurve);
			if (extrema.NbExtrema() > 0)
			{
				gp_Pnt P1, P2;
				extrema.NearestPoints(P1, P2);
				std::pair<double, double> weights = processPoints(P1, P2, isoInterpolatePoints[i], oppsiteInterpolatePoints[j]);
				// ��ȡȨ��ֵ
				double w1 = weights.first;
				double w2 = weights.second;
				gp_Pnt midPoint;
				double x = (P1.X() * w1 + P2.X() * w2);
				double y = (P1.Y() * w1 + P2.Y() * w2);
				double z = (P1.Z() * w1 + P2.Z() * w2);
				midPoint.SetX(x); midPoint.SetY(y);midPoint.SetZ(z);
				intersectionPoints.push_back(midPoint);
				interPoints.push_back(midPoint);
			}
		}

		interPoints.insert(interPoints.end(), isoInterpolatePoints[i].begin() + 1, isoInterpolatePoints[i].end() - 1);
		std::sort(interPoints.begin(), interPoints.end(),
			[&startPoint](const gp_Pnt& p1, const gp_Pnt& p2)
			{
				return p1.Distance(startPoint) < p2.Distance(startPoint);
			});
		std::sort(intersectionPoints.begin(), intersectionPoints.end(),
			[&startPoint](const gp_Pnt& p1, const gp_Pnt& p2)
			{
				return p1.Distance(startPoint) < p2.Distance(startPoint);
			});

		intersectionPoints.insert(intersectionPoints.begin(), startPoint);
		intersectionPoints.push_back(endPoint);

		for (Standard_Integer j = 1; j < intersectionPoints.size() - 1; j++)
		{
			if (intersectionPoints[j].Distance(intersectionPoints[j - 1]) < startPoint.Distance(endPoint) / (isoCount * 2) ||
				intersectionPoints[j].Distance(intersectionPoints[j + 1]) < startPoint.Distance(endPoint) / (isoCount * 2))
			{
				intersectionPoints.erase(intersectionPoints.begin() + j);
				--j;
			}
		}

		// ��������Լ��
	// U_Tangentָ�ڸõ��U������V_Tangentָ�ڸõ��V������
		gp_Vec U_StartTangent, V_StartTangent;
		gp_Vec U_EndTangent, V_EndTangent;
		// ���������ָ������ĽǶȣ�ȡС��
		double cosAngleStartU, cosAngleStartV, cosAngleEndU, cosAngleEndV;
		gp_Vec startDirection = (intersectionPoints[1].XYZ() - intersectionPoints[0].XYZ()).Normalized(); // ����ָ��
		gp_Vec endDirection = (intersectionPoints[intersectionPoints.size() - 1].XYZ() - intersectionPoints[intersectionPoints.size() - 2].XYZ()).Normalized(); // ����ָ��

		// startPoint���ߺ�endPoint������
		gp_Vec FirstD1, LastD1;
		if (!surfaceArr[0].IsNull())
		{

			Standard_Real uMin, uMax, vMin, vMax;
			surfaceArr[0]->Bounds(uMin, uMax, vMin, vMax); // ��ȡ������Χ
			GeomAPI_ProjectPointOnSurf projector(intersectionPoints.front(), surfaceArr[0]);
			if (projector.NbPoints() > 0)
			{
				gp_Pnt closestPoint = projector.NearestPoint();
				double uParam, vParam;
				projector.LowerDistanceParameters(uParam, vParam);

				// ��ȡU��V�����ƫ��������������
				surfaceArr[0]->D1(uParam, vParam, closestPoint, U_StartTangent, V_StartTangent);
				// ����ƫ��Ϊ������
				if (!U_StartTangent.Magnitude() || !V_StartTangent.Magnitude())
				{
					// ȷ�� uParam ����Ч������
					if (uParam - Precision::Confusion() < uMin)
					{
						uParam += Precision::Confusion();  // ��������
					}
					if (uParam + Precision::Confusion() > uMax)
					{
						uParam -= Precision::Confusion();  // ��С�������
					}
					else
					{
						uParam += Precision::Confusion();  // ��������
					}

					// ȷ�� vParam ����Ч������
					if (vParam - Precision::Confusion() < vMin)
					{
						vParam += Precision::Confusion();  // ��������
					}
					if (vParam + Precision::Confusion() > vMax)
					{
						vParam -= Precision::Confusion();  // ��С�������
					}
					else
					{
						vParam += Precision::Confusion();  // ��������
					}
					surfaceArr[0]->D1(uParam, vParam, closestPoint, U_StartTangent, V_StartTangent);
				}
				// ���㷨������U��Vƫ�����Ĳ����
				gp_Vec normalVector = U_StartTangent.Crossed(V_StartTangent);
				normalVector = normalVector.Normalized();
				// ȷ����������startDirection�н�Ϊ���
				if (startDirection.Dot(normalVector) < 0)
				{
					normalVector.Reverse(); // ��������
				}

				// ����startDirection�ڷ������ϵ�ͶӰ����
				double projectionMagnitude = startDirection.Dot(normalVector) / normalVector.Magnitude();
				gp_Vec projection = normalVector * projectionMagnitude; // ����ͶӰ����

				// ����startDirection��ƽ���ϵ�ͶӰ��ȥ���������ķ�����
				gp_Vec projectedDirection = startDirection - projection;

				// ��һ��ͶӰ����
				projectedDirection.Normalize();

				FirstD1 = projectedDirection; // ��ͶӰ��ķ�����Ϊ����
			}
		}
		else
		{
			// ����ǰ�����������
			if (intersectionPoints.size() >= 3)
			{
				std::vector<gp_Pnt> points(intersectionPoints.begin(), intersectionPoints.begin() + 3);
				FirstD1 = CalTangent(points, 1).Normalized();
			}
		}
		FirstD1.Multiply(CalPointsChordLen(intersectionPoints) / FirstD1.Magnitude());
		TangentArray.push_back(BRepBuilderAPI_MakeEdge(intersectionPoints.front(), intersectionPoints.front().Translated(FirstD1 * 0.1)).Edge());

		if (!surfaceArr[1].IsNull())
		{
			Standard_Real uMin, uMax, vMin, vMax;
			surfaceArr[1]->Bounds(uMin, uMax, vMin, vMax); // ��ȡ������Χ
			GeomAPI_ProjectPointOnSurf projector(intersectionPoints.back(), surfaceArr[1]);
			if (projector.NbPoints() > 0)
			{
				gp_Pnt closestPoint = projector.NearestPoint();
				double uParam, vParam;
				projector.LowerDistanceParameters(uParam, vParam);

				// ��ȡU��V�����ƫ��������������
				surfaceArr[1]->D1(uParam, vParam, closestPoint, U_EndTangent, V_EndTangent);
				if (!U_EndTangent.Magnitude() || !V_EndTangent.Magnitude())
				{
					// ȷ�� uParam ����Ч������
					if (uParam - Precision::Confusion() < uMin)
					{
						uParam += Precision::Confusion();  // ��������
					}
					if (uParam + Precision::Confusion() > uMax)
					{
						uParam -= Precision::Confusion();  // ��С�������
					}
					else
					{
						uParam += Precision::Confusion();  // ��������
					}

					// ȷ�� vParam ����Ч������
					if (vParam - Precision::Confusion() < vMin)
					{
						vParam += Precision::Confusion();  // ��������
					}
					if (vParam + Precision::Confusion() > vMax)
					{
						vParam -= Precision::Confusion();  // ��С�������
					}
					else
					{
						vParam += Precision::Confusion();  // ��������
					}
					surfaceArr[1]->D1(uParam, vParam, closestPoint, U_EndTangent, V_EndTangent);
				}
				// ���㷨������U��Vƫ�����Ĳ����
				gp_Vec normalVector = U_EndTangent.Crossed(V_EndTangent);
				normalVector = normalVector.Normalized();
				// ȷ����������startDirection�н�Ϊ���
				if (endDirection.Dot(normalVector) < 0)
				{
					normalVector.Reverse(); // ��������
				}
				// ����endDirection�ڷ������ϵ�ͶӰ����
				double projectionMagnitude = endDirection.Dot(normalVector) / normalVector.Magnitude();
				gp_Vec projection = normalVector * projectionMagnitude; // ����ͶӰ����

				// ����endDirection��ƽ���ϵ�ͶӰ��ȥ���������ķ�����
				gp_Vec projectedDirection = endDirection - projection;

				// ��һ��ͶӰ����
				projectedDirection.Normalize();

				LastD1 = projectedDirection; // ��ͶӰ��ķ�����Ϊ����
			}
		}
		else
		{
			// ����ǰ�����������
			if (intersectionPoints.size() >= 3)
			{
				std::vector<gp_Pnt> points(intersectionPoints.rbegin(), intersectionPoints.rbegin() + 3);
				LastD1 = CalTangent(points, 1).Normalized();
				LastD1.Reverse();
			}
		}


		LastD1.Multiply(CalPointsChordLen(intersectionPoints) / LastD1.Magnitude());
		TangentArray.push_back(BRepBuilderAPI_MakeEdge(intersectionPoints.back(),intersectionPoints.back().Translated(LastD1 * 0.1)).Edge());
		std::vector<double> params = ComputeUniformParam(intersectionPoints.size(), 0., 1.);
		std::vector<double> tempKnots = KnotGernerationByParams(params, 3, degree);
		std::vector<double> insertKnots;

		Handle(Geom_BSplineCurve) aBSplineCurve = IterateApproximate(insertKnots, intersectionPoints, FirstD1, LastD1, params, tempKnots, degree, 50, 0.01);
		knotsArray.push_back(GetKnotsSequence(aBSplineCurve));
		isoCurvesArray_Final.emplace_back(aBSplineCurve);

		gp_Vec TestFirstD1, TestLastD1;
		gp_Pnt pnt;
		aBSplineCurve->D1(aBSplineCurve->FirstParameter(),pnt, TestFirstD1);
		aBSplineCurve->D1(aBSplineCurve->LastParameter(),pnt, TestLastD1);
		auto areVectorsEqual = [](const gp_Vec& vec1, const gp_Vec& vec2, double tol = 1e-6)
			{
			return std::abs(vec1.X() - vec2.X()) <= tol &&
				std::abs(vec1.Y() - vec2.Y()) <= tol &&
				std::abs(vec1.Z() - vec2.Z()) <= tol;
			};
		if (!areVectorsEqual(TestFirstD1, FirstD1))
		{
			std::cout << "TestFirstD1 �� FirstD1 �����" << std::endl;
		}
		if (!areVectorsEqual(TestLastD1, LastD1))
		{
			std::cout << "TestLastD1 �� LastD1 �����" << std::endl;
		}


	}
}

// ����������������
Handle(Geom_BSplineSurface) FindClosestSurface(
	gp_Pnt point,
	std::vector<Handle(Geom_BSplineSurface)>& surfaceArr,
	double threshold)
{
	double minDistance = std::numeric_limits<double>::max();
	Standard_Integer closestSurfaceIndex = -1;

	for (size_t i = 0; i < surfaceArr.size(); ++i)
	{
		const auto& surface = surfaceArr[i];
		GeomAPI_ProjectPointOnSurf projector(point, surface);
		double distance = 0.0;

		if (projector.NbPoints() > 0)
		{
			gp_Pnt closestPoint = projector.NearestPoint();
			distance = point.Distance(closestPoint);
		}

		if (distance < minDistance)
		{
			minDistance = distance;
			closestSurfaceIndex = static_cast<Standard_Integer>(i);
		}
	}

	if (closestSurfaceIndex != -1 && minDistance < threshold)
	{
		return surfaceArr[closestSurfaceIndex];
	}
	return nullptr; // ��û���ҵ���Ч�ı��棬���ؿ�ָ��
}

void SurfaceModelingTool::CreateFinalISOCurves(
	const std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_New,
	const std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_New,
	std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Final,
	std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Final,
	std::vector<std::vector<gp_Pnt>>& uInterpolatePoints,
	std::vector<std::vector<gp_Pnt>>& vInterpolatePoints,
	std::vector<std::vector<double>>& uKnots,
	std::vector<std::vector<double>>& vKnots,
	std::vector<gp_Pnt>& boundaryPoints,
	std::vector<gp_Pnt>& interPoints,
	Standard_Integer isoCount,
	std::vector<TopoDS_Edge>& TangentArray,
	std::vector<Handle(Geom_BSplineSurface)>& surfaceArr)
{
	auto uCurve = uISOcurvesArray_New[0];
	Standard_Integer degree = uCurve->Degree();
	UniformCurve(uCurve);
	uKnots.push_back(GetKnotsSequence(uCurve));

	auto vCurve = vISOcurvesArray_New[0];
	UniformCurve(vCurve);
	vKnots.push_back(GetKnotsSequence(vCurve));

	// Surfacce[0] ���� startpoint���ڵ��������ڵ���
	// Surfacce[1] ���� endpoint ���ڵ��������ڵ���
	std::vector<Handle(Geom_BSplineSurface)> uTangentSurface(2);
	std::vector<Handle(Geom_BSplineSurface)> vTangentSurface(2);

	double minDistance = INT_MAX;  // ��ʼΪ���ֵ
	Standard_Integer closestSurfaceIndex = -1;  // ��ʼ��Ϊ��Ч����

	gp_Pnt startPoint = uISOcurvesArray_New[uISOcurvesArray_New.size() / 2]->StartPoint();
	gp_Pnt endPoint = uISOcurvesArray_New[uISOcurvesArray_New.size() / 2]->EndPoint();

	if (startPoint.Distance(endPoint) > 1e-6)
	{
		// ���� uTangentSurface[0] �� uTangentSurface[1]
		uTangentSurface[0] = FindClosestSurface(startPoint, surfaceArr, startPoint.Distance(endPoint) / 1000.0);
		uTangentSurface[1] = FindClosestSurface(endPoint, surfaceArr, startPoint.Distance(endPoint) / 1000.0);
	}


	// ʹ�� vISOcurvesArray_New ���м���������ȡ�µ������յ�
	gp_Pnt newStartPoint = vISOcurvesArray_New[vISOcurvesArray_New.size() / 2]->StartPoint();
	gp_Pnt newEndPoint = vISOcurvesArray_New[vISOcurvesArray_New.size() / 2]->EndPoint();

	if (newStartPoint.Distance(newEndPoint) > 1e-6)
	{
		// ���� vTangentSurface[0] �� vTangentSurface[1]
		vTangentSurface[0] = FindClosestSurface(newStartPoint, surfaceArr, newStartPoint.Distance(newEndPoint) / 1000);
		vTangentSurface[1] = FindClosestSurface(newEndPoint, surfaceArr, newStartPoint.Distance(newEndPoint) / 1000);
	}
	ProcessISOCurvesWithTangent(uISOcurvesArray_New, vISOcurvesArray_New,
		uInterpolatePoints, vInterpolatePoints,
		uISOcurvesArray_Final, 
		uKnots, boundaryPoints,
		interPoints, isoCount,
		TangentArray,
		uTangentSurface);
	ProcessISOCurvesWithTangent(vISOcurvesArray_New, uISOcurvesArray_New, 
		vInterpolatePoints, uInterpolatePoints,
		vISOcurvesArray_Final, 
		vKnots, boundaryPoints, 
		interPoints, isoCount, 
		TangentArray,
		vTangentSurface);

	uCurve = uISOcurvesArray_New.back();
	UniformCurve(uCurve);
	uKnots.push_back(GetKnotsSequence(uCurve));

	vCurve = vISOcurvesArray_New.back();
	UniformCurve(vCurve);
	vKnots.push_back(GetKnotsSequence(vCurve));
}
void SurfaceModelingTool::LoadBSplineCurves(const std::string& filePath, std::vector<Handle(Geom_BSplineCurve)>& curveArray)
{
	// ��ȡ�ļ���׺
	std::string extension = filePath.substr(filePath.find_last_of('.') + 1);

	TopoDS_Shape boundary;
	if (extension == "brep") 
	{
		// ��ʼ���߽�Shape
		BRep_Builder B1;
		// ���ļ���ȡBRep����
		BRepTools::Read(boundary, filePath.c_str(), B1);
	}
	else if(extension == "step" || extension == "stp")
	{
		// ���� STEP �ļ���ȡ��
		STEPControl_Reader reader;
		IFSelect_ReturnStatus status = reader.ReadFile(filePath.c_str());

		if (status == IFSelect_ReturnStatus::IFSelect_RetDone) 
		{
			// �����ȡ������
			reader.TransferRoots();
			boundary = reader.OneShape();
		}
	}
	else if (extension == "igs" || extension == "iges")
	{
		// ���� IGES �ļ���ȡ��
		IGESControl_Reader igesReader;
		IFSelect_ReturnStatus status = igesReader.ReadFile(filePath.c_str());

		if (status == IFSelect_ReturnStatus::IFSelect_RetDone)
		{
			// �����ȡ������
			igesReader.TransferRoots();
			boundary = igesReader.OneShape();
		}
	}

	// ����Shape�еı�
	TopExp_Explorer explorer(boundary, TopAbs_EDGE);
	for (; explorer.More(); explorer.Next()) 
	{
		TopoDS_Edge edge = TopoDS::Edge(explorer.Current());

		// ��ȡ�ߵļ��α�ʾ
		TopLoc_Location loc;
		Standard_Real first, last;
		Handle(Geom_Curve) gcurve = BRep_Tool::Curve(edge, loc, first, last);
		gcurve = Handle(Geom_Curve)::DownCast(gcurve->Copy());

		// �����������
		if (gcurve->DynamicType() == STANDARD_TYPE(Geom_Line)) 
		{
			// �����ֱ�ߣ�ת��ΪBSpline
			Handle(Geom_TrimmedCurve) aTrimmedLine = new Geom_TrimmedCurve(gcurve, first, last);
			Handle(Geom_BSplineCurve) aGeom_BSplineCurve = GeomConvert::CurveToBSplineCurve(aTrimmedLine);
			if (!aGeom_BSplineCurve.IsNull() && aGeom_BSplineCurve->IsKind(STANDARD_TYPE(Geom_BSplineCurve)))
			{
				curveArray.push_back(aGeom_BSplineCurve);
			}
		}
		else if (gcurve->DynamicType() == STANDARD_TYPE(Geom_BSplineCurve)) 
		{
			// ����Ѿ���BSpline��ֱ�Ӵ���
			Handle(Geom_BSplineCurve) aGeom_BSplineCurve = Handle(Geom_BSplineCurve)::DownCast(gcurve);
			if (!aGeom_BSplineCurve.IsNull()) {
				aGeom_BSplineCurve->Segment(first, last);
				if (!aGeom_BSplineCurve.IsNull() && aGeom_BSplineCurve->IsKind(STANDARD_TYPE(Geom_BSplineCurve))) 
				{
					curveArray.push_back(aGeom_BSplineCurve);
				}
			}
		}
	}
}
void SurfaceModelingTool::LoadBSplineSurfaces(const std::string& filePath, std::vector<Handle(Geom_BSplineSurface)>& surfaceArray)
{
	if (!std::filesystem::exists(filePath))
	{
		return;
	}
	try
	{
		// ��ȡ�ļ���չ��
		std::filesystem::path path(filePath);
		std::string extension = path.extension().string();
		std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

		TopoDS_Shape shape;

		// �����ļ�����ѡ����Ӧ�Ķ�ȡ��
		if (extension == ".stp" || extension == ".step")
		{
			STEPControl_Reader reader;
			IFSelect_ReturnStatus status = reader.ReadFile(filePath.c_str());

			if (status != IFSelect_RetDone)
			{
				throw Standard_Failure("Failed to read STEP file");
			}

			reader.TransferRoots();
			shape = reader.OneShape();
		}
		else if (extension == ".igs" || extension == ".iges")
		{
			IGESControl_Reader reader;
			IFSelect_ReturnStatus status = reader.ReadFile(filePath.c_str());

			if (status != IFSelect_RetDone)
			{
				throw Standard_Failure("Failed to read IGES file");
			}

			reader.TransferRoots();
			shape = reader.OneShape();
		}
		else if (extension == ".brep" || extension == ".brp")
		{
			BRep_Builder builder;
			Standard_Boolean result = BRepTools::Read(shape, filePath.c_str(), builder);

			if (!result)
			{
				throw Standard_Failure("Failed to read BREP file");
			}
		}
		else
		{
			throw Standard_Failure("Unsupported file format");
		}

		if (shape.IsNull())
		{
			throw Standard_Failure("No valid shape found in file");
		}

		// ���������沢��ȡB��������
		TopExp_Explorer explorer(shape, TopAbs_FACE);
		while (explorer.More())
		{
			const TopoDS_Face& face = TopoDS::Face(explorer.Current());
			Handle(Geom_Surface) surface = BRep_Tool::Surface(face);

			if (!surface.IsNull())
			{
				Handle(Geom_BSplineSurface) bsplineSurface;

				// ����Ѿ���BSplineSurface
				if (surface->IsKind(STANDARD_TYPE(Geom_BSplineSurface)))
				{
					bsplineSurface = Handle(Geom_BSplineSurface)::DownCast(surface);
				}
				// ����ת��ΪBSplineSurface
				else
				{
					try
					{
						bsplineSurface = GeomConvert::SurfaceToBSplineSurface(surface);
					}
					catch (const Standard_Failure& e)
					{
						std::cerr << "Warning: Failed to convert surface to BSpline: "
							<< e.GetMessageString() << std::endl;
						continue;
					}
				}

				if (!bsplineSurface.IsNull())
				{
					surfaceArray.push_back(bsplineSurface);
				}
			}
			explorer.Next();
		}

		if (surfaceArray.empty())
		{
			std::cout << "Warning: No B-Spline surfaces found in file." << std::endl;
		}
		else
		{
			std::cout << "Successfully loaded " << surfaceArray.size()
				<< " B-Spline surface(s) from " << path.filename().string() << std::endl;
		}
	}
	catch (const Standard_Failure& e)
	{
		std::cerr << "Error: " << e.GetMessageString() << std::endl;
		throw;
	}
}
void SurfaceModelingTool::GetISOCurveWithNormal(const Handle(Geom_BSplineSurface)& surfacecoons, 
	std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Initial, 
	std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Initial, std::vector<gp_Vec>& normalsOfUISOLines, std::vector<gp_Vec>& normalsOfVISOLines, Standard_Integer numIsoCurves)
{
	auto IsPointOnSurface = [](const gp_Pnt& point, const Handle(Geom_Surface)& surface)
	{
		GeomAPI_ProjectPointOnSurf proj(point, surface);

		if (proj.NbPoints() > 0) 
		{
			gp_Pnt closestPoint = proj.NearestPoint();

			return point.Distance(closestPoint) < 1e-6;
		}

		return false;
	};

	auto GetClosestPointOnSurface = [](const gp_Pnt& point, const Handle(Geom_Surface)& surface)
	{
		GeomAPI_ProjectPointOnSurf proj(point, surface);

		if (proj.NbPoints() > 0)
		{
			return proj.NearestPoint(); // ���������
		}

		return point;
	};

	std::vector<std::pair<gp_Vec, gp_Vec>> tangentOfUISOLines;
	std::vector<std::pair<gp_Vec, gp_Vec>> tangentOfVISOLines;
	const Standard_Integer numSamplePoints = 10;
	for (Standard_Integer i = 1; i < numIsoCurves; i++)
	{
		std::vector<gp_Vec> normalsU;
		std::vector<gp_Vec> normalsV;

		Handle(Geom_BSplineCurve) aUGeom_BSplineCurve = Handle(Geom_BSplineCurve)::DownCast(surfacecoons->UIso(((double)i / numIsoCurves) 
			* (surfacecoons->UKnot(surfacecoons->LastUKnotIndex()) - surfacecoons->UKnot(surfacecoons->FirstUKnotIndex())) + surfacecoons->UKnot(surfacecoons->FirstUKnotIndex())));
		uISOcurvesArray_Initial.emplace_back(aUGeom_BSplineCurve);

		// ��ȡu����BSpline���ߵĲ�����Χ
		double uStart = aUGeom_BSplineCurve->FirstParameter();
		double uEnd = aUGeom_BSplineCurve->LastParameter();

		gp_Vec startTangent, endTangent;
		for (Standard_Integer j = 0; j < numSamplePoints; j++)
		{
			double t = uStart + j * (uEnd - uStart) / (numSamplePoints - 1); // ��ʵ�ʲ�����Χ�ھ���ȡ��
			gp_Pnt p1 = aUGeom_BSplineCurve->Value(t);
			gp_Vec DXu, DXv, N;

			if (IsPointOnSurface(p1, surfacecoons)) 
			{
				surfacecoons->D1(((double)i / numIsoCurves) * (surfacecoons->UKnot(surfacecoons->LastUKnotIndex()) - surfacecoons->UKnot(surfacecoons->FirstUKnotIndex())) + surfacecoons->UKnot(surfacecoons->FirstUKnotIndex()), 
					((double)j / numSamplePoints) * (surfacecoons->VKnot(surfacecoons->LastVKnotIndex()) - surfacecoons->VKnot(surfacecoons->FirstVKnotIndex())) + surfacecoons->VKnot(surfacecoons->FirstVKnotIndex()), p1, DXu, DXv);
			}
			else
			{
				gp_Pnt closestPoint = GetClosestPointOnSurface(p1, surfacecoons);
				surfacecoons->D1(((double)i / numIsoCurves) * (surfacecoons->UKnot(surfacecoons->LastUKnotIndex()) - surfacecoons->UKnot(surfacecoons->FirstUKnotIndex())) + surfacecoons->UKnot(surfacecoons->FirstUKnotIndex()),
					((double)j / numSamplePoints) * (surfacecoons->VKnot(surfacecoons->LastVKnotIndex()) - surfacecoons->VKnot(surfacecoons->FirstVKnotIndex())) + surfacecoons->VKnot(surfacecoons->FirstVKnotIndex()), closestPoint, DXu, DXv);
			}

			DXu.Cross(DXv);
			N = DXu.Normalized();
			normalsU.push_back(N);
		}

		// ƽ��������
		gp_Vec avgNormalU = std::accumulate(normalsU.begin(), normalsU.end(), gp_Vec()) / normalsU.size();
		normalsOfUISOLines.emplace_back(avgNormalU);

		// ����v�����Iso����
		Handle(Geom_BSplineCurve) aVGeom_BSplineCurve = Handle(Geom_BSplineCurve)::DownCast(surfacecoons->VIso(((double)i / numIsoCurves)
			* (surfacecoons->VKnot(surfacecoons->LastVKnotIndex()) - surfacecoons->VKnot(surfacecoons->FirstVKnotIndex())) + surfacecoons->VKnot(surfacecoons->FirstVKnotIndex())));
		vISOcurvesArray_Initial.emplace_back(aVGeom_BSplineCurve);

		double vStart = aVGeom_BSplineCurve->FirstParameter();
		double vEnd = aVGeom_BSplineCurve->LastParameter();
		for (Standard_Integer j = 0; j < numSamplePoints; j++)
		{
			double t = vStart + j * (vEnd - vStart) / (numSamplePoints - 1);
			gp_Pnt p1 = aVGeom_BSplineCurve->Value(t);
			gp_Vec DXu, DXv, N;
			gp_Vec startTangent, endTangent;

			// ���㷨����
			if (IsPointOnSurface(p1, surfacecoons))
			{
				surfacecoons->D1(((double)j / numSamplePoints) * (surfacecoons->UKnot(surfacecoons->LastUKnotIndex()) - surfacecoons->UKnot(surfacecoons->FirstUKnotIndex())) + surfacecoons->UKnot(surfacecoons->FirstUKnotIndex()), 
					((double)i / numIsoCurves) * (surfacecoons->VKnot(surfacecoons->LastVKnotIndex()) - surfacecoons->VKnot(surfacecoons->FirstVKnotIndex())) + surfacecoons->VKnot(surfacecoons->FirstVKnotIndex()), p1, DXu, DXv);
			}
			else
			{
				gp_Pnt closestPoint = GetClosestPointOnSurface(p1, surfacecoons);
				surfacecoons->D1(((double)j / numSamplePoints) * (surfacecoons->UKnot(surfacecoons->LastUKnotIndex()) - surfacecoons->UKnot(surfacecoons->FirstUKnotIndex())) + surfacecoons->UKnot(surfacecoons->FirstUKnotIndex()), 
					((double)i / numIsoCurves) * (surfacecoons->VKnot(surfacecoons->LastVKnotIndex()) - surfacecoons->VKnot(surfacecoons->FirstVKnotIndex())) + surfacecoons->VKnot(surfacecoons->FirstVKnotIndex()), closestPoint, DXu, DXv);
			}

			//gp_Vec direction = aVGeom_BSplineCurve->EndPoint().XYZ() - aVGeom_BSplineCurve->StartPoint().XYZ(); // �õ���������
			//direction.Normalize();
			//if (j == 0)
			//{
			//	startTangent = DXv.Normalized();
			//	if (startTangent.Dot(direction) < 0)
			//	{
			//		startTangent = -startTangent; // ��ת�������ķ���
			//	}
			//}
			//if (j == numSamplePoints - 1)
			//{
			//	endTangent = DXv.Normalized();
			//	if (endTangent.Dot(direction) < 0)
			//	{
			//		endTangent = -endTangent; // ��ת�������ķ���
			//	}
			//}
			//tangentOfVISOLines.push_back(std::make_pair(startTangent.Normalized(), endTangent.Normalized()));

			DXu.Cross(DXv);
			N = DXu.Normalized();
			normalsV.push_back(N);
		}

		// ƽ��������
		gp_Vec avgNormalV = std::accumulate(normalsV.begin(), normalsV.end(), gp_Vec()) / normalsV.size();
		normalsOfVISOLines.emplace_back(avgNormalV);

	}
}


bool SurfaceModelingTool::ExportBSplineCurves(const std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_Final,
	const std::string& Filename)
{
	// ���� TopoDS_Compound ����
	TopoDS_Compound Result;
	BRep_Builder builder;
	builder.MakeCompound(Result);

	// �����������飬ת��Ϊ TopoDS_Edge ����ӵ� Result
	for (const auto& curve : ISOcurvesArray_Final)
	{
		if (curve.IsNull())
		{
			std::cerr << "���棺���ֿյ� BSpline ���ߣ�������" << std::endl;
			continue;
		}

		// ���� TopoDS_Edge
		TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);

		// �����Ƿ���Ч
		if (edge.IsNull())
		{
			std::cerr << "���棺�޷������ߣ�������������" << std::endl;
			continue;
		}

		// ������ӵ�������
		builder.Add(Result, edge);
	}

	// ����Ƿ�����Ч�ı���ӵ�������
	if (Result.IsNull())
	{
		std::cerr << "����û����Ч�����߱���ӵ������塣" << std::endl;
		return false;
	}

	// ���������嵽�ļ�
	if (!BRepTools::Write(Result, Filename.c_str()))
	{
		std::cerr << "�����޷������ߵ������ļ� " << Filename << std::endl;
		return false;
	}
	else
	{
		std::cout << "�ɹ��������ѵ������ļ� " << Filename << std::endl;
	}

	return true;
}
void SurfaceModelingTool::ApproximateBoundaryCurves(std::vector<Handle(Geom_BSplineCurve)>& curves, Standard_Integer samplingNum)
{
	for (auto& curve : curves) 
	{
		TColStd_Array1OfReal curveKnots(1, curve->NbKnots());
		curve->Knots(curveKnots);

		// ���²��������ߵĽڵ�
		if (!(curveKnots(curveKnots.Lower()) == 0 && curveKnots(curveKnots.Upper()) == 1))
		{
			BSplCLib::Reparametrize(0, 1, curveKnots);
			curve->SetKnots(curveKnots);
		}

		// ���������������
		std::vector<gp_Pnt> samplingPnts;
		std::vector<Standard_Real> samplingParams;
		Standard_Real vMin = curve->FirstParameter();
		Standard_Real vMax = curve->LastParameter();

		for (Standard_Integer j = 1; j <= samplingNum; j++) 
		{
			Standard_Real param = vMin + (vMax - vMin) * (j - 1) / (samplingNum - 1);
			gp_Pnt pnt = curve->Value(param);
			samplingParams.push_back(param);
			samplingPnts.push_back(pnt);
		}

		// ��ʼ���ڵ㲢�������
		std::vector<Standard_Real> init_knots = KnotGernerationByParams(samplingParams, 6, 3);
		std::vector<Standard_Real> insertKnots;
		curve = IterateApproximate(insertKnots, samplingPnts, samplingParams, init_knots, 3, 10, 1);
	}
}

void SurfaceModelingTool::UpdateFinalCurves(const std::vector<Handle(Geom_BSplineCurve)>& aBoundarycurveArray, 
	std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Final, 
	std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Final)
{
	Handle(Geom_BSplineCurve) uCurve = uISOcurvesArray_Final[0];
	Handle(Geom_BSplineCurve) vCurve = vISOcurvesArray_Final[0];
	std::vector<Handle(Geom_BSplineCurve)> uBoundaryCurve(2), vBoundaryCurve(2);

	bool isThreeBoundary = false;
	for (auto theBoundaryCurve : aBoundarycurveArray)
	{
		double distance = theBoundaryCurve->StartPoint().Distance(theBoundaryCurve->EndPoint());
		if (distance < 0.1)
		{
			// ����
			isThreeBoundary = true;
			break;
		}
	}
	if (isThreeBoundary)
	{
		GeomAPI_ExtremaCurveCurve extrema1(uCurve, aBoundarycurveArray[0]);
		GeomAPI_ExtremaCurveCurve extrema2(uCurve, aBoundarycurveArray[2]);

		gp_Pnt p1, p2, p3;
		extrema1.NearestPoints(p1,p2);
		extrema2.NearestPoints(p1,p3);

		double distance = p2.Distance(p3);
		if (distance < Precision::Confusion())
		{
			uBoundaryCurve[0] = aBoundarycurveArray[0];
			uBoundaryCurve[1] = aBoundarycurveArray[2];
			vBoundaryCurve[0] = aBoundarycurveArray[1];
			vBoundaryCurve[1] = aBoundarycurveArray[3];
		}
		else
		{
			uBoundaryCurve[0] = aBoundarycurveArray[1];
			uBoundaryCurve[1] = aBoundarycurveArray[3];
			vBoundaryCurve[0] = aBoundarycurveArray[0];
			vBoundaryCurve[1] = aBoundarycurveArray[2];
		}
	}
	else
	{
		bool isUBoundaryFirst = MathTool::ComputeCurveCurveDistance(aBoundarycurveArray[0], uCurve) >
			MathTool::ComputeCurveCurveDistance(aBoundarycurveArray[1], uCurve);

		uBoundaryCurve[0] = isUBoundaryFirst ? aBoundarycurveArray[0] : aBoundarycurveArray[1];
		uBoundaryCurve[1] = isUBoundaryFirst ? aBoundarycurveArray[2] : aBoundarycurveArray[3];
		vBoundaryCurve[0] = isUBoundaryFirst ? aBoundarycurveArray[1] : aBoundarycurveArray[0];
		vBoundaryCurve[1] = isUBoundaryFirst ? aBoundarycurveArray[3] : aBoundarycurveArray[2];
	}

	gp_Pnt uCurveSamplePoint = MathTool::ComputeAverageSamplePoint(uCurve, 10);
	bool isUFirstCloser = uCurveSamplePoint.Distance(MathTool::ComputeAverageSamplePoint(uBoundaryCurve[0], 10)) <
		uCurveSamplePoint.Distance(MathTool::ComputeAverageSamplePoint(uBoundaryCurve[1], 10));

	uISOcurvesArray_Final.insert(uISOcurvesArray_Final.begin(), isUFirstCloser ? uBoundaryCurve[0] : uBoundaryCurve[1]);
	uISOcurvesArray_Final.insert(uISOcurvesArray_Final.end(), isUFirstCloser ? uBoundaryCurve[1] : uBoundaryCurve[0]);

	gp_Pnt vCurveSamplePoint = MathTool::ComputeAverageSamplePoint(vCurve, 10);
	bool isVFirstCloser = vCurveSamplePoint.Distance(MathTool::ComputeAverageSamplePoint(vBoundaryCurve[0], 10)) <
		vCurveSamplePoint.Distance(MathTool::ComputeAverageSamplePoint(vBoundaryCurve[1], 10));

	vISOcurvesArray_Final.insert(vISOcurvesArray_Final.begin(), isVFirstCloser ? vBoundaryCurve[0] : vBoundaryCurve[1]);
	vISOcurvesArray_Final.insert(vISOcurvesArray_Final.end(), isVFirstCloser ? vBoundaryCurve[1] : vBoundaryCurve[0]);
	gp_Pnt origin(0, 0, 0);
	auto minDistanceToPoint = [&origin](const Handle(Geom_BSplineCurve)& curve)
		{
			return std::min(curve->StartPoint().Distance(origin), curve->EndPoint().Distance(origin));
		};

	// ������һ�����߱ȵ�һ�����߾���ԭ���������ת uISOcurvesArray_Final
	if (minDistanceToPoint(uISOcurvesArray_Final[uISOcurvesArray_Final.size() - 1]) < minDistanceToPoint(uISOcurvesArray_Final[0]))
	{
		std::reverse(uISOcurvesArray_Final.begin(), uISOcurvesArray_Final.end());
	}
	// �� vISO ���߽���ͬ���Ĵ���
	if (minDistanceToPoint(vISOcurvesArray_Final[vISOcurvesArray_Final.size() - 1]) < minDistanceToPoint(vISOcurvesArray_Final[0]))
	{
		std::reverse(vISOcurvesArray_Final.begin(), vISOcurvesArray_Final.end());
	}
	MathTool::ReverseIfNeeded(uISOcurvesArray_Final);
	MathTool::ReverseIfNeeded(vISOcurvesArray_Final);
}

Standard_Boolean SurfaceModelingTool::GetInternalCurves(
	std::vector<Handle(Geom_BSplineCurve)>& theBoundaryCurveArray,
	std::vector<Handle(Geom_BSplineCurve)>& theInternalBSplineCurves,
	std::vector<Handle(Geom_BSplineCurve)>& theUInternalCurve,
	std::vector<Handle(Geom_BSplineCurve)>& theVInternalCurve,
	Standard_Real& theUAngleSum,
	Standard_Real& theVAngleSum,
	Standard_Real theAngleTolerance) {

	if (theBoundaryCurveArray.size() != 4)
	{
		return Standard_False;
	}

	Handle(Geom_BSplineCurve) aBoundaryCurve1 = theBoundaryCurveArray[0];
	Handle(Geom_BSplineCurve) aBoundaryCurve2 = theBoundaryCurveArray[1];
	Handle(Geom_BSplineCurve) aBoundaryCurve3 = theBoundaryCurveArray[2];
	Handle(Geom_BSplineCurve) aBoundaryCurve4 = theBoundaryCurveArray[3];

	std::vector<PlanarCurve> aPlanarCurveArray;
	for (Standard_Integer i = 0; i < theBoundaryCurveArray.size(); ++i) 
	{
		aPlanarCurveArray.emplace_back(PlanarCurve(theBoundaryCurveArray[i]));
	}

	Standard_Boolean isPlanar = Standard_True;
	for (const PlanarCurve& aCurve : aPlanarCurveArray)
	{
		if (aCurve.GetCurveType() == CurveType::NOTPLANAR) 
		{
			isPlanar = Standard_False;
			break;
		}
	}

	if (!isPlanar) {
		return Standard_False;
	}

	theUInternalCurve.clear();
	theVInternalCurve.clear();

	for (auto& anInternalCurve : theInternalBSplineCurves)
	{
		PlanarCurve anInternalPlanarCurve(anInternalCurve);

		if (anInternalPlanarCurve.GetCurveType() == CurveType::NOTPLANAR) 
		{
			continue;
		}

		Standard_Real aDistance1 = MathTool::ComputeCurveCurveDistance(anInternalPlanarCurve.GetCurve(), aBoundaryCurve1);
		Standard_Real aDistance2 = MathTool::ComputeCurveCurveDistance(anInternalPlanarCurve.GetCurve(), aBoundaryCurve2);
		Standard_Real aDistance3 = MathTool::ComputeCurveCurveDistance(anInternalPlanarCurve.GetCurve(), aBoundaryCurve3);
		Standard_Real aDistance4 = MathTool::ComputeCurveCurveDistance(anInternalPlanarCurve.GetCurve(), aBoundaryCurve4);

		Standard_Real aSplitPointParams[2] = { 0.0 };

		if ((aDistance1 < 10.0 && aDistance3 < 10.0) || (aDistance2 < 10.0 && aDistance4 < 10.0))
		{
			GeomAPI_ExtremaCurveCurve anExtrema1(anInternalPlanarCurve.GetCurve(), aDistance1 < 10.0 ? aBoundaryCurve1 : aBoundaryCurve2);
			GeomAPI_ExtremaCurveCurve anExtrema2(anInternalPlanarCurve.GetCurve(), aDistance3 < 10.0 ? aBoundaryCurve3 : aBoundaryCurve4);
			gp_Pnt anInternalPoint;
			gp_Pnt aReplacePoint1, aReplacePoint2;

			Standard_Real aParameter = 0;
			if (anExtrema1.NbExtrema() > 0) 
			{
				anExtrema1.LowerDistanceParameters(aSplitPointParams[0], aParameter);
				anExtrema1.NearestPoints(anInternalPoint, aReplacePoint1);
			}

			if (anExtrema2.NbExtrema() > 0)
			{
				anExtrema2.LowerDistanceParameters(aSplitPointParams[1], aParameter);
				anExtrema2.NearestPoints(anInternalPoint, aReplacePoint2);
			}

			if (aSplitPointParams[0] > aSplitPointParams[1])
			{
				std::swap(aSplitPointParams[0], aSplitPointParams[1]);
				std::swap(aReplacePoint1, aReplacePoint2);
			}

			Handle(Geom_TrimmedCurve) aTrimmedCurve = new Geom_TrimmedCurve(anInternalPlanarCurve.GetCurve(), aSplitPointParams[0], aSplitPointParams[1]);
			Handle(Geom_BSplineCurve) aModifiedCurve = GeomConvert::CurveToBSplineCurve(aTrimmedCurve, Convert_TgtThetaOver2);

			Standard_Real anAngle1 = MathTool::ComputeAngleBetweenPlanarCurves(aPlanarCurveArray[0], anInternalPlanarCurve);
			Standard_Real anAngle3 = MathTool::ComputeAngleBetweenPlanarCurves(aPlanarCurveArray[2], anInternalPlanarCurve);
			Standard_Real anAngle2 = MathTool::ComputeAngleBetweenPlanarCurves(aPlanarCurveArray[1], anInternalPlanarCurve);
			Standard_Real anAngle4 = MathTool::ComputeAngleBetweenPlanarCurves(aPlanarCurveArray[3], anInternalPlanarCurve);

			aModifiedCurve->SetPole(1, aReplacePoint1);
			aModifiedCurve->SetPole(aModifiedCurve->NbPoles(), aReplacePoint2);
			anInternalPlanarCurve.SetCurve(aModifiedCurve);

			if (std::abs(anAngle1) < theAngleTolerance && std::abs(anAngle3) < theAngleTolerance &&
				std::abs(anAngle2) < theAngleTolerance && std::abs(anAngle4) < theAngleTolerance)
			{
				if (aDistance1 < 10.0 && aDistance3 < 10.0)
				{
					theVAngleSum += (std::abs(anAngle2) + std::abs(anAngle4)) / 2.0;
					theVInternalCurve.push_back(anInternalPlanarCurve.GetCurve());
				}
				else if(aDistance2 < 10.0 && aDistance4 < 10.0)
				{
					theUAngleSum += (std::abs(anAngle1) + std::abs(anAngle3)) / 2.0;
					theUInternalCurve.push_back(anInternalPlanarCurve.GetCurve());
				}
			}
			else if (std::abs(anAngle1) < theAngleTolerance && std::abs(anAngle3) < theAngleTolerance)
			{
				theUAngleSum += (std::abs(anAngle1) + std::abs(anAngle3)) / 2.0;
				theUInternalCurve.push_back(anInternalPlanarCurve.GetCurve());
			}
			else if (std::abs(anAngle2) < theAngleTolerance && std::abs(anAngle4) < theAngleTolerance) 
			{
				theVAngleSum += (std::abs(anAngle2) + std::abs(anAngle4)) / 2.0;
				theVInternalCurve.push_back(anInternalPlanarCurve.GetCurve());
			}
		}
	}

	theUInternalCurve.insert(theUInternalCurve.begin(), aBoundaryCurve1);
	theUInternalCurve.push_back(aBoundaryCurve3);
	theVInternalCurve.insert(theVInternalCurve.begin(), aBoundaryCurve2);
	theVInternalCurve.push_back(aBoundaryCurve4);

	MathTool::SortBSplineCurves(theUInternalCurve, theUInternalCurve[0]);
	MathTool::SortBSplineCurves(theVInternalCurve, theVInternalCurve[0]);
	MathTool::CheckSelfIntersect(theUInternalCurve);
	MathTool::CheckSelfIntersect(theVInternalCurve);

	return theUInternalCurve.size() >= 4 || theVInternalCurve.size() >= 4;
}


Handle(Geom_BSplineSurface) SurfaceModelingTool::GenerateReferSurface(
	std::vector<Handle(Geom_BSplineCurve)> theBoundaryCurveArray,
	std::vector<Handle(Geom_BSplineCurve)>& theUInternalCurve,
	std::vector<Handle(Geom_BSplineCurve)>& theVInternalCurve,
	Standard_Real theUAngleSum,
	Standard_Real theVAngleSum,
	Standard_Integer theIsoCount,
	ReferSurfaceType theReferSurfaceType)
{
	if (theReferSurfaceType == ReferSurfaceType::GORDEN_ONE_DIRECTION_GORDEN) {
		// ��ȡ����ı߽�����
		Handle(Geom_BSplineCurve) aBsplineCurve1 = theBoundaryCurveArray[0];
		Handle(Geom_BSplineCurve) aBsplineCurve2 = theBoundaryCurveArray[1];
		Handle(Geom_BSplineCurve) aBsplineCurve3 = theBoundaryCurveArray[2];
		Handle(Geom_BSplineCurve) aBsplineCurve4 = theBoundaryCurveArray[3];

		// �洢���ɵ�Gorden�Ȳ������ߺ�ʣ������
		std::vector<Handle(Geom_BSplineCurve)> aGordenISOCurves;
		std::vector<Handle(Geom_BSplineCurve)> aRemainCurves;

		// ʹ�����㷨�ı�־
		Standard_Boolean aUseNewAlgorithm = true;

		// �ж��ڲ����ߵ�������ѡ����Gorden����ķ�ʽ
		if (theUInternalCurve.size() > theVInternalCurve.size() && theUInternalCurve.size() >= 4) {
			// ѡ��u������ڲ��ߺͱ߽�������Gorden����
			aGordenISOCurves.insert(aGordenISOCurves.end(), theUInternalCurve.begin(), theUInternalCurve.end());
			aRemainCurves.push_back(aBsplineCurve2);
			aRemainCurves.push_back(aBsplineCurve4);
		} else if (theVInternalCurve.size() > theUInternalCurve.size() && theVInternalCurve.size() >= 4) {
			// ѡ��v������ڲ��ߺͱ߽�������Gorden����
			aGordenISOCurves.insert(aGordenISOCurves.end(), theVInternalCurve.begin(), theVInternalCurve.end());
			aRemainCurves.push_back(aBsplineCurve1);
			aRemainCurves.push_back(aBsplineCurve3);
		} else if (theUInternalCurve.size() == theVInternalCurve.size() && theUInternalCurve.size() >= 4) {
			// ���u�����v������ڲ�����������ȣ����ݽǶ�֮����ѡ��
			if (theUAngleSum < theVAngleSum) {
				aGordenISOCurves.insert(aGordenISOCurves.end(), theUInternalCurve.begin(), theUInternalCurve.end());
				aRemainCurves.push_back(aBsplineCurve2);
				aRemainCurves.push_back(aBsplineCurve4);
			} else {
				aGordenISOCurves.insert(aGordenISOCurves.end(), theVInternalCurve.begin(), theVInternalCurve.end());
				aRemainCurves.push_back(aBsplineCurve1);
				aRemainCurves.push_back(aBsplineCurve3);
			}
		} else {
			// ������������㣬���˵������㷨
			aUseNewAlgorithm = false;
			return nullptr;
		}

		// �洢���ɵĵȲ������ߺͷ���
		std::vector<Handle(Geom_BSplineCurve)> aUCreateGordenCurves, aVCreateGordenCurves;
		std::vector<gp_Vec> aNormalsOfUISOLines, aNormalsOfVISOLines;

		// �������ɵĲο�����
		Handle(Geom_BSplineSurface) aReferSurface;

		if (aUseNewAlgorithm) {
			// ��������������֮��ĽǶ�
			Standard_Real aAngleUwithG = MathTool::ComputeAngleBetweenCurves(aBsplineCurve1, aGordenISOCurves[0]);
			Standard_Real aAngleVwithG = MathTool::ComputeAngleBetweenCurves(aBsplineCurve2, aGordenISOCurves[0]);

			// ���ݽǶ�ѡ������
			if (aAngleUwithG > aAngleVwithG) {
				// ���u����ĽǶȸ��󣬵���u��v���������˳��
				aVCreateGordenCurves.clear();
				aVCreateGordenCurves.insert(aVCreateGordenCurves.begin(), aGordenISOCurves.begin(), aGordenISOCurves.end());
				aUCreateGordenCurves.insert(aUCreateGordenCurves.begin(), aRemainCurves[0]);
				aUCreateGordenCurves.insert(aUCreateGordenCurves.end(), aRemainCurves[1]);
			} else {
				// ����v���������˳��
				aUCreateGordenCurves.clear();
				aUCreateGordenCurves.insert(aUCreateGordenCurves.begin(), aGordenISOCurves.begin(), aGordenISOCurves.end());
				aVCreateGordenCurves.insert(aVCreateGordenCurves.begin(), aRemainCurves[0]);
				aVCreateGordenCurves.insert(aVCreateGordenCurves.end(), aRemainCurves[1]);
			}

			// �����ɵ����߽������򲢼�齻��
			MathTool::SortBSplineCurves(aUCreateGordenCurves, aUCreateGordenCurves[0]);
			MathTool::SortBSplineCurves(aVCreateGordenCurves, aVCreateGordenCurves[0]);
			MathTool::ReverseIfNeeded(aUCreateGordenCurves);
			MathTool::ReverseIfNeeded(aVCreateGordenCurves);
			TopoDS_Face aGordenFace;

			//MathTool::ReverseIfNeeded(aUCreateGordenCurves);
			//MathTool::ReverseIfNeeded(aVCreateGordenCurves);
			// ���ó��ε� Compatible
			std::for_each(aUCreateGordenCurves.begin(), aUCreateGordenCurves.end(), UniformCurve);
			std::for_each(aVCreateGordenCurves.begin(), aVCreateGordenCurves.end(), UniformCurve);

			CurveOperate::CompatibleWithInterPoints(aVCreateGordenCurves, aUCreateGordenCurves);
			CurveOperate::CompatibleWithInterPoints(aUCreateGordenCurves, aVCreateGordenCurves);

			std::vector<gp_Pnt> upoints, vpoints;
			std::vector<Standard_Real> uparams, vparams;
			std::tie(upoints, uparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(aUCreateGordenCurves, aVCreateGordenCurves[0]);
			std::tie(vpoints, vparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(aVCreateGordenCurves, aUCreateGordenCurves[0]);
			// �������
			std::vector<std::pair<double, Handle(Geom_BSplineCurve)>> combinedv;
			for (size_t i = 0; i < vparams.size(); ++i) {
				combinedv.emplace_back(vparams[i], aVCreateGordenCurves[i]);
			}
			std::sort(combinedv.begin(), combinedv.end(), [](const auto& a, const auto& b) {
				return a.first < b.first;
			});
			for (size_t i = 0; i < combinedv.size(); ++i) {
				vparams[i] = combinedv[i].first;
				aVCreateGordenCurves[i] = combinedv[i].second;
			}

			std::vector<std::pair<double, Handle(Geom_BSplineCurve)>> combinedu;
			for (size_t i = 0; i < uparams.size(); ++i) {
				combinedu.emplace_back(uparams[i], aUCreateGordenCurves[i]);
			}
			std::sort(combinedu.begin(), combinedu.end(), [](const auto& a, const auto& b) {
				return a.first < b.first;
			});
			for (size_t i = 0; i < combinedu.size(); ++i) {
				uparams[i] = combinedu[i].first;
				aUCreateGordenCurves[i] = combinedu[i].second;
			}

			GordenSurface::BuildMyGordonSurf(aUCreateGordenCurves, aVCreateGordenCurves, uparams, vparams, aGordenFace);
			Handle(Geom_Surface) aGeomSurface = BRep_Tool::Surface(aGordenFace);
			aReferSurface = Handle(Geom_BSplineSurface)::DownCast(aGeomSurface);
		}

		// �������ɵĲο�����
		return aReferSurface;
	}

	if (theReferSurfaceType == ReferSurfaceType::GORDEN_TWO_DIRECTION_GORDEN) {
		// �����ɵ����߽������򲢼�齻��
		MathTool::SortBSplineCurves(theUInternalCurve, theUInternalCurve[0]);
		MathTool::SortBSplineCurves(theVInternalCurve, theVInternalCurve[0]);
		MathTool::ReverseIfNeeded(theUInternalCurve);
		MathTool::ReverseIfNeeded(theVInternalCurve);
		TopoDS_Face aGordenFace;

		/*MathTool::ReverseIfNeeded(theUInternalCurve);
		MathTool::ReverseIfNeeded(theVInternalCurve);*/
		// ���ó��ε� Compatible
		std::for_each(theUInternalCurve.begin(), theUInternalCurve.end(), UniformCurve);
		std::for_each(theVInternalCurve.begin(), theVInternalCurve.end(), UniformCurve);

		CurveOperate::CompatibleWithInterPoints(theVInternalCurve, theUInternalCurve);
		CurveOperate::CompatibleWithInterPoints(theUInternalCurve, theVInternalCurve);

		std::vector<gp_Pnt> upoints, vpoints;
		std::vector<Standard_Real> uparams, vparams;
		std::tie(upoints, uparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(theUInternalCurve, theVInternalCurve[0]);
		std::tie(vpoints, vparams) = CurveOperate::CalCurvesInterPointsParamsToCurve(theVInternalCurve, theUInternalCurve[0]);
		// �������
		std::vector<std::pair<double, Handle(Geom_BSplineCurve)>> combinedv;
		for (size_t i = 0; i < vparams.size(); ++i) 
		{
			combinedv.emplace_back(vparams[i], theVInternalCurve[i]);
		}
		std::sort(combinedv.begin(), combinedv.end(), [](const auto& a, const auto& b)
			{
			return a.first < b.first;
		});
		for (size_t i = 0; i < combinedv.size(); ++i)
		{
			vparams[i] = combinedv[i].first;
			theVInternalCurve[i] = combinedv[i].second;
		}

		std::vector<std::pair<double, Handle(Geom_BSplineCurve)>> combinedu;
		for (size_t i = 0; i < uparams.size(); ++i) 
		{
			combinedu.emplace_back(uparams[i], theUInternalCurve[i]);
		}
		std::sort(combinedu.begin(), combinedu.end(), [](const auto& a, const auto& b) 
		{
			return a.first < b.first;
		});
		for (size_t i = 0; i < combinedu.size(); ++i)
		{
			uparams[i] = combinedu[i].first;
			theUInternalCurve[i] = combinedu[i].second;
		}
		GordenSurface::BuildMyGordonSurf(theUInternalCurve, theVInternalCurve, uparams, vparams, aGordenFace);


		// �����ɵ���ת��ΪBSplineSurface
		Handle(Geom_Surface) aGeomSurface = BRep_Tool::Surface(aGordenFace);
		return Handle(Geom_BSplineSurface)::DownCast(aGeomSurface);
	}

	return nullptr;
}

PlanarCurve::PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, Standard_Real theTolerance)
	: curveType(CurveType::NOTPLANAR), curve(theCurve), line(), plane()
{
	IsPlanarCurve(curve, theTolerance);
}

PlanarCurve::PlanarCurve()
	:curveType(CurveType::NOTPLANAR), curve(), line(), plane()
{
	curveType = CurveType::NOTPLANAR;
}

Standard_Boolean PlanarCurve::IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, Standard_Real theTolerance)
{
	// ��������Ƿ�Ϊ��
	if (theCurve.IsNull())
	{
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	Standard_Boolean isLinear = IsBSplineCurveLinear(theCurve);
	if (isLinear)
	{
		curveType = CurveType::LINEAR;
		line = gp_Lin(theCurve->StartPoint(), gp_Vec(theCurve->StartPoint(), theCurve->EndPoint()));
		return true;
	}

	Standard_Boolean isPoint = IsBSplineCurvePoint(theCurve);
	if (isPoint)
	{
		curveType = CurveType::POINT;
		return true;
	}

	// ���������ϵĵ�
	std::vector<gp_Pnt> aBoundarySampling;
	Standard_Integer aNumSamples = 100; // �������������ɸ�����Ҫ����
	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();
	Standard_Real step = (lastParam - firstParam) / (aNumSamples - 1);

	aBoundarySampling.reserve(aNumSamples);
	for (Standard_Integer i = 0; i < aNumSamples; ++i)
	{
		Standard_Real aParam = firstParam + i * step;
		gp_Pnt aPnt;
		theCurve->D0(aParam, aPnt); // ��ȡ�����ϵĵ�
		aBoundarySampling.push_back(aPnt);
	}

	// ����Ƿ�ɹ�����
	if (aBoundarySampling.empty())
	{
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	// ������������������
	gp_XYZ N, X, Y;
	// �������� P	
	gp_Pnt P;
	P.SetCoord(0, 0, 0);
	Standard_Integer num = aBoundarySampling.size();
	for (const auto& point : aBoundarySampling)
	{
		P.SetCoord(P.X() + point.X(), P.Y() + point.Y(), P.Z() + point.Z());
	}
	P.SetCoord(P.X() / num, P.Y() / num, P.Z() / num);

	// �滻Э������󹹽�������ֵ�ֽⲿ��
	Standard_Integer m = 3;
	Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(m, m);

	Standard_Real c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
	gp_Pnt point; // ���� 'point'

	for (Standard_Integer i = 0; i < num; i++)
	{
		point = aBoundarySampling[i];

		c00 += (point.X() - P.X()) * (point.X() - P.X());
		c01 += (point.X() - P.X()) * (point.Y() - P.Y());
		c02 += (point.X() - P.X()) * (point.Z() - P.Z());

		c11 += (point.Y() - P.Y()) * (point.Y() - P.Y());
		c12 += (point.Y() - P.Y()) * (point.Z() - P.Z());

		c22 += (point.Z() - P.Z()) * (point.Z() - P.Z());
	}

	A1(0, 0) = c00;
	A1(0, 1) = c01;
	A1(0, 2) = c02;

	A1(1, 0) = c01;
	A1(1, 1) = c11;
	A1(1, 2) = c12;

	A1(2, 0) = c02;
	A1(2, 1) = c12;
	A1(2, 2) = c22;


	Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A1);
	Eigen::VectorXcd E1 = eigensolver.eigenvalues();
	auto  E2 = eigensolver.eigenvectors();

	auto eigen1 = E1.col(0)[0];
	auto eigen2 = E1.col(0)[1];
	auto eigen3 = E1.col(0)[2];

	Standard_Real eigenvalue1 = E1.col(0)[0].real();
	Standard_Real eigenvalue2 = E1.col(0)[1].real();
	Standard_Real eigenvalue3 = E1.col(0)[2].real();

	// ������С����ֵȷ����������������
	if (eigenvalue1 < eigenvalue2 && eigenvalue1 < eigenvalue3)
	{
		N.SetCoord(E2.col(0)[0].real(), E2.col(0)[1].real(), E2.col(0)[2].real());
		X.SetCoord(E2.col(1)[0].real(), E2.col(1)[1].real(), E2.col(1)[2].real());
		Y.SetCoord(E2.col(2)[0].real(), E2.col(2)[1].real(), E2.col(2)[2].real());
	}
	else if (eigenvalue2 < eigenvalue1 && eigenvalue2 < eigenvalue3)
	{
		N.SetCoord(E2.col(1)[0].real(), E2.col(1)[1].real(), E2.col(1)[2].real());
		Y.SetCoord(E2.col(0)[0].real(), E2.col(0)[1].real(), E2.col(0)[2].real());
		X.SetCoord(E2.col(2)[0].real(), E2.col(2)[1].real(), E2.col(2)[2].real());
	}
	else
	{
		N.SetCoord(E2.col(2)[0].real(), E2.col(2)[1].real(), E2.col(2)[2].real());
		Y.SetCoord(E2.col(1)[0].real(), E2.col(1)[1].real(), E2.col(1)[2].real());
		X.SetCoord(E2.col(0)[0].real(), E2.col(0)[1].real(), E2.col(0)[2].real());
	}

	// ���� gp_Dir �������� gp_pln
	gp_Dir dirN(N);
	plane = gp_Pln(P, dirN);
	// ����ƽ����
	Standard_Real maxDistance = 0.0;
	Standard_Real sumDistance = 0.0;
	for (const auto& point : aBoundarySampling)
	{
		Standard_Real distance = MathTool::ComputeDistancePointToPlane(point, plane);
		sumDistance += distance;
		if (distance > maxDistance)
		{
			maxDistance = distance;
		}
	}

	Standard_Real averageDistance = sumDistance / aBoundarySampling.size();

	// ����ƽ������ֵ�����ݾ������������
	if (maxDistance < theTolerance)
	{
		curveType = CurveType::PLANAR;
		return true;
	}
	else
	{
		curveType = CurveType::NOTPLANAR;
		return false;
	}
}

Standard_Boolean PlanarCurve::IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve, Standard_Real theTolerance)
{
	if (theCurve.IsNull())
	{
		return false;
	}

	// ��ȡ���Ƶ�����
	Standard_Integer aNumPoles = theCurve->NbPoles();
	if (aNumPoles < 2)
	{
		// �����������Ƶ㣬�޷�ȷ��ֱ��
		return false;
	}

	// ��ȡ��һ���͵ڶ������Ƶ�
	gp_Pnt aP0 = theCurve->Pole(1);
	gp_Pnt aP1 = theCurve->Pole(2);

	// ���㷽������ v = P1 - P0
	gp_Vec aV(aP0, aP1);
	if (aV.Magnitude() < theTolerance)
	{
		// ǰ�������غϣ��޷����巽��
		return false;
	}

	// ����ÿһ�������Ŀ��Ƶ㣬����Ƿ��뷽����������
	for (Standard_Integer i = 3; i <= aNumPoles; ++i)
	{
		gp_Pnt aPi = theCurve->Pole(i);
		gp_Vec aU(aP0, aPi);

		// ������ v �� u
		gp_Vec aCross = aV.Crossed(aU);

		// �������ģ���Ƿ����ݲΧ��
		if (aCross.Magnitude() > theTolerance)
		{
			// ������
			return false;
		}
	}
	// ���п��Ƶ㹲��
	return true;
}


Standard_Boolean PlanarCurve::IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, Standard_Real theTolerance)
{
	// ��ȡ���Ƶ�����
	Standard_Integer aNumPoles = theCurve->NbPoles();
	// ��ȡ��һ���͵ڶ������Ƶ�
	gp_Pnt aP0 = theCurve->Pole(1);
	gp_Pnt aP1 = theCurve->Pole(2);

	// ���㷽������ v = P1 - P0
	gp_Vec aV(aP0, aP1);
	if (aNumPoles == 2 && aV.Magnitude() < theTolerance)
	{
		return true;
	}
	return false;
}

Standard_Real MathTool::ComputeCurveCurveDistance(const Handle(Geom_BSplineCurve)& theCurve, const Handle(Geom_BSplineCurve)& theBoundaryCurve)
{
	GeomAPI_ExtremaCurveCurve aExtrema(theCurve, theBoundaryCurve);
	// ����Ƿ��ҵ��˼�ֵ��
	if (aExtrema.NbExtrema() > 0)
	{
		// �������м�ֵ�㣬�ҵ���С����
		Standard_Real aMinDistance = RealLast();
		for (Standard_Integer i = 1; i <= aExtrema.NbExtrema(); ++i)
		{
			Standard_Real aDist = aExtrema.Distance(i);
			if (aDist < aMinDistance)
			{
				aMinDistance = aDist;
			}
		}
		return aMinDistance;
	}

	return INT_MAX;
}

// �������ߵĲ�����ƽ������
gp_Pnt MathTool::ComputeAverageSamplePoint(const Handle(Geom_BSplineCurve)& theCurve, Standard_Integer theNumSamples)
{
	GeomAdaptor_Curve adaptor(theCurve);
	Standard_Real aStartParam = adaptor.FirstParameter();
	Standard_Real aEndParam = adaptor.LastParameter();
	Standard_Real aDeltaParam = (aEndParam - aStartParam) / (theNumSamples - 1);
	Standard_Real X = 0, Y = 0, Z = 0;
	for (Standard_Integer i = 0; i < theNumSamples; ++i) {
		Standard_Real aParam = aStartParam + i * aDeltaParam;
		gp_Pnt sample = adaptor.Value(aParam);
		X += sample.X();
		Y += sample.Y();
		Z += sample.Z();
	}
	X /= theNumSamples;
	Y /= theNumSamples;
	Z /= theNumSamples;
	return gp_Pnt(X, Y, Z);
}

Standard_Real MathTool::ComputeAngleWithAxis(const gp_Vec& theVec, const gp_Vec& theAxis)
{
	Standard_Real aDotProduct = theVec.Dot(theAxis); // ���
	Standard_Real aMagnitudeVec = theVec.Magnitude(); // ������ģ
	Standard_Real aMagnitudeAxis = theAxis.Magnitude(); // ���ģ
	Standard_Real aCosAngle = aDotProduct / (aMagnitudeVec * aMagnitudeAxis); // ����ֵ
	// ��ֹ��ֵ���³��� [-1, 1] ��Χ
	aCosAngle = std::max(-1.0, std::min(1.0, aCosAngle));
	return std::acos(aCosAngle); // ���ؼн�
}

void MathTool::CheckSelfIntersect(std::vector<Handle(Geom_BSplineCurve)> theBSplineCurvesArray)
{
	for (Standard_Integer i = 0; i < theBSplineCurvesArray.size(); i++)
	{
		for (Standard_Integer j = i + 1; j < theBSplineCurvesArray.size(); j++)
		{
			Standard_Real aDistance = MathTool::ComputeCurveCurveDistance(theBSplineCurvesArray[i], theBSplineCurvesArray[j]);
			if (aDistance < 1e-3)
			{
				// �����Խ�������������ƽ�����������
				gp_Pnt aPnt1 = MathTool::ComputeAverageSamplePoint(theBSplineCurvesArray[i], 10);
				gp_Pnt aPnt2 = MathTool::ComputeAverageSamplePoint(theBSplineCurvesArray[j], 10);

				// �������ߵ�ԭ�������
				gp_Vec aVec1 = aPnt1.XYZ() - gp_Pnt(0, 0, 0).XYZ();
				gp_Vec aVec2 = aPnt2.XYZ() - gp_Pnt(0, 0, 0).XYZ();

				// ����x�ᡢy�ᡢz��
				gp_Vec xAxis(1, 0, 0);
				gp_Vec yAxis(0, 1, 0);
				gp_Vec zAxis(0, 0, 1);

				// ��������������x��y��z��ļн�
				Standard_Real angle1_x = ComputeAngleWithAxis(aVec1, xAxis);
				Standard_Real angle2_x = ComputeAngleWithAxis(aVec2, xAxis);
				Standard_Real angle1_y = ComputeAngleWithAxis(aVec1, yAxis);
				Standard_Real angle2_y = ComputeAngleWithAxis(aVec2, yAxis);
				Standard_Real angle1_z = ComputeAngleWithAxis(aVec1, zAxis);
				Standard_Real angle2_z = ComputeAngleWithAxis(aVec2, zAxis);

				// �Ƚ������ļнǣ������нǸ�С������
				if (angle1_x < angle2_x || angle1_y < angle2_y || angle1_z < angle2_z)
				{
					// ���� theBSplineCurvesArray[i]���Ƴ� theBSplineCurvesArray[j]
					theBSplineCurvesArray.erase(theBSplineCurvesArray.begin() + j);
					j--; // ��֤��������һ��Ԫ��
				}
				else
				{
					// ���� theBSplineCurvesArray[j]���Ƴ� theBSplineCurvesArray[i]
					theBSplineCurvesArray.erase(theBSplineCurvesArray.begin() + i);
					i--; // ��֤��������һ��Ԫ��
					break; // �˳��ڲ�ѭ������Ϊ i �Ѿ��ı�
				}
			}
		}
	}
}

gp_Dir MathTool::ComputeAverageTangent(const Handle(Geom_BSplineCurve)& theCurve, Standard_Integer theNumSamples)
{
	if (theCurve.IsNull())
	{
		throw std::invalid_argument("Curve is null.");
	}
	if (theNumSamples <= 0)
	{
		throw std::invalid_argument("Number of samples must be positive.");
	}

	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();
	Standard_Real step = (lastParam - firstParam) / (theNumSamples - 1);

	gp_Vec sumTangent(0.0, 0.0, 0.0);
	Standard_Integer validSamples = 0;

	for (Standard_Integer i = 0; i < theNumSamples; ++i)
	{
		Standard_Real param = firstParam + i * step;
		gp_Pnt pnt;
		gp_Vec tangent;
		theCurve->D1(param, pnt, tangent); // D1 ��ȡ���һ�׵���
		sumTangent += tangent;
		++validSamples;
	}

	if (validSamples == 0)
	{
		throw std::runtime_error("No valid samples were taken from the curve.");
	}

	gp_Vec averageTangent = sumTangent / validSamples;
	gp_Dir averageDir(averageTangent);

	return averageDir;
}

Standard_Real MathTool::ComputeAngleBetweenCurves(Handle(Geom_BSplineCurve)& theCurve1,
	Handle(Geom_BSplineCurve)& theCurve2,
	Standard_Integer theNumSamples)
{
	PlanarCurve aPlanarCurve1(theCurve1);
	PlanarCurve aPlanarCurve2(theCurve2);
	if (aPlanarCurve1.GetCurveType() != CurveType::NOTPLANAR && aPlanarCurve2.GetCurveType() != CurveType::NOTPLANAR)
	{
		return MathTool::ComputeAngleBetweenPlanarCurves(aPlanarCurve1, aPlanarCurve2);
	}

	gp_Dir aAvgDir1 = ComputeAverageTangent(theCurve1, theNumSamples);
	gp_Dir aAvgDir2 = ComputeAverageTangent(theCurve2, theNumSamples);

	Standard_Real aDotProduct = aAvgDir1.Dot(aAvgDir2);

	// ȷ������� [-1, 1] ��Χ�ڣ��Ա�����ֵ���
	aDotProduct = std::max(-1.0, std::min(1.0, aDotProduct));

	Standard_Real aAngleRad = std::acos(aDotProduct);
	Standard_Real aAngleDeg = aAngleRad * 180.0 / M_PI;

	return aAngleDeg;
}

void MathTool::SortBSplineCurves(std::vector<Handle(Geom_BSplineCurve)>& theCurves,
	Handle(Geom_BSplineCurve) theReferCurve)
{
	Standard_Integer aNumSamples = 10;
	// ������������Ƿ�Ϊ��
	if (theCurves.empty())
	{
		return;
	}

	// ����ο��㣬ʹ�������еĵ�һ������
	gp_Pnt aReferencePoint = MathTool::ComputeAverageSamplePoint(theReferCurve, aNumSamples);

	// ʹ�� std::sort �����������������
	std::sort(theCurves.begin(), theCurves.end(),
		[&](const Handle(Geom_BSplineCurve)& theCurve1, const Handle(Geom_BSplineCurve)& theCurve2) -> Standard_Boolean
		{
			// ����ÿ�����ߵ�ƽ��������
			gp_Pnt aAvg = MathTool::ComputeAverageSamplePoint(theCurve1, aNumSamples);
			gp_Pnt bAvg = MathTool::ComputeAverageSamplePoint(theCurve2, aNumSamples);

			// ����ƽ�������㵽�ο���ľ���
			Standard_Real distA = aReferencePoint.Distance(aAvg);
			Standard_Real distB = aReferencePoint.Distance(bAvg);

			// ���վ����С��������
			return distA < distB;
		}
	);

	std::cout << "���������Ѹ���ƽ�������㵽�ο���ľ���ɹ�����" << std::endl;
}

void MathTool::ReverseIfNeeded(std::vector<Handle(Geom_BSplineCurve)>& theCurves)
{
	// ������������Ƿ�Ϊ��
	if (theCurves.empty()) return;

	// ��ȡ��һ��������Ϊ��ʼ�ο�
	Handle(Geom_BSplineCurve) firstCurve = theCurves[0];
	gp_Pnt firstStart = firstCurve->StartPoint();
	gp_Pnt firstEnd = firstCurve->EndPoint();
	gp_Vec firstDirection(firstStart, firstEnd); // ����ο����ߵķ�������

	// �����һ��������һ���㣨��������Ϊ������������ѡ�����һ��������Ϊ�ο�
	if (firstDirection.Magnitude() == 0)
	{
		if (theCurves.size() < 2)
		{
			return;
		}
		firstCurve = theCurves.back();
		firstStart = firstCurve->StartPoint();
		firstEnd = firstCurve->EndPoint();
		firstDirection = gp_Vec(firstStart, firstEnd); // ���¼���ο����ߵķ�������

		// �ٴμ�鷽�������Ƿ�Ϊ������
		if (firstDirection.Magnitude() == 0)
		{
			return;
		}
	}

	// ����ÿ�����ߣ���鷽�򲢷�ת
	for (auto& aCurve : theCurves)
	{
		gp_Pnt start = aCurve->StartPoint();
		gp_Pnt end = aCurve->EndPoint();
		gp_Vec direction(start, end); // ���㵱ǰ���ߵķ�������

		// �����ǰ���ߵķ�����ο������෴����ת����
		if (direction.Dot(firstDirection) < 0)
		{
			aCurve->Reverse(); // ��ת���߷���
		}
	}
}


Standard_Real MathTool::ComputeAngleBetweenLineAndPlane(const gp_Lin& theLine, const gp_Pln& thePlane)
{
	// ��ȡֱ�ߵķ�������
	gp_Vec aLineDirection = theLine.Direction();  // ֱ�߷�������

	// ��ȡƽ��ķ�����
	gp_Vec aPlaneNormal = thePlane.Axis().Direction();  // ƽ�淨����

	// ����ֱ�߷���������ƽ�淨�����ĵ��
	Standard_Real aDotProduct = aLineDirection.Dot(aPlaneNormal);

	// ����ֱ�߷��������ͷ�������ģ��
	Standard_Real aLineMagnitude = aLineDirection.Magnitude();
	Standard_Real aPlaneNormalMagnitude = aPlaneNormal.Magnitude();

	// ����нǣ����ȣ�
	Standard_Real aAngle = std::acos(aDotProduct / (aLineMagnitude * aPlaneNormalMagnitude));

	return aAngle;
}

Standard_Real MathTool::ComputeAngleBetweenLines(const gp_Lin& theLine1, const gp_Lin& theLine2)
{
	// ��ȡ����ֱ�ߵķ�������
	gp_Vec direction1 = theLine1.Direction();
	gp_Vec direction2 = theLine2.Direction();

	// ���㷽�������ĵ��
	Standard_Real aDotProduct = direction1.Dot(direction2);

	// ���㷽��������ģ��
	Standard_Real magnitude1 = direction1.Magnitude();
	Standard_Real magnitude2 = direction2.Magnitude();

	// ����нǣ����ȣ�
	Standard_Real aAngle = std::acos(aDotProduct / (magnitude1 * magnitude2));

	return aAngle;
}

Standard_Real MathTool::ComputeAngleBetweenPlanes(const gp_Pln& thePlane1, const gp_Pln& thePlane2)
{
	// ��ȡƽ��ķ�����
	gp_Dir aNormal1 = thePlane1.Axis().Direction();
	gp_Dir aNormal2 = thePlane2.Axis().Direction();

	// ���㷨����֮��ĵ��
	Standard_Real aDotProduct = aNormal1.Dot(aNormal2);

	// ���� aDotProduct �ķ�Χ�� [-1, 1] ֮�䣬�Է�ֹ��ֵ���µ� acos �������
	aDotProduct = std::max(-1.0, std::min(1.0, aDotProduct));

	// ����нǵĻ���ֵ
	Standard_Real aAngleRad = std::acos(aDotProduct);

	// ������ת��Ϊ����
	Standard_Real aAngleDeg = aAngleRad * 180.0 / M_PI;

	// ��ȡ��ǣ�0�㵽90�㣩
	if (aAngleDeg > 90.0)
	{
		aAngleDeg = 180.0 - aAngleDeg;
	}

	return aAngleDeg;
}

Standard_Real MathTool::ComputeDistancePointToLine(const gp_Pnt& thePoint, const gp_Lin& theLine)
{
	// ֱ�ߵķ�������
	gp_Vec aLineDirection = theLine.Direction();

	// �㵽ֱ���ϵ�����һ�㣨����ѡ��ֱ���ϵ�һ������Ϊ�ο��㣩
	gp_Pnt aLinePoint = theLine.Location(); // ֱ���ϵ�һ����

	// ����ֱ���ϵĵ������
	gp_Vec pointToLine = thePoint.XYZ() - aLinePoint.XYZ();

	// ����㵽ֱ�ߵĴ�ֱ����
	gp_Vec aCrossProduct = pointToLine.Crossed(aLineDirection);

	// ���ؾ��룬ʹ�ò�˵�ģ�����Է���������ģ��
	return aCrossProduct.Magnitude() / aLineDirection.Magnitude();
}

Standard_Real MathTool::ComputeAngleBetweenPlanarCurves(const PlanarCurve& theCurve1, const PlanarCurve& theCurve2)
{
	Standard_Real aAngle = INT_MAX;
	Standard_Real aDistance = INT_MIN;
	if (theCurve1.GetCurveType() == CurveType::LINEAR && theCurve2.GetCurveType() == CurveType::LINEAR)
	{
		// ����ֱ��ֱ����
		aAngle = ComputeAngleBetweenLines(theCurve1.GetLine(), theCurve2.GetLine());
	}
	else if (theCurve1.GetCurveType() == CurveType::PLANAR && theCurve2.GetCurveType() == CurveType::PLANAR)
	{
		// ����ƽ��ƽ����
		aAngle = ComputeAngleBetweenPlanes(theCurve1.GetPlane(), theCurve2.GetPlane());
	}
	else if (theCurve1.GetCurveType() == CurveType::PLANAR && theCurve2.GetCurveType() == CurveType::LINEAR)
	{
		// ����ƽ���ֱ����
		aAngle = ComputeAngleBetweenLineAndPlane(theCurve2.GetLine(), theCurve1.GetPlane());
	}
	else if (theCurve1.GetCurveType() == CurveType::LINEAR && theCurve2.GetCurveType() == CurveType::PLANAR)
	{
		// ����ƽ���ֱ����
		aAngle = ComputeAngleBetweenLineAndPlane(theCurve1.GetLine(), theCurve2.GetPlane());
	}
	else if (theCurve1.GetCurveType() == CurveType::POINT && theCurve2.GetCurveType() == CurveType::POINT)
	{
		aDistance = theCurve1.GetPoint().Distance(theCurve2.GetPoint());
	}
	else if (theCurve1.GetCurveType() == CurveType::POINT && theCurve2.GetCurveType() == CurveType::PLANAR)
	{
		aDistance = ComputeDistancePointToPlane(theCurve1.GetPoint(), theCurve2.GetPlane());
	}
	else if (theCurve1.GetCurveType() == CurveType::POINT && theCurve2.GetCurveType() == CurveType::LINEAR)
	{
		aDistance = ComputeDistancePointToLine(theCurve1.GetPoint(), theCurve2.GetLine());
	}
	else if (theCurve1.GetCurveType() == CurveType::PLANAR && theCurve2.GetCurveType() == CurveType::POINT)
	{
		aDistance = ComputeDistancePointToPlane(theCurve2.GetPoint(), theCurve1.GetPlane());
	}
	else if (theCurve1.GetCurveType() == CurveType::LINEAR && theCurve2.GetCurveType() == CurveType::POINT)
	{
		aDistance = ComputeDistancePointToLine(theCurve2.GetPoint(), theCurve1.GetLine());
	}

	if (aDistance > 1) aAngle = 0.0;

	return aAngle;
}

std::vector<gp_Pnt> MathTool::GetSamplePointsOnCurve(const Handle(Geom_Curve)& theCurve, Standard_Integer theNumPoints)
{
	std::vector<gp_Pnt> aSampledPoints;

	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();
	Standard_Real delta = (lastParam - firstParam) / (theNumPoints - 1);

	for (Standard_Integer i = 0; i < theNumPoints; ++i)
	{
		Standard_Real param = firstParam + i * delta;
		gp_Pnt aP;
		theCurve->D0(param, aP);
		aSampledPoints.push_back(aP);
	}

	return aSampledPoints;
}

void MathTool::TrimInternalCurves(
	std::vector<Handle(Geom_BSplineCurve)>& theInternalBSplineCurves,
	const std::vector<Handle(Geom_BSplineCurve)>& theBoundaryCurveArray,
	Standard_Real theToleranceDistance) {
	for (auto& internalCurve : theInternalBSplineCurves)
	{
		// ���ڴ洢��������߽����ߵĽ���Ͳü�����
		gp_Pnt replacePoints[2];
		Standard_Real splitParams[2] = { 0 };
		Standard_Integer foundCount = 0;

		// �����ڲ�������ÿ���߽����ߵľ��벢����
		std::vector<std::pair<Standard_Real, Handle(Geom_BSplineCurve)>> aBoundaryCurves;
		for (auto& boundaryCurve : theBoundaryCurveArray)
		{
			Standard_Real distance = MathTool::ComputeCurveCurveDistance(internalCurve, boundaryCurve);
			aBoundaryCurves.emplace_back(distance, boundaryCurve);
		}

		// �������С��������
		std::sort(aBoundaryCurves.begin(), aBoundaryCurves.end(),
			[](const auto& a, const auto& b) { return a.first < b.first; });

		if (aBoundaryCurves[0].first >= theToleranceDistance) continue;

		// �ҵ����ڲ����߾�������������߽����ߵĽ���
		for (size_t i = 0; i < 2 && i < aBoundaryCurves.size(); ++i)
		{
			GeomAPI_ExtremaCurveCurve extrema(internalCurve, aBoundaryCurves[i].second);

			if (extrema.NbExtrema() > 0)
			{
				gp_Pnt internalPoint, boundaryPoint;
				Standard_Real paramOnCurve;
				extrema.NearestPoints(internalPoint, boundaryPoint);
				extrema.LowerDistanceParameters(splitParams[foundCount], paramOnCurve);

				replacePoints[foundCount] = boundaryPoint;
				foundCount++;
			}
		}

		// ���û���ҵ�������Ч�Ľ��㣬������
		if (foundCount != 2)
			continue;

		// ȷ���ü���������������
		if (splitParams[0] > splitParams[1])
		{
			std::swap(splitParams[0], splitParams[1]);
			std::swap(replacePoints[0], replacePoints[1]);
		}

		// �ü��ڲ����߲�����Ϊ�µ�B��������
		Handle(Geom_TrimmedCurve) trimmedCurve = new Geom_TrimmedCurve(internalCurve, splitParams[0], splitParams[1]);
		Handle(Geom_BSplineCurve) modifiedCurve = GeomConvert::CurveToBSplineCurve(trimmedCurve, Convert_TgtThetaOver2);

		// ���òü���Ķ˵�
		modifiedCurve->SetPole(1, replacePoints[0]);
		modifiedCurve->SetPole(modifiedCurve->NbPoles(), replacePoints[1]);

		// �����ڲ�����
		internalCurve = modifiedCurve;
	}
}
void MathTool::SortPoints(std::vector<gp_Pnt>& thePoints, const gp_Pnt& theReferPoint)
{
	// ���������Ƿ�Ϊ��
	if (thePoints.empty())
	{
		return;
	}

	// ʹ�� std::sort �Ե������������
	std::sort(thePoints.begin(), thePoints.end(),
		[&](const gp_Pnt& aPnt1, const gp_Pnt& aPnt2) -> Standard_Boolean
		{
			// ����ÿ���㵽�ο���ľ���
			Standard_Real distA = aPnt1.Distance(theReferPoint);
			Standard_Real distB = aPnt2.Distance(theReferPoint);

			// ���վ����С��������
			return distA < distB;
		}
	);
}

Standard_Real MathTool::ComputeDistancePointToPlane(const gp_Pnt& theP, const gp_Pln& thePlane)
{
	// ��ȡƽ��ķ������͵�
	gp_Dir aNormal = thePlane.Axis().Direction();
	gp_Pnt aPlanePoint = thePlane.Location();

	// ���� (p - aPlanePoint)
	gp_Vec vec(theP.X() - aPlanePoint.X(), theP.Y() - aPlanePoint.Y(), theP.Z() - aPlanePoint.Z());

	// �㵽ƽ��ľ��� = |vec . aNormal| / |aNormal|
	// ���� aNormal �ǵ�λ��������ĸΪ1
	Standard_Real aDistance = std::abs(vec.Dot(aNormal));
	return aDistance;
}


//�����Բ�
/*
* 1����������ĸ��д
* 2�������� ��the
* 3�������� Сд
* 4������֪��
* 5����������
* 6����������ע��
* 7�����������δʹ���β�
* 8��������δ���ñ���
* 9���޾���
*
*/

//�궨���쳣����
#define HANDLE_EXCEPTIONS_CONTINUE \
    catch (const std::exception & e) { \
        std::cerr << "Exception caught: " << e.what() << std::endl; \
        continue; \
    } catch (...) { \
        std::cerr << "Unknown exception caught." << std::endl; \
        continue; \
    }

#define HANDLE_EXCEPTIONS_RETURN_FALSE \
    catch (const std::exception & e) { \
        std::cerr << "Exception caught: " << e.what() << std::endl; \
        return false; \
    } catch (...) { \
        std::cerr << "Unknown exception caught." << std::endl; \
        return false; \
    }


// �����������ж��������Ƿ񼸺���ͬ������С�ڸ������ݲ
bool ArePointsEqual(const gp_Pnt& thePoint1, const gp_Pnt& thePoint2, Standard_Real theTolerance = 1e-6) {
	return thePoint1.Distance(thePoint2) <= theTolerance;
}

// �����������ж�����ֵ�Ƿ����ݲΧ�����
bool AreValuesEqual(Standard_Real theValue1, Standard_Real theValue2, Standard_Real theTolerance = 0.01) {
	return std::fabs(theValue1 - theValue2) <= theTolerance;
}

/// <summary>
/// �ϲ�theCurves�Ľڵ㡢��theTolerance����Ϊһ���ڵ�
/// </summary>
/// <param name="theCurves"></param>
/// <param name="theTolerance"></param>
/// <returns></returns> �ϲ���Ľڵ�
std::vector<Standard_Real> CalSameKnotFromCurves(std::vector< Handle(Geom_BSplineCurve) >& theCurves, Standard_Real theTolerance = 0.01) {
	if (theCurves.empty()) return {};

	// ʹ��map�Զ�ȥ�ز��ϲ�
	std::map<Standard_Real, Standard_Integer> knotMap;

	for (const auto& curve : theCurves) {
		if (curve.IsNull()) continue;

		Standard_Integer nbKnots = curve->NbKnots();
		for (Standard_Integer i = 1; i <= nbKnots; ++i) {
			Standard_Real knot = curve->Knot(i);
			Standard_Integer multiplicity = curve->Multiplicity(i);

			auto it = knotMap.lower_bound(knot - theTolerance);
			if (it != knotMap.end() && std::fabs(it->first - knot) <= theTolerance) {
				// �����ظ���
				it->second = std::max(it->second, multiplicity);
			}
			else {
				knotMap[knot] = multiplicity;
			}
		}
	}

	// չ�����
	std::vector<Standard_Real> result;
	for (const auto& [knot, multiplicity] : knotMap) {
		result.insert(result.end(), multiplicity, knot);
	}

	return result;
}



/// <summary>
/// ����theCurve��theCurves����Ĳ���ֵ������Ϊtuple����ʽ��һһ��Ӧ������theTolerance�ж�ֱ���Ƿ�ȷ�ཻ
/// ����������ߺ�һ�����ߵĽ����Լ�����
/// </summary>
/// <param name="theCurves"></param> �����������
/// <param name="theCurve"></param>	 ĳ������
/// <param name="theTolerance"></param>
/// <returns></returns>	�����Լ��������ֵ
std::tuple<std::vector<gp_Pnt>, std::vector<Standard_Real>> CurveOperate::CalCurvesInterPointsParamsToCurve(
	const std::vector<Handle(Geom_BSplineCurve)>& theCurves,
	const Handle(Geom_BSplineCurve)& theCurve,
	Standard_Real theTolerance) {
	std::vector<gp_Pnt> pointsOnTheCurve;
	std::vector<Standard_Real> paramsOnTheCurve;
	const Standard_Real intersectionToleranceSq = theTolerance * theTolerance; // �ݲ��ƽ�������ھ���Ƚ�

	// ��������Ƿ���Ч
	if (theCurve.IsNull()) {
		std::cerr << "Error: theCurve is null." << std::endl;
		return std::make_tuple(pointsOnTheCurve, paramsOnTheCurve);
	}

	// ����ÿ������
	for (const auto& curve : theCurves) {
		if (curve.IsNull()) {
			// ����������
			continue;
		}

		try {
			// ʹ�� GeomAPI_ExtremaCurveCurve ������������֮��ļ�ֵ��
			GeomAPI_ExtremaCurveCurve extrema(theCurve, curve);

			if (extrema.NbExtrema() > 0) {
				bool hasIntersection = false;
				std::vector<gp_Pnt> intersectionsForThisCurve;
				std::vector<Standard_Real> paramsForThisCurve;

				// �������м�ֵ�㣬���ҽ���
				for (Standard_Integer i = 1; i <= extrema.NbExtrema(); ++i) {
					Standard_Real distanceSq = extrema.Distance(i);
					if (distanceSq <= intersectionToleranceSq) { // �ж��Ƿ�Ϊ����
						Standard_Real U1, U2;
						extrema.Parameters(i, U1, U2); // ��ȡ theCurve �� curve �ϵĲ���ֵ

						gp_Pnt p1 = theCurve->Value(U1); // ��ȡ theCurve �ϵĵ�

						// ����Ƿ��Ѿ����ڼ�����ͬ�ĵ㣬�����ظ�
						bool alreadyExists = false;
						for (const auto& existingPnt : intersectionsForThisCurve) {
							if (ArePointsEqual(p1, existingPnt, theTolerance)) {
								alreadyExists = true;
								break;
							}
						}

						if (!alreadyExists) {
							intersectionsForThisCurve.emplace_back(p1);
							paramsForThisCurve.emplace_back(U1);
							hasIntersection = true;
						}
					}
				}

				if (hasIntersection) {
					// ��������ҵ��Ľ���
					pointsOnTheCurve.insert(pointsOnTheCurve.end(), intersectionsForThisCurve.begin(), intersectionsForThisCurve.end());
					paramsOnTheCurve.insert(paramsOnTheCurve.end(), paramsForThisCurve.begin(), paramsForThisCurve.end());
				}
				else {
					// û�н��㣬�ҵ������
					Standard_Real U1, U2;
					bool hasNearest = extrema.TotalLowerDistanceParameters(U1, U2);
					if (hasNearest) {
						gp_Pnt nearestP1 = theCurve->Value(U1);

						// ����Ƿ��Ѿ����ڼ�����ͬ�ĵ㣬�����ظ�
						bool alreadyExists = false;
						for (const auto& existingPnt : pointsOnTheCurve) {
							if (ArePointsEqual(nearestP1, existingPnt, theTolerance)) {
								alreadyExists = true;
								break;
							}
						}

						if (!alreadyExists) {
							pointsOnTheCurve.emplace_back(nearestP1);
							paramsOnTheCurve.emplace_back(U1);
						}
					}
					else {
						// �޷��ҵ������
						std::cerr << "Warning: Unable to find nearest point parameters between two curves." << std::endl;
					}
				}
			}
			else {
				// NbExtrema == 0��û�м�ֵ�������ҵ������
				Standard_Real U1, U2;
				bool hasNearest = extrema.TotalLowerDistanceParameters(U1, U2);
				if (hasNearest) {
					gp_Pnt nearestP1 = theCurve->Value(U1);

					// ����Ƿ��Ѿ����ڼ�����ͬ�ĵ㣬�����ظ�
					bool alreadyExists = false;
					for (const auto& existingPnt : pointsOnTheCurve) {
						if (ArePointsEqual(nearestP1, existingPnt, theTolerance)) {
							alreadyExists = true;
							break;
						}
					}

					if (!alreadyExists) {
						pointsOnTheCurve.emplace_back(nearestP1);
						paramsOnTheCurve.emplace_back(U1);
					}
				}
				else {
					// �޷��ҵ������
					std::cerr << "Warning: Unable to find nearest point parameters between two curves." << std::endl;
				}
			}
		}
		catch (Standard_Failure& failure) {
			// ��׽�쳣������������һ������
			std::cerr << "Exception: " << failure.GetMessageString() << std::endl;
			continue;
		}
	}

	return std::make_tuple(pointsOnTheCurve, paramsOnTheCurve);
}



/// <summary>
/// ��theCurve��theOriginParams�Ĳ��������ڵĵ���ܲ��������ܴ���ΪtheSamplingNum�����ҷ��ؼ��ܺ�ĵ��Լ�����
/// </summary>
/// <param name="theCurve"></param>	�����ܵ�����
/// <param name="theOriginParams"></param> ���ܵ�����
/// <param name="theSamplingNum"></param> ������Ŀ
/// <returns></returns> ���ܺ�ĵ��Լ�����ֵ
std::tuple<std::vector<gp_Pnt>, std::vector<Standard_Real>> DenseSampling(const Handle(Geom_BSplineCurve)& theCurve,
	std::vector<Standard_Real>& theOriginParams, Standard_Integer theSamplingNum) {
	std::vector<gp_Pnt> sampledPoints;
	std::vector<Standard_Real> sampledParams;

	// ��� theOriginParams �ĺϷ���
	if (theOriginParams.size() < 2) {
		Standard_Failure::Raise("DenseSampling Error: originParams must contain at least two parameters.");
	}

	// ��ȡ���ߵĲ�����
	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();

	// ��� theOriginParams �ĵ�һ�������һ�������Ƿ����ݲΧ��
	if (!AreValuesEqual(theOriginParams.front(), firstParam, 0.01)) {
		Standard_Failure::Raise("DenseSampling Error: The first element of originParams is not within tolerance of theCurve's first parameter.");
	}

	if (!AreValuesEqual(theOriginParams.back(), lastParam, 0.01)) {
		Standard_Failure::Raise("DenseSampling Error: The last element of originParams is not within tolerance of theCurve's last parameter.");
	}
	theOriginParams.front() = firstParam;
	theOriginParams.back() = lastParam;

	// ��� theOriginParams �Ƿ���������
	for (size_t i = 1; i < theOriginParams.size(); ++i) {
		if (theOriginParams[i] < theOriginParams[i - 1]) {
			Standard_Failure::Raise("DenseSampling Error: originParams must be in ascending order.");
		}
	}

	// ���м��ܲ���
	for (size_t i = 0; i < theOriginParams.size() - 1; ++i) {
		Standard_Real aStartParam = theOriginParams[i];
		Standard_Real aEndParam = theOriginParams[i + 1];
		Standard_Real interval = aEndParam - aStartParam;

		// ����ʼ��������������
		if (i == 0) { // ���ڵ�һ��ѭ��ʱ��ӵ�һ����ʼ��
			sampledParams.emplace_back(aStartParam);
			sampledPoints.emplace_back(theCurve->Value(aStartParam));
		}

		// ����ÿ�������ڵĲ�������
		Standard_Real step = interval / (theSamplingNum + 1);

		// ���ɲ�����
		for (Standard_Integer j = 1; j <= theSamplingNum; ++j) {
			Standard_Real aSamplingParam = aStartParam + j * step;
			gp_Pnt aSamplingPoint = theCurve->Value(aSamplingParam);
			sampledParams.emplace_back(aSamplingParam);
			sampledPoints.emplace_back(aSamplingPoint);
		}

		// ��������������������
		sampledParams.emplace_back(aEndParam);
		sampledPoints.emplace_back(theCurve->Value(aEndParam));
	}

	return std::make_tuple(sampledPoints, sampledParams);
}

/// <summary>
/// ����������任�����Ӧ��������ܵ����Ҳ��Ҫ�仯�����ӿ�Ϊ�����������ȱ����仯������ܲ���
/// </summary>
/// <param name="theBaseParams"></param> �ο��Ĳ���ֵ����׼����ֵ��ͨ��ƽ���������
/// <param name="theBaseIndex"></param>	 ��Ӧ���ܵĲ����������±�
/// <param name="theParams"></param> ���ܺ�Ĳ���ֵ
/// <returns></returns>	�任��Ĳ���ֵ
std::vector<Standard_Real> ScalingParamsByBaseParams(const std::vector<Standard_Real>& theBaseParams, Standard_Integer theBaseIndex, std::vector<Standard_Real>& theParams) {
	std::vector<Standard_Real> scaledParams;

	// ��� (theParams.size() + 1) �Ƿ��ܱ� theBaseIndex ����
	if ((theParams.size() - 1) % theBaseIndex != 0) {
		throw std::invalid_argument("ScalingParamsByBaseParams Error: (params.size() + 1) is not divisible by baseIndex.");
	}

	// ���������� theBaseParams ��С
	size_t expectedBaseSize = ((theParams.size() - 1) / theBaseIndex) + 1;
	if (theBaseParams.size() != expectedBaseSize) {
		throw std::invalid_argument("ScalingParamsByBaseParams Error: baseParams size does not match the expected number of intervals.");
	}

	// ��� theBaseParams �ĵ�һ��Ԫ���Ƿ��� theParams �ĵ�һ��Ԫ�����ݲΧ�����
	if (!AreValuesEqual(theBaseParams.front(), theParams.front(), 0.01)) {
		throw std::invalid_argument("ScalingParamsByBaseParams Error: The first element of baseParams is not within tolerance of params[0].");
	}

	// ��� theBaseParams �����һ��Ԫ���Ƿ��� theParams �����һ����׼Ԫ�����ݲΧ�����
	// ��׼Ԫ�ص�λ��Ϊ (theBaseParams.size() -1) * theBaseIndex
	size_t lastParamPos = (theBaseParams.size() - 1) * theBaseIndex;
	if (lastParamPos >= theParams.size()) {
		throw std::invalid_argument("ScalingParamsByBaseParams Error: baseIndex exceeds params size.");
	}
	if (!AreValuesEqual(theBaseParams.back(), theParams[lastParamPos], 0.01)) {
		throw std::invalid_argument("ScalingParamsByBaseParams Error: The last element of baseParams is not within tolerance of params[" + std::to_string(lastParamPos) + "].");
	}

	// �������Ŵ���
	for (size_t i = 0; i < theBaseParams.size() - 1; ++i) {
		Standard_Real baseStart = theBaseParams[i];
		Standard_Real baseEnd = theBaseParams[i + 1];
		Standard_Real baseInterval = baseEnd - baseStart;

		size_t paramStartPos = i * theBaseIndex;
		size_t paramEndPos = (i + 1) * theBaseIndex;
		Standard_Real paramStart = theParams[paramStartPos];
		Standard_Real paramEnd = theParams[paramEndPos];
		Standard_Real paramInterval = paramEnd - paramStart;

		if (paramInterval == 0) {
			throw std::invalid_argument("ScalingParamsByBaseParams Error: Zero interval in params between positions " +
				std::to_string(paramStartPos) + " and " + std::to_string(paramEndPos));
		}

		Standard_Integer startIndex = (i == 0) ? paramStartPos : (paramStartPos + 1);
		// ����ÿ�������ڵĲ���
		for (size_t j = startIndex; j <= paramEndPos; ++j) {
			Standard_Real scaledParam = baseStart + (theParams[j] - paramStart) * baseInterval / paramInterval;
			scaledParams.push_back(scaledParam);
		}
	}

	return scaledParams;
}

/// <summary>
/// ����compatible  ��compatible�󽻵����ֵ��ͬ���ڵ���ͬ
/// </summary>
/// <param name="theInterCurves"></param> �ཻ������
/// <param name="theCompatibleCurves"></param> ��Ҫcompatible������
/// <param name="theTolerance"></param> ����ݲ�
/// <returns></returns>	�Ƿ�compatible�ɹ�
Standard_Boolean CurveOperate::CompatibleWithInterPoints(const std::vector<Handle(Geom_BSplineCurve)>& theInterCurves, std::vector<Handle(Geom_BSplineCurve)>& theCompatibleCurves, Standard_Real theTolerance)
{
	//1.��ȡ�����Լ����������
	std::vector<std::vector<gp_Pnt>> interPoints;
	std::vector<std::vector<Standard_Real>> interPointOrgParams;
	interPoints.reserve(theCompatibleCurves.size());
	interPointOrgParams.reserve(theCompatibleCurves.size());

	for (const auto& curve : theCompatibleCurves) 
	{
		try {
			std::vector<gp_Pnt> pointsOnTheCurve;
			std::vector<Standard_Real> paramsOnTheCurve;
			std::tie(pointsOnTheCurve, paramsOnTheCurve) = [&]() 
				{
				auto result = CalCurvesInterPointsParamsToCurve(theInterCurves, curve);
				auto& points = std::get<0>(result);
				auto& params = std::get<1>(result);

				if (points.size() != params.size()) {
					std::cerr << "Error: The number of points and parameters do not match for a curve." << std::endl;
					std::cerr << "Points size: " << points.size()
						<< ", Parameters size: " << params.size() << std::endl;
				}

				// �������
				std::vector<std::pair<Standard_Real, gp_Pnt>> combined;
				for (size_t i = 0; i < params.size(); ++i)
				{
					combined.emplace_back(params[i], points[i]);
				}
				std::sort(combined.begin(), combined.end(), [](const auto& a, const auto& b) 
					{
					return a.first < b.first;
					});
				for (size_t i = 0; i < combined.size(); ++i)
				{
					params[i] = combined[i].first;
					points[i] = combined[i].second;
				}
				return result;
				}();
			interPoints.emplace_back(std::move(pointsOnTheCurve));
			interPointOrgParams.emplace_back(std::move(paramsOnTheCurve));
		} HANDLE_EXCEPTIONS_CONTINUE
	}

	//2.��compatibleCurves�ϼ��ܲ���
	Standard_Integer denseSamplingNum = 30;
	std::vector<std::vector<gp_Pnt>> denseSamplingPoints;
	std::vector<std::vector<Standard_Real>> denseSamplingPointsParams;
	denseSamplingPoints.reserve(theCompatibleCurves.size());
	denseSamplingPointsParams.reserve(theCompatibleCurves.size());

	for (Standard_Integer i = 0; i < theCompatibleCurves.size(); i++)
	{
		try {
			std::vector<gp_Pnt> densePoints;
			std::vector<Standard_Real> denseParams;
			std::tie(densePoints, denseParams) = DenseSampling(theCompatibleCurves[i], interPointOrgParams[i], denseSamplingNum);

			if (densePoints.size() != denseParams.size())
			{
				std::cerr << "Error: The number of points and parameters do not match for a curve." << std::endl;
				std::cerr << "Points size: " << densePoints.size()
					<< ", Parameters size: " << denseParams.size() << std::endl;
				continue;
			}
			denseSamplingPoints.emplace_back(std::move(densePoints));
			denseSamplingPointsParams.emplace_back(std::move(denseParams));
		} HANDLE_EXCEPTIONS_CONTINUE
	}

	//3.�����в�����ȡƽ��
	if (interPointOrgParams.empty())
	{
		throw std::runtime_error("interPointOrgParams is empty, cannot compute average parameters.");
	}
	size_t paramSize = interPointOrgParams[0].size();
	for (const auto& paramVec : interPointOrgParams)
	{
		if (paramVec.size() != paramSize) 
		{
			throw std::runtime_error("interPointOrgParams contains vectors of differing sizes.");
		}
	}
	std::vector<Standard_Real> avgParams(interPointOrgParams[0].size(), 0.0);
	for (const auto& paramVec : interPointOrgParams) 
	{
		for (size_t i = 0; i < paramVec.size(); ++i) 
		{
			avgParams[i] += paramVec[i];
		}
	}
	for (auto& val : avgParams) {
		val /= interPointOrgParams.size();
	}

	//4.�����в������²�����
	for (size_t i = 0; i < denseSamplingPointsParams.size(); ++i)
	{
		try {
			denseSamplingPointsParams[i] = ScalingParamsByBaseParams(avgParams, denseSamplingNum + 1, denseSamplingPointsParams[i]);
		} HANDLE_EXCEPTIONS_CONTINUE
	}

	//5.�ϲ��ڵ�����
	std::vector<Standard_Real> knots;
	try {
		knots = CalSameKnotFromCurves(theCompatibleCurves, 0.05);
	} HANDLE_EXCEPTIONS_RETURN_FALSE

		//6.�ҵ�������
		Standard_Integer aDegree = 0;
	try {
		auto maxCurve = std::max_element(theCompatibleCurves.begin(), theCompatibleCurves.end(),
			[](const Handle(Geom_BSplineCurve)& curve1, const Handle(Geom_BSplineCurve)& curve2) -> bool {
				if (curve1.IsNull()) return true;
				if (curve2.IsNull()) return false;
				return curve1->Degree() < curve2->Degree();
			}
		);

		if (maxCurve != theCompatibleCurves.end() && !(*maxCurve).IsNull()) {
			aDegree = (*maxCurve)->Degree();
		}
		else {
			std::cerr << "Error: No valid curves found in compatibleCurves." << std::endl;
			return false;
		}
	} HANDLE_EXCEPTIONS_RETURN_FALSE

		//7.�������
		Standard_Integer addNum = 0;
		if (std::fabs(knots.back() - 1.) > 0.01) {
			for (Standard_Integer i = knots.size() - 1; i > 0;) {
				if (std::fabs(knots[i] - knots[i - 1]) <= 0.01) {
					knots[i] = 1.;
					addNum++;
					i--;
					if (std::fabs(knots[i] - knots[i - 1]) > 0.01) {
						break;
					}
				}
			}
			knots.insert(knots.end(), aDegree + 1 - addNum, 1.);
		}
		for (size_t i = 0; i < theCompatibleCurves.size(); ++i) {
			try {
				theCompatibleCurves[i] = ApproximateC(denseSamplingPoints[i], denseSamplingPointsParams[i], knots, aDegree);
			} HANDLE_EXCEPTIONS_CONTINUE
		}

	return true;
}
