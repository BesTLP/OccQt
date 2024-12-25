#include "GordenSurface.h"
#include "GordenSurface.h"

// 标准库头文件
#include <vector>
#include <algorithm>
#include <iostream>
#include <utility> // for std::pair

#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TopoDS_Face.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_XYZ.hxx>
#include <Precision.hxx>
#include <math_Matrix.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include "InterPolate.h"
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <SurfaceModelingTool.h>

void GeomLib_ChangeUBounds(Handle(Geom_BSplineSurface)& aSurface,
	const Standard_Real newU1,
	const Standard_Real newU2)
{
	TColStd_Array1OfReal  knots(1, aSurface->NbUKnots());
	aSurface->UKnots(knots);
	BSplCLib::Reparametrize(newU1, newU2, knots);
	aSurface->SetUKnots(knots);
}

void GeomLib_ChangeVBounds(Handle(Geom_BSplineSurface)& aSurface,
	const Standard_Real newV1,
	const Standard_Real newV2)
{
	TColStd_Array1OfReal  knots(1, aSurface->NbVKnots());
	aSurface->VKnots(knots);
	BSplCLib::Reparametrize(newV1, newV2, knots);
	aSurface->SetVKnots(knots);
}

void GeomLib_ChangeCurveBounds(Handle(Geom_BSplineCurve)& aCurve,
	const Standard_Real newU1,
	const Standard_Real newU2)
{
	TColStd_Array1OfReal  knots(1, aCurve->NbKnots());
	aCurve->Knots(knots);
	BSplCLib::Reparametrize(newU1, newU2, knots);
	aCurve->SetKnots(knots);
}


void GordenSurface::BuildMyGordonSurf(std::vector<Handle(Geom_BSplineCurve)>& uCurves,
	std::vector<Handle(Geom_BSplineCurve)>& vCurves,
	TopoDS_Face& face)
{

	//--------------- 检查传入曲线是否空 ---------------
	if (uCurves.empty())
	{
		std::cout << "U向曲线为空！" << std::endl;
		return;
	}
	if (uCurves.empty())
	{
		std::cout << "V向曲线为空！" << std::endl;
		return;
	}


	//--------------- 提升曲线次数至相同 ---------------
	int uDegree = 3, vDegree = 3;
	for (auto c : uCurves)
	{
		uDegree = std::max(uDegree, c->Degree());
	}
	for (auto c : vCurves)
	{
		vDegree = std::max(uDegree, c->Degree());
	}

	for (auto c : uCurves)
	{
		c->IncreaseDegree(uDegree);
	}
	for (auto c : vCurves)
	{
		c->IncreaseDegree(vDegree);
	}


	//--------------- 将曲线方向调整为一致 ---------------
	int usize = uCurves.size();
	int vsize = vCurves.size();
	int i = 0;

	while (i < usize - 1)
	{
		gp_Pnt p1, p2;
		gp_Vec v1, v2;
		Standard_Real midPara1 = 0.5 * (uCurves[i]->FirstParameter() + uCurves[i]->LastParameter());
		Standard_Real midPara2 = 0.5 * (uCurves[i + 1]->FirstParameter() + uCurves[i + 1]->LastParameter());
		uCurves[i]->D1(midPara1, p1, v1);
		uCurves[i + 1]->D1(midPara2, p2, v2);

		if (v1.Dot(v2) < 0)
		{
			uCurves[i + 1]->Reverse();
		}

		i++;
	}

	i = 0;
	while (i < vsize - 1)
	{
		gp_Pnt p1, p2;
		gp_Vec v1, v2;
		Standard_Real midPara1 = 0.5 * (vCurves[i]->FirstParameter() + vCurves[i]->LastParameter());
		Standard_Real midPara2 = 0.5 * (vCurves[i + 1]->FirstParameter() + vCurves[i + 1]->LastParameter());
		vCurves[i]->D1(midPara1, p1, v1);
		vCurves[i + 1]->D1(midPara2, p2, v2);

		if (v1.Dot(v2) < 0)
		{
			vCurves[i + 1]->Reverse();
		}

		i++;
	}


	//--------------- 调整两个方向的曲线起点一致 ---------------
	std::vector< std::pair<unsigned, unsigned> >
		oris = { {0u, 0u},
				 {0u, 1u},
				 {1u, 0u},
				 {1u, 1u} };

	bool   syncStop = false;
	size_t syncAttempt = 0;
	do
	{
		std::vector<Handle(Geom_BSplineCurve)> _uCurves;
		std::vector<Handle(Geom_BSplineCurve)> _vCurves;

		if (oris[syncAttempt].first)
		{
			for (const auto& C : uCurves)
			{
				_uCurves.push_back(Handle(Geom_BSplineCurve)::DownCast(C->Reversed()));
			}
		}
		else
		{
			_uCurves = uCurves;
		}

		if (oris[syncAttempt].second)
		{
			for (const auto& C : vCurves)
			{
				_vCurves.push_back(Handle(Geom_BSplineCurve)::DownCast(C->Reversed()));
			}
		}
		else
		{
			_vCurves = vCurves;
		}

		const gp_Pnt OP = _uCurves[0]->Value(_uCurves[0]->FirstParameter());
		const gp_Pnt OG = _vCurves[0]->Value(_vCurves[0]->FirstParameter());

		if (OP.Distance(OG) < 1.e-2)
		{
			syncStop = true;
			uCurves = _uCurves;
			vCurves = _vCurves;
		}

		if (++syncAttempt > 3)
		{
			syncStop = true;
		}
	} while (!syncStop);


	//--------------- 将两组曲线分别设置为相容 ---------------

	for (auto c : uCurves)
	{
		GeomLib_ChangeCurveBounds(c, 0, 1);
	}
	for (auto c : vCurves)
	{
		GeomLib_ChangeCurveBounds(c, 0, 1);
	}
	for (int i = 0; i < uCurves.size() - 1; i++)
	{
		SurfaceModelingTool::SetSameDistribution(uCurves[i], uCurves[i + 1]);
	}
	for (int i = 0; i < vCurves.size() - 1; i++)
	{
		SurfaceModelingTool::SetSameDistribution(vCurves[i], vCurves[i + 1]);
	}
	for (int i = 0; i < uCurves.size() - 1; i++)
	{
		SurfaceModelingTool::SetSameDistribution(uCurves[i], uCurves[uCurves.size() - 1]);
	}
	for (int i = 0; i < vCurves.size() - 1; i++)
	{
		SurfaceModelingTool::SetSameDistribution(vCurves[i], vCurves[vCurves.size() - 1]);
	}


	//--------------- 计算曲线的交点参数 ---------------
	math_Matrix interParaMatrixU(0, usize - 1, 0, vsize - 1);
	math_Matrix interParaMatrixV(0, usize - 1, 0, vsize - 1);
	TColgp_Array2OfPnt interPoints(1, usize, 1, vsize);
	std::vector<gp_Pnt> Pnts;
	std::vector<gp_Pnt2d> PntParams;

	for (int i = 0; i < usize; i++)
	{
		for (int j = 0; j < vsize; j++)
		{
			// 考虑三边的情况
			// u向退化边
			if (uCurves[i]->StartPoint().IsEqual(uCurves[i]->EndPoint(), 1.e-2))
			{
				interPoints(i + 1, j + 1) = uCurves[i]->StartPoint();
				Pnts.push_back(uCurves[i]->StartPoint());
				gp_Pnt2d Pnt2d((j + 1.0) / vsize, 1);
				PntParams.push_back(Pnt2d);
				continue;
			}
			// v向退化边
			else if (vCurves[j]->StartPoint().IsEqual(vCurves[j]->EndPoint(), 1.e-2))
			{
				interPoints(i + 1, j + 1) = vCurves[j]->StartPoint();
				Pnts.push_back(vCurves[j]->StartPoint());
				gp_Pnt2d Pnt2d(1, (i + 1.0) / usize);
				PntParams.push_back(Pnt2d);
				continue;
			}

			GeomAPI_ExtremaCurveCurve extrema(uCurves[i], vCurves[j]);

			if (extrema.NbExtrema() != 1)
			{
				std::cout << "最近点对不唯一！" << std::endl;
			}
			int nbEx = extrema.NbExtrema();
			//最近点对如果出问题，可从内部线转等参线的工作中取

			Standard_Real para1, para2;
			extrema.Parameters(1, para1, para2);
			interParaMatrixU(i, j) = para1;
			interParaMatrixV(i, j) = para2;

			gp_Pnt p1, p2;
			extrema.NearestPoints(p1, p2);
			gp_Pnt interPnt((p1.XYZ() + p2.XYZ()) / 2.0);
			interPoints(i + 1, j + 1) = interPnt;

			Pnts.push_back(interPnt);
			gp_Pnt2d Pnt2d(para1, para2);
			PntParams.push_back(Pnt2d);
		}
	}

	if (interPoints.Size() != usize * vsize)
	{
		std::cout << "获取的等参线最近点对数量有误！" << std::endl;
	}



	//--------------- 构造三张曲面 ---------------
	Handle(Geom_BSplineSurface) L1, L2, T;

	int uDegree1 = usize <= 3 ? 1 : 3;
	int vDegree1 = vsize <= 3 ? 1 : 3;

	// chenxin's loft
	L1 = InterPolateTool::Loft(uCurves, uDegree1);
	L2 = InterPolateTool::LoftV(vCurves, vDegree1);

	if (L1.IsNull() || L2.IsNull())
	{
		std::cout << "放样失败！" << std::endl;
		return;
	}

	TColStd_Array1OfReal uKnotsTCol = L2->UKnots();
	TColStd_Array1OfInteger uMultsTCol = L2->UMultiplicities();
	TColStd_Array1OfReal vKnotsTCol = L1->VKnots();
	TColStd_Array1OfInteger vMultsTCol = L1->VMultiplicities();
	std::vector<double> uKnots;
	std::vector<double> vKnots;
	std::vector<int> uMults;
	std::vector<int> vMults;
	for (int i = 1; i <= uKnotsTCol.Size(); i++)
	{
		uKnots.push_back(uKnotsTCol(i));
	}
	for (int i = 1; i <= vKnotsTCol.Size(); i++)
	{
		vKnots.push_back(vKnotsTCol(i));
	}
	for (int i = 1; i <= uMultsTCol.Size(); i++)
	{
		uMults.push_back(uMultsTCol(i));
	}
	for (int i = 1; i <= vMultsTCol.Size(); i++)
	{
		vMults.push_back(vMultsTCol(i));
	}

	T = InterPolateTool::Interpolate(Pnts, PntParams, uKnots, vKnots, uMults, vMults, vDegree1, uDegree1);

	L1->IncreaseDegree(uCurves[0]->Degree(), vCurves[0]->Degree());
	L2->IncreaseDegree(uCurves[0]->Degree(), vCurves[0]->Degree());
	T->IncreaseDegree(uCurves[0]->Degree(), vCurves[0]->Degree());


	//--------------- 三张曲面compatible ---------------
	//将三张曲面的参数域都scale到[0,1]
	GeomLib_ChangeUBounds(L1, 0, 1);
	GeomLib_ChangeVBounds(L1, 0, 1);
	GeomLib_ChangeUBounds(L2, 0, 1);
	GeomLib_ChangeVBounds(L2, 0, 1);
	GeomLib_ChangeUBounds(T, 0, 1);
	GeomLib_ChangeVBounds(T, 0, 1);

	// Get the u knot vector
	Standard_Integer NbUKnot1 = L1->NbUKnots();
	TColStd_Array1OfReal    UKnots1(1, NbUKnot1);
	TColStd_Array1OfInteger UMults1(1, NbUKnot1);
	L1->UKnots(UKnots1);
	L1->UMultiplicities(UMults1);
	// Get the v knot vector
	Standard_Integer NbVKnot1 = L1->NbVKnots();
	TColStd_Array1OfReal    VKnots1(1, NbVKnot1);
	TColStd_Array1OfInteger VMults1(1, NbVKnot1);
	L1->VKnots(VKnots1);
	L1->VMultiplicities(VMults1);

	for (int i = 1; i <= NbUKnot1; i++)
	{
		L2->InsertUKnot(UKnots1(i), UMults1(i), 1.e-15, false);
	}
	for (int i = 1; i <= NbVKnot1; i++)
	{
		L2->InsertVKnot(VKnots1(i), VMults1(i), 1.e-15, false);
	}

	// Get the u knot vector
	Standard_Integer NbUKnot2 = L2->NbUKnots();
	TColStd_Array1OfReal    UKnots2(1, NbUKnot2);
	TColStd_Array1OfInteger UMults2(1, NbUKnot2);
	L2->UKnots(UKnots2);
	L2->UMultiplicities(UMults2);
	// Get the v knot vector
	Standard_Integer NbVKnot2 = L2->NbVKnots();
	TColStd_Array1OfReal    VKnots2(1, NbVKnot2);
	TColStd_Array1OfInteger VMults2(1, NbVKnot2);
	L2->VKnots(VKnots2);
	L2->VMultiplicities(VMults2);

	for (int i = 1; i <= NbUKnot2; i++)
	{
		T->InsertUKnot(UKnots2(i), UMults2(i), 1.e-15, false);
	}
	for (int i = 1; i <= NbVKnot2; i++)
	{
		T->InsertVKnot(VKnots2(i), VMults2(i), 1.e-15, false);
	}

	// Get the u knot vector
	Standard_Integer NbUKnot3 = T->NbUKnots();
	TColStd_Array1OfReal    UKnots3(1, NbUKnot3);
	TColStd_Array1OfInteger UMults3(1, NbUKnot3);
	T->UKnots(UKnots3);
	T->UMultiplicities(UMults3);
	// Get the v knot vector
	Standard_Integer NbVKnot3 = T->NbVKnots();
	TColStd_Array1OfReal    VKnots3(1, NbVKnot3);
	TColStd_Array1OfInteger VMults3(1, NbVKnot3);
	T->VKnots(VKnots3);
	T->VMultiplicities(VMults3);

	for (int i = 1; i <= NbUKnot3; i++)
	{
		L1->InsertUKnot(UKnots3(i), UMults3(i), 1.e-15, false);
		L2->InsertUKnot(UKnots3(i), UMults3(i), 1.e-15, false);
	}
	for (int i = 1; i <= NbVKnot3; i++)
	{
		L1->InsertVKnot(VKnots3(i), VMults3(i), 1.e-15, false);
		L2->InsertVKnot(VKnots3(i), VMults3(i), 1.e-15, false);
	}
	//至此，三张曲面的节点向量完全相同


	//--------------- 得到共同节点向量和次数 ---------------
	const TColStd_Array1OfReal knotsU = L1->UKnots();
	const TColStd_Array1OfReal knotsV = L1->VKnots();
	const TColStd_Array1OfReal knotsU2 = L2->UKnots();
	const TColStd_Array1OfReal knotsV2 = L2->VKnots();
	const TColStd_Array1OfReal knotsUT = T->UKnots();
	const TColStd_Array1OfReal knotsVT = T->VKnots();
	const TColStd_Array1OfInteger multsU = L1->UMultiplicities();
	const TColStd_Array1OfInteger multsV = L1->VMultiplicities();
	const TColStd_Array1OfInteger multsU2 = L2->UMultiplicities();
	const TColStd_Array1OfInteger multsV2 = L2->VMultiplicities();
	const TColStd_Array1OfInteger multsU3 = T->UMultiplicities();
	const TColStd_Array1OfInteger multsV3 = T->VMultiplicities();
	const int degreeU = L1->UDegree();
	const int degreeV = L1->VDegree();


	//--------------- 计算控制点 ---------------
	const TColgp_Array2OfPnt poles1 = L1->Poles();
	const TColgp_Array2OfPnt poles2 = L2->Poles();
	const TColgp_Array2OfPnt poles12 = T->Poles();

	const int nbPolesU = L1->NbUPoles();
	const int nbPolesV = L1->NbVPoles();
	const int nbPolesU2 = L2->NbUPoles();
	const int nbPolesV2 = L2->NbVPoles();
	const int nbPolesU12 = T->NbUPoles();
	const int nbPolesV12 = T->NbVPoles();

	TColgp_Array2OfPnt resPole(1, nbPolesU, 1, nbPolesV);

	for (int i = 1; i <= nbPolesU12; i++)
	{
		for (int j = 1; j <= nbPolesV12; j++)
		{
			gp_XYZ coord = poles1(i, j).Coord() + poles2(i, j).Coord() - poles12(i, j).Coord();
			resPole(i, j).SetCoord(coord.X(), coord.Y(), coord.Z());
		}
	}


	//--------------- 构造Gordon曲面 ---------------
	Handle(Geom_Surface) gordon = new Geom_BSplineSurface(resPole, knotsU,
		knotsV, multsU, multsV, degreeU, degreeV);

	face = BRepBuilderAPI_MakeFace(gordon, Precision::Confusion());
}

void GordenSurface::BuildMyGordonSurf(std::vector<Handle(Geom_BSplineCurve)>& uCurves,
	std::vector<Handle(Geom_BSplineCurve)>& vCurves,
	const std::vector<Standard_Real>& theUParams,
	const std::vector<Standard_Real>& theVParams,
	TopoDS_Face& face) {

	//--------------- 检查传入曲线是否空 ---------------
	if (uCurves.empty()) {
		std::cout << "U向曲线为空！" << std::endl;
		return;
	}
	if (uCurves.empty()) {
		std::cout << "V向曲线为空！" << std::endl;
		return;
	}


	//--------------- 提升曲线次数至相同 ---------------
	int uDegree = 3, vDegree = 3;
	for (auto c : uCurves) {
		uDegree = std::max(uDegree, c->Degree());
	}
	for (auto c : vCurves) {
		vDegree = std::max(uDegree, c->Degree());
	}

	for (auto c : uCurves) {
		c->IncreaseDegree(uDegree);
	}
	for (auto c : vCurves) {
		c->IncreaseDegree(vDegree);
	}


	//--------------- 将曲线方向调整为一致 ---------------
	int usize = uCurves.size();
	int vsize = vCurves.size();
	int i = 0;

	while (i < usize - 1) {
		gp_Pnt p1, p2;
		gp_Vec v1, v2;
		Standard_Real midPara1 = 0.5 * (uCurves[i]->FirstParameter() + uCurves[i]->LastParameter());
		Standard_Real midPara2 = 0.5 * (uCurves[i + 1]->FirstParameter() + uCurves[i + 1]->LastParameter());
		uCurves[i]->D1(midPara1, p1, v1);
		uCurves[i + 1]->D1(midPara2, p2, v2);

		if (v1.Dot(v2) < 0) {
			uCurves[i + 1]->Reverse();
		}

		i++;
	}

	i = 0;
	while (i < vsize - 1) {
		gp_Pnt p1, p2;
		gp_Vec v1, v2;
		Standard_Real midPara1 = 0.5 * (vCurves[i]->FirstParameter() + vCurves[i]->LastParameter());
		Standard_Real midPara2 = 0.5 * (vCurves[i + 1]->FirstParameter() + vCurves[i + 1]->LastParameter());
		vCurves[i]->D1(midPara1, p1, v1);
		vCurves[i + 1]->D1(midPara2, p2, v2);

		if (v1.Dot(v2) < 0) {
			vCurves[i + 1]->Reverse();
		}

		i++;
	}


	//--------------- 调整两个方向的曲线起点一致 ---------------
	std::vector< std::pair<unsigned, unsigned> >
		oris = { {0u, 0u},
				 {0u, 1u},
				 {1u, 0u},
				 {1u, 1u} };

	bool   syncStop = false;
	size_t syncAttempt = 0;
	do {
		std::vector<Handle(Geom_BSplineCurve)> _uCurves;
		std::vector<Handle(Geom_BSplineCurve)> _vCurves;

		if (oris[syncAttempt].first) {
			for (const auto& C : uCurves) {
				_uCurves.push_back(Handle(Geom_BSplineCurve)::DownCast(C->Reversed()));
			}
		} else {
			_uCurves = uCurves;
		}

		if (oris[syncAttempt].second) {
			for (const auto& C : vCurves) {
				_vCurves.push_back(Handle(Geom_BSplineCurve)::DownCast(C->Reversed()));
			}
		} else {
			_vCurves = vCurves;
		}

		const gp_Pnt OP = _uCurves[0]->Value(_uCurves[0]->FirstParameter());
		const gp_Pnt OG = _vCurves[0]->Value(_vCurves[0]->FirstParameter());

		if (OP.Distance(OG) < 1.e-2) {
			syncStop = true;
			uCurves = _uCurves;
			vCurves = _vCurves;
		}

		if (++syncAttempt > 3) {
			syncStop = true;
		}
	} while (!syncStop);


	//--------------- 将两组曲线分别设置为相容 ---------------

	for (auto c : uCurves) {
		GeomLib_ChangeCurveBounds(c, 0, 1);
	}
	for (auto c : vCurves) {
		GeomLib_ChangeCurveBounds(c, 0, 1);
	}
	for (int i = 0; i < uCurves.size() - 1; i++) {
		SurfaceModelingTool::SetSameDistribution(uCurves[i], uCurves[i + 1]);
	}
	for (int i = 0; i < vCurves.size() - 1; i++) {
		SurfaceModelingTool::SetSameDistribution(vCurves[i], vCurves[i + 1]);
	}
	for (int i = 0; i < uCurves.size() - 1; i++) {
		SurfaceModelingTool::SetSameDistribution(uCurves[i], uCurves[uCurves.size() - 1]);
	}
	for (int i = 0; i < vCurves.size() - 1; i++) {
		SurfaceModelingTool::SetSameDistribution(vCurves[i], vCurves[vCurves.size() - 1]);
	}


	//--------------- 计算曲线的交点参数 ---------------
	math_Matrix interParaMatrixU(0, usize - 1, 0, vsize - 1);
	math_Matrix interParaMatrixV(0, usize - 1, 0, vsize - 1);
	TColgp_Array2OfPnt interPoints(1, usize, 1, vsize);
	std::vector<gp_Pnt> Pnts;
	std::vector<gp_Pnt2d> PntParams;

	for (int i = 0; i < usize; i++) {
		for (int j = 0; j < vsize; j++) {
			// 考虑三边的情况
			// u向退化边
			if (uCurves[i]->StartPoint().IsEqual(uCurves[i]->EndPoint(), 1.e-2)) {
				interPoints(i + 1, j + 1) = uCurves[i]->StartPoint();
				Pnts.push_back(uCurves[i]->StartPoint());
				gp_Pnt2d Pnt2d(theVParams[j], theUParams[i]);
				PntParams.push_back(Pnt2d);
				continue;
			}
			// v向退化边
			else if (vCurves[j]->StartPoint().IsEqual(vCurves[j]->EndPoint(), 1.e-2)) {
				interPoints(i + 1, j + 1) = vCurves[j]->StartPoint();
				Pnts.push_back(vCurves[j]->StartPoint());
				gp_Pnt2d Pnt2d(theVParams[j], theUParams[i]);
				PntParams.push_back(Pnt2d);
				continue;
			}

			GeomAPI_ExtremaCurveCurve extrema(uCurves[i], vCurves[j]);

			if (extrema.NbExtrema() != 1) {
				std::cout << "最近点对不唯一！" << std::endl;
			}
			int nbEx = extrema.NbExtrema();

			gp_Pnt p1, p2;
			extrema.NearestPoints(p1, p2);
			gp_Pnt interPnt((p1.XYZ() + p2.XYZ()) / 2.0);
			interPoints(i + 1, j + 1) = interPnt;

			Pnts.push_back(interPnt);
			gp_Pnt2d Pnt2d(theVParams[j], theUParams[i]);
			PntParams.push_back(Pnt2d);
		}
	}

	if (interPoints.Size() != usize * vsize) {
		std::cout << "获取的等参线最近点对数量有误！" << std::endl;
	}



	//--------------- 构造三张曲面 ---------------
	Handle(Geom_BSplineSurface) L1, L2, T;

	int uDegree1 = usize <= 3 ? 1 : 3;
	int vDegree1 = vsize <= 3 ? 1 : 3;

	// chenxin's loft
	L1 = InterPolateTool::Loft(uCurves, theUParams, uDegree1);
	L2 = InterPolateTool::LoftV(vCurves, theVParams, vDegree1);

	if (L1.IsNull() || L2.IsNull()) {
		std::cout << "放样失败！" << std::endl;
		return;
	}

	TColStd_Array1OfReal uKnotsTCol = L2->UKnots();
	TColStd_Array1OfInteger uMultsTCol = L2->UMultiplicities();
	TColStd_Array1OfReal vKnotsTCol = L1->VKnots();
	TColStd_Array1OfInteger vMultsTCol = L1->VMultiplicities();
	std::vector<double> uKnots;
	std::vector<double> vKnots;
	std::vector<int> uMults;
	std::vector<int> vMults;
	for (int i = 1; i <= uKnotsTCol.Size(); i++) {
		uKnots.push_back(uKnotsTCol(i));
	}
	for (int i = 1; i <= vKnotsTCol.Size(); i++) {
		vKnots.push_back(vKnotsTCol(i));
	}
	for (int i = 1; i <= uMultsTCol.Size(); i++) {
		uMults.push_back(uMultsTCol(i));
	}
	for (int i = 1; i <= vMultsTCol.Size(); i++) {
		vMults.push_back(vMultsTCol(i));
	}

	T = InterPolateTool::Interpolate(Pnts, PntParams, uKnots, vKnots, uMults, vMults, vDegree1, uDegree1);

	L1->IncreaseDegree(uCurves[0]->Degree(), vCurves[0]->Degree());
	L2->IncreaseDegree(uCurves[0]->Degree(), vCurves[0]->Degree());
	T->IncreaseDegree(uCurves[0]->Degree(), vCurves[0]->Degree());


	//--------------- 三张曲面compatible ---------------
	//将三张曲面的参数域都scale到[0,1]
	GeomLib_ChangeUBounds(L1, 0, 1);
	GeomLib_ChangeVBounds(L1, 0, 1);
	GeomLib_ChangeUBounds(L2, 0, 1);
	GeomLib_ChangeVBounds(L2, 0, 1);
	GeomLib_ChangeUBounds(T, 0, 1);
	GeomLib_ChangeVBounds(T, 0, 1);

	// Get the u knot vector
	Standard_Integer NbUKnot1 = L1->NbUKnots();
	TColStd_Array1OfReal    UKnots1(1, NbUKnot1);
	TColStd_Array1OfInteger UMults1(1, NbUKnot1);
	L1->UKnots(UKnots1);
	L1->UMultiplicities(UMults1);
	// Get the v knot vector
	Standard_Integer NbVKnot1 = L1->NbVKnots();
	TColStd_Array1OfReal    VKnots1(1, NbVKnot1);
	TColStd_Array1OfInteger VMults1(1, NbVKnot1);
	L1->VKnots(VKnots1);
	L1->VMultiplicities(VMults1);

	for (int i = 1; i <= NbUKnot1; i++) {
		L2->InsertUKnot(UKnots1(i), UMults1(i), 1.e-15, false);
	}
	for (int i = 1; i <= NbVKnot1; i++) {
		L2->InsertVKnot(VKnots1(i), VMults1(i), 1.e-15, false);
	}

	// Get the u knot vector
	Standard_Integer NbUKnot2 = L2->NbUKnots();
	TColStd_Array1OfReal    UKnots2(1, NbUKnot2);
	TColStd_Array1OfInteger UMults2(1, NbUKnot2);
	L2->UKnots(UKnots2);
	L2->UMultiplicities(UMults2);
	// Get the v knot vector
	Standard_Integer NbVKnot2 = L2->NbVKnots();
	TColStd_Array1OfReal    VKnots2(1, NbVKnot2);
	TColStd_Array1OfInteger VMults2(1, NbVKnot2);
	L2->VKnots(VKnots2);
	L2->VMultiplicities(VMults2);

	for (int i = 1; i <= NbUKnot2; i++) {
		T->InsertUKnot(UKnots2(i), UMults2(i), 1.e-15, false);
	}
	for (int i = 1; i <= NbVKnot2; i++) {
		T->InsertVKnot(VKnots2(i), VMults2(i), 1.e-15, false);
	}

	// Get the u knot vector
	Standard_Integer NbUKnot3 = T->NbUKnots();
	TColStd_Array1OfReal    UKnots3(1, NbUKnot3);
	TColStd_Array1OfInteger UMults3(1, NbUKnot3);
	T->UKnots(UKnots3);
	T->UMultiplicities(UMults3);
	// Get the v knot vector
	Standard_Integer NbVKnot3 = T->NbVKnots();
	TColStd_Array1OfReal    VKnots3(1, NbVKnot3);
	TColStd_Array1OfInteger VMults3(1, NbVKnot3);
	T->VKnots(VKnots3);
	T->VMultiplicities(VMults3);

	for (int i = 1; i <= NbUKnot3; i++) {
		L1->InsertUKnot(UKnots3(i), UMults3(i), 1.e-15, false);
		L2->InsertUKnot(UKnots3(i), UMults3(i), 1.e-15, false);
	}
	for (int i = 1; i <= NbVKnot3; i++) {
		L1->InsertVKnot(VKnots3(i), VMults3(i), 1.e-15, false);
		L2->InsertVKnot(VKnots3(i), VMults3(i), 1.e-15, false);
	}
	//至此，三张曲面的节点向量完全相同


	//--------------- 得到共同节点向量和次数 ---------------
	const TColStd_Array1OfReal knotsU = L1->UKnots();
	const TColStd_Array1OfReal knotsV = L1->VKnots();
	const TColStd_Array1OfReal knotsU2 = L2->UKnots();
	const TColStd_Array1OfReal knotsV2 = L2->VKnots();
	const TColStd_Array1OfReal knotsUT = T->UKnots();
	const TColStd_Array1OfReal knotsVT = T->VKnots();
	const TColStd_Array1OfInteger multsU = L1->UMultiplicities();
	const TColStd_Array1OfInteger multsV = L1->VMultiplicities();
	const TColStd_Array1OfInteger multsU2 = L2->UMultiplicities();
	const TColStd_Array1OfInteger multsV2 = L2->VMultiplicities();
	const TColStd_Array1OfInteger multsU3 = T->UMultiplicities();
	const TColStd_Array1OfInteger multsV3 = T->VMultiplicities();
	const int degreeU = L1->UDegree();
	const int degreeV = L1->VDegree();


	//--------------- 计算控制点 ---------------
	const TColgp_Array2OfPnt poles1 = L1->Poles();
	const TColgp_Array2OfPnt poles2 = L2->Poles();
	const TColgp_Array2OfPnt poles12 = T->Poles();

	const int nbPolesU = L1->NbUPoles();
	const int nbPolesV = L1->NbVPoles();
	const int nbPolesU2 = L2->NbUPoles();
	const int nbPolesV2 = L2->NbVPoles();
	const int nbPolesU12 = T->NbUPoles();
	const int nbPolesV12 = T->NbVPoles();

	TColgp_Array2OfPnt resPole(1, nbPolesU, 1, nbPolesV);

	for (int i = 1; i <= nbPolesU12; i++) {
		for (int j = 1; j <= nbPolesV12; j++) {
			gp_XYZ coord = poles1(i, j).Coord() + poles2(i, j).Coord() - poles12(i, j).Coord();
			resPole(i, j).SetCoord(coord.X(), coord.Y(), coord.Z());
		}
	}


	//--------------- 构造Gordon曲面 ---------------
	Handle(Geom_Surface) gordon = new Geom_BSplineSurface(resPole, knotsU,
		knotsV, multsU, multsV, degreeU, degreeV);

	face = BRepBuilderAPI_MakeFace(gordon, Precision::Confusion());
}

