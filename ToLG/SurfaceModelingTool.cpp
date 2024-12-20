Standard_Integer SurfaceModelingTool::SetSameDistribution(Handle(
	Geom_BSplineCurve)& C1,
	Handle(Geom_BSplineCurve)& C2)
{

	Standard_Integer C1_Degree = C1->Degree();
	Standard_Integer C2_Degree = C2->Degree();

	// 如果Degree不相同的话，对低阶进行升阶
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
			// 三边
			isThreeBoundary = true;
			break;
		}
	}
	if (isThreeBoundary)
	{
		GeomAPI_ExtremaCurveCurve extrema1(uCurve, aBoundarycurveArray[0]);
		GeomAPI_ExtremaCurveCurve extrema2(uCurve, aBoundarycurveArray[2]);

		gp_Pnt p1, p2, p3;
		extrema1.NearestPoints(p1, p2);
		extrema2.NearestPoints(p1, p3);

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

	// 如果最后一条曲线比第一条曲线距离原点更近，反转 uISOcurvesArray_Final
	if (minDistanceToPoint(uISOcurvesArray_Final[uISOcurvesArray_Final.size() - 1]) < minDistanceToPoint(uISOcurvesArray_Final[0]))
	{
		std::reverse(uISOcurvesArray_Final.begin(), uISOcurvesArray_Final.end());
	}
	// 对 vISO 曲线进行同样的处理
	if (minDistanceToPoint(vISOcurvesArray_Final[vISOcurvesArray_Final.size() - 1]) < minDistanceToPoint(vISOcurvesArray_Final[0]))
	{
		std::reverse(vISOcurvesArray_Final.begin(), vISOcurvesArray_Final.end());
	}
	MathTool::ReverseIfNeeded(uISOcurvesArray_Final);
	MathTool::ReverseIfNeeded(vISOcurvesArray_Final);
}
bool SurfaceModelingTool::GetInternalCurves(
	std::vector<Handle(Geom_BSplineCurve)>& aBoundarycurveArray,
	std::vector<Handle(Geom_BSplineCurve)>& anInternalBSplineCurves,
	std::vector<Handle(Geom_BSplineCurve)>& uInternalCurve,
	std::vector<Handle(Geom_BSplineCurve)>& vInternalCurve,
	double& uAngleSum,
	double& vAngleSum,
	double AngleTolerance)
{
	Handle(Geom_BSplineCurve) bslpineCurve1 = aBoundarycurveArray[0];
	Handle(Geom_BSplineCurve) bslpineCurve2 = aBoundarycurveArray[1];
	Handle(Geom_BSplineCurve) bslpineCurve3 = aBoundarycurveArray[2];
	Handle(Geom_BSplineCurve) bslpineCurve4 = aBoundarycurveArray[3];

	// 步骤1：初始化PlanarCurveArray，包含所有边界曲线
	std::vector<PlanarCurve> PlanarCurveArray;
	for (int i = 0; i < aBoundarycurveArray.size(); i++)
	{
		PlanarCurveArray.emplace_back(PlanarCurve(aBoundarycurveArray[i]));
	}

	// 步骤2：检查所有边界曲线是否都是平面曲线
	bool canUseInternalLines = true;
	for (const PlanarCurve& curve : PlanarCurveArray)
	{
		if (curve.GetCurveType() == CurveType::NOTPLANAR)
		{
			canUseInternalLines = false; // 如果有非平面曲线，设置为 false
			break;
		}
	}

	if (!canUseInternalLines)
	{
		return false; // 如果任何一条边界曲线不是平面曲线，则直接返回false
	}

	// 清空uInternalCurve和vInternalCurve，准备存储结果
	uInternalCurve.clear();
	vInternalCurve.clear();

	// 遍历内部BSpline曲线
	for (auto& internalCurve : anInternalBSplineCurves)
	{
		PlanarCurve InternalPlanarCurve(internalCurve);

		// 如果内部曲线不是平面曲线，则跳过
		if (InternalPlanarCurve.GetCurveType() == CurveType::NOTPLANAR)
		{
			continue;
		}

		// 计算内部曲线和四条边界曲线之间的距离
		double distance1 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve1);
		double distance2 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve2);
		double distance3 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve3);
		double distance4 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve4);

		double SplitPointParameters[2] = { 0 };

		// 步骤5：检查曲线是否靠近边界曲线，如果靠近，计算分割点
		if ((distance1 < 10 && distance3 < 10) || (distance2 < 10 && distance4 < 10))
		{
			GeomAPI_ExtremaCurveCurve extrema1(InternalPlanarCurve.GetCurve(), distance1 < 10 ? bslpineCurve1 : bslpineCurve2);
			GeomAPI_ExtremaCurveCurve extrema2(InternalPlanarCurve.GetCurve(), distance3 < 10 ? bslpineCurve3 : bslpineCurve4);
			gp_Pnt internalPnt;
			gp_Pnt replacePnt1, replacePnt2;
			// 获取分割点参数
			if (extrema1.NbExtrema() > 0)
			{
				double U;
				extrema1.LowerDistanceParameters(SplitPointParameters[0], U);
				extrema1.NearestPoints(internalPnt, replacePnt1);
			}
			if (extrema2.NbExtrema() > 0)
			{
				double U;
				extrema2.LowerDistanceParameters(SplitPointParameters[1], U);
				extrema2.NearestPoints(internalPnt, replacePnt2);
			}

			// 确保分割点参数正确排序
			if (SplitPointParameters[0] > SplitPointParameters[1])
			{
				std::swap(SplitPointParameters[1], SplitPointParameters[0]);
				std::swap(replacePnt1, replacePnt2);
			}

			// 对内部曲线进行裁剪
			Handle(Geom_TrimmedCurve) trimmedCurve = new Geom_TrimmedCurve(InternalPlanarCurve.GetCurve(), SplitPointParameters[0], SplitPointParameters[1]);
			Handle(Geom_BSplineCurve) aBsplineCurve = GeomConvert::CurveToBSplineCurve(trimmedCurve, Convert_TgtThetaOver2);
			//UniformCurve(aBsplineCurve);
			aBsplineCurve->SetPole(1, replacePnt1);
			aBsplineCurve->SetPole(aBsplineCurve->NbPoles(), replacePnt2);
			InternalPlanarCurve.SetCurve(aBsplineCurve);

			//distance1 = aBsplineCurve->StartPoint().Distance(replacePnt1);
			//distance2 = aBsplineCurve->EndPoint().Distance(replacePnt2);
			//std::vector<gp_Pnt> ApproximatePoints = MathTool::GetSamplePointsOnCurve(aBsplineCurve, 100);
			//ApproximatePoints[0] = replacePnt1;
			//MathTool::SortPoints(ApproximatePoints, ApproximatePoints[0]);
			//ApproximatePoints[ApproximatePoints.size() - 1] = replacePnt2;
			//std::vector<double> params = ComputeUniformParam(ApproximatePoints.size(), 0., 1.);
			//std::vector<double> tempKnots = KnotGernerationByParams(params, 3, InternalPlanarCurve.GetCurve()->Degree());
			//std::vector<double> insertKnots;
			//Handle(Geom_BSplineCurve) aBSplineCurve = IterateApproximate(insertKnots, ApproximatePoints, params, tempKnots, InternalPlanarCurve.GetCurve()->Degree(), 50, 0.01);
			//InternalPlanarCurve.SetCurve(aBSplineCurve);

			// 计算角度并分类曲线
			double angle1, angle3, angle2, angle4;
			angle1 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[0], InternalPlanarCurve);
			angle2 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[1], InternalPlanarCurve);
			angle3 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[2], InternalPlanarCurve);
			angle4 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[3], InternalPlanarCurve);

			// 根据角度值判断该曲线属于u方向还是v方向
			if (std::abs(angle1) < AngleTolerance && std::abs(angle3) < AngleTolerance)
			{
				uAngleSum += (std::abs(angle1) + std::abs(angle3)) / 2;
				uInternalCurve.push_back(InternalPlanarCurve.GetCurve());
			}
			else if (std::abs(angle2) < AngleTolerance && std::abs(angle4) < AngleTolerance)
			{
				vAngleSum += (std::abs(angle2) + std::abs(angle4)) / 2;
				vInternalCurve.push_back(InternalPlanarCurve.GetCurve());
			}
		}
	}

	// 将边界曲线添加到uInternalCurve和vInternalCurve的头部和尾部
	uInternalCurve.insert(uInternalCurve.begin(), bslpineCurve1);
	uInternalCurve.insert(uInternalCurve.end(), bslpineCurve3);
	vInternalCurve.insert(vInternalCurve.begin(), bslpineCurve2);
	vInternalCurve.insert(vInternalCurve.end(), bslpineCurve4);

	// 排序曲线并检查自交
	MathTool::SortBSplineCurves(uInternalCurve, uInternalCurve[0]);
	MathTool::SortBSplineCurves(vInternalCurve, vInternalCurve[0]);

	MathTool::CheckSelfIntersect(uInternalCurve);
	MathTool::CheckSelfIntersect(vInternalCurve);

	//SurfaceModelingTool tool;
	//tool.CompatibleWithInterPoints(uInternalCurve, vInternalCurve);
	//tool.CompatibleWithInterPoints(vInternalCurve, uInternalCurve);
	// 如果uInternalCurve和vInternalCurve的曲线数有一个大于4，则返回true
	return uInternalCurve.size() > 4 || vInternalCurve.size() > 4;
}

Handle(Geom_BSplineSurface) SurfaceModelingTool::GenerateReferSurface(
	std::vector<Handle(Geom_BSplineCurve)> aBoundarycurveArray,
	std::vector<Handle(Geom_BSplineCurve)>& uInternalCurve,
	std::vector<Handle(Geom_BSplineCurve)>& vInternalCurve,
	double uAngleSum,
	double vAngleSum,
	int isoCount,
	ReferSurfaceType referSurfaceType)
{
	if (referSurfaceType == ReferSurfaceType::GORDEN_ONE_DIRECTION_GORDEN)
	{
		// 获取输入的边界曲线
		Handle(Geom_BSplineCurve) bslpineCurve1 = aBoundarycurveArray[0];
		Handle(Geom_BSplineCurve) bslpineCurve2 = aBoundarycurveArray[1];
		Handle(Geom_BSplineCurve) bslpineCurve3 = aBoundarycurveArray[2];
		Handle(Geom_BSplineCurve) bslpineCurve4 = aBoundarycurveArray[3];

		// 存储生成的Gorden等参线曲线和剩余曲线
		std::vector<Handle(Geom_BSplineCurve)> GordenISOCurves;
		std::vector<Handle(Geom_BSplineCurve)> remainCurves;

		// 使用新算法的标志
		bool useNewAlgorithm = true;

		// 判断内部曲线的数量来选择构造Gorden曲面的方式
		if (uInternalCurve.size() > vInternalCurve.size() && uInternalCurve.size() >= 4)
		{
			// 选择u方向的内部线和边界来构造Gorden曲面
			GordenISOCurves.insert(GordenISOCurves.end(), uInternalCurve.begin(), uInternalCurve.end());
			remainCurves.push_back(bslpineCurve2);
			remainCurves.push_back(bslpineCurve4);
		}
		else if (vInternalCurve.size() > uInternalCurve.size() && vInternalCurve.size() >= 4)
		{
			// 选择v方向的内部线和边界来构造Gorden曲面
			GordenISOCurves.insert(GordenISOCurves.end(), vInternalCurve.begin(), vInternalCurve.end());
			remainCurves.push_back(bslpineCurve1);
			remainCurves.push_back(bslpineCurve3);
		}
		else if (uInternalCurve.size() == vInternalCurve.size() && uInternalCurve.size() >= 4)
		{
			// 如果u方向和v方向的内部曲线数量相等，根据角度之和来选择
			if (uAngleSum < vAngleSum)
			{
				GordenISOCurves.insert(GordenISOCurves.end(), uInternalCurve.begin(), uInternalCurve.end());
				remainCurves.push_back(bslpineCurve2);
				remainCurves.push_back(bslpineCurve4);
			}
			else
			{
				GordenISOCurves.insert(GordenISOCurves.end(), vInternalCurve.begin(), vInternalCurve.end());
				remainCurves.push_back(bslpineCurve1);
				remainCurves.push_back(bslpineCurve3);
			}
		}
		else
		{
			// 如果条件不满足，回退到现有算法
			useNewAlgorithm = false;
			return nullptr;
		}

		// 存储生成的等参线曲线和法线
		std::vector<Handle(Geom_BSplineCurve)> uCreateGordenCurves, vCreateGordenCurves;
		std::vector<gp_Vec> normalsOfUISOLines, normalsOfVISOLines;

		// 最终生成的参考曲面
		Handle(Geom_BSplineSurface) referSurface;

		if (useNewAlgorithm)
		{
			// 计算曲线与曲线之间的角度
			double AngleUwithG = MathTool::ComputeAngleBetweenCurves(bslpineCurve1, GordenISOCurves[0]);
			double AngleVwithG = MathTool::ComputeAngleBetweenCurves(bslpineCurve2, GordenISOCurves[0]);

			// 根据角度选择曲线
			if (AngleUwithG > AngleVwithG)
			{
				// 如果u方向的角度更大，调整u和v方向的曲线顺序
				vCreateGordenCurves.clear();
				vCreateGordenCurves.insert(vCreateGordenCurves.begin(), GordenISOCurves.begin(), GordenISOCurves.end());
				uCreateGordenCurves.insert(uCreateGordenCurves.begin(), remainCurves[0]);
				uCreateGordenCurves.insert(uCreateGordenCurves.end(), remainCurves[1]);
			}
			else
			{
				// 调整v方向的曲线顺序
				uCreateGordenCurves.clear();
				uCreateGordenCurves.insert(uCreateGordenCurves.begin(), GordenISOCurves.begin(), GordenISOCurves.end());
				vCreateGordenCurves.insert(vCreateGordenCurves.begin(), remainCurves[0]);
				vCreateGordenCurves.insert(vCreateGordenCurves.end(), remainCurves[1]);
			}

			// 对生成的曲线进行排序并检查交点
			MathTool::SortBSplineCurves(uCreateGordenCurves, uCreateGordenCurves[0]);
			MathTool::SortBSplineCurves(vCreateGordenCurves, vCreateGordenCurves[0]);
			TopoDS_Face GordenFace;
			GordenSurface::BuildMyGordonSurf(uCreateGordenCurves, vCreateGordenCurves, GordenFace);
			Handle(Geom_Surface) geomSurface = BRep_Tool::Surface(GordenFace);
			referSurface = Handle(Geom_BSplineSurface)::DownCast(geomSurface);
			/*uInternalCurve.clear();
			vInternalCurve.clear();
			uInternalCurve = uCreateGordenCurves;
			vInternalCurve = vCreateGordenCurves;*/
			// 将生成的面转换为BSplineSurface

		}

		// 返回生成的参考曲面
		return referSurface;

	}
	if (referSurfaceType == ReferSurfaceType::GORDEN_TWO_DIRECTION_GORDEN)
	{
		TopoDS_Face GordenFace;
		GordenSurface::BuildMyGordonSurf(uInternalCurve, vInternalCurve, GordenFace);

		// 将生成的面转换为BSplineSurface
		Handle(Geom_Surface) geomSurface = BRep_Tool::Surface(GordenFace);
		return Handle(Geom_BSplineSurface)::DownCast(geomSurface);
	}
}

PlanarCurve::PlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance)
	: curveType(CurveType::NOTPLANAR), curve(theCurve), line(), plane()
{
	IsPlanarCurve(curve, 10);
}

PlanarCurve::PlanarCurve()
	:curveType(CurveType::NOTPLANAR), curve(), line(), plane()
{
	curveType = CurveType::NOTPLANAR;
}

bool PlanarCurve::IsPlanarCurve(Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	// 检查曲线是否为空
	if (theCurve.IsNull())
	{
		std::cerr << "输入的曲线为空。" << std::endl;
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	bool isLinear = IsBSplineCurveLinear(theCurve);
	if (isLinear)
	{
		curveType = CurveType::LINEAR;
		line = gp_Lin(theCurve->StartPoint(), gp_Vec(theCurve->StartPoint(), theCurve->EndPoint()));
		return true;
	}

	bool isPoint = IsBSplineCurvePoint(theCurve);
	if (isPoint)
	{
		curveType = CurveType::POINT;
		return true;
	}

	// 采样曲线上的点
	std::vector<gp_Pnt> boundary_Sampling;
	int numSamples = 100; // 采样点数量，可根据需要调整
	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();
	Standard_Real step = (lastParam - firstParam) / (numSamples - 1);

	boundary_Sampling.reserve(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		Standard_Real param = firstParam + i * step;
		gp_Pnt pnt;
		theCurve->D0(param, pnt); // 获取曲线上的点
		boundary_Sampling.push_back(pnt);
	}

	// 检查是否成功采样
	if (boundary_Sampling.empty())
	{
		std::cerr << "曲线采样失败，没有采样到任何点。" << std::endl;
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	// 声明法向量和坐标轴
	gp_XYZ N, X, Y;
	// 计算质心 P	
	gp_Pnt P;
	P.SetCoord(0, 0, 0);
	int num = boundary_Sampling.size();
	for (const auto& point : boundary_Sampling)
	{
		P.SetCoord(P.X() + point.X(), P.Y() + point.Y(), P.Z() + point.Z());
	}
	P.SetCoord(P.X() / num, P.Y() / num, P.Z() / num);

	// 替换协方差矩阵构建和特征值分解部分
	int m = 3;
	Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(m, m);

	double c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
	gp_Pnt point; // 声明 'point'

	for (int i = 0; i < num; i++)
	{
		point = boundary_Sampling[i];

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

	double eigenvalue1 = E1.col(0)[0].real();
	double eigenvalue2 = E1.col(0)[1].real();
	double eigenvalue3 = E1.col(0)[2].real();

	// 根据最小特征值确定法向量和坐标轴
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

	// 创建 gp_Dir 对象用于 gp_pln
	gp_Dir dirN(N);
	plane = gp_Pln(P, dirN);
	// 评估平面性
	double maxDistance = 0.0;
	double sumDistance = 0.0;
	for (const auto& point : boundary_Sampling)
	{
		double distance = MathTool::ComputeDistancePointToPlane(point, plane);
		sumDistance += distance;
		if (distance > maxDistance)
		{
			maxDistance = distance;
		}
	}

	double averageDistance = sumDistance / boundary_Sampling.size();

	// 设置平面性阈值（根据具体需求调整）
	if (maxDistance < tolerance)
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

bool PlanarCurve::IsBSplineCurveLinear(const Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	if (theCurve.IsNull())
	{
		std::cerr << "The B-Spline curve is null." << std::endl;
		return false;
	}

	// 获取控制点数量
	Standard_Integer numPoles = theCurve->NbPoles();
	if (numPoles < 2)
	{
		// 少于两个控制点，无法确定直线
		std::cerr << "The B-Spline curve has fewer than two control points." << std::endl;
		return false;
	}

	// 获取第一个和第二个控制点
	gp_Pnt P0 = theCurve->Pole(1); // 注意：OpenCASCADE的索引从1开始
	gp_Pnt P1 = theCurve->Pole(2);

	// 计算方向向量 v = P1 - P0
	gp_Vec v(P0, P1);
	if (v.Magnitude() < tolerance)
	{
		// 前两个点重合，无法定义方向
		std::cerr << "The first two control points are coincident; cannot define a direction." << std::endl;
		return false;
	}

	// 对于每一个后续的控制点，检查是否与方向向量共线
	for (Standard_Integer i = 3; i <= numPoles; ++i)
	{
		gp_Pnt Pi = theCurve->Pole(i);
		gp_Vec u(P0, Pi);

		// 计算叉积 v × u
		gp_Vec cross = v.Crossed(u);

		// 检查叉积的模长是否在容差范围内
		if (cross.Magnitude() > tolerance)
		{
			// 不共线
			return false;
		}
	}
	// 所有控制点共线
	return true;
}

bool PlanarCurve::IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	// 获取控制点数量
	Standard_Integer numPoles = theCurve->NbPoles();
	// 获取第一个和第二个控制点
	gp_Pnt P0 = theCurve->Pole(1); // 注意：OpenCASCADE的索引从1开始
	gp_Pnt P1 = theCurve->Pole(2);

	// 计算方向向量 v = P1 - P0
	gp_Vec v(P0, P1);
	if (numPoles == 2 && v.Magnitude() < tolerance)
	{
		return true;
	}
	return false;
}

double MathTool::ComputeCurveCurveDistance(const Handle(Geom_BSplineCurve)& curve, const Handle(Geom_BSplineCurve)& boundaryCurve)
{
	GeomAPI_ExtremaCurveCurve extrema(curve, boundaryCurve);
	// 检查是否找到了极值点
	if (extrema.NbExtrema() > 0)
	{
		// 遍历所有极值点，找到最小距离
		double minDistance = RealLast();
		for (int i = 1; i <= extrema.NbExtrema(); ++i)
		{
			Standard_Real dist = extrema.Distance(i);
			if (dist < minDistance)
			{
				minDistance = dist;
			}
		}
		return minDistance;
	}

	return INT_MAX;
}

// 计算曲线的采样点平均坐标
gp_Pnt MathTool::ComputeAverageSamplePoint(const Handle(Geom_BSplineCurve)& curve, int numSamples)
{
	GeomAdaptor_Curve adaptor(curve);
	double startParam = adaptor.FirstParameter();
	double endParam = adaptor.LastParameter();
	double deltaParam = (endParam - startParam) / (numSamples - 1);
	double X = 0, Y = 0, Z = 0;
	for (int i = 0; i < numSamples; ++i) {
		double param = startParam + i * deltaParam;
		gp_Pnt sample = adaptor.Value(param);
		X += sample.X();
		Y += sample.Y();
		Z += sample.Z();
	}
	X /= numSamples;
	Y /= numSamples;
	Z /= numSamples;
	return gp_Pnt(X, Y, Z);
}

double MathTool::ComputeAngleWithAxis(const gp_Vec& vec, const gp_Vec& axis)
{
	double dotProduct = vec.Dot(axis); // 点积
	double magnitudeVec = vec.Magnitude(); // 向量的模
	double magnitudeAxis = axis.Magnitude(); // 轴的模
	double cosAngle = dotProduct / (magnitudeVec * magnitudeAxis); // 余弦值
	// 防止数值误差导致超出 [-1, 1] 范围
	cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
	return std::acos(cosAngle); // 返回夹角
}

void MathTool::CheckSelfIntersect(std::vector<Handle(Geom_BSplineCurve)> theBSplineCurvesArray)
{
	for (int i = 0; i < theBSplineCurvesArray.size(); i++)
	{
		for (int j = i + 1; j < theBSplineCurvesArray.size(); j++)
		{
			double distance = MathTool::ComputeCurveCurveDistance(theBSplineCurvesArray[i], theBSplineCurvesArray[j]);
			if (distance < 1e-3)
			{
				// 曲线自交，保留距离主平面更近的曲线
				gp_Pnt pnt1 = MathTool::ComputeAverageSamplePoint(theBSplineCurvesArray[i], 10);
				gp_Pnt pnt2 = MathTool::ComputeAverageSamplePoint(theBSplineCurvesArray[j], 10);

				// 计算曲线到原点的向量
				gp_Vec vec1 = pnt1.XYZ() - gp_Pnt(0, 0, 0).XYZ();
				gp_Vec vec2 = pnt2.XYZ() - gp_Pnt(0, 0, 0).XYZ();

				// 定义x轴、y轴、z轴
				gp_Vec xAxis(1, 0, 0);
				gp_Vec yAxis(0, 1, 0);
				gp_Vec zAxis(0, 0, 1);

				// 计算两个向量与x、y、z轴的夹角
				double angle1_x = ComputeAngleWithAxis(vec1, xAxis);
				double angle2_x = ComputeAngleWithAxis(vec2, xAxis);
				double angle1_y = ComputeAngleWithAxis(vec1, yAxis);
				double angle2_y = ComputeAngleWithAxis(vec2, yAxis);
				double angle1_z = ComputeAngleWithAxis(vec1, zAxis);
				double angle2_z = ComputeAngleWithAxis(vec2, zAxis);

				// 比较与各轴的夹角，保留夹角更小的曲线
				if (angle1_x < angle2_x || angle1_y < angle2_y || angle1_z < angle2_z)
				{
					// 保留 theBSplineCurvesArray[i]，移除 theBSplineCurvesArray[j]
					theBSplineCurvesArray.erase(theBSplineCurvesArray.begin() + j);
					j--; // 保证不跳过下一个元素
				}
				else
				{
					// 保留 theBSplineCurvesArray[j]，移除 theBSplineCurvesArray[i]
					theBSplineCurvesArray.erase(theBSplineCurvesArray.begin() + i);
					i--; // 保证不跳过下一个元素
					break; // 退出内层循环，因为 i 已经改变
				}
			}
			// 还要处理如果distance没有达到阈值，但是两个曲线的夹角很大的情况
		}
	}
}

gp_Dir MathTool::ComputeAverageTangent(const Handle(Geom_BSplineCurve)& curve, int numSamples)
{
	if (curve.IsNull())
	{
		throw std::invalid_argument("Curve is null.");
	}
	if (numSamples <= 0)
	{
		throw std::invalid_argument("Number of samples must be positive.");
	}

	Standard_Real firstParam = curve->FirstParameter();
	Standard_Real lastParam = curve->LastParameter();
	Standard_Real step = (lastParam - firstParam) / (numSamples - 1);

	gp_Vec sumTangent(0.0, 0.0, 0.0);
	int validSamples = 0;

	for (int i = 0; i < numSamples; ++i)
	{
		Standard_Real param = firstParam + i * step;
		gp_Pnt pnt;
		gp_Vec tangent;
		curve->D1(param, pnt, tangent); // D1 获取点和一阶导数
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

double MathTool::ComputeAngleBetweenCurves(Handle(Geom_BSplineCurve)& curve1,
	Handle(Geom_BSplineCurve)& curve2,
	int numSamples)
{
	PlanarCurve PlanarCurve1(curve1);
	PlanarCurve PlanarCurve2(curve2);
	if (PlanarCurve1.GetCurveType() != CurveType::NOTPLANAR && PlanarCurve2.GetCurveType() != CurveType::NOTPLANAR)
	{
		return MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurve1, PlanarCurve2);
	}

	gp_Dir avgDir1 = ComputeAverageTangent(curve1, numSamples);
	gp_Dir avgDir2 = ComputeAverageTangent(curve2, numSamples);

	double dotProduct = avgDir1.Dot(avgDir2);

	// 确保点积在 [-1, 1] 范围内，以避免数值误差
	dotProduct = std::max(-1.0, std::min(1.0, dotProduct));

	double angleRad = std::acos(dotProduct);
	double angleDeg = angleRad * 180.0 / M_PI;

	return angleDeg;
}

void MathTool::SortBSplineCurves(std::vector<Handle(Geom_BSplineCurve)>& theCurves,
	Handle(Geom_BSplineCurve) referCurve)
{
	int numSamples = 10;
	// 检查曲线数组是否为空
	if (theCurves.empty())
	{
		std::cerr << "错误：传入的曲线数组为空，无法进行排序。" << std::endl;
		return;
	}

	// 计算参考点，使用数组中的第一个曲线
	gp_Pnt referencePoint = MathTool::ComputeAverageSamplePoint(referCurve, numSamples);

	// 使用 std::sort 对曲线数组进行排序
	std::sort(theCurves.begin(), theCurves.end(),
		[&](const Handle(Geom_BSplineCurve)& a, const Handle(Geom_BSplineCurve)& b) -> bool
		{
			// 计算每条曲线的平均采样点
			gp_Pnt aAvg = MathTool::ComputeAverageSamplePoint(a, numSamples);
			gp_Pnt bAvg = MathTool::ComputeAverageSamplePoint(b, numSamples);

			// 计算平均采样点到参考点的距离
			double distA = referencePoint.Distance(aAvg);
			double distB = referencePoint.Distance(bAvg);

			// 按照距离从小到大排序
			return distA < distB;
		}
	);

	std::cout << "曲线数组已根据平均采样点到参考点的距离成功排序。" << std::endl;
}

void MathTool::ReverseIfNeeded(std::vector<Handle(Geom_BSplineCurve)>& curves)
{
	// 检查曲线数组是否为空
	if (curves.empty()) return;

	// 获取第一条曲线作为初始参考
	Handle(Geom_BSplineCurve) firstCurve = curves[0];
	gp_Pnt firstStart = firstCurve->StartPoint();
	gp_Pnt firstEnd = firstCurve->EndPoint();
	gp_Vec firstDirection(firstStart, firstEnd); // 计算参考曲线的方向向量

	// 如果第一条曲线是一个点（方向向量为零向量），则选择最后一条曲线作为参考
	if (firstDirection.Magnitude() == 0)
	{
		if (curves.size() < 2)
		{
			return;
		}
		firstCurve = curves.back();
		firstStart = firstCurve->StartPoint();
		firstEnd = firstCurve->EndPoint();
		firstDirection = gp_Vec(firstStart, firstEnd); // 重新计算参考曲线的方向向量

		// 再次检查方向向量是否为零向量
		if (firstDirection.Magnitude() == 0)
		{
			return;
		}
	}

	// 遍历每条曲线，检查方向并反转
	for (auto& curve : curves)
	{
		gp_Pnt start = curve->StartPoint();
		gp_Pnt end = curve->EndPoint();
		gp_Vec direction(start, end); // 计算当前曲线的方向向量

		// 如果当前曲线的方向与参考方向相反，则反转曲线
		if (direction.Dot(firstDirection) < 0)
		{
			curve->Reverse(); // 反转曲线方向
		}
	}
}


double MathTool::ComputeAngleBetweenLineAndPlane(const gp_Lin& line, const gp_Pln& plane)
{
	// 获取直线的方向向量
	gp_Vec lineDirection = line.Direction();  // 直线方向向量

	// 获取平面的法向量
	gp_Vec planeNormal = plane.Axis().Direction();  // 平面法向量

	// 计算直线方向向量与平面法向量的点积
	double dotProduct = lineDirection.Dot(planeNormal);

	// 计算直线方向向量和法向量的模长
	double lineMagnitude = lineDirection.Magnitude();
	double planeNormalMagnitude = planeNormal.Magnitude();

	// 计算夹角（弧度）
	double angle = std::acos(dotProduct / (lineMagnitude * planeNormalMagnitude));

	return angle;
}

double MathTool::ComputeAngleBetweenLines(const gp_Lin& line1, const gp_Lin& line2)
{
	// 获取两条直线的方向向量
	gp_Vec direction1 = line1.Direction();
	gp_Vec direction2 = line2.Direction();

	// 计算方向向量的点积
	double dotProduct = direction1.Dot(direction2);

	// 计算方向向量的模长
	double magnitude1 = direction1.Magnitude();
	double magnitude2 = direction2.Magnitude();

	// 计算夹角（弧度）
	double angle = std::acos(dotProduct / (magnitude1 * magnitude2));

	return angle;
}

double MathTool::ComputeAngleBetweenPlanes(const gp_Pln& plane1, const gp_Pln& plane2)
{
	// 获取平面的法向量
	gp_Dir normal1 = plane1.Axis().Direction();
	gp_Dir normal2 = plane2.Axis().Direction();

	// 计算法向量之间的点积
	double dotProduct = normal1.Dot(normal2);

	// 限制 dotProduct 的范围在 [-1, 1] 之间，以防止数值误差导致的 acos 计算错误
	dotProduct = std::max(-1.0, std::min(1.0, dotProduct));

	// 计算夹角的弧度值
	double angleRad = std::acos(dotProduct);

	// 将弧度转换为度数
	double angleDeg = angleRad * 180.0 / M_PI;

	// 获取锐角（0°到90°）
	if (angleDeg > 90.0)
	{
		angleDeg = 180.0 - angleDeg;
	}

	return angleDeg;
}

double MathTool::ComputeDistancePointToLine(const gp_Pnt& point, const gp_Lin& line)
{
	// 直线的方向向量
	gp_Vec lineDirection = line.Direction();

	// 点到直线上的任意一点（我们选择直线上的一个点作为参考点）
	gp_Pnt linePoint = line.Location(); // 直线上的一个点

	// 点与直线上的点的向量
	gp_Vec pointToLine = point.XYZ() - linePoint.XYZ();

	// 计算点到直线的垂直距离
	gp_Vec crossProduct = pointToLine.Crossed(lineDirection);

	// 返回距离，使用叉乘的模长除以方向向量的模长
	return crossProduct.Magnitude() / lineDirection.Magnitude();
}

double MathTool::ComputeAngleBetweenPlanarCurves(const PlanarCurve& curve1, const PlanarCurve& curve2)
{
	double angle = INT_MAX;
	double distance = INT_MIN;
	if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::LINEAR)
	{
		// 调用直线直线求交
		angle = ComputeAngleBetweenLines(curve1.GetLine(), curve2.GetLine());
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::PLANAR)
	{
		// 调用平面平面求交
		angle = ComputeAngleBetweenPlanes(curve1.GetPlane(), curve2.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::LINEAR)
	{
		// 调用平面和直线求交
		angle = ComputeAngleBetweenLineAndPlane(curve2.GetLine(), curve1.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::PLANAR)
	{
		// 调用平面和直线求交
		angle = ComputeAngleBetweenLineAndPlane(curve1.GetLine(), curve2.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::POINT && curve2.GetCurveType() == CurveType::POINT)
	{
		distance = curve1.GetPoint().Distance(curve2.GetPoint());
	}
	else if (curve1.GetCurveType() == CurveType::POINT && curve2.GetCurveType() == CurveType::PLANAR)
	{
		distance = ComputeDistancePointToPlane(curve1.GetPoint(), curve2.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::POINT && curve2.GetCurveType() == CurveType::LINEAR)
	{
		distance = ComputeDistancePointToLine(curve1.GetPoint(), curve2.GetLine());
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::POINT)
	{
		distance = ComputeDistancePointToPlane(curve2.GetPoint(), curve1.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::POINT)
	{
		distance = ComputeDistancePointToLine(curve2.GetPoint(), curve1.GetLine());
	}

	if (distance > 1) angle = 0.0;

	return angle;
}

std::vector<gp_Pnt> MathTool::GetSamplePointsOnCurve(const Handle(Geom_Curve)& curve, int numPoints)
{
	std::vector<gp_Pnt> sampledPoints;

	double firstParam = curve->FirstParameter();
	double lastParam = curve->LastParameter();
	double delta = (lastParam - firstParam) / (numPoints - 1);

	for (int i = 0; i < numPoints; ++i)
	{
		double param = firstParam + i * delta;
		gp_Pnt p;
		curve->D0(param, p);
		sampledPoints.push_back(p);
	}

	return sampledPoints;
}

void MathTool::SortPoints(std::vector<gp_Pnt>& thePoints, const gp_Pnt& referPoint)
{
	// 检查点数组是否为空
	if (thePoints.empty())
	{
		return;
	}

	// 使用 std::sort 对点数组进行排序
	std::sort(thePoints.begin(), thePoints.end(),
		[&](const gp_Pnt& a, const gp_Pnt& b) -> bool
		{
			// 计算每个点到参考点的距离
			double distA = a.Distance(referPoint);
			double distB = b.Distance(referPoint);

			// 按照距离从小到大排序
			return distA < distB;
		}
	);
}

double MathTool::ComputeDistancePointToPlane(const gp_Pnt& p, const gp_Pln& plane)
{
	// 获取平面的法向量和点
	gp_Dir normal = plane.Axis().Direction();
	gp_Pnt planeP = plane.Location();

	// 向量 (p - planeP)
	gp_Vec vec(p.X() - planeP.X(), p.Y() - planeP.Y(), p.Z() - planeP.Z());

	// 点到平面的距离 = |vec . normal| / |normal|
	// 由于 normal 是单位向量，分母为1
	double distance = std::abs(vec.Dot(normal));
	return distance;
}
