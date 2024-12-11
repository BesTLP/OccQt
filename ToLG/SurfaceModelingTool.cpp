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
				extrema1.NearestPoints(internalPnt, replacePnt2);
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
			InternalPlanarCurve.SetCurve(aBsplineCurve);
			/*
			aBsplineCurve->SetPole(1, replacePnt1);
			aBsplineCurve->SetPole(InternalPlanarCurve.GetCurve()->NbPoles(), replacePnt2);
			distance1 = aBsplineCurve->StartPoint().Distance(replacePnt1);
			distance2 = aBsplineCurve->EndPoint().Distance(replacePnt2);
			std::vector<gp_Pnt> ApproximatePoints = MathTool::GetSamplePointsOnCurve(InternalPlanarCurve.GetCurve(), 100);
			ApproximatePoints[0] = replacePnt1;
			ApproximatePoints[ApproximatePoints.size() - 1] = replacePnt2;
			MathTool::SortPoints(ApproximatePoints, ApproximatePoints[0]);
			std::vector<double> params = ComputeUniformParam(ApproximatePoints.size(), 0., 1.);
			std::vector<double> tempKnots = KnotGernerationByParams(params, 3, InternalPlanarCurve.GetCurve()->Degree());
			std::vector<double> insertKnots;
			Handle(Geom_BSplineCurve) aBSplineCurve = IterateApproximate(insertKnots, ApproximatePoints, params, tempKnots, InternalPlanarCurve.GetCurve()->Degree(), 50, 0.01);
			InternalPlanarCurve.SetCurve(aBSplineCurve);
			*/

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
			//GordenSurface::BuildMyGordonSurf(uCreateGordenCurves, vCreateGordenCurves, GordenFace);
			//Handle(Geom_Surface) geomSurface = BRep_Tool::Surface(GordenFace);
			//referSurface = Handle(Geom_BSplineSurface)::DownCast(geomSurface);
			uInternalCurve.clear();
			vInternalCurve.clear();
			uInternalCurve = uCreateGordenCurves;
			vInternalCurve = vCreateGordenCurves;
			// 将生成的面转换为BSplineSurface

		}

		// 返回生成的参考曲面
		return referSurface;

	}
}