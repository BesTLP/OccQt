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

	// ����1����ʼ��PlanarCurveArray���������б߽�����
	std::vector<PlanarCurve> PlanarCurveArray;
	for (int i = 0; i < aBoundarycurveArray.size(); i++)
	{
		PlanarCurveArray.emplace_back(PlanarCurve(aBoundarycurveArray[i]));
	}

	// ����2��������б߽������Ƿ���ƽ������
	bool canUseInternalLines = true;
	for (const PlanarCurve& curve : PlanarCurveArray)
	{
		if (curve.GetCurveType() == CurveType::NOTPLANAR)
		{
			canUseInternalLines = false; // ����з�ƽ�����ߣ�����Ϊ false
			break;
		}
	}

	if (!canUseInternalLines)
	{
		return false; // ����κ�һ���߽����߲���ƽ�����ߣ���ֱ�ӷ���false
	}

	// ���uInternalCurve��vInternalCurve��׼���洢���
	uInternalCurve.clear();
	vInternalCurve.clear();

	// �����ڲ�BSpline����
	for (auto& internalCurve : anInternalBSplineCurves)
	{
		PlanarCurve InternalPlanarCurve(internalCurve);

		// ����ڲ����߲���ƽ�����ߣ�������
		if (InternalPlanarCurve.GetCurveType() == CurveType::NOTPLANAR)
		{
			continue;
		}

		// �����ڲ����ߺ������߽�����֮��ľ���
		double distance1 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve1);
		double distance2 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve2);
		double distance3 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve3);
		double distance4 = MathTool::ComputeCurveCurveDistance(InternalPlanarCurve.GetCurve(), bslpineCurve4);

		double SplitPointParameters[2] = { 0 };

		// ����5����������Ƿ񿿽��߽����ߣ��������������ָ��
		if ((distance1 < 10 && distance3 < 10) || (distance2 < 10 && distance4 < 10))
		{
			GeomAPI_ExtremaCurveCurve extrema1(InternalPlanarCurve.GetCurve(), distance1 < 10 ? bslpineCurve1 : bslpineCurve2);
			GeomAPI_ExtremaCurveCurve extrema2(InternalPlanarCurve.GetCurve(), distance3 < 10 ? bslpineCurve3 : bslpineCurve4);
			gp_Pnt internalPnt;
			gp_Pnt replacePnt1, replacePnt2;
			// ��ȡ�ָ�����
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

			// ȷ���ָ�������ȷ����
			if (SplitPointParameters[0] > SplitPointParameters[1])
			{
				std::swap(SplitPointParameters[1], SplitPointParameters[0]);
				std::swap(replacePnt1, replacePnt2);
			}

			// ���ڲ����߽��вü�
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

			// ����ǶȲ���������
			double angle1, angle3, angle2, angle4;
			angle1 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[0], InternalPlanarCurve);
			angle2 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[1], InternalPlanarCurve);
			angle3 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[2], InternalPlanarCurve);
			angle4 = MathTool::ComputeAngleBetweenPlanarCurves(PlanarCurveArray[3], InternalPlanarCurve);

			// ���ݽǶ�ֵ�жϸ���������u������v����
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

	// ���߽�������ӵ�uInternalCurve��vInternalCurve��ͷ����β��
	uInternalCurve.insert(uInternalCurve.begin(), bslpineCurve1);
	uInternalCurve.insert(uInternalCurve.end(), bslpineCurve3);
	vInternalCurve.insert(vInternalCurve.begin(), bslpineCurve2);
	vInternalCurve.insert(vInternalCurve.end(), bslpineCurve4);

	// �������߲�����Խ�
	MathTool::SortBSplineCurves(uInternalCurve, uInternalCurve[0]);
	MathTool::SortBSplineCurves(vInternalCurve, vInternalCurve[0]);

	MathTool::CheckSelfIntersect(uInternalCurve);
	MathTool::CheckSelfIntersect(vInternalCurve);

	//SurfaceModelingTool tool;
	//tool.CompatibleWithInterPoints(uInternalCurve, vInternalCurve);
	//tool.CompatibleWithInterPoints(vInternalCurve, uInternalCurve);
	// ���uInternalCurve��vInternalCurve����������һ������4���򷵻�true
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
		// ��ȡ����ı߽�����
		Handle(Geom_BSplineCurve) bslpineCurve1 = aBoundarycurveArray[0];
		Handle(Geom_BSplineCurve) bslpineCurve2 = aBoundarycurveArray[1];
		Handle(Geom_BSplineCurve) bslpineCurve3 = aBoundarycurveArray[2];
		Handle(Geom_BSplineCurve) bslpineCurve4 = aBoundarycurveArray[3];

		// �洢���ɵ�Gorden�Ȳ������ߺ�ʣ������
		std::vector<Handle(Geom_BSplineCurve)> GordenISOCurves;
		std::vector<Handle(Geom_BSplineCurve)> remainCurves;

		// ʹ�����㷨�ı�־
		bool useNewAlgorithm = true;

		// �ж��ڲ����ߵ�������ѡ����Gorden����ķ�ʽ
		if (uInternalCurve.size() > vInternalCurve.size() && uInternalCurve.size() >= 4)
		{
			// ѡ��u������ڲ��ߺͱ߽�������Gorden����
			GordenISOCurves.insert(GordenISOCurves.end(), uInternalCurve.begin(), uInternalCurve.end());
			remainCurves.push_back(bslpineCurve2);
			remainCurves.push_back(bslpineCurve4);
		}
		else if (vInternalCurve.size() > uInternalCurve.size() && vInternalCurve.size() >= 4)
		{
			// ѡ��v������ڲ��ߺͱ߽�������Gorden����
			GordenISOCurves.insert(GordenISOCurves.end(), vInternalCurve.begin(), vInternalCurve.end());
			remainCurves.push_back(bslpineCurve1);
			remainCurves.push_back(bslpineCurve3);
		}
		else if (uInternalCurve.size() == vInternalCurve.size() && uInternalCurve.size() >= 4)
		{
			// ���u�����v������ڲ�����������ȣ����ݽǶ�֮����ѡ��
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
			// ������������㣬���˵������㷨
			useNewAlgorithm = false;
			return nullptr;
		}

		// �洢���ɵĵȲ������ߺͷ���
		std::vector<Handle(Geom_BSplineCurve)> uCreateGordenCurves, vCreateGordenCurves;
		std::vector<gp_Vec> normalsOfUISOLines, normalsOfVISOLines;

		// �������ɵĲο�����
		Handle(Geom_BSplineSurface) referSurface;

		if (useNewAlgorithm)
		{
			// ��������������֮��ĽǶ�
			double AngleUwithG = MathTool::ComputeAngleBetweenCurves(bslpineCurve1, GordenISOCurves[0]);
			double AngleVwithG = MathTool::ComputeAngleBetweenCurves(bslpineCurve2, GordenISOCurves[0]);

			// ���ݽǶ�ѡ������
			if (AngleUwithG > AngleVwithG)
			{
				// ���u����ĽǶȸ��󣬵���u��v���������˳��
				vCreateGordenCurves.clear();
				vCreateGordenCurves.insert(vCreateGordenCurves.begin(), GordenISOCurves.begin(), GordenISOCurves.end());
				uCreateGordenCurves.insert(uCreateGordenCurves.begin(), remainCurves[0]);
				uCreateGordenCurves.insert(uCreateGordenCurves.end(), remainCurves[1]);
			}
			else
			{
				// ����v���������˳��
				uCreateGordenCurves.clear();
				uCreateGordenCurves.insert(uCreateGordenCurves.begin(), GordenISOCurves.begin(), GordenISOCurves.end());
				vCreateGordenCurves.insert(vCreateGordenCurves.begin(), remainCurves[0]);
				vCreateGordenCurves.insert(vCreateGordenCurves.end(), remainCurves[1]);
			}

			// �����ɵ����߽������򲢼�齻��
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
			// �����ɵ���ת��ΪBSplineSurface

		}

		// �������ɵĲο�����
		return referSurface;

	}
}