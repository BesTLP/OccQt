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
				extrema2.NearestPoints(internalPnt, replacePnt2);
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
			GordenSurface::BuildMyGordonSurf(uCreateGordenCurves, vCreateGordenCurves, GordenFace);
			Handle(Geom_Surface) geomSurface = BRep_Tool::Surface(GordenFace);
			referSurface = Handle(Geom_BSplineSurface)::DownCast(geomSurface);
			/*uInternalCurve.clear();
			vInternalCurve.clear();
			uInternalCurve = uCreateGordenCurves;
			vInternalCurve = vCreateGordenCurves;*/
			// �����ɵ���ת��ΪBSplineSurface

		}

		// �������ɵĲο�����
		return referSurface;

	}
	if (referSurfaceType == ReferSurfaceType::GORDEN_TWO_DIRECTION_GORDEN)
	{
		TopoDS_Face GordenFace;
		GordenSurface::BuildMyGordonSurf(uInternalCurve, vInternalCurve, GordenFace);

		// �����ɵ���ת��ΪBSplineSurface
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
	// ��������Ƿ�Ϊ��
	if (theCurve.IsNull())
	{
		std::cerr << "���������Ϊ�ա�" << std::endl;
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

	// ���������ϵĵ�
	std::vector<gp_Pnt> boundary_Sampling;
	int numSamples = 100; // �������������ɸ�����Ҫ����
	Standard_Real firstParam = theCurve->FirstParameter();
	Standard_Real lastParam = theCurve->LastParameter();
	Standard_Real step = (lastParam - firstParam) / (numSamples - 1);

	boundary_Sampling.reserve(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		Standard_Real param = firstParam + i * step;
		gp_Pnt pnt;
		theCurve->D0(param, pnt); // ��ȡ�����ϵĵ�
		boundary_Sampling.push_back(pnt);
	}

	// ����Ƿ�ɹ�����
	if (boundary_Sampling.empty())
	{
		std::cerr << "���߲���ʧ�ܣ�û�в������κε㡣" << std::endl;
		curveType = CurveType::NOTPLANAR;
		return false;
	}

	// ������������������
	gp_XYZ N, X, Y;
	// �������� P	
	gp_Pnt P;
	P.SetCoord(0, 0, 0);
	int num = boundary_Sampling.size();
	for (const auto& point : boundary_Sampling)
	{
		P.SetCoord(P.X() + point.X(), P.Y() + point.Y(), P.Z() + point.Z());
	}
	P.SetCoord(P.X() / num, P.Y() / num, P.Z() / num);

	// �滻Э������󹹽�������ֵ�ֽⲿ��
	int m = 3;
	Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(m, m);

	double c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
	gp_Pnt point; // ���� 'point'

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

	// ����ƽ������ֵ�����ݾ������������
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

	// ��ȡ���Ƶ�����
	Standard_Integer numPoles = theCurve->NbPoles();
	if (numPoles < 2)
	{
		// �����������Ƶ㣬�޷�ȷ��ֱ��
		std::cerr << "The B-Spline curve has fewer than two control points." << std::endl;
		return false;
	}

	// ��ȡ��һ���͵ڶ������Ƶ�
	gp_Pnt P0 = theCurve->Pole(1); // ע�⣺OpenCASCADE��������1��ʼ
	gp_Pnt P1 = theCurve->Pole(2);

	// ���㷽������ v = P1 - P0
	gp_Vec v(P0, P1);
	if (v.Magnitude() < tolerance)
	{
		// ǰ�������غϣ��޷����巽��
		std::cerr << "The first two control points are coincident; cannot define a direction." << std::endl;
		return false;
	}

	// ����ÿһ�������Ŀ��Ƶ㣬����Ƿ��뷽����������
	for (Standard_Integer i = 3; i <= numPoles; ++i)
	{
		gp_Pnt Pi = theCurve->Pole(i);
		gp_Vec u(P0, Pi);

		// ������ v �� u
		gp_Vec cross = v.Crossed(u);

		// �������ģ���Ƿ����ݲΧ��
		if (cross.Magnitude() > tolerance)
		{
			// ������
			return false;
		}
	}
	// ���п��Ƶ㹲��
	return true;
}

bool PlanarCurve::IsBSplineCurvePoint(const Handle(Geom_BSplineCurve)& theCurve, double tolerance)
{
	// ��ȡ���Ƶ�����
	Standard_Integer numPoles = theCurve->NbPoles();
	// ��ȡ��һ���͵ڶ������Ƶ�
	gp_Pnt P0 = theCurve->Pole(1); // ע�⣺OpenCASCADE��������1��ʼ
	gp_Pnt P1 = theCurve->Pole(2);

	// ���㷽������ v = P1 - P0
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
	// ����Ƿ��ҵ��˼�ֵ��
	if (extrema.NbExtrema() > 0)
	{
		// �������м�ֵ�㣬�ҵ���С����
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

// �������ߵĲ�����ƽ������
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
	double dotProduct = vec.Dot(axis); // ���
	double magnitudeVec = vec.Magnitude(); // ������ģ
	double magnitudeAxis = axis.Magnitude(); // ���ģ
	double cosAngle = dotProduct / (magnitudeVec * magnitudeAxis); // ����ֵ
	// ��ֹ��ֵ���³��� [-1, 1] ��Χ
	cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
	return std::acos(cosAngle); // ���ؼн�
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
				// �����Խ�������������ƽ�����������
				gp_Pnt pnt1 = MathTool::ComputeAverageSamplePoint(theBSplineCurvesArray[i], 10);
				gp_Pnt pnt2 = MathTool::ComputeAverageSamplePoint(theBSplineCurvesArray[j], 10);

				// �������ߵ�ԭ�������
				gp_Vec vec1 = pnt1.XYZ() - gp_Pnt(0, 0, 0).XYZ();
				gp_Vec vec2 = pnt2.XYZ() - gp_Pnt(0, 0, 0).XYZ();

				// ����x�ᡢy�ᡢz��
				gp_Vec xAxis(1, 0, 0);
				gp_Vec yAxis(0, 1, 0);
				gp_Vec zAxis(0, 0, 1);

				// ��������������x��y��z��ļн�
				double angle1_x = ComputeAngleWithAxis(vec1, xAxis);
				double angle2_x = ComputeAngleWithAxis(vec2, xAxis);
				double angle1_y = ComputeAngleWithAxis(vec1, yAxis);
				double angle2_y = ComputeAngleWithAxis(vec2, yAxis);
				double angle1_z = ComputeAngleWithAxis(vec1, zAxis);
				double angle2_z = ComputeAngleWithAxis(vec2, zAxis);

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
			// ��Ҫ�������distanceû�дﵽ��ֵ�������������ߵļнǺܴ�����
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
		curve->D1(param, pnt, tangent); // D1 ��ȡ���һ�׵���
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

	// ȷ������� [-1, 1] ��Χ�ڣ��Ա�����ֵ���
	dotProduct = std::max(-1.0, std::min(1.0, dotProduct));

	double angleRad = std::acos(dotProduct);
	double angleDeg = angleRad * 180.0 / M_PI;

	return angleDeg;
}

void MathTool::SortBSplineCurves(std::vector<Handle(Geom_BSplineCurve)>& theCurves,
	Handle(Geom_BSplineCurve) referCurve)
{
	int numSamples = 10;
	// ������������Ƿ�Ϊ��
	if (theCurves.empty())
	{
		std::cerr << "���󣺴������������Ϊ�գ��޷���������" << std::endl;
		return;
	}

	// ����ο��㣬ʹ�������еĵ�һ������
	gp_Pnt referencePoint = MathTool::ComputeAverageSamplePoint(referCurve, numSamples);

	// ʹ�� std::sort �����������������
	std::sort(theCurves.begin(), theCurves.end(),
		[&](const Handle(Geom_BSplineCurve)& a, const Handle(Geom_BSplineCurve)& b) -> bool
		{
			// ����ÿ�����ߵ�ƽ��������
			gp_Pnt aAvg = MathTool::ComputeAverageSamplePoint(a, numSamples);
			gp_Pnt bAvg = MathTool::ComputeAverageSamplePoint(b, numSamples);

			// ����ƽ�������㵽�ο���ľ���
			double distA = referencePoint.Distance(aAvg);
			double distB = referencePoint.Distance(bAvg);

			// ���վ����С��������
			return distA < distB;
		}
	);

	std::cout << "���������Ѹ���ƽ�������㵽�ο���ľ���ɹ�����" << std::endl;
}

void MathTool::ReverseIfNeeded(std::vector<Handle(Geom_BSplineCurve)>& curves)
{
	// ������������Ƿ�Ϊ��
	if (curves.empty()) return;

	// ��ȡ��һ��������Ϊ��ʼ�ο�
	Handle(Geom_BSplineCurve) firstCurve = curves[0];
	gp_Pnt firstStart = firstCurve->StartPoint();
	gp_Pnt firstEnd = firstCurve->EndPoint();
	gp_Vec firstDirection(firstStart, firstEnd); // ����ο����ߵķ�������

	// �����һ��������һ���㣨��������Ϊ������������ѡ�����һ��������Ϊ�ο�
	if (firstDirection.Magnitude() == 0)
	{
		if (curves.size() < 2)
		{
			return;
		}
		firstCurve = curves.back();
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
	for (auto& curve : curves)
	{
		gp_Pnt start = curve->StartPoint();
		gp_Pnt end = curve->EndPoint();
		gp_Vec direction(start, end); // ���㵱ǰ���ߵķ�������

		// �����ǰ���ߵķ�����ο������෴����ת����
		if (direction.Dot(firstDirection) < 0)
		{
			curve->Reverse(); // ��ת���߷���
		}
	}
}


double MathTool::ComputeAngleBetweenLineAndPlane(const gp_Lin& line, const gp_Pln& plane)
{
	// ��ȡֱ�ߵķ�������
	gp_Vec lineDirection = line.Direction();  // ֱ�߷�������

	// ��ȡƽ��ķ�����
	gp_Vec planeNormal = plane.Axis().Direction();  // ƽ�淨����

	// ����ֱ�߷���������ƽ�淨�����ĵ��
	double dotProduct = lineDirection.Dot(planeNormal);

	// ����ֱ�߷��������ͷ�������ģ��
	double lineMagnitude = lineDirection.Magnitude();
	double planeNormalMagnitude = planeNormal.Magnitude();

	// ����нǣ����ȣ�
	double angle = std::acos(dotProduct / (lineMagnitude * planeNormalMagnitude));

	return angle;
}

double MathTool::ComputeAngleBetweenLines(const gp_Lin& line1, const gp_Lin& line2)
{
	// ��ȡ����ֱ�ߵķ�������
	gp_Vec direction1 = line1.Direction();
	gp_Vec direction2 = line2.Direction();

	// ���㷽�������ĵ��
	double dotProduct = direction1.Dot(direction2);

	// ���㷽��������ģ��
	double magnitude1 = direction1.Magnitude();
	double magnitude2 = direction2.Magnitude();

	// ����нǣ����ȣ�
	double angle = std::acos(dotProduct / (magnitude1 * magnitude2));

	return angle;
}

double MathTool::ComputeAngleBetweenPlanes(const gp_Pln& plane1, const gp_Pln& plane2)
{
	// ��ȡƽ��ķ�����
	gp_Dir normal1 = plane1.Axis().Direction();
	gp_Dir normal2 = plane2.Axis().Direction();

	// ���㷨����֮��ĵ��
	double dotProduct = normal1.Dot(normal2);

	// ���� dotProduct �ķ�Χ�� [-1, 1] ֮�䣬�Է�ֹ��ֵ���µ� acos �������
	dotProduct = std::max(-1.0, std::min(1.0, dotProduct));

	// ����нǵĻ���ֵ
	double angleRad = std::acos(dotProduct);

	// ������ת��Ϊ����
	double angleDeg = angleRad * 180.0 / M_PI;

	// ��ȡ��ǣ�0�㵽90�㣩
	if (angleDeg > 90.0)
	{
		angleDeg = 180.0 - angleDeg;
	}

	return angleDeg;
}

double MathTool::ComputeDistancePointToLine(const gp_Pnt& point, const gp_Lin& line)
{
	// ֱ�ߵķ�������
	gp_Vec lineDirection = line.Direction();

	// �㵽ֱ���ϵ�����һ�㣨����ѡ��ֱ���ϵ�һ������Ϊ�ο��㣩
	gp_Pnt linePoint = line.Location(); // ֱ���ϵ�һ����

	// ����ֱ���ϵĵ������
	gp_Vec pointToLine = point.XYZ() - linePoint.XYZ();

	// ����㵽ֱ�ߵĴ�ֱ����
	gp_Vec crossProduct = pointToLine.Crossed(lineDirection);

	// ���ؾ��룬ʹ�ò�˵�ģ�����Է���������ģ��
	return crossProduct.Magnitude() / lineDirection.Magnitude();
}

double MathTool::ComputeAngleBetweenPlanarCurves(const PlanarCurve& curve1, const PlanarCurve& curve2)
{
	double angle = INT_MAX;
	double distance = INT_MIN;
	if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::LINEAR)
	{
		// ����ֱ��ֱ����
		angle = ComputeAngleBetweenLines(curve1.GetLine(), curve2.GetLine());
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::PLANAR)
	{
		// ����ƽ��ƽ����
		angle = ComputeAngleBetweenPlanes(curve1.GetPlane(), curve2.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::PLANAR && curve2.GetCurveType() == CurveType::LINEAR)
	{
		// ����ƽ���ֱ����
		angle = ComputeAngleBetweenLineAndPlane(curve2.GetLine(), curve1.GetPlane());
	}
	else if (curve1.GetCurveType() == CurveType::LINEAR && curve2.GetCurveType() == CurveType::PLANAR)
	{
		// ����ƽ���ֱ����
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
	// ���������Ƿ�Ϊ��
	if (thePoints.empty())
	{
		return;
	}

	// ʹ�� std::sort �Ե������������
	std::sort(thePoints.begin(), thePoints.end(),
		[&](const gp_Pnt& a, const gp_Pnt& b) -> bool
		{
			// ����ÿ���㵽�ο���ľ���
			double distA = a.Distance(referPoint);
			double distB = b.Distance(referPoint);

			// ���վ����С��������
			return distA < distB;
		}
	);
}

double MathTool::ComputeDistancePointToPlane(const gp_Pnt& p, const gp_Pln& plane)
{
	// ��ȡƽ��ķ������͵�
	gp_Dir normal = plane.Axis().Direction();
	gp_Pnt planeP = plane.Location();

	// ���� (p - planeP)
	gp_Vec vec(p.X() - planeP.X(), p.Y() - planeP.Y(), p.Z() - planeP.Z());

	// �㵽ƽ��ľ��� = |vec . normal| / |normal|
	// ���� normal �ǵ�λ��������ĸΪ1
	double distance = std::abs(vec.Dot(normal));
	return distance;
}
