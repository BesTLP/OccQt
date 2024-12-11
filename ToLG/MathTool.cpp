#include "MathTool.h"

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
