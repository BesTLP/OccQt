#include "MathTool.h"

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
