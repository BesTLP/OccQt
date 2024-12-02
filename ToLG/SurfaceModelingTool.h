#pragma once
#include "Geom_BSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Edge.hxx"
#include "vector"
class SurfaceModelingTool
{
public:
	//约定输入四条边满足如下的方向要求
	// 四条边连接和方向满足要求 G0连续的Coons曲面构造
	// 
	//                -------------curve1--->-------
	//                |                             |
	//              curve4                        curve 2
	//                v                             v
	//                |                             |
	//                |-----------curve3--->--------
	static void Coons_G0(Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4, Handle(Geom_BSplineSurface)& mySurface_coons);

	//              约定输入四条边和四条边的跨界导矢满足如下的方向要求
	//                           G1连续的Coons曲面构造
	// 
	//                ----------------------c1--->------------
	//                |                     |                |
	//               c4-->c4_d            c1_d              c2 -->-c2_d-->
	//                v                     v                v
	//                |                     |                |
	//                |--------------------c3--->-------------
	//                                      |
	//                                      v
	//                                    c3_d  
	//                                      | 
	static void Coons_G1(Handle(Geom_BSplineCurve)& c1, Handle(Geom_BSplineCurve)& c2, Handle(Geom_BSplineCurve)& c3, Handle(Geom_BSplineCurve)& c4, Handle(Geom_BSplineCurve)& c1_derivative, Handle(Geom_BSplineCurve)& curve2_derivative, Handle(Geom_BSplineCurve)& curve3_derivative, Handle(Geom_BSplineCurve)& curve4_derivative, Handle(Geom_BSplineSurface)& mySurface_coons);

	//make the curve same degree and knots
	static Standard_Integer SetSameDistribution(Handle(Geom_BSplineCurve)& C1, Handle(Geom_BSplineCurve)& C2);

	//make the curve arranged for the G0 construction
	// Input: curveArray, the curve array
	// Output: bslpineCurve1-bslpineCurve4 the arranged four boundary curves
	// IsModify: If the curve are not connected end to end, whether change the end point of the curve
	// Tol: the tolerance to check whether the four curves are connected or not

	static int Arrange_Coons_G0(std::vector<Handle(Geom_BSplineCurve)>& curveArray, Handle(Geom_BSplineCurve)& bslpineCurve1, Handle(Geom_BSplineCurve)& bslpineCurve2, Handle(Geom_BSplineCurve)& bslpineCurve3, Handle(Geom_BSplineCurve)& bslpineCurve4, double Tol, int IsModify);

	static void ClassifyAndSortISOcurves(const std::vector<Handle(Geom_BSplineCurve)>& anISOcurvesArray, std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray, std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray);

	static void CreateLoftingSurface(const std::vector<Handle(Geom_BSplineCurve)>& curvesArray,
		const std::vector<gp_Vec>& normals,
		std::vector<TopoDS_Shape>& loftingSurfaces);

	static void LoftSurfaceIntersectWithCurve(
		const std::vector<TopoDS_Shape>& LoftingSur,
		const std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_Initial,
		const std::vector<Handle(Geom_BSplineCurve)>& anInternalBSplineCurves,
		std::vector<Handle(Geom_BSplineCurve)>& ISOcurvesArray_New,
		Standard_Integer isoCount,
		std::vector<std::vector<gp_Pnt>>& InterpolatePoints,
		Handle(Geom_BSplineSurface) CoonsSurface);

	static void CreateFinalISOCurves(
		const std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_New,
		const std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_New,
		std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Final,
		std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Final,
		std::vector<std::vector<gp_Pnt>>& uInterpolatePoints,
		std::vector<std::vector<gp_Pnt>>& vInterpolatePoints,
		Standard_Integer isoCount,
		std::vector<Handle(Geom_BSplineSurface)>& surfaceArr);



	static void UpdateFinalCurves(const std::vector<Handle(Geom_BSplineCurve)>& aBoundarycurveArray,
		std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Final,
		std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Final);

	static void LoadBSplineCurves(const std::string& filePath, std::vector<Handle(Geom_BSplineCurve)>& curveArray);
	static void LoadBSplineSurfaces(const std::string& filePath, std::vector<Handle(Geom_BSplineSurface)>& surfaceArray);
	static void GetISOCurveWithNormal(
		const Handle(Geom_BSplineSurface)& surfacecoons,
		std::vector<Handle(Geom_BSplineCurve)>& uISOcurvesArray_Initial,
		std::vector<Handle(Geom_BSplineCurve)>& vISOcurvesArray_Initial,
		std::vector<gp_Vec>& normalsOfUISOLines,
		std::vector<gp_Vec>& normalsOfVISOLines,
		int numIsoCurves = 10);

	void setKnotsOutputPath(std::string knotsOuputPath)
	{
		this->knotsOutputPath = knotsOuputPath;
	}
	std::string getKnotsOuputPath()
	{
		return this->knotsOutputPath;
	}


	// 将节点数据输出到txt文件
	Standard_Boolean KnotsToTxt(const std::vector<double>& knots) const 
	{
		// 打开文件
		std::ofstream outFile(knotsOutputPath, std::ios::app);

		if (!outFile.is_open())
		{
			return Standard_False;
		}

		// 设置输出到小数点后三位
		outFile << std::fixed << std::setprecision(3); 

		// 写入数据到文件
		outFile << "[";
		for (size_t i = 0; i < knots.size(); ++i) 
		{
			outFile << knots[i];
			if (i < knots.size() - 1)
			{
				outFile << ",";
			}
		}
		outFile << "]\n";

		outFile.close();
		return Standard_True;
	}
	Standard_Boolean ContextToTxt(const std::string context) const
	{
		std::ofstream outFile(knotsOutputPath, std::ios::app);

		if (!outFile.is_open())
		{
			return Standard_False;
		}

		outFile << context << "\n";

		outFile.close();

		return Standard_True;
	}

	static void ApproximateBoundaryCurves(std::vector<Handle(Geom_BSplineCurve)>& curves, int samplingNum = 50);

private:

	std::string knotsOutputPath;
	
};

