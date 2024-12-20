#include "Interpolate.h"

#include <map>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <gp_Pln.hxx>
#include <gp_Cylinder.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <STEPControl_Writer.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <BRep_Tool.hxx>
#include <GeomConvert.hxx>
#include <Math_matrix.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <TopoDS_Face.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Eigen/Dense>


Handle(Geom_BSplineCurve) InterPolateTool::Interpolate(const std::vector<gp_Pnt>& DataPnts, const std::vector<double>& PntParams,
	const std::vector<double>& knotSeqs, int degree) {
	int ctrl_point_num = knotSeqs.size() - degree - 1;
	if (ctrl_point_num != DataPnts.size()) {
		std::cerr << "data points number must equal to control point number in interpolation" << std::endl;
	}
	Eigen::MatrixXd matN(ctrl_point_num, ctrl_point_num);
	for (size_t i = 0; i < ctrl_point_num; i++) {
		for (size_t j = 0; j < ctrl_point_num; j++) {
			double value = OneBasicFun(PntParams[i], j, degree, knotSeqs);
			matN(i, j) = value;
		}
	}
	Eigen::VectorXd VRx(ctrl_point_num);
	Eigen::VectorXd VRy(ctrl_point_num);
	Eigen::VectorXd VRz(ctrl_point_num);
	for (size_t i = 0; i < ctrl_point_num; i++) {
		double x, y, z;
		DataPnts[i].Coord(x, y, z);
		VRx(i) = x;
		VRy(i) = y;
		VRz(i) = z;
	}
	Eigen::VectorXd Sx = matN.lu().solve(VRx);
	Eigen::VectorXd Sy = matN.lu().solve(VRy);
	Eigen::VectorXd Sz = matN.lu().solve(VRz);
	TColgp_Array1OfPnt ctrlPnts(1, ctrl_point_num);
	for (size_t i = 1; i <= ctrl_point_num; i++) {
		gp_Pnt pntTemp(Sx(i - 1), Sy(i - 1), Sz(i - 1));
		ctrlPnts[i] = pntTemp;
	}
	std::vector<double> knots;
	std::vector<int> mutis;
	sequenceToKnots(knotSeqs, knots, mutis);
	TColStd_Array1OfReal Knots(1, knots.size());
	TColStd_Array1OfInteger Mutis(1, knots.size());
	for (int i = 0; i < knots.size(); i++) {
		Knots[i + 1] = knots[i];
		Mutis[i + 1] = mutis[i];
	}
	Handle(Geom_BSplineCurve) bspline = new Geom_BSplineCurve(ctrlPnts, Knots, Mutis, degree);
	return bspline;
}

Handle(Geom_BSplineSurface) InterPolateTool::Interpolate(const std::vector<gp_Pnt>& Pnts, const std::vector<gp_Pnt2d>& PntParams,
	std::vector<double>& UKnots, std::vector<double>& VKnots, std::vector<int>& UMuti, std::vector<int>& VMuti, int Udegree, int Vdegree) {
	std::vector<double> USeq = knotsToSequence(UKnots, UMuti);
	std::vector<double> VSeq = knotsToSequence(VKnots, VMuti);
	int u_direct_poles_num = USeq.size() - Udegree - 1;
	int v_direct_poles_num = VSeq.size() - Vdegree - 1;
	Eigen::MatrixXd matNi(Pnts.size(), u_direct_poles_num);
	Eigen::MatrixXd matNj(Pnts.size(), v_direct_poles_num);
	double valueTemp = 0.0;
	for (size_t i = 0; i < Pnts.size(); i++)
	{
		for (size_t j = 0; j < u_direct_poles_num; j++) 
		{
			valueTemp = OneBasicFun(PntParams[i].X(), j, Udegree, USeq);
			matNi(i, j) = valueTemp;
		}
		for (size_t j = 0; j < v_direct_poles_num; j++)
		{
			valueTemp = OneBasicFun(PntParams[i].Y(), j, Vdegree, VSeq);
			matNj(i, j) = valueTemp;
		}
	}
	//std::cout << matNi << std::endl;
	//std::cout << matNj << std::endl;
	auto N = composeTensowMat(matNi, matNj);
	//std::cout << N.cols()<<" x " <<N.rows()<< std::endl;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(N);
	double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	std::cout << "the matrix cond num : " << cond;
	Eigen::VectorXd VBx(Pnts.size());
	Eigen::VectorXd VBy(Pnts.size());
	Eigen::VectorXd VBz(Pnts.size());
	for (size_t i = 0; i < Pnts.size(); i++) {
		double x, y, z;
		Pnts[i].Coord(x, y, z);
		VBx(i) = x;
		VBy(i) = y;
		VBz(i) = z;
	}
#ifdef MATLIBSOLVE
	Eigen::VectorXd Sx;
	Eigen::VectorXd Sy;
	Eigen::VectorXd Sz;
	solveLeastSquaresInMATLAB(N, VBx, VBy, VBz, Sx, Sy, Sz);
#else
	//Eigen::VectorXd Sx = N.householderQr().solve(VBx);
	//Eigen::VectorXd Sy = N.householderQr().solve(VBy);
	//Eigen::VectorXd Sz = N.householderQr().solve(VBz);
	Eigen::VectorXd Sx = N.fullPivHouseholderQr().solve(VBx);
	Eigen::VectorXd Sy = N.fullPivHouseholderQr().solve(VBy);
	Eigen::VectorXd Sz = N.fullPivHouseholderQr().solve(VBz);
	//Eigen::VectorXd Sx = N.lu().solve(VBx);
	//Eigen::VectorXd Sy = N.lu().solve(VBy);
	//Eigen::VectorXd Sz = N.lu().solve(VBz);
#endif // MATLIBSOLVE

	std::vector<gp_Pnt> colPoles(Sx.size());
	for (size_t i = 0; i < Sx.size(); i++) {
		gp_Pnt pntTemp(Sx(i), Sy(i), Sz(i));
		colPoles[i] = pntTemp;
	}
	TColgp_Array2OfPnt Poles = vectorToOCCMatrix(colPoles, v_direct_poles_num);
	TColStd_Array1OfReal UKnots_occ(1, UKnots.size());
	TColStd_Array1OfInteger UMuti_occ(1, UKnots.size());
	for (int i = 0; i < UKnots.size(); i++) {
		UKnots_occ[i + 1] = UKnots[i];
		UMuti_occ[i + 1] = UMuti[i];
	}
	TColStd_Array1OfReal VKnots_occ(1, VKnots.size());
	TColStd_Array1OfInteger VMuti_occ(1, VKnots.size());
	for (int i = 0; i < VKnots.size(); i++) {
		VKnots_occ[i + 1] = VKnots[i];
		VMuti_occ[i + 1] = VMuti[i];
	}
	Handle(Geom_BSplineSurface) bspline = new Geom_BSplineSurface(Poles, UKnots_occ, VKnots_occ, UMuti_occ, VMuti_occ, Udegree, Vdegree);
	return bspline;
}

Handle(Geom_BSplineSurface) InterPolateTool::Loft(const std::vector<Handle(Geom_BSplineCurve)>& isoCurves, int perpendDegree) {
	auto ctrl_point_matrix = ExtractControlPoints(isoCurves);
	auto point_params = ChordLenParam(ctrl_point_matrix);
	auto perpend_knot_sequence = CalKnots(point_params, perpendDegree);
	std::vector<Handle(Geom_BSplineCurve)> interpolate_curves(ctrl_point_matrix.size());
	for (int i = 0; i < ctrl_point_matrix.size(); i++) {
		auto curve = Interpolate(ctrl_point_matrix[i], point_params, perpend_knot_sequence, perpendDegree);
		interpolate_curves[i] = curve;
	}
	TColgp_Array2OfPnt Poles(1, ctrl_point_matrix.size(), 1, ctrl_point_matrix[0].size());
	for (int i = 0; i < ctrl_point_matrix.size(); i++) {
		for (int j = 0; j < ctrl_point_matrix[0].size(); j++) {
			Poles(i + 1, j + 1) = interpolate_curves[i]->Pole(j + 1);
		}
	}
	auto UKnots = isoCurves[0]->Knots();
	auto UMulti = isoCurves[0]->Multiplicities();
	auto VKnots = interpolate_curves[0]->Knots();
	auto VMulti = interpolate_curves[0]->Multiplicities();

	Handle(Geom_BSplineSurface) bspline = new Geom_BSplineSurface(Poles, UKnots, VKnots, UMulti, VMulti, isoCurves[0]->Degree(), perpendDegree);
	return bspline;

}

std::vector<double> InterPolateTool::CalKnots(const std::vector<double>& params, int degree) {
	std::vector<double> knots(params.size() + degree + 1);
	for (int i = 0; i <= degree; i++) {
		knots[i] = 0.;
	}
	for (int j = 1; j <= params.size() - degree - 1; j++) {
		double knot_temp = 0.;
		for (int i = j; i < j + degree; i++) {
			knot_temp += params[i];
		}
		knots[degree + j] = (1. / degree) * knot_temp;
	}
	for (int i = params.size(); i < knots.size(); i++) {
		knots[i] = 1.;
	}
	return knots;
}

std::vector<double> InterPolateTool::ChordLenParam(const std::vector<std::vector<gp_Pnt>>& PointMatrix) {
	int iso_curve_ctrl_point_num = PointMatrix.size();
	int iso_curve_num = PointMatrix[0].size();
	std::vector<double> params(iso_curve_num);
	params[0] = 0.;
	for (int j = 0; j < iso_curve_num - 1; j++) {
		double chord_length = 0.;
		for (int i = 0; i < iso_curve_ctrl_point_num; i++) {
			chord_length += PointMatrix[i][j].Distance(PointMatrix[i][j + 1]);
		}
		chord_length /= iso_curve_ctrl_point_num;
		params[j + 1] = params[j] + chord_length;
	}
	for (int i = 1; i < iso_curve_num; i++) {
		params[i] /= params[iso_curve_num - 1];
	}
	return params;
}

std::vector<std::vector<gp_Pnt>> InterPolateTool::ExtractControlPoints(const std::vector<Handle(Geom_BSplineCurve)>& bsplineCurves) {
	if (bsplineCurves.empty()) {
		throw std::invalid_argument("输入的 B 样条曲线向量为空。");
	}

	Standard_Integer numPoles = bsplineCurves[0]->NbPoles();

	for (size_t i = 1; i < bsplineCurves.size(); ++i) {
		int n = bsplineCurves[i]->NbPoles();
		if (bsplineCurves[i]->NbPoles() != numPoles) {
			throw std::invalid_argument("所有 B 样条曲线必须具有相同数量的控制点。");
		}
	}

	std::vector<std::vector<gp_Pnt>> controlPointMatrix(numPoles);

	for (Standard_Integer i = 1; i <= numPoles; ++i) {
		std::vector<gp_Pnt> controlPointsRow;
		for (const auto& curve : bsplineCurves) {
			gp_Pnt point = curve->Pole(i);
			controlPointsRow.push_back(point);
		}
		controlPointMatrix[i - 1] = controlPointsRow;
	}

	return controlPointMatrix;
}

double InterPolateTool::OneBasicFun(double u, int i, int p, const std::vector<double>& Knots)
{
	double Nip, uleft, uright, saved, temp;
	int m = Knots.size() - 1;
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

void InterPolateTool::sequenceToKnots(const std::vector<double>& sequence, std::vector<double>& knots, std::vector<int>& multiplicities) {
	if (sequence.empty()) return;

	std::map<double, int> knotMap;

	// 使用map来统计每个节点的重复次数
	for (double value : sequence) {
		bool found = false;
		for (auto& knot : knotMap) {
			if (isEqual(value, knot.first)) {
				knot.second++;
				found = true;
				break;
			}
		}
		if (!found) {
			knotMap[value] = 1;
		}
	}

	// 将map的内容转移到knots和multiplicities向量
	for (const auto& knot : knotMap) {
		knots.push_back(knot.first);
		multiplicities.push_back(knot.second);
	}
}

std::vector<double> InterPolateTool::knotsToSequence(const std::vector<double>& knots, const std::vector<int>& multiplicities) {
	std::vector<double> sequence;
	for (size_t i = 0; i < knots.size(); i++) {
		for (int j = 0; j < multiplicities[i]; j++) {
			sequence.push_back(knots[i]);
		}
	}
	return sequence;
}

Eigen::MatrixXd InterPolateTool::composeTensowMat(const Eigen::MatrixXd& matA, const Eigen::MatrixXd& matB) {
	int col_num1, col_num2, row_num1, row_num2;
	col_num1 = matA.cols();
	row_num1 = matA.rows();
	col_num2 = matB.cols();
	row_num2 = matB.rows();
	if (row_num1 != row_num2) {
		std::cerr << "the row number of matA must equal to row number of matB!";
	}
	Eigen::MatrixXd desMat(row_num1, col_num1 * col_num2);
	for (int rowLoopIndex = 0; rowLoopIndex < row_num1; rowLoopIndex++) {
		int index = 0;
		for (int i = 0; i < col_num1; i++) {
			for (int j = 0; j < col_num2; j++) {
				desMat(rowLoopIndex, index++) = matA(rowLoopIndex, i) * matB(rowLoopIndex, j);
			}
		}
	}
	return desMat;
}

TColgp_Array2OfPnt InterPolateTool::vectorToOCCMatrix(const std::vector<gp_Pnt>& pnts, int m) {
	int nRows = (pnts.size() + m - 1) / m;  // 计算需要的行数
	TColgp_Array2OfPnt array2D(1, nRows, 1, m);  // 创建二维数组

	int row = 1, col = 1;
	for (const auto& point : pnts) {
		array2D(row, col) = point;
		col++;
		if (col > m) {
			col = 1;
			row++;
		}
	}

	return array2D;
}

Handle(Geom_BSplineSurface) InterPolateTool::LoftV(const std::vector<Handle(Geom_BSplineCurve)>& isoCurves, int perpendDegree)
{
	auto ctrl_point_matrix = ExtractControlPoints(isoCurves);
	auto point_params = ChordLenParam(ctrl_point_matrix);
	auto perpend_knot_sequence = CalKnots(point_params, perpendDegree);
	std::vector<Handle(Geom_BSplineCurve)> interpolate_curves(ctrl_point_matrix.size());
	for (int i = 0; i < ctrl_point_matrix.size(); i++) {
		auto curve = Interpolate(ctrl_point_matrix[i], point_params, perpend_knot_sequence, perpendDegree);
		interpolate_curves[i] = curve;
	}
	TColgp_Array2OfPnt Poles(1, ctrl_point_matrix[0].size(), 1, ctrl_point_matrix.size());
	for (int i = 0; i < ctrl_point_matrix.size(); i++) {
		for (int j = 0; j < ctrl_point_matrix[0].size(); j++) {
			Poles(j + 1, i + 1) = interpolate_curves[i]->Pole(j + 1);
		}
	}
	auto VKnots = isoCurves[0]->Knots();
	auto VMulti = isoCurves[0]->Multiplicities();
	auto UKnots = interpolate_curves[0]->Knots();
	auto UMulti = interpolate_curves[0]->Multiplicities();


	Handle(Geom_BSplineSurface) bspline = new Geom_BSplineSurface(Poles, UKnots, VKnots, UMulti, VMulti, perpendDegree, isoCurves[0]->Degree());
	return bspline;
}
// 比较两个浮点数是否相等
bool InterPolateTool::isEqual(double x, double y, double epsilon) {
	return std::fabs(x - y) < epsilon;
}

// 比较 x 是否大于 y
bool InterPolateTool::isGreaterThan(double x, double y, double epsilon) {
	return (x - y) > epsilon;
}

// 比较 x 是否小于 y
bool InterPolateTool::isLessThan(double x, double y, double epsilon) {
	return (y - x) > epsilon;
}

// 比较 x 是否大于等于 y
bool InterPolateTool::isGreaterThanOrEqual(double x, double y, double epsilon) {
	return (x - y) > -epsilon;
}

// 比较 x 是否小于等于 y
bool InterPolateTool::isLessThanOrEqual(double x, double y, double epsilon) {
	return (y - x) > -epsilon;
}