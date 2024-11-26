/*
*    Copyright (c) 2018 Shing Liu All Rights Reserved.
*
*           File : occQt.cpp
*         Author : Shing Liu(eryar@163.com)
*           Date : 2018-01-08 21:00
*        Version : OpenCASCADE7.2.0 & Qt5.7.1
*
*    Description : Qt main window for OpenCASCADE.
*/

#include "occQt.h"
#include "occView.h"

#include <QToolBar>
#include <QTreeView>
#include <QMessageBox>
#include <QDockWidget>
#include <QProgressBar>

#include <gp_Circ.hxx>
#include <gp_Elips.hxx>
#include <gp_Pln.hxx>

#include <gp_Lin2d.hxx>

#include <Geom_ConicalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>

#include <GCE2d_MakeSegment.hxx>

#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <Geom_BSplineCurve.hxx>

#include <BRepLib.hxx>

#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>

#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>

#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>

#include <iostream>
#include <chrono>
#include <thread>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\ModelTriangulation.h"
#include <AIS_Shape.hxx>
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qfiledialog.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/STEPControl/STEPControl_Reader.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/BRepAlgoAPI/BRepAlgoAPI_Section.hxx"
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\ModelImporter.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qprogressdialog.h"
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\RandomExport.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qinputdialog.h"
#include <BRepTools.hxx>
#include <Geom_Hyperbola.hxx>
#include <Geom_TrimmedCurve.hxx>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/BRepClass3d/BRepClass3d_SolidClassifier.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/BRepExtrema/BRepExtrema_DistShapeShape.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qboxlayout.h"
#include "D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\MakeShape.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/testPoly/ModelExporter.h"
#include <future>
#include <mutex>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomLib/GeomLib.hxx"
#include <GeomAPI_PointsToBSpline.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <IGESControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <vector>
#include <limits>
#include <iostream>
#include "../src/IGESCAFControl/IGESCAFControl_Reader.hxx"
#include "../src/BRepBuilderAPI/BRepBuilderAPI_MakeVertex.hxx"
#include "../inc/BRepAdaptor_Curve.hxx"
#include "../src/GCPnts/GCPnts_TangentialDeflection.hxx"
#include "../inc/BRepAdaptor_Surface.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "../inc/STEPControl_Reader.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "../src/GeomLib/GeomLib.hxx"
#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>

#include <math_DoubleTab.hxx>
#include <math_Vector.hxx>
#include <Standard_OStream.hxx>
#include <math_Matrix.hxx>
#include "../src/Convert/Convert_ParameterisationType.hxx"
#include "../src/TColgp/TColgp_Array1OfXYZ.hxx"
#include "../src/GeomConvert/GeomConvert_CompCurveToBSplineCurve.hxx"
#include "../src/GeomLib/GeomLib.cxx"
#include <TColgp_Array1OfXYZ.hxx>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomAPI/GeomAPI_ExtremaCurveCurve.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomAPI/GeomAPI_ProjectPointOnCurve.hxx"
#include "SurfaceModelingTool.h"
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GCPnts/GCPnts_QuasiUniformAbscissa.hxx"
#include <GCPnts_AbscissaPoint.hxx>

// Thread
#include <thread>
#include <chrono>
#include <atomic>
#include "../../OpenCASCADE-7.7.0-vc14-64/opencascade-7.7.0/src/GeomFill/GeomFill_Coons.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include <STEPControl_Writer.hxx>
#include <GeomFill_BSplineCurves.hxx>

template <typename T>
void Visualize(const std::vector<T>& objects, Handle(AIS_InteractiveContext) context, OccView* myOccView, const Quantity_Color& color = Quantity_NOC_BISQUE)
{
    for (const auto& obj : objects)
    {
        try
        {
            Handle(AIS_Shape) aisShape;

            // ���ݶ������ʹ�����Ӧ�� AIS_Shape
            if constexpr (std::is_same<T, gp_Pnt>::value) 
            {
                TopoDS_Vertex ver = BRepBuilderAPI_MakeVertex(obj);
                aisShape = new AIS_Shape(ver);
            }
            else if constexpr (std::is_same<T, Handle(Geom_BSplineCurve)>::value) 
            {
                TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(obj);
                aisShape = new AIS_Shape(edge);
            }
            else if constexpr (std::is_same<T, Handle(Geom_BSplineSurface)>::value)
            {
                Handle(Geom_Surface) genericSurface = Handle(Geom_Surface)::DownCast(obj);
                if (genericSurface.IsNull())
                {
                    continue;
                }
                TopoDS_Face face = BRepBuilderAPI_MakeFace(genericSurface, Precision::Confusion());
                aisShape = new AIS_Shape(face);
            }
            else if constexpr (std::is_same<T, TopoDS_Shape>::value) 
            {
                aisShape = new AIS_Shape(obj);
            }
            else if constexpr (std::is_same<T, TopoDS_Edge>::value) 
            {
                aisShape = new AIS_Shape(obj);
            }

            aisShape->SetColor(color); // ������ɫ
            context->Display(aisShape, Standard_True);
        }
        catch (Standard_Failure& e)
        {
            QMessageBox::critical(nullptr, "Error", QString("Failed to display object: %1").arg(e.GetMessageString()));
            return;
        }
    }

    // �����ӽ�
    myOccView->fitAll();
}


void ExportBSplineSurface(const Handle(Geom_BSplineSurface)& bsplineSurface, const std::string& filename)
{
    // ��ȡ����Ĳ�����Χ
    Standard_Real uMin, uMax, vMin, vMax;
    bsplineSurface->Bounds(uMin, uMax, vMin, vMax);

    // ʹ������Ͳ�����Χ������
    TopoDS_Face face = BRepBuilderAPI_MakeFace(bsplineSurface, uMin, uMax, vMin, vMax, 1e-7);

    // ������Ƿ���Ч
    if (face.IsNull())
    {
        std::cerr << "�洴��ʧ�ܣ�" << std::endl;
        return;
    }

    // ���ļ���ת��ΪСд�Խ��в����ִ�Сд�ıȽ�
    std::string filename_lower = filename;
    std::transform(filename_lower.begin(), filename_lower.end(), filename_lower.begin(), ::tolower);

    // �����ļ���չ��ѡ�������ʽ
    if (filename_lower.size() >= 5 && filename_lower.substr(filename_lower.size() - 5) == ".brep")
    {
        // ���汣�浽 BREP �ļ�
        if (BRepTools::Write(face, filename.c_str()))
        {
            std::cout << "�ɹ������� BREP �ļ�: " << filename << std::endl;
        }
        else
        {
            std::cerr << "���� BREP �ļ�ʧ�ܣ�" << std::endl;
        }
    }
    else if (filename_lower.size() >= 5 && filename_lower.substr(filename_lower.size() - 5) == ".step")
    {
        // ���汣�浽 STEP �ļ�
        STEPControl_Writer writer;
        IFSelect_ReturnStatus status = writer.Transfer(face, STEPControl_AsIs);
        if (status == IFSelect_RetDone)
        {
            status = writer.Write(filename.c_str());
            if (status == IFSelect_RetDone)
            {
                std::cout << "�ɹ������� STEP �ļ�: " << filename << std::endl;
            }
            else
            {
                std::cerr << "���� STEP �ļ�ʧ�ܣ�" << std::endl;
            }
        }
        else
        {
            std::cerr << "��ת��Ϊ STEP ��ʽʧ�ܣ�" << std::endl;
        }
    }
    else
    {
        std::cerr << "��֧�ֵ��ļ���չ������ʹ�� .brep �� .step��" << std::endl;
    }
}


void PrintMatrix(const math_Matrix& Mat) {
    for (Standard_Integer i = Mat.LowerRow(); i <= Mat.UpperRow(); ++i) {
        for (Standard_Integer j = Mat.LowerCol(); j <= Mat.UpperCol(); ++j) {
            std::cout << Mat(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "----------------------------" << std::endl;
}

void PrintArray1OfXYZ(const TColgp_Array1OfXYZ& array)
{
    for (Standard_Integer i = array.Lower(); i <= array.Upper(); ++i)
    {
        const gp_XYZ& point = array.Value(i);
        std::cout << "Point " << i << ": ("
            << point.X() << ", "
            << point.Y() << ", "
            << point.Z() << ")" << std::endl;
    }
    std::cout << "---------------------" << std::endl;
}
void PrintArray1OfPnt(const TColgp_Array1OfPnt& array)
{
    for (Standard_Integer i = array.Lower(); i <= array.Upper(); ++i)
    {
        const gp_Pnt& point = array.Value(i);
        std::cout << "Point " << i << ": ("
            << point.X() << ", "
            << point.Y() << ", "
            << point.Z() << ")" << std::endl;
    }
    std::cout << "---------------------" << std::endl;
}
void M_ExtendCurveToPoint(Handle(Geom_BoundedCurve)& Curve,
    const gp_Pnt& Point,
    const Standard_Integer Continuity,
    const Standard_Boolean After)
{
    if (Continuity < 1 || Continuity > 3) return;
    Standard_Integer size = Continuity + 2;
    Standard_Real Ubord, Tol = 1.e-6;
    math_Matrix  MatCoefs(1, size, 1, size);
    Standard_Real Lambda, L1;
    Standard_Integer ii, jj;
    gp_Vec d1, d2, d3;
    gp_Pnt p0;
    // il faut Convertir l'entree (en preservant si possible le parametrage)
    GeomConvert_CompCurveToBSplineCurve Concat(Curve, Convert_QuasiAngular);

    // Les contraintes de constructions
    TColgp_Array1OfXYZ Cont(1, size);
    if (After) {
        Ubord = Curve->LastParameter();

    }
    else {
        Ubord = Curve->FirstParameter();
    }
    PLib::HermiteCoefficients(0, 1,           // Les Bornes
        Continuity, 0,  // Les Ordres de contraintes
        MatCoefs);
    PrintMatrix(MatCoefs);

    Curve->D3(Ubord, p0, d1, d2, d3);
    if (!After) { // Inversion du parametrage
        d1 *= -1;
        d3 *= -1;
    }

    L1 = p0.Distance(Point);
    if (L1 > Tol) {
        // Lambda est le ratio qu'il faut appliquer a la derive de la courbe
        // pour obtenir la derive du prolongement (fixe arbitrairement a la
        // longueur du segment bout de la courbe - point cible.
        // On essai d'avoir sur le prolongement la vitesse moyenne que l'on
        // a sur la courbe.
        gp_Vec daux;
        gp_Pnt pp;
        Standard_Real f = Curve->FirstParameter(), t, dt, norm;
        dt = (Curve->LastParameter() - f) / 9;
        norm = d1.Magnitude();
        for (ii = 1, t = f + dt; ii <= 8; ii++, t += dt) {
            Curve->D1(t, pp, daux);
            norm += daux.Magnitude();
        }
        norm /= 9;
        dt = d1.Magnitude() / norm;
        if ((dt < 1.5) && (dt > 0.75)) { // Le bord est dans la moyenne on le garde
            Lambda = ((Standard_Real)1) / Max(d1.Magnitude() / L1, Tol);
        }
        else {
            Lambda = ((Standard_Real)1) / Max(norm / L1, Tol);
        }
    }
    else {
        return; // Pas d'extension
    }

    // Optimisation du Lambda
    math_Matrix Cons(1, 3, 1, size);
    Cons(1, 1) = p0.X();  Cons(2, 1) = p0.Y(); Cons(3, 1) = p0.Z();
    Cons(1, 2) = d1.X();  Cons(2, 2) = d1.Y(); Cons(3, 2) = d1.Z();
    Cons(1, size) = Point.X();  Cons(2, size) = Point.Y(); Cons(3, size) = Point.Z();
    if (Continuity >= 2) {
        Cons(1, 3) = d2.X();  Cons(2, 3) = d2.Y(); Cons(3, 3) = d2.Z();
    }
    if (Continuity >= 3) {
        Cons(1, 4) = d3.X();  Cons(2, 4) = d3.Y(); Cons(3, 4) = d3.Z();
    }
    PrintMatrix(MatCoefs);
    // Construction dans la Base Polynomiale
    Cont(1) = p0.XYZ();
    Cont(2) = d1.XYZ() * Lambda;
    if (Continuity >= 2) Cont(3) = d2.XYZ() * Pow(Lambda, 2);
    if (Continuity >= 3) Cont(4) = d3.XYZ() * Pow(Lambda, 3);
    Cont(size) = Point.XYZ();
    PrintArray1OfXYZ(Cont);

    TColgp_Array1OfPnt ExtrapPoles(1, size);
    TColgp_Array1OfPnt ExtraCoeffs(1, size);

    gp_Pnt PNull(0., 0., 0.);
    ExtraCoeffs.Init(PNull);
    for (ii = 1; ii <= size; ii++) {
        for (jj = 1; jj <= size; jj++) {
            ExtraCoeffs(jj).ChangeCoord() += MatCoefs(ii, jj) * Cont(ii);
        }
    }
    PrintArray1OfPnt(ExtraCoeffs);
    // Convertion Dans la Base de Bernstein
    PLib::CoefficientsPoles(ExtraCoeffs, PLib::NoWeights(),
        ExtrapPoles, PLib::NoWeights());


    PrintArray1OfPnt(ExtrapPoles);
    Handle(Geom_BezierCurve) Bezier = new (Geom_BezierCurve) (ExtrapPoles);

    Standard_Real dist = ExtrapPoles(1).Distance(p0);
    Standard_Boolean Ok;
    Tol += dist;

    // Concatenation
    Ok = Concat.Add(Bezier, Tol, After);
    if (!Ok) throw Standard_ConstructionError("ExtendCurveToPoint");

    Curve = Concat.BSplineCurve();
}
gp_Ax2 getCoordinateSystemFromUserInput(bool ok)
{
    double xOrigin = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for origin:", 0, -1000, 1000, 8, &ok);
    double yOrigin = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for origin:", 0, -1000, 1000, 8, &ok);
    double zOrigin = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for origin:", 0, -1000, 1000, 8, &ok);

    double zDirX = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for Z direction:", 0, -1000, 1000, 8, &ok);
    double zDirY = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for Z direction:", 0, -1000, 1000, 8, &ok);
    double zDirZ = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for Z direction:", 0, -1000, 1000, 8, &ok);

    double xDirX = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for X direction:", 0, -1000, 1000, 8, &ok);
    double xDirY = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for X direction:", 0, -1000, 1000, 8, &ok);
    double xDirZ = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for X direction:", 0, -1000, 1000, 8, &ok);

    // Create gp_Ax2 with user-provided values
    gp_Ax2 coordinateSystemAx2 = gp_Ax2(gp_Pnt(xOrigin, yOrigin, zOrigin), gp_Dir(zDirX, zDirY, zDirZ), gp_Dir(xDirX, xDirY, xDirZ));

    return coordinateSystemAx2;
}
// Function to create different types of curves
TopoDS_Shape createCurve(const QString& curveType)
{
    bool ok;
    if (curveType.toLower() == "circle")
    {
        // Create a circle here
        // ...
        // Create gp_Ax2 with user-provided values
        gp_Ax2 circleAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);

        gp_Circ cir(circleAx2, radius);
        BRepBuilderAPI_MakeEdge edgeMaker(cir);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }
    }
    else if (curveType.toLower() == "ellipse")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 ellipseAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double majorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Major Radius:", 0, -1000, 1000, 8, &ok);
        double minorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Minor Radius:", 0, -1000, 1000, 8, &ok);



        gp_Elips ellipse(ellipseAx2, majorRadius, minorRadius);
        BRepBuilderAPI_MakeEdge edgeMaker(ellipse);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }
    }
    else if (curveType.toLower() == "parabola")
    {
        // Create a parabola here
        // ...
                // Get user input for origin coordinates
        // Create gp_Ax2 with user-provided values
        gp_Ax2 parabolaAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double theFocal = QInputDialog::getDouble(nullptr, "Input", "Enter the Focal:", 0, -1000, 1000, 8, &ok);

        gp_Parab parabola(parabolaAx2, theFocal);
        BRepBuilderAPI_MakeEdge edgeMaker(parabola);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }

    }
    else if (curveType.toLower() == "hyperbola")
    {
        // Get user input for origin coordinates
        gp_Ax2 hyperbolaAx2 = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double majorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Major Radius:", 0, -1000, 1000, 8, &ok);
        double minorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Minor Radius:", 0, -1000, 1000, 8, &ok);

        // ���������ϴ���˫����
        gp_Hypr hyperbolaPositive(hyperbolaAx2, majorRadius, minorRadius);
        BRepBuilderAPI_MakeEdge edgeMakerPositive(hyperbolaPositive);

        // ���� y ����� 90�� ����ת���õ��������ϵ�˫����
        gp_Trsf rotationTrsf;
        gp_Ax1 axis(hyperbolaAx2.Location(), hyperbolaAx2.YDirection());

        rotationTrsf.SetRotation(axis, M_PI);  // M_PI ��Բ����
        gp_Hypr hyperbolaNegative = hyperbolaPositive.Transformed(rotationTrsf);
        BRepBuilderAPI_MakeEdge edgeMakerNegative(hyperbolaNegative);

        // �����������ϵ�˫���߷���һ��������״��
        TopoDS_Compound compound;
        BRep_Builder compoundBuilder;
        compoundBuilder.MakeCompound(compound);
        compoundBuilder.Add(compound, edgeMakerPositive.Shape());
        compoundBuilder.Add(compound, edgeMakerNegative.Shape());

        return compound;
    }
    else if (curveType.toLower() == "line")
    {
        // Create a line here
        // ...
        // Get user input for origin coordinates
        double xStart = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for Start Point:", 0, -1000, 1000, 8, &ok);
        double yStart = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for Start Point:", 0, -1000, 1000, 8, &ok);
        double zStart = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for Start Point:", 0, -1000, 1000, 8, &ok);

        // Get user input for Z direction coordinates
        double DirX = QInputDialog::getDouble(nullptr, "Input", "Enter X coordinate for direction:", 0, -1000, 1000, 8, &ok);
        double DirY = QInputDialog::getDouble(nullptr, "Input", "Enter Y coordinate for direction:", 0, -1000, 1000, 8, &ok);
        double DirZ = QInputDialog::getDouble(nullptr, "Input", "Enter Z coordinate for direction:", 0, -1000, 1000, 8, &ok);

        gp_Pnt startPoint(xStart, yStart, zStart);
        gp_Dir dir(DirX, DirY, DirZ);
        gp_Lin line(startPoint, dir);
        BRepBuilderAPI_MakeEdge edgeMaker(line);
        if (edgeMaker.IsDone())
        {
            return edgeMaker.Shape();
        }
    }

    // Return an empty shape if the curve type is not recognized
    return TopoDS_Shape();
}
TopoDS_Shape createSurface(const QString& surfaceType)
{
    bool ok;
    if (surfaceType.toLower() == "plane")
    {
        gp_Ax2 planeAxis = getCoordinateSystemFromUserInput(ok);
        double length = QInputDialog::getDouble(nullptr, "Input", "Enter length:", 0, -1000, 1000, 8, &ok);
        double width = QInputDialog::getDouble(nullptr, "Input", "Enter width:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakePlane(planeAxis.Location(), planeAxis.XDirection() ^ planeAxis.YDirection(), planeAxis.XDirection(), length, width);

    }
    else if (surfaceType.toLower() == "cylinder")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 cylinderAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);
        double height = QInputDialog::getDouble(nullptr, "Input", "Enter Height:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeCylinder(cylinderAxis.Location(), cylinderAxis.XDirection() ^ cylinderAxis.YDirection(), cylinderAxis.XDirection(), radius, height);
    }
    else if (surfaceType.toLower() == "sphere")
    {
        gp_Ax2 sphereAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeSphere(sphereAxis.Location(), sphereAxis.XDirection() ^ sphereAxis.YDirection(), sphereAxis.XDirection(), radius);

    }
    else if (surfaceType.toLower() == "cone")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 coneAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double radius = QInputDialog::getDouble(nullptr, "Input", "Enter Radius:", 0, -1000, 1000, 8, &ok);
        double height = QInputDialog::getDouble(nullptr, "Input", "Enter Height:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeCone(coneAxis.Location(), coneAxis.XDirection() ^ coneAxis.YDirection(), coneAxis.XDirection(), radius, height);
    }
    else if (surfaceType.toLower() == "torus")
    {
        // Create gp_Ax2 with user-provided values
        gp_Ax2 torusAxis = getCoordinateSystemFromUserInput(ok);
        // Get user input for major and minor radii
        double majorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Major Radius:", 0, -1000, 1000, 8, &ok);
        double minorRadius = QInputDialog::getDouble(nullptr, "Input", "Enter Minor Radius:", 0, -1000, 1000, 8, &ok);

        return MakeShape::MakeTorus(torusAxis.Location(), torusAxis.XDirection() ^ torusAxis.YDirection(), torusAxis.XDirection(), majorRadius, minorRadius);
    }

    // Return an empty shape if the curve type is not recognized
    return TopoDS_Shape();
}
occQt::occQt(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

    createActions();
    createMenus();
    createToolBars();

    myOccView = new OccView(this);
    setCentralWidget(myOccView);

}

occQt::~occQt()
{

}

void occQt::createActions( void )
{
    // File
    connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));


    // View
    connect(ui.actionZoom, SIGNAL(triggered()), myOccView, SLOT(zoom()));
    connect(ui.actionPan, SIGNAL(triggered()), myOccView, SLOT(pan()));
    connect(ui.actionRotate, SIGNAL(triggered()), myOccView, SLOT(rotate()));

    connect(ui.actionReset, SIGNAL(triggered()), myOccView, SLOT(reset()));
    connect(ui.actionFitAll, SIGNAL(triggered()), myOccView, SLOT(fitAll()));

    // Primitive
    connect(ui.actionBox, SIGNAL(triggered()), this, SLOT(importFile()));
    connect(ui.actionCone, SIGNAL(triggered()), this, SLOT(Triangulation()));
    connect(ui.actionSphere, SIGNAL(triggered()), this, SLOT(TriangulationIntersection()));
    connect(ui.actionCylinder, SIGNAL(triggered()), this, SLOT(ClearDisplay()));
    connect(ui.actionTorus, SIGNAL(triggered()), this, SLOT(PrintInfo()));

    // Modeling
    connect(ui.actionFillet, SIGNAL(triggered()), this, SLOT(RandomExport()));
    connect(ui.actionChamfer, SIGNAL(triggered()), this, SLOT(MakePoint()));
    connect(ui.actionExtrude, SIGNAL(triggered()), this, SLOT(MakeCurve()));
    connect(ui.actionRevolve, SIGNAL(triggered()), this, SLOT(MakeSurface()));
    connect(ui.actionLoft, SIGNAL(triggered()), this, SLOT(ExportFile()));
    connect(ui.actionExtendCurveToPoint, SIGNAL(triggered()), this, SLOT(GenerateIsoCurves()));

    // Help
    connect(ui.actionAbout, SIGNAL(triggered()), this, SLOT(about()));
}

void occQt::createMenus( void )
{
}

void occQt::createToolBars( void )
{
    QToolBar* aToolBar;
    //aToolBar = addToolBar(tr("&Navigate"));
    //aToolBar->addAction(ui.actionZoom);
    //aToolBar->addAction(ui.actionPan);
    //aToolBar->addAction(ui.actionRotate);

    //aToolBar = addToolBar(tr("&View"));
    //aToolBar->addAction(ui.actionReset);
    //aToolBar->addAction(ui.actionFitAll);

    aToolBar = addToolBar(tr("&Intersection"));
    aToolBar->addAction(ui.actionBox);
    aToolBar->addAction(ui.actionCone);
    aToolBar->addAction(ui.actionSphere);
    aToolBar->addAction(ui.actionCylinder);
    aToolBar->addAction(ui.actionTorus);

    aToolBar = addToolBar(tr("&Utils"));
    aToolBar->addAction(ui.actionFillet);
    aToolBar->addAction(ui.actionChamfer);
    aToolBar->addAction(ui.actionExtrude);
    aToolBar->addAction(ui.actionRevolve);
    aToolBar->addAction(ui.actionLoft);
    aToolBar->addAction(ui.actionExtendCurveToPoint);
    aToolBar->addSeparator();
    aToolBar->addSeparator();

    //aToolBar = addToolBar(tr("About"));
    //aToolBar->addAction(ui.actionAbout);
}

void occQt::about()
{
    QMessageBox::about(this, tr("About the Project"),
        tr("<h2>Opencascade Application 2.0</h2>"
        "<p>Copyright &copy; 281885807@qq.com"));
}

void occQt::importFile()
{
    myOccView->getContext()->EraseAll(Standard_True);
    TopoDS_Shape aTopoBox = BRepPrimAPI_MakeBox(3.0, 4.0, 5.0).Shape();
    TopoDS_Shape shape;
    STEPControl_Reader reader;
    QString filename = QFileDialog::getOpenFileName(this, "ѡ���ļ�", QDir::homePath(), "All Files (*.*);;Text Files (*.txt);;Image Files (*.png *.jpg)");
    if (reader.ReadFile(filename.toStdString().c_str()) == IFSelect_RetDone) {
        reader.TransferRoots();
        shape = reader.OneShape();
    }
    TopExp_Explorer exp;
    int cnt = 0;
    for (exp.Init(shape, TopAbs_FACE); exp.More(); exp.Next()) 
    {
        if (!cnt) {
            shape1 = exp.Current();
        }
        else
        {
            shape2 = exp.Current();
        }
        cnt++;
    }
    //Handle(AIS_Shape) anAisBox = new AIS_Shape(aTopoBox);
    Handle(AIS_Shape) aShape1 = new AIS_Shape(shape1);
    Handle(AIS_Shape) aShape2 = new AIS_Shape(shape2);

    //anAisBox->SetColor(Quantity_NOC_AZURE);
    aShape1->SetColor(Quantity_NOC_BISQUE);
    aShape2->SetColor(Quantity_NOC_BISQUE);

    myOccView->getContext()->Display(aShape1, Standard_True);
    myOccView->getContext()->Display(aShape2, Standard_True);
    intersectionDone = false;
}

void occQt::Triangulation()
{
    std::vector<TopoDS_Shape> shapes;
    int choice;

    if (shape1.IsNull() && shape2.IsNull())
    {
        QString filename = QFileDialog::getOpenFileName(this, "choose files", QDir::homePath(), "All Files (*.*);;Text Files (*.txt);;Image Files (*.png *.jpg)");
        ModelImporter importer;
        shapes = importer.loadStepShape(filename.toStdString());
        shape1 = shapes[0];
        shape2 = shapes[1];
    }
    else
    {
        shapes.push_back(shape1);
        shapes.push_back(shape2);
    }

    // Display shape1 and shape2
    myOccView->getContext()->RemoveAll(Standard_True);
    Handle(AIS_Shape) aShape1 = new AIS_Shape(shape1);
    aShape1->SetColor(Quantity_NOC_BISQUE);
    myOccView->getContext()->Display(aShape1, Standard_True);
    Handle(AIS_Shape) aShape2 = new AIS_Shape(shape2);
    aShape2->SetColor(Quantity_NOC_BISQUE);
    myOccView->getContext()->Display(aShape2, Standard_True);

    int cnt = 0;
    // Store triangulated results and build R-trees
    for (auto shape : shapes)
    {
        faceList[cnt] = ModelTriangulation::ModelTriangulate(shape, triPoints[cnt]);
        cnt++;
    }

    choice = QMessageBox::question(this, "Confirmation", "Display Results? ", QMessageBox::Yes | QMessageBox::No);
    if (choice == QMessageBox::Yes)
    {
        myOccView->getContext()->RemoveAll(Standard_True);
        std::vector<std::future<void>> futures;
        std::mutex displayMutex;

        for (int j = 0; j <= 1; j++)
        {
            for (size_t i = 0; i < faceList[j].size(); ++i)
            {
                futures.push_back(std::async(std::launch::async, [&, j, i]() 
                    {
                    Handle(AIS_Shape) triShape = new AIS_Shape(faceList[j][i]);
                    triShape->SetColor(Quantity_NOC_BISQUE);
                    // ������,����ʱ�Զ�����
                    std::lock_guard<std::mutex> lock(displayMutex);
                    myOccView->getContext()->Display(triShape, Standard_True);
                    }));
            }
        }

        // �ȴ������첽�������
        for (auto& f : futures) {
            f.get();
        }
    }
}


void occQt::TriangulationIntersection()
{
    if (faceList[0].size() && faceList[1].size())
    {
        auto start_time1 = std::chrono::high_resolution_clock::now();
        // ����R��
        R_Tree RTree1, RTree2;
        RTree1.BuildRtree(faceList[0], triPoints[0]);
        RTree2.BuildRtree(faceList[1], triPoints[1]);
        // ��ȡ����ʱ���

        auto end_time1 = std::chrono::high_resolution_clock::now();
        // �������ִ��ʱ�䣬��λΪ��
        std::chrono::duration<double> execution_time1 = end_time1 - start_time1;

        // ��ִ��ʱ��ת��Ϊ����
        double execution_time_ms1 = execution_time1.count() * 1000.0;

        std::cout << "����ִ��ʱ��Ϊ��" << execution_time_ms1 << " ����" << std::endl;
        // ��ʼ�����
        result.Nullify();
        BRep_Builder builder;
        builder.MakeCompound(result);

        std::vector<std::future<void>> futures;
        std::mutex sectionMutex; // ���ڱ���BRepAlgoAPI_Section����

        // ����һ���ݹ麯��������ÿ���ڵ�
        std::function<void(R_TreeNode*)> processNode;
        processNode = [&](R_TreeNode* currentNode) {
            // �ڵڶ��� R �����ҵ������ཻ�Ľڵ�
            std::vector<R_TreeNode*> potentialIntersectNodes = RTree2.Search(currentNode->Box, RTree2.Root);

            for (R_TreeNode* node2 : potentialIntersectNodes)
            {
                // ��ÿ���󽻲�����Ϊһ���첽����
                futures.push_back(std::async(std::launch::async, [&, currentNode, node2]() 
                    {
                    std::unique_lock<std::mutex> lock(sectionMutex);
                    BRepAlgoAPI_Section section(currentNode->Face, node2->Face, Standard_True);
                    section.ComputePCurveOn1(Standard_True);
                    section.Approximation(Standard_True);
                    section.Build();
                    lock.unlock(); // �ͷ���

                    if (section.IsDone())
                    {
                        TopoDS_Shape intersectionLine = section.Shape();
                        // �����ﲻ��Ҫ������Ϊbuilder.Add���ٶ�Ϊ�̰߳�ȫ��Խ����Ӱ�첻����ɳ�ͻ
                        builder.Add(result, intersectionLine);
                    }
                    }));
            }

            // ����ǰ�ڵ���ӽڵ����ݹ鴦��
            for (R_TreeNode* child : currentNode->Childs)
            {
                processNode(child);
            }
            };

        // ��ʼ�Ӹ��ڵ㴦��
        processNode(RTree1.Root);

        // �ȴ������첽�������
        for (auto& future : futures) {
            future.get();
        }

        Handle(AIS_Shape) aIntersectionLine = new AIS_Shape(result);
        aIntersectionLine->SetColor(Quantity_NOC_RED);
        myOccView->getContext()->Display(aIntersectionLine, Standard_True);
        faceList[0].clear();
        faceList[1].clear();
        intersectionDone = true;
    }
    else
    {
        Triangulation();
        TriangulationIntersection();
    }

}

void occQt::ClearDisplay()
{
    myOccView->getContext()->RemoveAll(Standard_True);
    shape1 = shape2 = TopoDS_Shape();
    curveShape = TopoDS_Shape();
    pointShape = TopoDS_Shape();
    faceList[0].clear();
    faceList[1].clear();
    result.Nullify();
    boundedCurve.Nullify();

}

void occQt::PrintInfo()
{
    // ���ı��������������Ϣ
    QString outputText;
    TopExp_Explorer exp;
    int cnt = 1;
    for (exp.Init(result, TopAbs_EDGE); exp.More(); exp.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(exp.Current());
        // ��ȡ�ߵ���ʼ�ͽ�������
        TopoDS_Vertex startVertex, endVertex;
        TopExp::Vertices(edge, startVertex, endVertex);

        // ��ȡ���������
        gp_Pnt startPoint = BRep_Tool::Pnt(startVertex);
        gp_Pnt endPoint = BRep_Tool::Pnt(endVertex);

        // ��������Ϣ��ӵ�����ı�
        outputText += QString("Edge%1: startV : (%2, %3, %4), endV : (%5, %6, %7)\n")
            .arg(cnt).arg(startPoint.X()).arg(startPoint.Y()).arg(startPoint.Z())
            .arg(endPoint.X()).arg(endPoint.Y()).arg(endPoint.Z());
        cnt++;
    }

    // �ڶԻ�������ʾ�ı���
    QDialog dialog(this);
    QVBoxLayout* layout = new QVBoxLayout(&dialog);


    QTextEdit* outputTextEdit = new QTextEdit(&dialog);
    outputTextEdit->setPlainText(outputText);
    // �����ı����������
    layout->addWidget(outputTextEdit);

    // �����ı�����չ�����ռ�
    outputTextEdit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    dialog.setWindowTitle("Edge Info");

    // ��ʾ�Ի���
    dialog.exec();

}

void occQt::RandomExport()
{

    // ����һ���ļ��Ի���
    QFileDialog dialog;

    // ���öԻ���ģʽΪѡ���ļ���
    dialog.setFileMode(QFileDialog::Directory);
    QString folderPath;
    // ��ʾ�Ի���
    if (dialog.exec()) {
        // ��ȡ�û�ѡ����ļ���·��
        folderPath = dialog.selectedFiles().first();
    }
    // ��ȡ�������ַ���
    bool ok;
    QString string1 = QInputDialog::getText(nullptr, "Input", "Enter String 1: (Basic Model Type: plane, cone, cylinder, sphere, torus)", QLineEdit::Normal, "", &ok);
    QString string2 = QInputDialog::getText(nullptr, "Input", "Enter String 2: (Basic Model Type: plane, cone, cylinder, sphere, torus)", QLineEdit::Normal, "", &ok);
    int value = QInputDialog::getInt(nullptr, "Input", "Enter an Integer:", 0, 0, 100, 1, &ok);
     

    // ʹ���û��������RandomExport::randomRotateAndExport()
    RandomExport::randomRotateAndExport(folderPath.toStdString().c_str(), string1.toStdString().c_str(), string2.toStdString().c_str(), value);
}


void occQt::MakePoint()
{
    bool ok;

    double x = QInputDialog::getDouble(nullptr, "Input", "Enter X for point:", 0, -1000, 1000, 8, &ok);
    double y = QInputDialog::getDouble(nullptr, "Input", "Enter Y for point:", 0, -1000, 1000, 8, &ok);
    double z = QInputDialog::getDouble(nullptr, "Input", "Enter Z for point:", 0, -1000, 1000, 8, &ok);

    TopoDS_Shape point = MakeShape::MakePoint(x, y, z);
    if (!point.IsNull())
    {
        Handle(AIS_Shape) aPoint = new AIS_Shape(point);
        aPoint->SetColor(Quantity_NOC_BISQUE);
        myOccView->getContext()->Display(aPoint, Standard_True);
    }
}


void occQt::MakeCurve()
{
    bool ok;
    QString curveType = QInputDialog::getText(nullptr, "Input", "Enter Curve Type: (circle, ellipse, parabola, hyperbola, line)", QLineEdit::Normal, "", &ok);
    if (ok)
        curveShape = createCurve(curveType);

    if (!curveShape.IsNull())
    {
        Handle(AIS_Shape) aCurve = new AIS_Shape(curveShape);
        aCurve->SetColor(Quantity_NOC_BISQUE);
        myOccView->getContext()->Display(aCurve, Standard_True);
    }
}

void occQt::MakeSurface()
{
    bool ok;
    TopoDS_Shape shape;
    QString surfaceType = QInputDialog::getText(nullptr, "Input", "Enter Surface Type: (plane, cylinder, sphere, cone, torus)", QLineEdit::Normal, "", &ok);
    if (ok)
    {
        shape = createSurface(surfaceType);
        if (changeCount % 2 == 0)
        {
            shape1 = shape;
        }
        else 
        {
            shape2 = shape;
        }
        changeCount++;
    }

    if (!shape.IsNull())
    {
        Handle(AIS_Shape) aShape = new AIS_Shape(shape);
        aShape->SetColor(Quantity_NOC_BISQUE);
        myOccView->getContext()->Display(aShape, Standard_True);
    }
    intersectionDone = false;
}

void occQt::ExportFile()
{
    // ����һ���ļ��Ի���
    QFileDialog dialog;

    // ���öԻ���ģʽΪѡ���ļ���
    dialog.setFileMode(QFileDialog::Directory);
    QString folderPath;
    // ��ʾ�Ի���
    if (dialog.exec()) {
        // ��ȡ�û�ѡ����ļ���·��
        folderPath = dialog.selectedFiles().first();
    }

    // ָ���ļ���
    QString filePath = folderPath + "/exportModel.step";

    // ��ȡ�ļ��Ļ�����Ϣ
    QFileInfo fileInfo(filePath);

    // ����ļ��Ѵ��ڣ����Զ���Ӻ�׺
    int suffix = 1;
    while (fileInfo.exists()) {
        filePath = folderPath + QString("/exportModel_%1.step").arg(suffix);
        fileInfo.setFile(filePath);
        suffix++;
    }

    BRep_Builder builder;
    TopoDS_Compound exportResult;
    builder.MakeCompound(exportResult);
    if (intersectionDone)
    {
        builder.Add(exportResult, result);
    }
    builder.Add(exportResult, shape1);
    builder.Add(exportResult, shape2);

    ModelExporter exporter;
    exporter.exportModel(exportResult, filePath.toStdString());

    exportResult.Nullify();

}

// �Ե����Ȳ��߽��еȾ���������㻡��
std::vector<std::pair<gp_Pnt, double>> SampleIsoCurveWithArcLength(const Handle(Geom_BSplineCurve)& bsplineCurve, int numSamples) 
{
    Handle(GeomAdaptor_Curve) curve = new GeomAdaptor_Curve(bsplineCurve);
    GCPnts_QuasiUniformAbscissa sampler(*curve, numSamples);
    std::vector<std::pair<gp_Pnt, double>> sampledPointsWithArcLength;

    for (int i = 1; i <= sampler.NbPoints(); ++i)
    {
        Standard_Real param = sampler.Parameter(i);
        gp_Pnt point = curve->Value(param); // ��ȡ������

        // ���㵱ǰ�����һ����֮��Ļ���
        Standard_Real arcLength = CPnts_AbscissaPoint::Length(*curve, sampler.Parameter(1), param);

        // ���µ�Ͷ�Ӧ�����洢
        sampledPointsWithArcLength.emplace_back(point, arcLength);
    }

    return sampledPointsWithArcLength;
}
// ���������Ͷ�Ӧ�����ĺ���
void PrintSampledPointsWithArcLength(const std::vector<std::pair<gp_Pnt, double>>& sampledPoints) {
    for (const auto& pair : sampledPoints)
    {
        const gp_Pnt& point = pair.first;
        double arcLength = pair.second;
        std::cout << "Point: (" << point.X() << ", " << point.Y() << ", " << point.Z() << "), Arc Length: " << arcLength << std::endl;
    }
}

Handle(Geom_BSplineSurface) GenerateCoonsSurface(
    Handle(Geom_BSplineCurve)& curve1, Handle(Geom_BSplineCurve)& curve2, Handle(Geom_BSplineCurve)& curve3, Handle(Geom_BSplineCurve)& curve4
) {

    // ���� GeomFill_BSplineCurves ����ʹ�� Coons �����ʽ
    GeomFill_BSplineCurves fillCoonsStyle(
        curve1, // U=0
        curve2, // V=1
        curve3, // U=1
        curve4, // V=0
        GeomFill_CoonsStyle
    );

    // ��ȡ���ɵ�����
    Handle(Geom_Surface) surface = fillCoonsStyle.Surface();

    // ��������Ƿ����ɳɹ�
    if (surface.IsNull()) {
        throw std::runtime_error("ʹ�� CoonsStyle ��������ʧ�ܣ�");
    }

    std::cout << "ʹ�� CoonsStyle �ɹ��������棡" << std::endl;

    // ���Խ� Geom_Surface ת��Ϊ Geom_BSplineSurface
    Handle(Geom_BSplineSurface) bsplineSurface = Handle(Geom_BSplineSurface)::DownCast(surface);
    if (bsplineSurface.IsNull()) 
    {
        throw std::runtime_error("���ɵ����治�� B-Spline ���棡");
    }

    return bsplineSurface;
}

void occQt::GenerateIsoCurves(void)
{
    for (int i = 12; i <= 12; i++)
    {
        myOccView->getContext()->RemoveAll(Standard_True);
        // ����߽���
        std::vector<Handle(Geom_BSplineCurve)> tempArray;
        tempArray.clear();
        std::string filename = "E:/Models/Constraint/test";
        filename += std::to_string(i);
        filename += "/";
        std::string boundaryPath = filename + "boundary.brep";
        SurfaceModelingTool::LoadBSplineCurves(boundaryPath.c_str(), tempArray);
        Visualize(tempArray, myOccView->getContext(), myOccView);

        std::vector<std::vector<double>> uKnots;
        std::vector<std::vector<double>> vKnots;
        std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_Final, vISOcurvesArray_Final;

        SurfaceModelingTool::ApproximateBoundaryCurves(tempArray);
        if (tempArray.size() == 3)
        {
            gp_Pnt pnt1 = tempArray[0]->StartPoint(), pnt2 = tempArray[0]->EndPoint(), pnt3 = tempArray[1]->StartPoint();
            gp_Pnt pnt4 = tempArray[1]->EndPoint(), pnt5 = tempArray[2]->StartPoint(), pnt6 = tempArray[2]->EndPoint();

            double tol = tempArray[0]->EndPoint().Distance(tempArray[0]->StartPoint()) / 1000;
            if (tempArray[1]->EndPoint().Distance(tempArray[0]->EndPoint()) < tol)
            {
                tempArray[1]->Reverse();
            }
            else if (tempArray[2]->StartPoint().Distance(tempArray[0]->EndPoint()) < tol)
            {
                std::swap(tempArray[1], tempArray[2]);
            }
            else if (tempArray[2]->EndPoint().Distance(tempArray[0]->EndPoint()) < tol)
            {
                std::swap(tempArray[1], tempArray[2]);
                tempArray[1]->Reverse();
            }

            if (tempArray[2]->EndPoint().Distance(tempArray[1]->EndPoint()) < tol)
            {
                tempArray[2]->Reverse();
            }

            pnt1 = tempArray[0]->StartPoint(); pnt2 = tempArray[0]->EndPoint(); pnt3 = tempArray[1]->StartPoint();
            pnt4 = tempArray[1]->EndPoint(); pnt5 = tempArray[2]->StartPoint(); pnt6 = tempArray[2]->EndPoint();

            // ��������������˻��߹����ı�
            std::vector<gp_Pnt> boundaryPoints = { tempArray[0]->StartPoint(), tempArray[1]->StartPoint(), tempArray[2]->StartPoint() };

            // �����
            gp_Vec line_01(boundaryPoints[1].XYZ() - boundaryPoints[0].XYZ());
            gp_Vec line_12(boundaryPoints[2].XYZ() - boundaryPoints[1].XYZ());
            gp_Vec line_20(boundaryPoints[0].XYZ() - boundaryPoints[2].XYZ());

            auto calculateAngle = [](const gp_Vec& v1, const gp_Vec& v2)
            {
                double dotProduct = v1.Dot(v2);
                double magnitudes = v1.Magnitude() * v2.Magnitude();
                return std::acos(dotProduct / magnitudes);  // ���ؽǶ�
            };

            // ���������ǵļн�
            double angleAtPoint0 = calculateAngle(-line_20, line_01);  // ��0�ļн�
            double angleAtPoint1 = calculateAngle(-line_01, line_12);  // ��1�ļн�
            double angleAtPoint2 = calculateAngle(-line_12, line_20);  // ��2�ļн�

            // �ҳ����Ƕ�
            double maxAngle = std::max({ angleAtPoint0, angleAtPoint1, angleAtPoint2 });
            int maxAngleIndex = 0;
            if (maxAngle == angleAtPoint1) maxAngleIndex = 1;
            else if (maxAngle == angleAtPoint2) maxAngleIndex = 2;

            // �����˻���
            auto CreateDegenerateEdge = [](const gp_Pnt & point)
            {
                TColgp_Array1OfPnt poles(1, 2);
                poles.SetValue(1, point);
                poles.SetValue(2, point);

                TColStd_Array1OfReal knots(1, 2);
                knots.SetValue(1, 0.0);
                knots.SetValue(2, 1.0);

                TColStd_Array1OfInteger multiplicities(1, 2);
                multiplicities.SetValue(1, 2);
                multiplicities.SetValue(2, 2);

                return new Geom_BSplineCurve(poles, knots, multiplicities, 1);
            };

            //maxAngleIndex = 0;
            //maxAngleIndex = 1;
            //maxAngleIndex = 2;
            if (maxAngleIndex == 0)
            {
                tempArray.push_back(CreateDegenerateEdge(boundaryPoints[maxAngleIndex]));
            }
            else
            {
                tempArray.insert(tempArray.begin() + maxAngleIndex, CreateDegenerateEdge(boundaryPoints[maxAngleIndex]));
            }
        }

        Handle(Geom_BSplineCurve) bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4;
        SurfaceModelingTool::Arrange_Coons_G0(tempArray, bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4,10, Standard_True);
        // �洢�߽�����
        std::vector<Handle(Geom_BSplineCurve)> aBoundarycurveArray = { bslpineCurve1 , bslpineCurve2, bslpineCurve3, bslpineCurve4 };

        // ��ʾ�߽�����
        myOccView->getContext()->RemoveAll(Standard_True);
        Visualize(aBoundarycurveArray, myOccView->getContext(), myOccView);

        // ��ȡCoons����
        Handle(Geom_BSplineSurface) surfacecoons;
        surfacecoons = GenerateCoonsSurface(bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4);
        std::string SurfaceCoonsFilename = filename + "SurfaceCoons_OCC.step";
        ExportBSplineSurface(surfacecoons, SurfaceCoonsFilename);

        SurfaceModelingTool::Coons_G0(bslpineCurve1, bslpineCurve2, bslpineCurve3, bslpineCurve4, surfacecoons);
        SurfaceCoonsFilename = filename + "SurfaceCoons_y.step";
        ExportBSplineSurface(surfacecoons, SurfaceCoonsFilename);  
        // ��Coons�����ȡ��ʼ�Ȳ��ߣ����Ҽ���ÿ���Ȳ�������Ӧ�ķ���
        std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_Initial, vISOcurvesArray_Initial;
        std::vector<gp_Vec> normalsOfUISOLines, normalsOfVISOLines;
        int isoCount = 20;
        SurfaceModelingTool::GetISOCurveWithNormal(surfacecoons, uISOcurvesArray_Initial, vISOcurvesArray_Initial, normalsOfUISOLines, normalsOfVISOLines,isoCount);

        // ���ӻ��׶ν��
        //std::vector<Handle(Geom_BSplineSurface)> visualSurfaceArray = { surfacecoons };
        //VisualizeBSplineSurface(visualSurfaceArray, myOccView->getContext(), myOccView);
        Visualize(uISOcurvesArray_Initial, myOccView->getContext(), myOccView);
        Visualize(vISOcurvesArray_Initial, myOccView->getContext(), myOccView);

        //����Lofting����
        std::vector<TopoDS_Shape>  uLoftingSur, vLoftingSur;
        SurfaceModelingTool::CreateLoftingSurface(uISOcurvesArray_Initial, normalsOfUISOLines, uLoftingSur);
        SurfaceModelingTool::CreateLoftingSurface(vISOcurvesArray_Initial, normalsOfVISOLines, vLoftingSur);
        //VisualizeShapes(uLoftingSur, myOccView->getContext(), myOccView);
        //VisualizeShapes(vLoftingSur, myOccView->getContext(), myOccView);

        //�����ڲ���
        std::vector<Handle(Geom_BSplineCurve)> anInternalBSplineCurves;
        std::string internalPath = filename + "internal.brep";
        SurfaceModelingTool::LoadBSplineCurves(internalPath, anInternalBSplineCurves);
        //VisualizeBSplineCurves(anInternalBSplineCurves, myOccView->getContext(), myOccView, Quantity_NOC_RED);
        myOccView->getContext()->RemoveAll(Standard_True);

        // ���������ĵȲ���
        std::vector<Handle(Geom_BSplineCurve)> uISOcurvesArray_New, vISOcurvesArray_New;
        std::vector<gp_Pnt> interPoints;
        std::vector<std::vector<gp_Pnt>> uInterpolatePoints;
        std::vector<std::vector<gp_Pnt>> vInterpolatePoints;
        std::vector<TopoDS_Edge> uInterpoalteTangentArray;
        std::vector<TopoDS_Edge> uInterpoalteTangentArray2;
        SurfaceModelingTool::LoftSurfaceIntersectWithCurve(uLoftingSur, uISOcurvesArray_Initial, anInternalBSplineCurves, uISOcurvesArray_New, isoCount, uInterpolatePoints, uInterpoalteTangentArray, uInterpoalteTangentArray2,surfacecoons);
        Visualize(uInterpoalteTangentArray, myOccView->getContext(), myOccView, Quantity_NOC_GOLD);
        Visualize(uInterpoalteTangentArray2, myOccView->getContext(), myOccView, Quantity_NOC_GOLD);
        for (auto interPoints : uInterpolatePoints)
            Visualize(interPoints, myOccView->getContext(), myOccView, Quantity_NOC_RED);

        std::vector<TopoDS_Edge> vInterpoalteTangentArray;
        std::vector<TopoDS_Edge> vInterpoalteTangentArray2;
        SurfaceModelingTool::LoftSurfaceIntersectWithCurve(vLoftingSur, vISOcurvesArray_Initial, anInternalBSplineCurves, vISOcurvesArray_New, isoCount, vInterpolatePoints, vInterpoalteTangentArray, vInterpoalteTangentArray2,surfacecoons);
        Visualize(vInterpoalteTangentArray, myOccView->getContext(), myOccView, Quantity_NOC_GOLD);
        Visualize(vInterpoalteTangentArray2, myOccView->getContext(), myOccView, Quantity_NOC_GOLD);
        //VisualizeBSplineSurface({surfacecoons}, myOccView->getContext(), myOccView, Quantity_NOC_BLUE1);
        for(auto interPoints : vInterpolatePoints)
            Visualize(interPoints, myOccView->getContext(), myOccView, Quantity_NOC_RED);

        Visualize(aBoundarycurveArray, myOccView->getContext(), myOccView);
        Visualize(vISOcurvesArray_New, myOccView->getContext(), myOccView);
        Visualize(uISOcurvesArray_New, myOccView->getContext(), myOccView);

        // �����µȲ���
        {
            //std::this_thread::sleep_for(std::chrono::milliseconds(5000));
            myOccView->getContext()->RemoveAll(Standard_True);

            // ����u��v�Ȳ���֮��Ľ��㣬�������յȲ���
            interPoints.clear();
            std::vector<gp_Pnt> boundaryPoints;
            std::vector<Handle(Geom_BSplineSurface)> surfaceArray;
            std::vector<TopoDS_Edge> TangentArray;
            std::string gordenSurf1 = filename + "gordonSurf1.step";
            std::string gordenSurf2 = filename + "gordonSurf2.step";
            std::string gordenSurf3 = filename + "gordonSurf3.step";
            std::string gordenSurf4 = filename + "gordonSurf4.step";
            SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf1, surfaceArray);
            SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf2, surfaceArray);
            SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf3, surfaceArray);
            SurfaceModelingTool::LoadBSplineSurfaces(gordenSurf4, surfaceArray);
            Visualize(surfaceArray, myOccView->getContext(), myOccView);
            SurfaceModelingTool::CreateFinalISOCurves(uISOcurvesArray_New, vISOcurvesArray_New, uISOcurvesArray_Final, vISOcurvesArray_Final, uInterpolatePoints, vInterpolatePoints,uKnots, vKnots, boundaryPoints, interPoints, isoCount, TangentArray, surfaceArray);
            Visualize(TangentArray, myOccView->getContext(), myOccView, Quantity_NOC_RED);

            SurfaceModelingTool::UpdateFinalCurves(aBoundarycurveArray, uISOcurvesArray_Final, vISOcurvesArray_Final);
            
            for (auto boundaryPoint : boundaryPoints)
            {
                interPoints.push_back(boundaryPoint);
            }
            interPoints.push_back(uISOcurvesArray_Final[0]->StartPoint());
            interPoints.push_back(uISOcurvesArray_Final[0]->EndPoint());
            interPoints.push_back(uISOcurvesArray_Final[uISOcurvesArray_Final.size() - 1]->StartPoint());
            interPoints.push_back(uISOcurvesArray_Final[uISOcurvesArray_Final.size() - 1]->EndPoint());

            //VisualizePoints(interPoints, myOccView->getContext(), myOccView);
            // ���� u(v)ISOcurvesArray_Final ���п��ӻ�
            Visualize(uISOcurvesArray_Final, myOccView->getContext(), myOccView);
            Visualize(vISOcurvesArray_Final, myOccView->getContext(), myOccView);
            auto ExportPointsToBREP = [](const std::vector<gp_Pnt>& boundaryPoints, const std::string& filename)
                {
                    // ����һ���������Դ洢���ж���
                    TopoDS_Compound compound;
                    BRep_Builder builder;
                    builder.MakeCompound(compound);

                    // �� gp_Pnt ת��Ϊ TopoDS_Vertex ����ӵ�������
                    for (const auto& point : boundaryPoints)
                    {
                        TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(point);
                        builder.Add(compound, vertex);
                    }

                    // �������屣�浽 BREP �ļ�
                    if (BRepTools::Write(compound, filename.c_str()))
                    {
                        std::cout << "�ɹ������� BREP �ļ�: " << filename << std::endl;
                    }
                    else {
                        std::cerr << "���� BREP �ļ�ʧ�ܣ�" << std::endl;
                    }
                };

            ExportPointsToBREP(interPoints, filename + std::string("points.brep"));

            TopoDS_Compound UResult;
            BRep_Builder builder1;
            builder1.MakeCompound(UResult);
            for (auto curve : uISOcurvesArray_Final)
            {
                TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
                builder1.Add(UResult, edge);
            }
            TopoDS_Compound VResult;
            BRep_Builder builder;
            builder.MakeCompound(VResult);
            for (auto curve : vISOcurvesArray_Final)
            {
                TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
                builder.Add(VResult, edge);
            }

            std::string UresultPath = filename + "UResult.brep"; std::string VresultPath = filename + "VResult.brep";
            BRepTools::Write(UResult, UresultPath.c_str()); BRepTools::Write(VResult, VresultPath.c_str());

            SurfaceModelingTool tool;
            std::string knotsPath = filename + "knots.txt";
            tool.setKnotsOutputPath(knotsPath.c_str());
            // ����ļ��Ƿ���ڣ�������ڣ�����ļ�����
            std::ifstream checkFile(tool.getKnotsOuputPath());
            if (checkFile.is_open())
            {
                // �رռ���ļ���������
                checkFile.close();
                // ����ļ����ݣ�����д
                std::ofstream clearFile(tool.getKnotsOuputPath(), std::ios::trunc);
                clearFile.close();
            }

            tool.ContextToTxt("U:");
            for (auto debugKnots : uKnots)
                tool.KnotsToTxt(debugKnots);

            tool.ContextToTxt("------------------------------\nV:");
            for (auto debugKnots : vKnots)
                tool.KnotsToTxt(debugKnots);
        }
    }
}
