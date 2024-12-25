/*
*    Copyright (c) 2018 Shing Liu All Rights Reserved.
*
*           File : occQt.h
*         Author : Shing Liu(eryar@163.com)
*           Date : 2018-01-08 20:00
*        Version : OpenCASCADE7.2.0 & Qt5.7.1
*
*    Description : OpenCASCADE in Qt.
*/

#ifndef OCCQT_H
#define OCCQT_H

#include "ui_occQt.h"

#include <AIS_InteractiveContext.hxx>
#include "TopoDS_Shape.hxx"
#include "TopoDS_Face.hxx"\
#include <vector>
#include <D:\OpenCASCADE-7.7.0-vc14-64\opencascade-7.7.0\testPoly\RTree.h>
#include "TopoDS_Compound.hxx"  
#include "TopoDS_Edge.hxx"
#include "Geom_BoundedCurve.hxx"
#include "AIS_Shape.hxx"
#include "../../OpenCASCADE-7.7.0-vc14-64/qt5.11.2-vc14-64/include/QtWidgets/qtextbrowser.h"
#include <mutex>
class OccView;

//! Qt main window which include OpenCASCADE for its central widget.
class occQt : public QMainWindow
{
    Q_OBJECT

public:
    //! constructor/destructor.
    occQt(QWidget *parent = nullptr);
    ~occQt();
    template <typename T>
    void Visualize(const T& object, const Quantity_Color& color = Quantity_NOC_BISQUE);
    template <typename T>
    void Visualize(const std::vector<T>& object, const Quantity_Color& color = Quantity_NOC_BISQUE);
    void SetVisualize(bool flag) 
    {
        isVisualize = flag;
    }
protected:
    //! create all the actions.
    void createActions(void);

    //! create all the menus.
    void createMenus(void);

    //! create the toolbar.
    void createToolBars(void);

private slots:
    //! show about box.
    void about(void);

    //! make box test.
    void importFile(void);

    //! make cone test.
    void Triangulation(void);

    //! make sphere test.
    void TriangulationIntersection(void);

    //! make cylinder test.
    void ClearDisplay(void);

    //! make torus test.
    void PrintInfo(void);

    //! fillet test.
    void RandomExport(void);

    //! chamfer test.
    void MakePoint(void);

    //! test extrude algorithm.
    void MakeCurve(void);

    //! test revol algorithm.
    void MakeSurface(void);

    //! test loft algorithm.
    void ExportFile(void);

    void GenerateIsoCurves(void);

private:
    Ui::occQtClass ui;
    Handle(Geom_BoundedCurve) boundedCurve;
    // wrapped the widget for occ.
    OccView* myOccView;
    TopoDS_Shape shape1, shape2;
    std::vector<TopoDS_Face> faceList[2];
    std::vector<TriangleVertex> triPoints[2];
    TopoDS_Compound result;
    TopoDS_Shape curveShape, pointShape;
    std::mutex coneMutex;  // 互斥锁，用于同步 Triangulation 和 TriangulationIntersection
    std::condition_variable coneCondition;  // 条件变量，用于通知 TriangulationIntersection 等待结束
    bool coneProcessingComplete{ false };  // 标志 Triangulation 是否完成
    QTextEdit* outputTextEdit;
    int changeCount;
    bool intersectionDone;

    bool isVisualize;
};

#endif // OCCQT_H
