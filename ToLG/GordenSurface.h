#ifndef GORDENSURFACE_H
#define GORDENSURFACE_H

#pragma once

#include <vector>
#include <Geom_BSplineCurve.hxx>
#include <TopoDS_Face.hxx>

// ������Ҫ����������ص�ͷ�ļ����缸�δ�����

class GordenSurface
{
public:
    // ���ݸ����� u �� v �������߹��� Gordon ����
    static void BuildMyGordonSurf(
        std::vector<Handle(Geom_BSplineCurve)>& uCurves,
        std::vector<Handle(Geom_BSplineCurve)>& vCurves,
        TopoDS_Face& face);
};

#endif