#ifndef GORDENSURFACE_H
#define GORDENSURFACE_H

#pragma once

#include <vector>
#include <Geom_BSplineCurve.hxx>
#include <TopoDS_Face.hxx>

// 可能需要包含其他相关的头文件，如几何处理库等

class GordenSurface
{
public:
    // 根据给定的 u 和 v 方向曲线构建 Gordon 曲面
    static void BuildMyGordonSurf(
        std::vector<Handle(Geom_BSplineCurve)>& uCurves,
        std::vector<Handle(Geom_BSplineCurve)>& vCurves,
        TopoDS_Face& face);
};

#endif