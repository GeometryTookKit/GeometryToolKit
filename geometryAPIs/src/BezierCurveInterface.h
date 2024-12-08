#pragma once

#include <vector>

//Forward declarations
class Point3D;

class BezierCurveInterface
{

public:

    //Function to get points along bazier curve
    virtual void GetPoinstAlongBazierCurve(const int numPoints, std::vector<Point3D>& pointsAlongCurve)const = 0;

    //Function to get approx point of projection on Bezier curve
    virtual Point3D GetApproxPointOfProjectionOnBezierCurve(const Point3D& pointToProject)const = 0;

    //Default destructor
    virtual ~BezierCurveInterface() = default;
};