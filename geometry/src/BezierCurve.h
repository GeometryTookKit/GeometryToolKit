#pragma once

//Geom header
#include "BezierCurveInterface.h"

//Forward declarations
class Point3D;

class BezierCurve : public BezierCurveInterface
{
private:
    int m_degree;                                  //Degree of Bazier curve
    std::vector<Point3D> m_pointsGivenByUser;     // Vector to store control points given by user
    std::vector<Point3D> m_controlPoints;         // Vector to store control points

    //Function to get factorial of the number
    int factorial(int number)const;

public:

    //Parameterized constructor
    BezierCurve(const std::vector<Point3D>& pointsFromUser);

    //Function to get points along bazier curve
    void GetPoinstAlongBazierCurve(std::vector<Point3D>& pointsAlongCurve)const override;

    //Function to get approx point of projection on Bezier curve
    Point3D GetApproxPointOfProjectionOnBezierCurve(const Point3D& pointToProject)const override;

    //Default destructor
    ~BezierCurve() ;
};