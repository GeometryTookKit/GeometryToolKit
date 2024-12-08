#pragma once

//Geom header
#include "BezierCurveInterface.h"

//Forward declarations
class Point3D;

class BezierCurve : public BezierCurveInterface
{
private:
    int m_degree;                                  //Degree of Bazier curve
    std::vector<Point3D> m_controlPoints;         // Vector to store control points

    //Function to get factorial of the number
    int factorial(int number)const;

    //De Casteljau's Algorithm
    Point3D deCasteljauAlgorithm(const double param, const int i, const int j)const;

public
    :
    //Default constructor
    BezierCurve() = default;

    //Parameterized constructor
    BezierCurve(const std::vector<Point3D>& pointsFromUser);

    //Function to get points along bazier curve
    void GetPoinstAlongBazierCurve(const int numPoints, std::vector<Point3D>& pointsAlongCurve)const override;

    //Function to get approx point of projection on Bezier curve
    Point3D GetApproxPointOfProjectionOnBezierCurve(const Point3D& pointToProject)const override;

    //Function to insert a control point at given index
    void AddControlPoint(const size_t index, const Point3D& controlPoint);

    //Function to remove control point
    void RemoveControlPoint(const size_t index);

    //Function to update the control point at given index
    void UpdadeControlPoint(const size_t index, const Point3D& controlPoint);

    //Function to get control point at given index
    void GetControlPoint(const size_t index, Point3D& controlPoint)const;

    //Function to get all control points of Bezier curve
    void GetControlPoints(std::vector<Point3D>& controlPoints)const;

    //Function to set control points of the curve
    void SetControlPoints(const std::vector<Point3D>& controlPoints);

    //Function to get degree of the curve
    int GetDegree()const;

    //Function to get number of control points
    size_t GetNumControlPoints()const;

    //Function to get point along curve for given parametric value
    void GetPointOnCurve(const double param, Point3D& pointOnCurve)const;

    //Function to get point along curve using De Casteljau's algorithm
    void GetPointOnCurveUsingDeCasteljau(const double param, Point3D& point)const;

    //Default destructor
    ~BezierCurve() ;
};