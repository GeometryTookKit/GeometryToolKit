#pragma once
#include "Point2DInterface.h"

class Point2D : public Point2DInterface
{
private:
	double m_x;    //x co-ordinate of point
	double m_y;	   //y co-ordinate of point

public:

	//Default constructor
	Point2D();
	
	//Parameterized constructor
	Point2D(double x, double y);

	//Function to get x co-ordinate of the point 
	double GetX()const override;

	//Function to get y co-ordinate of the point
	double GetY()const override;

	//Function to get x and y co-ordinate of the point
	void  GetXY(double& x, double& y)const override;

	//Function to set x co-ordinate of the point 
	void SetX(double x) override;

	//Function to set y co-ordinate of the point 
	void SetY(double y) override;

	//Function to set x and y co-ordinate of the point 
	void SetXY(double x, double y) override;

	//This function calculate the distance between two points
	double GetDistance(Point2D secondPoint)const override; 

	//Destructor
	~Point2D() override;
};