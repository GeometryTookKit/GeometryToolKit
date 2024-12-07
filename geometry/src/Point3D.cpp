#include "Point3D.h"
#include "Vector3D.h"
#include <cmath>

//-----------------------------------------------------------------------------

Point3D::Point3D()
{
	m_x = 0.0;
	m_y = 0.0;
	m_z = 0.0;
}

//-----------------------------------------------------------------------------

Point3D::Point3D(double x, double y, double z)
{
	m_x = x;
	m_y = y;
	m_z = z;
}

//-----------------------------------------------------------------------------

double Point3D::GetX()const
{
	return m_x;
}

//-----------------------------------------------------------------------------

double Point3D::GetY()const
{
	return m_y;
}

//-----------------------------------------------------------------------------

double Point3D::GetZ()const
{
	return m_z;
}

//-----------------------------------------------------------------------------

void Point3D::GetXYZ(double& x, double& y, double& z)const
{
	x = m_x;
	y = m_y;
	z = m_z;
}

//-----------------------------------------------------------------------------

void Point3D::SetX(double x)
{
	m_x = x;
}

//-----------------------------------------------------------------------------

void Point3D::SetY(double y)
{
	m_y = y;
}

//-----------------------------------------------------------------------------

void Point3D::SetZ(double z)
{
	m_z = z;
}

//-----------------------------------------------------------------------------

void Point3D::SetXYZ(double x, double y, double z)
{
	m_x = x;
	m_y = y;
	m_z = z;
}

//-----------------------------------------------------------------------------

double Point3D::GetDistance(Point3D secondPoint)const
{
	return std::sqrt(std::pow((m_x - secondPoint.m_x), 2)
					+ std::pow((m_y - secondPoint.m_y), 2)
	                + std::pow((m_z - secondPoint.m_z), 2));
}

//-----------------------------------------------------------------------------

Point3D Point3D::GetPointAlongVector(const Vector3D& vector, double distance)
{
	return Point3D((m_x + distance*vector.GetI()),(m_y + distance * vector.GetJ()),(m_z + distance * vector.GetK()));
}

//-----------------------------------------------------------------------------

Point3D::~Point3D()
{
}

//-----------------------------------------------------------------------------
