#include "BezierCurve.h"
#include "Point3D.h"
#include "Matrix.h"

//-----------------------------------------------------------------------------

BezierCurve::BezierCurve(const std::vector<Point3D>& pointsFromUser)
{
	m_degree = (int)pointsFromUser.size() - 1;
	m_pointsGivenByUser = pointsFromUser;
	m_controlPoints.resize(pointsFromUser.size());

	double param = (double)1 / pointsFromUser.size();

	std::vector<double> coefficientVector((m_degree + 1) * (m_degree + 1), 0);

	std::vector<double> rightVectorForX(m_degree + 1, 0);
	std::vector<double> rightVectorForY(m_degree + 1, 0);
	std::vector<double> rightVectorForZ(m_degree + 1, 0);

	for (int iRow = 0; iRow <= m_degree; ++iRow)
	{
		double paramValue = 0;
		for (int iColumn = 0; iColumn <= m_degree; ++iColumn)
		{
			double basisFunction = (factorial(m_degree) / (factorial(m_degree - iColumn) * factorial(iColumn)))
				* pow(paramValue, iColumn) * pow((1 - paramValue), m_degree - iColumn);
			coefficientVector[(iRow * (m_degree + 1))+ iColumn] = basisFunction;
			paramValue += param;
		}
		rightVectorForX[iRow] = m_pointsGivenByUser[iRow].GetX();
		rightVectorForY[iRow] = m_pointsGivenByUser[iRow].GetY();
		rightVectorForZ[iRow] = m_pointsGivenByUser[iRow].GetZ();
	}

	Matrix coefficientMatrix(m_degree + 1, m_degree + 1, coefficientVector);
	Matrix inversCoefficientMatrix(coefficientMatrix.GetNumRows(), coefficientMatrix.GetNumColumns());
	coefficientMatrix.Inverse(inversCoefficientMatrix);

	Matrix righMatrixForX(m_degree + 1, 1, rightVectorForX);
	Matrix righMatrixForY(m_degree + 1, 1, rightVectorForY);
	Matrix righMatrixForZ(m_degree + 1, 1, rightVectorForZ);

	Matrix xCoordinatesOfControlPoints = inversCoefficientMatrix.GetMatrixMultiplication(righMatrixForX);
	Matrix yCoordinatesOfControlPoints = inversCoefficientMatrix.GetMatrixMultiplication(righMatrixForY);
	Matrix zCoordinatesOfControlPoints = inversCoefficientMatrix.GetMatrixMultiplication(righMatrixForY);

	for (int iSize =0; iSize < m_controlPoints.size(); ++iSize)
	{
		m_controlPoints[iSize].SetX(xCoordinatesOfControlPoints.GetMatrixElement(iSize,0));
		m_controlPoints[iSize].SetY(yCoordinatesOfControlPoints.GetMatrixElement(iSize, 0));
		m_controlPoints[iSize].SetZ(zCoordinatesOfControlPoints.GetMatrixElement(iSize, 0));
	}
}

//-----------------------------------------------------------------------------

int BezierCurve::factorial(int number) const
{
	if (number == 0 || number == 1)
		return 1;

	return number * factorial(number - 1);
}

//-----------------------------------------------------------------------------

void BezierCurve::GetPoinstAlongBazierCurve(std::vector<Point3D>& pointsAlongCurve) const
{
	double param = (double)1 / 1000;        //To get 1000 points along the curve
	for (double iParam = 0; iParam <= 1; iParam = iParam + param)
	{
		double xCoordinate = 0;
		double yCoordinate = 0;
		double zCoordinate = 0;

		for (int iDegree = 0; iDegree <= m_degree; ++iDegree)
		{
			double basisFunction = (factorial(m_degree) / (factorial(m_degree - iDegree) * factorial(iDegree)))
				* pow(iParam, iDegree) * pow((1 - iParam), m_degree - iDegree);

			xCoordinate += basisFunction * m_pointsGivenByUser[iDegree].GetX();
			yCoordinate += basisFunction * m_pointsGivenByUser[iDegree].GetY();
			zCoordinate += basisFunction * m_pointsGivenByUser[iDegree].GetZ();
		}
		pointsAlongCurve.emplace_back(xCoordinate, yCoordinate, zCoordinate);
	}
}

//-----------------------------------------------------------------------------

Point3D BezierCurve::GetApproxPointOfProjectionOnBezierCurve(const Point3D& pointToProject) const
{
	std::vector<Point3D> pointsAlongCurve;
	GetPoinstAlongBazierCurve(pointsAlongCurve);

	double distance = std::numeric_limits<double>::max();
	Point3D projectionPoint;

	for (int iPoint = 0; iPoint < pointsAlongCurve.size(); ++iPoint)
	{
		double tempDistance = pointsAlongCurve[iPoint].GetDistance(pointToProject);
		if (tempDistance < distance)
		{
			distance = tempDistance;
			projectionPoint = pointsAlongCurve[iPoint];
		}
	}

	return projectionPoint;
}

//-----------------------------------------------------------------------------


BezierCurve::~BezierCurve()
{

}

//-----------------------------------------------------------------------------
