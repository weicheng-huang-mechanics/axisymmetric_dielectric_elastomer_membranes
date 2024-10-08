#include "electricActuation.h"
#include <iostream>

electricActuation::electricActuation(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
  stepper = &m_stepper;

	Jss.setZero(4, 4);
	flocal = VectorXd::Zero(4);

  thickness = plate->thickness;

	localDOF = VectorXi::Zero(4);

  epsilon = 1.0;

  phi = 0.0;
}

electricActuation::~electricActuation()
{
	;
}

void electricActuation::computeFs()
{
	totalForce = VectorXd::Zero(plate->ndof);

  totalEnergyS = 0.0;

	for (int k = 1; k < plate->edgeNum - 1; k++)
	{
		flocal = VectorXd::Zero(4);

		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		p = plate->getVertex(ind1);
		p1 = plate->getVertex(ind2);

		xk = p[0];
		yk = p[1];
		xkp1 = p1[0];
		ykp1 = p1[1];

		l_k = plate->v_edgeElement[k].refLength;

		rBar1 = plate->v_edgeElement[k].x_1_start(0);
		rBar2 = plate->v_edgeElement[k].x_2_start(0);
    r_k = (rBar1 + rBar2) / 2;

    //cout << (xk+xkp1) / (2 * r_k) << " " << (p1 - p).norm() / l_k << endl;

		flocal = computeStretchingForce(xk, yk, xkp1, ykp1, l_k, r_k);

		flocal = (- 2 * l_k * M_PI * r_k * thickness ) * flocal;

		localDOF = VectorXi::Zero(4);

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			stepper->addForce(localDOF(i), - flocal(i));

			totalForce(localDOF(i)) = totalForce(localDOF(i)) + flocal(i);
		}
	}
}

void electricActuation::computeJs()
{
	for (int k = 1; k < plate->edgeNum - 1; k++)
	{
		Jss.setZero(4, 4);

		ind1 = plate->v_edgeElement[k].nv_1;
    ind2 = plate->v_edgeElement[k].nv_2;

    p = plate->getVertex(ind1);
    p1 = plate->getVertex(ind2);

    xk = p[0];
    yk = p[1];
    xkp1 = p1[0];
    ykp1 = p1[1];

    l_k = plate->v_edgeElement[k].refLength;

    rBar1 = plate->v_edgeElement[k].x_1_start(0);
    rBar2 = plate->v_edgeElement[k].x_2_start(0);
    r_k = (rBar1 + rBar2) / 2;

		Jss = computeStretchingJacobian(xk, yk, xkp1, ykp1, l_k, r_k);
		Jss = ( 2 * l_k * M_PI * r_k * thickness ) * Jss; // scale with stiffness

		localDOF = VectorXi::Zero(4);

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				stepper->addJacobian(localDOF(i), localDOF(j), Jss(i,j));
			}
		}
	}
}

void electricActuation::setFirstJacobian()
{
	for (int k = 0; k < plate->edgeNum; k++)
	{
		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		localDOF = VectorXi::Zero(4);

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				stepper->addJacobian(localDOF(i), localDOF(j), 1);
			}
		}
	}
}


VectorXd electricActuation::computeStretchingForce(double xa, double ya, double xb, double yb, double lBar, double rBar)
{
  VectorXd vecResult;

  vecResult = ListVec((0.25*epsilon*pow(phi,2)*(-xa + xb)*pow(xa + xb,2))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
    (0.25*epsilon*pow(phi,2)*(xa + xb)*(pow(-xa + xb,2) + pow(-ya + yb,2)))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
   (0.25*epsilon*pow(phi,2)*pow(xa + xb,2)*(-ya + yb))/
    (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
   (-0.25*epsilon*pow(phi,2)*(-xa + xb)*pow(xa + xb,2))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
    (0.25*epsilon*pow(phi,2)*(xa + xb)*(pow(-xa + xb,2) + pow(-ya + yb,2)))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
   (-0.25*epsilon*pow(phi,2)*pow(xa + xb,2)*(-ya + yb))/
    (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)));

  return vecResult;
}


MatrixXd electricActuation::computeStretchingJacobian(double xa, double ya, double xb, double yb, double lBar, double rBar)
{
  MatrixXd matResult;

  matResult = ListMat(ListVec((1.*epsilon*pow(phi,2)*(-xa + xb)*(xa + xb))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
     (0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
     (0.25*epsilon*pow(phi,2)*(pow(-xa + xb,2) + pow(-ya + yb,2)))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    0. + (0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
     (0.25*epsilon*pow(phi,2)*(pow(-xa + xb,2) + pow(-ya + yb,2)))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (-0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2))),
   ListVec((0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (-0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2))),
   ListVec(0. + (0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
     (0.25*epsilon*pow(phi,2)*(pow(-xa + xb,2) + pow(-ya + yb,2)))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (-1.*epsilon*pow(phi,2)*(-xa + xb)*(xa + xb))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
     (0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)) - 
     (0.25*epsilon*pow(phi,2)*(pow(-xa + xb,2) + pow(-ya + yb,2)))/
      (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (-0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2))),
   ListVec((-0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (-0.5*epsilon*pow(phi,2)*(xa + xb)*(-ya + yb))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2)),
    (-0.25*epsilon*pow(phi,2)*pow(xa + xb,2))/
     (pow(lBar,2)*pow(rBar,2)*pow(thickness,2))));

  return matResult;
}

VectorXd electricActuation::ListVec(double a1, double a2, double a3, double a4)
{
  VectorXd vecResult;

  vecResult.setZero(4, 1);

  vecResult(0) = a1;
  vecResult(1) = a2;
  vecResult(2) = a3;
  vecResult(3) = a4;

  return vecResult;
}

MatrixXd electricActuation::ListMat(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4)
{
  MatrixXd matResult;

  matResult.setZero(4, 4);

  matResult.col(0) = a1;
  matResult.col(1) = a2;
  matResult.col(2) = a3;
  matResult.col(3) = a4;

  return matResult;
}