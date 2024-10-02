#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include input file and option
#include "setInput.h"

// include elastic Plate class
#include "elasticPlate.h"

// include time stepper
#include "timeStepper.h"

// include force
#include "inertialForce.h"
#include "dampingForce.h"
#include "externalPressureForce.h"
#include "hyperElasticMM.h"
#include "electricActuation.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	void setPlateStepper();

	void updateTimeStep();

	int simulationRunning();

	int numStretchingPair();
	int getNv();
	Vector2d getScaledCoordinate(int i, int j);
	Vector2d getScaledPos(int i);
		
private:

	// physical parameters
	bool render;
	bool saveData;
	double deltaTime;
	double totalTime;
	double YoungM;
	double density;
	double Possion;
	double stol;
	double forceTol;
	double scaleRendering;
	int maxIter;
	Vector3d gVector;
	double viscosity;
	double thickness;
	int nv;
	double inputPressure;
	double speed;
	double inputPhi;

	int Nstep;
	int timeStep;

	double characteristicForce;

	double currentTime;

	void plateBoundaryCondition();

	// Plate
	elasticPlate *plate;

	// stepper
	timeStepper *stepper;

	// force
	inertialForce *m_inertialForce;
	dampingForce *m_dampingForce;
	externalPressureForce *m_externalPressureForce;
	hyperElasticMM *m_hyperElasticMM;
	electricActuation *m_electricActuation;

	void updateEachStep();

	void computeReactionForce();
	VectorXd reactionForce;

	Vector2d xStart;
};

#endif
