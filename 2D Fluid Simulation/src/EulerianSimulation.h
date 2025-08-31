#pragma once
#include <random>
#include <iostream>
#include "primitives.h"

#define EPSILON 1e-7


// For Cell https://www.youtube.com/watch?v=iKAVRgIrUOU
//	Eulerian - grid
// The we try Lagrange (particle based)



// Assume incompressibility - (free gass + water)
struct Cell2D {
	Vec2 velocities[4];		// 0 = left, 1 = right, 2 = down, 3 = up
	double pressure = 0;
	double viscocity = 0;
	int s[4] = { 1,1,1,1 };	//each side of the cell - 1=fluid // 0 = left, 1 = right, 2 = down, 3 = up
	Vec2 forces;
};


class EulerianFluidSim2D {



public:

	EulerianFluidSim2D(Vec2 dim, double _gridSpacing, double _density, double _overrelaxationO = 1.9);

	void step(double dt);

	Cell2D* getState();
	Cell2D getCell(int x, int y);
	Vec2 getDimensions();


private:
	Cell2D* grid;
	Vec2 dimensions;
	int gridLen;
	double gridSpacing;
	double density;
	double overrelaxationO;


	int index(int x, int y);

	int compute_divergence(Cell2D cell);

	int compute_s(Cell2D cell);


	// 1.update velocity based on forces
	void applyForces(double dt);

	//2. Make fluid imcompressible
	//	- divergance, d = total velocity out
	void makeIncompressible(double dt);

	// Semi-Lagrangian advection for velocity field
	void applyAdvection(double dt);



	// Bilinear interpolation of velocity at arbitrary position
	Vec2 sampleVelocity(const Vec2& pos);



	// Average all four face velocities for a given cell
	Vec2 averageVelocitiesAtCell(int x, int y);
};


