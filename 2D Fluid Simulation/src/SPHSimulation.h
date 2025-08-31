#pragma once
#include <random>
#include <iostream>
#include <array>
#include "primitives.h"

#define EPSILON 1e-7
#define PI 3.14159265358979323846

struct Particle
{
	Particle(Vec2 _pos = Vec2(0,0)) : pos(_pos.x, _pos.y), velocity(0.f, 0.f), force(0.f, 0.f), density(0), pressure(0.f), mass(1), forceExt(0,0) {}
	Vec2 pos, velocity, force, forceExt;
	float density, pressure, mass;
};

struct Cell {
	std::vector<int> particleIdxs;					// Indices of particles in this cell
	Vec2 position;									// Position of the cell in the grid
	std::array<int, 9> neighbourIdxs;				// Pointers to neighbouring cells
};



class SPHFluidSimulation2D {

public:
	SPHFluidSimulation2D(int _nParticles, Vec2 dimentions);
	void step(double dt);

	Particle* getState() { return particles; }
	Particle* getParticle(int index) { return &(particles[index]); }

	Cell* getCells() { return cells; }
	Cell* getCell(int index) { return &(cells[index]); }

	Vec2 getCellDimensions() { return cellsDimentions; }
	Vec2 getSimulationDimensions() { return simulationDimensions; }

	int nParticles = 0;
	int nCells;
	double cellSize;

private:
	Particle* particles;
	Cell* cells;

	Vec2 cellsDimentions;
	Vec2 simulationDimensions = Vec2(100, 100);

	double kernelRadius = 16.0f;
	double viscosity = 200.0f;


	// Simulation Constant	
	double RestDensity = 300.f;  // rest density
	double GasConstant = 2000.f; // const for equation of state
	float BoundDamping = -0.5f;


	// smoothing kernels defined in Müller and their gradients
	// adapted to 2D per "SPH Based Shallow Water Simulation" by Solenthaler et al.
	float POLY6 = 4.0f / (PI * pow(kernelRadius, 8.f));
	float SPIKY_GRAD = -10.f / (PI * pow(kernelRadius, 5.f));
	float VISC_LAP = 40.f / (PI * pow(kernelRadius, 5.f));



	void computeDensityPressure();
	void computeForces();

	int cell_index(int x, int y);

	std::array<int, 9> calculate_cell_neighbours(int x, int y);
	bool is_valid_cell(int x, int y);

	void updateCells();

};