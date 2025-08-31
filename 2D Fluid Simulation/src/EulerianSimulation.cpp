#include <random>
#include <iostream>
#include "EulerianSimulation.h"

#define EPSILON 1e-7


EulerianFluidSim2D::EulerianFluidSim2D(Vec2 dim, double _gridSpacing, double _density, double _overrelaxationO) {
	dimensions = dim;
	grid = new Cell2D[dim.x * dim.y];
	gridLen = dim.x * dim.y;
	gridSpacing = _gridSpacing;
	density = _density;
	overrelaxationO = _overrelaxationO;
}


void EulerianFluidSim2D::step(double dt) {
	applyForces(dt);
	makeIncompressible(dt);
	applyAdvection(dt);
}

Cell2D* EulerianFluidSim2D::getState() {
	return grid;
}

Cell2D EulerianFluidSim2D::getCell(int x, int y) {
	return grid[index(x, y)];
}

Vec2 EulerianFluidSim2D::getDimensions() {
	return dimensions;
}



int EulerianFluidSim2D::index(int x, int y) {
	return x * dimensions.y + y;
}

int EulerianFluidSim2D::compute_divergence(Cell2D cell) {
	return overrelaxationO * (cell.velocities[0].x - cell.velocities[1].x) + (cell.velocities[2].y - cell.velocities[3].y);
}

int EulerianFluidSim2D::compute_s(Cell2D cell) {
	return cell.s[0] + cell.s[1] + cell.s[2] + cell.s[3];
}


// 1.update velocity based on forces
void EulerianFluidSim2D::applyForces(double dt) {
	for (int i = 0; i < gridLen; i++)
	{
		for (int n = 0; n < 4; n++)
		{
			grid[i].velocities[n] += grid[i].forces * dt;
		}
	}

}

//2. Make fluid imcompressible
//	- divergance, d = total velocity out
void EulerianFluidSim2D::makeIncompressible(double dt) {

	for (int x = 0; x < dimensions.x; x++)
	{
		for (int y = 0; y < dimensions.y; y++)
		{
			Cell2D* cell = &grid[index(x, y)]; // err

			double d = compute_divergence(*cell);
			double s = compute_s(*cell);

			// force incompressibility
			cell->velocities[0].x += d * (cell->s[0] / s); // left
			cell->velocities[1].x -= d * (cell->s[1] / s); // right
			cell->velocities[2].y += d * (cell->s[2] / s); // down
			cell->velocities[3].y -= d * (cell->s[3] / s); // up

			// update pressure (not necessary for simulation)
			cell->pressure = d / s * (density * gridSpacing / dt);
		}
	}
}

// Semi-Lagrangian advection for velocity field
void EulerianFluidSim2D::applyAdvection(double dt) {
	Cell2D* newGrid = new Cell2D[gridLen];
	memcpy(newGrid, grid, sizeof(Cell2D) * gridLen); // copy current states

	for (int j = 1; j < dimensions.y - 1; j++) {
		for (int i = 1; i < dimensions.x - 1; i++) {
			int idx = index(i, j);

			// Skip solid cells (assuming s[] per side indicates solid boundaries)
			// Here, a cell is fluid if all sides are non-solid (0)
			bool isSolid = false;
			for (int side = 0; side < 4; side++) {
				if (grid[idx].s[side] == 1) {
					isSolid = true;
					break;
				}
			}
			if (isSolid)
				continue;



			for (int face = 0; face < 4; face++) {
				// Define sample position for velocity on cell face
				Vec2 pos;
				switch (face) {
				case 0: pos = Vec2(i, j + 0.5); break;     // left face center
				case 1: pos = Vec2(i + 1, j + 0.5); break; // right face center
				case 2: pos = Vec2(i + 0.5, j); break;     // down face center
				case 3: pos = Vec2(i + 0.5, j + 1); break; // up face center
				}

				// Get velocity at pos
				Vec2 vel = sampleVelocity(pos);

				// Backtrace position
				Vec2 prevPos = pos - vel * dt;

				// Clamp to domain boundaries
				prevPos.x = std::max(0.5, std::min(prevPos.x, dimensions.x - 1.5));
				prevPos.y = std::max(0.5, std::min(prevPos.y, dimensions.y - 1.5));

				// Sample velocity at previous position
				Vec2 newVel = sampleVelocity(prevPos);

				// Store advected velocity in new grid
				newGrid[idx].velocities[face] = newVel;
			}
		}
	}

	// Copy new velocities back to main grid
	for (int i = 0; i < gridLen; i++) {
		for (int face = 0; face < 4; face++) {
			grid[i].velocities[face] = newGrid[i].velocities[face];
		}
	}
	delete[] newGrid;
}



// Bilinear interpolation of velocity at arbitrary position
Vec2 EulerianFluidSim2D::sampleVelocity(const Vec2& pos) {
	int x0 = (int)floor(pos.x);
	int y0 = (int)floor(pos.y);
	int x1 = x0 + 1;
	int y1 = y0 + 1;

	double sx = pos.x - x0;
	double sy = pos.y - y0;

	Vec2 v00 = averageVelocitiesAtCell(x0, y0);
	Vec2 v10 = averageVelocitiesAtCell(x1, y0);
	Vec2 v01 = averageVelocitiesAtCell(x0, y1);
	Vec2 v11 = averageVelocitiesAtCell(x1, y1);

	return v00 * (1 - sx) * (1 - sy) +
		v10 * sx * (1 - sy) +
		v01 * (1 - sx) * sy +
		v11 * sx * sy;
}



// Average all four face velocities for a given cell
Vec2 EulerianFluidSim2D::averageVelocitiesAtCell(int x, int y) {
	Vec2 sum(0, 0);

	for (int i = 0; i < 4; i++) {
		sum += grid[index(x, y)].velocities[i];
	}
	return sum * 0.25;
}