// SETUP HERE:: https://www.youtube.com/watch?v=uO__ntYT-2Q
#include <iostream>
#include <random>
#include "UI.h"
#include "EulerianSimulation.h"
#include "SPHSimulation.h"

#define EPSILON 1e-7





class EulerianSimulation2D : public Window {
public:


	int itterationsPerFrame = 100;

	double zoomRate = 1E-3;
	double panRate = 2E-2;
	Vec3 initalScalefactor = Vec3(1E-2, 1E-2, 1);

	Vec2 gridDimensions = Vec2(100, 100);
	double dt = 0.001;
	EulerianFluidSim2D fluidSim = EulerianFluidSim2D(gridDimensions, 1, 1, 1);



	// Correctly calling the base class constructor using the initializer list
	EulerianSimulation2D(unsigned int height, unsigned int width)
		: Window(height, width) {
		wireframe = false;
	}



	bool OnUserCreate() override {

		SetScaleFactor(initalScalefactor);


		Cell2D* grid = fluidSim.getState(); // Assuming this returns pointer to your grid array
		Vec2 dim = fluidSim.getDimensions(); // Assuming you have a way to get dimensions

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(-1.0, 1.0);

		for (int j = 0; j < dim.y; j++) {
			for (int i = 0; i < dim.x; i++) {
				int idx = i + j * dim.x;

				grid[idx].forces = Vec2(0, -9.81);

				// Set side cells' s to 0 on left and right edges
				for (int j = 0; j < dim.y; j++) {
					int leftIdx = 0 + j * dim.x;        // Left edge cell index
					int rightIdx = (dim.x - 1) + j * dim.x; // Right edge cell index

					grid[leftIdx].s[0] = 0;  // Left side 0 for left cells
					//grid[leftIdx].s[1] = 0;  // Right side 0 for left cells
					//grid[leftIdx].s[2] = 0;  // Down side 0 for left cells
					//grid[leftIdx].s[3] = 0;  // Up side 0 for left cells

					//grid[rightIdx].s[0] = 0; // Left side 0 for right cells
					grid[rightIdx].s[1] = 0; // Right side 0 for right cells
					//grid[rightIdx].s[2] = 0; // Down side 0 for right cells
					//grid[rightIdx].s[3] = 0; // Up side 0 for right cells
				}

				// Set side cells' s to 0 on top and bottom edges
				for (int i = 0; i < dim.x; i++) {
					int bottomIdx = i + 0 * dim.x;       // Bottom edge cell index
					int topIdx = i + (dim.y - 1) * dim.x; // Top edge cell index

					//grid[bottomIdx].s[0] = 0;  // Left side 0 for bottom cells
					//grid[bottomIdx].s[1] = 0;  // Right side 0 for bottom cells
					grid[bottomIdx].s[2] = 0;  // Down side 0 for bottom cells
					//grid[bottomIdx].s[3] = 0;  // Up side 0 for bottom cells

					//grid[topIdx].s[0] = 0;     // Left side 0 for top cells
					//grid[topIdx].s[1] = 0;     // Right side 0 for top cells
					//grid[topIdx].s[2] = 0;     // Down side 0 for top cells
					grid[topIdx].s[3] = 0;     // Up side 0 for top cells
				}


				// Assign random initial velocities for all faces
				grid[idx].velocities[0] = Vec2(dis(gen), dis(gen));  // left face velocity random
				grid[idx].velocities[1] = Vec2(dis(gen), dis(gen));  // right face velocity random
				grid[idx].velocities[2] = Vec2(dis(gen), dis(gen));  // up face velocity random
				grid[idx].velocities[3] = Vec2(dis(gen), dis(gen));  // down face velocity random
				
				//// Create a pressure gradient: low pressure on left, high on right
				//grid[idx].pressure = (double)i / dim.x * 10.0;
			}
		}




		return true;
	}



	bool OnUserUpdate() override {
		

		// Control Logic
		if (GetMousePosition().x >= 0 && GetMousePosition().y >= 0) {
			if (GetMousePosition().x <= 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(-panRate, 0, 0));
			}
			else if (GetMousePosition().x >= GetScreenSize().x - 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(panRate, 0, 0));
			}
			else if (GetMousePosition().y <= 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(0, panRate, 0));
			}
			else if (GetMousePosition().y >= GetScreenSize().y - 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(0, -panRate, 0));
			}
		}



		if (IsScrollUp()) {
			SetScaleFactor(GetScaleFactor() + Vec3(zoomRate, zoomRate, 1));
		}

		if (IsScrollDown()) {
			SetScaleFactor(GetScaleFactor() - Vec3(zoomRate, zoomRate, 1));
		}




		// Graphics Logic


		//step simulation
		for (int i = 0; i < itterationsPerFrame; i++)
		{
			fluidSim.step(dt);
		}

		// Retrieve grid state
		Cell2D* state = fluidSim.getState();


		// Find pressure range min and max for normalization first
		double minPressure = state[0].pressure;
		double maxPressure = state[0].pressure;

		for (int i = 0; i < gridDimensions.x; i++) {
			for (int j = 0; j < gridDimensions.y; j++) {
				double p = abs(fluidSim.getCell(i, j).pressure);
				if (p < minPressure) minPressure = p;
				if (p > maxPressure) maxPressure = p;
			}
		}

		// Draw grid with pressure-based colors
		for (int i = 0; i < gridDimensions.x; i++) {
			for (int j = 0; j < gridDimensions.y; j++) {
				Quad cellQuad = Quad(
					Vec3(i, j, 0),
					Vec3(i + 1, j, 0),
					Vec3(i + 1, j + 1, 0),
					Vec3(i, j + 1, 0)
				);

				double p = fluidSim.getCell(i, j).pressure;

				float normalizedP = (float)((p - minPressure) / (maxPressure - minPressure + EPSILON));


				// Map pressure normalized [0,1] to color (e.g., blue to red)
				Vec3 cellColor = Vec3(normalizedP, 0.0f, 1.0f - normalizedP);
				DrawQuad(cellQuad, cellColor);
			}
		}



		return true;

	}
};






























class SPHSimulation2D : public Window {
public:


	int itterationsPerFrame = 1;

	double zoomRate = 1E-3;
	double panRate = 2E-2;
	Vec3 initalScalefactor = Vec3(1E-2, 1E-2, 1);

	int nParticles = 5000;
	Vec2 simulationDimentions = Vec2(256, 128);
	double dt = 0.001;

	SPHFluidSimulation2D fluidSim = SPHFluidSimulation2D(nParticles, simulationDimentions);

	FPSCounter fpsCounter;

	SPHSimulation2D(unsigned int height, unsigned int width)
		: Window(height, width) {

		wireframe = false;
		fpsCounter = FPSCounter();
	}



	bool OnUserCreate() override {

		SetScaleFactor(initalScalefactor);
		
		// Pressure gradient along y-axis
		float particleDiameter = 1.5;
		int perRow = floor(simulationDimentions.x / particleDiameter);
		int numRows = ceil(float(nParticles) / perRow);

		float minSpacing = particleDiameter * 0.5f; // Most dense region
		float maxSpacing = particleDiameter * 2.0f; // Least dense region

		for (int i = 0; i < nParticles; i++)
		{
			int row = i / perRow;
			int col = i % perRow;

			// Interpolate spacing for this row to create gradient (0: bottom, 1: top)
			float t = float(row) / float(numRows - 1);
			float rowSpacing = minSpacing * (1 - t) + maxSpacing * t;

			fluidSim.getParticle(i)->pos = Vec2(
				particleDiameter / 2 + col * particleDiameter,
				rowSpacing / 2 + row * rowSpacing
			);
			fluidSim.getParticle(i)->forceExt = Vec2(0, -9.81); // gravity
		}
	
		return true;
	}



	bool OnUserUpdate() override {
		
		fpsCounter.step();

		SetTitle("2D SPH Fluid Simulation	[FPS:" + std::to_string(fpsCounter.getFPS()) + "]");

		// Control Logic
		if (GetMousePosition().x >= 0 && GetMousePosition().y >= 0) {
			if (GetMousePosition().x <= 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(-panRate, 0, 0));
			}
			else if (GetMousePosition().x >= GetScreenSize().x - 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(panRate, 0, 0));
			}
			else if (GetMousePosition().y <= 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(0, panRate, 0));
			}
			else if (GetMousePosition().y >= GetScreenSize().y - 25) {
				SetScreenOrigin(GetScreenOrigin() + Vec3(0, -panRate, 0));
			}
		}



		if (IsScrollUp()) {
			SetScaleFactor(GetScaleFactor() + Vec3(zoomRate, zoomRate, 1));
		}

		if (IsScrollDown()) {
			SetScaleFactor(GetScaleFactor() - Vec3(zoomRate, zoomRate, 1));
		}



		// ------------------- Update Solver State -------------------

		//step simulation
		for (int i = 0; i < itterationsPerFrame; i++)
		{
			fluidSim.step(dt);
		}

		// Retrieve state
		Particle* state = fluidSim.getState();


		// ==============================================================
		//						Graphics Logic
		// ==============================================================




		Vec3 boundaryColor = Vec3(0.f, 0.f, 0.f);
		Vec3 cellBoundaryColor = Vec3(0.8f, 0.2f, 0.2f);


		// ------------------ Draw Cells Boundary ------------------ 

		for (int i = 0; i < fluidSim.nCells; i++)
		{
			DrawWireframeQuad(Quad(
				Vec3(fluidSim.getCell(i)->position.x, fluidSim.getCell(i)->position.y, 0),
				Vec3(fluidSim.getCell(i)->position.x + fluidSim.cellSize, fluidSim.getCell(i)->position.y, 0),
				Vec3(fluidSim.getCell(i)->position.x + fluidSim.cellSize, fluidSim.getCell(i)->position.y + fluidSim.cellSize, 0),
				Vec3(fluidSim.getCell(i)->position.x, fluidSim.getCell(i)->position.y + fluidSim.cellSize, 0)
			), cellBoundaryColor, 0.375f);
		}



		// ------------------ Draw Simulation Boundary ------------------ 

		Vec2 simMin = Vec2(0, 0);
		Vec2 simMax = simulationDimentions;

		DrawWireframeQuad(Quad(
			Vec3(simMin.x, simMin.y, 0),
			Vec3(simMax.x, simMin.y, 0),
			Vec3(simMax.x, simMax.y, 0),
			Vec3(simMin.x, simMax.y, 0)
		), boundaryColor, 0.5f);



		// ------------------ Draw Particles ------------------

		// normalise pressure for color mapping
		float minPressure = FLT_MAX;
		float maxPressure = -FLT_MAX;

		// Find minimum and maximum pressure first
//#pragma omp parallel for
		for (int i = 0; i < fluidSim.nParticles; ++i) {
			float p = state[i].pressure;
			if (p < minPressure) minPressure = p;
			if (p > maxPressure) maxPressure = p;
		}



		// Plot particles with color based on pressure
		for (int i = 0; i < fluidSim.nParticles; ++i) {
			Particle particle = state[i];

			Circle particleCirc = Circle(1.0f, Vec3(particle.pos.x, particle.pos.y, 0), 6);

			float normPressure = (particle.pressure - minPressure) / (maxPressure - minPressure);
			Vec3 particleColor = Vec3(normPressure, 0.0f, 1.0f - normPressure);			// Interpolate color: Blue (low) to Red (high)

			DrawCircle(particleCirc, particleColor);
		}



		return true;
	}
};










int main() {

	//EulerianSimulation2D sim = EulerianSimulation2D(800,1000);
	//sim.Init();

	SPHSimulation2D sim = SPHSimulation2D(800,1000);
	sim.Init();



}