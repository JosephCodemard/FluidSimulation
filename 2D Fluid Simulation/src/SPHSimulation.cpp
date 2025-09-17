#include "SPHSimulation.h"
#include <algorithm>

// https://lucasschuermann.com/writing/implementing-sph-in-2d


SPHFluidSimulation2D::SPHFluidSimulation2D(int _nParticles, Vec2 dimentions) {
	nParticles = _nParticles;
    simulationDimensions = dimentions;

	particles = new Particle[nParticles];


    // ----------------- generate cells ----------------- 
	cellSize = kernelRadius;

    cellsDimentions = Vec2(
        std::ceil(double(simulationDimensions.x) / kernelRadius),
        std::ceil(double(simulationDimensions.y) / kernelRadius)
    );
	nCells = cellsDimentions.x * cellsDimentions.y;
    cells = new Cell[nCells];

    // Initialize cells
    for (int x = 0; x < cellsDimentions.x; ++x) {
        for (int y = 0; y < cellsDimentions.y; ++y) {
            cells[cell_index(x, y)].position = x * kernelRadius;
            cells[cell_index(x, y)].position.y = y * kernelRadius;
        }
    }

	// calculate neighbour indices for each cell
    for (int x = 0; x < cellsDimentions.x; ++x) {
        for (int y = 0; y < cellsDimentions.y; ++y) {
            cells[cell_index(x, y)].neighbourIdxs = calculate_cell_neighbours(x, y);
        }
    }
}


bool SPHFluidSimulation2D::is_valid_cell(int x, int y) {
    return (x >= 0 && x < cellsDimentions.x && y >= 0 && y < cellsDimentions.y);
}


int SPHFluidSimulation2D::cell_index(int x, int y) {
    return x + y * cellsDimentions.x;
}



std::array<int,9> SPHFluidSimulation2D::calculate_cell_neighbours(int x, int y) {
    std::array<int, 9> neighbourIdxs;
    int i = 0;

    for (int nx = -1; nx <= 1; ++nx) {
        for (int ny = -1; ny <= 1; ++ny) {

            if (nx == 0 && ny == 0) continue;
            
            int neighborX = x + nx;
            int neighborY = y + ny;
            
            if (is_valid_cell(neighborX, neighborY)) {
                neighbourIdxs[i] = cell_index(neighborX, neighborY);
            }
            else {
                neighbourIdxs[i] = -1;
            }
            i++;
        }
    }
    neighbourIdxs[8] = cell_index(x, y); // self at end
    return neighbourIdxs;
}



void SPHFluidSimulation2D::updateCells() {
    for (int i = 0; i < nCells; i++)
        cells[i].particleIdxs.clear(); // clear for new step

    for (int i = 0; i < nParticles; i++) {
        Particle& particle = particles[i];
        
        int cellX = std::max(0, std::min(int(particle.pos.x / cellSize), (int)(cellsDimentions.x - 1)));
        int cellY = std::max(0, std::min(int(particle.pos.y / cellSize), (int)(cellsDimentions.y - 1)));

        cells[cell_index(cellX, cellY)].particleIdxs.push_back(i);
    }
}











void SPHFluidSimulation2D::step(double dt) {

    updateCells();
	computeDensityPressure();
    computeForces();



    for (int i = 0; i < nParticles; i++)
    {
        Particle* particle = &particles[i];

        particle->velocity += particle->force * (dt / particle->density);
        particle->pos += particle->velocity * dt;



        // enforce boundary conditions
        if (particle->pos.x < EPSILON)
        {
            particle->velocity.x = -BoundDamping * particle->velocity.x;
            particle->pos.x = EPSILON;
        }
        if (particle->pos.x > simulationDimensions.x)
        {
            particle->velocity.x = -BoundDamping * particle->velocity.x;
            particle->pos.x = simulationDimensions.x - EPSILON;
        }

        if (particle->pos.y < EPSILON)
        {
            particle->velocity.y = -BoundDamping * particle->velocity.y;
            particle->pos.y = EPSILON;
        }
        if (particle->pos.y > simulationDimensions.y)
        {
            particle->velocity.y = -BoundDamping * particle->velocity.y;
            particle->pos.y = simulationDimensions.y - EPSILON;
        }
    }
}



//void SPHFluidSimulation2D::computeDensityPressure() {
//
//#pragma omp parallel for
//    for (int i = 0; i < nParticles; i++)
//    {
//        Particle* particleA = &particles[i];
//
//        particleA->density = 0.0f;
//
//        for (int j = 0; j < nParticles; j++)
//        {
//            Particle particleB = particles[j];
//
//            double d = (particleA->pos - particleB.pos).magnitude();
//
//            if (d < kernelRadius)
//            {
//                particleA->density += particleB.mass * POLY6 * pow(smoothingRadius - d, 3.0f);
//            }
//        }
//
//        particleA->pressure = GasConstant * (particleA->density - RestDensity);
//    }
//}
//


//
//void SPHFluidSimulation2D::computeForces() {
//
//    #pragma omp parallel for
//    for (int i = 0; i < nParticles; i++)
//    {
//        Particle* particleA = &particles[i];
//        
//        Vec2 fpress(0.f, 0.f);
//        Vec2 fvisc(0.f, 0.f);
//
//        for (int j = 0; j < nParticles; j++)
//        {
//            Particle* particleB = &particles[j];
//            
//            if (i == j)   continue;
//
//            Vec2 d = (particleA->pos - particleB->pos);
//			double r = d.magnitude();
//
//
//            if (r < kernelRadius)
//            {
//				Vec2 deltaVel = particleB->velocity - particleA->velocity;
//
//                fpress += d.unit() * (particleB->mass * (particleA->pressure + particleB->pressure) / (2.f * particleB->density) * SPIKY_GRAD * pow(kernelRadius - r, 3.f)); // * -1 (removed??)
//                fvisc += deltaVel * ((viscosity * particleB->mass) / particleB->density * VISC_LAP * (kernelRadius - r));
//            }
//        }
//    
//        particleA->force = fpress + fvisc + particleA->forceExt;
//
//    }
//}







void SPHFluidSimulation2D::computeForces() {

#pragma omp parallel for
    for (int c = 0; c < nCells; c++) {

        Cell* cell = &cells[c];
        
        for (int idx = 0; idx < cell->particleIdxs.size(); ++idx) {
        
            int i = cell->particleIdxs[idx];
            Particle* particleA = &particles[i];

            Vec2 fpress(0.f, 0.f);
            Vec2 fvisc(0.f, 0.f);

            // Loop over neighbor cells
            for (int k = 0; k < 9; k++) {

                int neighborCellIdx = cell->neighbourIdxs[k];
                if (neighborCellIdx == -1) continue;
                
                Cell* neighborCell = &cells[neighborCellIdx];

                for (int jdx = 0; jdx < neighborCell->particleIdxs.size(); ++jdx) {

                    int j = neighborCell->particleIdxs[jdx];
                    if (i == j) continue;
                    
                    Particle* particleB = &particles[j];
                    Vec2 d = (particleA->pos - particleB->pos);
                    
                    double rSq = d.sqrMagnitude();  // only perform expensive sqrt-op if in kernal radius
                    if (rSq < kernalRadiusSq) {
                        double r = sqrt(rSq);
                        double diff = kernelRadius - r;

                        Vec2 deltaVel = particleB->velocity - particleA->velocity;
                        fpress += d.unit() * (particleB->mass * (particleA->pressure + particleB->pressure) / (2.f * particleB->density) * SPIKY_GRAD * (diff * diff * diff));
                        fvisc += deltaVel * ((viscosity * particleB->mass) / particleB->density * VISC_LAP * diff);
                    }
                }
            }
            particleA->force = fpress + fvisc + particleA->forceExt;
        }
    }
}



void SPHFluidSimulation2D::computeDensityPressure() {

#pragma omp parallel for
    for (int c = 0; c < nCells; c++)    // loop through each cell
    {
		Cell* cell = &cells[c];


        for (int i = 0; i < cell->particleIdxs.size(); i++) // loop through each particle in current cell
        {

			int particleAIdx = cell->particleIdxs[i];
            Particle* particleA = &particles[ particleAIdx ];
            particleA->density = 0.0f;


            // loop through all particles in current + neighbouring cells
            
            for (int k = 0; k < 9; k++) // no neighbour idxs
            {
        
				int neighbourCellIdx = cell->neighbourIdxs[k];
                if (neighbourCellIdx == -1) { continue; } // skip if no neighbour

				Cell* neighbourCell = &cells[neighbourCellIdx];

                for (int j = 0; j < neighbourCell->particleIdxs.size(); j++)
                {					
                    int particleBIdx = neighbourCell->particleIdxs[j];
                    Particle* particleB = &particles[ particleBIdx ];

                    double d = (particleA->pos - particleB->pos).sqrMagnitude();
                    double diff = kernalRadiusSq - d;

                    if (diff > 0) {
                        particleA->density += particleB->mass * POLY6 * (diff * diff * diff);
                    }
                }
            }

            particleA->pressure = GasConstant * (particleA->density - RestDensity);
        }
    }
}
