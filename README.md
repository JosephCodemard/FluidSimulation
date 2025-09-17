# 2D Fluid Simulation

A particle-based fluid dynamics simulation
This project explores computational fluid dynamics (CFD) by implementing a particle based solver, using cell based optimiations, as well as 
OpenMP to parallelise simple loops in order to increase simulation speed. 
Currently, it implements a particle-based solver, while an Eulerian grid-based version is under active development.

***

### Features

- **Particle-Based Fluid Solver** — real-time particle interaction with fluid-like motion.  
- **Custom Rendering Loop** — built using OpenGL (and GLFW) with a custom, extendable rendering/game loop
- **Cell-Based Optimization** — spatial partitioning to efficiently handle particle-particle interactions.  
- **Multithreading Support** — uses openMP for better scalability across CPU cores for increased simulation speed.  

***

### Future Pipeline

This project is still in active development. Planned improvements include:

- [ ] **Interactive Controls**  
  Configure parameters (viscosity, density, pressure, gravity, etc.) live during simulation
- [ ] **Freetype UI Integration**  
  On-screen text rendering for real-time statistics and adjustable settings.  
- [ ] **Performance Enhancements**  
  Improved multi-core CPU performance, and GPU acceleration via CUDA/OpenMP.
- [ ] **Support for Irregular Solid Boundaries**  
  Collisions with non-uniform obstacles inside simulation bounds for more realism.
- [ ] **Eulerian Solver Integration**  
  An alternate grid-based fluid model for stability and accuracy.  
- [ ] **WebAssembly Export**  
  Build cross-compiled simulations for the web using Emscripten.
- [ ] **3D Support**  
  Support 3D simulations, while mainting a stable simulation speed.
