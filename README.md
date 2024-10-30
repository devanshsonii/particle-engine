
# Particle Simulation with Energy Conservation

This project is a simple particle simulation built using [raylib](https://www.raylib.com/) and C++. It simulates particles under gravitational forces with collisions and boundary checks while testing for energy conservation. Additionally, "wall particles" can be added, which act as immovable obstacles that interact with other particles.

## Features

- **Basic Physics Simulation**: Particles experience gravity and can collide with each other elastically.
- **Energy Conservation Check**: Calculates and checks for conservation of kinetic and potential energy within a specified tolerance.
- **Boundary Collisions**: Particles are constrained within screen boundaries.
- **Wall Particles**: Particles that act as static obstacles, causing other particles to bounce off them.
- **RK4 Integration**: Uses the Runge-Kutta 4th order method for stable and accurate physics integration.

## Project Structure

- **`main.cpp`**: Contains all definitions, including physics integration, collision handling, drawing functions, and an energy conservation test.

## Installation and Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/devanshsonii/particle-engine
   cd particle-engine
   ```

2. **Install raylib**: Make sure you have [raylib](https://www.raylib.com/) installed on your system.

3. **Compile and Run**:
   ```bash
   g++ -o particle_simulation main.cpp -lraylib -std=c++17 -O2
   ./particle_simulation
   ```

## Controls

- **Left Mouse Button**: Adds a new particle at the mouse position.
- **Key 'W'**: Adds a wall particle at the mouse position, creating a static obstacle.

## Code Overview

- **Classes**:
  - `Particle`: Represents a dynamic particle with position, velocity, and mass. Handles updates and collisions with other particles and boundaries.
  - `WallParticle`: Represents a static particle that interacts with dynamic particles but does not move.
- **Functions**:
  - `Update(float dt)`: Updates all particles' positions and resolves collisions.
  - `Draw()`: Draws all particles and wall particles on the screen.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

Thanks to the raylib team for providing an easy-to-use library for 2D graphics in C++.
