#include <iostream>
#include <vector>
#include <cmath>
#include "raylib.h"

const float g = 9.81f;

struct State {
    Vector2 pos; // position
    Vector2 velocity; // velocity
};

Vector2 operator+(const Vector2& v1, const Vector2& v2) {
    return {v1.x + v2.x, v1.y + v2.y};
}

Vector2 operator*(const Vector2& v, float scalar) {
    return {v.x * scalar, v.y * scalar};
}

Vector2 operator*(float scalar, const Vector2& v) {
    return {v.x * scalar, v.y * scalar};
}

Vector2 computeAcceleration(const State &s) {
    return {0, g}; // Gravity points downward
}

State derivative(const State &s) {
    Vector2 acc = computeAcceleration(s);
    return {s.velocity, acc};
}

State rk4(const State& s, float dt) {
    State k1 = derivative(s);
    State k2 = derivative({
        s.pos + k1.pos * (dt * 0.5f),
        s.velocity + k1.velocity * (dt * 0.5f)
    });
    State k3 = derivative({
        s.pos + k2.pos * (dt * 0.5f),
        s.velocity + k2.velocity * (dt * 0.5f)
    });
    State k4 = derivative({
        s.pos + k3.pos * dt,
        s.velocity + k3.velocity * dt
    });
    return {
        s.pos + (k1.pos + 2.0f*k2.pos + 2.0f*k3.pos + k4.pos) * (dt / 6.0f),
        s.velocity + (k1.velocity + 2.0f*k2.velocity + 2.0f*k3.velocity + k4.velocity) * (dt / 6.0f)
    };
}

class Particle {
public:
    float x, y, radius;
    float mass;
    Color color;
    State state;
    Particle(float x, float y, float radius, Color color)
        : x(x), y(y), radius(radius), color(color) {
            state.pos = {x,y};
            state.velocity = {0,0};
            mass = radius * radius * 3.14159f; // Mass proportional to area
    }
    void DrawParticle(){
        DrawCircle(x, y, radius, color);
    }
    void updatePos(Vector2 mouse);
    void applyForce(float forceX, float forceY, float dt);
    void resolveCollision(Particle& other);
    void checkBounds(float dt);
    void update(float dt) {
        State newState = rk4(state, dt);
        state = newState;
        x = state.pos.x;
        y = state.pos.y;
        checkBounds(dt);
    }
};

void Particle::updatePos(Vector2 mouse) {
    x = mouse.x;
    y = mouse.y;
    state.pos = {x, y};
    state.velocity = {0, 0};
}

void Particle::applyForce(float forceX, float forceY, float dt) {
    state.velocity.x += forceX * dt / mass;
    state.velocity.y += forceY * dt / mass;
}

void Particle::resolveCollision(Particle &other) {
    Vector2 delta = {other.x - x, other.y - y};
    float distance = sqrt(delta.x * delta.x + delta.y * delta.y);
    float minDist = radius + other.radius;

    if (distance < minDist) {
        Vector2 normal = {delta.x / distance, delta.y / distance};
        Vector2 relativeVelocity = {
            other.state.velocity.x - state.velocity.x,
            other.state.velocity.y - state.velocity.y
        };

        float normalVelocity = relativeVelocity.x * normal.x + relativeVelocity.y * normal.y;

        if (normalVelocity > 0) return;

        float restitution = 1.0f; // Perfectly elastic collisions
        float impulseScalar = -(1 + restitution) * normalVelocity;
        impulseScalar /= 1/mass + 1/other.mass;

        Vector2 impulse = {impulseScalar * normal.x, impulseScalar * normal.y};

        state.velocity.x -= impulse.x / mass;
        state.velocity.y -= impulse.y / mass;
        other.state.velocity.x += impulse.x / other.mass;
        other.state.velocity.y += impulse.y / other.mass;

        // Separate particles to prevent overlap
        float correction = (minDist - distance) / 2.0f;
        Vector2 correctionVector = {normal.x * correction, normal.y * correction};
        x -= correctionVector.x;
        y -= correctionVector.y;
        other.x += correctionVector.x;
        other.y += correctionVector.y;

        state.pos = {x, y};
        other.state.pos = {other.x, other.y};
    }
}

void Particle::checkBounds(float dt) {
    if (y + radius >= GetScreenHeight()) {
        y = GetScreenHeight() - radius;
        state.velocity.y *= -1.0f; 
    }
    if (y - radius <= 0) {
        y = radius;
        state.velocity.y *= -1.0f; 
    }
    if (x + radius >= GetScreenWidth()) {
        x = GetScreenWidth() - radius;
        state.velocity.x *= -1.0f; 
    }
    if (x - radius <= 0) {
        x = radius;
        state.velocity.x *= -1.0f; 
    }
    state.pos = {x, y};
}

std::vector<Particle> particles;

class WallParticle {
public:
    float x, y, radius;
    float mass;
    Color color;
    WallParticle(float x, float y, float radius, Color color)
        : x(x), y(y), radius(radius), color(color) {
            mass = radius * radius * 3.14159f; // Mass proportional to area
    }
    void DrawParticle(){
        DrawCircle(x, y, radius, color);
    }
    // all other particles bounce off it, this particle is a wall and will not move 
    void resolveCollision(Particle& other);

};

void WallParticle::resolveCollision(Particle &other) {
    Vector2 delta = {other.x - x, other.y - y};
    float distance = sqrt(delta.x * delta.x + delta.y * delta.y);
    float minDist = radius + other.radius;

    if (distance < minDist) {
        Vector2 normal = {delta.x / distance, delta.y / distance};
        Vector2 relativeVelocity = other.state.velocity;

        float normalVelocity = relativeVelocity.x * normal.x + relativeVelocity.y * normal.y;

        if (normalVelocity > 0) return;

        float restitution = 1.0f; // Perfectly elastic collisions
        float impulseScalar = -(1 + restitution) * normalVelocity;
        impulseScalar /= 1 / other.mass;

        Vector2 impulse = {impulseScalar * normal.x, impulseScalar * normal.y};

        other.state.velocity.x += impulse.x / other.mass;
        other.state.velocity.y += impulse.y / other.mass;

        // Separate particles to prevent overlap
        float correction = (minDist - distance) / 2.0f;
        Vector2 correctionVector = {normal.x * correction, normal.y * correction};
        other.x += correctionVector.x;
        other.y += correctionVector.y;

        other.state.pos = {other.x, other.y};
    }
}

std::vector<WallParticle> wall_particles;



float calculateTotalEnergy(const std::vector<Particle>& particles) {
    float totalEnergy = 0.0f;
    for (const auto& p : particles) {
        float kineticEnergy = 0.5f * p.mass * (p.state.velocity.x * p.state.velocity.x + 
                                               p.state.velocity.y * p.state.velocity.y);
        float potentialEnergy = p.mass * g * (GetScreenHeight() - p.y);
        totalEnergy += kineticEnergy + potentialEnergy;
    }
    return totalEnergy;
}

void testEnergyConservation(int numSteps, float dt) {
    float initialEnergy = calculateTotalEnergy(particles);
    float maxEnergyDifference = 0.0f;

    for (int i = 0; i < numSteps; i++) {
        for (auto& pen : particles) {
            pen.update(dt);
        }
        if (particles.size() >= 2) {
            for (size_t i = 0; i < particles.size(); i++) {
                for (size_t j = i + 1; j < particles.size(); j++) {
                    particles[i].resolveCollision(particles[j]);
                }
            }
        }
        float currentEnergy = calculateTotalEnergy(particles);
        float energyDifference = std::abs(currentEnergy - initialEnergy) / initialEnergy;
        maxEnergyDifference = std::max(maxEnergyDifference, energyDifference);
    }

    std::cout << "Maximum energy difference: " << maxEnergyDifference * 100 << "%" << std::endl;
    if (maxEnergyDifference < 0.01) {
        std::cout << "Energy is well conserved (less than 1% difference)" << std::endl;
    } else {
        std::cout << "Energy conservation error exceeds 1%" << std::endl;
    }
}

void Draw() {
    for (auto& pen : particles) {
        pen.DrawParticle();
    }
    for (auto& wall_pen : wall_particles) {
        wall_pen.DrawParticle();
    }
}

void Update(float dt) {
    Vector2 mouse = GetMousePosition();
    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
        particles.emplace_back(mouse.x, mouse.y, 25, RED);
    }
    if (IsKeyDown(KEY_W)){
        wall_particles.emplace_back(mouse.x, mouse.y, 5, BLUE);
    }
    for (auto& pen : particles) {
        pen.update(dt);
    }
    for(auto &wall_pen : wall_particles){
        for(auto &pen : particles){
            wall_pen.resolveCollision(pen);
        }
    }
    if (particles.size() >= 2) {
        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = i + 1; j < particles.size(); j++) {
                particles[i].resolveCollision(particles[j]);
            }
        }
    }
    // std::cout << calculateTotalEnergy(particles) << "\n";
}

int main() {
    const int screenHeight = 800;
    const int screenWidth = 800;
    const float dt = 1.0f / 60.0f;
    SetTargetFPS(90);
    InitWindow(screenWidth, screenHeight, "Particle");
    // testEnergyConservation(10000, dt);
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);
        Update(dt);
        Draw();
        EndDrawing();
    }
    CloseWindow();
    return 0;
}