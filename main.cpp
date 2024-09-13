#include <iostream>
#include <vector>
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
    return {0, -g};
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
        s.pos + (k1.pos + 2*k2.pos + 2*k3.pos + k4.pos) * (dt / 6.0f),
        s.velocity + (k1.velocity + 2*k2.velocity + 2*k3.velocity + k4.velocity) * (dt / 6.0f)
    };
}


class Pendulum {
public:
    float x, y, length, radius;
    Color color;
    State state;
    Pendulum(float x, float y, float radius, float length, Color color)
        : x(x), y(y), radius(radius), length(length), color(color) {
            state.pos = {x,y};
            state.velocity = {0,0};
    }
    void DrawPendulum(){
        DrawCircle(x, y, radius, color);
    }
    void updatePos(Vector2 mouse);
    void applyGravity(float dt);
    void applyForce(float forceX, float forceY, float dt);
    void resolveCollision(Pendulum& other);
        void checkBounds(float dt) {
        // Apply gravity
        state.velocity.y += g * dt;
        state.pos.y += state.velocity.y * dt;
        state.pos.x += state.velocity.x * dt; 
        x = state.pos.x;
        y = state.pos.y;

        // Check for collision with screen bounds
        if (y + radius >= GetScreenHeight()) {
            y = GetScreenHeight() - radius;
            state.velocity.y = -(state.velocity.y * 0.8); // bounce back + dampening 
        }
        if (y - radius <= 0) {
            y = radius;
            state.velocity.y = -(state.velocity.y * 0.8); // bounce back + dampening 
        }
        if (x + radius >= GetScreenWidth()) {
            x = GetScreenWidth() - radius;
            state.velocity.x = -(state.velocity.x * 0.8); // bounce back + dampening 
        }
        if (x - radius <= 0) {
            x = radius;
            state.velocity.x = -(state.velocity.x * 0.8); // bounce back + dampening 
        }

        // Update positions in state
        state.pos = {x, y};
    }

};
void Pendulum :: updatePos(Vector2 mouse) {
    x = mouse.x;
    y = mouse.y;
    state.pos.x = x;
    state.pos.y = y;
    state.velocity.x = 0;
    state.velocity.y = 0;
}

void Pendulum :: applyGravity(float dt) {
    if(y + radius >= GetScreenHeight()){
        state.velocity.y = -(state.velocity.y * 0.8); // bounce back + dampening 
    }
    state.velocity.y += g * dt;
    state.pos.y += state.velocity.y * dt;
    state.pos.x += state.velocity.x * dt; 
    x = state.pos.x;
    y = state.pos.y;
}

void Pendulum :: applyForce(float forceX, float forceY, float dt){
    state.velocity.x += forceX * dt;
    state.pos.x += state.velocity.x * dt;
    state.velocity.y += forceY * dt;
    state.pos.y += state.velocity.y * dt;
    x = state.pos.x;
    y = state.pos.y;
}

void Pendulum :: resolveCollision(Pendulum &other){
    Vector2 delta = {other.x - x, other.y - y};
    float distance = sqrt(delta.x * delta.x + delta.y * delta.y);
    float minDist = radius + other.radius;

    if (distance < minDist) {
        float overlap = 0.5f * (distance - minDist);

        // Displace current pendulum
        x -= overlap * (x - other.x) / distance;
        y -= overlap * (y - other.y) / distance;

        // Displace other pendulum
        other.x += overlap * (x - other.x) / distance;
        other.y += overlap * (y - other.y) / distance;

        // Update positions in state
        state.pos = {x, y};
        other.state.pos = {other.x, other.y};

        // Resolve velocities
        Vector2 normal = {delta.x / distance, delta.y / distance};
        Vector2 relativeVelocity = {state.velocity.x - other.state.velocity.x, state.velocity.y - other.state.velocity.y};
        float velocityAlongNormal = relativeVelocity.x * normal.x + relativeVelocity.y * normal.y;

        if (velocityAlongNormal > 0) return;

        float restitution = 0.8f; // coefficient of restitution
        float j = -(1 + restitution) * velocityAlongNormal;
        j /= 1 / radius + 1 / other.radius;

        Vector2 impulse = {j * normal.x, j * normal.y};
        state.velocity.x -= 1 / radius * impulse.x;
        state.velocity.y -= 1 / radius * impulse.y;
        other.state.velocity.x += 1 / other.radius * impulse.x;
        other.state.velocity.y += 1 / other.radius * impulse.y;
    }
}




std::vector<Pendulum> pendulums;

void Draw() {
    for (auto& pen : pendulums) {
        pen.DrawPendulum();
    }
}

void Update(float dt) {
    Vector2 mouse = GetMousePosition();
    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
        pendulums.emplace_back(mouse.x, mouse.y, 25, 100, RED);
    }
    for (auto& pen : pendulums) {
        // pen.applyGravity(dt);
        pen.checkBounds(dt);
    }
    if(pendulums.size() >= 2){
        for(int i = 0; i < pendulums.size(); i++){
            for(int j = i+1; j < pendulums.size(); j++){
                pendulums[i].resolveCollision(pendulums[j]);
            }
        }
    }
}

int main() {
    const int screenHeight = 800;
    const int screenWidth = 800;
    const float dt = 1.0f / 60.0f;
    SetTargetFPS(60);
    InitWindow(screenWidth, screenHeight, "Pendulum");
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