#include <iostream>
#include "raylib.h"

const float g = 9.81f;
struct State {
    float y; // height 
    float v; // velocity
};

State derivative(const State &s){
    return {
        s.v, // dy/dt = v
        g // dv/dt = g
    };
}


State rk4(const State& s, float dt) {
    State k1 = derivative(s);
    State k2 = derivative({
        s.y + k1.y * dt * 0.5f,
        s.v + k1.v * dt * 0.5f
    });
    State k3 = derivative({
        s.y + k2.y * dt * 0.5f,
        s.v + k2.v * dt * 0.5f
    });
    State k4 = derivative({
        s.y + k3.y * dt,
        s.v + k3.v * dt
    });

    return {
        s.y + (k1.y + 2*k2.y + 2*k3.y + k4.y) * dt / 6.0f,
        s.v + (k1.v + 2*k2.v + 2*k3.v + k4.v) * dt / 6.0f
    };
}

class Circle {
public: 
    int x, y;
    float radius;
    Color color;
    State state;
    Circle(float x, float y, float radius, Color color)
        : x(x), y(y), radius(radius), color(color) {
        state.y = y;
        state.v = 0;  // Initial velocity
    }
};

Circle one(200, 200, 25, RED);


void Draw(){
    DrawCircle(one.x, one.y, one.radius, one.color);
}


void Update(float dt){
    one.state = rk4(one.state, dt);
    one.y = one.state.y;
    if(one.y > GetScreenHeight() - one.radius){
        one.y = GetScreenHeight() - one.radius;
        one.state.y = one.y;
        one.state.v = -one.state.v;
    }
    
}


int main(){
    const int screenHeight = 800;
    const int screenWidth = 800;
    const float dt = 1.0f/60.0f;
    SetTargetFPS(60);
    InitWindow(screenWidth, screenHeight, "Pendulum");
    while(!WindowShouldClose()){
        BeginDrawing();
        ClearBackground(BLACK);
        Update(dt);
        Draw();
        EndDrawing();
    }
    CloseWindow();
}
