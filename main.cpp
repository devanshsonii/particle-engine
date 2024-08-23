#include <iostream>
#include "raylib.h"


class Circle {
public: 
    int x, y;
    float radius;
    Color color;
    
};

void Draw(){
    DrawCircle(200, 200, 25, WHITE); // center point

}

void Update(){

}


int main(){
    const int screenHeight = 800;
    const int screenWidth = 800;

    SetTargetFPS(60);
    InitWindow(screenWidth, screenHeight, "Pendulum");
    while(!WindowShouldClose()){
        BeginDrawing();
        ClearBackground(BLACK);
        Update();
        Draw();
        EndDrawing();
    }
    CloseWindow();
}
