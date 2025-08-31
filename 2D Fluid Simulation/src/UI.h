#pragma once


#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <ft2build.h>
#include <iostream>
#include <map>
#include <chrono>


#include "UI.h"
#include "primitives.h"
#include "shader.h"

#include FT_FREETYPE_H  


struct Character {
    unsigned int TextureID; // ID handle of the glyph texture
    Vec3   Size;      // Size of glyph
    Vec3   Bearing;   // Offset from baseline to left/top of glyph
    unsigned int Advance;   // Horizontal offset to advance to next glyph
};


class FPSCounter {
private:
    int frames;
    int currentFrame;
    double fps;
    std::chrono::steady_clock::time_point lastTime;

public:
    FPSCounter();
    void step();
    double getFPS() const;
	double getFrameTime() const;
	int getCurrentFrame() const;
};


class Window {

public:

    bool wireframe = false;

    Window(unsigned int height, unsigned int width);
    ~Window();

    bool Init();

    virtual bool OnUserCreate();
    virtual bool OnUserUpdate();

    void SetScaleFactor(Vec3 scaleFactor);
    Vec3 GetScaleFactor();

    void DrawTriangle(Triangle tri, Vec3 col=Vec3(1,0,0));
    void DrawCircle(Circle circ, Vec3 col = Vec3(1, 0, 0));
    void DrawQuad(Quad quad, Vec3 col = Vec3(1, 0, 0));
	void DrawWireframeQuad(Quad quad, Vec3 col = Vec3(1, 0, 0), float edgeThickness = 1.0f);

    void DrawText(std::string text, float x, float y, float scale, Vec3 color);

    bool IsScrollDown();
    bool IsScrollUp();

    Vec3 GetMousePosition();
    Vec3 GetScreenSize();

    Vec3 GetScreenOrigin();
    void SetScreenOrigin(Vec3 origin);

	void SetTitle(std::string title);

private:

    GLFWwindow* window;

    unsigned int SCR_WIDTH;
    unsigned int SCR_HEIGHT;

    Vec3 windowScaleFactor;
    Vec3 userScaleFactor;
    Vec3 sceneOrigin;

    int scrollDirection = 0;

    Shader solidShader = Shader("./src/shaders/solid_shader.vs", "./src/shaders/solid_shader.fs");
    Shader fontShader = Shader("./src/shaders/font_shader.vs", "./src/shaders/font_shader.fs");


    std::map<GLchar, Character> Characters;
    unsigned int VAO, VBO = 0;



    unsigned int fragmentShader = 0;
    unsigned int vertexShader = 0;
    unsigned int shaderProgram = 0;

    GLint uniform_varColour;
    Vec3 mousePos = Vec3();


    bool initGlfw();
    bool initGlad();
    bool initFreetype();


    static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);

    bool terminate();

    void processInput(GLFWwindow* window);
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
};