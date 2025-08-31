#pragma once

#include <glad/glad.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "primitives.h"


// Taken from https://learnopengl.com/code_viewer_gh.php?code=includes/learnopengl/shader_s.h,
// because why re-invent the wheel...

class Shader
{
public:
    unsigned int ID;
    Shader(const char* vertexPath, const char* fragmentPath);

    void compile();
    void use();

    void setBool(const std::string& name, bool value) const;
    void setInt(const std::string& name, int value) const;
    void setFloat(const std::string& name, float value) const;
    void setVec3(const std::string& name, Vec3 value) const;

private:

    const char* vertexPath;
    const char* fragmentPath;

    void checkCompileErrors(unsigned int shader, std::string type);
};