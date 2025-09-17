#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>



class Vec3
{
public:
    double x;
    double y;
    double z;
       
    Vec3(double _x=0, double _y=0, double _z=0) : x(_x), y(_y), z(_z) {}

    Vec3 abs();
    Vec3 unit();
    double dot(const Vec3& v) const;


    void _printInfo();

    Vec3& operator+=(const Vec3& rhs);
    Vec3& operator-=(const Vec3& rhs);


    // inline these functions
    float magnitude()  const { return std::sqrt(x * x + y * y + z * z); }
    float sqrMagnitude()  const { return x * x + y * y + z * z; }

    Vec3 operator+(const Vec3& rhs) const { return Vec3(x + rhs.x, y + rhs.y, z + rhs.z); }
    Vec3 operator-(const Vec3& rhs) const { return Vec3(x - rhs.x, y - rhs.y, z - rhs.z); }
    Vec3 operator/(const double& rhs) const { return Vec3(x / rhs, y / rhs, z / rhs); }
    Vec3 operator*(const double& rhs) const { return Vec3(x * rhs, y * rhs, z * rhs); }

};


class Vec2
{
public:
    double x;
    double y;

    Vec2(double _x=0, double _y=0) : x(_x), y(_y) {}

    Vec2 abs();
    Vec2 unit();
    
    double dot(const Vec2& v) const;
    void _printInfo();

    Vec2& operator+=(const Vec2& rhs);
    Vec2& operator-=(const Vec2& rhs);


    // inline these functions
    float magnitude() const { return std::sqrt(x * x + y * y); }
    float sqrMagnitude() const { return x * x + y * y; }

    Vec2 operator+(const Vec2& rhs) const { return Vec2(x + rhs.x, y + rhs.y); }
    Vec2 operator-(const Vec2& rhs) const { return Vec2(x - rhs.x, y - rhs.y); }
    Vec2 operator/(double rhs) const { return Vec2(x / rhs, y / rhs); }
    Vec2 operator*(double rhs) const { return Vec2(x * rhs, y * rhs); }
};





class Triangle
{
public:
    Vec3 vertex1;
    Vec3 vertex2;
    Vec3 vertex3;

    Triangle(Vec3 v1 = Vec3(0,0,0), Vec3 v2 = Vec3(0, 0, 0), Vec3 v3 = Vec3(0, 0, 0));
};




class Quad
{
public:
    Vec3 vertex1;
    Vec3 vertex2;
    Vec3 vertex3;
    Vec3 vertex4;

    Triangle triangles[2];

    Quad(Vec3 v1 = Vec3(1, 0, 0), Vec3 v2 = Vec3(1, 1, 0), Vec3 v3 = Vec3(0, 1, 0), Vec3 v4 = Vec3(0, 0, 0));
};





class Circle
{
public:
    double radius;
    int n_gon;
    Vec3 position;


    std::vector<Triangle> triangles;

    Circle(double rad, Vec3 pos, int n);
};