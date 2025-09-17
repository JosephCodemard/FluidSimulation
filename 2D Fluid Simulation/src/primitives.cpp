#include "primitives.h"


Vec3 Vec3::abs() {
	return Vec3(std::abs(x), std::abs(y), std::abs(z));
}

void Vec3::_printInfo() {
	std::cout << "Vec3(" << x << "," << y << ", " << z << ")" << std::endl;
}


Vec3 Vec3::unit() {
	return Vec3(x,y,z) / magnitude();
}


double Vec3::dot(const Vec3& v) const {
	return x * v.x + y * v.y + z * v.z; 
}






Vec3& Vec3::operator+=(const Vec3& rhs) {

	this->x += rhs.x;
	this->y += rhs.y;
	this->z += rhs.z;
	return *this;
}

Vec3& Vec3::operator-=(const Vec3& rhs) {

	this->x -= rhs.x;
	this->y -= rhs.y;
	this->z -= rhs.z;
	return *this;
}











Vec2 Vec2::abs() {
	return Vec2(std::abs(x), std::abs(y));
}

Vec2 Vec2::unit() {
	double mag = magnitude();
	return (mag != 0) ? (*this / mag) : Vec2(0, 0);
}

double Vec2::dot(const Vec2& v) const {
	return x * v.x + y * v.y;
}

void Vec2::_printInfo() {
	std::cout << "Vec2(" << x << ", " << y << ")" << std::endl;
}


Vec2& Vec2::operator+=(const Vec2& rhs) {
	this->x += rhs.x;
	this->y += rhs.y;
	return *this;
}

Vec2& Vec2::operator-=(const Vec2& rhs) {
	this->x -= rhs.x;
	this->y -= rhs.y;
	return *this;
}








Triangle::Triangle(Vec3 v1, Vec3 v2, Vec3 v3) {
	vertex1 = v1;
	vertex2 = v2;
	vertex3 = v3;
}



Quad::Quad(Vec3 v1, Vec3 v2, Vec3 v3, Vec3 v4) {
	vertex1 = v1;
	vertex2 = v2;
	vertex3 = v3;
	vertex4 = v4;

	triangles[0] = Triangle(vertex1, vertex2, vertex4);
	//triangles[1] = Triangle(vertex4, vertex3, vertex1);
	triangles[1] = Triangle(vertex2, vertex3, vertex4);
}




Circle::Circle(double rad, Vec3 pos, int n) {

	radius = rad;
	position = pos;
	n_gon = n;
	triangles = {};


	if (n < 3) {
		throw std::invalid_argument("An n-agon must have at least 3 sides.");
	}

	double angleIncrement = 2 * 3.1415926535 / n_gon;

	for (int i = 0; i < n; ++i) {
		double angle1 = i * angleIncrement;
		double angle2 = (i + 1) * angleIncrement;


		Vec3 vert1 = Vec3(
			position.x + radius * cos(angle1),
			position.y + radius * sin(angle1),
			position.z
		);

		Vec3 vert2 = Vec3(
			position.x + radius * cos(angle2),
			position.y + radius * sin(angle2),
			position.z
		);

		triangles.emplace_back(Triangle(pos, vert1, vert2));
	}
}