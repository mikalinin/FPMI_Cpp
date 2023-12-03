#include <iostream>
#include <vector>
#include <initializer_list>
#include <cmath>

struct Point {
  double x;
  double y;

  Point(double x, double y) : x(x), y(y) {}

  friend bool operator==(const Point&, const Point&);
};

bool operator==(const Point& point1, const Point& point2) {
  return point1.x == point2.x && point1.y == point2.y;
}

bool operator!=(const Point& point1, const Point& point2) {
  return !(point1 == point2);
}

class Line {
// Ax + By + C = 0
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;

  Line() = default;

  Line(const Point& point1, const Point& point2) {
    a = point1.y - point2.y;
    b = point2.x - point1.x;
    c = point1.x * point2.y - point2.x * point1.y;
  }

  Line(const Point& point, double k) : a(-k), b(1) {
    c = point.x * k - point.y;
  }

  Line(double k, double b) : a(-k), b(1), c(-b) {}

  friend bool operator==(const Line&, const Line&);
};

bool operator==(const Line& line1, const Line& line2) {
  return line1.a == line2.a && line1.b == line2.b && line1.c == line2.c;
}

bool operator!=(const Line& line1, const Line& line2) {
  return !(line1 == line2);
}

class Vector {
  Point a = Point(0.0, 0.0);
  Point b = Point(0.0, 0.0);
public:
  Vector(const Point& a, const Point& b) : a(a), b(b) {}

  double len() const { return pow((a.x + b.x) * (a.x + b.x) + (a.y + b.y) * (a.y + b.y), 0.5); }
};

class Shape {

};

class Polygon : Shape {
  std::vector<Point> vertices_;
  int vertices_count_;
public:
  Polygon(std::vector<Point>& vertices) : vertices_(vertices), vertices_count_(static_cast<int>(vertices.size())) {}

  Polygon(std::initializer_list<Point> vertices) : vertices_(vertices),
                                                   vertices_count_(static_cast<int>(vertices.size())) {}

  int verticesCount() const { return vertices_count_; };

  std::vector<Point> getVertices() { return vertices_; };

  bool isConvex() {

  }
};