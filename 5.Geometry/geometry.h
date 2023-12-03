#include <iostream>
#include <vector>
#include <initializer_list>
#include <cmath>
#define PI 3.14

struct Point {
  double x;
  double y;

  Point(double x, double y) : x(x), y(y) {}

  friend bool operator==(const Point&, const Point&);

  friend double points_distance(const Point&, const Point&);
};

bool operator==(const Point& point1, const Point& point2) {
  return point1.x == point2.x && point1.y == point2.y;
}

bool operator!=(const Point& point1, const Point& point2) {
  return !(point1 == point2);
}

double points_distance(const Point& point1, const Point& point2) {
  return pow((point1.x - point2.x) * (point1.x - point2.x) + (point1.y - point2.y) * (point1.y - point2.y), 0.5);
}

class Line {
// Ax + By + C = 0
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
public:
  Line() = default;

  Line(double a, double b, double c) : a(a), b(b), c(c) {}

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
  double x = 0.0;
  double y = 0.0;
public:
  Vector(double x, double y) : x(x), y(y) {}

  Vector(const Point& a, const Point& b) : x(a.x - b.x), y(a.y - b.y) {}

  friend double scalarProduct(const Vector&, const Vector&);
};

double scalarProduct(const Vector& vec1, const Vector& vec2) {
  return vec1.x * vec2.x + vec1.y * vec2.y;
}

class Shape {

};

class Polygon : Shape {
protected:
  std::vector<Point> vertices_;
  int vertices_count_;
public:
  Polygon(std::vector<Point>& vertices) : vertices_(vertices), vertices_count_(static_cast<int>(vertices.size())) {}

  Polygon(std::initializer_list<Point> vertices) : vertices_(vertices),
                                                   vertices_count_(static_cast<int>(vertices.size())) {}

  int verticesCount() const { return vertices_count_; };

  std::vector<Point> getVertices() { return vertices_; };

  bool isConvex() {
    for (int i = 0; i < vertices_count_; ++i) {
      Vector vec1 = Vector(vertices_[i], vertices_[(i + 1) % vertices_count_]);
      Vector vec2 = Vector(vertices_[(i + 1) % vertices_count_], vertices_[(i + 2) % vertices_count_]);
      if (scalarProduct(vec1, vec2) < 0) {
        return false;
      }
    }
    return true;
  }
};

class Ellipse : Shape {
  // x^2/a^2 + y^2/b^2 = 1;
protected:
  double a = 0.0;
  double b = 0.0;
  std::pair<Point, Point> focus = std::pair(Point(0, 0), Point(0, 0));
public:
  Ellipse(const Point& point1, const Point& point2, double dist) {
    focus = std::make_pair(point1, point2);
    a = dist / 2;
    b = pow((a * a - (points_distance(point1, point2) / 2) * (points_distance(point1, point2) / 2)), 0.5);
  }

  std::pair<Point, Point> focuses() const { return focus; }

  std::pair<Line, Line> directrices() const {
    Line line1 = Line(1, 0, a / eccentricity());
    Line line2 = Line(1, 0, (-1 * a) / eccentricity());
    return std::make_pair(line1, line2);
  }

  double eccentricity() const { return pow((1 - (a * a / b * b)), 0.5); }

  Point center() const {
    return Point((focus.first.x + focus.second.x) / 2, (focus.first.y + focus.second.y) / 2);
  }
};

class Circle : public Ellipse {
  Circle(const Point& center, double radius) : Ellipse(center, center, 2 * radius) {}

  double radius() const {return a; }
};

class Rectangle : public Polygon {
  Rectangle(const Point& point1, const Point& point2, double k) : Polygon{point1, Point(0, 0), point2, Point(0, 0)} {

    double angle = 180 - 2 * atan(k) * 180 / PI;
    Point middle = middle_point(point1, point2);
    vertices[1] = point1.rotate(middle, angle);
    vertices[3] = point2.rotate(middle, angle);
  }
};