#include <iostream>
#include <vector>
#include <initializer_list>
#include <cmath>

#define NUM_PI 3.14159265

class Line;

struct Point {
  double x;
  double y;

  Point() = default;

  Point(double x, double y) : x(x), y(y) {}

  Point rotate(const Point& center, double angle) const {
    angle = (angle / 180) * NUM_PI;
    double x = center.x + cos(angle) * (this->x - center.x) - sin(angle) * (this->y - center.y);
    double y = center.y + sin(angle) * (this->x - center.x) + cos(angle) * (this->y - center.y);
    return Point(x, y);
  }

  void reflect(const Point& center) {
    x = 2 * center.x - x;
    y = 2 * center.y - y;
  }

  void reflect(const Line& axis);

  void scale(const Point& center, double coefficient) {
    x = (x - center.x) * coefficient + center.x;
    y = (y - center.y) * coefficient + center.y;
  }
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

Point middle_point(const Point& point1, const Point& point2) {
  return Point((point1.x + point2.x) / 2, (point1.y + point2.y) / 2);
}


class Line {
// Ax + By + C = 0
public:
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;

  Line() = default;

  Line(double a, double b, double c) : a(a), b(b), c(c) {}

  Line(const Point& point1, const Point& point2) : a(point1.y - point2.y), b(point2.x - point1.x),
                                                   c(point1.x * point2.y - point2.x * point1.y) {}

  Line(const Point& point, double k) : a(-k), b(1), c(point.x * k - point.y) {}

  Line(double k, double b) : a(-k), b(1), c(-b) {}

};

bool operator==(const Line& line1, const Line& line2) {
  return (line1.a * line2.b == line1.b * line2.a && line1.b * line2.c == line1.c * line2.b);
}

bool operator!=(const Line& line1, const Line& line2) {
  return !(line1 == line2);
}

Line perpendicular(const Line& line) {
  return Line(line.b, -line.a, 0);
}

Point cross(const Line& line1, const Line& line2) {
  Point tmp;
  tmp.x = (line1.b * line2.c - line2.b * line1.c) / (line1.a * line2.b - line2.a * line1.b);
  tmp.y = (line2.a * line1.c - line1.a * line2.c) / (line1.a * line2.b - line2.a * line1.b);
  return tmp;
}

Line middle_perpendicular(const Point& point1, const Point& point2) {
  Point middle = middle_point(point1, point2);
  Line axis(point1, point2);
  Line result = perpendicular(axis);
  result.c = -result.a * middle.x - result.b * middle.y;
  return result;
}

void Point::reflect(const Line& axis) {
  Line perpendicular1 = perpendicular(axis);
  perpendicular1.c = -perpendicular1.a * x - perpendicular1.b * y;
  scale(cross(perpendicular1, axis), -1);
}

class Vector {
public:
  double x = 0.0;
  double y = 0.0;

  Vector(double x, double y) : x(x), y(y) {}

  Vector(const Point& a, const Point& b) : x(a.x - b.x), y(a.y - b.y) {}

};

double ScalarProduct(const Vector& vec1, const Vector& vec2) {
  return vec1.x * vec2.x + vec1.y * vec2.y;
}

double VectorMultiply(const Vector& vec1, const Vector& vec2) {
  return vec1.x * vec2.y - vec2.x * vec1.y;
}

class Shape {
public:
  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool operator==(const Shape& another) const = 0;

  virtual bool operator!=(const Shape& another) const = 0;

  virtual bool isCongruentTo(const Shape&) const = 0;

  virtual bool isSimilarTo(const Shape&) const = 0;

  virtual bool containsPoint(Point point) const = 0;

  virtual void rotate(Point center, double angle) = 0;

  virtual void reflect(Point center) = 0;

  virtual void reflect(Line axis) = 0;

  virtual void scale(Point center, double k) = 0;

  virtual ~Shape() = 0;
};

Shape::~Shape() {}

class Segment {
public:
  Point p1;
  Point p2;
  Segment() = default;

  Segment(double x1, double y1, double x2, double y2) {
    p1 = Point(x1, y1);
    p2 = Point(x2, y2);
  }
  Segment(const Point &p1, const Point &p2) : p1(p1), p2(p2) {}
  bool ContainsPoint(const Point& point) const {
    Vector point_to_first(p1, point);
    Vector point_to_second(p2, point);
    return (VectorMultiply(point_to_first, point_to_second) == 0
            && ScalarProduct(point_to_first, point_to_second) <= 0);
  }
};

class Polygon : public Shape {
protected:
  std::vector<Point> vertices_;
public:
  Polygon() = default;

  Polygon(std::vector<Point>& vertices) : vertices_(vertices) {}

  template<class... Types>
  Polygon(Types... args) {

  }

  int verticesCount() const { return static_cast<int>(vertices_.size()); };

  std::vector<Point> getVertices() { return vertices_; };

  bool isConvex() {
    for (int i = 0; i < verticesCount(); ++i) {
      Vector vec1 = Vector(vertices_[i], vertices_[(i + 1) % verticesCount()]);
      Vector vec2 = Vector(vertices_[(i + 1) % verticesCount()], vertices_[(i + 2) % verticesCount()]);
      if (ScalarProduct(vec1, vec2) < 0) {
        return false;
      }
    }
    return true;
  }

  double perimeter() const override {
    double perimeter = 0;
    int size = verticesCount();
    for (int i = 0; i < size; ++i) {
      perimeter += points_distance(vertices_[i], vertices_[(i + 1) % size]);
    }
    return perimeter;
  }

  double area() const override {
    double area = 0;
    int size = verticesCount();
    for (int i = 0; i < size; ++i) {
      area += vertices_[i].x * vertices_[(i + 1) % size].y - vertices_[i].y * vertices_[(i + 1) % size].x;
    }
    return std::abs(area / 2);
  }

  bool operator==(const Shape& another) const override {
    auto polygon_ptr = dynamic_cast<const Polygon*>(&another);
    if (!polygon_ptr) return false;
    if (verticesCount() != polygon_ptr->verticesCount()) return false;
    Point first = polygon_ptr->vertices_[0];
    for (int i = 0; i < verticesCount(); ++i) {
      if (vertices_[i] == first) {
        bool is_equal = true;
        for (int j = 0; j < verticesCount(); ++j) {
          if (polygon_ptr->vertices_[j] != vertices_[(i + j) % verticesCount()]) {
            is_equal = false;
            break;
          }
        }
        if (is_equal) return true;
        is_equal = true;
        for (int j = 0; j < verticesCount(); ++j) {
          int k;
          if (i - j < 0) k = verticesCount() + (i - j);
          else k = i - j;
          if (polygon_ptr->vertices_[j] != vertices_[k]) {
            is_equal = false;
            break;
          }
        }
        if (is_equal) return true;
      }
    }
    return false;
  }

  bool operator!=(const Shape& another) const override {
    return !(*this == another);
  }

  bool isSimilarTo(const Shape& another) const override {
    auto ptr = dynamic_cast<const Polygon*>(&another);
    if (!ptr) return false;
    double k = perimeter() / ptr->perimeter();
    return area() / ptr->area() == k * k;
  }

  bool isCongruentTo(const Shape& another) const override {
    return this->isSimilarTo(another) && this->area() == another.area();
  }

  bool containsPoint(Point point) const override {
    for (int i = 0; i < verticesCount(); ++i) {
      auto seg = Segment(vertices_[i], vertices_[(i + 1) % verticesCount()]);
      if (seg.ContainsPoint(point)) {
        return true;
      }
    }
    bool is_inside = false;
    int j = verticesCount() - 1;
    for (int i = 0; i < verticesCount(); ++i) {
      // false if there is an even number of sides of poly on the left and true with an odd one.
      bool one_direction = vertices_[i].y < point.y && vertices_[j].y >= point.y;
      bool other_direction = vertices_[j].y < point.y && vertices_[i].y >= point.y;
      if (one_direction || other_direction) {
        double ray_poly_cross_x = vertices_[i].x + (point.y - vertices_[i].y) /
                                                   (vertices_[j].y - vertices_[i].y) * (vertices_[j].x - vertices_[i].x);
        if (ray_poly_cross_x < point.x) {
          is_inside = !is_inside;
        }
      }
      j = i;
    }
    return is_inside;
  }

  void rotate(Point center, double angle) override {
    angle = (angle / 180) * NUM_PI;
    for (int i = 0; i < verticesCount(); ++i) {
      vertices_[i] = vertices_[i].rotate(center, angle);
    }
  }

  void reflect(Point center) override {
    scale(center, -1);
  }

  void reflect(Line axis) override {
    for (int i = 0; i < verticesCount(); ++i) {
      vertices_[i].reflect(axis);
    }
  }

  void scale(Point center, double coefficient) override {
    for (int i = 0; i < verticesCount(); ++i) {
      vertices_[i].scale(center, coefficient);
    }

  }
};

class Ellipse : public Shape {
  // x^2/a^2 + y^2/b^2 = 1;
protected:
  double a = 0.0;
  double b = 0.0;
  std::pair<Point, Point> focus = std::pair(Point(0, 0), Point(0, 0));
public:
  Ellipse() = default;

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

  double eccentricity() const { return sqrt(1 - (b * b / a * a)); }

  Point center() const {
    return Point((focus.first.x + focus.second.x) / 2, (focus.first.y + focus.second.y) / 2);
  }

  double perimeter() const override {
    return NUM_PI * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)));
  }

  double area() const override {
    return NUM_PI * a * b;
  }

  bool operator==(const Shape& another) const override {
    auto ellipse_ptr = dynamic_cast<const Ellipse*>(&another);
    if (!ellipse_ptr) return false;
    return (((focus.first == ellipse_ptr->focus.first && focus.second == ellipse_ptr->focus.second) ||
             (focus.first == ellipse_ptr->focus.second && focus.second == ellipse_ptr->focus.first)) &&
            a == ellipse_ptr->a);
  }

  bool operator!=(const Shape& another) const override {
    return !(*this == another);
  }

  bool isCongruentTo(const Shape& another) const override {
    auto ellipse_ptr = dynamic_cast<const Ellipse*>(&another);
    if (!ellipse_ptr) return false;
    return points_distance(focus.first, focus.second) ==
           points_distance(ellipse_ptr->focus.first, ellipse_ptr->focus.first) &&
           a == ellipse_ptr->a;
  }

  bool isSimilarTo(const Shape& another) const override {
    auto ellipse_ptr = dynamic_cast<const Ellipse*>(&another);
    if (!ellipse_ptr) return false;
    return eccentricity() == ellipse_ptr->eccentricity();
  }

  bool containsPoint(Point point) const override {
    double s = points_distance(focus.first, point) + points_distance(focus.second, point);
    return s < a * 2 || s == 2 * a;
  }

  void rotate(Point center, double angle) override {
    angle *= NUM_PI / 180;
    focus.first = focus.first.rotate(center, angle);
    focus.second = focus.second.rotate(center, angle);
  }

  void reflect(Point center) override {
    scale(center, -1);
  }

  void reflect(Line axis) override {
    focus.first.reflect(axis);
    focus.second.reflect(axis);
  }

  void scale(Point center, double k) override {
    focus.first.scale(center, k);
    focus.second.scale(center, k);
    a *= std::abs(k);
    b *= std::abs(k);
  }
};

class Circle : public Ellipse {
public:
  Circle() : Ellipse() {};

  Circle(const Point& center, double radius) : Ellipse(center, center, 2 * radius) {}

  double radius() const { return a; }
};

class Rectangle : public Polygon {
public:
  Rectangle() : Polygon() {};

  Rectangle(const Point& point1, const Point& point2, double k) : Polygon{point1, Point(0, 0), point2, Point(0, 0)} {
    double angle = 180 - 2 * atan(k) * 180 / NUM_PI;
    Point middle = Point((point1.x + point2.x) / 2, (point1.y + point2.y) / 2);
    vertices_[1] = point1.rotate(middle, angle);
    vertices_[3] = point2.rotate(middle, angle);
  }

  Point center() const {
    return middle_point(vertices_[0], vertices_[2]);
  }

  std::pair<Line, Line> diagonals() const {
    return std::pair(Line(vertices_[0], vertices_[2]), Line(vertices_[1], vertices_[3]));
  }
};

class Square : public Rectangle {
public:
  Square() : Rectangle() {}

  Square(const Point& point1, const Point& point2) : Rectangle(point1, point2, 1) {}

  Circle circumscribedCircle() const {
    return Circle(middle_point(vertices_[0], vertices_[2]), points_distance(vertices_[0], vertices_[2]) / 2);
  }

  Circle inscribedCircle() const {
    return Circle(middle_point(vertices_[0], vertices_[2]), points_distance(vertices_[0], vertices_[1]) / 2);
  }
};

class Triangle : public Polygon {
public:
  Triangle() : Polygon() {}

  Triangle(const Point& point1, const Point& point2, const Point& point3) : Polygon{point1, point2, point3} {}

  Circle circumscribedCircle() const {
    Line mp1 = middle_perpendicular(vertices_[0], vertices_[1]);
    Line mp2 = middle_perpendicular(vertices_[1], vertices_[2]);
    return Circle(cross(mp1, mp2), points_distance(cross(mp1, mp2), vertices_[0]));
  }

  Circle inscribedCircle() const {
    double a = points_distance(vertices_[1], vertices_[2]);
    double b = points_distance(vertices_[2], vertices_[0]);
    double c = points_distance(vertices_[0], vertices_[1]);
    Point tmp;
    tmp.x = (vertices_[0].x * a + vertices_[1].x * b + vertices_[2].x * c) / (a + b + c);
    tmp.y = (vertices_[0].y * a + vertices_[1].y * b + vertices_[2].y * c) / (a + b + c);
    return Circle(tmp, 2 * area() / perimeter());
  }

  Point centroid() const {
    double x = (vertices_[0].x + vertices_[1].x + vertices_[2].x) / 3;
    double y = (vertices_[0].y + vertices_[1].y + vertices_[2].y) / 3;
    return Point(x, y);
  }

  Point orthocenter() const {
    Line ha = perpendicular(Line(vertices_[0], vertices_[1]));
    ha.c = -ha.a * vertices_[2].x - ha.b * vertices_[2].y;
    Line hb = perpendicular(Line(vertices_[1], vertices_[2]));
    hb.c = -hb.a * vertices_[0].x - hb.b * vertices_[0].y;
    return cross(ha, hb);
  }

  Line EulerLine() const {
    return Line(circumscribedCircle().center(), orthocenter());
  }

  Circle ninePointsCircle() const {
    Point center = middle_point(circumscribedCircle().center(), orthocenter());
    return Circle(center, circumscribedCircle().radius() / 2);
  }
};


