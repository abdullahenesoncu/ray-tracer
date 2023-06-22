#include <cmath>

using namespace std;

class Point {
public:
    double x, y, z;
    Point() {}
    Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    double distanceTo(Point other) {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return sqrt( dx * dx  + dy * dy + dz * dz );
    }

    double magnitude() {
        return sqrt( x * x + y * y + z * z );
    }

    Point normalized() {
        double l = magnitude();
        return Point( x / l, y / l, z / l );
    }

    double angleBetween( Point p ) {
        double dotProduct = x * p.x + y * p.y + z * p.z;
        double cosine = dotProduct / ( magnitude() * p.magnitude() );
        return acos(cosine);
    }

    Point cross( Point other ) {
        double x = y * other.z - z * other.y;
        double y = z * other.x - x * other.z;
        double z = x * other.y - y * other.x;
        return Point( x, y, z );
    }

    double dot( Point other ) {
        return x * other.x + y * other.y + z * other.z;
    }
};

Point operator-( Point a, Point b ) {
    return Point( a.x - b.x, a.y - b.y, a.z - b.z );
}

Point operator+( Point a, Point b ) {
    return Point( a.x + b.x, a.y + b.y, a.z + b.z );
}

Point operator*( int a, Point b ) {
    return Point( a * b.x, a * b.y, a * b.z );
}

class Ray {
public:
    Point origin, direction;
    Ray() {}
    Ray(Point origin_, Point direction_) : origin(origin_), direction(direction_) {
        double length = direction.magnitude();
        direction.x /= length;
        direction.y /= length;
        direction.z /= length;
    }
    
    Point getPointAt(double t) {
        double x = origin.x + t * direction.x;
        double y = origin.y + t * direction.y;
        double z = origin.z + t * direction.z;
        return Point(x, y, z);
    }

    double angleBetween( Ray r ) {
        return direction.angleBetween( r.direction );
    }

};

class Color {
public:
    unsigned char b, g, r;
    double ka, kd, ks, kr;
    Color() {}
    Color( unsigned char b_, unsigned char g_, unsigned char r_, double ka_, double kd_, double ks_, double kr_ ) 
        : b( b_ ), g( g_ ), r( r_ ), ka( ka_ ), kd( kd_ ), ks( ks_ ), kr( kr_ ) {}
    Color( unsigned char b_, unsigned char g_, unsigned char r_ ) 
        : b( b_ ), g( g_ ), r( r_ ) {}
    
    Color getColorUnderLight( double lightIntensity ) {
        double multiplier = 1 / ( 1 + exp( - lightIntensity / 100 + 2 ) );
        return Color( b * multiplier, g * multiplier, r * multiplier, ka, kd, ks, kr );
    }
};

class Object {
public:
    Color color;
    Object() {}
    Object( Color color_ ) : color( color_ ) {}
    virtual Point * intersectsWithRayOrNull( Ray ray ) = 0;
    virtual Point getNormal(Point p) = 0;
};

class Sphere : public Object {
public:
    Point center;
    double radius;
    Sphere() {}
    Sphere( Point center_, double radius_, Color color_ ) : center( center_ ), radius( radius_ ), Object( color_ ) {}

    Point * intersectsWithRayOrNull( Ray ray ) {
        Point origin = ray.origin;
        Point dir = ray.direction;
        double a = dir.magnitude() * dir.magnitude();

        double b = 2 * ( dir.x * ( origin.x - center.x ) + dir.y * ( origin.y - center.y ) + dir.z * ( origin.z - center.z ) );

        double c = center.magnitude() * center.magnitude() + origin.magnitude() * origin.magnitude()
                    - 2 * ( center.x * origin.x + center.y * origin.y + center.z * origin.z ) - radius * radius;

        double disc = b * b - 4 * a * c;

        if ( disc < 0 ) {
            return NULL;
        }

        double t1 = ( - b - sqrt( disc ) ) / ( 2 * a );
        double t2 = ( - b + sqrt( disc ) ) / ( 2 * a );

        if ( t1 < 1e-5 && t2 < 1e-5 ) {
            return nullptr;
        }

        double t;
        if ( t1 < 1e-5 )   t = t2;
        else if ( t2 < 1e-5 )  t = t1;
        else    t = min( t1, t2 );

        return new Point( ray.getPointAt( t ) );
    }

    Point getNormal( Point p ) {
        return ( p - center ).normalized();
    }

    bool intersectsSphere( Sphere other ) {
        return radius + other.radius >= center.distanceTo( other.center );
    }

    double getMaxX() {
        return center.x + radius;
    }
};

class Quadrangle : public Object {
public:
    Point p1, p2, p3, p4;
    Quadrangle() {}
    Quadrangle( Point p1_, Point p2_, Point p3_, Point p4_, Color color_ ) : p1( p1_ ), p2( p2_ ), p3( p3_ ), p4( p4_ ), Object( color_ ) {}

    Point * intersectsWithRayOrNull( Ray ray ) {
        Point p = ray.origin;
        Point d = ray.direction;

        Point n = getNormal( Point() );

        double t = -( n.x * p.x + n.y * p.y + n.z * p.z - n.x * p1.x - n.y * p1.y - n.z * p1.z ) 
                    / ( n.x * d.x + n.y * d.y + n.z * d.z );

        if ( t < 1e-5 ) {
            return NULL;
        }

        Point intersection = Point( p.x + t * d.x, p.y + t * d.y, p.z + t * d.z );

        if ( isInside( intersection ) ) {
            return new Point( intersection );
        }

        return NULL;
    }

    Point getNormal( Point p ) {
        Point edge1 = Point( p2.x - p1.x, p2.y - p1.y, p2.z - p1.z );
        Point edge2 = Point( p3.x - p1.x, p3.y - p1.y, p3.z - p1.z );
        Point normal = Point( edge1.y * edge2.z - edge1.z * edge2.y, edge1.z * edge2.x - edge1.x * edge2.z, edge1.x * edge2.y - edge1.y * edge2.x );

        return normal.normalized();
    }

    bool isInside( Point p ) {
        Point v1 = p1 - p;
        Point v2 = p2 - p;
        Point v3 = p3 - p;
        Point v4 = p4 - p;

        // Calculate the cross products of adjacent vectors
        Point n1 = v1.cross(v2);
        Point n2 = v2.cross(v3);
        Point n3 = v3.cross(v4);
        Point n4 = v4.cross(v1);

        Point refNormal = getNormal( Point() );

        // Check if the dot product of the normal vectors and the reference normal
        // are all positive or all negative, indicating that the point is inside the quadrangle
        if (n1.dot(refNormal) > 0 && n2.dot(refNormal) > 0 &&
            n3.dot(refNormal) > 0 && n4.dot(refNormal) > 0) {
            return true;
        }
        if (n1.dot(refNormal) < 0 && n2.dot(refNormal) < 0 &&
            n3.dot(refNormal) < 0 && n4.dot(refNormal) < 0) {
            return true;
        }

        return false;
    }
};
