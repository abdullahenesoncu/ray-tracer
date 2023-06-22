#include <iostream>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include "maths.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;
namespace fs = std::filesystem;

int maxNumberOfSpheres = 15;
double minZ = 200;
double maxZ = 1000;
double rMin = 50;
double rMax = 100;
int width = 1000;
int height = 1000;
bool video = false;

double randomDouble( double rangeStart, double rangeEnd ){
	if ( rangeStart == rangeEnd ) return rangeStart;
	double ratio = abs(rand()*rand()%1000000000/1000000000.);
	return rangeStart + ( rangeEnd - rangeStart ) * ratio;
}

class Args {
public:
	int scenarioNumber = 0;
	int nSpheres = 1;
	string filename;
	bool video = false;
	bool hasAmbientTerms = false;
	bool hasDiffuseTerms = false;
	bool hasSpecularTerms = false;
	bool hasReflectionTerms = false;

	Args() {}
	Args( int argc, char *argv[] ){
		vector< string > v;
		for ( int i=1; i<argc; i++ ){
			v.push_back( string( argv[i] ) );
		}
		for ( int i=0; i<v.size(); i++ ) {
			if( v[i] == "-scenario" ){
				scenarioNumber = atoi( v[ i + 1 ].c_str() );
				i++;
			}
			else if( v[i] == "-n" ){
				nSpheres = atoi( v[ i + 1 ].c_str() );
				i++;
			}
			else if( v[i] == "-f" ){
				filename = v[ i + 1 ];
				i++;
			}
			else if ( v[i] == "-v" )	video = true;
			else if ( v[i] == "-a" )	hasAmbientTerms = true;
			else if ( v[i] == "-d" )	hasDiffuseTerms = true;
			else if ( v[i] == "-s" )	hasSpecularTerms = true;
			else if ( v[i] == "-r" )	hasReflectionTerms = true;
		}
	}
};

class Scenario {
public:
	
	double defaultKa = 500;
	double defaultKd = 500;
	double defaultKs = 500;
	double defaultKr = 0.5;

	vector< Object* > objects;
	vector< Sphere* > spheres;
	Args args;
	int seed;

	Scenario() {}
	Scenario( Args args_, int seed_, vector< Sphere* > spheres_ ) {
		args = args_;
		seed = seed_;
		
		for ( auto sphere: spheres_ ) {
			spheres.push_back( sphere );
		}
		
		if ( args.filename.size() ) {
			ifstream in( args.filename );
			int n, _;
			in >> _ >> n;
			while ( n-- ) {
				double ka = args.hasAmbientTerms ? 500 : 0;
				double kd = args.hasDiffuseTerms ? 500 : 0;
				double ks = args.hasSpecularTerms ? 500 : 0;
				double kr = args.hasReflectionTerms ? 0.2 : 0;
				double x, y, z, radius;
				int r, g, b;
				in >> x >> y >> z >> radius;
				in >> r >> g >> b;;
				Sphere *s = new Sphere(
					Point( x, y, z ),
					radius,
					Color( r, g, b, ka, kd, ks, kr )
				);
				spheres.push_back( s );
			}
		}

		double maxX = 50;
		for ( auto sphere: spheres ) {
			maxX = max( maxX, sphere->getMaxX() );
			objects.push_back( sphere );
		}

		addGround( maxX );
	}

	void addGround( double x ) {
		for ( int i = -1000; i < 1000; i += 100 ) {
			for ( int j = -1000; j < 1000; j += 100 ) {
				unsigned char c = ( i + j ) / 100 % 2 * 255;
				Color color = Color( c, c, c, 
									args.hasAmbientTerms ? defaultKa : 0,
									args.hasDiffuseTerms ? defaultKd : 0,
									args.hasSpecularTerms ? defaultKs : 0,
									args.hasReflectionTerms ? defaultKr : 0 );
				Quadrangle *q = new Quadrangle(
					Point( x, i, j ),
					Point( x, i, j+100 ),
					Point( x, i+100, j+100 ),
					Point( x, i+100, j ),
					color
				);
				objects.push_back( q );
			}
		}
	}

	string dump() {
		string seedStr = to_string( seed );
		ofstream out( "data/"+seedStr+".txt" );
		out << seed << endl << spheres.size() << endl;
		for ( auto sphere: spheres ) {
			out << sphere->center.x << " " << sphere->center.y << " " << sphere->center.z << " " << sphere->radius << endl;
			Color color = sphere->color;
			out << int(color.r) << " " << int(color.g) << " " << int(color.b) << endl;
		}
		return seedStr;
	}

};

class Runner{
public:
	Scenario scenario;
	Point lightPosition = Point( -500, 500, 500 );
	Point eyePosition = Point( 0, 0, 0 );
	vector< Object* > objects;
	double lightLuminance = 100000;

	Runner() {}
	Runner( Scenario scenario_ ) {
		scenario = scenario_;
		objects = scenario.objects;
	}

	pair< Object*, Point* > getCollidedObject( Ray ray ) {
		Point * intersectionPoint = NULL;
		Object * collidedObject;
		for ( auto object: objects ) {
			Point * p = object->intersectsWithRayOrNull( ray );
			if ( p && ( intersectionPoint == NULL || 
							p->distanceTo( ray.origin ) < intersectionPoint->distanceTo( ray.origin ) ) ){
				intersectionPoint = p;
				collidedObject = object;
			}
		}
		return { collidedObject, intersectionPoint };
	}

	Color* getColorOnRay( Ray ray, double reflectionMultiplier=1.0 ) {
		if ( reflectionMultiplier < 0.01 ) return NULL;
		
		auto [ collidedObject, intersectionPoint ] = getCollidedObject( ray );
		if ( !intersectionPoint )	return NULL;
		
		Color* resColor = new Color(collidedObject->color);

		Ray lightRay = Ray( lightPosition, *intersectionPoint - lightPosition );
		auto [ _, lightRayIntersectionPoint ] = getCollidedObject( lightRay );

		double totalI = 0;
		if ( lightRayIntersectionPoint->distanceTo( lightPosition ) + 1e-5 < intersectionPoint->distanceTo( lightPosition ) ) {
		}
		else {
			double d = intersectionPoint->distanceTo( lightPosition );
			// Handle ambient term
			totalI += resColor->ka * lightLuminance / ( d * d );
			// Handle diffuse term
			double theta = ( lightPosition - *lightRayIntersectionPoint ).angleBetween( collidedObject->getNormal( *lightRayIntersectionPoint ) );
			totalI += resColor->kd * lightLuminance * max( 0.0, cos( theta ) ) / ( d * d );
			// Handle specular term
			Point bis = ( ( lightPosition - *lightRayIntersectionPoint ) + ( ray.origin - *lightRayIntersectionPoint ) ).normalized();
			double alpha = bis.angleBetween( collidedObject->getNormal( *lightRayIntersectionPoint ) );
			totalI += resColor->ks * lightLuminance * max( 0.0, pow( cos( alpha ), 2 ) ) / ( d * d );
		}
		resColor = new Color( resColor->getColorUnderLight( totalI ) );
		Point r = ( ray.origin - *intersectionPoint ).normalized();
		Point normal = collidedObject->getNormal( *intersectionPoint );
		Point reflectionDir = ( 2 * normal - r );
		double kr = collidedObject->color.kr;
		Color* reflectedColor = getColorOnRay( Ray( *intersectionPoint, reflectionDir ), reflectionMultiplier * kr );
		if ( reflectedColor ) {
			resColor->b = resColor->b * ( 1.0 - kr ) + reflectedColor->b * kr;
			resColor->g = resColor->g * ( 1.0 - kr ) + reflectedColor->g * kr;
			resColor->r = resColor->r * ( 1.0 - kr ) + reflectedColor->r * kr;
		}
	
		return resColor;
	}

	vector< unsigned char > render( Point eyePosition, Point lightPosition ) {
		vector< unsigned char > pixels;
		for ( int row = 0; row < height; row ++ ) {
			for ( int col = 0; col < width; col ++ ) {
				double x = ( row + 0.5 ) / height * 100.0 - 50.0;
				double y = ( col + 0.5 ) / height * 100.0 - 50.0;
				auto pixelPoint = Point( x, y, 100.0 );
				auto ray = Ray( eyePosition, pixelPoint - eyePosition );
				Color* color = getColorOnRay( ray );
				if ( !color )	color = new Color( 0, 0, 0 );
				pixels.push_back( color->r );
				pixels.push_back( color->g );
				pixels.push_back( color->b );
			}
		}
		return pixels;
	}

	void run() {
		if ( !scenario.args.video ) {
			string seedStr = scenario.dump();
			auto pixels = render( eyePosition, lightPosition );
			stbi_write_png( ("data/"+seedStr+".png").c_str(), width, height, 3, &(pixels[0]), 3*width);
		}
		else {
			string seedStr = scenario.dump();
			fs::create_directories( "video/"+seedStr );
			for ( auto theta=200; theta>=0; theta-- ) {
				lightPosition = Point( 500*cos(theta/180.0*acos(-1)), 0, 500*sin(theta/180.0*acos(-1)) );
				auto pixels = render( eyePosition, lightPosition );
				string thetaStr = to_string( theta );
				stbi_write_png( ("video/"+seedStr+"/"+thetaStr+".png").c_str(), width, height, 3, &(pixels[0]), 3*width);
			}
		}
	}
};

Scenario getScenario1( Args args, int seed ) {
	double ka = args.hasAmbientTerms ? 500 : 0;
	double kd = args.hasDiffuseTerms ? 500 : 0;
	double ks = args.hasSpecularTerms ? 500 : 0;
	double kr = args.hasReflectionTerms ? 0.2 : 0;
	vector< Sphere* > spheres = {
		new Sphere( Point( 0, 0, 400 ), 50, Color( 0, 0, 255, ka, kd, ks, kr ) ),
		new Sphere( Point( 0, -100, 400 ), 50, Color( 255, 0, 0, ka, kd, ks, kr ) ),
	};
	
	return Scenario( args, seed, spheres );
}

Scenario getRandomScenario( Args args, int seed ) {
	vector< Sphere* > spheres;
	for ( int i = 0; i < args.nSpheres; i ++ ) {
		double ka = args.hasAmbientTerms ? 500 : 0;
		double kd = args.hasDiffuseTerms ? 500 : 0;
		double ks = args.hasSpecularTerms ? 500 : 0;
		double kr = args.hasReflectionTerms ? 0.2 : 0;
		Color color = Color( rand() % 255, rand() % 255, rand() % 255, ka, kd, ks, kr );

		do{
			double z = randomDouble( minZ, maxZ );
			double x = randomDouble( - z / 2, + z / 2 );
			double y = randomDouble( - z / 2, + z / 2 );
			double r = 50;
			Sphere* sphere = new Sphere( Point( x, y, z ), r, color );
			bool ok = true;
			for ( auto sphere2: spheres ) {
				if ( sphere2->intersectsSphere( *sphere ) ){
					ok = false;
					break;
				}
			}
			if ( ok ) {
				spheres.push_back( sphere );
				break;
			}
		} while ( true );
	}
	return Scenario( args, seed, spheres );
}

int main( int argc, char *argv[] ) {
	int seed = time( 0 );
	srand( seed );
	
	Args args = Args( argc, argv );
	Scenario scenario;
	
	if ( args.scenarioNumber == 1 ) {
		scenario = getScenario1( args, 1 );
	}
	else if ( args.filename.size() ) {
		scenario = Scenario( args, seed, {} );
	}
	else {
		scenario = getRandomScenario( args, seed );
	}
	
	Runner runner( scenario );
	runner.run();
}
