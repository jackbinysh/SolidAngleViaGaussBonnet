//  FN_knot_code.h
//
//
//  Created by Carl Whitfield on 17/05/2016.
//
//  Last modified 3/11/16
#include "FN_Constants.h"
#include "TriCubicInterpolator.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
using namespace std;

const double sixth = 1.0/6.0;
const double ONETHIRD = 1.0/3.0;
const double oneoverepsilon = 1.0/epsilon;
const double oneoverhsq = 1.0/(h*h);

struct griddata 
{
	double Nx,Ny,Nz;	
};
struct parameters
{
	gsl_vector *v,*f,*b;
    likely::TriCubicInterpolator* ucvmag;
    griddata mygriddata;
};

struct triangle
{
    double xvertex[3];   //x components of vertices
    double yvertex[3];   //y components of vertices
    double zvertex[3];   //z components of vertices
    double normal[3];   //normal vector
    double area;      //surface area
    double centre[3];  //centre position vector
};

struct knotpoint
{
    double xcoord;   //position vector x coord
    double ycoord;   //position vector y coord
    double zcoord;   //position vector z coord
    double modxcoord;   //position vector x coord, modded out by the lattice
    double modycoord;   //position vector y coord, ""
    double modzcoord;   //position vector z coord, ""
    double ax;       //grad vector x coord
    double ay;       //grad vector y coord
    double az;       //grad vector z coord
    double tx;       //grad vector x coord
    double ty;       //grad vector y coord
    double tz;       //grad vector z coord
    double nx;       //grad vector x coord
    double ny;       //grad vector y coord
    double nz;       //grad vector z coord
    double bx;       //grad vector x coord
    double by;       //grad vector y coord
    double bz;       //grad vector z coord
    double vx;       //grad vector x coord
    double vy;       //grad vector y coord
    double vz;       //grad vector z coord
    double vdotnx;       //grad vector x coord
    double vdotny;       //grad vector x coord
    double vdotnz;       //grad vector y coord
    double vdotbx;       //grad vector x coord
    double vdotby;       //grad vector x coord
    double vdotbz;       //grad vector y coord
    double curvature;        // curvature
    double torsion;        // torsion
    double twist;    //local twist value
    double writhe;   //local writhe value
    double length;   //length of line
};

struct knotcurve
{
    std::vector<knotpoint> knotcurve; // the actual data of the curve
    // global data for the knot component
    double twist;    //total twist value
    double writhe;   //total  writhe value
    double length;   //total lengthh of line
    double xavgpos;  // average position of the knot
    double yavgpos;
    double zavgpos;
};
/*************************General maths and integer functions*****************************/

// little inline guys
inline double x(int i, const griddata& griddata);
inline double y(int i,const griddata& griddata);
inline double z(int i,const griddata& griddata);
inline int sign(int i);
inline  int pt( int i,  int j,  int k,const griddata& griddata);       //convert i,j,k to single index
inline int circularmod(int i, int N);    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
inline int incw(int i, int p, int N);    //increment with reflecting boundary between -1 and 0 and N-1 and N
inline int gridinc(int i, int p, int N, int direction );    //increment with reflecting boundary between -1 and 0 and N-1 and N

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product);
inline void Cross(float* a, float* b, float* result);
inline float dot(float* a, float* b);
inline void Norm(float* a);
inline double my_f(const gsl_vector* minimum, void* params);
void rotatedisplace(double& xcoord, double& ycoord, double& zcoord, const double theta, const double dispx,const double dispy,const double dispz);

// spherical geometry functions
inline float GeodesicDistance(float* a, float* b);
inline bool GeodesicIntersection(float* a1,float* a2,float* b1,float* b2);
/*************************Functions for knot initialisation*****************************/

double initialise_knot(std::vector<triangle>& knotsurface);

double init_from_surface_file(std::vector<triangle>& knotsurface);

void scalefunction(double *scale, double *midpoint, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin);

/*************************Functions for B and Phi calcs*****************************/

void phi_calc( vector<double>&phi,std::vector<triangle>& knotsurface,const griddata& griddata);
void DualConePhiCalc(vector<double>&phi, const griddata& griddata);
void phi_calc_manual( vector<double>&phi,const griddata& griddata);

//FitzHugh Nagumo functions
void uv_initialise(vector<double>&phi, vector<double>&u, vector<double>&v,const griddata& griddata);
void crossgrad_calc( vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>&ucvmag,const griddata& griddata);
void find_knot_properties( vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>& ucvmag,vector<double>&u,vector<knotcurve>& knotcurves,double t, gsl_multimin_fminimizer* minimizerstate, const griddata& griddata);
void find_knot_velocity(const vector<knotcurve>& knotcurves,vector<knotcurve>& knotcurvesold,const griddata& griddata,const double deltatime);
void uv_update(vector<double>&u, vector<double>&v,  vector<double>&ku, vector<double>&kv,const griddata& griddata);
// 3d geometry functions
int intersect3D_SegmentPlane( knotpoint SegmentStart, knotpoint SegmentEnd, knotpoint PlaneSegmentStart, knotpoint PlaneSegmentEnd, double& IntersectionFraction, std::vector<double>& IntersectionPoint );
void resizebox(vector<double>&u,vector<double>&v,vector<double>&ucvx,vector<double>&ucvy,vector<double>&ucvz,vector<knotcurve>&knotcurves,vector<double>&ku,vector<double>&kv,griddata& oldgriddata);
void overlayknots(vector<knotcurve>& knotcurves,const vector<knotcurve>& knotcurvesold,const griddata& griddata);

/*************************File reading and writing*****************************/



void print_marked( vector<int>&marked,int shelllabel, const griddata& griddata);

void print_B_phi( vector<double>&phi,const griddata& griddata);
void print_uv( vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz,vector<double>&ucvmag, double t, const griddata& griddata);
int phi_file_read(vector<double>&phi,const griddata& griddata);
void print_knot( double t, vector<knotcurve>& knotcurves,const griddata& griddata,string seriesname="");
int uvfile_read(vector<double>&u, vector<double>&v, vector<double>& ku, vector<double>& kv, vector<double>& ucvx, vector<double>& ucvy,vector<double>& ucvz,griddata& griddata);
int uvfile_read_ASCII(vector<double>&u, vector<double>&v,const griddata& griddata); // for legacy purposes
int uvfile_read_BINARY(vector<double>&u, vector<double>&v,const griddata& griddata);
float FloatSwap( float f );
void ByteSwap(const char* TobeSwapped, char* swapped );


// things for the grown function

inline int incabsorb(int i, int p, int N);
void growshell(vector<double>&u,vector<int>& marked,double ucrit, const griddata& griddata);
void grow(const vector<double>&u,vector<int>&marked,double ucrit,const griddata& griddata);



