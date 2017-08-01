#include "Constants.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
using namespace std;

struct knotpoint
{
    double xcoord;   //position vector x coord
    double ycoord;   //position vector y coord
    double zcoord;   //position vector z coord
    double tx;       //grad vector x coord
    double ty;       //grad vector y coord
    double tz;       //grad vector z coord
    double length;   //length of line
};

struct knotcurve
{
    std::vector<knotpoint> knotcurve; // the actual data of the curve
    double length;   //total lengthh of line
};

struct griddata 
{
	double Nx,Ny,Nz;	
};

/*************************General maths and integer functions*****************************/

// little inline guys
inline double x(int i, const griddata& griddata);
inline double y(int i,const griddata& griddata);
inline double z(int i,const griddata& griddata);
inline int sign(int i);
inline  int pt( int i,  int j,  int k,const griddata& griddata);       //convert i,j,k to single index
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
inline int gridinc(int i, int p, int N, int direction );    //increment with reflecting boundary between -1 and 0 and N-1 and N

inline void Cross(float* a, float* b, float* result);
inline float dot(float* a, float* b);
inline void Norm(float* a);

// spherical geometry functions
inline float GeodesicDistance(float* a, float* b);
inline bool GeodesicIntersection(float* a1,float* a2,float* b1,float* b2);
/*************************Functions for knot initialisation*****************************/

void scalefunction(double *scale, double *midpoint, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin);

/*************************Functions for B and Phi calcs*****************************/

void DualConePhiCalc(vector<double>&phi, const griddata& griddata);
void print_phi( vector<double>&phi,const griddata& griddata);

