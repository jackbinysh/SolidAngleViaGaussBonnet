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
    double tx;       //tangent vector x component
    double ty;       //tangent vector y component
    double tz;       //tangent vector z component
    double length;   //length of line
    double kappaNx;  //curvature vector x component
    double kappaNy;  //curvature vector x component
    double kappaNz;  //curvature vector x component
    double kappaBx;  //integrand we need, x component
    double kappaBy;  //integrand we need, y component
    double kappaBz;  //integrand we need, z component
};

struct knotcurve
{
    std::vector<knotpoint> knotcurve; //the actual data of the curve
    double length;   //total length of line
    double writhe;
    double twist;
};

struct Link
{
    std::vector<knotcurve> Components;
    int NumComponents;
    int NumPoints;
    double length;   //total length of line
    double writhe;
    double twist;
}; 

struct viewpoint
{
    double xcoord;
    double ycoord;
    double zcoord;
};

struct griddata 
{
    int Nx;
    int Ny;
    int Nz;
};

/*************************Functions for knot initialisation*****************************/

void InitialiseFromFile(struct Link& Curve);
void ScaleFunction(double *scale, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin);
void InitialiseCurve(struct Link& Curve);
void RefineCurve(struct Link& Curve);

/**********************Functions for curve geometry************************/

void ComputeLengths(struct Link& Curve);
void ComputeTangent(struct Link& Curve);
void ComputeKappaB(struct Link& Curve);
void Writhe(struct Link& Curve);
void Twist(struct Link& Curve, const struct viewpoint& View);

/***********************Functions for outputting the solid angle*************************/

double SolidAngleCalc(const struct Link& Curve, const struct viewpoint& View);
double SolidAngleCalc2(const struct Link& Curve, const struct viewpoint& View);
void OutputSolidAngle(const struct Link& Curve);

/*************************General maths and integer functions*****************************/

void CheckThreshold(const struct Link& Curve, const struct viewpoint& Point,double& ndotTmin,double& ndotTmax);
// little inline guys
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
inline double x(int i,const griddata& griddata);
inline double y(int i,const griddata& griddata);
inline double z(int i,const griddata& griddata);
