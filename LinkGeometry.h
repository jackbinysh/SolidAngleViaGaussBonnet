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
  double kappaBx;  //integrand we need, x component
  double kappaBy;  //integrand we need, y component
  double kappaBz;  //integrand we need, z component
};

struct knotcurve
{
  std::vector<knotpoint> knotcurve; //the actual data of the curve
  double length;   //total length of line
};

struct viewpoint
{
  double xcoord;
  double ycoord;
  double zcoord;
};

/*************************Functions for knot initialisation*****************************/

void InitialiseLink(vector<struct knotcurve>& Link);

/**********************Functions for curve geometry************************/

void ComputeLengths(struct knotcurve& Curve);
void ComputeTangent(struct knotcurve& Curve);
void ComputeKappaB(struct knotcurve& Curve);
double Writhe(const struct knotcurve& Curve);
double Twist(const struct knotcurve& Curve, const struct viewpoint& View);

/***********************Functions for outputting the solid angle*************************/

double SolidAngleCalc(const struct knotcurve& Curve, const struct viewpoint& View);
double SolidAngleCalc2(const struct knotcurve& Curve, const struct viewpoint& View);
void OutputSolidAngle(const vector<struct knotcurve>& Link, const double *writhe);

/*************************General maths and integer functions*****************************/

double CheckThreshold(const struct knotcurve& Curve, const struct viewpoint& View);
// little inline guys
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
