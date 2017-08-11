#include "LinkGeometry.h"    
#include <math.h>
#include <string.h>

int main (void)
{
  knotcurve Curve1,Curve2;
  InitialiseCurve(Curve1,Curve2);
  // compute the writhe of each component
  double writhe1 = Writhe(Curve1);
  double writhe2 = Writhe(Curve2);

  // it would be good to construct a set of points to compute the twist from
  viewpoint View;
  View.xcoord = 0.0;
  View.ycoord = 5.0;
  View.zcoord = 2.02;

  // compute the twist of each component with the projective framing
  double twist1 = Twist(Curve1,View);
  double twist2 = Twist(Curve2,View);
  // compute the solid angle from the view point
  double omega1 = SolidAngleCalc(Curve1,Curve2,writhe1,writhe2,View);
  double omega2 = SolidAngleCalc2(Curve1,Curve2,writhe1,writhe2,View);

  cout << "Wr1 = " << writhe1 << "  Tw1 = " << twist1 << "  SL1 = " << writhe1+twist1 << endl; 
  cout << "Wr2 = " << writhe2 << "  Tw2 = " << twist2 << "  SL2 = " << writhe2+twist2 << endl; 
  cout << "omega1 = " << omega1 << "  omega2 = " << omega2 << endl;

  // output the solid angle on a grid
  cout << "calculating the solid angle..." << endl;
  OutputSolidAngle(Curve1,Curve2,writhe1,writhe2);

  return 0;
}

// initialise the link -- torus links only
void InitialiseCurve(knotcurve& Curve1, knotcurve& Curve2)
{
  // torus knot -- ( (1 + r cos qu) cos pu , (1 + r cos qu) sin pu , r sin qu )
  int p = 1;
  int q = 2;
  double r = 0.57;
  int NP = 2000;
  for (int s=0; s<NP; s++)
    {
      knotpoint Point;
      Point.xcoord = (1.24+r*cos(2.0*M_PI*q*s/NP))*cos(2.0*M_PI*p*s/NP);
      Point.ycoord = (1.24+r*cos(2.0*M_PI*q*s/NP))*sin(2.0*M_PI*p*s/NP);
      Point.zcoord = r*sin(2.0*M_PI*q*s/NP);
      Curve1.knotcurve.push_back(Point);
      // second link component
      //      Point.xcoord = (1.24-r*cos(2.0*M_PI*q*s/NP))*cos(2.0*M_PI*p*s/NP);
      //      Point.ycoord = (1.24-r*cos(2.0*M_PI*q*s/NP))*sin(2.0*M_PI*p*s/NP);
      //      Point.zcoord = -r*sin(2.0*M_PI*q*s/NP);
      // second link component -- reverse orientation
      Point.xcoord = (1.24-r*cos(2.0*M_PI*q*s/NP))*cos(2.0*M_PI*p*s/NP);
      Point.ycoord = -(1.24-r*cos(2.0*M_PI*q*s/NP))*sin(2.0*M_PI*p*s/NP);
      Point.zcoord = r*sin(2.0*M_PI*q*s/NP);
      Curve2.knotcurve.push_back(Point);
    }

  // first component
  ComputeLengths(Curve1);
  ComputeTangent(Curve1);
  ComputeKappaB(Curve1);
  // second component
  ComputeLengths(Curve2);
  ComputeTangent(Curve2);
  ComputeKappaB(Curve2);
}

// computes the lengths of each curve segment
void ComputeLengths(knotcurve& Curve)
{
  int NP = Curve.knotcurve.size();

  for(int s=0; s<NP; s++)        
    {
      double dx = (Curve.knotcurve[incp(s,1,NP)].xcoord - Curve.knotcurve[s].xcoord);  
      double dy = (Curve.knotcurve[incp(s,1,NP)].ycoord - Curve.knotcurve[s].ycoord);
      double dz = (Curve.knotcurve[incp(s,1,NP)].zcoord - Curve.knotcurve[s].zcoord);
      double deltas = sqrt(dx*dx+dy*dy+dz*dz);
      Curve.knotcurve[s].length = deltas;
      Curve.length += deltas;
    }
}

// computes the tangent vector to the curve
void ComputeTangent(knotcurve& Curve)
{
  int NP = Curve.knotcurve.size();

  for(int s=0; s<NP; s++)        // new scheme -- central difference
    {
      double dsp = Curve.knotcurve[s].length;
      double dsm = Curve.knotcurve[incp(s,-1,NP)].length;
      Curve.knotcurve[s].tx = (dsm/(dsp*(dsp+dsm)))*Curve.knotcurve[incp(s,1,NP)].xcoord + ((dsp-dsm)/(dsp*dsm))*Curve.knotcurve[s].xcoord - (dsp/(dsm*(dsp+dsm)))*Curve.knotcurve[incp(s,-1,NP)].xcoord;
      Curve.knotcurve[s].ty = (dsm/(dsp*(dsp+dsm)))*Curve.knotcurve[incp(s,1,NP)].ycoord + ((dsp-dsm)/(dsp*dsm))*Curve.knotcurve[s].ycoord - (dsp/(dsm*(dsp+dsm)))*Curve.knotcurve[incp(s,-1,NP)].ycoord;
      Curve.knotcurve[s].tz = (dsm/(dsp*(dsp+dsm)))*Curve.knotcurve[incp(s,1,NP)].zcoord + ((dsp-dsm)/(dsp*dsm))*Curve.knotcurve[s].zcoord - (dsp/(dsm*(dsp+dsm)))*Curve.knotcurve[incp(s,-1,NP)].zcoord;
    }
}

// computes the integrand we need
void ComputeKappaB(knotcurve& Curve)
{
  int NP = Curve.knotcurve.size();

  for(int s=0; s<NP; s++)    // new scheme -- central difference
    {
      double dsp = Curve.knotcurve[s].length;
      double dsm = Curve.knotcurve[incp(s,-1,NP)].length;
      double kappaNx = (dsm/(dsp*(dsp+dsm)))*Curve.knotcurve[incp(s,1,NP)].tx + ((dsp-dsm)/(dsp*dsm))*Curve.knotcurve[s].tx - (dsp/(dsm*(dsp+dsm)))*Curve.knotcurve[incp(s,-1,NP)].tx;
      double kappaNy = (dsm/(dsp*(dsp+dsm)))*Curve.knotcurve[incp(s,1,NP)].ty + ((dsp-dsm)/(dsp*dsm))*Curve.knotcurve[s].ty - (dsp/(dsm*(dsp+dsm)))*Curve.knotcurve[incp(s,-1,NP)].ty;
      double kappaNz = (dsm/(dsp*(dsp+dsm)))*Curve.knotcurve[incp(s,1,NP)].tz + ((dsp-dsm)/(dsp*dsm))*Curve.knotcurve[s].tz - (dsp/(dsm*(dsp+dsm)))*Curve.knotcurve[incp(s,-1,NP)].tz;
      Curve.knotcurve[s].kappaBx = Curve.knotcurve[s].ty*kappaNz - Curve.knotcurve[s].tz*kappaNy;
      Curve.knotcurve[s].kappaBy = Curve.knotcurve[s].tz*kappaNx - Curve.knotcurve[s].tx*kappaNz;
      Curve.knotcurve[s].kappaBz = Curve.knotcurve[s].tx*kappaNy - Curve.knotcurve[s].ty*kappaNx;
    }
}

// computes the writhe of a curve
double Writhe(const knotcurve& Curve)
{
  double Wr = 0.0;
  int NP = Curve.knotcurve.size();
  
  for (int s=0; s<NP; s++) // running over the knot
    {
      double a1x = Curve.knotcurve[s].xcoord;
      double a1y = Curve.knotcurve[s].ycoord;
      double a1z = Curve.knotcurve[s].zcoord;
      double t1x = Curve.knotcurve[s].tx;
      double t1y = Curve.knotcurve[s].ty;
      double t1z = Curve.knotcurve[s].tz;
      double ds = Curve.knotcurve[s].length;

      for (int t=s+1; t<NP; t++) // run over all points ahead of s
	{
	  double a2x = Curve.knotcurve[t].xcoord;
	  double a2y = Curve.knotcurve[t].ycoord;
	  double a2z = Curve.knotcurve[t].zcoord;
	  double t2x = Curve.knotcurve[t].tx;
	  double t2y = Curve.knotcurve[t].ty;
	  double t2z = Curve.knotcurve[t].tz;
	  double dt = Curve.knotcurve[t].length;

	  double dist = sqrt((a1x-a2x)*(a1x-a2x)+(a1y-a2y)*(a1y-a2y)+(a1z-a2z)*(a1z-a2z));
	  Wr += ds*dt*((a1x-a2x)*(t1y*t2z-t1z*t2y)+(a1y-a2y)*(t1z*t2x-t1x*t2z)+(a1z-a2z)*(t1x*t2y-t1y*t2x))/(dist*dist*dist);
	}
    }

  Wr /= 2.0*M_PI; 
  return Wr;
}

// computes the twist of a curve with projective framing from a view point
double Twist(const knotcurve& Curve, const viewpoint& View)
{
  double Angle = 0.0; 
  int NP = Curve.knotcurve.size();
  
  for (int s=0; s<NP; s++) // running over the knot
    {
      double curvex = Curve.knotcurve[s].xcoord;
      double curvey = Curve.knotcurve[s].ycoord;
      double curvez = Curve.knotcurve[s].zcoord;
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = curvex - View.xcoord;
      double viewy = curvey - View.ycoord;
      double viewz = curvez - View.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
      //      viewx /= dist;
      //      viewy /= dist;
      //      viewz /= dist;
      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve.knotcurve[s].tx;
      double ty = Curve.knotcurve[s].ty;
      double tz = Curve.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve.knotcurve[s].kappaBx;
      double kappaBy = Curve.knotcurve[s].kappaBy;
      double kappaBz = Curve.knotcurve[s].kappaBz;
      double ds = Curve.knotcurve[s].length;
      // and here's the integrand
      //      Angle += ds*ndotT*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(1.0 - ndotT*ndotT);
      Angle += ds*ndotT*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist*dist - ndotT*ndotT);
      // new version -- this approach is significantly worse!!
      /*      double dkappaNx = Curve.knotcurve[s].tx - Curve.knotcurve[incp(s,-1,NP)].tx;
      double dkappaNy = Curve.knotcurve[s].ty - Curve.knotcurve[incp(s,-1,NP)].ty;
      double dkappaNz = Curve.knotcurve[s].ty - Curve.knotcurve[incp(s,-1,NP)].tz;
      double dkappaBx = ty*dkappaNz - tz*dkappaNy;
      double dkappaBy = tz*dkappaNx - tx*dkappaNz;
      double dkappaBz = tx*dkappaNy - ty*dkappaNx;
      // and here's the integrand
      Angle += ndotT*(viewx*dkappaBx + viewy*dkappaBy + viewz*dkappaBz)/(1.0 - ndotT*ndotT);*/
    }

  Angle /= 2.0*M_PI;  
  return Angle;
}

double SolidAngleCalc(const knotcurve& Curve1, const knotcurve& Curve2, const double writhe1, const double writhe2, const viewpoint& View)
{
  double omega;
  int NP1 = Curve1.knotcurve.size();
  int NP2 = Curve2.knotcurve.size();

  double Angle = 0.0;
  // first link component
  for (int s=0; s<NP1; s++) 
    {
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = Curve1.knotcurve[s].xcoord - View.xcoord;
      double viewy = Curve1.knotcurve[s].ycoord - View.ycoord;
      double viewz = Curve1.knotcurve[s].zcoord - View.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve1.knotcurve[s].tx;
      double ty = Curve1.knotcurve[s].ty;
      double tz = Curve1.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve1.knotcurve[s].kappaBx;
      double kappaBy = Curve1.knotcurve[s].kappaBy;
      double kappaBz = Curve1.knotcurve[s].kappaBz;
      double ds = Curve1.knotcurve[s].length;
      // and here's the integrand
      Angle += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
    }
  // second link component
  for (int s=0; s<NP2; s++) 
    {
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = Curve2.knotcurve[s].xcoord - View.xcoord;
      double viewy = Curve2.knotcurve[s].ycoord - View.ycoord;
      double viewz = Curve2.knotcurve[s].zcoord - View.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve2.knotcurve[s].tx;
      double ty = Curve2.knotcurve[s].ty;
      double tz = Curve2.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve2.knotcurve[s].kappaBx;
      double kappaBy = Curve2.knotcurve[s].kappaBy;
      double kappaBz = Curve2.knotcurve[s].kappaBz;
      double ds = Curve2.knotcurve[s].length;
      // and here's the integrand
      Angle += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
    }
    
  //  omega = 2.0*M_PI*(1.0+writhe1+writhe2) - Angle;  // 1 replaced by |L|=2
  omega = 2.0*M_PI*(writhe1+writhe2) - Angle;  
  // put in the interval [-2pi,2pi] or [0,4pi]
  //  while(omega>2*M_PI) omega -= 4*M_PI;
  while(omega>4*M_PI) omega -= 4*M_PI;
  //  while(omega<-2*M_PI) omega += 4*M_PI;
  while(omega<0) omega += 4*M_PI;

  return omega;
}

double SolidAngleCalc2(const knotcurve& Curve1, const knotcurve& Curve2, const double writhe1, const double writhe2, const viewpoint& View)
{
  double omega;
  int NP1 = Curve1.knotcurve.size();
  int NP2 = Curve2.knotcurve.size();

  double Angle = 0.0;
  // first link component
  for (int s=0; s<NP1; s++) 
    {
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = Curve1.knotcurve[s].xcoord - View.xcoord;
      double viewy = Curve1.knotcurve[s].ycoord - View.ycoord;
      double viewz = Curve1.knotcurve[s].zcoord - View.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve1.knotcurve[s].tx;
      double ty = Curve1.knotcurve[s].ty;
      double tz = Curve1.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve1.knotcurve[s].kappaBx;
      double kappaBy = Curve1.knotcurve[s].kappaBy;
      double kappaBz = Curve1.knotcurve[s].kappaBz;
      double ds = Curve1.knotcurve[s].length;
      // and here's the integrand
      Angle += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist - ndotT);
    }
  // second link component
  for (int s=0; s<NP2; s++) 
    {
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = Curve2.knotcurve[s].xcoord - View.xcoord;
      double viewy = Curve2.knotcurve[s].ycoord - View.ycoord;
      double viewz = Curve2.knotcurve[s].zcoord - View.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve2.knotcurve[s].tx;
      double ty = Curve2.knotcurve[s].ty;
      double tz = Curve2.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve2.knotcurve[s].kappaBx;
      double kappaBy = Curve2.knotcurve[s].kappaBy;
      double kappaBz = Curve2.knotcurve[s].kappaBz;
      double ds = Curve2.knotcurve[s].length;
      // and here's the integrand
      Angle += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist - ndotT);
    }
    
  //  omega = 2.0*M_PI*(1.0-writhe1-writhe2) - Angle;  // 1 replaced by |L|=2
  omega = -2.0*M_PI*(writhe1+writhe2) - Angle;  
  // put in the interval [-2pi,2pi] or [0,4pi]
  //  while(omega>2*M_PI) omega -= 4*M_PI;
  while(omega>4*M_PI) omega -= 4*M_PI;
  //  while(omega<-2*M_PI) omega += 4*M_PI;
  while(omega<0) omega += 4*M_PI;

  return omega;
}

void OutputSolidAngle(const knotcurve& Curve1, const knotcurve& Curve2, const double writhe1, const double writhe2)
{
  double SolidAngle;
  viewpoint Point;
  int NP1 = Curve1.knotcurve.size();
  int NP2 = Curve2.knotcurve.size();

  // run over a grid of points of size Nx*Ny*Nz -- actual size Lx*Ly*Lz
  int Nx,Ny,Nz;
  Nx = 101; Ny = 101; Nz = 101;
  double Lx,Ly,Lz;
  Lx = 5.0; Ly = 5.0; Lz = 5.0;

  // header stuff for the vtk file
  string fn = "solid_angle.vtk";
  ofstream Aout (fn.c_str());
  
  Aout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
  Aout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
  Aout << "ORIGIN " << -0.5*Lx << ' ' << -0.5*Ly << ' ' << -0.5*Lz << '\n';
  Aout << "SPACING " << Lx/(Nx-1) << ' ' << Ly/(Ny-1) << ' ' << Lz/(Nz-1) << '\n';
  Aout << "POINT_DATA " << Nx*Ny*Nz << '\n';
  Aout << "SCALARS Phi float\nLOOKUP_TABLE default\n";

  //  for (int i=0; i<Nx; i++) 
  for (int k=0; k<Nz; k++) 
    {
      //      Point.xcoord = -0.5*Lx + Lx*i/(Nx-1);
      Point.zcoord = -0.5*Lz + Lz*k/(Nz-1);
      for (int j=0; j<Ny; j++) 
	{  
	  Point.ycoord = -0.5*Ly + Ly*j/(Ny-1);
	  //	  for (int k=0; k<Nz; k++) 
	  for (int i=0; i<Nx; i++) 
	    {
	      //	      Point.zcoord = -0.5*Lz + Lz*k/(Nz-1);
	      Point.xcoord = -0.5*Lx + Lx*i/(Nx-1);

	      double threshold = -0.985;
	      double ndotTmin1 = CheckThreshold(Curve1,Point);
	      //	      double ndotTmax1 = 0.0;
	      double ndotTmin2 = CheckThreshold(Curve2,Point);
	      //	      double ndotTmax2 = 0.0;
	      
	      if ((ndotTmin1 > threshold)&&(ndotTmin2 > threshold)) // more to it than this now ...
		{
		  SolidAngle = SolidAngleCalc(Curve1,Curve2,writhe1,writhe2,Point);
		}
	      else 
		{
		  SolidAngle = SolidAngleCalc2(Curve1,Curve2,writhe1,writhe2,Point);
		}
	      
	      Aout << SolidAngle << '\n';	      
	    }
	}
    }
  Aout.close();
}

// function to improve behaviour near tangent developable surface
double CheckThreshold(const knotcurve& Curve, const viewpoint& Point)
{
  int NP = Curve.knotcurve.size();
  double ndotTmin = 0.0;
  //	      double ndotTmax = 0.0;

  for (int s=0; s<NP; s++) 
    {
      double curvex = Curve.knotcurve[s].xcoord;
      double curvey = Curve.knotcurve[s].ycoord;
      double curvez = Curve.knotcurve[s].zcoord;
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = curvex - Point.xcoord;
      double viewy = curvey - Point.ycoord;
      double viewz = curvez - Point.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);		  
      double tx = Curve.knotcurve[s].tx;
      double ty = Curve.knotcurve[s].ty;
      double tz = Curve.knotcurve[s].tz;
      double ndotT = (viewx*tx + viewy*ty + viewz*tz)/dist;
      if (ndotT<ndotTmin) ndotTmin = ndotT;
      //		  if (ndotT>ndotTmax) ndotTmax = ndotT;
    }

  return ndotTmin;	      
}

inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}
