#include "LinkGeometry.h"    
#include <math.h>
#include <string.h>

int main (void)
{
  int L = 2; // number of link components
  vector<knotcurve> Link (L);
  InitialiseLink(Link);

  // compute the writhe of each component
  double writhe[L];
  for (int i=0; i<L; i++)
    {
      writhe[i] = Writhe(Link[i]);
    }

  // it would be good to construct a set of points to compute the twist from
  viewpoint View;
  View.xcoord = 5.0;
  View.ycoord = -5.0;
  View.zcoord = 5.0;

  // compute the twist of each component with the projective framing
  double twist[L];
  for (int i=0; i<L; i++)
    {
      twist[i] = Twist(Link[i],View);
    }

  // compute the solid angle from the view point
  double omega1 = 0.0;
  for (int i=0; i<L; i++)
    {
      omega1 += 2.0*M_PI*(1.0 + writhe[i]) - SolidAngleCalc(Link[i],View);
    }
  // put in interval [0,4pi]
  while(omega1>4*M_PI) omega1 -= 4*M_PI;
  while(omega1<0) omega1 += 4*M_PI;

  // there are really NL^2 choices not just 2
  double omega2 = 0.0;
  for (int i=0; i<L; i++)
    {
      omega2 += 2.0*M_PI*(1.0 - writhe[i]) - SolidAngleCalc2(Link[i],View);
    }
  // put in interval [0,4pi]
  while(omega2>4*M_PI) omega2 -= 4*M_PI;
  while(omega2<0) omega2 += 4*M_PI;

  for (int i=0; i<L; i++)
    {
      cout << "curve " << i << ":  Wr = " << writhe[i] << "  Tw = " << twist[i] << "  SL = " << writhe[i]+twist[i] << endl; 
    }
  cout << "omega1 = " << omega1 << "  omega2 = " << omega2 << endl;

  // output the solid angle on a grid
  cout << "calculating the solid angle..." << endl;
  OutputSolidAngle(Link,writhe);

  return 0;
}

// initialise the link 
void InitialiseLink(vector<knotcurve>& Link)
{
  int L = Link.size(); // number of link components
  double x0 = 0.024; // small shift of origin;
  double y0 = 0.036;
  double z0 = -0.022;
  // torus knot -- ( (1 + r cos qu) cos pu , (1 + r cos qu) sin pu , r sin qu )
  int p = 1;
  int q = 2;
  double r = 0.57;
  int NP = 2000;
  for (int s=0; s<NP; s++)
    {
      knotpoint Point;
      Point.xcoord = x0 + (1.24+r*cos(2.0*M_PI*q*s/NP))*cos(2.0*M_PI*p*s/NP);
      Point.ycoord = y0 + (1.24+r*cos(2.0*M_PI*q*s/NP))*sin(2.0*M_PI*p*s/NP);
      Point.zcoord = z0 + r*sin(2.0*M_PI*q*s/NP);
      Link[0].knotcurve.push_back(Point);
      // second link component
      //      Point.xcoord = x0 + (1.24-r*cos(2.0*M_PI*q*s/NP))*cos(2.0*M_PI*p*s/NP);
      //      Point.ycoord = y0 + (1.24-r*cos(2.0*M_PI*q*s/NP))*sin(2.0*M_PI*p*s/NP);
      //      Point.zcoord = z0 - r*sin(2.0*M_PI*q*s/NP);
      // second link component -- reverse orientation
      Point.xcoord = x0 + (1.24-r*cos(2.0*M_PI*q*s/NP))*cos(2.0*M_PI*p*s/NP);
      Point.ycoord = y0 - (1.24-r*cos(2.0*M_PI*q*s/NP))*sin(2.0*M_PI*p*s/NP);
      Point.zcoord = z0 + r*sin(2.0*M_PI*q*s/NP);
      Link[1].knotcurve.push_back(Point);
    }

  // calculate geometry
  for (int i=0; i<L; i++) 
    {
      ComputeLengths(Link[i]);
      ComputeTangent(Link[i]);
      ComputeKappaB(Link[i]);
    }
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
      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve.knotcurve[s].tx;
      double ty = Curve.knotcurve[s].ty;
      double tz = Curve.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve.knotcurve[s].kappaBx;
      double kappaBy = Curve.knotcurve[s].kappaBy;
      double kappaBz = Curve.knotcurve[s].kappaBz;
      double ds = 0.5*(Curve.knotcurve[s].length+Curve.knotcurve[incp(s,-1,NP)].length);
      // and here's the integrand
      Angle += ds*ndotT*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist*dist - ndotT*ndotT);
    }

  Angle /= 2.0*M_PI;  
  return Angle;
}

double SolidAngleCalc(const knotcurve& Curve, const viewpoint& View)
{
  double Iplus = 0.0;
  int NP = Curve.knotcurve.size();

 // run over the knot 
  for (int s=0; s<NP; s++) 
    {
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = Curve.knotcurve[s].xcoord - View.xcoord;
      double viewy = Curve.knotcurve[s].ycoord - View.ycoord;
      double viewz = Curve.knotcurve[s].zcoord - View.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve.knotcurve[s].tx;
      double ty = Curve.knotcurve[s].ty;
      double tz = Curve.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve.knotcurve[s].kappaBx;
      double kappaBy = Curve.knotcurve[s].kappaBy;
      double kappaBz = Curve.knotcurve[s].kappaBz;
      double ds = 0.5*(Curve.knotcurve[s].length+Curve.knotcurve[incp(s,-1,NP)].length);
      // and here's the integrand
      Iplus += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
    }
    
  return Iplus;
}

double SolidAngleCalc2(const knotcurve& Curve, const viewpoint& View)
{
  double Iminus = 0.0;
  int NP = Curve.knotcurve.size();

 // run over the knot 
  for (int s=0; s<NP; s++) 
    {
      // define the view vector -- n = (Curve - View)/|Curve - View|
      double viewx = Curve.knotcurve[s].xcoord - View.xcoord;
      double viewy = Curve.knotcurve[s].ycoord - View.ycoord;
      double viewz = Curve.knotcurve[s].zcoord - View.zcoord;
      double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
      // now do the Riemann integral -- I think I can rewrite this better
      double tx = Curve.knotcurve[s].tx;
      double ty = Curve.knotcurve[s].ty;
      double tz = Curve.knotcurve[s].tz;
      double ndotT = viewx*tx + viewy*ty + viewz*tz;
      double kappaBx = Curve.knotcurve[s].kappaBx;
      double kappaBy = Curve.knotcurve[s].kappaBy;
      double kappaBz = Curve.knotcurve[s].kappaBz;
      double ds = 0.5*(Curve.knotcurve[s].length+Curve.knotcurve[incp(s,-1,NP)].length);
      // and here's the integrand
      Iminus += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist - ndotT);
    }
    
  return Iminus;
}

void OutputSolidAngle(const vector<knotcurve>& Link, const double *writhe)
{
  int L = Link.size(); // number of link components
  double SolidAngle;
  viewpoint Point;

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

  for (int k=0; k<Nz; k++) 
    {
      Point.zcoord = -0.5*Lz + Lz*k/(Nz-1);
      for (int j=0; j<Ny; j++) 
	{  
	  Point.ycoord = -0.5*Ly + Ly*j/(Ny-1);
	  for (int i=0; i<Nx; i++) 
	    {
	      Point.xcoord = -0.5*Lx + Lx*i/(Nx-1);

	      double threshold = -0.985;
	      double SolidAngle = 0.0;
	      for (int l=0; l<L; l++)
		{
		  double ndotTmin = CheckThreshold(Link[l],Point);
		  //	      double ndotTmax = 0.0;	      
		  if (ndotTmin > threshold) 
		    {
		      SolidAngle += 2.0*M_PI*(1.0 + writhe[l]) - SolidAngleCalc(Link[l],Point);
		    }
		  else 
		    {
		      SolidAngle += 2.0*M_PI*(1.0 - writhe[l]) - SolidAngleCalc2(Link[l],Point);
		    }
		}

	      // put in the interval [0,4pi]
	      while(SolidAngle>4*M_PI) SolidAngle -= 4*M_PI;
	      while(SolidAngle<0) SolidAngle += 4*M_PI;
	      
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
