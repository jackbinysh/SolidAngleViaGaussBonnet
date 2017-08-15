#include "CurveGeometry.h"    
#include "Constants.h"    
#include <math.h>
#include <string.h>

int main (void)
{
    Link Curve;
    //  InitialiseCurve(Curve);
    InitialiseFromFile(Curve);
    // compute the writhe of the curve
    Writhe(Curve);

    // it would be good to construct a set of points to compute the twist from
    viewpoint View;
    View.xcoord = 1.0;
    View.ycoord = 0.0;
    View.zcoord = -5.0;
    // compute the twist with the view induced framing
    Twist(Curve,View);
    // compute the solid angle from the view point
    double omega = SolidAngleCalc(Curve,View);
    // put in interval [0,4pi]
    while(omega>4*M_PI) omega -= 4*M_PI;
    while(omega<0) omega += 4*M_PI;
    // 2nd version
    double omega2 = SolidAngleCalc2(Curve,View);
    // put in interval [0,4pi]
    while(omega2>4*M_PI) omega2 -= 4*M_PI;
    while(omega2<0) omega2 += 4*M_PI;

    cout << "  omega = " << omega << "  omega2 = " << omega2 << endl;

    // output the solid angle on a grid
    cout << "calculating the solid angle..." << endl;
    OutputSolidAngle(Curve);

    return 0;
}

// initialise the curve from explicit formula
//
//void InitialiseCurve(knotcurve& Curve)
//{
//    // torus knot -- ( (1 + r cos qu) cos pu , (1 + r cos qu) sin pu , r sin qu )
//    int p = 2;
//    int q = 3;
//    double r = 0.505;
//    int NP = 6000;
//    for (int s=0; s<NP; s++)
//    {
//        knotpoint Point;
//        Point.xcoord = (1.11+r*cos(2.0*M_PI*q*s/NP))*cos(2.0*M_PI*p*s/NP);
//        Point.ycoord = (1.11+r*cos(2.0*M_PI*q*s/NP))*sin(2.0*M_PI*p*s/NP);
//        Point.zcoord = r*sin(2.0*M_PI*q*s/NP);
//        Curve.knotcurve.push_back(Point);
//    }
//    /*
//    // twisted unknot -- ( cos u , sin u , r sin pu )
//    int p = 2;
//    double r = 0.3;
//    int NP = 6000;
//    for (int s=0; s<NP; s++)
//    {
//    knotpoint Point;
//    Point.xcoord = 1.03*cos(2.0*M_PI*s/NP);
//    Point.ycoord = 1.03*sin(2.0*M_PI*s/NP);
//    Point.zcoord = r*sin(2.0*M_PI*p*s/NP);
//    Curve.knotcurve.push_back(Point);
//    }
//    *//*    
//    // 7_4 knot -- ( cos 3u , sin 2u , r sin 7u )
//    double r = 0.61;
//    int NP = 12000;
//    for (int s=0; s<NP; s++)
//    {
//    knotpoint Point;
//    Point.xcoord = 0.02+0.91*cos(6.0*M_PI*s/NP);
//    Point.ycoord = 0.06+1.12*sin(4.0*M_PI*s/NP);
//    Point.zcoord = -0.0135+r*sin(14.0*M_PI*s/NP);
//    Curve.knotcurve.push_back(Point);
//    }
//    */
//    // compute the basic geometry
//    ComputeLengths(Curve);
//    ComputeTangent(Curve);
//    ComputeKappaB(Curve);
//}

// initialise the curve from an input file
void InitialiseFromFile(Link& Curve)
{

    // first up, how many components are there?
    Curve.NumComponents = NumComponents; 
    Curve.Components.resize(NumComponents);

    /*  For recording max and min input values*/
    double maxxin = 0;
    double maxyin = 0;
    double maxzin = 0;
    double minxin = 0;
    double minyin = 0;
    double minzin = 0;

    for(int i=0;i<Curve.NumComponents;i++)
    {
        stringstream ss;
        string buff,filename;

        ss.clear();
        ss.str("");
        if (Curve.NumComponents==1){ ss << knot_filename << ".txt";}
        ss << knot_filename <<"_"<< i <<  ".txt";

        filename = ss.str();
        ifstream CurveInputStream;   //knot file(s)
        CurveInputStream.open(filename.c_str());
        while(CurveInputStream.good())   //read in points for knot
        {
            double xcoord,ycoord,zcoord;    //temporary variables
            if(getline(CurveInputStream,buff))
            {
                ss.clear();
                ss.str("");
                ss << buff;
                ss >> xcoord >> ycoord >> zcoord;
            }
            else break;
            // construct a point and put it on the curve
            knotpoint Point;
            Point.xcoord = xcoord;
            Point.ycoord = ycoord;
            Point.zcoord = zcoord;
            Curve.Components[i].knotcurve.push_back(Point);
            // track max and min input values
            if(xcoord>maxxin) maxxin = xcoord;
            if(ycoord>maxyin) maxyin = ycoord;
            if(zcoord>maxzin) maxzin = zcoord;
            if(xcoord<minxin) minxin = xcoord;
            if(ycoord<minyin) minyin = ycoord;
            if(zcoord<minzin) minzin = zcoord;
        }
        CurveInputStream.close();
        // keep track of how many total points are added to the link
        Curve.NumPoints += Curve.Components[i].knotcurve.size();
    }

    // now centre and scale to a standard size
    double midpoint[3];
    midpoint[0] = 0.5*(maxxin+minxin);
    midpoint[1] = 0.5*(maxyin+minyin);
    midpoint[2] = 0.5*(maxzin+minzin);
    double scale[3];
    ScaleFunction(scale,maxxin,minxin,maxyin,minyin,maxzin,minzin);
    for(int i=0; i<Curve.NumComponents; i++)
    {
        for(int s=0; s<Curve.Components[i].knotcurve.size(); s++)
        {
            Curve.Components[i].knotcurve[s].xcoord = scale[0]*(Curve.Components[i].knotcurve[s].xcoord - midpoint[0]);
            Curve.Components[i].knotcurve[s].ycoord = scale[1]*(Curve.Components[i].knotcurve[s].ycoord - midpoint[1]);
            Curve.Components[i].knotcurve[s].zcoord = scale[2]*(Curve.Components[i].knotcurve[s].zcoord - midpoint[2]);
        }
    }

    ComputeLengths(Curve);
    ComputeTangent(Curve);
    ComputeKappaB(Curve);

    // refine the curve to increase number of points
    //  while(Curve.NumPoints < 5000*Curve.NumComponents)
    //  {
    //      RefineCurve(Curve);
    //      ComputeLengths(Curve);
    //      ComputeTangent(Curve);
    //      ComputeKappaB(Curve);

    cout << "curve has size " << Curve.NumPoints << endl;
    //  }
}

void ScaleFunction(double *scale, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin)
{
    int i;
    bool nonzeroheight[3];  //marker: true if this dimension has non zero height in stl file
    if(maxxin-minxin>0) { scale[0] = xmax/(maxxin-minxin); nonzeroheight[0] = true; }
    else { scale[0] = 1;  nonzeroheight[0] = false; }
    if(maxyin-minyin>0) { scale[1] = ymax/(maxyin-minyin); nonzeroheight[1] = true; }
    else { scale[1] = 1;  nonzeroheight[1] = false; }
    if(maxzin-minzin>0) { scale[2] = zmax/(maxzin-minzin); nonzeroheight[2] = true; }
    else { scale[2] = 1;  nonzeroheight[2] = false; }
#if PRESERVE_RATIOS
    double minscale=1000000000;
    int imin=3;
    for(i=0; i<3; i++)   //find minimum scale factor
    {
        if(scale[i] < minscale && nonzeroheight[i])
        {
            imin = i;
            minscale = scale[i];
        }
    }
    if(imin < 3)      //scale x,y, and z directions by same scale factor
    {
        for(i=0; i<3; i++) scale[i] = scale[imin];
    }
#endif
}

void RefineCurve(Link& Curve)
{
    Link NewCurve;
    NewCurve.NumComponents=Curve.NumComponents;
    NewCurve.Components.resize(Curve.NumComponents);

    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++)
        {
            knotpoint Point;
            Point.xcoord = Curve.Components[i].knotcurve[s].xcoord;
            Point.ycoord = Curve.Components[i].knotcurve[s].ycoord;
            Point.zcoord = Curve.Components[i].knotcurve[s].zcoord;
            NewCurve.Components[i].knotcurve.push_back(Point);

            knotpoint NewPoint;
            double ds = 0.5*Curve.Components[i].knotcurve[s].length;
            double x1 = Curve.Components[i].knotcurve[s].xcoord + ds*Curve.Components[i].knotcurve[s].tx + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNx;
            double x2 = Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].tx + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNx;
            NewPoint.xcoord = 0.5*(x1+x2);
            double y1 = Curve.Components[i].knotcurve[s].ycoord + ds*Curve.Components[i].knotcurve[s].ty + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNy;
            double y2 = Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].ty + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNy;
            NewPoint.ycoord = 0.5*(y1+y2);
            double z1 = Curve.Components[i].knotcurve[s].zcoord + ds*Curve.Components[i].knotcurve[s].tz + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNz;
            double z2 = Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].tz + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNz;
            NewPoint.zcoord = 0.5*(z1+z2);

            NewCurve.Components[i].knotcurve.push_back(NewPoint);
        }
        NewCurve.NumPoints += NewCurve.Components[i].knotcurve.size();
    }
    Curve = NewCurve;
}

// computes the lengths of each segment of the knot
void ComputeLengths(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)        
        {
            double dx = (Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - Curve.Components[i].knotcurve[s].xcoord);  
            double dy = (Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - Curve.Components[i].knotcurve[s].ycoord);
            double dz = (Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - Curve.Components[i].knotcurve[s].zcoord);
            double deltas = sqrt(dx*dx+dy*dy+dz*dz);
            Curve.Components[i].knotcurve[s].length = deltas;
            Curve.Components[i].length += deltas;
        }
    }
}

// computes the unit tangent of the knot
void ComputeTangent(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)        // new scheme -- central difference
        {
            double dsp = Curve.Components[i].knotcurve[s].length;
            double dsm = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            Curve.Components[i].knotcurve[s].tx = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].xcoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].xcoord;
            Curve.Components[i].knotcurve[s].ty = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].ycoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].ycoord;
            Curve.Components[i].knotcurve[s].tz = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].zcoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].zcoord;
        }
    }
}

// computes the integrand we need -- kappa B
void ComputeKappaB(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)    // new scheme -- central difference
        {
            double dsp = Curve.Components[i].knotcurve[s].length;
            double dsm = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            double kappaNx = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].tx + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].tx - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].tx;
            double kappaNy = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].ty + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].ty - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].ty;
            double kappaNz = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].tz + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].tz - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].tz;
            Curve.Components[i].knotcurve[s].kappaNx = kappaNx;
            Curve.Components[i].knotcurve[s].kappaNy = kappaNy;
            Curve.Components[i].knotcurve[s].kappaNz = kappaNz;
            Curve.Components[i].knotcurve[s].kappaBx = Curve.Components[i].knotcurve[s].ty*kappaNz - Curve.Components[i].knotcurve[s].tz*kappaNy;
            Curve.Components[i].knotcurve[s].kappaBy = Curve.Components[i].knotcurve[s].tz*kappaNx - Curve.Components[i].knotcurve[s].tx*kappaNz;
            Curve.Components[i].knotcurve[s].kappaBz = Curve.Components[i].knotcurve[s].tx*kappaNy - Curve.Components[i].knotcurve[s].ty*kappaNx;
        }
    }
}

// computes the writhe of each link component 
void Writhe(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double Wr = 0.0;
        int NP = Curve.Components[i].knotcurve.size();

        for (int s=0; s<NP; s++) // running over the knot
        {
            double a1x = Curve.Components[i].knotcurve[s].xcoord;
            double a1y = Curve.Components[i].knotcurve[s].ycoord;
            double a1z = Curve.Components[i].knotcurve[s].zcoord;
            double t1x = Curve.Components[i].knotcurve[s].tx;
            double t1y = Curve.Components[i].knotcurve[s].ty;
            double t1z = Curve.Components[i].knotcurve[s].tz;
            double ds = Curve.Components[i].knotcurve[s].length;

            for (int t=s+1; t<NP; t++) // run over all points ahead of s
            {
                double a2x = Curve.Components[i].knotcurve[t].xcoord;
                double a2y = Curve.Components[i].knotcurve[t].ycoord;
                double a2z = Curve.Components[i].knotcurve[t].zcoord;
                double t2x = Curve.Components[i].knotcurve[t].tx;
                double t2y = Curve.Components[i].knotcurve[t].ty;
                double t2z = Curve.Components[i].knotcurve[t].tz;
                double dt = Curve.Components[i].knotcurve[t].length;

                double dist = sqrt((a1x-a2x)*(a1x-a2x)+(a1y-a2y)*(a1y-a2y)+(a1z-a2z)*(a1z-a2z));
                Wr += ds*dt*((a1x-a2x)*(t1y*t2z-t1z*t2y)+(a1y-a2y)*(t1z*t2x-t1x*t2z)+(a1z-a2z)*(t1x*t2y-t1y*t2x))/(dist*dist*dist);
            }
        }
        Wr /= 2.0*M_PI; 
        Curve.Components[i].writhe=Wr;
    }
}

// computes the twist of a curve with an induced framing
void Twist(Link& Curve, const viewpoint& View)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double Angle = 0.0; 
        int NP = Curve.Components[i].knotcurve.size();

        for (int s=0; s<NP; s++) // running over the knot
        {
            double curvex = Curve.Components[i].knotcurve[s].xcoord;
            double curvey = Curve.Components[i].knotcurve[s].ycoord;
            double curvez = Curve.Components[i].knotcurve[s].zcoord;
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = curvex - View.xcoord;
            double viewy = curvey - View.ycoord;
            double viewz = curvez - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);

            // now do the Riemann integral -- I think I can rewrite this better
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ndotT = viewx*tx + viewy*ty + viewz*tz;
            double kappaBx = Curve.Components[i].knotcurve[s].kappaBx;
            double kappaBy = Curve.Components[i].knotcurve[s].kappaBy;
            double kappaBz = Curve.Components[i].knotcurve[s].kappaBz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            // and here's the integrand
            Angle += ds*ndotT*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist*dist - ndotT*ndotT);
        }

        Angle /= 2.0*M_PI;  
        Curve.Components[i].twist=Angle;
    }
}

double SolidAngleCalc(const Link& Curve, const viewpoint& View)
{
    double totalomega =0.0;
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double componentomega =0.0;
        // run over the knot 
        double Iplus = 0.0;
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++) 
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
            // now do the Riemann integral -- I think I can rewrite this better
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ndotT = viewx*tx + viewy*ty + viewz*tz;
            double kappaBx = Curve.Components[i].knotcurve[s].kappaBx;
            double kappaBy = Curve.Components[i].knotcurve[s].kappaBy;
            double kappaBz = Curve.Components[i].knotcurve[s].kappaBz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            // and here's the integrand
            Iplus += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
        }
        componentomega = 2.0*M_PI*(1.0 + Curve.Components[i].writhe)-Iplus;
        totalomega += componentomega;
    }
    return totalomega;
}

double SolidAngleCalc2(const Link& Curve, const viewpoint& View)
{
    double totalomega =0.0;
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double componentomega =0.0;
        // run over the knot 
        double Iminus = 0.0;
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++) 
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
            // now do the Riemann integral -- I think I can rewrite this better
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ndotT = viewx*tx + viewy*ty + viewz*tz;
            double kappaBx = Curve.Components[i].knotcurve[s].kappaBx;
            double kappaBy = Curve.Components[i].knotcurve[s].kappaBy;
            double kappaBz = Curve.Components[i].knotcurve[s].kappaBz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            // and here's the integrand
            Iminus += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist - ndotT);
        }
        componentomega = 2.0*M_PI*(1.0 + Curve.Components[i].writhe)-Iminus;
        totalomega += componentomega;
    }
    return totalomega;
}

void OutputSolidAngle(const Link& Curve)
{
    double SolidAngle;
    viewpoint Point;

    // header stuff for the vtk file
    string fn = "solid_angle.vtk";
    ofstream Aout (fn.c_str());

    griddata griddata;
    griddata.Nx = Nx;
    griddata.Ny = Ny;
    griddata.Nz = Nz;
    Aout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Aout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    Aout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    Aout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    Aout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    Aout << "SCALARS Phi float\nLOOKUP_TABLE default\n";

    for (int k=0; k<Nz; k++) 
    {
        Point.zcoord = z(k,griddata);
        for (int j=0; j<Ny; j++) 
        {  
            Point.ycoord = y(j,griddata) ;
            for (int i=0; i<Nx; i++) 
            {
                Point.xcoord = x(i,griddata);

                double ndotTmin,ndotTmax;
                CheckThreshold(Curve,Point,ndotTmin,ndotTmax);

                // is ndotT very close to 1? if so, better us 1/(1+ndotT)
                if(fabs(1 - ndotTmax) < fabs(-1-ndotTmin)) 
                {
                    SolidAngle = SolidAngleCalc(Curve,Point);
                }
                else 
                {
                    SolidAngle = SolidAngleCalc2(Curve,Point);
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
void CheckThreshold(const Link& Curve, const viewpoint& Point, double& ndotTmin,double& ndotTmax)
{
    // set these guys at impossible bounds
    ndotTmin = 2.0;
    ndotTmax = -2.0;
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();

        for (int s=0; s<NP; s++) 
        {
            double curvex = Curve.Components[i].knotcurve[s].xcoord;
            double curvey = Curve.Components[i].knotcurve[s].ycoord;
            double curvez = Curve.Components[i].knotcurve[s].zcoord;
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = curvex - Point.xcoord;
            double viewy = curvey - Point.ycoord;
            double viewz = curvez - Point.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);		  
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ndotT = (viewx*tx + viewy*ty + viewz*tz)/dist;
            if (ndotT<ndotTmin) ndotTmin = ndotT;
            if (ndotT>ndotTmax) ndotTmax = ndotT;
        }
    }
}

inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}
inline double x(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Nx/2.0)*h;
}
inline double y(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Ny/2.0)*h;
}
inline double z(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Nz/2.0)*h;
}
