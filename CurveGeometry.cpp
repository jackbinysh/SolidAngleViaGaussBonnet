#include "CurveGeometry.h"    
#include "Constants.h"    
#include <math.h>
#include <string.h>


int main (void)
{
    Link Curve;
    InitialiseFromFile(Curve);
    cout << "calculating the solid angle..." << endl;
    OutputSolidAngle(Curve);

    return 0;
}

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
        if (Curve.NumComponents==1)
        {
            ss << knot_filename << ".txt";
        }
        else
        {
            ss << knot_filename <<"_"<< i <<  ".txt";
        }

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
    ComputeWrithe(Curve);
    cout << "curve has size " << Curve.NumPoints << endl;
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
            Curve.Components[i].knotcurve[s].curvature = sqrt(kappaNx*kappaNx + kappaNy*kappaNy + kappaNz*kappaNz);
            double curvature = Curve.Components[i].knotcurve[s].curvature;
        }
    }
}

// computes the writhe of each link component 
void ComputeWrithe(Link& Curve)
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

// the regularised one
double SolidAngleCalcR(const Link& Curve, const viewpoint& View)
{

    double totalomega =0;
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double Integral= 0.0;
        double correctedpartIntegral = 0.0;
        double uncorrectedpartIntegral = 0.0;

        int M = 8; // even integer for defining a window around the cusp point
        int NP = Curve.Components[i].knotcurve.size();
        int window[NP]; // window of integration

        for (int s=0; s<NP; s++) {window[s]=1;} // there must be something smarter than this!!

        for (int s=0; s<NP; s++)
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
            double ndotT =  (viewx*Curve.Components[i].knotcurve[s].tx + viewy*Curve.Components[i].knotcurve[s].ty + viewz*Curve.Components[i].knotcurve[s].tz)/dist;
            // calculate dot product with the Frenet normal
            double ndotN = viewx*Curve.Components[i].knotcurve[s].kappaNx + viewy*Curve.Components[i].knotcurve[s].kappaNy + viewz*Curve.Components[i].knotcurve[s].kappaNz;
            ndotN /= (dist*Curve.Components[i].knotcurve[s].curvature); // shouldn't really be necessary
            // repeat for the next point along
            double viewx2 = Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - View.xcoord;
            double viewy2 = Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - View.ycoord;
            double viewz2 = Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - View.zcoord;
            double dist2 = sqrt(viewx2*viewx2 + viewy2*viewy2 + viewz2*viewz2);      
            double ndotN2 = viewx2*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNx + viewy2*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNy + viewz2*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNz;
            ndotN2 /= (dist2*Curve.Components[i].knotcurve[incp(s,1,NP)].curvature); // shouldn't really be necessary
            // only continue if ndotN changes sign and ndotT is negative
            if ((ndotN*ndotN2<0) && (ndotT<0))
            {
                // tag the point of the curve
                int sp = s;
                //if ((-ndotN2/ndotN)<1) sp = incp(s,1,NP);
                // compute ndotT and ndotB
                viewx = Curve.Components[i].knotcurve[sp].xcoord - View.xcoord;
                viewy = Curve.Components[i].knotcurve[sp].ycoord - View.ycoord;
                viewz = Curve.Components[i].knotcurve[sp].zcoord - View.zcoord;
                dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
                // note the minus sign introduced here -- makes the quantity positive
                ndotT = -(viewx*Curve.Components[i].knotcurve[sp].tx + viewy*Curve.Components[i].knotcurve[sp].ty + viewz*Curve.Components[i].knotcurve[sp].tz)/dist;
                double viewx2 = Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - View.xcoord;
                double viewy2 = Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - View.ycoord;
                double viewz2 = Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - View.zcoord;
                double dist2 = sqrt(viewx2*viewx2 + viewy2*viewy2 + viewz2*viewz2);      

                double ndotB = -(viewx*Curve.Components[i].knotcurve[sp].kappaBx + viewy*Curve.Components[i].knotcurve[sp].kappaBy + viewz*Curve.Components[i].knotcurve[sp].kappaBz)/(dist*Curve.Components[i].knotcurve[sp].curvature);
                double ndotB2 = -(viewx2*Curve.Components[i].knotcurve[incp(sp,1,NP)].kappaBx + viewy2*Curve.Components[i].knotcurve[incp(sp,1,NP)].kappaBy + viewz2*Curve.Components[i].knotcurve[incp(sp,1,NP)].kappaBz)/(dist2*Curve.Components[i].knotcurve[incp(sp,1,NP)].curvature);
                // check if a threshold is exceeded
                double threshold = Curve.Components[i].knotcurve[sp].curvature*M*(Curve.Components[i].length/NP)/(2.0*sqrt(2.0)*sqrt(1.0/ndotT - 1.0));
                if (threshold > 1)
                {
                    // the lambda which, when put into a linear interpolatio of s n, would give the location of the zero in ndotN.
                    double lambda = -(ndotN/(ndotN2 - ndotN));
                    double s0 =sp +lambda;

                    // an interpolation of the quanities needed in the integral:
                    double k0 = Curve.Components[i].knotcurve[sp].curvature + lambda*(Curve.Components[i].knotcurve[incp(sp,1,NP)].curvature-Curve.Components[i].knotcurve[sp].curvature);
                    double ndotB0 = ndotB + lambda*(ndotB2-ndotB);

                    // value of the integral over the narrow Lorentzian peak

                    // ok - this is a hack . if you check, you will find that it matters that the arclengths are not constant, and are not the same as the total average, when it comes to evaluting the below integral. rather than properly implementing getting x0 and x1 for varying arc length, im just setting the "average" to the value at sp, assuming its constant over the lorentzian peak - but note that the value at sp may be quite different from the overall average.
                    double avgds = Curve.Components[i].knotcurve[sp].length;

                    double x0 = (sp - M/2 +1-s0)*avgds;
                    double x1 = (sp + M/2  -s0)*avgds;
                    double theta = asin(ndotB0);
                    double Iplus0 = 2*(atan(k0*x0/theta) -atan(k0*x1/theta));
                    correctedpartIntegral +=Iplus0;
                    // set an excluded window about sp
                    // do the ones above
                    for (int q=0; q<=M/2; q++) 
                    {
                        double y1 = (sp + q -s0)*avgds;
                        window[incp(sp,q,NP)] = 0;
                    }
                    // do the ones below 
                    for (int q=-1; q>=(-M/2)+1; q--) 
                    {
                        double y0 = (sp + q -s0)*avgds;
                        window[incp(sp,q,NP)] = 0;
                    }
                }
            }
        }

        // now do the remaining integral as usual
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

            double factor =1;
            if (window[s]==0)
            {
                factor = 0;
                // if we are are the edge of the excluded region
                if(window[incp(s,1,NP)]==1||window[incp(s,-1,NP)]==1)
                {
                    factor = 0.5;
                }
            }
            uncorrectedpartIntegral += ds*factor*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
        }

        Integral = correctedpartIntegral+uncorrectedpartIntegral;
        double componentomega = 2.0*M_PI*(1.0 + Curve.Components[i].writhe)-Integral;
        totalomega += componentomega;
    }
    return totalomega;
}

void OutputSolidAngle(const Link& Curve)
{
    vector<double>phi(Nx*Ny*Nz);
    griddata griddata;
    griddata.Nx = Nx;
    griddata.Ny = Ny;
    griddata.Nz = Nz;

#pragma omp parallel default(none) shared(phi,Curve,griddata)
    {
        double SolidAngle;
        viewpoint Point;

#pragma omp for
        for (int i=0; i<Nx; i++) 
        {
            for (int j=0; j<Ny; j++) 
            {
                for (int k=0; k<Nz; k++) 
                {
                    Point.xcoord = x(i,griddata);
                    Point.ycoord = y(j,griddata) ;
                    Point.zcoord = z(k,griddata);
                    int n = pt(i,j,k,griddata);

                    SolidAngle = SolidAngleCalcR(Curve,Point);

                    // put in the interval [0,4pi]
                    while(SolidAngle>4*M_PI) SolidAngle -= 4*M_PI;
                    while(SolidAngle<0) SolidAngle += 4*M_PI;

                    phi[n]=SolidAngle;
                }
            }
        }

    }

    // header stuff for the vtk file
    string fn = "solid_angle.vtk";
    ofstream Aout (fn.c_str());
    Aout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Aout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    Aout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    Aout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    Aout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    Aout << "SCALARS Phi float\nLOOKUP_TABLE default\n";
    for(int k=0; k<Nz; k++)
    {
        for(int j=0; j<Ny; j++)
        {
            for(int i=0; i<Nx; i++)
            {
                int n = pt(i,j,k,griddata);
                Aout << phi[n] << '\n';
            }
        }
    }
    Aout.close();
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
inline  int pt( int i,  int j,  int k,const griddata& griddata)       //convert i,j,k to single index
{
    return (i*griddata.Ny*griddata.Nz+j*griddata.Nz+k);
}
