#include "SolidAngleGaussBonnet.h"    //contains user defined variables for the simulation, and the parameters used 
#include <omp.h>
#include <math.h>
#include <string.h>

int main (void)
{
    griddata griddata;
    griddata.Nx = initialNx;
    griddata.Ny = initialNy;
    griddata.Nz = initialNz;
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    vector<double>phi(Nx*Ny*Nz);  //scalar potential

    DualConePhiCalc(phi,griddata); 
    return 0;
}

void DualConePhiCalc(vector<double>&phi, const griddata& griddata)
{
    //first up, initialise the knotcurve
    knotcurve InitialisationCurve;

    stringstream ss;
    string buff,filename;

    ss.clear();
    ss.str("");
    ss << knot_filename << ".txt";

    filename = ss.str();
    ifstream CurveInputStream;   //knot file(s)
    CurveInputStream.open(filename.c_str());
    /*  For recording max and min input values*/
    double maxxin = 0;
    double maxyin = 0;
    double maxzin = 0;
    double minxin = 0;
    double minyin = 0;
    double minzin = 0;
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
        InitialisationCurve.knotcurve.push_back(Point);
        if(xcoord>maxxin) maxxin = xcoord;
        if(ycoord>maxyin) maxyin = ycoord;
        if(zcoord>maxzin) maxzin = zcoord;
        if(xcoord<minxin) minxin = xcoord;
        if(ycoord<minyin) minyin = ycoord;
        if(zcoord<minzin) minzin = zcoord;
    }
    CurveInputStream.close();

    // now scale it to fit
    double scale[3];
    double midpoint[3];
    scalefunction(scale,midpoint,maxxin,minxin,maxyin,minyin,maxzin,minzin);
    for(int s=0;s<InitialisationCurve.knotcurve.size();s++)
    {
        InitialisationCurve.knotcurve[s].xcoord = scale[0]*(InitialisationCurve.knotcurve[s].xcoord - midpoint[0]);
        InitialisationCurve.knotcurve[s].ycoord = scale[1]*(InitialisationCurve.knotcurve[s].ycoord - midpoint[1]);
        InitialisationCurve.knotcurve[s].zcoord = scale[2]*(InitialisationCurve.knotcurve[s].zcoord - midpoint[2]);
    }

#pragma omp parallel default(none) shared(InitialisationCurve,phi,griddata)
    {
        // now for the guts of the thing - construct the dual curve
        knotcurve ProjectedCurve;
        knotcurve DualCurve;
        int intersectioncount = 0;

        int Nx = griddata.Nx;
        int Ny = griddata.Ny;
        int Nz = griddata.Nz;
#pragma omp for 
        for(int i=0;i<Nx;i++)
        {
            for(int j=0; j<Ny; j++)
            {
                for(int k=0; k<Nz; k++)
                {

                    // first ensure the curves are clear from the list iteration
                    ProjectedCurve.knotcurve.clear();
                    DualCurve.knotcurve.clear();
                    DualCurve.length = 0 ;

                    // first, take the curve and project it onto the unit sphere around our observation point
                    for(int s=0;s<InitialisationCurve.knotcurve.size();s++)
                    {

                        double initialx = InitialisationCurve.knotcurve[s].xcoord;
                        double initialy = InitialisationCurve.knotcurve[s].ycoord;
                        double initialz = InitialisationCurve.knotcurve[s].zcoord;

                        double projxcoord = initialx - x(i,griddata);
                        double projycoord = initialy - y(j,griddata);
                        double projzcoord = initialz - z(k,griddata);

                        double mag = sqrt(projxcoord*projxcoord +projycoord*projycoord +projzcoord*projzcoord);
                        projxcoord /= mag;
                        projycoord /= mag;
                        projzcoord /= mag;

                        knotpoint Point;
                        Point.xcoord = projxcoord;
                        Point.ycoord = projycoord;
                        Point.zcoord = projzcoord;
                        ProjectedCurve.knotcurve.push_back(Point);
                    }

                    // get segment lengths, tangent vectors, for the projected curve
                    float topdeltas = 0;
                    for(int s=0;s<ProjectedCurve.knotcurve.size();s++)
                    {
                        int NP = ProjectedCurve.knotcurve.size();
                        double dx = (ProjectedCurve.knotcurve[incp(s,1,NP)].xcoord - ProjectedCurve.knotcurve[incp(s,0,NP)].xcoord);   //central diff as a is defined on the points
                        double dy = (ProjectedCurve.knotcurve[incp(s,1,NP)].ycoord - ProjectedCurve.knotcurve[incp(s,0,NP)].ycoord);
                        double dz = (ProjectedCurve.knotcurve[incp(s,1,NP)].zcoord - ProjectedCurve.knotcurve[incp(s,0,NP)].zcoord);
                        double deltas = sqrt(dx*dx+dy*dy+dz*dz);
                        ProjectedCurve.knotcurve[s].tx= dx/(deltas);
                        ProjectedCurve.knotcurve[s].ty= dy/(deltas);
                        ProjectedCurve.knotcurve[s].tz= dz/(deltas);
                        ProjectedCurve.knotcurve[s].length = deltas;
                        if (deltas > topdeltas){ topdeltas = deltas;}
                    }

                    int intersectioncount = 0;
                    for(int s=0;s<ProjectedCurve.knotcurve.size();s++)
                    {
                        int NP = ProjectedCurve.knotcurve.size();
                        float a1[3];
                        float a2[3];
                        a1[0]= ProjectedCurve.knotcurve[s].xcoord;
                        a1[1]= ProjectedCurve.knotcurve[s].ycoord;
                        a1[2]= ProjectedCurve.knotcurve[s].zcoord;
                        a2[0]= ProjectedCurve.knotcurve[incp(s,1,NP)].xcoord;
                        a2[1]= ProjectedCurve.knotcurve[incp(s,1,NP)].ycoord;
                        a2[2]= ProjectedCurve.knotcurve[incp(s,1,NP)].zcoord;

                        ProjectedCurve.knotcurve[s].length = 0;
                        for(int t=s+1;t<ProjectedCurve.knotcurve.size();t++)
                        {
                            float b1[3];
                            float b2[3];
                            b1[0]= ProjectedCurve.knotcurve[t].xcoord;
                            b1[1]= ProjectedCurve.knotcurve[t].ycoord;
                            b1[2]= ProjectedCurve.knotcurve[t].zcoord;
                            b2[0]= ProjectedCurve.knotcurve[incp(t,1,NP)].xcoord;
                            b2[1]= ProjectedCurve.knotcurve[incp(t,1,NP)].ycoord;
                            b2[2]= ProjectedCurve.knotcurve[incp(t,1,NP)].zcoord;

                            bool intersection = false;
                            if(s!=t && s!=incp(t,1,NP) && s!=incp(t,-1,NP))
                            {
                                // first up, a quick and dirty test:
                                if( ((a1[0] - b1[0])*(a1[0] - b1[0])+(a1[1] - b1[1])*(a1[1] - b1[1])+(a1[2] - b1[2])*(a1[2] - b1[2])) < 4*topdeltas*topdeltas)
                                {
                                    intersection = GeodesicIntersection(a1,a2,b1,b2);

                                }
                            }
                            intersectioncount += intersection;
                        }
                    }

                    // now from the projected curve, construct the dual. In the notation below, I follow Mark Levi's convention of using n to denote position on the sphere
                    for(int s=0;s<ProjectedCurve.knotcurve.size();s++)
                    {

                        double tx = ProjectedCurve.knotcurve[s].tx ;
                        double ty = ProjectedCurve.knotcurve[s].ty ;
                        double tz = ProjectedCurve.knotcurve[s].tz ;

                        double nx = ProjectedCurve.knotcurve[s].xcoord;
                        double ny = ProjectedCurve.knotcurve[s].ycoord;
                        double nz = ProjectedCurve.knotcurve[s].zcoord;
                        // get dual curve positions, with cross product
                        double nstarx = ty*nz - tz*ny;
                        double nstary = tz*nx - tx*nz;
                        double nstarz = tx*ny - ty*nx;

                        knotpoint Point;
                        Point.xcoord = nstarx; 
                        Point.ycoord = nstary;
                        Point.zcoord = nstarz;
                        Point.tx = -tx; 
                        Point.ty = -ty; 
                        Point.tz = -tz; 
                        DualCurve.knotcurve.push_back(Point);
                    }

                    // compute the dual curves length. We do this using the exact formula in Mark Levi's paper - got to be careful with signs here
                    for(int s=0;s<DualCurve.knotcurve.size();s++)
                    {
                        int NP = DualCurve.knotcurve.size();

                        // computing the quantities in Levi's eqn (1.1)
                        double dx = (DualCurve.knotcurve[incp(s,1,NP)].xcoord - DualCurve.knotcurve[incp(s,0,NP)].xcoord);   //central diff as a is defined on the points
                        double dy = (DualCurve.knotcurve[incp(s,1,NP)].ycoord - DualCurve.knotcurve[incp(s,0,NP)].ycoord);
                        double dz = (DualCurve.knotcurve[incp(s,1,NP)].zcoord - DualCurve.knotcurve[incp(s,0,NP)].zcoord);
                        double deltas = sqrt(dx*dx+dy*dy+dz*dz);

                        double dndsstarx = dx/deltas; 
                        double dndsstary = dy/deltas;  
                        double dndsstarz = dz/deltas; 

                        double tstarx =  DualCurve.knotcurve[s].tx;   
                        double tstary =  DualCurve.knotcurve[s].ty;
                        double tstarz =  DualCurve.knotcurve[s].tz;

                        int triadsign =(tstarx*dndsstarx+tstary*dndsstary+tstarz*dndsstarz)<0? -1 : 1;

                        // I want to try using the actual geodesic length:
                        double n1x = DualCurve.knotcurve[incp(s,0,NP)].xcoord;
                        double n1y = DualCurve.knotcurve[incp(s,0,NP)].ycoord;
                        double n1z = DualCurve.knotcurve[incp(s,0,NP)].zcoord;
                        double n2x = DualCurve.knotcurve[incp(s,1,NP)].xcoord;
                        double n2y = DualCurve.knotcurve[incp(s,1,NP)].ycoord;
                        double n2z = DualCurve.knotcurve[incp(s,1,NP)].zcoord;

                        double ax = n1y*n2z - n1z*n2y;
                        double ay = n1z*n2x - n1x*n2z;
                        double az = n1x*n2y - n1y*n2x;

                        double geodesicdistance = atan2(sqrt(ax*ax+ay*ay+az*az),(n1x*n2x+n1y*n2y+n1z*n2z));

                        DualCurve.length += triadsign*geodesicdistance;
                    }

                    // clean up, and set phi
                    int parity = intersectioncount%2;
                    double SolidAngle = DualCurve.length + 2*M_PI*parity;
                    while(SolidAngle>4*M_PI) SolidAngle-= 4*M_PI;
                    while(SolidAngle<0) SolidAngle += 4*M_PI;

                    int n = pt(i,j,k,griddata);
                    phi[n] = SolidAngle;
                }
            }
        }
    }
    cout << "Printing phi...\n";
    print_phi(phi,griddata);
    return;
}

void scalefunction(double *scale, double *midpoint, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin)
{
    int i;
    bool nonzeroheight[3];  //marker: true if this dimension has non zero height in stl file
    if(maxxin-minxin>0) { scale[0] = xmax/(maxxin-minxin); nonzeroheight[0] = true; }
    else { scale[0] = 1;  nonzeroheight[0] = false; }
    if(maxyin-minyin>0) { scale[1] = ymax/(maxyin-minyin); nonzeroheight[1] = true; }
    else { scale[1] = 1;  nonzeroheight[1] = false; }
    if(maxzin-minzin>0) { scale[2] = zmax/(maxzin-minzin); nonzeroheight[2] = true; }
    else { scale[2] = 1;  nonzeroheight[2] = false; }
    //double p1x,p1y,p1z,p2x,p2y,p2z,nx,ny,nz;
    midpoint[0] = 0.5*(maxxin+minxin);
    midpoint[1] = 0.5*(maxyin+minyin);
    midpoint[2] = 0.5*(maxzin+minzin);
#if PRESERVE_RATIOS
    double minscale=1000000000;
    int imin=3;
    for(i = 0;i<3;i++)   //find minimum scale factor
    {
        if(scale[i] < minscale && nonzeroheight[i])
        {
            imin = i;
            minscale = scale[i];
        }
    }
    if(imin < 3)      //scale x,y, and z directions by same scale factor
    {
        for(i = 0;i<3;i++) scale[i] = scale[imin];
    }
#endif
}

void print_phi( vector<double>&phi, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n;
    string fn = "phi.vtk";

    ofstream Bout (fn.c_str());

    Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Bout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    Bout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    Bout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    Bout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    Bout << "SCALARS Phi float\nLOOKUP_TABLE default\n";
    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                Bout << phi[n] << '\n';
            }
        }
    }
    Bout.close();
}

inline bool GeodesicIntersection(float* a1,float* a2,float* b1,float* b2)
{
    float epsilon = 0.00001f;
    bool intersection = 0;

    float a[3];
    Cross(a1,a2,a);
    float b[3];
    Cross(b1,b2,b);
    float line[3];
    Cross(a,b,line);
    Norm(line);
    // the vector line intersects the sphere twice - at the points {line, -line}. we now test if these points lie on the arcs
    float test1 = fabs(GeodesicDistance(a1,a2)-GeodesicDistance(a1,line)-GeodesicDistance(a2,line));
    float test2 = fabs(GeodesicDistance(b1,b2)-GeodesicDistance(b1,line)-GeodesicDistance(b2,line));
    if(test1 < epsilon && test2 <epsilon)
    {
        intersection=true;
    }
    else
    {
        line[0] = -line[0];
        line[1] = -line[1];
        line[2] = -line[2];

        test1 = fabs(GeodesicDistance(a1,a2)-GeodesicDistance(a1,line)-GeodesicDistance(a2,line));
        test2 = fabs(GeodesicDistance(b1,b2)-GeodesicDistance(b1,line)-GeodesicDistance(b2,line));

        if(test1 < epsilon && test2 <epsilon){intersection=true;}
    }
    //
    return intersection;
}

inline float GeodesicDistance(float* a, float* b)
{
    float result[3];
    Cross(a,b,result);
    float distance = atan2(sqrt(dot(result,result)),dot(a,b));
    return distance;
}

inline void Cross(float* a, float* b, float* result)
{
    result[0]= a[1]*b[2] - a[2]*b[1];
    result[1]= a[2]*b[0] - a[0]*b[2];
    result[2]= a[0]*b[1] - a[1]*b[0];
}

inline float dot(float* a, float* b)
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline void Norm(float* a)
{
    float mag = sqrt(dot(a,a));
    a[0] /= mag;
    a[1] /= mag;
    a[2] /= mag;
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
inline int sign(int i)
{
    if(i==0) return 0;
    else return i/abs(i);
}
