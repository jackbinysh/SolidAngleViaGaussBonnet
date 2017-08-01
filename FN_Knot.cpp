/* Fitzhugh-Nagumo reaction diffusion simulation with arbitrary vortex lines
   OPENMP VERSION
   Created by Carl Whitfield
   Last modified 03/01/17

   Operational order of the code:
   1) The code takes as input an stl file (defined in knot_filename) which defines an orientable surface with a boundary.
   2) This surface is scaled to fill a box of size xmax x ymax x zmax.
   3) A nunmerical integral is performed to calculate a phase field (phi_calc) on the 3D grid which winds around the boundary of the surface.
   4) This is then used to initialise the Fitzhugh-Nagumo set of partial differential equations such that on initialisation (uv_initialise):
   u = 2cos(phi) - 0.4        v = sin(phi) - 0.4
   The pde's used are
   dudt = (u - u^3/3 - v)/epsilon + Del^2 u
   dvdt = epsilon*(u + beta - gam v)
   5) The update method is Runge-Kutta fourth order (update_uv) unless RK4 is set to 0, otherwise Euler forward method is used.
   6) A parametric curve for the knot is found at each unit T


   The parameters epsilon, beta and gam are set to give rise to scroll waves (see Sutcliffe, Winfree, PRE 2003 and Maucher Sutcliffe PRL 2016) which eminate from a closed curve, on initialisation this curve corresponds to the boundary of the surface in the stl file.

   See below for various options to start the code from previous outputted data.*/
#include "FN_Knot.h"    //contains user defined variables for the simulation, and the parameters used 
#include "TriCubicInterpolator.h"    //contains user defined variables for the simulation, and the parameters used 
#include <omp.h>
#include <math.h>
#include <string.h>
//includes for the signal processing
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

int main (void)
{
    griddata griddata;
    griddata.Nx = initialNx;
    griddata.Ny = initialNy;
    griddata.Nz = initialNz;
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    // all major allocations are here
    // the main data storage arrays, contain info associated with the grid
    vector<double>phi(Nx*Ny*Nz);  //scalar potential
    vector<double>u(Nx*Ny*Nz);
    vector<double>v(Nx*Ny*Nz);
    vector<double>ucvx(Nx*Ny*Nz);
    vector<double>ucvy(Nx*Ny*Nz);
    vector<double>ucvz(Nx*Ny*Nz);
    vector<double>ucvmag(Nx*Ny*Nz);// mod(grad u cross grad v)
    vector<double>ku(4*Nx*Ny*Nz);
    vector<double>kv(4*Nx*Ny*Nz);
    // objects to hold information about the knotcurve we find, andthe surface we read in
    vector<knotcurve > knotcurves; // a structure containing some number of knot curves, each curve a list of knotpoints
    vector<knotcurve > knotcurvesold; // a structure containing some number of knot curves, each curve a list of knotpoints
    vector<triangle> knotsurface;    //structure for storing knot surface coordinates
    // GSL initialization
    const gsl_multimin_fminimizer_type *Type;
    gsl_multimin_fminimizer *minimizerstate;
    Type = gsl_multimin_fminimizer_nmsimplex2;
    minimizerstate = gsl_multimin_fminimizer_alloc (Type,2);


    if (option == FROM_PHI_FILE)
    {
        cout << "Reading input file...\n";
        phi_file_read(phi,griddata);
    }
    else
    {
        if(option == FROM_UV_FILE)
        {
            cout << "Reading input file...\n";
            if(uvfile_read(u,v,ku,kv, ucvx,ucvy,ucvz,griddata)) return 1;
        }
        else
        {
            if(option == FROM_FUNCTION)
            { 
                phi_calc_manual(phi,griddata);
            }
            else
            {
                if(option == FROM_SURFACE_FILE) 
                {
                    //Initialise knot
                    if(initialise_knot(knotsurface)==0)
                    {
                        cout << "Error reading input option. Aborting...\n";
                        return 1;
                    }

                    if(option == FROM_SURFACE_FILE) cout << "Total no. of surface points: ";
                    cout << knotsurface.size() << '\n';

                    //Calculate phi for initial conditions
                    phi_calc(phi,knotsurface,griddata);
                }
                else
                {
                    if(option == FROM_CURVE_FILE)
                    {
                        DualConePhiCalc(phi,griddata); 
                    }
                }
            }
        }
    }
    vector<triangle> ().swap(knotsurface);   //empty knotsurface memory
    if(option!=FROM_UV_FILE)
    {
        cout << "Calculating u and v...\n";
        uv_initialise(phi,u,v,griddata);
    }

    cout << "Updating u and v...\n";
    // initialising the start time to 0. it gets overwritten if we use a uv file
    int starttime = 0;
    if(option == FROM_UV_FILE)
    {
        // we hack this together as so: the filename looks like uv_plotxxx.vtk, we want the xxx. so we find the t, find the ., and grab everyting between
        int startpos = B_filename.find('t');
        int endpos = B_filename.find('.');
        string number = B_filename.substr(B_filename.find('t')+1,B_filename.find('.')-B_filename.find('t')-1);
        starttime = atoi(number.c_str()); 
    }
    // initialise the time to the starttime
    double CurrentTime = starttime;
    int iterationcounter = 0;

    // initialising timers
    time_t then = time(NULL);
    time_t rawtime;
    time (&rawtime);
    struct tm * timeinfo;
#pragma omp parallel default(none) shared (u,v,ku,kv,ucvx, iterationcounter,ucvy, ucvz,ucvmag,cout, rawtime, starttime, timeinfo,CurrentTime, knotcurves,knotcurvesold,minimizerstate,griddata)
    {
        while(CurrentTime <= TTime)
        {
#pragma omp single
            {
                // in this section we do all the on the fly analysis. A few things happen


                // if we want to do the box resizing, it happens here
                if(BoxResizeFlag && ( abs(CurrentTime-BoxResizeTime) < (dtime/2) ))   
                {
                    cout << "resizingbox";
                    resizebox(u,v,ucvx,ucvy,ucvz,knotcurves,ku,kv,griddata);
                }

                // its useful to have an oppurtunity to print the knotcurve, without doing the velocity tracking, whihc doesnt work too well if we go more frequenclty
                // than a cycle
                if( ( CurrentTime > InitialSkipTime ) && ( fabs( CurrentTime - FrequentKnotplotPrintTime*round(CurrentTime/FrequentKnotplotPrintTime))<(dtime/2) ) )  
                {
                    crossgrad_calc(u,v,ucvx,ucvy,ucvz,ucvmag,griddata); //find Grad u cross Grad v
                    find_knot_properties(ucvx,ucvy,ucvz,ucvmag,u,knotcurves,CurrentTime,minimizerstate ,griddata);      //find knot curve and twist and writhe
                    print_knot(CurrentTime, knotcurves, griddata);
                }

                // run the curve tracing, and find the velocity of the one we previously stored, then print that previous one 
                if( ( CurrentTime > InitialSkipTime ) && ( fabs( CurrentTime - VelocityKnotplotPrintTime*round(CurrentTime/VelocityKnotplotPrintTime))<(dtime/2) ) )  
                {
                    crossgrad_calc(u,v,ucvx,ucvy,ucvz,ucvmag,griddata); //find Grad u cross Grad v

                    find_knot_properties(ucvx,ucvy,ucvz,ucvmag,u,knotcurves,CurrentTime,minimizerstate ,griddata);      //find knot curve and twist and writhe
                    if(!knotcurvesold.empty())
                    {
                        overlayknots(knotcurves,knotcurvesold, griddata);
                        find_knot_velocity(knotcurves,knotcurvesold,griddata,VelocityKnotplotPrintTime);
                        print_knot(CurrentTime - VelocityKnotplotPrintTime , knotcurvesold, griddata);
                    }
                    knotcurvesold = knotcurves;

                    // at this point, let people know how things are going
                    cout << "T = " << CurrentTime << endl;
                    time (&rawtime);
                    timeinfo = localtime (&rawtime);
                    cout << "current time \t" << asctime(timeinfo) << "\n";
                }

                // print the UV, and ucrossv data
                if(fmod(CurrentTime,UVPrintTime)<(dtime/2))  
                {
                    crossgrad_calc(u,v,ucvx,ucvy,ucvz,ucvmag,griddata); //find Grad u cross Grad v
                    print_uv(u,v,ucvx,ucvy,ucvz,ucvmag,CurrentTime,griddata);
                }
                //though its useful to have a double time, we want to be careful to avoid double round off accumulation in the timer
                iterationcounter++;
                CurrentTime  = starttime + ((double)(iterationcounter) * dtime);
            }
            uv_update(u,v,ku,kv,griddata);
        }
    }
    return 0;
}

/*************************Functions for knot initialisation*****************************/
double initialise_knot(vector<triangle>& knotsurface)
{
    double L;
    switch (option)
    {
        case FROM_SURFACE_FILE: L = init_from_surface_file(knotsurface);
                                break;

        default: L=0;
                 break;
    }

    return L;
}

double init_from_surface_file(vector<triangle>& knotsurface)
{
    string filename, buff;
    stringstream ss;
    double A = 0;   //total area
    int i=0;
    int j;
    double r10,r20,r21,s,xcoord,ycoord,zcoord;
    string temp;
    ifstream knotin;
    /*  For recording max and min input values*/
    double maxxin = 0;
    double maxyin = 0;
    double maxzin = 0;
    double minxin = 0;
    double minyin = 0;
    double minzin = 0;

    ss.clear();
    ss.str("");
    ss << knot_filename << ".stl";

    filename = ss.str();
    knotin.open(filename.c_str());
    if(knotin.good())
    {
        if(getline(knotin,buff)) temp = buff;
    }
    else cout << "Error reading file\n";
    while(knotin.good())   //read in points for knot
    {
        if(getline(knotin,buff))  //read in surface normal
        {
            ss.clear();
            ss.str("");
            ss << buff;
            ss >> temp;
            if(temp.compare("endsolid") == 0) break;
            knotsurface.push_back(triangle());
            ss >> temp >> knotsurface[i].normal[0] >> knotsurface[i].normal[1] >> knotsurface[i].normal[2];
        }

        if(getline(knotin,buff)) temp = buff;   //read in "outer loop"
        knotsurface[i].centre[0] = 0;
        knotsurface[i].centre[1] = 0;
        knotsurface[i].centre[2] = 0;
        for(j=0;j<3;j++)
        {
            if(getline(knotin,buff))  //read in vertices
            {
                ss.clear();
                ss.str("");
                ss << buff;
                ss >> temp >> xcoord >> ycoord >> zcoord;

                if(xcoord>maxxin) maxxin = xcoord;
                if(ycoord>maxyin) maxyin = ycoord;
                if(zcoord>maxzin) maxzin = zcoord;
                if(xcoord<minxin) minxin = xcoord;
                if(ycoord<minyin) minyin = ycoord;
                if(zcoord<minzin) minzin = zcoord;

                knotsurface[i].xvertex[j] = xcoord;
                knotsurface[i].yvertex[j] = ycoord;
                knotsurface[i].zvertex[j] = zcoord;
                knotsurface[i].centre[0] += knotsurface[i].xvertex[j]/3.0;
                knotsurface[i].centre[1] += knotsurface[i].yvertex[j]/3.0;
                knotsurface[i].centre[2] += knotsurface[i].zvertex[j]/3.0;
            }
        }
        //cout << i << " (" << knotsurface[i].centre[0] << ',' << knotsurface[i].centre[1] << ',' << knotsurface[i].centre[2] << ") , (" << knotsurface[i].normal[0] << ',' << knotsurface[i].normal[1] << ',' << knotsurface[i].normal[2] << ") \n";

        if(getline(knotin,buff)) temp = buff;   //read in "outer loop"
        if(getline(knotin,buff)) temp = buff;   //read in "outer loop"

        i++;
    }


    /* Work out space scaling for knot surface */
    double scale[3];
    double midpoint[3];
    double norm;
    scalefunction(scale,midpoint,maxxin,minxin,maxyin,minyin,maxzin,minzin);

    /*Rescale points and normals to fit grid properly*/
    for(i=0;i<knotsurface.size();i++)
    {
        for(j=0;j<3;j++)
        {
            knotsurface[i].xvertex[j] = scale[0]*(knotsurface[i].xvertex[j] - midpoint[0]);
            knotsurface[i].yvertex[j] = scale[1]*(knotsurface[i].yvertex[j] - midpoint[1]);
            knotsurface[i].zvertex[j] = scale[2]*(knotsurface[i].zvertex[j] - midpoint[2]);
            knotsurface[i].centre[j] = scale[j]*(knotsurface[i].centre[j] - midpoint[j]);
        }

        norm = sqrt(scale[1]*scale[1]*scale[2]*scale[2]*knotsurface[i].normal[0]*knotsurface[i].normal[0] +
                scale[0]*scale[0]*scale[2]*scale[2]*knotsurface[i].normal[1]*knotsurface[i].normal[1] +
                scale[0]*scale[0]*scale[1]*scale[1]*knotsurface[i].normal[2]*knotsurface[i].normal[2]);

        knotsurface[i].normal[0] *= scale[1]*scale[2]/norm;
        knotsurface[i].normal[1] *= scale[0]*scale[2]/norm;
        knotsurface[i].normal[2] *= scale[0]*scale[1]/norm;

        /*Check surface normal is correct
          p1x = knotsurface[i].xvertex[1] - knotsurface[i].xvertex[0];
          p1y = knotsurface[i].yvertex[1] - knotsurface[i].yvertex[0];
          p1z = knotsurface[i].zvertex[1] - knotsurface[i].zvertex[0];
          p2x = knotsurface[i].xvertex[2] - knotsurface[i].xvertex[0];
          p2y = knotsurface[i].yvertex[2] - knotsurface[i].yvertex[0];
          p2z = knotsurface[i].zvertex[2] - knotsurface[i].zvertex[0];
          nx = p1y*p2z - p2y*p1z;
          ny = p1z*p2x - p2z*p1x;
          nz = p1x*p2y - p2x*p1y;
          norm = sqrt(nx*nx+ny*ny+nz*nz);
          nx = nx/norm;
          ny = ny/norm;
          nz = nz/norm;
          cout << nx*knotsurface[i].normal[0] + ny*knotsurface[i].normal[1] + nz*knotsurface[i].normal[2] << '\n';
          */

        r10 = sqrt((knotsurface[i].xvertex[1]-knotsurface[i].xvertex[0])*(knotsurface[i].xvertex[1]-knotsurface[i].xvertex[0]) + (knotsurface[i].yvertex[1]-knotsurface[i].yvertex[0])*(knotsurface[i].yvertex[1]-knotsurface[i].yvertex[0]) + (knotsurface[i].zvertex[1]-knotsurface[i].zvertex[0])*(knotsurface[i].zvertex[1]-knotsurface[i].zvertex[0]));
        r20 = sqrt((knotsurface[i].xvertex[2]-knotsurface[i].xvertex[0])*(knotsurface[i].xvertex[2]-knotsurface[i].xvertex[0]) + (knotsurface[i].yvertex[2]-knotsurface[i].yvertex[0])*(knotsurface[i].yvertex[2]-knotsurface[i].yvertex[0]) + (knotsurface[i].zvertex[2]-knotsurface[i].zvertex[0])*(knotsurface[i].zvertex[2]-knotsurface[i].zvertex[0]));
        r21 = sqrt((knotsurface[i].xvertex[2]-knotsurface[i].xvertex[1])*(knotsurface[i].xvertex[2]-knotsurface[i].xvertex[1]) + (knotsurface[i].yvertex[2]-knotsurface[i].yvertex[1])*(knotsurface[i].yvertex[2]-knotsurface[i].yvertex[1]) + (knotsurface[i].zvertex[2]-knotsurface[i].zvertex[1])*(knotsurface[i].zvertex[2]-knotsurface[i].zvertex[1]));
        s = (r10+r20+r21)/2;
        knotsurface[i].area = sqrt(s*(s-r10)*(s-r20)*(s-r21));
        A += knotsurface[i].area;

        // apply any rotations and displacements  of the initial coniditions the user has specified
        for(j=0;j<3;j++) rotatedisplace(knotsurface[i].xvertex[j],knotsurface[i].yvertex[j],knotsurface[i].zvertex[j],initialthetarotation,initialxdisplacement,initialydisplacement,initialzdisplacement);
        rotatedisplace(knotsurface[i].normal[0],knotsurface[i].normal[1],knotsurface[i].normal[2],initialthetarotation,initialxdisplacement,initialydisplacement,initialzdisplacement);
        rotatedisplace(knotsurface[i].centre[0],knotsurface[i].centre[1],knotsurface[i].centre[2],initialthetarotation,initialxdisplacement,initialydisplacement,initialzdisplacement);
    }

    cout << "Input scaled by: " << scale[0] << ' ' << scale[1] << ' ' << scale[2] << " in x,y and z\n";

    return A;
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

/*************************Functions for B and Phi calcs*****************************/

void phi_calc(vector<double>&phi,vector<triangle>& knotsurface, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n,s;
    double rx,ry,rz,r;
    cout << "Calculating scalar potential...\n";
#pragma omp parallel default(none) shared (Nx,Ny,Nz,griddata, knotsurface, phi ) private ( i, j, k, n, s, rx, ry, rz , r)
    {
#pragma omp for
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)
                {
                    n = pt(i,j,k,griddata);
                    phi[n] = 0;
                    for(s=0;s<knotsurface.size();s++)
                    {
                        rx = knotsurface[s].centre[0]-x(i,griddata);
                        ry = knotsurface[s].centre[1]-y(j,griddata);
                        rz = knotsurface[s].centre[2]-z(k,griddata);
                        r = sqrt(rx*rx+ry*ry+rz*rz);
                        if(r>0) phi[n] += (rx*knotsurface[s].normal[0] + ry*knotsurface[s].normal[1] + rz*knotsurface[s].normal[2])*knotsurface[s].area/(2*r*r*r);
                    }
                    while(phi[n]>M_PI) phi[n] -= 2*M_PI;
                    while(phi[n]<-M_PI) phi[n] += 2*M_PI;
                }
            }
        }
    }
    cout << "Printing B and phi...\n";
    print_B_phi(phi,griddata);

}
void DualConePhiCalc(vector<double>&phi, const griddata& griddata)
{
    //first up, initialise the knotcurve
    knotcurve InitialisationCurve;

    stringstream ss;
    string buff,filename;

    ss.clear();
    ss.str("");
    ss << knot_filename << ".fe";

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

#pragma omp parallel default(none) shared(InitialisationCurve,phi,griddata,cout)
    {
        // now for the guts of the thing - construct the dual curve
        knotcurve ProjectedCurve;
        knotcurve DualCurve;
        //    vector<knotcurve> Curvetoprint;

        int intersectioncount = 0;

        int Nx = griddata.Nx;
        int Ny = griddata.Ny;
        int Nz = griddata.Nz;
#pragma omp for 
        for(int i=0;i<Nx;i++)
        {
            cout << i;
            time_t rawtime;
            struct tm * timeinfo;
            time (&rawtime);
            timeinfo = localtime (&rawtime);
            cout << "current time \t" << asctime(timeinfo) << "\n";
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

                    // get segment lengths
                    float topdeltas = 0;
                    for(int s=0;s<ProjectedCurve.knotcurve.size();s++)
                    {
                        int NP = ProjectedCurve.knotcurve.size();
                        double dx = (ProjectedCurve.knotcurve[incp(s,1,NP)].xcoord - ProjectedCurve.knotcurve[incp(s,0,NP)].xcoord);   //central diff as a is defined on the points
                        double dy = (ProjectedCurve.knotcurve[incp(s,1,NP)].ycoord - ProjectedCurve.knotcurve[incp(s,0,NP)].ycoord);
                        double dz = (ProjectedCurve.knotcurve[incp(s,1,NP)].zcoord - ProjectedCurve.knotcurve[incp(s,0,NP)].zcoord);
                        double deltas = sqrt(dx*dx+dy*dy+dz*dz);
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
                        int NP = ProjectedCurve.knotcurve.size();
                        // forward difference on the tangents
                        double dx = (ProjectedCurve.knotcurve[incp(s,1,NP)].xcoord - ProjectedCurve.knotcurve[incp(s,0,NP)].xcoord);   //central diff as a is defined on the points
                        double dy = (ProjectedCurve.knotcurve[incp(s,1,NP)].ycoord - ProjectedCurve.knotcurve[incp(s,0,NP)].ycoord);
                        double dz = (ProjectedCurve.knotcurve[incp(s,1,NP)].zcoord - ProjectedCurve.knotcurve[incp(s,0,NP)].zcoord);
                        double deltas = sqrt(dx*dx+dy*dy+dz*dz);
                        double tx = dx/(deltas);
                        double ty = dy/(deltas);
                        double tz = dz/(deltas);

                        // a central differencing scheme
                        //double tx = (ProjectedCurve.knotcurve[incp(s,1,NP)].xcoord - ProjectedCurve.knotcurve[incp(s,-1,NP)].xcoord)/(ProjectedCurve.knotcurve[incp(s,-1,NP)].length+ProjectedCurve.knotcurve[incp(s,0,NP)].length) ;
                        //double ty = (ProjectedCurve.knotcurve[incp(s,1,NP)].ycoord - ProjectedCurve.knotcurve[incp(s,-1,NP)].ycoord)/(ProjectedCurve.knotcurve[incp(s,-1,NP)].length+ProjectedCurve.knotcurve[incp(s,0,NP)].length);
                        // double tz = (ProjectedCurve.knotcurve[incp(s,1,NP)].zcoord - ProjectedCurve.knotcurve[incp(s,-1,NP)].zcoord)/(ProjectedCurve.knotcurve[incp(s,-1,NP)].length+ProjectedCurve.knotcurve[incp(s,0,NP)].length) ;

                        ProjectedCurve.knotcurve[s].tx = tx;
                        ProjectedCurve.knotcurve[s].ty = ty;
                        ProjectedCurve.knotcurve[s].tz = tz;

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

                        // data to output to the dual
                        DualCurve.knotcurve[s].bx = dndsstarx;
                        DualCurve.knotcurve[s].by = dndsstary;
                        DualCurve.knotcurve[s].bz = dndsstarz;

                        DualCurve.knotcurve[s].length = triadsign*geodesicdistance;
                        DualCurve.knotcurve[s].torsion = triadsign;

                        DualCurve.length += triadsign*geodesicdistance;
                    }


                    //            if(i==55 && j == 95)
                    //            {
                    //                Curvetoprint.push_back(ProjectedCurve);
                    //                print_knot(k, Curvetoprint, griddata,"ProjectedCurve");
                    //                Curvetoprint.clear();
                    //                Curvetoprint.push_back(DualCurve);
                    //                print_knot(k, Curvetoprint, griddata,"DualCurve");
                    //                Curvetoprint.clear();
                    //                Curvetoprint.push_back(InitialisationCurve);
                    //                print_knot(0, Curvetoprint, griddata,"InitialisationCurve");
                    //                Curvetoprint.clear();
                    //            }
                    // clean up, and set phi
                    int parity = intersectioncount%2;
                    double SolidAngle = DualCurve.length + 2*M_PI*parity;
                    while(SolidAngle>4*M_PI) SolidAngle-= 4*M_PI;
                    while(SolidAngle<0) SolidAngle += 4*M_PI;

                    int n = pt(i,j,k,griddata);
                    phi[n] = SolidAngle/2;
                }
            }
        }
    }
    cout << "Printing B and phi...\n";
    print_B_phi(phi,griddata);
    return;
}

void phi_calc_manual(vector<double>&phi, griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n;
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)
            {
                n = pt(i,j,k,griddata);
                phi[n] = 0;
                double theta = 0.5;
                phi[n] = atan2(y(j,griddata)-lambda,x(i,griddata)-lambda)- atan2(y(j,griddata),-sin(theta)*z(k,griddata) +cos(theta)*x(i,griddata));
                while(phi[n]>M_PI) phi[n] -= 2*M_PI;
                while(phi[n]<-M_PI) phi[n] += 2*M_PI;
            }
        }
    }
    cout << "Printing B and phi...\n";
    print_B_phi(phi,griddata);
}

/*************************Functions for FN dynamics*****************************/

void uv_initialise(vector<double>&phi, vector<double>&u, vector<double>&v, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int n;

    for(n=0; n<Nx*Ny*Nz; n++)
    {
        u[n] = (2*cos(phi[n]) - 0.4);
        v[n] = (sin(phi[n]) - 0.4);
    }
}

void crossgrad_calc( vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>&ucvmag,const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n,kup,kdown;
    double dxu,dyu,dzu,dxv,dyv,dzv;
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)   //Central difference
            {
                kup = gridinc(k,1,Nz,2);
                kdown = gridinc(k,-1,Nz,2);
                dxu = 0.5*(u[pt(gridinc(i,1,Nx,0),j,k,griddata)]-u[pt(gridinc(i,-1,Nx,0),j,k,griddata)])/h;
                dxv = 0.5*(v[pt(gridinc(i,1,Nx,0),j,k,griddata)]-v[pt(gridinc(i,-1,Nx,0),j,k,griddata)])/h;
                dyu = 0.5*(u[pt(i,gridinc(j,1,Ny,1),k,griddata)]-u[pt(i,gridinc(j,-1,Ny,1),k,griddata)])/h;
                dyv = 0.5*(v[pt(i,gridinc(j,1,Ny,1),k,griddata)]-v[pt(i,gridinc(j,-1,Ny,1),k,griddata)])/h;
                dzu = 0.5*(u[pt(i,j,kup,griddata)]-u[pt(i,j,kdown,griddata)])/h;
                dzv = 0.5*(v[pt(i,j,kup,griddata)]-v[pt(i,j,kdown,griddata)])/h;
                //          dxu =(-u[pt(gridinc(i,2,Nx,0),j,k,griddata)]+8*u[pt(gridinc(i,1,Nx,0),j,k,griddata)]-8*u[pt(gridinc(i,-1,Nx,0),j,k,griddata)]+u[pt(gridinc(i,-2,Nx,0),j,k,griddata)])/(12*h);
                //          dxv =(-v[pt(gridinc(i,2,Nx,0),j,k,griddata)]+8*v[pt(gridinc(i,1,Nx,0),j,k,griddata)]-8*v[pt(gridinc(i,-1,Nx,0),j,k,griddata)]+v[pt(gridinc(i,-2,Nx,0),j,k,griddata)])/(12*h);
                //          dyu =(-u[pt(gridinc(j,2,Ny,1),j,k,griddata)]+8*u[pt(gridinc(j,1,Ny,1),j,k,griddata)]-8*u[pt(gridinc(j,-1,Ny,1),j,k,griddata)]+u[pt(gridinc(j,-2,Ny,1),j,k,griddata)])/(12*h);
                //          dyv =(-v[pt(gridinc(j,2,Ny,1),j,k,griddata)]+8*v[pt(gridinc(j,1,Ny,1),j,k,griddata)]-8*v[pt(gridinc(j,-1,Ny,1),j,k,griddata)]+v[pt(gridinc(j,-2,Ny,1),j,k,griddata)])/(12*h);
                //          dzu =(-u[pt(gridinc(k,2,Nz,2),j,k,griddata)]+8*u[pt(gridinc(k,1,Nz,2),j,k,griddata)]-8*u[pt(gridinc(k,-1,Nz,2),j,k,griddata)]+u[pt(gridinc(k,-2,Nz,2),j,k,griddata)])/(12*h);
                //          dzv =(-v[pt(gridinc(k,2,Nz,2),j,k,griddata)]+8*v[pt(gridinc(k,1,Nz,2),j,k,griddata)]-8*v[pt(gridinc(k,-1,Nz,2),j,k,griddata)]+v[pt(gridinc(k,-2,Nz,2),j,k,griddata)])/(12*h);
                n = pt(i,j,k,griddata);
                ucvx[n] = dyu*dzv - dzu*dyv;
                ucvy[n] = dzu*dxv - dxu*dzv;    //Grad u cross Grad v
                ucvz[n] = dxu*dyv - dyu*dxv;
                ucvmag[n] = sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]);
            }
        }
    }
}

void find_knot_properties( vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz, vector<double>& ucvmag,vector<double>&u,vector<knotcurve>& knotcurves,double t, gsl_multimin_fminimizer* minimizerstate, const griddata& griddata)
{
    // first thing, clear the knotcurve object before we begin writing a new one 
    knotcurves.clear(); //empty vector with knot curve points

    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;

    // initialise the tricubic interpolator for ucvmag
    likely::TriCubicInterpolator interpolateducvmag(ucvmag, h, Nx,Ny,Nz);

    int c =0;
    bool knotexists = true;

    double   ucvmax = -1.0; // should always be +ve, so setting it to an initially -ve # means it always gets written to once.
    int n,i,j,k,imax,jmax,kmax;
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)   //Central difference
            {
                n = pt(i,j,k,griddata);
                if( ucvmag[n] > ucvmax)
                {
                    ucvmax = ucvmag[n];
                    imax = i;
                    jmax = j;
                    kmax=k;
                }
            }
        }
    }
    if(ucvmax<0.45) knotexists = false;
    else
    {
        knotexists = true; 
    }
    if(knotexists)
    {
        knotcurves.push_back(knotcurve() );
        knotcurves[c].knotcurve.push_back(knotpoint());
        knotcurves[c].knotcurve[0].xcoord=x(imax,griddata);
        knotcurves[c].knotcurve[0].ycoord=y(jmax,griddata);
        knotcurves[c].knotcurve[0].zcoord=z(kmax,griddata);

        int idwn,jdwn,kdwn, modidwn, modjdwn, modkdwn,m,pts,iinc,jinc,kinc,attempts;
        double ucvxs, ucvys, ucvzs, graducvx, graducvy, graducvz, prefactor, xd, yd ,zd, fx, fy, fz, xdiff, ydiff, zdiff;
        int s=1;
        bool finish=false;
        // we will discard the first few points from the knot, using this flag
        bool burnin=true;
        /*calculate local direction of grad u x grad v (the tangent to the knot curve) at point s-1, then move to point s by moving along tangent + unit confinement force*/
        while (finish==false)
        {

            /**Find nearest gridpoint**/
            idwn = (int) ((knotcurves[c].knotcurve[s-1].xcoord/h) - 0.5 + Nx/2.0);
            jdwn = (int) ((knotcurves[c].knotcurve[s-1].ycoord/h) - 0.5 + Ny/2.0);
            kdwn = (int) ((knotcurves[c].knotcurve[s-1].zcoord/h) - 0.5 + Nz/2.0);
            // idwn etc can be off the actual grid , into "ghost" grids around the real one. this is useful for knotcurve tracing over periodic boundaries
            // but we also need the corresponding real grid positions!
            modidwn = circularmod(idwn,Nx);
            modjdwn = circularmod(jdwn,Ny);
            modkdwn = circularmod(kdwn,Nz);
            if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
            if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
            pts=0;
            ucvxs=0;
            ucvys=0;
            ucvzs=0;
            /*curve to gridpoint down distance*/
            xd = (knotcurves[c].knotcurve[s-1].xcoord - x(idwn,griddata))/h;
            yd = (knotcurves[c].knotcurve[s-1].ycoord - y(jdwn,griddata))/h;
            zd = (knotcurves[c].knotcurve[s-1].zcoord - z(kdwn,griddata))/h;
            for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
            {
                /* Work out increments*/
                iinc = m%2;
                jinc = (m/2)%2;
                kinc = (m/4)%2;
                /*Loop over nearest points*/
                i = gridinc(modidwn, iinc, Nx,0);
                j = gridinc(modjdwn, jinc, Ny,1);
                k = gridinc(modkdwn,kinc, Nz,2);
                prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                /*interpolate grad u x grad v over nearest points*/
                ucvxs += prefactor*ucvx[pt(i,j,k,griddata)];
                ucvys += prefactor*ucvy[pt(i,j,k,griddata)];
                ucvzs += prefactor*ucvz[pt(i,j,k,griddata)];
            }
            double norm = sqrt(ucvxs*ucvxs + ucvys*ucvys + ucvzs*ucvzs);
            ucvxs = ucvxs/norm; //normalise
            ucvys = ucvys/norm; //normalise
            ucvzs = ucvzs/norm; //normalise

            // okay we have our first guess, move forward in this direction
            double testx = knotcurves[c].knotcurve[s-1].xcoord + 0.5*ucvxs*lambda/(2*M_PI);
            double testy = knotcurves[c].knotcurve[s-1].ycoord + 0.5*ucvys*lambda/(2*M_PI);
            double testz = knotcurves[c].knotcurve[s-1].zcoord + 0.5*ucvzs*lambda/(2*M_PI);

            // now get the grad at this point
            idwn = (int) ((testx/h) - 0.5 + Nx/2.0);
            jdwn = (int) ((testy/h) - 0.5 + Ny/2.0);
            kdwn = (int) ((testz/h) - 0.5 + Nz/2.0);
            modidwn = circularmod(idwn,Nx);
            modjdwn = circularmod(jdwn,Ny);
            modkdwn = circularmod(kdwn,Nz);
            // again, bear in mind these numbers can be into the "ghost" grids
            if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
            if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
            pts=0;
            graducvx=0;
            graducvy=0;
            graducvz=0;
            /*curve to gridpoint down distance*/
            xd = (testx - x(idwn,griddata))/h;
            yd = (testy - y(jdwn,griddata))/h;
            zd = (testz - z(kdwn,griddata))/h;
            for(m=0;m<8;m++)  //linear interpolation from 8 nearest neighbours
            {
                /* Work out increments*/
                iinc = m%2;
                jinc = (m/2)%2;
                kinc = (m/4)%2;
                /*Loop over nearest points*/
                i = gridinc(modidwn, iinc, Nx,0);
                j = gridinc(modjdwn, jinc, Ny,1);
                k = gridinc(modkdwn,kinc, Nz,2);
                prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);
                /*interpolate gradients of |grad u x grad v|*/
                graducvx += prefactor*(sqrt(ucvx[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvx[pt(gridinc(i,1,Nx,0),j,k,griddata)] + ucvy[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvy[pt(gridinc(i,1,Nx,0),j,k,griddata)] + ucvz[pt(gridinc(i,1,Nx,0),j,k,griddata)]*ucvz[pt(gridinc(i,1,Nx,0),j,k,griddata)]) - sqrt(ucvx[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvx[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + ucvy[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvy[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + ucvz[pt(gridinc(i,-1,Nx,0),j,k,griddata)]*ucvz[pt(gridinc(i,-1,Nx,0),j,k,griddata)]))/(2*h);
                graducvy += prefactor*(sqrt(ucvx[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvx[pt(i,gridinc(j,1,Ny,1),k,griddata)] + ucvy[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvy[pt(i,gridinc(j,1,Ny,1),k,griddata)] + ucvz[pt(i,gridinc(j,1,Ny,1),k,griddata)]*ucvz[pt(i,gridinc(j,1,Ny,1),k,griddata)]) - sqrt(ucvx[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvx[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + ucvy[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvy[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + ucvz[pt(i,gridinc(j,-1,Ny,1),k,griddata)]*ucvz[pt(i,gridinc(j,-1,Ny,1),k,griddata)]))/(2*h);
                graducvz += prefactor*(sqrt(ucvx[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvx[pt(i,j,gridinc(k,1,Nz,2),griddata)] + ucvy[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvy[pt(i,j,gridinc(k,1,Nz,2),griddata)] + ucvz[pt(i,j,gridinc(k,1,Nz,2),griddata)]*ucvz[pt(i,j,gridinc(k,1,Nz,2),griddata)]) - sqrt(ucvx[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvx[pt(i,j,gridinc(k,-1,Nz,2),griddata)] + ucvy[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvy[pt(i,j,gridinc(k,-1,Nz,2),griddata)] + ucvz[pt(i,j,gridinc(k,-1,Nz,2),griddata)]*ucvz[pt(i,j,gridinc(k,-1,Nz,2),griddata)]))/(2*h);

            }
            knotcurves[c].knotcurve.push_back(knotpoint());
            // one of the vectors in the plane we wish to perfrom our minimisation in
            fx = (graducvx - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvxs); 
            fy = (graducvy - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvys);
            fz = (graducvz - (graducvx*ucvxs + graducvy*ucvys + graducvz*ucvzs)*ucvzs);
            norm = sqrt(fx*fx + fy*fy + fz*fz);
            fx = fx/norm;
            fy = fy/norm;
            fz = fz/norm;

            // okay we have our direction to perfrom the line minimisation in
            // the point
            gsl_vector* v = gsl_vector_alloc (3);
            gsl_vector_set (v, 0, testx);
            gsl_vector_set (v, 1, testy);
            gsl_vector_set (v, 2, testz);
            // one vector in the plane we with to minimize in
            gsl_vector* f = gsl_vector_alloc (3);
            gsl_vector_set (f, 0, fx);
            gsl_vector_set (f, 1, fy);
            gsl_vector_set (f, 2, fz);
            // the ucv vector
            gsl_vector* ucv = gsl_vector_alloc (3);
            gsl_vector_set (ucv, 0, ucvxs);
            gsl_vector_set (ucv, 1, ucvys);
            gsl_vector_set (ucv, 2, ucvzs);
            // take a cross product to get the other vector in the plane 
            gsl_vector* b = gsl_vector_alloc (3);
            cross_product(f,ucv,b); 
            // initial conditions
            gsl_vector* minimum = gsl_vector_alloc (2);
            gsl_vector_set (minimum, 0, 0);
            gsl_vector_set (minimum, 1, 0);
            struct parameters params; struct parameters* pparams = &params;
            pparams->ucvmag=&interpolateducvmag;
            pparams->v = v; pparams->f = f;pparams->b=b;
            pparams->mygriddata = griddata;
            // some initial values
            gsl_multimin_function F;
            F.n=2;
            F.f = &my_f;
            F.params = (void*) pparams;
            gsl_vector* stepsize = gsl_vector_alloc (2);
            gsl_vector_set (stepsize, 0, lambda/(8*M_PI));
            gsl_vector_set (stepsize, 1, lambda/(8*M_PI));
            gsl_multimin_fminimizer_set (minimizerstate, &F, minimum, stepsize);

            int iter=0;
            int status =0;
            double minimizersize=0;
            do
            {
                iter++;
                status = gsl_multimin_fminimizer_iterate(minimizerstate);

                if (status) 
                    break;

                minimizersize = gsl_multimin_fminimizer_size (minimizerstate);
                status = gsl_multimin_test_size (minimizersize, 1e-2);

            }
            while (status == GSL_CONTINUE && iter < 500);


            gsl_vector_scale(f,gsl_vector_get(minimizerstate->x, 0));
            gsl_vector_scale(b,gsl_vector_get(minimizerstate->x, 1));
            gsl_vector_add(f,b);
            gsl_vector_add(v,f);
            knotcurves[c].knotcurve[s].xcoord = gsl_vector_get(v, 0);
            knotcurves[c].knotcurve[s].ycoord= gsl_vector_get(v, 1);
            knotcurves[c].knotcurve[s].zcoord= gsl_vector_get(v, 2);

            gsl_vector_free(v);
            gsl_vector_free(f);
            gsl_vector_free(b);
            gsl_vector_free(ucv);
            gsl_vector_free(stepsize);

            xdiff = knotcurves[c].knotcurve[0].xcoord - knotcurves[c].knotcurve[s].xcoord;     //distance from start/end point
            ydiff = knotcurves[c].knotcurve[0].ycoord - knotcurves[c].knotcurve[s].ycoord;
            zdiff = knotcurves[c].knotcurve[0].zcoord - knotcurves[c].knotcurve[s].zcoord;
            if(sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff) <3*h  && s > 10) finish = true;
            if(s>50000) finish = true;

            // okay, we just added a point in position s in the vector
            // if we have a few points in the vector, discard the first few and restart the whole thing - burn it in
            int newstartingposition =5;
            if(s==newstartingposition && burnin)
            {
                knotcurves[c].knotcurve.erase(knotcurves[c].knotcurve.begin(),knotcurves[c].knotcurve.begin()+newstartingposition);
                s =0;
                burnin =false;
            }

            s++;
        }

        int NP = knotcurves[c].knotcurve.size();  //store number of points in knot curve


        /*******Vertex averaging*********/

        double totlength, dl, dx,dy,dz;
        for(i=0;i<3;i++)   //repeat a couple of times because of end point
        {
            totlength=0;
            for(s=0; s<NP; s++)   //Work out total length of curve
            {
                dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
                dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
                dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
                totlength += sqrt(dx*dx + dy*dy + dz*dz);
            }
            dl = totlength/NP;
            for(s=0; s<NP; s++)    //Move points to have spacing dl
            {
                dx = knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[s].xcoord;
                dy = knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[s].ycoord;
                dz = knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[s].zcoord;
                double norm = sqrt(dx*dx + dy*dy + dz*dz);
                knotcurves[c].knotcurve[incp(s,1,NP)].xcoord = knotcurves[c].knotcurve[s].xcoord + dl*dx/norm;
                knotcurves[c].knotcurve[incp(s,1,NP)].ycoord = knotcurves[c].knotcurve[s].ycoord + dl*dy/norm;
                knotcurves[c].knotcurve[incp(s,1,NP)].zcoord = knotcurves[c].knotcurve[s].zcoord + dl*dz/norm;
            }
        }

        /*************Curve Smoothing*******************/
        vector<double> coord(NP);
        gsl_fft_real_wavetable * real;
        gsl_fft_halfcomplex_wavetable * hc;
        gsl_fft_real_workspace * work;
        work = gsl_fft_real_workspace_alloc (NP);
        real = gsl_fft_real_wavetable_alloc (NP);
        hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
        for(j=1; j<4; j++)
        {
            switch(j)
            {
                case 1 :
                    for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].xcoord ; break;
                case 2 :
                    for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ycoord ; break;
                case 3 :
                    for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].zcoord ; break;
            }
            double* data = coord.data();
            // take the fft
            gsl_fft_real_transform (data, 1, NP, real, work);
            // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
            // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
            // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
            // plug in above. im doing a rough length calc below, this might be overkill.
            // at the moment its just a hard filter, we can choose others though.
            // compute a rough length to set scale
            double filter;
            const double cutoff = 2*M_PI*(totlength/(6*lambda));
            for (i = 0; i < NP; ++i)
            {
                filter = 1/sqrt(1+pow((i/cutoff),8));
                data[i] *= filter;
            };
            // transform back
            gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
            switch(j)
            {
                case 1 :
                    for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].xcoord = coord[i] ; break;
                case 2 :
                    for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ycoord = coord[i] ; break;
                case 3 :
                    for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].zcoord = coord[i] ; break;
            }
        }



        /******************Interpolate direction of grad u for twist calc*******/
        /**Find nearest gridpoint**/
        double dxu, dyu, dzu, dxup, dyup, dzup;
        for(s=0; s<NP; s++)
        {
            idwn = (int) ((knotcurves[c].knotcurve[s].xcoord/h) - 0.5 + Nx/2.0);
            jdwn = (int) ((knotcurves[c].knotcurve[s].ycoord/h) - 0.5 + Ny/2.0);
            kdwn = (int) ((knotcurves[c].knotcurve[s].zcoord/h) - 0.5 + Nz/2.0);
            modidwn = circularmod(idwn,Nx);
            modjdwn = circularmod(jdwn,Ny);
            modkdwn = circularmod(kdwn,Nz);
            if((BoundaryType==ALLREFLECTING) && (idwn<0 || jdwn<0 || kdwn<0 || idwn > Nx-1 || jdwn > Ny-1 || kdwn > Nz-1)) break;
            if((BoundaryType==ZPERIODIC) && (idwn<0 || jdwn<0 || idwn > Nx-1 || jdwn > Ny-1 )) break;
            dxu=0;
            dyu=0;
            dzu=0;
            /*curve to gridpoint down distance*/
            xd = (knotcurves[c].knotcurve[s].xcoord - x(idwn,griddata))/h;
            yd = (knotcurves[c].knotcurve[s].ycoord - y(jdwn,griddata))/h;
            zd = (knotcurves[c].knotcurve[s].zcoord - z(kdwn,griddata))/h;
            for(m=0;m<8;m++)  //linear interpolation of 8 NNs
            {
                /* Work out increments*/
                iinc = m%2;
                jinc = (m/2)%2;
                kinc = (m/4)%2;
                /*Loop over nearest points*/
                i = gridinc(modidwn, iinc, Nx,0);
                j = gridinc(modjdwn, jinc, Ny,1);
                k = gridinc(modkdwn,kinc, Nz,2);
                prefactor = (1-iinc + pow(-1,1+iinc)*xd)*(1-jinc + pow(-1,1+jinc)*yd)*(1-kinc + pow(-1,1+kinc)*zd);   //terms of the form (1-xd)(1-yd)zd etc. (interpolation coefficient)
                /*interpolate grad u over nearest points*/
                dxu += prefactor*0.5*(u[pt(gridinc(i,1,Nx,0),j,k,griddata)] -  u[pt(gridinc(i,-1,Nx,0),j,k,griddata)])/h;  //central diff
                dyu += prefactor*0.5*(u[pt(i,gridinc(j,1,Ny,1),k,griddata)] -  u[pt(i,gridinc(j,-1,Ny,1),k,griddata)])/h;
                dzu += prefactor*0.5*(u[pt(i,j,gridinc(k,1,Nz,2),griddata)] -  u[pt(i,j,gridinc(k,-1,Nz,2),griddata)])/h;
            }
            //project du onto perp of tangent direction first
            dx = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].xcoord);   //central diff as a is defined on the points
            dy = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,-1,NP)].ycoord);
            dz = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,-1,NP)].zcoord);
            dxup = dxu - (dxu*dx + dyu*dy + dzu*dz)*dx/(dx*dx+dy*dy+dz*dz);               //Grad u_j * (delta_ij - t_i t_j)
            dyup = dyu - (dxu*dx + dyu*dy + dzu*dz)*dy/(dx*dx+dy*dy+dz*dz);
            dzup = dzu - (dxu*dx + dyu*dy + dzu*dz)*dz/(dx*dx+dy*dy+dz*dz);
            /*Vector a is the normalised gradient of u, should point in direction of max u perp to t*/
            double norm = sqrt(dxup*dxup+dyup*dyup+dzup*dzup);
            knotcurves[c].knotcurve[s].ax = dxup/norm;
            knotcurves[c].knotcurve[s].ay = dyup/norm;
            knotcurves[c].knotcurve[s].az = dzup/norm;
        }

        for(j=1; j<4; j++)
        {
            switch(j)
            {
                case 1 :
                    for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ax ; break;
                case 2 :
                    for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].ay ; break;
                case 3 :
                    for(i=0; i<NP; i++) coord[i] =  knotcurves[c].knotcurve[i].az ; break;
            }
            double* data = coord.data();
            // take the fft
            gsl_fft_real_transform (data, 1, NP, real, work);
            // 21/11/2016: make our low pass filter. To apply our filter. we should sample frequencies fn = n/Delta N , n = -N/2 ... N/2
            // this is discretizing the nyquist interval, with extreme frequency ~1/2Delta.
            // to cut out the frequencies of grid fluctuation size and larger we need a lengthscale Delta to
            // plug in above. im doing a rough length calc below, this might be overkill.
            // at the moment its just a hard filter, we can choose others though.
            // compute a rough length to set scale
            double filter;
            const double cutoff = 2*M_PI*(totlength/(1*lambda));
            for (i = 0; i < NP; ++i)
            {
                filter = 1/sqrt(1+pow((i/cutoff),8));
                data[i] *= filter;
            };
            // transform back
            gsl_fft_halfcomplex_inverse (data, 1, NP, hc, work);
            switch(j)
            {
                case 1 :
                    for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ax= coord[i] ; break;
                case 2 :
                    for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].ay= coord[i] ; break;
                case 3 :
                    for(i=0; i<NP; i++)  knotcurves[c].knotcurve[i].az = coord[i] ; break;
            }
        }
        gsl_fft_real_wavetable_free (real);
        gsl_fft_halfcomplex_wavetable_free (hc);
        gsl_fft_real_workspace_free (work);


        // CURVE GEOMETRY - get curvatures, torsions, frennet serret frame


        NP = knotcurves[c].knotcurve.size();
        for(s=0; s<NP; s++)   
        {
            // forward difference on the tangents
            double dx = (knotcurves[c].knotcurve[incp(s,1,NP)].xcoord - knotcurves[c].knotcurve[incp(s,0,NP)].xcoord);   
            double dy = (knotcurves[c].knotcurve[incp(s,1,NP)].ycoord - knotcurves[c].knotcurve[incp(s,0,NP)].ycoord);
            double dz = (knotcurves[c].knotcurve[incp(s,1,NP)].zcoord - knotcurves[c].knotcurve[incp(s,0,NP)].zcoord);
            double deltas = sqrt(dx*dx+dy*dy+dz*dz);
            knotcurves[c].knotcurve[s].tx = dx/(deltas);
            knotcurves[c].knotcurve[s].ty = dy/(deltas);
            knotcurves[c].knotcurve[s].tz = dz/(deltas);
            knotcurves[c].knotcurve[s].length = deltas;
            knotcurves[c].length +=deltas;
        }
        for(s=0; s<NP; s++)  
        {
            // backwards diff for the normals, amounting to a central diff overall
            double nx = 2.0*(knotcurves[c].knotcurve[s].tx-knotcurves[c].knotcurve[incp(s,-1,NP)].tx)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
            double ny = 2.0*(knotcurves[c].knotcurve[s].ty-knotcurves[c].knotcurve[incp(s,-1,NP)].ty)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
            double nz = 2.0*(knotcurves[c].knotcurve[s].tz-knotcurves[c].knotcurve[incp(s,-1,NP)].tz)/(knotcurves[c].knotcurve[s].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
            double curvature = sqrt(nx*nx+ny*ny+nz*nz);
            nx /=curvature;
            ny /=curvature;
            nz /=curvature;
            double tx = knotcurves[c].knotcurve[s].tx ;
            double ty =  knotcurves[c].knotcurve[s].ty ;
            double tz = knotcurves[c].knotcurve[s].tz ;
            double bx = ty*nz - tz*ny;
            double by = tz*nx - tx*nz;
            double bz = tx*ny - ty*nx;
            knotcurves[c].knotcurve[s].nx = nx ;
            knotcurves[c].knotcurve[s].ny = ny ;
            knotcurves[c].knotcurve[s].nz = nz ;
            knotcurves[c].knotcurve[s].bx = bx ;
            knotcurves[c].knotcurve[s].by = by ;
            knotcurves[c].knotcurve[s].bz = bz ;
            knotcurves[c].knotcurve[s].curvature = curvature ;
        }
        // torsions with a central difference
        for(s=0; s<NP; s++)   
        {
            double bx = knotcurves[c].knotcurve[s].bx;
            double by =  knotcurves[c].knotcurve[s].by;
            double bz = knotcurves[c].knotcurve[s].bz;

            double dnxds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].nx-knotcurves[c].knotcurve[incp(s,-1,NP)].nx)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
            double dnyds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].ny-knotcurves[c].knotcurve[incp(s,-1,NP)].ny)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);
            double dnzds = 2.0*(knotcurves[c].knotcurve[incp(s,1,NP)].nz-knotcurves[c].knotcurve[incp(s,-1,NP)].nz)/(knotcurves[c].knotcurve[incp(s,1,NP)].length+knotcurves[c].knotcurve[incp(s,-1,NP)].length);

            double torsion = bx*dnxds+by*dnyds+bz*dnzds;
            knotcurves[c].knotcurve[s].torsion = torsion ;
        }


        // RIBBON TWIST AND WRITHE

        for(s=0; s<NP; s++)   
        {

            // twist of this segment
            double ds = knotcurves[c].knotcurve[s].length;
            double dxds = knotcurves[c].knotcurve[s].tx;
            double dyds = knotcurves[c].knotcurve[s].ty;
            double dzds = knotcurves[c].knotcurve[s].tz;
            double bx = (knotcurves[c].knotcurve[incp(s,1,NP)].ax - knotcurves[c].knotcurve[s].ax)/ds;
            double by = (knotcurves[c].knotcurve[incp(s,1,NP)].ay - knotcurves[c].knotcurve[s].ay)/ds;
            double bz = (knotcurves[c].knotcurve[incp(s,1,NP)].az - knotcurves[c].knotcurve[s].az)/ds;
            knotcurves[c].knotcurve[s].twist = (dxds*(knotcurves[c].knotcurve[s].ay*bz - knotcurves[c].knotcurve[s].az*by) + dyds*(knotcurves[c].knotcurve[s].az*bx - knotcurves[c].knotcurve[s].ax*bz) + dzds*(knotcurves[c].knotcurve[s].ax*by - knotcurves[c].knotcurve[s].ay*bx))/(2*M_PI*sqrt(dxds*dxds + dyds*dyds + dzds*dzds));

            // "writhe" of this segment. writhe is nonlocal, this is the thing in the integrand over s
            knotcurves[c].knotcurve[s].writhe = 0;
            for(m=0; m<NP; m++)
            {
                if(s != m)
                {
                    xdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].xcoord + knotcurves[c].knotcurve[s].xcoord - knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord);   //interpolate, consistent with fwd diff
                    ydiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].ycoord + knotcurves[c].knotcurve[s].ycoord - knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord);
                    zdiff = 0.5*(knotcurves[c].knotcurve[incp(s,1,NP)].zcoord + knotcurves[c].knotcurve[s].zcoord - knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord);
                    double dxdm = (knotcurves[c].knotcurve[incp(m,1,NP)].xcoord - knotcurves[c].knotcurve[m].xcoord)/(ds);
                    double dydm = (knotcurves[c].knotcurve[incp(m,1,NP)].ycoord - knotcurves[c].knotcurve[m].ycoord)/(ds);
                    double dzdm = (knotcurves[c].knotcurve[incp(m,1,NP)].zcoord - knotcurves[c].knotcurve[m].zcoord)/(ds);
                    knotcurves[c].knotcurve[s].writhe += ds*(xdiff*(dyds*dzdm - dzds*dydm) + ydiff*(dzds*dxdm - dxds*dzdm) + zdiff*(dxds*dydm - dyds*dxdm))/(4*M_PI*(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)*sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff));
                }
            }

            //Add on writhe, twist 
            knotcurves[c].writhe += knotcurves[c].knotcurve[s].writhe*ds;
            knotcurves[c].twist  += knotcurves[c].knotcurve[s].twist*ds;
            // while we are computing the global quantites, get the average position too
            knotcurves[c].xavgpos += knotcurves[c].knotcurve[s].xcoord/NP;
            knotcurves[c].yavgpos += knotcurves[c].knotcurve[s].ycoord/NP;
            knotcurves[c].zavgpos += knotcurves[c].knotcurve[s].zcoord/NP;
        }

        // the ghost grid has been useful for painlessly computing all the above quantities, without worrying about the periodic bc's
        // but for storage and display, we should put it all in the box
        double xmax = x(griddata.Nx -1 ,griddata);
        double xmin = x(0,griddata);
        double deltax = griddata.Nx * h;
        double ymax = y(griddata.Ny -1,griddata);
        double ymin = y(0,griddata);
        double deltay = griddata.Ny * h;
        double zmax = z(griddata.Nz -1 ,griddata);
        double zmin = z(0,griddata);
        double deltaz = griddata.Nz * h;

        for(s=0; s<NP; s++)
        {
            knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord;
            knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord;
            knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord;
            if(knotcurves[c].knotcurve[s].xcoord > xmax) {
                knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord-deltax;
            } ;
            if(knotcurves[c].knotcurve[s].xcoord < xmin) {
                knotcurves[c].knotcurve[s].modxcoord = knotcurves[c].knotcurve[s].xcoord+deltax;
            };
            if(knotcurves[c].knotcurve[s].ycoord > ymax) {
                knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord-deltay;
            };
            if(knotcurves[c].knotcurve[s].ycoord < ymin) {
                knotcurves[c].knotcurve[s].modycoord = knotcurves[c].knotcurve[s].ycoord+deltay;
            };
            if(knotcurves[c].knotcurve[s].zcoord > zmax) {
                knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord-deltaz;
            };
            if(knotcurves[c].knotcurve[s].zcoord < zmin) 
            {
                knotcurves[c].knotcurve[s].modzcoord = knotcurves[c].knotcurve[s].zcoord+deltaz;
            };
        }
    }
}

void find_knot_velocity(const vector<knotcurve>& knotcurves,vector<knotcurve>& knotcurvesold,const griddata& griddata,const double deltatime)
{
    for(int c=0;c<knotcurvesold.size();c++)
    {

        int NP = knotcurves[c].knotcurve.size();
        int NPold = knotcurvesold[c].knotcurve.size();

        for(int s = 0; s< knotcurvesold[c].knotcurve.size(); s++)
        {
            double IntersectionFraction =-1;
            std::vector<double> IntersectionPoint(3);
            std::vector<double> ClosestIntersection(3);
            double closestdistancesquare = knotcurvesold[c].length;
            for(int t = 0 ; t<knotcurves[c].knotcurve.size();t++)
            {
                int intersection = 0;
                intersection = intersect3D_SegmentPlane( knotcurves[c].knotcurve[t%NP], knotcurves[c].knotcurve[(t+1)%NP], knotcurvesold[c].knotcurve[s%NPold], knotcurvesold[c].knotcurve[(s+1)%NPold], IntersectionFraction, IntersectionPoint );
                if(intersection ==1)
                { 
                    double intersectiondistancesquare = (IntersectionPoint[0] - knotcurvesold[c].knotcurve[s].xcoord )*(IntersectionPoint[0] - knotcurvesold[c].knotcurve[s].xcoord )+ (IntersectionPoint[1] - knotcurvesold[c].knotcurve[s].ycoord )*(IntersectionPoint[1] - knotcurvesold[c].knotcurve[s].ycoord )+ (IntersectionPoint[2] - knotcurvesold[c].knotcurve[s].zcoord )*(IntersectionPoint[2] - knotcurvesold[c].knotcurve[s].zcoord );
                    if(intersectiondistancesquare < closestdistancesquare)
                    {
                        closestdistancesquare = intersectiondistancesquare;
                        ClosestIntersection[0] = IntersectionPoint[0];
                        ClosestIntersection[1] = IntersectionPoint[1];
                        ClosestIntersection[2] = IntersectionPoint[2];
                    }

                }
            }
            // work out velocity and twist rate
            knotcurvesold[c].knotcurve[s].vx = (ClosestIntersection[0] - knotcurvesold[c].knotcurve[s].xcoord )/ deltatime;
            knotcurvesold[c].knotcurve[s].vy = (ClosestIntersection[1] - knotcurvesold[c].knotcurve[s].ycoord )/ deltatime;
            knotcurvesold[c].knotcurve[s].vz = (ClosestIntersection[2] - knotcurvesold[c].knotcurve[s].zcoord )/ deltatime;
            // for convenience, lets also output the decomposition into normal and binormal
            double vdotn = knotcurvesold[c].knotcurve[s].nx*knotcurvesold[c].knotcurve[s].vx+knotcurvesold[c].knotcurve[s].ny*knotcurvesold[c].knotcurve[s].vy+knotcurvesold[c].knotcurve[s].nz*knotcurvesold[c].knotcurve[s].vz;
            double vdotb = knotcurvesold[c].knotcurve[s].bx*knotcurvesold[c].knotcurve[s].vx+knotcurvesold[c].knotcurve[s].by*knotcurvesold[c].knotcurve[s].vy+knotcurvesold[c].knotcurve[s].bz*knotcurvesold[c].knotcurve[s].vz;

            knotcurvesold[c].knotcurve[s].vdotnx = vdotn * knotcurvesold[c].knotcurve[s].nx ;
            knotcurvesold[c].knotcurve[s].vdotny = vdotn * knotcurvesold[c].knotcurve[s].ny ;
            knotcurvesold[c].knotcurve[s].vdotnz = vdotn * knotcurvesold[c].knotcurve[s].nz ;
            knotcurvesold[c].knotcurve[s].vdotbx = vdotb * knotcurvesold[c].knotcurve[s].bx ;
            knotcurvesold[c].knotcurve[s].vdotby = vdotb * knotcurvesold[c].knotcurve[s].by ;
            knotcurvesold[c].knotcurve[s].vdotbz = vdotb * knotcurvesold[c].knotcurve[s].bz ;
        }
    }
}
void uv_update(vector<double>&u, vector<double>&v,  vector<double>&ku, vector<double>&kv,const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,l,n,kup,kdown,iup,idown,jup,jdown;
    double D2u;
    const int arraysize = Nx*Ny*Nz;
    // first loop. get k1, store (in testun] and testv[n], the value u[n]+h/2k1)
#pragma omp for 
    for(i=0;i<Nx;i++)
    {
        for(j=0; j<Ny; j++)
        {
            for(k=0; k<Nz; k++)   //Central difference
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k,griddata);
                    kup = gridinc(k,1,Nz,2);
                    kdown = gridinc(k,-1,Nz,2);
                    D2u = oneoverhsq*(u[pt(gridinc(i,1,Nx,0),j,k,griddata)] + u[pt(gridinc(i,-1,Nx,0),j,k,griddata)] + u[pt(i,gridinc(j,1,Ny,1),k,griddata)] + u[pt(i,gridinc(j,-1,Ny,1),k,griddata)] + u[pt(i,j,kup,griddata)] + u[pt(i,j,kdown,griddata)] - 6.0*u[n]);
                    ku[n] = oneoverepsilon*(u[n] - (ONETHIRD*u[n])*(u[n]*u[n]) - v[n]) + D2u;
                    kv[n] = epsilon*(u[n] + beta - gam*v[n]);
                }
            }
        }
    }
    // 2nd and 3rd loops
    double inc ;   
    for(l=1;l<=3;l++)  //u and v update for each fractional time step
    {
        switch (l)
        {
            case 1:
                {
                    inc=0.5;   //add k1 to uv and add to total k
                }
                break;

            case 2:
                {
                    inc=0.5 ;   //add k1 to uv and add to total k
                }
                break;
            case 3:
                {
                    inc=1 ;   //add k1 to uv and add to total k
                }
                break;
        }   
#pragma omp for 
        for(i=0;i<Nx;i++)
        {
            for(j=0; j<Ny; j++)
            {
                for(k=0; k<Nz; k++)   //Central difference
                {
                    n = pt(i,j,k,griddata);

                    iup = pt(gridinc(i,1,Nx,0),j,k,griddata);
                    idown =pt(gridinc(i,-1,Nx,0),j,k,griddata);
                    jup = pt(i,gridinc(j,1,Ny,1),k,griddata);
                    jdown =pt(i,gridinc(j,-1,Ny,1),k,griddata);
                    kup = pt(i,j,gridinc(k,1,Nz,2),griddata);
                    kdown = pt(i,j,gridinc(k,-1,Nz,2),griddata);
                    double currentu = u[n] + dtime*inc*ku[(l-1)*arraysize+n];
                    double currentv = v[n] + dtime*inc*kv[(l-1)*arraysize+n];

                    D2u = oneoverhsq*((u[iup]+dtime*inc*ku[(l-1)*arraysize+iup]) + (u[idown]+dtime*inc*ku[(l-1)*arraysize+idown]) +(u[jup]+dtime*inc*ku[(l-1)*arraysize+jup]) +(u[jdown]+dtime*inc*ku[(l-1)*arraysize+jdown]) + (u[kup]+dtime*inc*ku[(l-1)*arraysize+kup]) + (u[kdown]+dtime*inc*ku[(l-1)*arraysize+kdown])- 6.0*(currentu));


                    ku[arraysize*l+n] = oneoverepsilon*(currentu - (ONETHIRD*currentu)*(currentu*currentu) - currentv) + D2u;
                    kv[arraysize*l+n] = epsilon*(currentu + beta - gam*currentv);
                }
            }
        }
    }
#pragma omp for 
    for(n=0;n<Nx*Ny*Nz;n++)
    {

        u[n] = u[n] + dtime*sixth*(ku[n]+2*ku[arraysize+n]+2*ku[2*arraysize+n]+ku[3*arraysize+n]);
        v[n] = v[n] + dtime*sixth*(kv[n]+2*kv[arraysize+n]+2*kv[2*arraysize+n]+kv[3*arraysize+n]);
    }

}

/*************************File reading and writing*****************************/

void print_marked( vector<int>&marked,int shelllabel, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n;
    stringstream ss;
    ss << shelllabel<< "marked.vtk";
    ofstream uvout (ss.str().c_str());

    uvout << "# vtk DataFile Version 3.0\nUV fields\nASCII\nDATASET STRUCTURED_POINTS\n";
    uvout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    uvout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    uvout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    uvout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    uvout << "SCALARS marked float\nLOOKUP_TABLE default\n";


    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                uvout << marked[n] << '\n';
            }
        }
    }
    uvout.close();
}
void print_uv( vector<double>&u, vector<double>&v, vector<double>&ucvx, vector<double>&ucvy, vector<double>&ucvz,vector<double>&ucvmag, double t, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n;
    stringstream ss;
    ss << "uv_plot" << t << ".vtk";
    ofstream uvout (ss.str().c_str(),std::ios::binary | std::ios::out);

    uvout << "# vtk DataFile Version 3.0\nUV fields\nBINARY\nDATASET STRUCTURED_POINTS\n";
    uvout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    uvout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    uvout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    uvout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    uvout << "SCALARS u float\nLOOKUP_TABLE default\n";


    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                float val =  FloatSwap(u[n]);
                uvout.write((char*) &val, sizeof(float));
            }
        }
    }

    uvout << "\n" << "SCALARS v float\nLOOKUP_TABLE default\n";


    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                float val =  FloatSwap(v[n]);
                uvout.write( (char*) &val, sizeof(float));
            }
        }
    }

    uvout << "\n" << "SCALARS ucrossv float\nLOOKUP_TABLE default\n";

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                //float val = FloatSwap(sqrt(ucvx[n]*ucvx[n] + ucvy[n]*ucvy[n] + ucvz[n]*ucvz[n]));
                float val = FloatSwap(ucvmag[n]);
                uvout.write( (char*) &val, sizeof(float));
            }
        }
    }

    uvout.close();
}

void print_B_phi( vector<double>&phi, const griddata& griddata)
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


void print_knot( double t, vector<knotcurve>& knotcurves,const griddata& griddata, string seriesname)
{
    for( int c=0; c < (knotcurves.size()) ; c++)
    {

        /***Write values to file*******/
        stringstream ss;
        ss << "globaldata" << "_" << c <<  ".txt";
        ofstream wrout (ss.str().c_str(), std::ofstream::app);
        wrout << t << '\t' << knotcurves[c].writhe << '\t' << knotcurves[c].twist << '\t' << knotcurves[c].length << '\n';
        wrout.close();

        ss.str("");
        ss.clear();       

        ss << seriesname << c << "_" << t <<  ".vtk";
        ofstream knotout (ss.str().c_str());

        int i;
        int n = knotcurves[c].knotcurve.size();

        knotout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        knotout << "POINTS " << n << " float\n";

        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].xcoord << ' ' << knotcurves[c].knotcurve[i].ycoord << ' ' << knotcurves[c].knotcurve[i].zcoord << '\n';
        }

        knotout << "\n\nCELLS " << n << ' ' << 3*n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << 2 << ' ' << i << ' ' << incp(i,1,n) << '\n';
        }

        knotout << "\n\nCELL_TYPES " << n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << "3\n";
        }

        knotout << "\n\nPOINT_DATA " << n << "\n\n";

        //   knotout << "\nSCALARS Curvature float\nLOOKUP_TABLE default\n";
        //   for(i=0; i<n; i++)
        //   {
        //       knotout << knotcurves[c].knotcurve[i].curvature << '\n'; }

        knotout << "\nSCALARS Torsion float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].torsion << '\n';
        }

        //        knotout << "\nVECTORS A float\n";
        //        for(i=0; i<n; i++)
        //     {
        //          knotout << knotcurves[c].knotcurve[i].ax << ' ' << knotcurves[c].knotcurve[i].ay << ' ' << knotcurves[c].knotcurve[i].az << '\n';
        //      }

        //   knotout << "\nVECTORS V float\n";
        //   for(i=0; i<n; i++)
        //   {
        //       knotout << knotcurves[c].knotcurve[i].vx << ' ' << knotcurves[c].knotcurve[i].vy << ' ' << knotcurves[c].knotcurve[i].vz << '\n';
        //   }
        //            knotout << "\nVECTORS t float\n";
        //        for(i=0; i<n; i++)
        //        {
        //            knotout << knotcurves[c].knotcurve[i].tx << ' ' << knotcurves[c].knotcurve[i].ty << ' ' << knotcurves[c].knotcurve[i].tz << '\n';
        //        }
        //    knotout << "\nVECTORS n float\n";
        //    for(i=0; i<n; i++)
        //    {
        //        knotout << knotcurves[c].knotcurve[i].nx << ' ' << knotcurves[c].knotcurve[i].ny << ' ' << knotcurves[c].knotcurve[i].nz << '\n';
        //    }
        //        knotout << "\nVECTORS b float\n";
        //        for(i=0; i<n; i++)
        //        {
        //       knotout << knotcurves[c].knotcurve[i].bx << ' ' << knotcurves[c].knotcurve[i].by << ' ' << knotcurves[c].knotcurve[i].bz << '\n';
        //        }
        //    knotout << "\nVECTORS vdotn float\n";
        //    for(i=0; i<n; i++)
        //    {
        //        knotout << knotcurves[c].knotcurve[i].vdotnx << ' ' << knotcurves[c].knotcurve[i].vdotny << ' ' << knotcurves[c].knotcurve[i].vdotnz << '\n';
        //    }
        //    knotout << "\nVECTORS vdotb float\n";
        //    for(i=0; i<n; i++)
        //    {
        //        knotout << knotcurves[c].knotcurve[i].vdotbx << ' ' << knotcurves[c].knotcurve[i].vdotby << ' ' << knotcurves[c].knotcurve[i].vdotbz << '\n';
        //    }
        //    knotout << "\n\nCELL_DATA " << n << "\n\n";
        //    knotout << "\nSCALARS Writhe float\nLOOKUP_TABLE default\n";
        //    for(i=0; i<n; i++)
        //    {
        //        knotout << knotcurves[c].knotcurve[i].writhe << '\n';
        //    }

        //    knotout << "\nSCALARS Twist float\nLOOKUP_TABLE default\n";
        //    for(i=0; i<n; i++)
        //    {
        //        knotout << knotcurves[c].knotcurve[i].twist << '\n';
        //    }

        knotout << "\nSCALARS Length float\nLOOKUP_TABLE default\n";
        for(i=0; i<n; i++)
        {
            knotout << knotcurves[c].knotcurve[i].length << '\n';
        }
        knotout.close();
    }
}

int phi_file_read(vector<double>&phi,const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    string temp,buff;
    stringstream ss;
    ifstream fin (B_filename.c_str());
    int i,j,k,n;

    for(i=0;i<10;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                ss.clear();
                ss.str("");
                if(fin.good())
                {
                    if(getline(fin,buff))
                    {
                        ss << buff;
                        ss >> phi[n];
                    }
                }
                else
                {
                    cout << "Something went wrong!\n";
                    return 1;
                }
            }
        }
    }

    fin.close();

    return 0;
}

int uvfile_read(vector<double>&u, vector<double>&v, vector<double>& ku, vector<double>& kv, vector<double>& ucvx, vector<double>& ucvy,vector<double>& ucvz,griddata& griddata)
{
    string buff,datatype,dimensions,xdim,ydim,zdim;
    ifstream fin (B_filename.c_str());
    for(int i=0;i<4;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff) &&(i==2)) datatype = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }
    if(fin.good())
    {
        getline(fin,buff,' ' ); 
        if(getline(fin,buff,' ')) xdim = buff;
        if(getline(fin,buff,' ')) ydim = buff;
        if(getline(fin,buff,'\n')) zdim = buff;
    } 
    fin.close();
    int x = atoi(xdim.c_str()); 
    int y= atoi(ydim.c_str()); 
    int z = atoi(zdim.c_str()); 

    if(x!=griddata.Nx || y!=griddata.Ny ||z!=griddata.Nz)
    {
        cout << "CAREFUL! the gridsize you read in from the uv file isnt equal to the one you set! resizing to the uv file values read in \n";
        ucvx.resize(x*y*z);
        ucvy.resize(x*y*z);
        ucvz.resize(x*y*z);
        // better resize our scratchpad too
        ku.resize(4*x*y*z);
        kv.resize(4*x*y*z);
        u.resize(x*y*z);
        v.resize(x*y*z);

        griddata.Nx = x;
        griddata.Ny = y;
        griddata.Nz = z;
    }

    // grab the dimensions read in, resize u and v, and warn the user if they dont match!

    if(datatype.compare("ASCII")==0)
    {
        uvfile_read_ASCII(u,v,griddata);
    }
    else if(datatype.compare("BINARY")==0)
    {
        uvfile_read_BINARY(u,v,griddata);
    }
    return 0;
}

int uvfile_read_ASCII(vector<double>&u, vector<double>&v,const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    string temp,buff;
    stringstream ss;
    ifstream fin (B_filename.c_str());
    int i,j,k,n;

    for(i=0;i<10;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                ss.clear();
                ss.str("");
                if(fin.good())
                {
                    if(getline(fin,buff))
                    {
                        ss << buff;
                        ss >> u[n];
                    }
                }
                else
                {
                    cout << "Something went wrong!\n";
                    return 1;
                }
            }
        }
    }

    for(i=0;i<2;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                ss.clear();
                ss.str("");
                if(fin.good())
                {
                    if(getline(fin,buff)) ss << buff;
                    ss >> v[n];
                }
                else
                {
                    cout << "Something went wrong!\n";
                    return 1;
                }
            }
        }
    }

    fin.close();

    return 0;
}

int uvfile_read_BINARY(vector<double>&u, vector<double>&v,const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    string temp,buff;
    stringstream ss;
    ifstream fin (B_filename.c_str(), std::ios::in  | std::ios::binary);
    int i,j,k,n;

    for(i=0;i<10;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                char* memblock;
                char* swapped;
                memblock = new char [sizeof(float)];
                swapped = new char [sizeof(float)];
                fin.read(memblock,sizeof(float));
                ByteSwap(memblock, swapped);
                float value = 12;
                memcpy(&value, swapped, 4);
                u[n] = value;
                delete[] memblock;
                delete[] swapped;
            }
        }
    }

    for(i=0;i<3;i++)
    {
        if(fin.good())
        {
            if(getline(fin,buff)) temp = buff;
        }
        else
        {
            cout << "Something went wrong!\n";
            return 1;
        }
    }

    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n=pt(i,j,k,griddata);
                char* memblock;
                char* swapped;
                memblock = new char [sizeof(float)];
                swapped = new char [sizeof(float)];
                fin.read(memblock,sizeof(float));
                ByteSwap(memblock, swapped);
                float value = 12;
                memcpy(&value, swapped, 4);
                v[n] = value;
                delete[] memblock;
                delete[] swapped;
            }
        }
    }

    fin.close();

    return 0;
}
void resizebox(vector<double>&u,vector<double>&v,vector<double>&ucvx,vector<double>&ucvy,vector<double>&ucvz,vector<knotcurve>&knotcurves,vector<double>&ku,vector<double>&kv,griddata& oldgriddata)
{
    cout << "resizing box \n";
    int Nx = oldgriddata.Nx;
    int Ny = oldgriddata.Ny;
    int Nz = oldgriddata.Nz;
    double ucrit = -1.2;
    // first of all, take off the boundary; we set up the marked array to have 1's on the boudary of the box 
    std::vector<int>marked(u.size(),0);
    int shelllabel=1;
    for(int i=0;i<Nx;i++)
    {
        for(int j=0; j<Ny; j++)
        {
            for(int k=0; k<Nz; k++)   //Central difference
            {
                int n = pt(i,j,k,oldgriddata);
                if(i==0||i==Nx-1||j==0||j==Ny-1||k==0||k==Nz-1 && u[n]>ucrit) marked[n] =-1;
            }
        }
    }
    // okay , grow the shell
    growshell(u,marked,ucrit, oldgriddata);
    bool dontresize = false;
    for(int n = 0; n<u.size();n++)
    {
        if(marked[n]==-2)
        {
            marked[n]=shelllabel; 
            if(ucvx[n]*ucvx[n]+ucvy[n]*ucvy[n]+ucvz[n]*ucvz[n]>0.1) dontresize = true;
        }
    }
    shelllabel++;

    if (!dontresize)
    {
        bool hitinnershell = false;
        while(!hitinnershell)
        {
            int imax,jmax,kmax;
            imax = -1;
            jmax = -1;
            kmax = -1;
            // now we have no shells intersecting the boundary, but there may still be multiple shells before the knot; lets remove them one by one
            // to begin with , just grab some point on the outer shell
            for(int i=0;i<Nx;i++)
            {
                for(int j=0; j<Ny; j++)
                {
                    for(int k=0; k<Nz; k++)   //Central difference
                    {
                        int n = pt(i,j,k,oldgriddata);
                        if(u[n]>ucrit &&marked[n]==0 && i>imax && j>jmax && k> kmax) {imax = i ; jmax = j; kmax = k;}
                    }
                }
            }
            marked[pt(imax,jmax,kmax,oldgriddata)] = -1;
            // now grow the shell from here
            growshell(u,marked,ucrit, oldgriddata);
            for(int n = 0; n<u.size();n++)
            {
                if(marked[n]==-2)
                {
                    marked[n]=shelllabel;
                    if(ucvx[n]*ucvx[n]+ucvy[n]*ucvy[n]+ucvz[n]*ucvz[n]>0.1) hitinnershell = true;
                }
            }
            if(!hitinnershell)shelllabel++;
        }
        // at this point we have an array, marked, marked with integers increasing from the boudary, denoting shell numbers
        // we want to stip off all but the innermost shell
        // set everything outside this inner shell to the fixed point values 
        for(int n = 0; n<u.size();n++)
        {
            if(marked[n]>0 &&marked[n]<shelllabel){ u[n] = -1.03; v[n] = -0.66;}
        }
        // find the hull of the inner shell
        int imax = 0;
        int jmax = 0; 
        int kmax =0;
        int imin = Nx;
        int jmin =Ny; 
        int kmin =Nz; 
        for(int i = 0; i<Nx;i++)
        {
            for(int j = 0; j<Ny;j++)
            {
                for(int k = 0; k<Nz;k++)
                {
                    if(marked[pt(i,j,k,oldgriddata)] == shelllabel)
                    {
                        if(i>imax) imax = i;
                        if(j>jmax) jmax = j;
                        if(k>kmax) kmax = k;
                        if(i<imin) imin = i;
                        if(j<jmin) jmin = j;
                        if(k<kmin) kmin = k;
                    }
                }
            }
        }
        // we have our box dimensions in the ijk max min values already
        int deltai = imax - imin;
        int deltaj = jmax - jmin;
        int deltak = kmax - kmin;
        int N = (deltai<deltaj) ? deltaj:deltai;
        N = (N < deltak) ? deltak:N;

        griddata newgriddata;
        newgriddata.Nx = newgriddata.Ny = newgriddata.Nz = N;
        vector<double>utemp(N*N*N);
        vector<double>vtemp(N*N*N);
        for(int i = 0; i<N;i++)
        {
            for(int j = 0; j<N;j++)
            {
                for(int k = 0; k<N;k++)
                {
                    utemp[pt(i,j,k,newgriddata)] = u[pt(imin+i,jmin+j,kmin+k,oldgriddata)] ;
                    vtemp[pt(i,j,k,newgriddata)] = v[pt(imin+i,jmin+j,kmin+k,oldgriddata)] ;
                }
            }
        }
        // first of all, we can simply resize the ucvx data, since it gets recalculated anyhow
        ucvx.resize(N*N*N);
        ucvy.resize(N*N*N);
        ucvz.resize(N*N*N);
        // better resize our scratchpad too
        ku.resize(4*N*N*N);
        kv.resize(4*N*N*N);
        // the data is safely stored in the temp arrays, lets trash u and v
        u.resize(N*N*N);
        v.resize(N*N*N);
        u = utemp;
        v = vtemp;
        // finally, reset the grid data to the new griddata
        oldgriddata = newgriddata;
    }
    if(dontresize)
    {
        cout << "the inner shell is touching the boundary. Either the knot is spanning the whole box, or its across/very close to the box  boundary. For now, just aborting the resize \n" ;
    }
}
void growshell(vector<double>&u,vector<int>& marked,double ucrit, const griddata& griddata)
{
    bool stillboundaryleft = true;
    while(stillboundaryleft)
    {
        grow(u,marked,ucrit,griddata);
        stillboundaryleft = false;
        for(int n = 0; n<u.size();n++)
        {
            if(marked[n]==-1) stillboundaryleft =true;
        }

    }
    // okay we have our marked points - they are marked with a 2 in the marked array. lets set all the uv values we find their to the resting state values
}
void grow(const vector<double>&u,vector<int>&marked,double ucrit,const griddata& griddata)
{
    // the marked array has the following values
    // 0 - not evaluated 
    // -1 - a boundary, to be grown 
    // -2 - the interrior, already grown
    // -3 - a temporary state, marked as a boundary during the update
    // positive numbers - layers of shells already marked
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    for(int i=0;i<Nx;i++)
    {
        for(int j=0; j<Ny; j++)
        {
            for(int k=0; k<Nz; k++)   //Central difference
            {
                int n = pt(i,j,k,griddata);
                if(marked[n] ==-1)
                {

                    for(int iinc=-1;iinc<=1;iinc++)
                    {
                        for(int jinc=-1; jinc<=1; jinc++)
                        {
                            for(int kinc=-1; kinc<=1; kinc++)   //Central difference
                            {
                                int neighboringn = pt(incabsorb(i,iinc,Nx),incabsorb(j,jinc,Ny),incabsorb(k,kinc,Nz),griddata);
                                if(marked[neighboringn] == 0 && u[neighboringn] > ucrit) marked[neighboringn] = -3;
                            }
                        }
                    }
                }
            }
        }
    }
    for(int n = 0; n<u.size();n++)
    {
        if(marked[n]==-1){marked[n] =-2;}
        if(marked[n]==-3){marked[n] =-1;}
    }
}
int intersect3D_SegmentPlane( knotpoint SegmentStart, knotpoint SegmentEnd, knotpoint PlaneSegmentStart, knotpoint PlaneSegmentEnd, double& IntersectionFraction, std::vector<double>& IntersectionPoint )
{
    double ux = SegmentEnd.xcoord - SegmentStart.xcoord ;
    double uy = SegmentEnd.ycoord - SegmentStart.ycoord ;
    double uz = SegmentEnd.zcoord - SegmentStart.zcoord ;

    double wx= SegmentStart.xcoord - PlaneSegmentStart.xcoord ;
    double wy = SegmentStart.ycoord - PlaneSegmentStart.ycoord ;
    double wz = SegmentStart.zcoord - PlaneSegmentStart.zcoord ;

    double nx= PlaneSegmentEnd.xcoord  - PlaneSegmentStart.xcoord ;
    double ny = PlaneSegmentEnd.ycoord  - PlaneSegmentStart.ycoord ;
    double nz = PlaneSegmentEnd.zcoord  - PlaneSegmentStart.zcoord ;

    double D = nx*ux+ ny*uy + nz*uz;
    double N = - (nx*wx+ ny*wy + nz*wz);

    if (fabs(D) < 0.01)
    {           // segment is parallel to plane
        if (N == 0)                      // segment lies in plane
            return 2;
        else
            return 0;                    // no intersection
    }

    double sI = N / D;
    if (sI < 0 || sI > 1)
        return 0;                        // no intersection


    IntersectionFraction = sI;
    IntersectionPoint[0] = SegmentStart.xcoord + sI * ux;
    IntersectionPoint[1] = SegmentStart.ycoord + sI * uy;
    IntersectionPoint[2] = SegmentStart.zcoord + sI * uz;
    return 1;
}

double my_f(const gsl_vector* minimum, void* params)
{

    int i,j,k,idwn,jdwn,kdwn,modidwn,modjdwn,modkdwn,m,pts,iinc,jinc,kinc;
    double ucvxs, ucvys, ucvzs,  xd, yd ,zd, xdiff, ydiff, zdiff, prefactor;
    struct parameters* myparameters = (struct parameters *) params;
    griddata griddata = myparameters->mygriddata;
    likely::TriCubicInterpolator* interpolateducvmag = myparameters->ucvmag;
    double Nx = myparameters->mygriddata.Nx;
    double Ny = myparameters->mygriddata.Ny;
    double Nz = myparameters->mygriddata.Nz;
    gsl_vector* tempf = gsl_vector_alloc (3);
    gsl_vector* tempv = gsl_vector_alloc (3);
    gsl_vector* tempb = gsl_vector_alloc (3);
    gsl_vector_memcpy (tempf,myparameters->f);
    gsl_vector_memcpy (tempv,myparameters->v);
    gsl_vector_memcpy (tempb,myparameters->b);

    // s gives us how much of f to add to p

    gsl_vector_scale(tempf,gsl_vector_get (minimum, 0));
    gsl_vector_scale(tempb,gsl_vector_get (minimum, 1));
    gsl_vector_add(tempf,tempb);
    gsl_vector_add(tempv,tempf);
    double px = gsl_vector_get(tempv, 0);
    double py = gsl_vector_get(tempv, 1);
    double pz = gsl_vector_get(tempv, 2);
    gsl_vector_free(tempf);
    gsl_vector_free(tempv);
    gsl_vector_free(tempb);

    double value = -1*((*interpolateducvmag)(px,py,pz));
    return value;
}
void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
    double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
        - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

    double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
        - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

    double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
        - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

    gsl_vector_set(product, 0, p1);
    gsl_vector_set(product, 1, p2);
    gsl_vector_set(product, 2, p3);
}
void rotatedisplace(double& xcoord, double& ycoord, double& zcoord, const double theta, const double ux,const double uy,const double uz)
{

    double xprime = ( cos(theta) + (ux*ux)*(1-cos(theta)) )*xcoord + ( ux*uy*(1-cos(theta)) - uz*sin(theta) )*ycoord + ( ux*uz*(1-cos(theta)) + uy*sin(theta) )*zcoord;
    double yprime =(uy*ux*(1-cos(theta)) + uz*sin(theta) )*xcoord + ( cos(theta) + (uy*uy)*(1-cos(theta)) )*ycoord + ( uy*uz*(1-cos(theta)) - ux*sin(theta)  )*zcoord;
    double zprime = (uz*ux*(1-cos(theta)) - uy*sin(theta) )*xcoord + ( uz*uy*(1-cos(theta)) + ux*sin(theta)  )*ycoord + ( cos(theta) + (uz*uz)*(1-cos(theta)) )*zcoord;

    xcoord = xprime;
    ycoord = yprime;
    zcoord = zprime; 

}
inline int circularmod(int i, int N)    // mod i by N in a cirucler fashion, ie wrapping around both in the +ve and -ve directions
{
    if(i<0) return N - ((-i)%N);
    else return i%N;
}
// inlined functions for incrementing things respecting boundaries
inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}

inline int incw(int i, int p, int N)    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    if(i+p<0) return (-(i+p+1));
    if(i+p>N-1) return (2*N-(i+p+1));
    return (i+p);
}

inline int incabsorb(int i, int p, int N)    //increment with reflecting boundary between -1 and 0 and N-1 and N
{
    if(i+p<0) return (0);
    if(i+p>N-1) return (N-1);
    return (i+p);
}
// this function is specifically designed to incremenet, in the direction specified, respecting the boundary conditions, which are global enums
inline int gridinc(int i, int p, int N, int direction )    //increment with reflecting boundary between -1 and 0 and N-1 and N
{

    if(BoundaryType == ALLREFLECTING)
    {
        return incw(i,p,N);
    }

    if(BoundaryType == ALLPERIODIC)
    {
        return incp(i,p,N);
    }

    if(BoundaryType == ZPERIODIC)
    {
        if(direction ==2) return incp(i,p,N);
        else return incw(i,p,N);
    }
    return 0;
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
float FloatSwap( float f )
{    
    union
    {
        float f;
        char b[4];
    } dat1, dat2;

    dat1.f = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;
}

void ByteSwap(const char* TobeSwapped, char* swapped )
{    
    swapped[0] = TobeSwapped[3];
    swapped[1] = TobeSwapped[2];
    swapped[2] = TobeSwapped[1];
    swapped[3] = TobeSwapped[0];
    return; 
}
// this function takes two knots, and shifts the first one by grid spacing multiples until it literally lies over the second
void overlayknots(vector<knotcurve>& knotcurves,const vector<knotcurve>& knotcurvesold,const griddata& griddata)
{
    for(int c = 0; c <knotcurves.size();c++)
    {
        // how the average position is displaced between the two
        double deltax = knotcurves[c].xavgpos - knotcurvesold[c].xavgpos;
        double deltay = knotcurves[c].yavgpos - knotcurvesold[c].yavgpos;
        double deltaz = knotcurves[c].zavgpos - knotcurvesold[c].zavgpos;

        // we mod out by the lattice spacing in all of these
        // how many lattice spacings go into these deltas?
        int xlatticeshift = (int) (round(deltax/(griddata.Nx *h)));
        int ylatticeshift = (int) (round(deltay/(griddata.Ny *h)));
        int zlatticeshift = (int) (round(deltaz/(griddata.Nz *h)));

        // we should shift the knotcurve points by this much , to properly be able to compare the curves

        for(int s=0; s<knotcurves[c].knotcurve.size(); s++)
        {
            knotcurves[c].knotcurve[s].xcoord -= (double)(xlatticeshift) * (griddata.Nx *h);
            knotcurves[c].knotcurve[s].ycoord -= (double)(ylatticeshift) * (griddata.Ny *h);
            knotcurves[c].knotcurve[s].zcoord -= (double)(zlatticeshift) * (griddata.Nz *h);
        }
        // now we've done these shifts, we'd better move the knotcurve average position too.
        knotcurves[c].xavgpos = 0;
        knotcurves[c].yavgpos = 0;
        knotcurves[c].zavgpos = 0;
        double NP = knotcurves[c].knotcurve.size();
        for(int s=0; s<knotcurves[c].knotcurve.size(); s++)
        {
            knotcurves[c].xavgpos += knotcurves[c].knotcurve[s].xcoord/NP;
            knotcurves[c].yavgpos += knotcurves[c].knotcurve[s].ycoord/NP;
            knotcurves[c].zavgpos += knotcurves[c].knotcurve[s].zcoord/NP;
        }
    }
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
