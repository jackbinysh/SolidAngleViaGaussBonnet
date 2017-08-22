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

    cout << "Outputting Tangent Developable" << endl;
    OutputTangentDevelopable(Curve);
    // output the solid angle on a grid
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

    // refine the curve to increase number of points
    // while(Curve.NumPoints < 5000*Curve.NumComponents)
    {
        cout << "curve has size " << Curve.NumPoints << endl;
        //     RefineCurve(Curve);
        //     ComputeLengths(Curve);
        //     ComputeTangent(Curve);
        //     ComputeKappaB(Curve);
    }
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
    NewCurve.NumPoints=0;
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
            Curve.Components[i].knotcurve[s].curvature = sqrt(kappaNx*kappaNx + kappaNy*kappaNy + kappaNz*kappaNz);
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
            // Simpson's rule quadrature
            double ds0 = Curve.Components[i].knotcurve[s].length;
            double dsp1 = Curve.Components[i].knotcurve[incp(s,1,NP)].length;
            double dsm1 = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            double dsm2 = Curve.Components[i].knotcurve[incp(s,-2,NP)].length;
            double ds;
            if (s%2==0) ds = 0.5*(ds0+dsm1) + dsm1*dsm1/(6.0*ds0) + ds0*ds0/(6.0*dsm1);
            else ds = 0.5*(ds0+dsm1) + (dsm2-dsm1-ds0+dsp1)/6.0 - dsm2*dsm2/(6.0*dsm1) - dsp1*dsp1/(6.0*ds0); 
            // and here's the integrand
            Angle += ds*ndotT*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist*dist - ndotT*ndotT);
        }

        Angle /= 2.0*M_PI;  
        Curve.Components[i].twist=Angle;
    }
}

// a function which will output the entire integrand at a given point
void OutputIntegrands(const Link& Curve, const viewpoint& View, int Pointi, int  Pointj, int  Pointk)
{

    for(int i=0; i<Curve.NumComponents; i++)
    {
        std::ostringstream oss;
        oss << "PlusIntegrand_Component" << i << "_"  << Pointk << ".txt";
        string name = oss.str();
        ofstream Aout (name.c_str());
        std::ostringstream oss2;
        oss2 << "MinusIntegrand_Component" << i << "_"  << Pointk << ".txt";
        name = oss2.str();
        ofstream A2out (name.c_str());
        // run over the knot 
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
            // and here's the integrand
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            double PlusIntegrand = (viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
            double MinusIntegrand = (viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist - ndotT);
            Aout << PlusIntegrand*ds << '\n';
            A2out << MinusIntegrand*ds << '\n';
        }
        Aout.close();
        A2out.close();
    }

    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();

        // check proximity to tangent developable surface
        int M = 16; // even integer for defining a window around the cusp point
        int window[NP]; // window of integration
        double Integrand[NP]; // window of integration
        for (int s=0; s<NP; s++) {window[s]=1;} // there must be something smarter than this!!
        int count = 0; // keep track of the number of cusp points we encounter
        for (int s=0; s<NP; s++)
        {

            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
            double ndotT = viewx*Curve.Components[i].knotcurve[s].tx + viewy*Curve.Components[i].knotcurve[s].ty + viewz*Curve.Components[i].knotcurve[s].tz;
            // calculate dot product with the Frenet normal
            double ndotN = viewx*Curve.Components[i].knotcurve[s].kappaNx + viewy*Curve.Components[i].knotcurve[s].kappaNy + viewz*Curve.Components[i].knotcurve[s].kappaNz;
            ndotN /= (dist*Curve.Components[i].knotcurve[s].curvature); // shouldn't really be necessary
            // repeat for the next point along
            double viewx2 = Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - View.xcoord;
            double viewy2 = Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - View.ycoord;
            double viewz2 = Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - View.zcoord;
            double dist2 = sqrt(viewx2*viewx2 + viewy2*viewy2 + viewz2*viewz2);      
            double ndotN2 = viewx*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNx + viewy*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNy + viewz*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNz;
            ndotN2 /= (dist2*Curve.Components[i].knotcurve[incp(s,1,NP)].curvature); // shouldn't really be necessary
            // only continue if ndotN changes sign and ndotT is negative
            if ((ndotN*ndotN2<0) && (ndotT<0))
            {
                // tag the point of the curve
                int sp = s;
                if ((-ndotN2/ndotN)<1) sp = incp(s,1,NP);
                // compute ndotT and ndotB
                viewx = Curve.Components[i].knotcurve[sp].xcoord - View.xcoord;
                viewy = Curve.Components[i].knotcurve[sp].ycoord - View.ycoord;
                viewz = Curve.Components[i].knotcurve[sp].zcoord - View.zcoord;
                dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
                // note the minus sign introduced here -- makes the quantity positive
                ndotT = -(viewx*Curve.Components[i].knotcurve[sp].tx + viewy*Curve.Components[i].knotcurve[sp].ty + viewz*Curve.Components[i].knotcurve[sp].tz)/dist;
                // check if a threshold is exceeded
                double threshold = Curve.Components[i].knotcurve[sp].curvature*M*(Curve.Components[i].length/NP)/(2.0*sqrt(2.0)*sqrt(1.0/ndotT - 1.0));
                cout << threshold << '\n';
                if (threshold > 1)
                {
                    // note the minus sign introduced here -- makes the quantity positive
                    double ndotB = -(viewx*Curve.Components[i].knotcurve[sp].kappaBx + viewy*Curve.Components[i].knotcurve[sp].kappaBy + viewz*Curve.Components[i].knotcurve[sp].kappaBz)/(dist*Curve.Components[i].knotcurve[sp].curvature);
                    // value of the integral over the narrow Lorentzian peak
                    //Integrand[sp] = (-2.0*sqrt(2.0)*(ndotB/ndotT)/sqrt(1.0/ndotT - 1.0))*atan(threshold);
                    double theta = ( ndotB>0? 1 :-1 )*acos(ndotT);
                    double k =Curve.Components[i].knotcurve[sp].curvature ;
                    double avgds = Curve.Components[i].length/(double)NP;
                    double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
                    Integrand[sp] = -2 * ( ( k*theta )/(theta*theta)  )*ds;
                    double Iplus = (-2.0*sqrt(2.0)*(ndotB/ndotT)/sqrt(1.0/ndotT - 1.0))*atan(threshold);
                    cout << "value" << Iplus << '\n';
                    // set an excluded window about sp
                    window[sp] = 0;
                    for (int m=1; m<=M/2; m++) 
                    {
                        window[incp(sp,m,NP)] = 0;
                        window[incp(sp,-m,NP)] = 0;
                        Integrand[incp(sp,m,NP)] = -2 * ( ( k*theta )/(theta*theta +((double)(m)*avgds*k)*((double)(m)*avgds*k) )  )*ds;
                        Integrand[incp(sp,-m,NP)] = -2 * ( ( k*theta )/(theta*theta +( (double)(m)*avgds*k )*( (double)(m)*avgds*k) )*ds  );
                    }
                    count++;
                    cout << "cusp found at s = " << s << endl;
                }
            }
        }

        // now do the remaining integral as usual
        for (int s=0; s<NP; s++) 
        {
            if (window[s]==1)
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
                // Simpson's rule quadrature
               // double ds0 = Curve.Components[i].knotcurve[s].length;
               // double dsp1 = Curve.Components[i].knotcurve[incp(s,1,NP)].length;
               // double dsm1 = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
               // double dsm2 = Curve.Components[i].knotcurve[incp(s,-2,NP)].length;
               // double ds;
               // if (s%2==0) ds = 0.5*(ds0+dsm1) + dsm1*dsm1/(6.0*ds0) + ds0*ds0/(6.0*dsm1);
               // else ds = 0.5*(ds0+dsm1) + (dsm2-dsm1-ds0+dsp1)/6.0 - dsm2*dsm2/(6.0*dsm1) - dsp1*dsp1/(6.0*ds0); 
                double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
                // and here's the integrand
                double tempIntegrand = (viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
                // Iplus += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
                // cheap trick
                if ((window[incp(s,1,NP)]==0)||(window[incp(s,-1,NP)]==0))
                {
                    //  tempIntegrand += 0.5*(M-1)*ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
                    // Iplus += 0.5*(M-1)*ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
                }
                Integrand[s]=tempIntegrand*ds;
            }
        }

        std::ostringstream oss;
        oss << "RegularisedPlusIntegrand_Component" << i << "_"  << Pointk << ".txt";
        string name = oss.str();
        ofstream A3out (name.c_str());

        for (int s=0; s<NP; s++){A3out << Integrand[s] << '\n';};
        A3out.close();
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
            // Simpson's rule quadrature
          //  double ds0 = Curve.Components[i].knotcurve[s].length;
          //  double dsp1 = Curve.Components[i].knotcurve[incp(s,1,NP)].length;
          //  double dsm1 = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
          //  double dsm2 = Curve.Components[i].knotcurve[incp(s,-2,NP)].length;
          //  double ds;
          //  if (s%2==0) ds = 0.5*(ds0+dsm1) + dsm1*dsm1/(6.0*ds0) + ds0*ds0/(6.0*dsm1);
          //  else ds = 0.5*(ds0+dsm1) + (dsm2-dsm1-ds0+dsp1)/6.0 - dsm2*dsm2/(6.0*dsm1) - dsp1*dsp1/(6.0*ds0); 
            // and here's the integrand
                double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
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
            // Simpson's rule quadrature
            double ds0 = Curve.Components[i].knotcurve[s].length;
            double dsp1 = Curve.Components[i].knotcurve[incp(s,1,NP)].length;
            double dsm1 = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            double dsm2 = Curve.Components[i].knotcurve[incp(s,-2,NP)].length;
            double ds;
            if (s%2==0) ds = 0.5*(ds0+dsm1) + dsm1*dsm1/(6.0*ds0) + ds0*ds0/(6.0*dsm1);
            else ds = 0.5*(ds0+dsm1) + (dsm2-dsm1-ds0+dsp1)/6.0 - dsm2*dsm2/(6.0*dsm1) - dsp1*dsp1/(6.0*ds0); 
            // and here's the integrand
            Iminus += ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist - ndotT);
        }
        componentomega = 2.0*M_PI*(1.0 + Curve.Components[i].writhe)-Iminus;
        totalomega += componentomega;
    }
    return totalomega;
}


// the regularised one
double SolidAngleCalcR(const Link& Curve, const viewpoint& View)
{
    double totalomega =0;
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double Iplus = 0.0;
        int NP = Curve.Components[i].knotcurve.size();

        // check proximity to tangent developable surface
        int M = 4; // even integer for defining a window around the cusp point
        int window[NP]; // window of integration
        for (int s=0; s<NP; s++) {window[s]=1;} // there must be something smarter than this!!
        int count = 0; // keep track of the number of cusp points we encounter
        for (int s=0; s<NP; s++)
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);      
            double ndotT = viewx*Curve.Components[i].knotcurve[s].tx + viewy*Curve.Components[i].knotcurve[s].ty + viewz*Curve.Components[i].knotcurve[s].tz;
            // calculate dot product with the Frenet normal
            double ndotN = viewx*Curve.Components[i].knotcurve[s].kappaNx + viewy*Curve.Components[i].knotcurve[s].kappaNy + viewz*Curve.Components[i].knotcurve[s].kappaNz;
            ndotN /= (dist*Curve.Components[i].knotcurve[s].curvature); // shouldn't really be necessary
            // repeat for the next point along
            double viewx2 = Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - View.xcoord;
            double viewy2 = Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - View.ycoord;
            double viewz2 = Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - View.zcoord;
            double dist2 = sqrt(viewx2*viewx2 + viewy2*viewy2 + viewz2*viewz2);      
            double ndotN2 = viewx*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNx + viewy*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNy + viewz*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNz;
            ndotN2 /= (dist2*Curve.Components[i].knotcurve[incp(s,1,NP)].curvature); // shouldn't really be necessary
            // only continue if ndotN changes sign and ndotT is negative
            if ((ndotN*ndotN2<0) && (ndotT<0))
            {
                // tag the point of the curve
                int sp = s;
                if ((-ndotN2/ndotN)<1) sp = incp(s,1,NP);
                // compute ndotT and ndotB
                viewx = Curve.Components[i].knotcurve[sp].xcoord - View.xcoord;
                viewy = Curve.Components[i].knotcurve[sp].ycoord - View.ycoord;
                viewz = Curve.Components[i].knotcurve[sp].zcoord - View.zcoord;
                dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
                // note the minus sign introduced here -- makes the quantity positive
                ndotT = -(viewx*Curve.Components[i].knotcurve[sp].tx + viewy*Curve.Components[i].knotcurve[sp].ty + viewz*Curve.Components[i].knotcurve[sp].tz)/dist;
                // check if a threshold is exceeded
                double threshold = Curve.Components[i].knotcurve[sp].curvature*M*(Curve.Components[i].length/NP)/(2.0*sqrt(2.0)*sqrt(1.0/ndotT - 1.0));
                if (threshold > 1)
                {
                    // note the minus sign introduced here -- makes the quantity positive
                    double ndotB = (-viewx*Curve.Components[i].knotcurve[sp].kappaBx + viewy*Curve.Components[i].knotcurve[sp].kappaBy + viewz*Curve.Components[i].knotcurve[sp].kappaBz)/(dist*Curve.Components[i].knotcurve[sp].curvature);
                    // value of the integral over the narrow Lorentzian peak
                    Iplus += (-2.0*sqrt(2.0)*(ndotB/ndotT)/sqrt(1.0/ndotT - 1.0))*atan(threshold);
                    // set an excluded window about sp
                    window[sp] = 0;
                    for (int m=1; m<=M/2; m++) 
                    {
                        window[incp(sp,m,NP)] = 0;
                        window[incp(sp,-m,NP)] = 0;
                    }
                    count++;
                   // cout << "cusp found at s = " << s << endl;
                }
            }
        }

        // now do the remaining integral as usual
        for (int s=0; s<NP; s++) 
        {
            if (window[s]==1)
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
                if(window[incp(s,1,NP)]==0 || window[incp(s,-1,NP)]==0){factor = 0.5;}
                // Simpson's rule quadrature
                //   double ds0 = Curve.Components[i].knotcurve[s].length;
                //   double dsp1 = Curve.Components[i].knotcurve[incp(s,1,NP)].length;
                //   double dsm1 = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
                //   double dsm2 = Curve.Components[i].knotcurve[incp(s,-2,NP)].length;
                //   double ds;
                //   if (s%2==0) ds = 0.5*(ds0+dsm1) + dsm1*dsm1/(6.0*ds0) + ds0*ds0/(6.0*dsm1);
                //   else ds = 0.5*(ds0+dsm1) + (dsm2-dsm1-ds0+dsp1)/6.0 - dsm2*dsm2/(6.0*dsm1) - dsp1*dsp1/(6.0*ds0); 
                // trapezium rule quadrature
                //      ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
                // and here's the integrand
                Iplus += ds*factor*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
                // cheap trick
                if ((window[incp(s,1,NP)]==0)||(window[incp(s,-1,NP)]==0))
                {
                //    Iplus += 0.5*(M-1)*ds*(viewx*kappaBx + viewy*kappaBy + viewz*kappaBz)/(dist + ndotT);
                }
            }
        }
        double componentomega = 2.0*M_PI*(1.0 + Curve.Components[i].writhe)-Iplus;
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

    // header stuff for the vtk file
    fn = "solid_angle_regular.vtk";
    ofstream Bout (fn.c_str());

    Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Bout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    Bout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    Bout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    Bout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    Bout << "SCALARS Phi float\nLOOKUP_TABLE default\n";
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

                if( i == Nx/2 && j == Ny-3)
                {
                    OutputIntegrands(Curve, Point,i,j,k);
                }

                // is ndotT very close to 1? if so, better us 1/(1+ndotT)
                // if(fabs(1 - ndotTmax) < fabs(-1-ndotTmin)) 
                {
                    SolidAngle = SolidAngleCalc(Curve,Point);
                }
                // else 
                {
                    //      SolidAngle = SolidAngleCalc2(Curve,Point);
                }

                // put in the interval [0,4pi]
                while(SolidAngle>4*M_PI) SolidAngle -= 4*M_PI;
                while(SolidAngle<0) SolidAngle += 4*M_PI;

                Aout << SolidAngle << '\n';	      

                SolidAngle = SolidAngleCalcR(Curve,Point);

                // put in the interval [0,4pi]
                while(SolidAngle>4*M_PI) SolidAngle -= 4*M_PI;
                while(SolidAngle<0) SolidAngle += 4*M_PI;

                Bout << SolidAngle << '\n';	      
            }
        }
    }
    Aout.close();
    Bout.close();
}

void OutputTangentDevelopable(const Link& Curve)
{
    // header stuff for the vtk file
    string fn = "Tangent_Developable.txt";
    ofstream Aout (fn.c_str());

    griddata griddata;
    griddata.Nx = Nx;
    griddata.Ny = Ny;
    griddata.Nz = Nz;

    double spacing = 5;

    // the forward branch
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++) 
        {

            double xcoord = 0;
            double ycoord = 0;
            double zcoord = 0;
            double lambda = 0;
            while( xcoord < x(Nx,griddata) && xcoord > x(0,griddata) && ycoord < y(Ny,griddata) && ycoord > y(0,griddata) && zcoord < z(Nz,griddata) && zcoord > z(0,griddata) )
            {
                xcoord = Curve.Components[i].knotcurve[s].xcoord + lambda*Curve.Components[i].knotcurve[s].tx;
                ycoord = Curve.Components[i].knotcurve[s].ycoord + lambda*Curve.Components[i].knotcurve[s].ty;
                zcoord = Curve.Components[i].knotcurve[s].zcoord + lambda*Curve.Components[i].knotcurve[s].tz;
                Aout << xcoord <<' '<<ycoord<<' '<<zcoord<< ' '<< 0 << '\n'; 	      
                lambda += spacing;
            }
        }
    }
    // the backward branch
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++) 
        {

            double xcoord = 0;
            double ycoord = 0;
            double zcoord = 0;
            double lambda = -spacing;
            while( xcoord < x(Nx,griddata) && xcoord > x(0,griddata) && ycoord < y(Ny,griddata) && ycoord > y(0,griddata) && zcoord < z(Nz,griddata) && zcoord > z(0,griddata) )
            {
                xcoord = Curve.Components[i].knotcurve[s].xcoord + lambda*Curve.Components[i].knotcurve[s].tx;
                ycoord = Curve.Components[i].knotcurve[s].ycoord + lambda*Curve.Components[i].knotcurve[s].ty;
                zcoord = Curve.Components[i].knotcurve[s].zcoord + lambda*Curve.Components[i].knotcurve[s].tz;
                Aout << xcoord <<' '<<ycoord<<' '<<zcoord<<' '<< 1 <<'\n'; 	      
                lambda -= spacing;
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
