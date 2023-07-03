#include <math.h>
// This is a collection of helpfunctions that were used to run convergence tests and benchmarks for the ellipse/ellipsoid
// in the particle closest point draft. Computing the theoretical closest points and distances from a given query point
// is done using the first four functions which were adopted from David Eberly "Distance from a Point to an Ellipse, an 
// Ellipsoid, or a Hyperellipsoid", 2013.
//
// Created by lschulze
double GetRoot (double r0, double z0, double z1, double g)
    {
	const int maxIter = 100;
        double n0 = r0*z0;
        double s0 = z1 - 1;
        double s1 = ( g < 0 ? 0 : sqrt(n0*n0+z1*z1) - 1 ) ;
        double s = 0;
        for ( int i = 0; i < maxIter; ++i ){
	    if (i == (maxIter - 1)) std::cout<<"distance point ellipse algorithm did not converge."<<std::endl;
            s = ( s0 + s1 ) / 2 ;
            if ( s == s0 || s == s1 ) {break; }
            double ratio0 = n0 /( s + r0 );
            double ratio1 = z1 /( s + 1 );
            g = ratio0*ratio0 + ratio1*ratio1 - 1 ;
            if (g > 0) {s0 = s;} else if (g < 0) {s1 = s ;} else {break ;}
        }
        return s;
    }

double GetRoot(double r0, double r1, double z0, double z1, double z2, double g)
{	
	const int maxIter = 100;
	double n0 = r0*z0;
	double n1 = r1*z1;
	double s0 = z2 - 1;
	double s1 = (g < 0 ? 0 : sqrt(n0*n0 + n1*n1 + z2*z2) - 1) ;
	double s = s;
	for(int i = 0 ; i < maxIter ; ++i )
	{
		if (i == (maxIter - 1)) std::cout<<"distance point ellipse algorithm did not converge."<<std::endl;
		s = ( s0 + s1 ) / 2 ;
		if ( s == s0 || s == s1 ) {break; }
		double ratio0 = n0 / ( s + r0 ); 
		double ratio1 = n1 / ( s + r1 );
		double ratio2 = z2 / ( s + 1 );
		g = ratio0*ratio0 + ratio1*ratio1 +ratio2*ratio2 - 1;
		if ( g > 0 ) { s0 = s ;} 
		else if ( g < 0 ) { s1 = s ; }
		else {break;}
	}
	return (s);
}

double DistancePointEllipse(double e0, double e1, double y0, double y1, double& x0, double& x1)
    {
        double distance;
        if ( y1 > 0){
            if ( y0 > 0){
                double z0 = y0 / e0;
                double z1 = y1 / e1;
                double g = z0*z0+z1*z1 - 1;
                if ( g != 0){
                    double r0 = (e0/e1)*(e0/e1);
                    double sbar = GetRoot(r0 , z0 , z1 , g);
                    x0 = r0 * y0 /( sbar + r0 );
                    x1 = y1 /( sbar + 1 );
                    distance = sqrt( (x0-y0)*(x0-y0) + (x1-y1)*(x1-y1) );
                    }else{
                        x0 = y0;
                        x1 = y1;
                        distance = 0;
                    }
                }
                else // y0 == 0
                    {x0 = 0 ; x1 = e1 ; distance = abs( y1 - e1 );}
        }else{ // y1 == 0
            double numer0 = e0*y0 , denom0 = e0*e0 - e1*e1;
            if ( numer0 < denom0 ){
                    double xde0 = numer0/denom0;
                    x0 = e0*xde0 ; x1 = e1*sqrt(1 - xde0*xde0 );
                    distance = sqrt( (x0-y0)*(x0-y0) + x1*x1 );
                }else{
                    x0 = e0;
                    x1 = 0;
                    distance = abs( y0 - e0 );
            }
        }
        return distance;
    }

double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2, double& x0, double& x1, double& x2)
{
	double distance;
	if( y2 > 0 )
	{
		if( y1 > 0 )
		{
			if( y0 > 0 )
			{
				double z0 = y0 / e0;
				double z1 = y1 / e1;
				double z2 = y2 / e2;
				double g = z0*z0 + z1*z1 + z2*z2 - 1 ;
				if( g != 0 )
				{
					double r0 = (e0/e2)*(e0/e2);
					double r1 = (e1/e2)*(e1/e2);
					double sbar = GetRoot ( r0 , r1 , z0 , z1 , z2 , g );
					x0 = r0 *y0 / ( sbar + r0 );
					x1 = r1 *y1 / ( sbar + r1 );
					x2 = y2 / ( sbar + 1 );
					distance = sqrt( (x0 - y0)*(x0 - y0) + (x1 - y1)*(x1 - y1) + (x2 - y2)*(x2 - y2));
				}
				else
				{
					x0 = y0;
					x1 = y1;
					x2 = y2;
					distance = 0;
				}
			}
			else // y0 == 0
			{
				x0 = 0;
				distance = DistancePointEllipse( e1 , e2 , y1 , y2, x1, x2);
			}
		}
		else // y1 == 0
		{
			if( y0 > 0 )
			{
				x1 = 0;
				distance = DistancePointEllipse( e0 , e2 , y0 , y2, x0, x2);
			}
			else // y0 == 0
			{
				x0 = 0;
				x1 = 0;
				x2 = e2;
				distance = abs(y2 - e2);
			}
		}
	}
	else // y2 == 0
	{
		double denom0 = e0*e0 - e2*e2;
		double denom1 = e1*e1 - e2*e2;
		double numer0 = e0*y0;
		double numer1 = e1*y1;
		bool computed = false;
		if((numer0 < denom0) && (numer1 < denom1))
		{
			double xde0 = numer0/denom0;
			double xde1 = numer1/denom1 ;
			double xde0sqr = xde0 *xde0;
			double xde1sqr = xde1 * xde1 ;
			double discr = 1 - xde0sqr - xde1sqr;
			if( discr > 0 )
			{
				x0 = e0*xde0;
				x1 = e1*xde1;
				x2 = e2*sqrt(discr);
				distance = sqrt((x0 - y0)*(x0 - y0) + (x1 - y1)*(x1 - y1) + x2*x2);
				computed = true;
			}
		}
		if( !computed )
		{
			x2 = 0;
			distance = DistancePointEllipse(e0 , e1 , y0 , y1, x0, x1);
		}
	}
	return distance;
}

constexpr unsigned int factorial(unsigned int x)
{
    unsigned int fact = 1;
    for(int i = 1; i < (x + 1); i++) fact = fact*i;
    return(fact);
}
constexpr unsigned int minter_lp_degree_one_num_coeffs(unsigned int dims, unsigned int poly_degree)
{
    return(factorial(dims + poly_degree)/(factorial(dims)*factorial(poly_degree)));
}

int return_sign(double phi)
{
	if (phi > 0) return 1;
	if (phi < 0) return -1;
	return 0;
}

double randMinusOneToOne()
{
    double temp = rand() / (RAND_MAX + 1.0);
    //std::cout<<(2.0*temp - 1.0)<<std::endl;
    return(2.0*temp - 1.0);
}
