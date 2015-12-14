/*
 * Copyrigtht 2014 Georgios Karagiannis
 *
 * This file is part of PISAA_ProtAB3D.
 *
 * PISAA_ProtAB3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 2 of the License.
 *
 * PISAA_ProtAB3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISAA_ProtAB3D.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Georgios Karagiannis 
 * Postdoctoral research associate
 * Department of Mathematics, Purdue University
 * 150 N. University Street
 * West Lafayette, IN 47907-2067, USA
 *
 * Telephone: +1 (765) 496-1007
 *
 * Email: gkaragia@purdue.edu
 *
 * Contact email: georgios.stats@gmail.com
*/


/*Spherical Coordinates
 *
 * Spherical coordinates, also called spherical polar coordinates
 * (Walton 1967, Arfken 1985), are a system of curvilinear coordinates that are
 * natural for describing positions on a sphere or spheroid. Define theta to be
 * the azimuthal angle in the xy-plane from the x-axis with 0<=theta<2pi
 * (denoted lambda when referred to as the longitude), phi to be the polar
 * angle (also known as the zenith angle and colatitude, with phi=90
 * degrees-delta where delta is the latitude) from the positive z-axis with
 * 0<=phi<=pi, and r to be distance (radius) from a point to the origin. This
 * is the convention commonly used in mathematics.
 * http://mathworld.wolfram.com/SphericalCoordinates.html
 * (radial, azimuthal, polar)	(r,theta,phi)	this work
 *
 *  The spherical coordinates (r,theta,phi) are related to the
 *  Cartesian coordinates (x,y,z) by
 *  r = sqrt(x^2+y^2+z^2)
 *  theta = tan^(-1)(y/x)
 *  phi	= cos^(-1)(z/r),
 *  where r in [0,infty), theta in [0,2pi), and phi in [0,pi],
 *  and the inverse tangent must be suitably defined to take
 *  the correct quadrant of (x,y) into account.
 *
 *  In terms of Cartesian coordinates,
 *  x = r*cos(theta)*sin(phi)
 *  y = r*sin(theta)*sin(phi)
 *  z = r*cos(phi).
*/

/* In this work, following the mathematics convention, the symbols for the
 * radial, azimuth, and zenith angle coordinates are taken as r, theta,
 * and phi, respectively. Note that this definition provides a logical
 * extension of the usual polar coordinates notation, with theta remaining the
 * angle in the xy-plane and phi becoming the angle out of that plane.
 *
 * */

#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include <stdio.h>

#ifndef INFINITY
	#include <float.h>
	#define INFINITY DBL_MAX
#endif

static int *seq = NULL ;

void get_data(char file_name[], int N_monomer, int *N_dimension){

	int i ;
	FILE *ins;

	/*N_monomer = (N_dimension+5)/2 ;*/

	*N_dimension = 2*N_monomer -5 ;

	seq = ivector(1,N_monomer) ;

	ins = fopen(file_name, "r") ;

	if (ins==NULL) {
		printf("No data set\n") ;
		abort() ;
	}

	for(i=1; i<=N_monomer; i++)
		fscanf(ins, " %d", &seq[i]);

	fclose(ins);
}

void cost_bounds(double *z_min, double *z_max, int i){
	if (i > 0){
		*z_min = (double) -INFINITY ;
		*z_max = (double) INFINITY ;
	}
}

double cost(double *angles, int N_dimension){

	int N_monomer ;
	int i ;
	int j ;

	double *theta ;
	double *phi ;

	double *x ;
	double *y ;
	double *z ;

	double V_theta ;
	double V_phi ;
	double V_LJ ;

	double r = 1.0 ; /* adjustent distance */

	double C_ij ;
	double d_ij ;

	double mypi = 3.14159265358979 ;

	N_monomer = (int) (N_dimension+5)/2 ;

	x = dvector(1, N_monomer) ;
	y = dvector(1, N_monomer) ;
	z = dvector(1, N_monomer) ;

	/* discriminate the azimuthial (theta) and polar (phi) coordinates*/

	/* theta[3:N_monomer] = angles[1:(N_monomer-2)] */
	theta = angles + 1 -3 ;
	/* wrap the azimuthial coordinates to [0,2*pi) */
	for (i=3; i<=N_monomer; i++)
		theta[i] -= ( floor( theta[i]/(2.0*mypi) ) * (2.0*mypi) ) ;

	/* phi[4:N_monomer] = angles[(N_monomer-1):(2*N_monomer-5)] */
	phi = angles + N_monomer-1 -4 ;
	/* wrap the polar coordinates to [0,pi] */
	for (i=4; i<=N_monomer; i++)
		phi[i] -= ( floor( phi[i]/(mypi) ) * (mypi) ) ;

	/* compute the bond vectors
	 * (... unit vectors represented by Cartecian coordinates)
	 * unit vector u(i) = (x[1],y[i],z[i]) connects monomer M(i)
	 * to monomer M(i+1) .
	 * */

	/* theta_1 = 0.0 ; phi_1 = pi/2  */
	x[2] = r * 1.0 ;
	y[2] = r * 0.0 ;
	z[2] = r * 0.0 ;

	/* phi_2 = pi/2 */
	x[3] = r * cos(theta[3]) ;
	y[3] = r * sin(theta[3]) ;
	z[3] = r * 0.0 ;

	for (i=4; i<=N_monomer; i++){
		x[i] = r * cos(theta[i]) * sin(phi[i]) ;
		y[i] = r * sin(theta[i]) * sin(phi[i]) ;
		z[i] = r * cos(phi[i]) ;
	}

	/* compute the contribution of the bond angles */

	V_theta = 0.0 ;
	for (i=1; i<=N_monomer-2; i++)
		V_theta += ( x[(i) +1]*x[(i+1) +1]
					+y[(i) +1]*y[(i+1) +1]
					+z[(i) +1]*z[(i+1) +1] ) ;

	/* compute the contribution of the torsional angles */

	V_phi = 0.0 ;
	for (i=1; i<=N_monomer-3; i++)
		V_phi -= 0.5*( x[(i) +1]*x[(i+2) +1]
			+y[(i) +1]*y[(i+2) +1]
			+z[(i) +1]*z[(i+2) +1] ) ;

	/* compute the Cartesian coordinates of each monomer
	 * (x[i],y[i],z[i]) stores the Cartecian coordinates of the M(i)-th monomer
	 * or equivalently, (x[i],y[i],z[i]) is the vector that connects (0,0,0)
	 * to the M(i)-th monomer .
	 * */

	x[1] = 0.0 ;
	y[1] = 0.0 ;
	z[1] = 0.0 ;
	for (i=2; i<=N_monomer; i++){
		x[i] = x[i-1] + x[i] ;
		y[i] = y[i-1] + y[i] ;
		z[i] = z[i-1] + z[i] ;
	}

	/* compute the contribution of the Lennard-Johnes */

	V_LJ = 0.0 ;
	for(i=1; i<=N_monomer-2; i++)
		for(j=i+2; j<=N_monomer; j++){
			C_ij = (seq[i]==1 && seq[j]==1) ? 1.0 : 0.5 ;
			d_ij = (x[i]-x[j])*(x[i]-x[j])
					+(y[i]-y[j])*(y[i]-y[j])
					+(z[i]-z[j])*(z[i]-z[j]) ;
			d_ij = sqrt( d_ij ) ;
			V_LJ += 4.0 * C_ij * ( pow(d_ij,-12.0) - pow(d_ij,-6.0) ) ;
		}

	free_dvector(x, 1, N_monomer) ;
	free_dvector(y, 1, N_monomer) ;
	free_dvector(z, 1, N_monomer) ;

	return V_theta + V_phi + V_LJ ;

}
