/*************************************************************************

	compile	: cc anb.c vec.c -g -lm -o anb
	input   : 'name'.in		for parameters check comments in example
	output	: 'name'.vec	readable by PAW (ascii file)
	  contents of each row:
	    photon energy, total coherent Intensity, diff coherent I = perp-para,
	    incoherent crystal I, incoherent amorphous I,
	    electron contrib crystal I, electron amorphous I
	  use PAW to calculate:
	    Irelexp = (Isum+Iinc)/Iamo  for experimental comparison
	    Irel    = (Isum+Iinc)/Iinc
	    Pol     = Idif/(Isum+Iinc)

*************************************************************************/


/*******************************************************
		INCLUDES
********************************************************/

#include <stdio.h>
#include <string.h>
#include<stdlib.h> // for exit function PenkovGS(c)
#include <math.h>
#include <time.h>
#include "vec.h"



/*******************************************************
		DEFINES
********************************************************/

#define TRUE 			1
#define FALSE			0

#define BEAMX_STEPS		20													/* no of steps for crystal angle variation - beam divergence */
#define BEAMY_STEPS		20
#define BEAMXY_REGION	2
#define MAX_LATTICEVEC  1000												/* maximum no of latticevectors */
#define MASS_E			0.511003373282										/* elektron mass in MeV as energy unit */
#define A				923.7												/* Diamond spacing, lattice constants */
#define M_2PI 			(2 * M_PI)
#define GUNIT			M_2PI/A
#define DEBYE_ON		TRUE 												/* for amorphous (amo) calc no debye - factor is needed */
#define DEBYE_OFF		FALSE												/* but for the incoherent crystal contrib (inc) */
#define SCREEN_DEBOF	111.												/* hubbell screening constant */
#define SCREEN_DEBON	135.6
#define TNULL			273.15
#define M_EulerC		0.5772156649
#define X_RANGE			1000												/* no of points in energy range 0<x<1 */
#define USTEPS			100													/* no of entries in Collifct table u=0..umax */
#define UMAX			5													/* umax: max angle for collifct-integral (<pi)*/
#define RSTEPS			5000												/* no of steps for the integral in Collifct */
#define THETARANGE		16.73												/* MCB: max polar photon angle : 0.01 rad * E0  MCB*/
#define DQ_STEPS		1000 												/* no of steps for q-integral in BH-intensity */
#define ESTEPS 			9													/* no of steps for primary beam spread folding (needs to beodd ! */
#define A_12C			12.
#define ALPHA2			(1/(137.0*137.))
#define DENS_DIAMOND	3.513

#define LOGTERM 		1
#define LOGFILE 		2													/* defines for logging, do not alter bit representation */
#define LOGERR			4
#define LOGBOTH 		LOGTERM|LOGFILE
#define LOGERROR		LOGERR |LOGFILE



/*******************************************************
		PROTOTYPES
********************************************************/

void	calc_collifct( double ), calc_latticevec( unsigned int ),
        sort( unsigned int ), calc_debyered( void );
double  collifunction( double, double ),
        collifctr( double, double, double ),
        bdiv_distr( double, double ), solveeq( double ),
        moliere_scattering(), sig_moliere (double, double ),
        energy_spread( double, double, double * );
void	cohbrems( double, double, double *, double * ),
        int_coh( double *, double * ),
        (*int_inc)( double *, double *, double, int ),
        int_BH( double *, double *, double, int ),
        int_Hub( double *, double *, double, int ),
        int_HubBd( double *, double *, double, int ),
        logging( int ), get_par ( FILE *, char, void * );
void input_parameter();



/*******************************************************
		GLOBAL VARS
********************************************************/

VEC 	lv[MAX_LATTICEVEC];													/* Lattice Vecs and structure factors */
double	S_2[MAX_LATTICEVEC];

double	g_2[MAX_LATTICEVEC], aps[MAX_LATTICEVEC], debye[MAX_LATTICEVEC], ff[MAX_LATTICEVEC], Irel[MAX_LATTICEVEC];
double	sig_colx2, sig_coly2, sig_colr2=0, sig_colas=0, sig_colli=0, collifct[USTEPS+1];
double	sigx, sigy, sigpx, sigpy, sig_mol=0, sig_bx, del_bx, sig_by, del_by;
double	umax, alpha, theta, phi, di_thickness;
double	E0MeV, sigE0MeV, E0, dist, collength, rkoll, C;
double	temperature, Zamo, Zcry, A3, Atot, Adebye;
int		type, max_latticevec, NMax;

FILE	*file_log;
char	input_file[50], logname[50], outname[50], comment[100], Deb_cmnt[100], _logbuf[200], date[100];




/**************************************************
	Main:
	Read input, calc parameters, call
	incoherent/coherent routines and
	write result in output file, Format of output:
	1000 (X_RANGE) rows with 7 coloumns
	photon energy [MeV] followed by 6 intensities
	total coherent, delta coherent = perp-para,
	incoherent crystal, incoherent amorphous,
	electron contrib crystal, electron amorphous
**************************************************/

int main(int argc, char **argv) {
	FILE	*fhdl;
	double	intsum[X_RANGE+1], intdif[X_RANGE+1], intinc[X_RANGE+1],
	        intamo[X_RANGE+1], intiel[X_RANGE+1], intael[X_RANGE+1];
	char	*s;
	int		i;


	time_t current_time;
	char date[26];


	/* get time, convert to ascii and strip trailing '\n' */
	current_time = time(NULL);
	strcpy( date, ctime(&current_time) );

	size_t len = strlen(date);
	if (len > 0 && date[len - 1] == '\n') {
		date[len - 1] = '\0';
	}


	for( i=1; i<=X_RANGE; i++ ) {
		intsum[i] = 0;
		intdif[i] = 0;
	}

	/**********************************************
	Set up files and read input parameter and
	calc several variables and collimation function
	**********************************************/

	if (argc == 2)
		strcpy (input_file, argv[1]);
	else {
		printf( "\nPlease enter the name of the input file (without the ending 'in')?\n" );
		scanf("%s", input_file);
	}
	if( (s=strstr(input_file,".in")) != NULL )								/* omit possible file-ending '.ein' */
		*s=0;																/* set filenames and open log file */
	strcpy( outname, input_file );
	strcat( outname, ".vec" );
	strcpy( logname, input_file );
	strcat( logname, ".log" );
	strcat( input_file, ".in" );
	if( (file_log = fopen( logname, "w" ) ) == NULL ) {
		fprintf( stderr, "\n\nFatal Error: couldn't open log file\n\n" );
		exit( -1 );
	}

	input_parameter();														/* read user input from parameter file */
	sprintf( _logbuf, "anb > Input file read, please wait approx 5 minutes...\nanb > Output is written in file %s and log file is %s\n",
	         outname, logname );
	logging( LOGBOTH );

	E0 = E0MeV/MASS_E;														/* calc some constants */
	if( max_latticevec > MAX_LATTICEVEC ) {
		max_latticevec = MAX_LATTICEVEC;
		sprintf( _logbuf, "anb > Warning: highest value for 'MAX_LATTICEVEC' = %i using\n", MAX_LATTICEVEC );
		logging( LOGFILE );
	}
	if( di_thickness > 0. )
		sig_mol = moliere_scattering();										/* add moliere angle deviation if requested, remind: sig_mol='plane-sigma' */
	sig_bx = sqrt( sig_mol*sig_mol + sigpx*sigpx*1E-6 );
	sig_by = sqrt( sig_mol*sig_mol + sigpy*sigpy*1E-6 );
	sprintf( _logbuf, "anb > Mean electron divergence: sig_bx/y: %lf\t%lf mrad\n", sig_bx*1E3, sig_by*1E3 );
	logging( LOGFILE );

	if( dist == 0. ) {													/* col angl var = beam div + spot var */
		C = E0 * M_2PI;														/* no collimation if coll-dist == 0 */
		umax = C;															/* for comparison with MCB use C = THETARANGE */
		if( type == 2 ) type = 1;
	}										/* if divergent hubbell choosen for uncoll case use 'normal' hubbell !! */
	else {
		C = rkoll*1E-3/(dist+collength) * E0;								/* eff colli => add colli length to distance */
		umax = UMAX*C;														/* eff spherical divergence sigma, remind: sig_mol='plane-sigma' */
		sig_colx2 = ( sigx*sigx*1E-6/(dist*dist) + sigpx*sigpx*1E-6 + sig_mol*sig_mol ) * E0*E0;
		sig_coly2 = ( sigy*sigy*1E-6/(dist*dist) + sigpy*sigpy*1E-6 + sig_mol*sig_mol ) * E0*E0;
		sig_colr2 = sig_colx2;
		sig_colas = 1./(2.*sig_coly2) - 1./(2.*sig_colx2);
		sig_colli = ( sig_colx2 + sig_coly2 );
	}

	if( type == 0 ) {
		int_inc = int_BH;
		sprintf( _logbuf, "anb > Calculating incoherent part via Bethe Heitler and trivial Electron contrib\n" );
	}
	if( type == 1 ) {
		if( sig_colli==0 )	int_inc = int_Hub;
		else				int_inc = int_HubBd;
		sprintf( _logbuf, "anb > Calculating incoherent part via Hubbell and Electron contrib from [Owens]\n" );
	}
	if( type > 1 ) {
		sprintf( _logbuf, "anb > FATAL ERROR: process type doesn't exist\n" );
		logging( LOGERROR );
		exit(-1);
	}
	logging( LOGFILE );

	sprintf( _logbuf, "anb > Collimation angle (E0*theta_c): %lg\tcolli-divergence: %lg\n", C, sig_colli );
	logging( LOGFILE );
	A3 = pow( A, -3. );
	del_bx  = sig_bx>0. ? 2*BEAMXY_REGION*sig_bx/BEAMX_STEPS : 100;			/* beam divergence -> crystal orientation variation */
	del_by  = sig_by>0. ? 2*BEAMXY_REGION*sig_by/BEAMY_STEPS : 100;			/* sig_px/py = 0 -> no beam div or beamx/y variation */

	calc_collifct( sig_colli*M_SQRT1_2 );									/* calc collifunction dependent on sig_colli */
	calc_latticevec( 10000 );												/* calc and sort contributing lattice vecs lv[], S_2[] */

	int_coh( intsum, intdif );												/* calc coherent and incoherent intensities */
	int_inc( intinc, intiel, Zcry, DEBYE_ON  );								/* inc->incoherent crystal radiator */
	int_inc( intamo, intael, Zamo, DEBYE_OFF );								/* inc->incoherent amorphous radiator */

	if( sigE0MeV > 0.0 ) {
		energy_spread( E0, sigE0MeV, intsum );								/* fold spectra with spreaded beam energy E0 */
		energy_spread( E0, sigE0MeV, intdif );
		energy_spread( E0, sigE0MeV, intinc );
		energy_spread( E0, sigE0MeV, intamo );
	}

	/*****************************************************
	Write result in output file, contents:
	energy(MeV) int_sum int_dif int_inc(hub) int_inc(elec)
	*****************************************************/
	if( (fhdl=fopen(outname,"w")) == NULL ) {
		sprintf( _logbuf, "Fatal Error: couldn't open output file\n" );
		logging( LOGERROR );
		exit( -1 );
	}

	for( i=1; i<=X_RANGE; i++ )
		fprintf( fhdl, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", (float)((float)i/1000.*E0MeV),
		         (float)intsum[i], (float)intdif[i], (float)intinc[i], (float)intamo[i], (float)intiel[i], (float)intael[i] );
	fclose(fhdl);

	/* get time, convert to ascii and strip trailing '\n' */
	current_time = time(NULL);
	strcpy( date, ctime(&current_time) );
	date[strlen(date)-1] = 0;
	sprintf( _logbuf, "anb > intensities I written in file %s\n\
contents: photon energy(MeV) Isum Idif Iinccrystal Iincamorphous Ieleccry Ielecamo\n\
use paw, xmgr or plotdata to calc Iinc,amo=(Iinc+Ielec)(cry,amo), Irel_exp=(Isum+Iinc)/Iamo, Pol=Idif/(Isum+Iinc)\n\
anb > Program finished (%s)\n", outname, date );
	logging( LOGBOTH );
	fclose( file_log );

	return 0;
}



/*************************************************
	read parameters from input file
*************************************************/

void input_parameter() {
	FILE	 *con, *file_debye;
	int 	 j = 0;
	char	 fnam_debye[100];

	if( (con = fopen(input_file,"r")) == NULL) {
		sprintf( _logbuf, "Fatal Error > unable to open input file: '%s'.\n", input_file );
		logging( LOGERROR );
		exit( -1 );
	}

	sprintf( _logbuf, "\n\n********* ANB Version 2.0 (%s) ********* Alexander Natter/11/98\n\n\
input > Reading input file %s\n", date, input_file );
	logging( LOGBOTH );

	get_par( con, 's', comment ); 												/* comment line for log-file */
	sprintf ( _logbuf, "\ninput > Comment:\n%s\n", comment);
	logging( LOGBOTH );
	get_par( con, 'd', &E0MeV );												/* energy of the incident electron beam */
	fprintf ( file_log, "energy = %lf MeV\n", E0MeV );
	get_par( con, 'd',& sigE0MeV ); 											/* energy spread of incident beam */
	fprintf ( file_log, "energy spread = %lf MeV\n", sigE0MeV );
	get_par( con, 'd', &theta );												/* lattice orientation angle theta,*/
	fprintf( file_log, "Theta = %lf rad\n", theta);
	get_par( con, 'd', &alpha );												/* alpha in [rad] */
	fprintf( file_log, "Alpha = %lf rad\n", alpha);
	get_par( con, 'd', &phi );  												/* phi in [rad], typically Phi=Pi/4 */
	fprintf( file_log, "phi = %lf rad\n", phi);
	get_par( con, 'd', &sigx ); 												/* h. beam divergence in [mm] */
	fprintf( file_log, "sigma_x = %lf mm\n", sigx);
	get_par( con, 'd', &sigy ); 												/* v. beam divergence in [mm] */
	fprintf( file_log, "sigma_y = %lf mm\n", sigy);
	get_par( con, 'd', &sigpx ); 												/* h. angle divergence of the beam in [mrad]  */
	fprintf( file_log, "sigma_px = %lf mrad\n", sigpx);
	get_par( con, 'd', &sigpy );
	fprintf( file_log, "sigma_py = %lf mrad\n", sigpy);
	get_par( con, 'd', &di_thickness ); 										/* thickness of the diamond in [mm] */
	fprintf( file_log, "Di-thickness = %lf mm\n", di_thickness );
	get_par( con, 'd', &dist );												/* distance target-collimator in [m] */
	fprintf(file_log,"Distance target - collimator = %lf m\n", dist );
	get_par( con, 'd', &collength ); 											/* length of collimator in [m] */
	fprintf(file_log,"Length of collimator = %lf m\n", collength );
	get_par( con, 'd', &rkoll );
	fprintf(file_log,"collimator radius = %lf mm,\n", rkoll ); 			 	/* radius of collimator [mm] */
	get_par( con, 'i', &type );												/* BH or Hub inc/amo Bremsstrahlung */
	fprintf (file_log, "%scoherente Bremsstrahlung (type:%i)\n", type>=9 ? "" : "in", type );
	get_par( con, 'i', &max_latticevec );  									/* Max Lattice vector index in [HKL] */
	fprintf (file_log, "Max Lattice vector index in Lattice vector table i=%i in [0..%i]\n",
	         max_latticevec,  MAX_LATTICEVEC );
	get_par( con, 'd', &Zcry ); 												 	/* Z of cyrstal radiator */
	fprintf( file_log, "using Z=%3.0lf for 'crystal' radiator\n", Zcry );
	get_par( con, 'd', &Zamo ); 												 	/* Z of amorphous radiator*/
	fprintf( file_log, "using Z=%3.0lf for 'amorphous' radiator\n", Zamo );
	fclose( con );

	fflush( file_log );
}



/*************************************************
	approximate treatment of beam divergence,
	multiple scattering, finite beam spot size
	in case of collimation on the incoherent part,
	here only on the colli-cutoff: x_c
*************************************************/


void calc_collifct( double sig ) {										/* calc collifct table sig=sig_angle*E0 */
	double	a, g, r, rstep, r1, r2, sum=0,
	                                u, ustep=umax/USTEPS;
	int		i=0;

	sprintf( _logbuf, "calc_collifct > calculating collimation function\n" );
	logging( LOGBOTH );

	if( sig  == 0. )														/* no beamdiv -> sharp colli */
		for( u=0.; u<umax; u+=ustep, i++ )
			collifct[i] = u<C ? 1 : 0;
	else
		for( u=0.; u<umax; u+=ustep, i++ ) {								/* loop over all theta angles */
			if( u>0. ) {
				rstep = ( C+u - fabs(C-u) )/RSTEPS;
				r1	 = fabs(C-u) + rstep;
				r2	 = C+u - rstep;

				sum   = collifctr( u, r1, sig ) * rstep/2;					/* perform integral */
				for( r = r1; r <= r2; r += rstep )							/* extended closed trapezoidal rule */
					sum += collifctr( u, r, sig ) * rstep;
				sum  += collifctr( u, r2, sig ) * rstep/2;
			}
			if( C > u )
				sum += 1-exp( -(C-u)*(C-u)/(2*sig*sig));					/* step fct */
			collifct[i] = sum;
		}

	sprintf( _logbuf, "calc_collifct > done\n" );
	logging( LOGFILE );
}


double collifctr( double u, double r, double s ) {
	return  r * exp( r*r/(-2.*s*s))/(s*s) * acos((r*r+u*u-C*C)/(2.*r*u+1e-10))/M_PI;
}


double collifunction( double x, double xd ) {
	double	u, ustep=umax/USTEPS, f1, f2;
	int	i;

	if( x==0. || x>xd )  return 0;
	u  = sqrt( (xd/x-1.)/(1.-xd) );										 	/* calc u(x), use table and interpolate */
	if( u>= umax ) return 0;

	i  = (int)(u/ustep);
	f1 = collifct[i];
	f2 = collifct[i+1];

	return	f1 + (f2-f1)/ustep * (u-i*ustep);
}



/*************************************************
	approximate treatment of beam divergence,
	multiple scattering, finite beam spot size
	in case of collimation on the incoherent part
*************************************************/

void int_HubBd( double *iinc, double *iiel, double Z, int debye_switch ) {
	double	norm=0, weight, U0=C, Cm, U, f=0;
	double	rho, drho=sig_colli/11., sg2 = 1/(2*sig_colli*sig_colli), ph, dph=M_2PI/31., phc;
	double	ii[X_RANGE+1], ie[X_RANGE+1], ir[X_RANGE+1], er[X_RANGE+1];
	int		i;

	for( i=1; i<X_RANGE; i++ )	iinc[i] = iiel[i] = 0;
	sprintf( _logbuf, "int_HubBd > folding electron divergence and hubbel intensity %s debye-factor\n",
	         debye_switch ? "with" : "without" );
	logging( LOGBOTH );

	for( rho = drho/2.;  rho < U0;  rho+=drho )								/* perform integration over beam divergence */
		for( ph = dph/2.;  ph < M_2PI;	ph += dph ) {							/* 'inside' the collimator */
			C  = sqrt( U0*U0 - pow(rho*sin(ph),2.) ) - rho*cos(ph);
			weight = dph*drho/M_PI * exp( -rho*rho*sg2 ) * rho / (M_2PI*sig_colli*sig_colli);
			norm  += weight;
			int_Hub( ii, ie, Z, debye_switch ); 								/* calc Hubbel intensity with Colli C=U(rho,ph)*/
			for( i=1; i<X_RANGE; i++ ) {								  		/* add contributions */
				iinc[i] += ii[i]*weight;
				iiel[i] += ie[i]*weight;
			}
		}

	C = U0;
	drho = sig_colli/6;
	for( rho = U0+drho/2.; rho < 3*U0; rho+=drho ) {						/* perform integration over beam divergence */
		/* 'outside' the collimator */
		phc = M_PI-asin(U0/rho);
		dph = (M_2PI-2.*phc)/10.;
		for( ph = phc+dph/2.; ph <= M_2PI-phc; ph += dph  ) {
			C  =  sqrt( U0*U0 - pow(rho*sin(ph),2.) ) - rho*cos(ph);
			Cm = -sqrt( U0*U0 - pow(rho*sin(ph),2.) ) - rho*cos(ph);
			weight = dph*drho/M_PI * exp( -rho*rho*sg2 ) * rho / (M_2PI*sig_colli*sig_colli);
			norm += weight;
			int_Hub( ii, ie, Z, debye_switch );
			for( i=1; i<X_RANGE; i++ ) {
				iinc[i] += ii[i]*weight;
				iiel[i] += ie[i]*weight;
			}
			C = Cm;
			int_Hub( ii, ie, Z, debye_switch );
			for( i=1; i<X_RANGE; i++ ) {
				iinc[i] -= ii[i]*weight;
				iiel[i] -= ie[i]*weight;
			}
		}
	}

	for( i=1; i<X_RANGE; i++ ) {
		iinc[i] /= norm;
		iiel[i] /= norm;
	}
	C = U0;																	/* restaurate global variable */
	sprintf( _logbuf, "int_HubBd > done\n" );
	logging( LOGBOTH );
}



/*************************************************
	calc the collimated incoherent contribution
	via the Hubbel Xsec and the electron brems
	part based on Matthews/Owens -- 'ideal' case
*************************************************/

void int_Hub( double *iinc, double *iiel, double Z, int debye_switch ) {	/* calc int_inc from HUBBEL and shell electrons */
	double	x, emx, dmy, fc, Z3, b, at, eps, psieps,
	        Mv, M1, psi1, psi2, v, v2, screening;
	int		i;

	fc = C*C/(1+C*C);														/* assume same angle dep for electron as nucl brem */
	Z3 = pow( Z, 1./3. );
	v  = 1./(1.+C*C);
	v2 = v*v;
	if( debye_switch==DEBYE_ON )	screening = SCREEN_DEBON;
	else							screening = SCREEN_DEBOF;
	sprintf( _logbuf, "int_inc > calculating inc hubbel + electron intensities\
int_inc > Z: %lg\tDebye Factor `mode`:\tColli: %lf\t", Z, debye_switch, C );
	logging( 0/*LOGFILE*/ );
	for( i=1; i<X_RANGE; i++ ) {											/* loop over photon energies */
		x  = (double)i/X_RANGE;
		emx= 1. - x;
		b  = 2.*Z3*E0/screening * emx/x;									/* Hubbel contrib */
		at = atan(C*C*b/(1+C*C+b*b)) * 2./b;
		dmy= pow( (x/(2.*E0*emx)), 2. );
		Mv = -log( dmy + pow( (v*Z3/screening), 2. ) );
		M1 = -log( dmy + pow( (  Z3/screening), 2. ) );
		psi1 = 2. * ( 1 + M1 - (1+Mv)*v - at );
		psi2 = -40./3.*v*v2 + 18.*v2 - (8./(b*b)+6.)*v + 8./(b*b) + 4./3.
		       + (4.*v*v2-6.*v2)*Mv + 2.*M1 - 6./(b*b)*(Mv-M1+2./3.*at);
		iinc[i] = ( (1+emx*emx)*psi1  - 2./3*emx*psi2 ); 					/* hubbell intensity */
		eps     = 100./(E0*Z3*Z3) * x/emx;									/* see Matthews and owens */
		if( eps<=0.88 ) {													/* for eps<0.88, hence about x<0.95  */
			eps     = 0.88 - eps;
			psieps  = 19.70 + eps * ( 4.177 + eps * ( -3.806
			                          + eps*( 31.84 + eps*( -58.63 + eps*40.77 ))));
		} else
			psieps  = 19.9 - 4.*log(eps);		   							/* for eps>0.88 */
		psi1    = psieps - 4 - 8/3 * log(Z);
		psi2    = psi1 + 0.6666666;
		iiel[i] = ( (1+emx*emx)*psi1  - 2./3*emx*psi2 ) * fc/Z;				/* shell electron intensity */
	}
	iinc[X_RANGE] = 0;
	iiel[X_RANGE] = 0;
	sprintf( _logbuf, "int_inc > done\n" );
	logging( 0/*LOGFILE*/ );
}



/*************************************************
	calc the collimated incoherent contribution
	via the bethe-heitler xsec
	(collimation effect treated with the
	asymptotic angular BH-dependence
*************************************************/

void int_BH( double *iinc, double *iiel, double Z, int debye_switch ) {	/* bethe heitler intensity */
	const double	psi1c = 13.79,  psi2c = 13.12,							/* approx of q-integral, s.below */
	                psi1e=   4.05,  psi2e =  3.94;							/* shell electron contrib */
	double 	x, emx, del, q, q2, q3, dq, psi1, psi2,
	        debfact = 1., ff, dmy, fc=C*C/(1+C*C);
	int		i;

	sprintf( _logbuf, "int_BH > starting to calculate the incoherent (Bethe Heitler) intensities\n\
int_inc > Z of radiator: %lg\t%s using Debye Factor\tColli: %lf\n", Z, debye_switch ? "" : "not", C );
	logging( LOGBOTH );
	for( i=1; i<X_RANGE; i++ ) {											/* loop over photon energies */
		x  = (double)i/1000.;
		emx= 1. - x;
		del= x/(2.*E0*emx);
		dq = (1.-del)/(double)DQ_STEPS;
		psi1 = 0.;
		psi2 = 0.;
		for( q=del; q<1.; q+=dq ) {										/* perfom q-integral */
			q2 = q*q;
			q3 = q*q2;
			if( debye_switch==DEBYE_ON )
				debfact = 1. - exp(-Adebye*q2);
			ff = 1 - ( 0.2283+1.8359*exp(-10528*q2)+1.8119*exp(-4678*q2)
			           +1.5809*exp(-239*q2)+0.5426*exp(-27116*q2) )/6.;
			dmy= debfact*ff*ff/q3;
			psi1 += dmy*(q-del)*(q-del);
			psi2 += dmy*(q2+del*del*(3.-6.*log(q/del)-4.*del/q));
		}
		psi1 = 4.    + 4.*psi1*dq;
		psi2 = 10./3.+ 4.*psi2*dq;
		iinc[i] = ( (1+emx*emx)*psi1  - 2./3*emx*psi2 )  * fc; 				/* Bethe Heitler intensity */
		iiel[i] = ( (1+emx*emx)*psi1e - 2./3*emx*psi2e ) * fc; 				/* assume same angle dep as nucl brem */
	}
	iinc[X_RANGE] = 0;
	iiel[X_RANGE] = 0;
	sprintf( _logbuf, "int_BH > done\n" );
	logging( LOGFILE );
}



/*************************************************
	approximate treatment of beam divergence,
	multiple scattering, finite beam spot size
	in case of collimation on the coherent part,
	here effect of crystal angles only
*************************************************/

void int_coh( double *intsum, double *intdif ) {
	VEC		g;
	double	al, th, beam_x, beam_y,
	        phi_beam, th_beam, discr, weight, norm=0;
	double	isum[X_RANGE+1], idif[X_RANGE+1];
	int		i, l;

	sprintf( _logbuf, "int_coh > calculating the coherent (sum, dif) intensities\n" );
	logging( LOGBOTH );
	for( l=0; l<max_latticevec; l++ ) {									/* loop over all lattice-vecs and add coh contrib*/
		g       = lv[l];
		g_2[l]  = v_sprd( g, g ) * GUNIT*GUNIT;								/* calc kin factors and mom transfers etc per lattice vec */
		aps[l]  = ( (g.y*g.y - g.z*g.z)*cos(2*phi) + 2*g.y*g.z*sin(2*phi) ) * GUNIT*GUNIT;
		debye[l]= exp(-Adebye*g_2[l]);
		ff[l]   = 0.2283 + 1.8359*exp(-10528*g_2[l]) + 1.8119*exp(-4678*g_2[l])
		          + 1.5809*exp(-239*g_2[l]) + 0.5426*exp(-27116*g_2[l]);
		ff[l]   = 1.-ff[l]/6.;
	}
	/* variate crystal orientation due to beam divergence */
	for( beam_x = -BEAMXY_REGION*sig_bx; beam_x <= BEAMXY_REGION*sig_bx; beam_x+=del_bx )/* loop over all possible beam directions */
		for( beam_y = -BEAMXY_REGION*sig_by; beam_y <= BEAMXY_REGION*sig_by; beam_y+=del_by ) {
			th_beam = sqrt( beam_x*beam_x + beam_y*beam_y );
			phi_beam= ( beam_x !=0.0 ? atan( beam_y/beam_x ) : 0.0 );
			discr   = theta*theta + th_beam*th_beam + 2.*theta*th_beam*cos(phi_beam);
			if( discr <= 0. ) {
				sprintf( _logbuf, "inc_coh > FATAL ERROR: beam divergence too large: phi:%lg\tth:%lg\tdiscr:%lg\n",
				         phi_beam, th_beam, discr );
				logging( LOGERROR );
				continue;
			}
			th = sqrt( discr );												/* calc 'rotated' crystal (al,th) due to th/phi_beam */
			al = alpha + asin( sin(phi_beam)*th_beam/th );					/* weight intensity beam divergence distribution */
			weight =	( sig_bx>0. ? exp( -(beam_x*beam_x)/(2*sig_bx*sig_bx) ) : 1 ) *
			            ( sig_by>0. ? exp( -(beam_y*beam_y)/(2*sig_by*sig_by) ) : 1 );
			norm  += weight;
			sprintf( _logbuf, "int_coh > calculating intensities for alpha:%lg\ttheta:%lg\tphib:%lg\tthb:%lg\tweight:%lg\n",
			         al, th, phi_beam, th_beam, weight );
			logging( LOGFILE );
			cohbrems( al, th, isum, idif );									/* calc intensities (para+-perp) with collimation and beamdiv */
			for( i=1; i<=X_RANGE; i++ ) {									/* and add contributions */
				intsum[i] += isum[i]*weight;
				intdif[i] += idif[i]*weight;
			}
		}
	sprintf( _logbuf, "int_coh > sum of theta/phi-crystal weight: %lg\n", norm );
	logging( LOGFILE );
	for( i=1; i<=X_RANGE; i++ ) {
		intsum[i] /= norm;
		intdif[i] /= norm;
	}
	sprintf( _logbuf, "int_coh > done\n" );
	logging( LOGFILE );
}



/**********************************************
	calc 'ideal' coherent contribution for
	all energies and lattice vectors
**********************************************/

void cohbrems( double alpha, double theta, double *isum, double *idif ) {
	int		i, l;
	double	gl, gt_2, x, del, xd, xc, Gg, dmy,
	        psi1, psi2, psi3;
	VEC		g;

	for( i=1; i<=X_RANGE; i++ ) {
		isum[i] = 0;
		idif[i] = 0;
	}

	for( l=0; l<max_latticevec; l++ ) {									/* loop over all lattice vectors */
		g   = lv[l];
		dmy = sin(theta) * ( g.y*cos(alpha) + g.z*sin(alpha) );
		gl  = ( cos(theta)*g.x + dmy ) * GUNIT;
		gt_2 = ( g.y*g.y + g.z*g.z + pow(sin(theta)*g.x,2.) - dmy*dmy ) * GUNIT*GUNIT;
		xd  = 2.*E0*gl / ( 1. + 2.*E0*gl );
		xc  = xd/(1.+C*C*(1.-xd));
		dmy = M_2PI*ff[l] / (g_2[l]*gl*gl);
		Gg  = A3/8.  * dmy*dmy * S_2[l] * debye[l];

		for( i=1; i<X_RANGE; i++ ) {										/* loop over all photon energies */
			x   = (double)i/1000.;
			if( x>xd ) {
				i=X_RANGE;
				continue;
			}
			del = x/(2.*E0*(1.-x));
			psi1 =  4 * Gg * del	* gt_2 * gl*gl;							/* calc psi - fcts */
			psi2 = 24 * Gg * del*del * gt_2 * (gl-del);
			psi3 = -4 * Gg * del*del*del * aps[l];							/* calc intensities (para+-perp) and use collifct */
			isum[i] += ( (1.+(1.-x)*(1.-x))*psi1 - 2./3.*(1.-x)*psi2 ) * collifunction( x, xd );
			idif[i] -= 2*(1.-x)*psi3 * collifunction( x, xd );
		}
	}
}



/**********************************************
	perform a 'convolution' of the spectra
	with a beam energy distribution (spread)
**********************************************/

double energy_spread( double E0, double sigE0MeV, double *I ) {
	double	Ie[X_RANGE+1], weight[ESTEPS+2];
	double	x, x0, x1, sig2E0;
	double	Emin, Emax, Estep, Es, sumweight=0;
	int	i, i0, Eidx;

	sprintf( _logbuf, "energy_spread > fold intensity with primary beam energy spread \n" );
	logging( LOGBOTH );
	for( i=1; i<X_RANGE; i++ )	Ie[i]=0.;
	Emin   = E0 - 2*sigE0MeV/MASS_E;
	Emax   = E0 + 2*sigE0MeV/MASS_E;
	Estep  = (Emax-Emin)/ESTEPS;
	sig2E0 = 2 * sigE0MeV*sigE0MeV/(MASS_E*MASS_E);
	sprintf( _logbuf, "energy_spread > E0:%lf\tsigE0MeV:%lf\tEstep:%lf\tEmin:%lf\tEmax:%lf\n",
	         E0, sigE0MeV, Estep, Emin, Emax );
	logging( LOGFILE );

	for( Eidx = 0; Eidx<=ESTEPS; Eidx++ ) {									/* calc weights of beam spreaded energy Es */
		Es = Emin + Eidx*Estep;
		weight[Eidx] = exp( -(Es-E0)*(Es-E0)/sig2E0 );
		sumweight += weight[Eidx];
		sprintf( _logbuf, "Eidx:%i\tEs:%lf\tweight:%lf\n", Eidx, Es, weight[Eidx] );
		logging( LOGFILE );
	}

	for( i=1; i<X_RANGE; i++ ) {
		x = (double)i/X_RANGE;
		for( Eidx=0;	Eidx<=ESTEPS; Eidx++ ) {
			x1 = x*E0/Es;
			i0 = (int)(X_RANGE*x1);
			x0 = i0/X_RANGE;
			Ie[i] += ( I[i0] + (I[i0+1]-I[i0]) * (x1-x0) ) * weight[Eidx];
		}
	}

	for( i=1; i<X_RANGE; i++ )
		I[i] = Ie[i]/sumweight;
}



/**********************************************
	calc, sort and store lattice vectors
	and dependend variables -> faster
**********************************************/

void calc_latticevec( unsigned int lvmax ) {
	int	mod, S2, lvidx=0;
	VEC	g;
	double	xd, gl, g_2, gt_2, ff, I, I1, I2, hklmax;
	double	Gg, del, psi1, dmy;

	hklmax = (double)((int)((pow((double)lvmax,1./3.)-1.)/2.));
	lvmax  = pow(2*hklmax+1,3.);
	sprintf( _logbuf, "calc_latticevec > calculating relative contributions of %i lattice vectors\n", lvmax );
	logging( LOGBOTH );														/* check for approximate expected ("contributing") lattice */
	sprintf( _logbuf, "calc_latticevec > checking %i lattice vectors, maximal miller index:%lf\n", lvmax, hklmax );
	logging( LOGFILE );
	if( 3*(int)pow(hklmax-.5,5./2.) > MAX_LATTICEVEC ) {					/* vectors lv_contrib = 3*(hklmax-1/2)^(5/2)*/
		sprintf( _logbuf, "calc_latticevec > Error: Too many lattice vects requested\nMaximal %i possible\n", MAX_LATTICEVEC );
		logging( LOGERROR );
		exit(-1);
	}

	for( g.x=-hklmax; g.x<=hklmax; g.x++ )
		for( g.y=-hklmax; g.y<=hklmax; g.y++ )
			for( g.z=-hklmax; g.z<=hklmax; g.z++ ) {
				mod = abs((int)g.x)%2+abs((int)g.y)%2+abs((int)g.z)%2;
				if( mod==0 && (int)(g.x+g.y+g.z)%4==0 )  S2 = 64; 						/*Structurefactor*/
				else if( mod==3 )  S2 = 32;
				else S2 = 0;
				if( S2 > 0 ) {
					dmy = sin(theta) * ( g.y*cos(alpha) + g.z*sin(alpha) );
					gl  = ( cos(theta)*g.x + dmy ) * GUNIT;
					xd  = 2.*E0*gl / ( 1. + 2.*E0*gl );
					if( xd>0. && xd<1. ) {
						g_2  = v_sprd( g, g ) * GUNIT*GUNIT;
						gt_2 = ( g.y*g.y + g.z*g.z + pow(sin(theta)*g.x,2.) - dmy*dmy ) * GUNIT*GUNIT;
						ff	= 0.2283 + 1.8359*exp(-10528*g_2) + 1.8119*exp(-4678*g_2)
						      + 1.5809*exp(-239*g_2) + 0.5426*exp(-27116*g_2);
						ff	= 1.-ff/6.;
						del  = xd/(2.*E0*(1.-xd));										/* @ x=xd : psi2 always zero !!! */
						psi1 =  A3*M_2PI*M_2PI/2. * ff*ff * S2 * exp(-Adebye*g_2)/del/g_2;
						I	=  (1.+(1.-xd)*(1.-xd)) * psi1;
						Irel[lvidx] = I;
						S_2[lvidx]  = S2;
						lv[lvidx++] = g;
						if( lvidx>MAX_LATTICEVEC ) {
							sprintf( _logbuf, "\n\ncalc_latticevec > Error: MAX_LATTICEVEC exceeded (lvidx:%i)\n\n", lvidx );
							logging( LOGERROR );
							exit(-1);
						}
					}
				}
			}

	sort( lvmax=lvidx );
	sprintf( _logbuf, "calc_latticevec > sorting %i contributing lattice vectors\n", lvmax );
	logging( LOGFILE );
	for( lvidx=0; lvidx<lvmax; lvidx++ ) {
		gl = ( cos(theta)*lv[lvidx].x + sin(theta)*(lv[lvidx].y*cos(alpha)+lv[lvidx].z*sin(alpha)) ) * GUNIT;
		xd = 2.*E0*gl / ( 1. + 2.*E0*gl );
		sprintf( _logbuf, "%lf\t%lf\t%lf\tg:", Irel[lvidx], S_2[lvidx], xd );
		v_out( file_log, _logbuf, lv[lvidx] );
	}
	sprintf( _logbuf, "calc_latticevec > using the %i most contributing lattice vectors out of %i for further calculation\n",
	         max_latticevec, lvmax );
	logging( LOGFILE );
}



/**********************************************
	integrate (average) moliere angle
	over the radiator-thickness
**********************************************/

double moliere_scattering() {
#define MSSTEPS	1E2
	double _sig_moliere (double s, double pe);
	double s, sm=0, step=di_thickness/MSSTEPS;

	sprintf( _logbuf, "moliere > starting to calculate moliere variance, sigmoliere(0,thickness): %lg (mrad)\n",
	         sig_moliere(di_thickness,E0MeV) * 1E3 );
	logging( LOGFILE );
	for( s=step;  s<di_thickness;  s+=step )									/* average moliere angles */
		sm += sig_moliere( s, E0MeV );
	sprintf( _logbuf, "moliere > mean sigma_moliere: %lg mrad\n",sm/MSSTEPS * 1E3 );
	logging( LOGFILE );
	return( sm/MSSTEPS );
}



/**********************************************
	mean angle divergence due to
	multiple scattering (moliere theorie)
	s : distance travelled in medium, unit: mm
	pe: electron momentum, unit: MeV
	approx used : beta(electron) == 1
**********************************************/

double sig_moliere (double s, double pe) {
	double	om, b, B, thsq;

	s /= 10.;
	om = 7800. * (1.+Zcry)*pow(Zcry,1./3.) * DENS_DIAMOND*s
	     / ( A_12C * (1.+3.35*(Zcry*Zcry*ALPHA2) ) );
	b    = log(om) + 1. - 2.*M_EulerC;
	B    = solveeq( b );
	thsq = 0.157 * Zcry*(1.+Zcry) * DENS_DIAMOND*s / (A_12C*pe*pe);

	return	sqrt( thsq*(B-1.2) );
}



/**********************************************
	solve B - log(B) = b
	using regula falsi method
**********************************************/

double solveeq ( double b ) {
	double	f_r, f_l, f_next, re=40.0, li=1.2, next;

	f_r = re - log(re) - b;
	f_l = li - log(li) - b;
	while( re-li > 1E-9 ) {
		next   = (re+li)/2;
		f_next = next - log(next) - b;
		if(f_next<=0.0) {
			li =next;
			f_l=f_next;
		} else {
			re =next;
			f_r=f_next;
		}
	}

	return	next;
}



/**********************************************
	sort lattice vectors depending on the
	relative contribution
	straight insertion method,
	see numerical recipes: piksrt
**********************************************/

void sort( int unsigned n ) {
	int		i, j;
	double	Ii, Si;
	VEC		gi;

	for( j=1; j<n; j++ ) {
		Ii = Irel[j];
		Si = S_2[j];
		gi = lv[j];
		i  = j-1;
		while( i >= 0 && Irel[i] < Ii ) {
			Irel[i+1] = Irel[i];
			S_2[i+1]  = S_2[i];
			lv[i+1]   = lv[i];
			i--;
		}
		Irel[i+1] = Ii;
		S_2[i+1]  = Si;
		lv[i+1]   = gi;
	}
}



/************************************************
	read one par from file, skip comment lines
	types as in `scanf`, but without `%` char
	e.g.: s=string, whole line
	d:double, f:float, c:char, i:int
************************************************/

void	get_par ( FILE *file, char type, void *data ) {
#define	BUFSIZE		500
	char	tmp[BUFSIZE], *form;
	memset (tmp, 0, BUFSIZE);

	do {
		if(fgets (tmp, BUFSIZE, file) == NULL)
			break;
		if( *tmp!='#' && *tmp!='\n' && *tmp!= ' ' && *tmp!='\t')
			break;
	} while( !feof(file) );
	if( *tmp == 0 ) {
		sprintf( _logbuf, "get_par > FATAL ERROR: couldn't get parameter" );
		logging( LOGERROR );
		exit( -1 );
	}

	switch( type ) {
		case 's':
			strcpy( data, tmp );
			break;
		case 'c':
			*(char * )data = *tmp;
			break;
		case 'f':
			sscanf( tmp, "%f",  data );
			break;
		case 'd':
			sscanf( tmp, "%lf", data );
			break;
		case 'i':
			sscanf( tmp, "%i",  data );
			break;
	}
}



/**********************************************
	print and log program messages
**********************************************/

void logging( int what ) {
	if( LOGTERM == (LOGTERM&what) ) {
		printf( _logbuf );
		fflush( stdout );
	}
	if( LOGFILE == (LOGFILE&what) )
		fprintf( file_log, _logbuf );
	if( LOGERR  == (LOGERR&what) )
		fprintf( stderr, _logbuf );
}

