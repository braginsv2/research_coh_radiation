
/*****************************************
 Primitive 3D double vektor - routines
*****************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define V_ZERO     1E-8
#define V_EQU      1E-10


typedef struct
{
  double x,y,z;
}
VEC;

typedef struct
{
  VEC x,y,z;
}
MATRIX;


/* Bem: die Fkt-Parameter werden der Reihe nach mit a,b,c... */

VEC v_set(double, double, double);	/*v=(a,b,c)*/
VEC v_setpol(double,double,double);/*v=(r,th,ph)*/
VEC v_add(VEC, VEC);		/*v=a+b    */
VEC v_sub(VEC, VEC);		/*analog   */
VEC v_smul(double, VEC);		/*v=a*b    */
VEC v_kreuz(VEC, VEC);		/*v=a^b    */
VEC v_lgs(MATRIX,VEC);		/*Av=b ->v */
VEC v_scale(VEC, VEC);		/*vi=ai*bi */
VEC v_scadd(VEC, VEC, VEC);	/*vi=ai*bi+ci*/
VEC v_mmul(MATRIX, VEC);		/*v=Ab     */
VEC v_normal(VEC);		/*v=a/|a|  */
VEC v_minus(VEC);			/*v=-a     */
VEC v_perpend(VEC, VEC);		/*senkr-Komponente*/
VEC v_parallel(VEC, VEC);		/*paralell-Komp.*/
VEC v_zrotate(VEC,VEC);		/*koordtrafo:b->b' wo a=z-axis*/
VEC v_polrot(double,double,VEC);	/*rotiere vec um theta,phi in sph. koord*/
VEC v_polrotinv(double,double,VEC);/*dito, but inverse operation*/


double v_sprd(VEC,VEC);		/*w=a.b    */
double v_sqr(VEC);			/*w=a.a	   */
double v_det(MATRIX);			/*w=|A|    */
double v_norm(VEC);		/*w=|v|    */
double v_vdet(VEC,VEC, VEC);		/*w=|a,b,c|*/
double v_transverse(VEC);		/*sqrt(x^2+y^2)*/


void v_out(FILE*, char* , VEC);	/*print '%s %vec' in file. Bei file=0->stdout und file neg warten auf Taste*/
void v_cart2pol( VEC, double*, double* , double*);	/* cart vec a -> r, theta, phi (pol coord) */


int  v_equ(VEC,VEC);			/*(a-b)^2<V_ZERO  */
int  v_lin_rect();			/*v_lin_rect( VEC p,l, q,m,n &r )  Schnitt gerade,ebene*/ 

VEC  v_null;
