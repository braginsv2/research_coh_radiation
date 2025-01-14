/*****************************************
 Primitive 3D double vektor - routines
 (C) Frank A. Natter, Uni. Tuebingen
 Version: 6 Apr 1998
*****************************************/
#include "vec.h"

VEC  v_null = { 0.0, 0.0, 0.0 };


VEC v_set(double x, double y, double z){
	
	VEC c;
	c.y=y;
	c.z=z; 
	return( c ); 
}

VEC v_setpol(double r,double th, double ph){ 
                                             
	VEC c = v_smul(r,v_set(sin(th)*cos(ph),sin(th)*sin(ph),cos(th))); 
	return(c); 
}

void v_out(FILE* file, char* txt, VEC a){ 
     
	char s[100];
	                                      
	sprintf( s,"%s %10.5lf %10.5lf %10.5lf\n",txt,a.x,a.y,a.z );
  	if( file == NULL){
		printf("%s", s );
		getchar();
	}
  	else
	  fprintf( file, "%s", s);
	
}

VEC v_add(VEC a,VEC b){ 
    a.x+=b.x;                   
	a.y+=b.y;
	a.z+=b.z;
	return( a ); 
}

VEC v_sub(VEC a, VEC b){ 
    a.x-=b.x;                     
	a.y-=b.y; 
	a.z-=b.z; 
	return( a ); 
}

VEC v_minus(VEC a){ 
    a.x= -a.x;                
	a.y= -a.y;
	a.z= -a.z;
	return( a ); 
}

VEC v_smul(double s, VEC a){ 
    a.x*=s;                         
	a.y*=s;
	a.z*=s;
	return( a );
}

double v_sprd(VEC a, VEC b){ 
    return( a.x*b.x+a.y*b.y+a.z*b.z );                        
}

double v_sqr(VEC a){ 

    return( a.x*a.x+a.y*a.y+a.z*a.z );                 
}

VEC v_kreuz(VEC a,VEC b){ 

    VEC c;                      
	c.x=a.y*b.z-a.z*b.y; 
	c.y=a.z*b.x-a.x*b.z; 
	c.z=a.x*b.y-a.y*b.x;
	return( c ); 
}

double v_norm(VEC a){
	
	return( sqrt(v_sqr(a)) );                    
}

double v_vdet(VEC a1, VEC a2, VEC a3){ 

    /*Berechnet det A=det(a1,a2,a3) nach Sarrus*/                            
	return( a1.x*a2.y*a3.z+a1.y*a2.z*a3.x+a1.z*a2.x*a3.y
          -a1.z*a2.y*a3.x-a1.x*a2.z*a3.y-a1.y*a2.x*a3.z ); 
}

double v_det(MATRIX A ){ 
                         
    return( A.x.x*A.y.y*A.z.z + A.x.y*A.y.z*A.z.x + A.x.z*A.y.x*A.z.y
          - A.x.z*A.y.y*A.z.x - A.x.x*A.y.z*A.z.y - A.x.y*A.y.x*A.z.z );
}

VEC v_scale(VEC scale, VEC a){
	a.x*=scale.x;                          
	a.y*=scale.y; 
	a.z*=scale.z;  
	return( a ); 	
}

VEC v_scadd(VEC scale, VEC a, VEC pos){ 
    a.x*=scale.x+pos.x;                                   
	a.y*=scale.y+pos.y;
	a.z*=scale.z+pos.z;
	return( a );
}

VEC v_mmul(MATRIX A, VEC x){
	
	VEC c;                           
	c.x=A.x.x*x.x+A.y.x*x.y+A.z.x*x.z; 
	c.y=A.x.y*x.x+A.y.y*x.y+A.z.y*x.z;
	c.z=A.x.z*x.x+A.y.z*x.y+A.z.z*x.z;
	return( c );
}

VEC v_lgs(MATRIX A, VEC b){ 
    
	VEC c;                        
	c.x=v_vdet(b,A.y,A.z);
	c.y=v_vdet(A.x,b,A.z);
	c.z=v_vdet(A.x,A.y,b);
  	return( v_smul( 1/v_det(A),c));
}

int v_equ(VEC a, VEC b){ 
    
	VEC c;                     
	c=v_sub( a,b );
  	return( v_sprd(c,c)<1.0E-5 );
}

VEC v_normal(VEC a){ 
	return( v_smul( 1/v_norm(a), a ) );               
}

VEC v_perpend(VEC a, VEC v){ 
	return( v_sub( a, v_parallel(a,v) ) );                           
}

VEC v_parallel(VEC a,VEC v){ 
	v = v_normal( v ); 
	return( v_smul( v_sprd(a,v), v ) ); 
}

VEC v_zrotate( VEC z, VEC v ){

 double	r,th,ph;

 v_cart2pol( z, &r, &th, &ph );
 return v_polrot( th, ph, v );
}

VEC v_polrot( double th, double ph, VEC v ){
	
 VEC    c;
 double r, cth, sth, cph, sph;

 if( th<0.0 ){      
 	fprintf( stderr, "vec.v_thphrot >  ERROR: polar angle (%lg) negative!!!\n\n", th );
    exit(-1); 
}

 while( ph<0.0 )        ph+=2*M_PI;
 cth = cos(th); sth = sin(th);
 cph = cos(ph); sph = sin(ph);
 c.x = cth*cph*v.x - cth*sph*v.y - sth*v.z;
 c.y = sph    *v.x + cph    *v.y;
 c.z = sth*cph*v.x - sth*sph*v.y + cth*v.z;

 return c;
}


VEC v_polrotinv( double th, double ph, VEC v )
{
 VEC    c;
 double r, cth, sth, cph, sph;

 if( th<0.0 ){
    fprintf( stderr, "vec.v_thphrot >  ERROR: polar angle (%lg) negative!!!\n\n", th );
    exit(-1);       
}
 while( ph<0.0 )   ph+=2*M_PI;
 th  = 2*M_PI-th;  ph  = 2*M_PI-ph;
 cth = cos(th);    sth = sin(th);
 cph = cos(ph);    sph = sin(ph);
 c.x = cph*cth*v.x - sph*v.y - cph*sth*v.z;
 c.y = sph*cth*v.x + cph*v.y - sph*sth*v.z;
 c.z = sth    *v.x           + cth    *v.z;

 return c;
}


double v_transverse(VEC a ){ 

	return( sqrt(a.x*a.x+a.y*a.y) );
}

void v_cart2pol(VEC a, double* r, double* t, double* p){
	
  double	rt2 = a.x*a.x + a.y*a.y;
  if( (*r=sqrt(rt2+a.z*a.z)) < V_ZERO )
  {	*r = *t = *p = 0.;	return;	}
  *t = acos(a.z/(*r));
  if( rt2==0 ) {	*p = 0.;  return;	}
  if( a.y>0. )	*p = acos(a.x/sqrt(rt2));
  else		*p = 2.*M_PI - acos(a.x/sqrt(rt2));
}
