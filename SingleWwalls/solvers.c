#define MAX(X1,X2) ( (X1>X2) ? X1 : X2 )

/******************************************************
**               QUADRATICSOLVER
** To always have accurate roots, even when a or c are
     small enough to incur in cancellation errors:
       q = -0.5 * (b + sign(b)*sqrt(b*b-4*a*c)) ;
       x1 = q / a ;
       x2 = c / q ;
******************************************************/


isfinite(blablabla)
/******************************************************
**               CUBICSOLVER
** Returns the real solutions of the equation 
**   a * x*x*x + b * x*x + c * x + d = 0
** It returns an array sol[] wherein sol[0] is the
**    number of real roots, followed by the roots 
**    ordered from the one with the largest absolute value
******************************************************/
void efficientcubicsolver(double a, double b, double c, double d, double *sol)
{
  if (a==0) return NAN ;
  b = b / a ;
  c = c / a ;
  d = d / a ;
  double Q = b*b/9. - c/3. ;
  double R = ( 2.*b*b*b - 9.*b*c + 27.*d ) / 54. ;
  // associated depressed cubic : y^3 + 3Qy + 2R = 0

  if (R == 0) {
    if (Q > 0) {
      sol[0] = 1 ;
      sol[1] = - b / 3 ;
    } else {
      sol[0] = 3 ;
      sol[2] = - b / 3 ;
      if (b < 0) {
	sol[1] = - b / 3 + sqrt(-3*Q) ;
	sol[3] = - b / 3 - sqrt(-3*Q) ;
      } else {
	sol[1] = - b / 3 - sqrt(-3*Q) ;
	sol[3] = - b / 3 + sqrt(-3*Q) ;
      }
    }
  } else {
    if (fabs(Q) < fabs(R)) {
      double ratio = Q / R ;
      double K = 1 - Q * ratio ;
      if (K < 0) {
	double theta = acos(1. / ratio / sqrt(Q)) ;
	sol[0] = 3 ;
	sol[3] = -2. * sqrt(Q) * cos((theta - 2*M_PI) / 3) - b / 3 ;
	if (theta < M_PI/3) {
	  sol[1] = -2. * sqrt(Q) * cos(theta / 3) - b / 3 ;
	  sol[2] = -2. * sqrt(Q) * cos((theta + 2*M_PI) / 3) - b / 3 ;
	} else {
	  sol[2] = -2. * sqrt(Q) * cos(theta / 3) - b / 3 ;
	  sol[1] = -2. * sqrt(Q) * cos((theta + 2*M_PI) / 3) - b / 3 ;
	}
      } else {
	double A = -sign(R) * cbrt(fabs(R) * (1 + sqrt(K))) ;
	double B = 0 ;
	if (A != 0) B = Q / A ;
	sol[0] = 1 ;
	sol[3] = A + B - b / 3 ;
      }

    } else {
      double ratio = R / Q ;
      double K = sign(Q)*(ratio * ratio / Q - 1) ;
      if (K < 0) {
	double theta = acos(ratio / sqrt(Q)) ;
	sol[0] = 3 ;
	sol[3] = -2. * sqrt(Q) * cos((theta - 2*M_PI) / 3) - b / 3 ;
	if (theta < M_PI/3) {
	  sol[1] = -2. * sqrt(Q) * cos(theta / 3) - b / 3 ;
	  sol[2] = -2. * sqrt(Q) * cos((theta + 2*M_PI) / 3) - b / 3 ;
	} else {
	  sol[2] = -2. * sqrt(Q) * cos(theta / 3) - b / 3 ;
	  sol[1] = -2. * sqrt(Q) * cos((theta + 2*M_PI) / 3) - b / 3 ;
	}
      } else {
	double A = -sign(R) * cbrt(fabs(R) + sqrt(Q) * fabs(Q) * sqrt(K)) ;
	double B = 0 ;
	if (A != 0) B = Q / A ;
	sol[0] = 1 ;
	sol[3] = A + B - b / 3 ;
      }
    }
  }

  // refinement of the solution with the Newton-Raphson method
  for (i=1; i<=sol[0]; i++) {
    float x = sol[i] , x0 ;
    float f = x*x*x + b*x*x + c*x + d , df , f0 ;
    int iter = 0 , iterMAX = 10 ;
    while( fabs(f) > 2.22045E−16*MAX( fabs(x*x*x) , MAX( fabs(b*x*x) , MAX( fabs(c*x) , d ) ) ) && iter < iterMAX ) {
      df = 3.*x*x + 2.*b*x + c ;
      if (df == 0) iter = 20 ; //terminate
      else {
	x0 = x ;
	f0 = f ;
	x = x - f / df ;
	f = x*x*x + b*x*x + c*x + d ;
	if (f == 0 ) iter = 20 ; //terminate
	if (fabs(f) > fabs(f0)) {
	  x = x0 ;
	  iter = 20 ; //terminate
	} else iter ++ ;
      }
    }
    sol[i] = x ;
  }

  return sol ;
}



/******************************************************
**               QUARTICSOLVER
** Returns the real solutions of the equation 
**   e * x*x*x*x + a * x*x*x + b * x*x + c * x + d = 0
** It returns an array sol[] wherein sol[0] is the
**    number of real roots, followed by the roots
******************************************************/
void efficientquarticsolver(double e, double a, double b, double c, double d, double *sol)
{
  double sol[5] ;
  sol[0] = 0 ;

  // degenerate case
  if (e==0) {
    efficientcubicsolver(a, b, c, d, sol) ;
    return ;
  }

  double a = a / e ;
  double b = b / e ;
  double c = c / e ;
  double d = d / e ;

  // biquadratic case
  double disc = (a*a*a - 4.*a*b + 8.*c) / 8 ;
  if (disc == 0 && (a*a-4.*b)/a/a > 0.02) {
    double a1 = 0.5 * (b - 3.*a*a/8) ;
    double a2 = -3.*a*a*a*a/256 + a*a*b/16 - a*c/4 + d ;
    double disc1 = a1*a1 - a2 ;
    if(disc1 >= 0) {
      double squared_sol = -a1 + sign(a1)*sqrt(disc1) ;
      if (squared_sol > 0) {
	sol[0] ++ ;
	sol[ sol[0] ] = -a/4 + sqrt(squared_sol) ;
	sol[0] ++ ;
	sol[ sol[0] ] = -a/4 - sqrt(squared_sol) ;
      }
      squared_sol = a2 / squared_sol ;
      if (squared_sol > 0) {
	sol[0] ++ ;
	sol[ sol[0] ] = -a/4 + sqrt(squared_sol) ;
	sol[0] ++ ;
	sol[ sol[0] ] = -a/4 - sqrt(squared_sol) ;
      }
    }
    return ;
  }

  //write according to the article of Alberto and Cristiano...
}
