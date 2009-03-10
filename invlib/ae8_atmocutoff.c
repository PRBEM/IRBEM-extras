#include <gsl/gsl_math.h>

double AE8_AtmoCutoff(const double Lm) {
  /* returns atmospheric cutoff in equatorial pitch angle, degrees */
  /* From Vette, AE-8, equation 5.4 */
  /* returns NaN for Lm<1 */
  double BB0, ac;
  if (Lm<1) {
    return(GSL_NAN);
  } else if (Lm < 2.4) {
    BB0 = 0.6572*pow(Lm,3.452);
  } else if (Lm <= 3.0) {
    BB0 = 0.196*pow(Lm,4.878);
  } else {
    BB0 = 1.4567*pow(Lm,3.050);
  }
  BB0 = GSL_MAX(BB0,1.0);
  ac = 180.0/M_PI*asin(1.0/sqrt(BB0));
  return(ac);
}

