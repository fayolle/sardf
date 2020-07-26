#ifndef SARDF_H
#define SARDF_H

#ifdef __cplusplus
extern "C" {
#endif

double uni_sardf_1(const double xt, const double yt, const double R, const double alpha);
double uni_sardf_2(double xt, double yt, double R1, double R2);
double uni_sardf_and_gradient_1(const double xt, const double yt, const double R, const double alpha, const double gxt[3], const double gyt[3], double grad[3]);
double uni_sardf_and_gradient_2(const double xt, const double yt, const double R1, const double R2, const double gxt[3], const double gyt[3], double grad[3]);
double int_sardf_1 (double xt, double yt, double R, double alpha);
double int_sardf_2 (double xt, double yt, double R1, double R2);
double int_sardf_and_gradient_1 (const double xt, const double yt, const double R, const double alpha, const double gxt[3], const double gyt[3], double grad[3]);
double int_sardf_and_gradient_2 (const double xt, const double yt, const double R1, const double R2, const double gxt[3], const double gyt[3], double grad[3]);

#ifdef __cplusplus
}
#endif
  
#endif
