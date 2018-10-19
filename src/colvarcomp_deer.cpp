// -*- c++ -*-

#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvarcomp_deer.h"


/* fresnel_c(x) - Fresnel Cosine Integral
 *  * C(x)=fresnel_c(x)=\dint\limits_{0}^{x}\cos (\frac{\pi}{2}t^{2})dt
 *   */
double fresnel_c(double x);
/* fresnel_s(x) - Fresnel Sine Integral
 *  * S(x)=fresnel_s(x)=\dint\limits_{0}^{x}\sin (\frac{\pi}{2}t^{2})dt
 *   */
double fresnel_s(double x);
/* Additional functions*/
/* fresnel_c1(x)
 *  * fresnel_c1(x)=fresnel_c(x*sqrt(2/pi))=
 *   * = \sqrt{\frac{2}{\pi }}\dint\limits_{0}^{x}\cos (t^{2})dt
 *    */
double fresnel_c2(double x);
/* fresnel_s1(x)
 *  * fresnel_s1(x)=fresnel_s(x*sqrt(2/pi))=
 *   * = \sqrt{\frac{2}{\pi }}\dint\limits_{0}^{x}\sin (t^{2})dt
 *    */
double fresnel_s2(double x);



// deer keernel routine

double colvar::deer_kernel::kdeer(cvm::real const &r, cvm::real const &t)
{
  double rc=r; //units of Ångström
  double tc=t; //units of nanoseconds
  double const g=-2.00231930436182;
  double const ub=9.27400968;
  double const pi=4*atan(1.);
  double const u0=4*pi;
  double const ht=1.054571800;
  double omed=(g*g)*(ub*ub)*u0/(4*pi*ht*(rc*rc*rc));
  double k=sqrt(6*omed*tc/pi);
  double c=fresnel_c(k);
  double s=fresnel_s(k);
  double gdeer=0.;
  if(tc<0) tc=-tc;
  if(rc>0 && tc!=0){
      gdeer=(sqrt(pi*((s*s)+(c*c))/(6*omed*tc)))*(cos(omed*tc-(atan2(s,c))));
  }
  if(rc==0) gdeer=0.;
  if(tc==0) gdeer=0.;
  return gdeer;
}

// end deer keernel routine

colvar::deer_kernel::deer_kernel(std::string const &conf)
  : cvc(conf)
{
  function_type = "deer_kernel";
  get_keyval(conf, "deertimefile", deer_time_file);
  int nlines=0;
  int i;
  int t;
  // read deertimefile
  std::string line;
  std::ifstream deerfile (deer_time_file);
  if (deerfile.is_open()){
    while ( getline (deerfile,line) )
    {
       if (line.size() == 0)
         continue;
       nlines++;
    }
    deersize=nlines;
    // allocations
    timesdeer.resize(deersize);
    deerexpvalues.resize(deersize);
    // read file again
    deerfile.clear();
    deerfile.seekg(0);
    // assign times and experimental values
    for (t=0; t<deersize; t++){
        deerfile >> timesdeer[t];
        deerfile >> deerexpvalues[t];
        getline (deerfile,line);
    }
    deerfile.close();
  }
  // define deer kernel grid
  get_keyval(conf, "deerWidth", deerwidth, 0.1);
  get_keyval(conf, "deerLower", deerlower, 0);
  get_keyval(conf, "deerUpper", deerupper, 100);
  rpoints=std::floor( (deerupper-deerlower) / deerwidth );

  deerk.resize(rpoints,deersize);

  // assign deer kernel grid

  for (i=0; i<rpoints; i++){
     cvm::real const rval = deerlower+(i+0.5)*deerwidth;
     for (t=0; t<deersize; t++){
        deerk[i][t]=kdeer(rval,timesdeer[t])-deerexpvalues[t];
     }
  }

  if (get_keyval(conf, "forceNoPBC", b_no_PBC, false)) {
    cvm::log("Computing distance using absolute positions (not minimal-image)");
  }

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  x.type(colvarvalue::type_vector);
  x.vector1d_value.resize(deersize);
}

colvar::deer_kernel::deer_kernel()
{
  function_type = "deer_kernel";
  x.type(colvarvalue::type_vector);
}

void colvar::deer_kernel::calc_value()
{
  int t;
  x.vector1d_value.resize(deersize);
  deerder.resize(deersize);

  if (b_no_PBC) {
    dist_v = group2->center_of_mass() - group1->center_of_mass();
  } else {
    dist_v = cvm::position_distance(group1->center_of_mass(),
                                     group2->center_of_mass());
  }

  int deerbin=floor( (dist_v.norm()-deerlower) / deerwidth );

  if(deerbin<0 || deerbin >=rpoints ){
    cvm::error("distance out of boundaries, expand deerLower or deerUpper!");
    return;
  }

  for (t=0; t<deersize; t++){
     x.vector1d_value[t]=deerk[deerbin][t]; // value calculated
     // calculate gradients
     cvm::rvector const u = dist_v.unit();
     if(deerbin==0) deerder[t] = (deerk[deerbin+1][t]-deerk[deerbin][t])/(deerwidth);
     if(deerbin==rpoints-1) deerder[t] = (deerk[deerbin][t]-deerk[deerbin-1][t])/(deerwidth);
     if(deerbin>0 && deerbin<rpoints-1) {
       deerder[t] = (deerk[deerbin+1][t]-deerk[deerbin-1][t])/(2.*deerwidth);
     }
     group1->set_weighted_gradient(-1.0 * u * deerder[t]);
     group2->set_weighted_gradient(       u * deerder[t]);
  }
}

void colvar::deer_kernel::calc_gradients()
{
  // will be calculated on the fly in apply_force()
}

void colvar::deer_kernel::apply_force(colvarvalue const &force)
{
  cvm::rvector const u = dist_v.unit();
  int t;
  for (t=0; t<deersize; t++){
     if (!group1->noforce)
       group1->apply_force(-1.0 * u * force[t] * deerder[t]);

     if (!group2->noforce)
       group2->apply_force(       u * force[t] * deerder[t]);
  }
}

cvm::real colvar::deer_kernel::dist2(colvarvalue const &x1,
                                        colvarvalue const &x2) const
{
  return (x1.vector1d_value - x2.vector1d_value).norm2();
}


colvarvalue colvar::deer_kernel::dist2_lgrad(colvarvalue const &x1,
                                                colvarvalue const &x2) const
{
  return colvarvalue((x1.vector1d_value - x2.vector1d_value), colvarvalue::type_vector);
}


colvarvalue colvar::deer_kernel::dist2_rgrad(colvarvalue const &x1,
                                                colvarvalue const &x2) const
{
  return colvarvalue((x2.vector1d_value - x1.vector1d_value), colvarvalue::type_vector);
}



// fresnel integrals routines (by expansion to Chebyshev series)

static double fresnel_cos_0_8(double x)
{
 double x_8 = x/8.0;
 double xx = 2.0*x_8*x_8 - 1.0;

 double t0 = 1.0;
 double t1 = xx;
 double sumC = f_data_fres_a[0] + f_data_fres_a[1]*t1;
 double t2;
 int n;
 for (n=2; n < 18; n++)
 {
  t2 = 2.0*xx*t1 - t0;
  sumC += f_data_fres_a[n]*t2;
  t0 = t1; t1 = t2;
 }
 return _1_sqrt_2pi_fres*sqrt(x)*sumC;
}

static double fresnel_sin_0_8(double x)
{
 double x_8 = x/8.0;
 double xx = 2.0*x_8*x_8 - 1.0;
 double t0 = 1.;
 double t1 = xx;
 double ot1 = x_8;
 double ot2 = 2.0*x_8*t1 - ot1;
 double sumS = f_data_fres_b[0]*ot1 + f_data_fres_b[1]*ot2;
 int n;
 double t2;
 for (n=2; n < 17; n++)
 {
  t2 = 2.0*xx*t1 - t0;
  ot1 = ot2;
  ot2 = 2.0*x_8*t2 - ot1;
  sumS += f_data_fres_b[n]*ot2;
  t0 = t1; t1 = t2;
 }
 return _1_sqrt_2pi_fres*sqrt(x)*sumS;
}

static double fresnel_cos_8_inf(double x)
{
 double xx = 128.0/(x*x) - 1.0;   /* 2.0*(8/x)^2 - 1 */
 double t0 = 1.0;
 double t1 = xx;
 double sumP = f_data_fres_e[0] + f_data_fres_e[1]*t1;
 double sumQ = f_data_fres_f[0] + f_data_fres_f[1]*t1;
 double t2;
 int n;
 for(n = 2; n < 35; n++)
 {
   t2 = 2.0*xx*t1 - t0;
   sumP += f_data_fres_e[n]*t2; /*  sumP += f_data_fres_e[n]*ChebyshevT(n,xx) */
   sumQ += f_data_fres_f[n]*t2; /*  sumQ += f_data_fres_f[n]*ChebyshevT(n,xx) */
   t0 = t1; t1 = t2;
 }
 for(n = 35; n < 41; n++)
 {
   t2 = 2.0*xx*t1 - t0;
   sumP += f_data_fres_e[n]*t2; /*  sumP += f_data_fres_e[n]*ChebyshevT(n,xx) */
   t0 = t1; t1 = t2;
 }
 return 0.5 - _1_sqrt_2pi_fres*(0.5*sumP*cos(x)/x - sumQ*sin(x))/sqrt(x);
}

static double fresnel_sin_8_inf(double x)
{
 double xx = 128.0/(x*x) - 1.0;   /* 2.0*(8/x)^2 - 1 */
 double t0 = 1.0;
 double t1 = xx;
 double sumP = f_data_fres_e[0] + f_data_fres_e[1]*t1;
 double sumQ = f_data_fres_f[0] + f_data_fres_f[1]*t1;
 double t2;
 int n;
 for(n = 2; n < 35; n++)
 {
   t2 = 2.0*xx*t1 - t0;
   sumP += f_data_fres_e[n]*t2; /*  sumP += f_data_fres_e[n]*ChebyshevT(n,xx) */
   sumQ += f_data_fres_f[n]*t2; /*  sumQ += f_data_fres_f[n]*ChebyshevT(n,xx) */
   t0 = t1; t1 = t2;
 }
 for(n = 35; n < 41; n++)
 {
   t2 = 2.0*xx*t1 - t0;
   sumP += f_data_fres_e[n]*t2; /*  sumQ += f_data_fres_f[n]*ChebyshevT(n,xx) */
   t0 = t1; t1 = t2;
 }
 return 0.5 - _1_sqrt_2pi_fres*(0.5*sumP*sin(x)/x + sumQ*cos(x))/sqrt(x);
}


double fresnel_c(double x)
{
  double xx = x*x*pi_2_fres;
  double ret_val;
  if(xx<=8.0)
   ret_val = fresnel_cos_0_8(xx);
  else
   ret_val = fresnel_cos_8_inf(xx);
  return (x<0.0) ? -ret_val : ret_val;
}

double fresnel_s(double x)
{
  double xx = x*x*pi_2_fres;
  double ret_val;
  if(xx<=8.0)
   ret_val = fresnel_sin_0_8(xx);
  else
   ret_val = fresnel_sin_8_inf(xx);
  return (x<0.0) ? -ret_val : ret_val;
}

double fresnel_c1(double x)
{
  return fresnel_c(x*sqrt_2_pi_fres);
}

double fresnel_s1(double x)
{
  return fresnel_s(x*sqrt_2_pi_fres);
}

// end fresnel integral routines
