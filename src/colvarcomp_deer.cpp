// -*- c++ -*-

#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvarcomp_deer.h"


namespace DEER_Kernel {

/// \brief fresnel integrals (by expansion to Chebyshev series)

static const double sqrt_pi_2_fres   = 1.2533141373155002512078826424; /* sqrt(pi/2) */
static const double sqrt_2_pi_fres   = 0.7978845608028653558798921199; /* sqrt(2/pi) */
static const double _1_sqrt_2pi_fres = 0.3989422804014326779399460599; /* 1/sqrt(2*pi) */
static const double pi_2_fres        = 1.5707963267948966192313216916; /* pi/2 */

static double f_data_fres_a[18] =
                {
                  0.76435138664186000189,
                 -0.43135547547660179313,
                  0.43288199979726653054,
                 -0.26973310338387111029,
                  0.08416045320876935378,
                 -0.01546524484461381958,
                  0.00187855423439822018,
                 -0.00016264977618887547,
                  0.00001057397656383260,
                 -0.00000053609339889243,
                  0.00000002181658454933,
                 -0.00000000072901621186,
                  0.00000000002037332546,
                 -0.00000000000048344033,
                  0.00000000000000986533,
                 -0.00000000000000017502,
                  0.00000000000000000272,
                 -0.00000000000000000004
                };

static double f_data_fres_b[17] =
                {
                  0.63041404314570539241,
                 -0.42344511405705333544,
                  0.37617172643343656625,
                 -0.16249489154509567415,
                  0.03822255778633008694,
                 -0.00564563477132190899,
                  0.00057454951976897367,
                 -0.00004287071532102004,
                  0.00000245120749923299,
                 -0.00000011098841840868,
                  0.00000000408249731696,
                 -0.00000000012449830219,
                  0.00000000000320048425,
                 -0.00000000000007032416,
                  0.00000000000000133638,
                 -0.00000000000000002219,
                  0.00000000000000000032
                };

static double f_data_fres_e[41] =
                {
                    0.97462779093296822410,
                   -0.02424701873969321371,
                    0.00103400906842977317,
                   -0.00008052450246908016,
                    0.00000905962481966582,
                   -0.00000131016996757743,
                    0.00000022770820391497,
                   -0.00000004558623552026,
                    0.00000001021567537083,
                   -0.00000000251114508133,
                    0.00000000066704761275,
                   -0.00000000018931512852,
                    0.00000000005689898935,
                   -0.00000000001798219359,
                    0.00000000000594162963,
                   -0.00000000000204285065,
                    0.00000000000072797580,
                   -0.00000000000026797428,
                    0.00000000000010160694,
                   -0.00000000000003958559,
                    0.00000000000001581262,
                   -0.00000000000000646411,
                    0.00000000000000269981,
                   -0.00000000000000115038,
                    0.00000000000000049942,
                   -0.00000000000000022064,
                    0.00000000000000009910,
                   -0.00000000000000004520,
                    0.00000000000000002092,
                   -0.00000000000000000982,
                    0.00000000000000000467,
                   -0.00000000000000000225,
                    0.00000000000000000110,
                   -0.00000000000000000054,
                    0.00000000000000000027,
                   -0.00000000000000000014,
                    0.00000000000000000007,
                   -0.00000000000000000004,
                    0.00000000000000000002,
                   -0.00000000000000000001,
                    0.00000000000000000001
        };

static double f_data_fres_f[35] =
                {
                    0.99461545179407928910,
                   -0.00524276766084297210,
                    0.00013325864229883909,
                   -0.00000770856452642713,
                    0.00000070848077032045,
                   -0.00000008812517411602,
                    0.00000001359784717148,
                   -0.00000000246858295747,
                    0.00000000050925789921,
                   -0.00000000011653400634,
                    0.00000000002906578309,
                   -0.00000000000779847361,
                    0.00000000000222802542,
                   -0.00000000000067239338,
                    0.00000000000021296411,
                   -0.00000000000007041482,
                    0.00000000000002419805,
                   -0.00000000000000861080,
                    0.00000000000000316287,
                   -0.00000000000000119596,
                    0.00000000000000046444,
                   -0.00000000000000018485,
                    0.00000000000000007527,
                   -0.00000000000000003131,
                    0.00000000000000001328,
                   -0.00000000000000000574,
                    0.00000000000000000252,
                   -0.00000000000000000113,
                    0.00000000000000000051,
                   -0.00000000000000000024,
                    0.00000000000000000011,
                   -0.00000000000000000005,
                    0.00000000000000000002,
                   -0.00000000000000000001,
                    0.00000000000000000001
                };

/* fresnel_c(x) - Fresnel Cosine Integral
 *  * C(x)=fresnel_c(x)=\dint\limits_{0}^{x}\cos (\frac{\pi}{2}t^{2})dt
 *   */
static double fresnel_c(double x);
/* fresnel_s(x) - Fresnel Sine Integral
 *  * S(x)=fresnel_s(x)=\dint\limits_{0}^{x}\sin (\frac{\pi}{2}t^{2})dt
 *   */
static double fresnel_s(double x);
/* Additional functions*/
/* fresnel_c1(x)
 *  * fresnel_c1(x)=fresnel_c(x*sqrt(2/pi))=
 *   * = \sqrt{\frac{2}{\pi }}\dint\limits_{0}^{x}\cos (t^{2})dt
 *    */
static double fresnel_c2(double x);
/* fresnel_s1(x)
 *  * fresnel_s1(x)=fresnel_s(x*sqrt(2/pi))=
 *   * = \sqrt{\frac{2}{\pi }}\dint\limits_{0}^{x}\sin (t^{2})dt
 *    */
static double fresnel_s2(double x);

}


// DEER kernel routine

double colvar::deer_kernel::kdeer(cvm::real const &r, cvm::real const &t)
{ 
  double rc=r; //units of Ångström
  double tc=t; //units of nanoseconds
  tc=sqrt(tc*tc); // deer_kernel is a symmetric function of tc
  double const g=-2.00231930436182;
  double const ub=9.27400968;
  double const pi=4*atan(1.0);
  double const u0=4*pi;
  double const ht=1.054571800;
  double omed=(g*g)*(ub*ub)*u0/(4*pi*ht*(rc*rc*rc));
  double k=sqrt(6*omed*tc/pi);
  double c=fresnel_c(k);
  double s=fresnel_s(k);
  double gdeer=0.;
  if(rc>0 && tc!=0){
      gdeer=(sqrt(pi*((s*s)+(c*c))/(6*omed*tc)))*(cos(omed*tc-(atan2(s,c))));
  }
  if(rc==0) gdeer=0.0; //checked
  if(tc==0) gdeer=1.0; //ckecked
  return gdeer;
}

// end DEER keernel routine

// DEER keernel analytical derivatives routine

double colvar::deer_kernel::kdeer_der(cvm::real const &r, cvm::real const &t)
{
  double rc=r; //units of Ångström
  double tc=t; //units of nanoseconds
  tc=sqrt(tc*tc); // deer_kernel derivative is a symmetric function of tc
  double const g=-2.00231930436182;
  double const ub=9.27400968;
  double const pi=4*atan(1.0);
  double const u0=4*pi;
  double const ht=1.054571800;
  double omed=(g*g)*(ub*ub)*u0/(4*pi*ht*(rc*rc*rc));
  double omedp=-3*(g*g)*(ub*ub)*u0/(4*pi*ht*(rc*rc*rc*rc));
  double k=sqrt(6*omed*tc/pi);
  double c=fresnel_c(k);
  double s=fresnel_s(k);
  double sp=0.5*sin(pi*k*k/2)*(6*omedp*tc/pi)/k;
  double cp=0.5*cos(pi*k*k/2)*(6*omedp*tc/pi)/k;
  double gdeerp=0.;
  if(rc>0 && tc!=0){
    gdeerp=(0.5*(pi*(2*s*sp+2*c*cp)*6*omed*tc-pi*6*omedp*tc*((s*s)+(c*c)))*
     (cos(omed*tc-(atan2(s,c))))/(((6*omed*tc)*(6*omed*tc))*sqrt(pi*((s*s)+
     (c*c))/(6*omed*tc))))-(sqrt(pi*((s*s)+(c*c))/(6*omed*tc)))*sin(omed*tc-
     (atan2(s,c)))*(omedp*tc-(1/(1+((s/c)*(s/c))))*((sp*c-s*cp)/(c*c)));
  }
  if(rc==0) gdeerp=0.0; //checked
  if(tc==0) gdeerp=0.0; //ckecked
  return gdeerp;
}

// end DEER keernel analytical derivatives routine


colvar::deer_kernel::deer_kernel(std::string const &conf)
  : cvc(conf)
{
  function_type = "deer_kernel";
  x.type(colvarvalue::type_vector);

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
       deerk[i][t] = static_cast<cvm::real>(kdeer(rval,timesdeer[t])) -
         deerexpvalues[t];
     }
  }

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

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

  if (is_enabled(f_cvc_pbc_minimum_image)) {
    dist_v = cvm::position_distance(group1->center_of_mass(),
                                    group2->center_of_mass());
  } else {
    dist_v = group2->center_of_mass() - group1->center_of_mass();
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


namespace DEER_Kernel {

/// \brief Fresnel integrals (by expansion to Chebyshev series)

static const double sqrt_pi_2_fres   = 1.2533141373155002512078826424; /* sqrt(pi/2) */
static const double sqrt_2_pi_fres   = 0.7978845608028653558798921199; /* sqrt(2/pi) */
static const double _1_sqrt_2pi_fres = 0.3989422804014326779399460599; /* 1/sqrt(2*pi) */
static const double pi_2_fres        = 1.5707963267948966192313216916; /* pi/2 */

static const double f_data_fres_a[18] =
                {
                  0.76435138664186000189,
                 -0.43135547547660179313,
                  0.43288199979726653054,
                 -0.26973310338387111029,
                  0.08416045320876935378,
                 -0.01546524484461381958,
                  0.00187855423439822018,
                 -0.00016264977618887547,
                  0.00001057397656383260,
                 -0.00000053609339889243,
                  0.00000002181658454933,
                 -0.00000000072901621186,
                  0.00000000002037332546,
                 -0.00000000000048344033,
                  0.00000000000000986533,
                 -0.00000000000000017502,
                  0.00000000000000000272,
                 -0.00000000000000000004
                };

static const double f_data_fres_b[17] =
                {
                  0.63041404314570539241,
                 -0.42344511405705333544,
                  0.37617172643343656625,
                 -0.16249489154509567415,
                  0.03822255778633008694,
                 -0.00564563477132190899,
                  0.00057454951976897367,
                 -0.00004287071532102004,
                  0.00000245120749923299,
                 -0.00000011098841840868,
                  0.00000000408249731696,
                 -0.00000000012449830219,
                  0.00000000000320048425,
                 -0.00000000000007032416,
                  0.00000000000000133638,
                 -0.00000000000000002219,
                  0.00000000000000000032
                };

static const double f_data_fres_e[41] =
                {
                    0.97462779093296822410,
                   -0.02424701873969321371,
                    0.00103400906842977317,
                   -0.00008052450246908016,
                    0.00000905962481966582,
                   -0.00000131016996757743,
                    0.00000022770820391497,
                   -0.00000004558623552026,
                    0.00000001021567537083,
                   -0.00000000251114508133,
                    0.00000000066704761275,
                   -0.00000000018931512852,
                    0.00000000005689898935,
                   -0.00000000001798219359,
                    0.00000000000594162963,
                   -0.00000000000204285065,
                    0.00000000000072797580,
                   -0.00000000000026797428,
                    0.00000000000010160694,
                   -0.00000000000003958559,
                    0.00000000000001581262,
                   -0.00000000000000646411,
                    0.00000000000000269981,
                   -0.00000000000000115038,
                    0.00000000000000049942,
                   -0.00000000000000022064,
                    0.00000000000000009910,
                   -0.00000000000000004520,
                    0.00000000000000002092,
                   -0.00000000000000000982,
                    0.00000000000000000467,
                   -0.00000000000000000225,
                    0.00000000000000000110,
                   -0.00000000000000000054,
                    0.00000000000000000027,
                   -0.00000000000000000014,
                    0.00000000000000000007,
                   -0.00000000000000000004,
                    0.00000000000000000002,
                   -0.00000000000000000001,
                    0.00000000000000000001
        };

static const double f_data_fres_f[35] =
                {
                    0.99461545179407928910,
                   -0.00524276766084297210,
                    0.00013325864229883909,
                   -0.00000770856452642713,
                    0.00000070848077032045,
                   -0.00000008812517411602,
                    0.00000001359784717148,
                   -0.00000000246858295747,
                    0.00000000050925789921,
                   -0.00000000011653400634,
                    0.00000000002906578309,
                   -0.00000000000779847361,
                    0.00000000000222802542,
                   -0.00000000000067239338,
                    0.00000000000021296411,
                   -0.00000000000007041482,
                    0.00000000000002419805,
                   -0.00000000000000861080,
                    0.00000000000000316287,
                   -0.00000000000000119596,
                    0.00000000000000046444,
                   -0.00000000000000018485,
                    0.00000000000000007527,
                   -0.00000000000000003131,
                    0.00000000000000001328,
                   -0.00000000000000000574,
                    0.00000000000000000252,
                   -0.00000000000000000113,
                    0.00000000000000000051,
                   -0.00000000000000000024,
                    0.00000000000000000011,
                   -0.00000000000000000005,
                    0.00000000000000000002,
                   -0.00000000000000000001,
                    0.00000000000000000001
                };


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

// end Fresnel integral routines
}
