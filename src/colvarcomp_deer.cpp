// -*- c++ -*-

#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvarcomp_deer.h"


namespace DEER_Kernel {

/// \brief fresnel integrals (by expansion to Chebyshev series)

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
  double c=DEER_Kernel::fresnel_c(k);
  double s=DEER_Kernel::fresnel_s(k);
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
  double c=DEER_Kernel::fresnel_c(k);
  double s=DEER_Kernel::fresnel_s(k);
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
  : cvc()
{
  function_type = "deer_kernel";
  enable(f_cvc_implicit_gradient);
  x.type(colvarvalue::type_vector);
  init(conf);
}


int colvar::deer_kernel::init(std::string const &conf)
{
  int error_code = COLVARS_OK;

  error_code |= cvc::init(conf);

  if (get_keyval(conf, "deertimefile", deer_time_file)) {
    int nlines=0;
    int i;
    int t;
    // read deertimefile
    std::string line;
    std::ifstream deerfile (deer_time_file);
    if (deerfile.is_open()){
      // count the lines
      while ( cvm::getline(deerfile, line) )
        {
          if (line.size() == 0)
            continue;
          nlines++;
        }
      // allocations
      timesdeer.resize(nlines);
      deerexpvalues.resize(nlines);
      // read file again
      deerfile.clear();
      deerfile.seekg(0);
      // assign times and experimental values
      for (t=0; t< nlines; t++){
        deerfile >> timesdeer[t];
        deerfile >> deerexpvalues[t];
        getline (deerfile,line);
      }
      deerfile.close();
    } else {
      return cvm::error("Unable to open deerTimeFile", INPUT_ERROR);
    }
  }

  if (timesdeer.size() == 0) {
    return cvm::error("Must provide a deerTimeFile with at least one "
                      "non-empty line.\m", INPUT_ERROR);
  }

  get_keyval(conf, "useDEERGrid", deer_grid, false);
  if (!deer_grid) deer_anal_der=true;

  if (deer_grid){
    // define deer kernel grid
    get_keyval(conf, "deerWidth", deerwidth, 0.1);
    get_keyval(conf, "deerLower", deerlower, 0);
    get_keyval(conf, "deerUpper", deerupper, 160);
    get_keyval(conf, "useAnalyDer", deer_anal_der, false);
    rpoints=std::floor( (deerupper-deerlower) / deerwidth );

    deerk.resize(rpoints, times.size());
    if (deer_anal_der) deerk_der.resize(rpoints, timesdeer.size());
    // assign deer kernel grid

    for (i=0; i<rpoints; i++){
       cvm::real const rval = deerlower+(i+0.5)*deerwidth;
       for (t=0; t<timesdeer.size(); t++){
          deerk[i][t]=kdeer(rval,timesdeer[t]);
       }
    }
    if (deer_anal_der){
      for (i=0; i<rpoints; i++){
         cvm::real const rval = deerlower+(i+0.5)*deerwidth;
         for (t=0; t<timesdeer.size(); t++){
            deerk_der[i][t]=kdeer_der(rval,timesdeer[t]);
         }
      }
    }
  }

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  x.vector1d_value.resize(timesdeer.size());

  return error_code;
}


colvar::deer_kernel::deer_kernel()
{
  function_type = "deer_kernel";
  enable(f_cvc_implicit_gradient);
  x.type(colvarvalue::type_vector);
}


template<bool gradients>
int colvar::deer_kernel::compute_deer_kernel(cvm::vector1d<cvm::real> &kernel,
                                             cvm::vector1d<cvm::real> &kernel_deriv,
                                             cvm::vector1d<cvm::real> const &times)
{
  int t;

  size_t const deersize = times.size();

  cvm::real const r = dist_v.norm();

  if (!deer_grid) {

    // Compute k(t) explicitly for each t
    for (t = 0; t < deersize; t++){
      kernel[t] = kdeer(r, times[t]);
      if (gradients) kernel_deriv[t] = kdeer_der(r, times[t]);
    }

  } else {

    int const deerbin = std::floor( (r-deerlower) / deerwidth );

    if (deerbin < 0 || deerbin >= rpoints ){
      return cvm::error("distance out of boundaries for DEER component \""+name+
                        "\", please expand deerLower and/or deerUpper.",
                        INPUT_ERROR);
    }

    if (deer_anal_der) {
      for (t=0; t<deersize; t++){
        kernel[t]=deerk[deerbin][t]; // value calculated
        if (gradients) kernel_deriv[t]=deerk_der[deerbin][t];
      }
    } else {
      for (t=0; t<deersize; t++){
        kernel[t]=deerk[deerbin][t]; // value calculated
        if (gradients) {
          if (deerbin==0) kernel_deriv[t] = (deerk[deerbin+1][t]-deerk[deerbin][t])/(deerwidth);
          if (deerbin==rpoints-1) kernel_deriv[t] = (deerk[deerbin][t]-deerk[deerbin-1][t])/(deerwidth);
          if (deerbin>0 && deerbin<rpoints-1) {
            kernel_deriv[t] = (deerk[deerbin+1][t]-deerk[deerbin-1][t])/(2.*deerwidth);
          }
        }
      }
    }
  }
  return COLVARS_OK;
}


void colvar::deer_kernel::calc_value()
{
  size_t const deersize = timesdeer.size();
  kernel_deriv.resize(deersize);

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    dist_v = group2->center_of_mass() - group1->center_of_mass();
  } else {
    dist_v = cvm::position_distance(group1->center_of_mass(),
                                    group2->center_of_mass());
  }

  if (is_enabled(f_cvc_gradient)) {
    compute_deer_kernel<true>(x.vector1d_value, deer_deriv, timesdeer);
  } else {
    compute_deer_kernel<false>(x.vector1d_value, deer_deriv, timesdeer);
  }
}


void colvar::deer_kernel::calc_gradients()
{
  // calculated on the fly in apply_force() and not stored
}


void colvar::deer_kernel::apply_force(colvarvalue const &force)
{
  size_t const deersize = times.size();
  cvm::rvector const u = dist_v.unit();
  cvm::vector1d<cvm::real> const &fv = force.vector1d_value;

  cvm::real totdeerforce = 0.0;
  for (int t = 0; t < deersize; t++){
     totdeerforce += fv[t] * deer_deriv[t];
  }

  if (!group1->noforce)
    group1->apply_force(-1.0 * u * totdeerforce);

  if (!group2->noforce)
    group2->apply_force(       u * totdeerforce);
}

cvm::real colvar::deer_kernel::dist2(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return (x1.vector1d_value - x2.vector1d_value).norm2();
}


colvarvalue colvar::deer_kernel::dist2_lgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return 2.0*colvarvalue((x1.vector1d_value - x2.vector1d_value), colvarvalue::type_vector);
}


colvarvalue colvar::deer_kernel::dist2_rgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return 2.0*colvarvalue((x2.vector1d_value - x1.vector1d_value), colvarvalue::type_vector);
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



colvar::deer::deer(std::string const &conf)
  : deer_kernel()
{
  function_type = "deer";
  mdepth = 1.0;
  alpha = 0.0;
  sample_dimensionality = 3;
  init(conf);
}


colvar::deer::deer()
  : deer_kernel()
{
  function_type = "deer";
  sample_dimensionality = 3;
  mdepth=1.0; //default value
  alpha=0.0; //default value
}


int colvar::deer::init(std::string const &conf)
{
  int error_code = COLVARS_OK;

  error_code |= deer_kernel::init(conf);

  get_keyval(conf, "deerMdepth", mdepth, mdepth);
  get_keyval(conf, "deerBackAlpha", alpha, alpha);
  get_keyval(conf, "sampleDimensionality", sample_dimensionality, sample_dimensionality);
  get_keyval(conf, "useDeerKernelWidth", use_kdeer_width, false);

  mdepth_unopt=mdepth; // save unoptimized parameters
  alpha_unopt=alpha; // save unoptimized parameters

  if (is_enabled(f_cvc_gradient)) {
    deriv_signal.resize(timesdeer.size());
    deer_signal_force.type(colvarvalue::type_vector);
    deer_signal_force.resize(timesdeer.size());
  }

  return error_code;
}


template<bool gradients, size_t deer_dim>
int colvar::deer::compute_exp_signal(cvm::vector1d<cvm::real> &kernel,
                                     cvm::vector1d<cvm::real> &kernel_deriv,
                                     cvm::vector1d<cvm::real> const &times,
                                     cvm::real deer_mdepth,
                                     cvm::real deer_alpha)
{
  size_t const nt = times.size();
  for (size_t it = 0; it < nt; it++) {
    cvm::real const t = times[it];
    cvm::real const exp_background = (deer_dim == 3) ?
      std::exp(-1.0 * deer_alpha*std::fabs(t)) :
      std::exp(-1.0 * std::pow(deer_alpha*std::fabs(t),
                               static_cast<cvm::real>(deer_dim)/3.0));
    cvm::real const k_t = kernel[it];
    if (gradients) {
      cvm::real const dk_t = kernel_deriv[it];
      cvm::real &dF_t = kernel_deriv[it];
      dF_t = deer_mdepth * exp_background * dk_t;
    }
    cvm::real &F_t = kernel[it];
    F_t = ((1.0 - deer_mdepth) + deer_mdepth*k_t) * exp_background;
  }
  return COLVARS_OK;
}

void colvar::deer::calc_value()
{
  deer_kernel::calc_value();
  size_t const dim = sample_dimensionality;
  // Compute in-place the experimental signal from the kernel
  if (is_enabled(f_cvc_gradient)) {
    compute_exp_signal<true, dim>(x.vector1d_value, deer_deriv, timesdeer,
                                  mdepth, alpha);
  } else {
    compute_exp_signal<false, dim>(x.vector1d_value, deer_deriv, timesdeer,
                                   mdepth, alpha);
  }
}

void colvar::deer::calc_gradients()
{
  // calculated on the fly in apply_force() and not stored
}


void colvar::deer::apply_force(colvarvalue const &force)
{
  // The signal's derivative has been already computed in deer_deriv
  deer_kernel::apply_force(force);
}

cvm::real colvar::deer_kernel::default_width() const
{
  return 0.1;
}

cvm::real colvar::deer::default_width() const
{
  use_kdeer_width=true;
  return 0.1;
}

cvm::real colvar::deer::scale_width(cvm::real const &refwidth) const
{
  if (!use_kdeer_width) {
    return refwidth;
  } else {
    cvm::real result = refwidth;
    size_t const varsize = timesdeer.size();
    cvm::real scalevalue = 0.0;
    for (int it = 0; it < varsize; it++){
      cvm::real const t = timesdeer[it];
      cvm::real const exp_background = (sample_dimensionality == 3) ?
        std::exp(-1.0 * alpha_unopt*std::fabs(t)) :
        std::exp(-1.0 * std::pow(alpha_unopt*std::fabs(t),
                                 static_cast<cvm::real>(sample_dimensionality)/3.0));
       scalevalue += (exp_background*exp_background*mdepth_unopt*mdepth_unopt);
    }
    scalevalue = std::sqrt(scalevalue/varsize);
    result = refwidth * scalevalue;
    return result;   
  }
}

cvm::real colvar::deer::rescale_width(cvm::real const &inputvalue) const
{
  cvm::real result = inputvalue;
  size_t const varsize = timesdeer.size();
  cvm::real scalevalue = 0.0;
  for (int it = 0; it < varsize; it++){
    cvm::real const t = timesdeer[it];
    cvm::real const exp_background = (sample_dimensionality == 3) ?
      std::exp(-1.0 * alpha_unopt*std::fabs(t)) :
      std::exp(-1.0 * std::pow(alpha_unopt*std::fabs(t),
                               static_cast<cvm::real>(sample_dimensionality)/3.0));
     scalevalue += (exp_background*exp_background*mdepth_unopt*mdepth_unopt);
  }
  scalevalue = std::sqrt(scalevalue/varsize);
  result = inputvalue/scalevalue;
  return result;
}

colvarvalue colvar::deer::rad_paramscale(colvarvalue const &inputvector) const
{
  size_t const varsize = inputvector.vector1d_value.size();
  cvm::vector1d<cvm::real> result = inputvector.vector1d_value;
  for (int it = 0; it < varsize; it++){
    cvm::real const t = timesdeer[it];
    cvm::real const exp_background = (sample_dimensionality == 3) ?
      std::exp(-1.0 * alpha*std::fabs(t)) :
      std::exp(-1.0 * std::pow(alpha*std::fabs(t),
                               static_cast<cvm::real>(sample_dimensionality)/3.0));
    result[it]=inputvector.vector1d_value[it]/(exp_background*mdepth);
  }
  return result;
}

void colvar::deer::get_params(vector1d<cvm::real> &vectorparams) const
{
  vectorparams.resize(2);
  vectorparams[0]=mdepth;
  vectorparams[1]=alpha;
  return;
}

void colvar::deer::set_params(vector1d<cvm::real> const &vectorparams)
{
  mdepth=vectorparams[0];
  alpha=vectorparams[1];
  return;
}

void colvar::deer_kernel::get_exp_val(colvarvalue &vectorexpval) const
{
  size_t const varsize = vectorexpval.vector1d_value.size();
  for (size_t it = 0; it < varsize; it++){
     vectorexpval.vector1d_value[it]=deerexpvalues[it]
  }
  return;
}

int colvar::deer::calc_params_deriv(colvarvalue const &lambdavector, colvarvalue const &centersvector,
                                     cvm::real const &coupling_time, cvm::real const &wt, cvm::real const &us,
                                     cvm::real const &width)
{
  size_t const varsize = centersvector.vector1d_value.size();
  cvm::real lambda2mod=0.0;
  cvm::real coef_mdepth=0.0;
  cvm::real grad_mdepth=0.0;
  cvm::real coef_alpha=0.0;
  cvm::real grad_alpha=0.0;
  for (size_t t=0;t<varsize;t++) {
     cvm::real const time = timesdeer[t];
     cvm::real const exp_background = (sample_dimensionality == 3) ?
       std::exp(-1.0 * alpha*std::fabs(time)) :
       std::exp(-1.0 * std::pow(alpha*std::fabs(time),
                                static_cast<cvm::real>(sample_dimensionality)/3.0));
     cvm::real const k_centers = 1.0+((centersvector.vector1d_value[t]-exp_background)/(exp_background*mdepth));
     cvm::real const lambda2=lambdavector.vector1d_value[t];
     lambda2mod+=std::fabs(exp_background*mdepth*lambda2);
     cvm::real const der_mdepth=(1.0-k_centers);
     coef_param[0]=coef_param[0]+(der_mdepth*der_mdepth);
     //grad_mdepth=grad_mdepth+lambda2*der_mdepth;
     deriv_param[t,0]=-exp_background*der_mdepth;
     cvm::real const alphaf = (sample_dimensionality == 3) ?
       std::fabs(time):
       (static_cast<cvm::real>(sample_dimensionality)/3.0)*(std::pow(std::fabs(time),
       static_cast<cvm::real>(sample_dimensionality)/3.0))*(std::pow(alpha,(static_cast<cvm::real>(sample_dimensionality)-3.0)/3.0));
     cvm::real const der_alpha = ((1-mdepth_deer[ii])+mdepth_deer[ii]*k_centers)*alphaf
     coef_param[1]=coef_param[1]+(der_alpha*der_alpha);
     //grad_alpha=grad_alpha+lambda2*der_alpha;
     deriv_param[t,1]=-exp_background*der_alpha;
  }
  lambda2mod=lambda2mod/varsize;

  cvm::real deltainv=std::sqrt(width)*lambda2mod; // scale factor so that coupling_time~one
  if(deltainv==0) deltainv=1.0;

  coef_param[0]=(mdepth*mdepth)/std::sqrt(coef_mdepth*coef_mdepth);
  coef_param[1]=(mdepth*mdepth)/std::sqrt(coef_alpha*coef_alpha);



  //mdepth=mdepth+coupling_time * coef_param[0] * grad_mdepth * wt  * cvm::dt()/(deltainv*us);
  //alpha=alpha+coupling_time * coef_param[1] * grad_alpha  * wt  * cvm::dt() /(deltainv*us);

  return COLVARS_OK;

}


int colvar::deer::update_params_rad(colvarvalue const &lambdavector, colvarvalue const &centersvector,
                                     cvm::real const &coupling_time, cvm::real const &wt, cvm::real const &us,
                                     cvm::real const &width)
{
  size_t const varsize = centersvector.vector1d_value.size();
  cvm::real lambda2mod=0.0;
  cvm::real coef_mdepth=0.0;
  cvm::real grad_mdepth=0.0;
  cvm::real coef_alpha=0.0;
  cvm::real grad_alpha=0.0;
  for (size_t t=0;t<varsize;t++) {
     cvm::real const time = timesdeer[t];
     cvm::real const exp_background = (sample_dimensionality == 3) ?
       std::exp(-1.0 * alpha*std::fabs(time)) :
       std::exp(-1.0 * std::pow(alpha*std::fabs(time),
                                static_cast<cvm::real>(sample_dimensionality)/3.0));
     cvm::real const k_centers = 1.0+((centersvector.vector1d_value[t]-exp_background)/(exp_background*mdepth));
     cvm::real const lambda2=exp_background*mdepth*lambdavector.vector1d_value[t];
     lambda2mod+=std::fabs(lambda2);
     //calculate derivatives of alpha_deer and mdepth_deer
     cvm::real const der_mdepth=(1.0-k_centers);
     coef_mdepth=coef_mdepth+(der_mdepth*der_mdepth);
     grad_mdepth=grad_mdepth+lambda2*der_mdepth;
     cvm::real const alphaf = (sample_dimensionality == 3) ?
       std::fabs(time):
       (static_cast<cvm::real>(sample_dimensionality)/3.0)*(std::pow(std::fabs(time),
       static_cast<cvm::real>(sample_dimensionality)/3.0))*(std::pow(alpha,(static_cast<cvm::real>(sample_dimensionality)-3.0)/3.0));
     cvm::real const der_alpha = ((1-mdepth_deer[ii])+mdepth_deer[ii]*k_centers)*alphaf
     coef_alpha=coef_alpha+(der_alpha*der_alpha);
     grad_alpha=grad_alpha+lambda2*der_alpha;
  }
  lambda2mod=lambda2mod/varsize;

  cvm::real deltainv=std::sqrt(width)*lambda2mod; // scale factor so that coupling_time~one
  if(deltainv==0) deltainv=1.0;

  coef_mdepth=std::sqrt(coef_mdepth*coef_mdepth);
  coef_alpha=std::sqrt(coef_alpha*coef_alpha);



  mdepth=mdepth-coupling_time * mdepth * grad_mdepth * wt  * cvm::dt()/(deltainv*us*coef_mdepth);
  alpha=alpha-coupling_time * mdepth *  grad_alpha  * wt  * cvm::dt() /(deltainv*us*coef_alpha);

  return COLVARS_OK;

}

int colvar::deer::update_params_rad_chis(colvarvalue const &aver_dev, colvarvalue const &exp_centers,
                                          size_t nsteps, cvm::real toll)
{
  size_t const varsize = exp_centers.vector1d_value.size();
  cvm::real grad_mdepth = 0;
  cvm::real grad_alpha = 0;
  cvm::real hess_mdepth = 0;
  cvm::real hess_alpha = 0;
  for (size_t it = 0; it < nsteps; it++){
     grad_mdepth = 0;
     grad_alpha = 0;
     hess_mdepth = 0;
     hess_alpha = 0;
     for (size_t t=0;t<varsize;t++) {
        cvm::real const time = timesdeer[t];
        cvm::real const exp_background = (sample_dimensionality == 3) ?
          std::exp(-1.0 * alpha*std::fabs(time)) :
          std::exp(-1.0 * std::pow(alpha*std::fabs(time),
                                   static_cast<cvm::real>(sample_dimensionality)/3.0));
        cvm::real const k_exp = 1.0+((exp_centers.vector1d_value[t]-exp_background)/(exp_background*mdepth));
        cvm::real const mean = aver_dev.vector1d_value[t]+k_exp;
        cvm::real const alphaf = (sample_dimensionality == 3) ?
          std::fabs(time):
          (static_cast<cvm::real>(sample_dimensionality)/3.0)*(std::pow(std::fabs(time),
          static_cast<cvm::real>(sample_dimensionality)/3.0))*(std::pow(alpha,(static_cast<cvm::real>(sample_dimensionality)-3.0)/3.0));

        cvm::real const der_fmean_mdepth=(1-mean)*exp_background;
        cvm::real const fmean = exp_background*(1-mdepth+mdepth*mean);
        cvm::real const der_fmean = -fmean*alphaf;
        cvm::real const sder_fmean =
         ((1/fmean)*der_fmean*der_fmean)+((1/alpha)*((static_cast<cvm::real>(sample_dimensionality)-3.0)/3.0)*der_fmean);
        grad_mdepth = grad_mdepth - aver_dev.vector1d_value[t]*
                                    der_fmean_mdepth*exp_background;
        grad_alpha = grad_alpha + aver_dev.vector1d_value[t]*
                                  der_fmean*exp_background;
        hess_mdepth = hess_mdepth + der_fmean_mdepth*der_fmean_mdepth;
        hess_alpha = hess_alpha + (der_fmean*der_fmean)+(aver_dev.vector1d_value[t]*sder_fmean*mdept*exp_background);
     }
     grad_mdepth = mdepth*grad_mdepth;
     grad_alpha = mdepth*grad_alpha;
     mdepth = mdepth - grad_mdepth/hess_mdepth;
     alpha = alpha - grad_alpha/hess_alpha;
     cvm::real const deltamdepth=sqrt(grad_mdepth*grad_mdepth)/hess_mdepth;
     cvm::real const deltamalpha=sqrt(grad_alphah*grad_alpha)/hess_alpha;
     if(deltamdepth/mdepth<toll&&deltalalpha/lalpha<toll) break;
  }
  return COLVARS_OK;
}

void colvar::deer::write_params_rad(std::ostream &os){

  os << " "
     << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
     << mdepth
     << " "
     << std::setprecision(cvm::cv_prec) << std::setw(cvm::cv_width)
     << alpha;

}

