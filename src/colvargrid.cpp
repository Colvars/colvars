// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvargrid.h"


colvar_grid_count::colvar_grid_count()
  : colvar_grid<size_t>()
{
  mult = 1;
}

colvar_grid_count::colvar_grid_count(std::vector<int> const &nx_i,
                                     size_t const           &def_count)
  : colvar_grid<size_t>(nx_i, def_count, 1)
{}

colvar_grid_count::colvar_grid_count(std::vector<colvar *>  &colvars,
                                     size_t const           &def_count)
  : colvar_grid<size_t>(colvars, def_count, 1)
{}

colvar_grid_scalar::colvar_grid_scalar()
  : colvar_grid<cvm::real>(), samples(NULL), grad(NULL)
{}

colvar_grid_scalar::colvar_grid_scalar(colvar_grid_scalar const &g)
  : colvar_grid<cvm::real>(g), samples(NULL), grad(NULL)
{
  grad = new cvm::real[nd];
}

colvar_grid_scalar::colvar_grid_scalar(std::vector<int> const &nx_i)
  : colvar_grid<cvm::real>(nx_i, 0.0, 1), samples(NULL), grad(NULL)
{
  grad = new cvm::real[nd];
}

colvar_grid_scalar::colvar_grid_scalar(std::vector<colvar *> &colvars, bool margin)
  : colvar_grid<cvm::real>(colvars, 0.0, 1, margin), samples(NULL), grad(NULL)
{
  grad = new cvm::real[nd];
}

colvar_grid_scalar::~colvar_grid_scalar()
{
  if (grad) {
    delete [] grad;
    grad = NULL;
  }
}

cvm::real colvar_grid_scalar::maximum_value() const
{
  cvm::real max = data[0];
  for (size_t i = 0; i < nt; i++) {
    if (data[i] > max) max = data[i];
  }
  return max;
}


cvm::real colvar_grid_scalar::minimum_value() const
{
  cvm::real min = data[0];
  for (size_t i = 0; i < nt; i++) {
    if (data[i] < min) min = data[i];
  }
  return min;
}

cvm::real colvar_grid_scalar::minimum_pos_value() const
{
  cvm::real minpos = data[0];
  size_t i;
  for (i = 0; i < nt; i++) {
    if(data[i] > 0) {
      minpos = data[i];
      break;
    }
  }
  for (i = 0; i < nt; i++) {
    if (data[i] > 0 && data[i] < minpos) minpos = data[i];
  }
  return minpos;
}

cvm::real colvar_grid_scalar::integral() const
{
  cvm::real sum = 0.0;
  for (size_t i = 0; i < nt; i++) {
    sum += data[i];
  }
  cvm::real bin_volume = 1.0;
  for (size_t id = 0; id < widths.size(); id++) {
    bin_volume *= widths[id];
  }
  return bin_volume * sum;
}


cvm::real colvar_grid_scalar::entropy() const
{
  cvm::real sum = 0.0;
  for (size_t i = 0; i < nt; i++) {
    sum += -1.0 * data[i] * std::log(data[i]);
  }
  cvm::real bin_volume = 1.0;
  for (size_t id = 0; id < widths.size(); id++) {
    bin_volume *= widths[id];
  }
  return bin_volume * sum;
}


colvar_grid_gradient::colvar_grid_gradient()
  : colvar_grid<cvm::real>(), samples(NULL)
{}

colvar_grid_gradient::colvar_grid_gradient(std::vector<int> const &nx_i)
  : colvar_grid<cvm::real>(nx_i, 0.0, nx_i.size()), samples(NULL)
{}

colvar_grid_gradient::colvar_grid_gradient(std::vector<colvar *> &colvars)
  : colvar_grid<cvm::real>(colvars, 0.0, colvars.size()), samples(NULL)
{}

void colvar_grid_gradient::write_1D_integral(std::ostream &os)
{
  cvm::real bin, min, integral;
  std::vector<cvm::real> int_vals;

  os << "#       xi            A(xi)\n";

  if ( cv.size() != 1 ) {
    cvm::error("Cannot write integral for multi-dimensional gradient grids.");
    return;
  }

  integral = 0.0;
  int_vals.push_back( 0.0 );
  min = 0.0;

  // correction for periodic colvars, so that the PMF is periodic
  cvm::real corr;
  if ( periodic[0] ) {
    corr = average();
  } else {
    corr = 0.0;
  }

  for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {

    if (samples) {
      size_t const samples_here = samples->value(ix);
      if (samples_here)
        integral += (value(ix) / cvm::real(samples_here) - corr) * cv[0]->width;
    } else {
      integral += (value(ix) - corr) * cv[0]->width;
    }

    if ( integral < min ) min = integral;
    int_vals.push_back( integral );
  }

  bin = 0.0;
  for ( int i = 0; i < nx[0]; i++, bin += 1.0 ) {
    os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
       << std::setw(cvm::cv_width)
       << std::setprecision(cvm::cv_prec)
       << int_vals[i] - min << "\n";
  }

  os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
     << std::setw(cvm::cv_width)
     << std::setprecision(cvm::cv_prec)
     << int_vals[nx[0]] - min << "\n";

  return;
}






// Parameters:
// b (divergence + BC): member of class integrate_cg; updated locally every ts
// x (solution PMF): reference to pmf object? or copy of the vector if more efficient
// atimes, asolve: member functions of class integrate_cg, relying on
// laplacian: member data (vector) of integrate_cg; sparse matrix representation of
// finite diff. Laplacian, defined once and for all at construction time.
// NOTE: most of the data needs complete updates if the grid size changes...

integrate_potential::integrate_potential(std::vector<colvar *> &colvars)
  : colvar_grid_scalar (colvars, true)
{

  // length of "vectors" is number of points in potential (data) and divergence grids
  p.reserve(nt);
  p.assign(nt, 0.0);
  pp.reserve(nt);
  pp.assign(nt, 0.0);
  r.reserve(nt);
  r.assign(nt, 0.0);
  rr.reserve(nt);
  rr.assign(nt, 0.0);
  z.reserve(nt);
  z.assign(nt, 0.0);
  zz.reserve(nt);
  zz.assign(nt, 0.0);

  divergence.resize(nt);
}


int integrate_potential::integrate(const int itmax, const cvm::real tol, cvm::real & err)
{
  int iter;

  nr_linbcg(divergence, data, 1, tol, itmax, iter, err);
  cvm::log ("Completed integration in " + cvm::to_str (iter) + " steps with"
     + " error " + cvm::to_str (err));

  // TODO remove this test of laplacian calcs
  std::ofstream p("pmf.dat");
  write_multicol(p);
  std::vector<cvm::real> lap = std::vector<cvm::real>(data.size());
  atimes(data, lap, 0);
  data = lap;
  std::ofstream l("laplacian.dat");
  write_multicol(l);
  data = divergence;
  std::ofstream d("divergence.dat");
  write_multicol(d);
  // TODO TODO TODO

  return iter;
}


void integrate_potential::get_local_grads(colvar_grid_gradient * gradient, const std::vector<int> & ix0)
{
  size_t      count;
  cvm::real   fact;
  const cvm::real   * g;
  std::vector<int> ix = ix0;

  gradient->wrap(ix);
  count = gradient->samples->value (ix);
  if (count) {
    g = &(gradient->value (ix));
    fact = 1.0 / count;
    g11[0] = fact * g[0];
    g11[1] = fact * g[1];
  } else
    g11[0] = g11[1] = 0.0;

  ix[0] = ix0[0] - 1;
  gradient->wrap(ix);
  count = gradient->samples->value(ix);
  if (count) {
    g = & (gradient->value(ix));
    fact = 1.0 / count;
    g01[0] = fact * g[0];
    g01[1] = fact * g[1];
  } else
    g01[0] = g01[1] = 0.0;

  ix[1] = ix0[1] - 1;
  gradient->wrap(ix);
  count = gradient->samples->value(ix);
  if (count) {
    g = & (gradient->value(ix));
    fact = 1.0 / count;
    g00[0] = fact * g[0];
    g00[1] = fact * g[1];
  } else
    g00[0] = g00[1] = 0.0;

  ix[0] = ix0[0];
  gradient->wrap(ix);
  count = gradient->samples->value(ix);
  if (count) {
    g = & (gradient->value(ix));
    fact = 1.0 / count;
    g10[0] = fact * g[0];
    g10[1] = fact * g[1];
  } else
    g10[0] = g10[1] = 0.0;
}

void integrate_potential::set_div(colvar_grid_gradient * gradient)
{
  for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix)) {
    update_div_local(gradient, ix);
  }

  int index = 0;
  if (cvm::debug()) {
    for (int i=0; i<nx[0]; i++) {
      for (int j=0; j<nx[1]; j++) {
        printf("%i %i %8.3f\n", i, j, divergence[index++]);
      }
      printf("\n");
    }
  }
}

void integrate_potential::update_div(colvar_grid_gradient * gradient, const std::vector<int> &ix0)
{
  std::vector<int> ix (ix0);

  update_div_local(gradient, ix);
  ix[0]++; wrap(ix);
  update_div_local(gradient, ix);
  ix[1]++; wrap(ix);
  update_div_local(gradient, ix);
  ix[0]--; wrap(ix);
  update_div_local(gradient, ix);
}

void integrate_potential::update_div_local(colvar_grid_gradient * gradient, const std::vector<int> &ix0)
{
  std::vector<int> ix (2);
  const int linear_index = nx[1] * ix0[0] + ix0[1];

  // 4 corners in non-periodic grids have divergence set to 0
  if ( ! (periodic[0] || periodic[1])
      && (ix0[0] == 0 || ix0[0] == nx[0]-1)
      && (ix0[1] == 0 || ix0[1] == nx[1]-1) ) {
    divergence[linear_index] = 0.0;
    return;
  }

  // if ix[i] = 0 or max, update edge...
  if ((ix0[0] == 0 || ix0[0] == nx[0]-1) && !periodic[0]) {
    // NOTE gradient grid is smaller than divergence grid by 1
    ix[0] = ix0[0] == 0 ? 0 : nx[0] - 2;
    ix[1] = ix0[1]-1;
    g00[1] = gradient->value_output(ix, 1);
    ix[1] = ix0[1];
    g01[1] = gradient->value_output(ix, 1);
    divergence[linear_index] = (g01[1]-g00[1]) / widths[1];
    return;
  }

  if ((ix0[1] == 0 || ix0[1] == nx[1]-1) && !periodic[1])  {
    ix[0] = ix0[0]-1;
    // NOTE gradient grid is smaller than divergence grid by 1
    ix[1] = ix0[1] == 0 ? 0 : nx[1] - 2;
    g00[0] = gradient->value_output(ix, 0);
    ix[0] = ix0[0];
    g10[0] = gradient->value_output(ix, 0);
    divergence[linear_index] = (g10[0]-g00[0]) / widths[0];
    return;
  }

  // FIXME: missing case of edges in semi-periodic case

  // otherwise update "center" (if periodic, everything is in the center)
  get_local_grads(gradient, ix0);
  divergence[linear_index] = (g10[0]-g00[0] + g11[0]-g01[0]) * 0.5 / widths[0]
                           + (g01[1]-g00[1] + g11[1]-g10[1]) * 0.5 / widths[1];
 }



/// Multiplication by sparse matrix representing Laplacian (or its transpose)
void integrate_potential::atimes (const std::vector<cvm::real> &x, std::vector<cvm::real> &r, const int transpose)
{
  size_t index, index2;

  // Note: center of the matrix is symmetric - we only worry about transposing the edges in non-PBC.

  const cvm::real fx = 1.0/widths[0];
  const cvm::real fy = 1.0/widths[1];
  const cvm::real ffx = fx * fx;
  const cvm::real ffy = fy * fy;
  // offsets for 4 reference points of the Laplacian
  int xm = -nx[1];
  int xp =  nx[1];
  int ym = -1;
  int yp =  1;

  index = nx[1] + 1;
  for (int i=1; i<nx[0]-1; i++) {
    for (int j=1; j<nx[1]-1; j++) {
      r[index] = ffx * (x[index + xm] + x[index + xp] - 2.0 * x[index])
               + ffy * (x[index + ym] + x[index + yp] - 2.0 * x[index]);
      index++;
    }
    index += 2;
  }

  // then, edges depending on BC
  if (periodic[0]) {
    // i = 0 and i = nx[0] are periodic images
    xm =  nx[1] * (nx[0] - 1);
    xp =  nx[1];
    ym = -1;
    yp =  1;
    index = 1;
    index2 = nx[1] * (nx[0] - 1) + 1;
    for (int j=1; j<nx[1]-1; j++) {
      r[index] = ffx * (x[index + xm] + x[index + xp] - 2.0 * x[index])
               + ffy * (x[index + ym] + x[index + yp] - 2.0 * x[index]);
      r[index2] = ffx * (x[index2 - xp] + x[index2 - xm] - 2.0 * x[index2])
                + ffy * (x[index2 + ym] + x[index2 + yp] - 2.0 * x[index2]);
      index++;
      index2++;
    }
  } else {
    //TODO
  }
  if (periodic[1]) {
    // j = 0 and j = nx[1] are periodic images
    xm = -nx[1];
    xp =  nx[1];
    ym =  nx[1] - 1;
    yp =  1;
    index = nx[1];
    index2 = 2 * nx[1] - 1;
    for (int i=1; i<nx[0]-1; i++) {
      r[index] = ffx * (x[index + xm] + x[index + xp] - 2.0 * x[index])
               + ffy * (x[index + ym] + x[index + yp] - 2.0 * x[index]);
      r[index2] = ffx * (x[index2 + xm] + x[index2 + xp] - 2.0 * x[index2])
                + ffy * (x[index2 - yp] + x[index2 - ym] - 2.0 * x[index2]);
      index  += nx[1];
      index2 += nx[1];
    }
  } else {
    //TODO
  }

  //FIXME
  if (! (periodic[0] || periodic[1])) {
    r[0] = 0.0;
    r[nx[1]-1] = 0.0;
    r[nx[1] * (nx[0]-1)] = 0.0;
    r[nx[1] * nx[0] - 1] = 0.0;
  }
  // TODO: corners in semi-PBC
  if (periodic[0] && periodic[1]) {
    xm = nx[1];
    xp = nx[1] * (nx[0] - 1);
    ym = 1;
    yp = nx[1] - 1;
    index = 0;
    r[index] = ffx * (x[xp] + x[xm] - 2.0 * x[0])
             + ffy * (x[yp] + x[ym] - 2.0 * x[0]);
    index = nx[1]-1;
    r[index] = ffx * (x[index + xp] + x[index + xm] - 2.0 * x[index])
             + ffy * (x[index - ym] + x[index - yp] - 2.0 * x[index]);
    index = nx[1] * (nx[0]-1);
    r[index] = ffx * (x[index - xm] + x[index - xp] - 2.0 * x[index])
             + ffy * (x[index + yp] + x[index + ym] - 2.0 * x[index]);
    index = nx[1] * nx[0] - 1;
    r[index] = ffx * (x[index - xm] + x[index - xp] - 2.0 * x[index])
             + ffy * (x[index - ym] + x[index - yp] - 2.0 * x[index]);
  }
}

/// Inversion of preconditioner matrix (e.g. diagonal of the Laplacian)
void integrate_potential::asolve(const std::vector<cvm::real> &b, std::vector<cvm::real> &x, const int itrnsp)
{
  //cvm::real inv_diag = - 1.0 / (2.0 / (widths[0] * widths[0]) + 2.0 / (widths[1] * widths[1]));
  // Identity works, so we leave it at that for now.
  for (size_t i=0; i<x.size(); i++) {
    //x[i] = inv_diag * b[i];
    x[i] = b[i];
  }
  return;
}


// b : RHS of equation
// x : initial guess for the solution; output is solution
// itol : convergence criterion
void integrate_potential::nr_linbcg(const std::vector<cvm::real> &b, std::vector<cvm::real> &x, const int itol, const cvm::real tol,
  const int itmax, int &iter, cvm::real &err)
{
  cvm::real ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
  const cvm::real EPS=1.0e-14;
  int j;

  // Making it symmetric

  iter=0;
  atimes(x,r,0);
  for (j=0;j<nt;j++) {
    r[j]=b[j]-r[j];
    // SYM rr[j]=r[j];
  }
  //atimes(r,rr,0);   // uncomment for minimum residual algorithm
  if (itol == 1) {
    bnrm=nr_snrm(b,itol);
    if (bnrm < EPS) {// Target is zero - somehow this is a problem
      return;
    }
    asolve(r,z,0);
  }
  else if (itol == 2) {
    asolve(b,z,0);
    bnrm=nr_snrm(z,itol);
    asolve(r,z,0);
  }
  else if (itol == 3 || itol == 4) {
    asolve(b,z,0);
    bnrm=nr_snrm(z,itol);
    asolve(r,z,0);
    znrm=nr_snrm(z,itol);
  } else {
    cvm::fatal_error("Illegal itol in linbcg");
  }
  //std::cout << std::fixed << std::setprecision(6);
  while (iter < itmax) {
    ++iter;
    // SYM asolve(rr,zz,1);
    for (bknum=0.0,j=0;j<nt;j++) {
      // SYM bknum += z[j]*rr[j];
      bknum += z[j]*r[j];
    }
    if (iter == 1) {
      for (j=0;j<nt;j++) {
        p[j]  = z[j];
        // SYM pp[j] = zz[j];
      }
    } else {
      bk=bknum/bkden;
      for (j=0;j<nt;j++) {
        p[j]  = bk*p[j]  + z[j];
        // SYM pp[j] = bk*pp[j] + zz[j];
      }
    }
    bkden = bknum;
    atimes(p,z,0);
    for (akden=0.0,j=0;j<nt;j++) {
      // SYM akden += z[j]*pp[j];
      akden += z[j]*p[j];
    }
    ak = bknum/akden;
    // SYM atimes(pp,zz,1);
    for (j=0;j<nt;j++) {
      x[j]  += ak*p[j];
      r[j]  -= ak*z[j];
      // SYM rr[j] -= ak*zz[j];
    }
    asolve(r,z,0);
    if (itol == 1)
      err = nr_snrm(r,itol)/bnrm;
    else if (itol == 2)
      err = nr_snrm(z,itol)/bnrm;
    else if (itol == 3 || itol == 4) {
      zm1nrm = znrm;
      znrm = nr_snrm(z,itol);
      if (fabs(zm1nrm-znrm) > EPS*znrm) {
        dxnrm = fabs(ak)*nr_snrm(p,itol);
        err = znrm/fabs(zm1nrm-znrm)*dxnrm;
      } else {
        err = znrm/bnrm;
        continue;
      }
      xnrm = nr_snrm(x,itol);
      if (err <= 0.5*xnrm) err /= xnrm;
      else {
        err = znrm/bnrm;
        continue;
      }
    }
    std::cout << "iter=" << std::setw(4) << iter+1 << std::setw(12) << err << std::endl;
    if (err <= tol)
      break;
  }
}

cvm::real integrate_potential::nr_snrm(const std::vector<cvm::real> &sx, const int itol)
{
  int i,isamax;
  cvm::real ans;

  int n=sx.size();
  if (itol <= 3) {
    ans = 0.0;
    for (i=0;i<n;i++) ans += sx[i]*sx[i];
    return ::sqrt(ans);
  } else {
    isamax=0;
    for (i=0;i<n;i++) {
      if (::fabs(sx[i]) > ::fabs(sx[isamax])) isamax=i;
    }
    return ::fabs(sx[isamax]);
  }
}