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
  // parent class colvar_grid_scalar is constructed with margin option set to true
  // hence PMF grid is wider than gradient grid if non-PBC

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
/*
  if (nx[0]*nx[1] <= 100) {
    // Write explicit Laplacian FIXME debug output
    std::ofstream lap_out("lap_op.dat");
    std::vector<cvm::real> id(nx[0]*nx[1]), lap_col(nx[0]*nx[1]);
    for (int i = 0; i < nx[0] * nx[1]; i++) {
      id[i] = 1.;
      atimes(id, lap_col);
      id[i] = 0.;
      for (int j = 0; j < nx[0] * nx[1]; j++) {
        lap_out << cvm::to_str(i) + " " + cvm::to_str(j) + " " + cvm::to_str(lap_col[j]) << std::endl;
      }
      lap_out << std::endl;
    }
  }
  // FIXME*/

  nr_linbcg_sym(divergence, data, 1, tol, itmax, iter, err);
  cvm::log ("Completed integration in " + cvm::to_str (iter) + " steps with"
     + " error " + cvm::to_str (err));
/*
  // TODO remove this test of laplacian calcs
  std::vector<cvm::real> backup (data);
  std::ofstream p("pmf.dat");
  add_constant(-1.0 * minimum_value());
  write_multicol(p);
  std::vector<cvm::real> lap = std::vector<cvm::real>(data.size());
  atimes(data, lap);
  data = lap;
  std::ofstream l("laplacian.dat");
  write_multicol(l);
  data = divergence;
  std::ofstream d("divergence.dat");
  write_multicol(d);
  data = backup;
  // TODO TODO TODO*/

  return iter;
}


void integrate_potential::get_local_grads(colvar_grid_gradient * gradient, const std::vector<int> & ix0)
{
  int count;
  cvm::real   fact;
  const cvm::real   * g;
  std::vector<int> ix = ix0;
  bool edge;

  edge = gradient->wrap_edge(ix);
  if (!edge && (count = gradient->samples->value (ix))) {
    g = &(gradient->value (ix));
    fact = 1.0 / count;
    g11[0] = fact * g[0];
    g11[1] = fact * g[1];
  } else
    g11[0] = g11[1] = 0.0;

  ix[0] = ix0[0] - 1;
  edge = gradient->wrap_edge(ix);
  if (!edge && (count = gradient->samples->value (ix))) {
    g = & (gradient->value(ix));
    fact = 1.0 / count;
    g01[0] = fact * g[0];
    g01[1] = fact * g[1];
  } else
    g01[0] = g01[1] = 0.0;

  ix[1] = ix0[1] - 1;
  edge = gradient->wrap_edge(ix);
  if (!edge && (count = gradient->samples->value (ix))) {
    g = & (gradient->value(ix));
    fact = 1.0 / count;
    g00[0] = fact * g[0];
    g00[1] = fact * g[1];
  } else
    g00[0] = g00[1] = 0.0;

  ix[0] = ix0[0];
  edge = gradient->wrap_edge(ix);
  if (!edge && (count = gradient->samples->value (ix))) {
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
}

void integrate_potential::update_div(colvar_grid_gradient * gradient, const std::vector<int> &ix0)
{
  std::vector<int> ix (ix0);

  // If not periodic, expanded grid ensures that neighbors of ix0 are valid grid points
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

  get_local_grads(gradient, ix0);
  // Special case of corners: there is only one value of the gradient to average
  cvm::real fact_corner = 0.5;
  if (!periodic[0] && !periodic[1] && (ix0[0] == 0 || ix0[0] == nx[0]-1) && (ix0[1] == 0 || ix0[1] == nx[1]-1)) {
    fact_corner = 1.0;
  }
  divergence[linear_index] = (g10[0]-g00[0] + g11[0]-g01[0]) * fact_corner / widths[0]
                           + (g01[1]-g00[1] + g11[1]-g10[1]) * fact_corner / widths[1];
 }


/// Multiplication by sparse matrix representing Laplacian
void integrate_potential::atimes(const std::vector<cvm::real> &p, std::vector<cvm::real> &l)
{
  size_t index, index2;

  const cvm::real fx = 1.0/widths[0];
  const cvm::real fy = 1.0/widths[1];
  const cvm::real ffx = fx * fx;
  const cvm::real ffy = fy * fy;
  // offsets for 4 reference points of the Laplacian stencil

  const int h = nx[1];
  const int w = nx[0];
  int xm = -h;
  int xp =  h;
  int ym = -1;
  int yp =  1;

  index = h + 1;
  for (int i=1; i<w-1; i++) {
    for (int j=1; j<h-1; j++) {
      l[index] = ffx * (p[index + xm] + p[index + xp] - 2.0 * p[index])
               + ffy * (p[index + ym] + p[index + yp] - 2.0 * p[index]);
      index++;
    }
    index += 2; // skip the edges and move to next column
  }

  // then, edges depending on BC
  if (periodic[0]) {
    // i = 0 and i = w are periodic images
    xm =  h * (w - 1);
    xp =  h;
    ym = -1;
    yp =  1;
    index = 1; // Follows left edge
    index2 = h * (w - 1) + 1; // Follows right edge
    for (int j=1; j<h-1; j++) {
      l[index] = ffx * (p[index + xm] + p[index + xp] - 2.0 * p[index])
               + ffy * (p[index + ym] + p[index + yp] - 2.0 * p[index]);
      l[index2] = ffx * (p[index2 - xp] + p[index2 - xm] - 2.0 * p[index2])
                + ffy * (p[index2 + ym] + p[index2 + yp] - 2.0 * p[index2]);
      index++;
      index2++;
    }
  } else {
    xm =  -h;
    xp =  h;
    ym = -1;
    yp =  1;
    index = 1; // Follows left edge
    index2 = h * (w - 1) + 1; // Follows right edge
    for (int j=1; j<h-1; j++) {
      // x gradient beyond the edge is taken to be zero
      // alternate: x gradient + y term of laplacian
      l[index] = ffx * (p[index + xp] - p[index])
               + ffy * (p[index + ym] + p[index + yp] - 2.0 * p[index]);
      l[index2] = ffx * (p[index2 + xm] - p[index2])
                + ffy * (p[index2 + ym] + p[index2 + yp] - 2.0 * p[index2]);
      index++;
      index2++;
    }
  }

  if (periodic[1]) {
    // j = 0 and j = h are periodic images
    xm = -h;
    xp =  h;
    ym =  h - 1;
    yp =  1;
    index = h; // Follows bottom edge
    index2 = 2 * h - 1; // Follows top edge
    for (int i=1; i<w-1; i++) {
      l[index] = ffx * (p[index + xm] + p[index + xp] - 2.0 * p[index])
               + ffy * (p[index + ym] + p[index + yp] - 2.0 * p[index]);
      l[index2] = ffx * (p[index2 + xm] + p[index2 + xp] - 2.0 * p[index2])
                + ffy * (p[index2 - yp] + p[index2 - ym] - 2.0 * p[index2]);
      index  += h;
      index2 += h;
    }
  } else {
    xm = -h;
    xp =  h;
    ym =  -1;
    yp =  1;
    index = h; // Follows bottom edge
    index2 = 2 * h - 1; // Follows top edge
    for (int i=1; i<w-1; i++) {
      // alternate: y gradient + x term of laplacian
      l[index] = ffx * (p[index + xm] + p[index + xp] - 2.0 * p[index])
             +ffy * (p[index + yp] - p[index]);
      l[index2] = ffx * (p[index2 + xm] + p[index2 + xp] - 2.0 * p[index2])
             +ffy * (p[index2 + ym] - p[index2]);
      index  += h;
      index2 += h;
    }
  }

  // 4 corners
  xm = h;
  xp = h * (w - 1);
  ym = 1;
  yp = h - 1;
  cvm::real lx, ly;

  index = 0;
  lx = periodic[0] ? (p[xp] + p[xm] - 2.0 * p[0]) : (p[h] - p[0]);
  ly = periodic[1] ?  (p[yp] + p[ym] - 2.0 * p[0]) : (p[1] - p[0]);
  l[index] = ffx * lx + ffy * ly;

  index = h-1;
  lx = periodic[0] ? (p[index + xp] + p[index + xm] - 2.0 * p[index]) : (p[index + h] - p[index]);
  ly = periodic[1] ? (p[index - ym] + p[index - yp] - 2.0 * p[index]) : (p[index - 1] - p[index]);
  l[index] = ffx * lx + ffy * ly;

  index = h * (w-1);
  lx = periodic[0] ? (p[index - xm] + p[index - xp] - 2.0 * p[index]) : (p[index - h] - p[index]);
  ly = periodic[1] ? (p[index + yp] + p[index + ym] - 2.0 * p[index]) : (p[index + 1] - p[index]);
  l[index] = ffx * lx + ffy * ly;

  index = h * w - 1;
  lx = periodic[0] ? (p[index - xm] + p[index - xp] - 2.0 * p[index]) : (p[index - h] - p[index]);
  ly = periodic[1] ? (p[index - ym] + p[index - yp] - 2.0 * p[index]) : (p[index - 1] - p[index]);
  l[index] = ffx * lx + ffy * ly;
}

// TODO remove if unused
/// Inversion of preconditioner matrix (e.g. diagonal of the Laplacian)
void integrate_potential::asolve(const std::vector<cvm::real> &b, std::vector<cvm::real> &x, const int itrnsp)
{
//   cvm::real inv_diag = - 1.0 / (2.0 / (widths[0] * widths[0]) + 2.0 / (widths[1] * widths[1]));
  // Identity works, so we leave it at that for now.
  for (size_t i=0; i<x.size(); i++) {
//     x[i] = inv_diag * b[i];
    x[i] = b[i];
  }
  return;
}


// b : RHS of equation
// x : initial guess for the solution; output is solution
// itol : convergence criterion
void integrate_potential::nr_linbcg_sym(const std::vector<cvm::real> &b, std::vector<cvm::real> &x, const int itol, const cvm::real tol,
  const int itmax, int &iter, cvm::real &err)
{
  cvm::real ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
  const cvm::real EPS=1.0e-14;
  int j;

  // Making it symmetric

  iter=0;
  atimes(x,r);
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
    atimes(p,z);
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
//     std::cout << "iter=" << std::setw(4) << iter+1 << std::setw(12) << err << std::endl;
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
