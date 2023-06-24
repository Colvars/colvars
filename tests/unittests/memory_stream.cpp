#include <iostream>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarvalue.h"
#include "colvars_memstream.h"


template <typename T> std::ostream &operator<<(std::ostream &os, std::vector<T> const &v)
{
  os << "{";
  for (auto &vi : v)
    os << " " << vi;
  os << " }";
  return os;
}

template <typename T> void read_and_print(cvm::memory_stream &is, T &t)
{
  if (is >> t) {
    std::cout << typeid(T).name() << ": " << t << std::endl;
  }
}


int main(int argc, char const *argv[])
{
  cvm::real x = -101.0;
  size_t count = 1240566;
  std::vector<cvm::real> const a{1, 2, 3, 4, 5, 6, 7, 8};
  cvm::rvector const v = (cvm::rvector(-1.0, 1.0, 0.5)).unit();
  std::string const text("Vabbé / Está bien / Ça va / 好吧");
  std::vector<std::vector<double>> const a2(1, a);

  cvm::vector1d<cvm::real> v_from_a(a.size(), a.data());

  colvarvalue cv(v, colvarvalue::type_unit3vector);

  size_t n = (1L << 26);
  if (argc > 1) {
    n = atoi(argv[1]);
  }

  cvm::memory_stream buf(n);

  // Use standard functions and operators interchangeably
  buf.write_object(x);
  buf.write_object(count);
  buf.write_object(v);
  buf << text;
  buf << a;

  buf << v_from_a;

  buf << cv;

  // // The following will raise a compile error
  // buf << a2;

  if (buf) {

    cvm::real new_x = 0.0;
    read_and_print(buf, new_x);

    size_t new_count = 0;
    read_and_print(buf, new_count);

    cvm::rvector new_v;
    read_and_print(buf, new_v);

    std::string new_text;
    read_and_print(buf, new_text);

    std::vector<double> new_a;
    read_and_print(buf, new_a);

    cvm::vector1d<cvm::real> new_v_from_a;
    read_and_print(buf, new_v_from_a);

    colvarvalue new_cv(colvarvalue::type_unit3vector);
    read_and_print(buf, new_cv);

  }

  return buf ? 0 : 1;
}
