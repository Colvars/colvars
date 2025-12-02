#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvarscript.h"
#include "colvartypes.h"

int main(int argc, char *argv[]) {

  colvarproxy_stub *proxy = new colvarproxy_stub();
  proxy->set_unit_system("real", false);
  proxy->set_output_prefix("test.out");
  proxy->cvmodule->setup_input();
  proxy->cvmodule->setup_output();

  // Hard-coded for decaalanine system
  const int natoms = 104;
  for (int ai = 0; ai < natoms; ai++) {
    proxy->init_atom(ai+1);
  }

  int err = proxy->cvmodule->read_config_file("test_error.in");

  delete proxy;

  if (err != COLVARS_OK) {
    return 0;
  } else {
    return 1;
  }
}
