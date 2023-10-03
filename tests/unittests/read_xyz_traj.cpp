#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvarscript.h"
#include "colvartypes.h"

extern "C" int main(int argc, char *argv[]) {

  colvarproxy_stub *proxy = new colvarproxy_stub();
  proxy->set_unit_system("real", false);
  proxy->set_output_prefix("test.out");
  proxy->colvars->setup_input();
  proxy->colvars->setup_output();

  // Hard-coded for decaalanine system
  const int natoms = 104;
  for (int ai = 0; ai < natoms; ai++) {
    proxy->init_atom(ai+1);
  }

  proxy->colvars->read_config_file("test.in");

  int err = 0;
  while (!err) {
    cvm::log("Frame " + cvm::to_str(cvm::step_absolute()));
    err = proxy->read_frame_xyz("da-traj.xyz");
  }
  proxy->post_run();
  cvm::log("Done");

  return 0;
}
