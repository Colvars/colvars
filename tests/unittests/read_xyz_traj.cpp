#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvarscript.h"
#include "colvartypes.h"

extern "C" int main(int argc, char *argv[]) {

  colvarproxy *proxy = new colvarproxy_stub();
  proxy->set_unit_system("real", false);

  const int natoms = 104;
  std::vector<std::vector<cvm::rvector>> traj;
  int err = 0, nframes = -1;

  while (!err) {
    nframes++;
    traj.push_back(std::vector<cvm::rvector>(natoms));
    err = proxy->colvars->load_coords_xyz("da-traj.xyz", &(traj[nframes]), nullptr, true);
  }
  traj.resize(nframes);

  for (int i = 0; i < nframes; i++) {
    cvm::log("Frame " + cvm::to_str(i) + " atom 0 coords: " + cvm::to_str(traj[i][0]));
  }
  return 0;
}
