// Check that an error raised without a module pointer (cvm::error_static(nullptr, ...))
// is caught at a module entry point and forwarded to the proxy, instead of terminating
// the program; the module must remain usable after a reset (this is what VMD relies on)

#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_stub.h"


/// Stub proxy that raises an error without a module pointer when a given atom is requested
class colvarproxy_stub_fatal : public colvarproxy_stub {
public:
  int init_atom(int atom_number) override
  {
    if (atom_number == 999) {
      return cvm::error_static(nullptr, "Error: test of a fatal error without a module pointer.\n",
                               COLVARS_BUG_ERROR);
    }
    return colvarproxy_stub::init_atom(atom_number);
  }
};


int main(int argc, char *argv[]) {

  colvarproxy_stub *proxy = new colvarproxy_stub_fatal();
  proxy->set_unit_system("real", false);
  colvarmodule *cvmodule = proxy->cvmodule;

  int retval = 0;

  // Parsing this configuration calls init_atom(999), which raises the error
  int err = cvmodule->read_config_string(R"(
colvar {
  name d
  distance {
    group1 {
      atomNumbers 1
    }
    group2 {
      atomNumbers 999
    }
  }
}
)");

  if (!(err & COLVARS_BUG_ERROR)) {
    std::cerr << "Error: the fatal error was not reported as COLVARS_BUG_ERROR (got "
              << err << ").\n";
    retval = 1;
  }

  if (!(cvmodule->get_error() & COLVARS_BUG_ERROR)) {
    std::cerr << "Error: the module error state was not set.\n";
    retval = 1;
  }

  // The module may be left in an inconsistent state, but must be usable after a reset
  cvmodule->reset();
  err = cvmodule->read_config_string("colvarsTrajFrequency 1\n");
  if (err != COLVARS_OK) {
    std::cerr << "Error: the module could not be used after a fatal error and a reset.\n";
    retval = 1;
  }

  delete proxy;

  return retval;
}
