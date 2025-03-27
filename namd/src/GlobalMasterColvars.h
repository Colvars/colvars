#ifndef GLOBALMASTERCOLVARS_H
#define GLOBALMASTERCOLVARS_H

#include <memory>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_namd.h"


class GlobalMasterColvars : public GlobalMaster {
public:

  GlobalMasterColvars() : cp(new colvarproxy_namd(this)) {}

  void calculate() override {
    cp->calculate();
  }

protected:
  
  std::unique_ptr<colvarproxy_namd> cp;
};


#endif
