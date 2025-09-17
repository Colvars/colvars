#ifndef GLOBALMASTERCOLVARS_H
#define GLOBALMASTERCOLVARS_H

#include <memory>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_namd.h"


class GlobalMasterColvars : public GlobalMaster {
public:

  GlobalMasterColvars() : proxy(new colvarproxy_namd(this)) {}

  void calculate() override {
    proxy->calculate();
  }

protected:

  std::unique_ptr<colvarproxy_namd> cp;
};


#endif
