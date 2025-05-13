#ifndef GLOBALMASTERCOLVARS_H
#define GLOBALMASTERCOLVARS_H

#include <memory>

#include "GlobalMaster.h"


class colvarproxy_namd;

class GlobalMasterColvars : public GlobalMaster {
public:

  GlobalMasterColvars();

  void calculate() override;

protected:

  std::unique_ptr<colvarproxy_namd> proxy;
};


#endif
