// Test linking Colvars through a C interface

#include "colvarproxy_C.h"
#include "colvarproxy_C_interface.h"
#include "colvarscript.h"

// This to avoid having to pass around proxy pointers
static colvarproxy* unique_colvarproxy_object;

colvarproxy_C::colvarproxy_C()
{
    std::cerr << "This is the colvarproxy_C constructor at address " << this << std::endl;
    colvars = new colvarmodule(this);
    colvars->log("This is the Module speaking.");
}

colvarproxy_C::~colvarproxy_C()
{
}

void *allocate_Colvars(char *t)
{
    std::cout << "A message from our caller: " << t << std::endl;
    colvarproxy_C *proxy = new colvarproxy_C();
    unique_colvarproxy_object = proxy;
    return proxy;
}

// This prototype C function can call any colvarproxy member by using the global static pointer
void call_proxy_member_static()
{
    unique_colvarproxy_object->log("Called through static pointer.");
}

// This prototype C function can call any colvarproxy member by using the global static pointer
void call_proxy_member(void *proxy_in)
{
    colvarproxy *proxy = static_cast<colvarproxy *>(proxy_in);
    proxy->log("Called through pointer passed as argument.");
}
