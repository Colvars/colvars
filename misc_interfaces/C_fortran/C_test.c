// Test linking Colvars through a C interface

#include <stdio.h>
#include "colvarproxy_C_interface.h"

int main() {
    void *proxy = allocate_Colvars("This is a message from C");
    printf("This is a proxy pointer: %p\n", proxy);
    call_proxy_member_static();
    call_proxy_member(proxy);

    // Cannot yet call colvarscript without language interface
    unsigned char cmd[][256] = {"", "version" };
    run_colvarscript_command(1, (unsigned char **)(cmd));
    const char *res = get_colvarscript_result();
    printf("Colvarscript claims it runs version %s\n", res);
    return 0;
}
