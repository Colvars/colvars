// Test linking Colvars through a C interface

#ifdef __cplusplus
extern "C" {
#endif

// This will be callable from C or Fortran
// Returns proxy pointer and stores it to static global unique_colvarproxy_object
void *allocate_Colvars(char *);

// Skeleton of a function to access proxy members
void call_proxy_member_static();
void call_proxy_member(void *proxy);

// Functions of this type should define a complete interface
// This could be based on colvarscript

int run_colvarscript_command(int objc, unsigned char *const objv[]);
const char * get_colvarscript_result();

#ifdef __cplusplus
}
#endif
