
#include <string.h>
#include <iostream>

extern "C" void mcfm_();

// use these to get the command line arguments into
// fortran 

static int    argc_;
static char** argv_;

extern "C" int cargc_() { return argc_; }

extern "C" void cargv_(int &i, char* s) {
  if ( i>=argc_ ) return;  
  strcpy(s, argv_[i]);
} 
extern "C" int clenargv_(int &i) {
  if ( i>=argc_ ) return 0;  
  return strlen(argv_[i]);
} 

// main starts here

int main(int argc, char** argv) { 

  argc_ = argc;
  argv_ = argv;
  
  mcfm_();

  return 0;
}
