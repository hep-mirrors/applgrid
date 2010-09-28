
#include <string.h>
#include <iostream>

extern "C" void mcfm_();

// use these to get the command line arguments into
// fortran 

static int    argc_ = 0;
static char** argv_ = 0;

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
  
  for ( int i=0 ; i<argc_ ; i++ ) { 
    std::cout << "argv[" << i << "] = " << argv_[i] << std::endl; 
  }

  mcfm_();

  return 0;
}
