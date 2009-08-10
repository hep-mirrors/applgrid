

#include "TFileString.h"
#include "TFile.h"

#include "appl_grid/appl_grid.h"
using namespace appl;



class git : public appl_pdf { 
  
public:
  git() : appl_pdf("git") { } 
  
  void evaluate(const double* f1, const double* f2, double* H) { } 
};




#include <sstream>
using std::ostringstream;
using std::stringstream;

void sgit(ostringstream& s) { 
  cerr << s << endl;
}

void sgit(ostream& s) { 
  cerr << s << endl;
}

void sgit(stringstream& s) { 
  cerr << s << endl;
}



int main() { 


  try { 

    sgit( ostringstream("hello") << 5 );

    //    appl_pdf* m = appl_pdf::getpdf("mcfm_pdf");
    appl_pdf* m = appl_pdf::getpdf("mcfm-w");
    appl_pdf* j = appl_pdf::getpdf("jetrad");
    appl_pdf* n = appl_pdf::getpdf("nlojet");

    double f1[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double f2[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double H[7];
    
    {
    git g;
    
    cout << m->Nproc() << endl;
    cout << j->Nproc() << endl;
    cout << n->Nproc() << endl;
    cout << g.Nproc() << endl;

    appl_pdf* _g = appl_pdf::getpdf("mcfm-z");
    _g->evaluate(f1, f2, H);

    }

    m->evaluate(f1, f2, H);

    appl_pdf* _g = appl_pdf::getpdf("mcfm-z");
    _g->evaluate(f1, f2, H);
    
  }
  catch ( appl_pdf::exception e) {   
    cout << "naughty, naughty" << endl;
  }
  

#if 0
  TFile f("git.root");

  TFileString g = *(TFileString*)f.Get("Transform");
  
  cout << g << endl;

  f.Close();
#endif


#if 0
  TFile f("git.root", "recreate");

  TFileString g("Transform");
  
  g.add("hello");

  g.Write();

  f.Close();
#endif

#if 0
  appl::grid g(6);
  
  g.Write("git.root");

  appl::grid g2("git.root");

  g2.untrim();
  cout << "g2 size=" << g2.size() << endl;

  g2.trim();
  cout << "g2 size=" << g2.size() << endl;
#endif
  

  return 0;
}
