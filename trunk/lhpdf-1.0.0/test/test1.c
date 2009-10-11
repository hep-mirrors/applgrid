/* Standard C includes */
#include <stdio.h>
#include <stdlib.h>


/* lhpdf includes */
#include "lhpdf.h"


#include <lhpdf_cteq6.h>
#include <lhpdf_mrst2001.h>
#include <lhpdf_fermi2002.h>
#include <lhpdf_alekhin.h>
#include <lhpdf_botje.h>



int main(int argc, char ** argv) 
{
  int j, p = 0, N = 41;
  double x, g1, g2, g3, g[13];
  lhpdf_pdf_t *pdf[41];
  
  if(argc > 1) N = atoi(argv[1]);

  for(j = 0; j < N; j++)
    pdf[j] = lhpdf_pdf_alloc(lhpdf_pdfset_botje_1000, j);
  
  for(j = 0; j < N; j++) {
    printf("\nPDF number : %d\n", j);
    printf("alpha_S(M_Z) = %.8f\n",lhpdf_pdf_evolveas(pdf[j], 61.71));
    
    for(p = 0; p <= 0; p++) {
      printf("    Parton label : %d\n", p);
      
      for(x = 0.1; x <= 1.0; x += 0.1) {
	lhpdf_pdf_evolvepdf(pdf[j], x, 10, g+6);
	g1 = g[p+6];
	
	lhpdf_pdf_evolvepdf(pdf[j], x, 100, g+6);
	g2 = g[p+6];
	
	lhpdf_pdf_evolvepdf(pdf[j], x, 1000.0, g+6);
	g3 = g[p+6];
	
	printf("%.2f  %.8f  %.8f  %.8f\n", x, g1, g2, g3);
      }
    }
  }
  
  for(j = 0; j < N; j++)
    lhpdf_pdf_free(pdf[j]);
  
  return 0;
}
