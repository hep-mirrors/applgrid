#include <iostream>
#include <string>
#include "TH2D.h"
#include "TFile.h"
#include "grid.h"
// #include "pdf-cteq6.h"

using std::cout;
using std::cerr;
using std::endl;



void make_pdf_combinations(const double* fA, const double* fB, double* H) {

  double GA=fA[6];
  double GB=fB[6];
  double QA=0; double QB=0; double QbA=0; double QbB=0; double D=0; double Db=0;
  for(int i=1;i<7;i++) {
    QA+=fA[i+6];
    QB+=fB[i+6];
  }
  for(int i=-6;i<0;i++) {
    QbA+=fA[i+6];
    QbB+=fB[i+6];
  }
  for(int i=-6;i<7;i++) {
    if(i!=0) {
      D+=fA[i+6]*fB[i+6];
      Db+=fA[i+6]*fB[-i+6];
    }
  }
  H[0]=GA*GB;    // gluon a gluon b 
  H[1]=(QA+QbA)*GB; // q+qbar a + gluon b
  H[2]=GA*(QB+QbB);  // gluon a + a+abar b ( H[1](x) == H[2](x) 
  H[3]=QA*QB+QbA*QbB-D;
  H[4]=D;
  H[5]=Db;
  H[6]=QA*QbB+QbA*QB-Db;
}



int main(int argc, char **argv) {
  //  grid* mygrid = new grid("weightgrid.root");
  grid* mygrid = new grid("grid2.root");


  TH1D* grid[2][7];
  char histname[30];
  for(int j=0;j<2;j++) {
    for(int i=0;i<7;i++) {
      sprintf(histname,"grid_%i_%i",j,i);
      grid[j][i] = new TH1D(histname,histname,mygrid->Nobs(),mygrid->obsmin(),mygrid->obsmax());
    }
  }
  TH1D* grid_total = new TH1D("grid_total","grid_total",mygrid->Nobs(),mygrid->obsmin(),mygrid->obsmax());

  TFile myfile2("weightgrid.root");
  TH1D* run = (TH1D*) myfile2.Get("run");
  run->SetDirectory(0);

  TH1D* reference[2][7];
  for(int j=0; j<2; j++) {
    for(int i=0; i<7; i++) {
      sprintf(histname,"reference_%i_%i",j,i);
      reference[j][i] = (TH1D*) myfile2.Get(histname);
      reference[j][i]->SetDirectory(0);
    }
  }

  myfile2.Close();
  double runs=run->GetBinContent(1);



  //  pdf_cteq6* nlojetpdf = new pdf_cteq6(0);

  double alphas, alphas2, alphas3, Q2,Q, born_weight, nlo_weight, x, x1, x2, xA, xB;
  double H[7];
  TH1D* m_pdf[13];
  double fA[13]; double fB[13]; double f[13];

  double pb_fac=3.89379656e8;
  //double pb_fac=1.0;

  //  cout<<"Cross section after "<<int(runs+0.5)<<" runs..."<<endl;


  for(int iobs=0;iobs<mygrid->Nobs();iobs++) {
    
    //    cout << "iobs=" << iobs << " (Nobs=" << mygrid->Nobs() << ")" << endl;

    for(int iinitial=0;iinitial<7;iinitial++) {
      //      cout << iinitial << "  " << mygrid->m_weight[0][iinitial][iobs] << endl;
      mygrid->m_weight[0][iinitial][iobs]->Scale(pb_fac);
      mygrid->m_weight[1][iinitial][iobs]->Scale(pb_fac);
    }

    double obssum=0;
    double sum[2][7];
    for(int j=0;j<2;j++) {
      for(int iinitial=0;iinitial<7;iinitial++) {
        sum[j][iinitial]=0;
      }
    }


    for(int itau=0;itau<mygrid->Ntau(iobs);itau++) {
      Q2=mygrid->fQ2(mygrid->m_weight[0][0][iobs]->GetZaxis()->GetBinCenter(itau+1));

      //     alphas=nlojetpdf->alpha_qcd(5,Q2);              // alphas/(2*pi)
      alphas2=alphas*alphas;
      alphas3=alphas2*alphas;

      // fill pdf grid

      for(int iflav=0;iflav<13;iflav++) {
        sprintf(histname,"pdf_%i_%i_%i",iflav, iobs, itau);
        m_pdf[iflav] = new TH1D(histname,histname,mygrid->Ny(iobs),mygrid->ymin(iobs),mygrid->ymax(iobs));
      }

      for(int iy=0;iy<=mygrid->Ny(iobs);iy++) {
        x=mygrid->fx(mygrid->fy(iobs,iy));
        // nlojetpdf->hadronA(x, Q2, 0, 0, f+6);
        for(int iflav=0;iflav<13;iflav++) {
          m_pdf[iflav]->SetBinContent(iy+1,x*f[iflav]);
        }
      }


      // calculate cross section from grid
      for(int iy1=0;iy1<mygrid->Ny(iobs);iy1++) {
        x1=mygrid->fx(mygrid->fy(iobs,iy1));
        for(int iflav=0;iflav<13;iflav++) {
          fA[iflav]=m_pdf[iflav]->GetBinContent(iy1+1)/x1;
        }
        for(int iy2=0;iy2<mygrid->Ny(iobs);iy2++) {
          x2=mygrid->fx(mygrid->fy(iobs,iy2));
          for(int iflav=0;iflav<13;iflav++) {
            fB[iflav]=m_pdf[iflav]->GetBinContent(iy2+1)/x2;
          }

          make_pdf_combinations(fA, fB, H);

          for(int iinitial=0;iinitial<7;iinitial++) {
            born_weight = mygrid->m_weight[0][iinitial][iobs]->GetBinContent(iy1+1,iy2+1,itau+1);
            nlo_weight  = mygrid->m_weight[1][iinitial][iobs]->GetBinContent(iy1+1,iy2+1,itau+1);

            obssum+=(born_weight*alphas2+nlo_weight*alphas3) * H[iinitial];

            sum[0][iinitial]+=born_weight*alphas2*H[iinitial];
            sum[1][iinitial]+=nlo_weight*alphas3*H[iinitial];

            //cout<<"x1="<<x1<<"\tx2="<<x2<<"\tQ2="<<Q2<<"\tobsbin="<<iobs<<endl;
            //cout<<"\tborn_weight="<<born_weight<<"\tnlo_weight="<<nlo_weight<<endl;
            //cout<<"weight="<<weight<<"\tH["<<iinitial<<"]="<<H[iinitial]<<"\talphas2="<<alphas2<<endl;
          }
        }
      }

      for(int iflav=0;iflav<13;iflav++) {             // delete pdf grid for next run
        delete m_pdf[iflav]; m_pdf[iflav]=NULL;
      }
    }


    double grid_value=obssum/mygrid->deltaobs();
    double reference_value=mygrid->m_obs_bins->GetBinContent(iobs+1)*pb_fac/runs;
    // double reference_value=mygrid->m_obs_bins->GetBinContent(iobs+1)*pb_fac;
    //    grid_total->SetBinContent(iobs+1,grid_value);
    //  cout<<"\tgrid="<<grid_value<<"\treference: "<<reference_value<<"\tgrid/reference="<<grid_value/reference_value<<endl;

    for(int j=0;j<2;j++) {
      //cout<<"\tiorder["<<j<<"]:"<<endl;
      for(int i=0;i<7;i++) {
        grid_value=sum[j][i]/mygrid->deltaobs();
        grid[j][i]->SetBinContent(iobs+1,grid_value);

        //reference_value=reference[j][i]->GetBinContent(iobs+1)*pb_fac/runs;
        //cout<<"\t\tgrid="<<grid_value<<"\t reference="<<reference_value<<"\t grid/reference="<<grid_value/reference_value<<endl;
      }
    }

  }

  TFile myfile("standalone.root","RECREATE");
  grid_total->Write();
  for(int j=0;j<2;j++) {
    for(int i=0;i<7;i++) {
       grid[j][i]->Write();
    }
  }
  myfile.Close();

  delete mygrid;

  return 0;
}

