

void check() {

 TFile *f1 = new TFile("out.root");



 TCanvas *c1= new TCanvas("c1");
 ratio->Draw();
 
//  TCanvas *c2= new TCanvas("corder");
//  c2->Divide(1,2);
//  order->cd();
//  c2->cd(1);
//  ratio_lo->Draw();
//  c2->cd(2);
//  ratio_nlo->Draw();

 TCanvas *c3= new TCanvas("csubprocess");
 c3->Divide(2,3);

 f1->cd();
 subProcess->cd();

 for (Int_t i=0; i<7; i++) {
  c3->cd(i+1);
  char hname[100];
  sprintf(hname,"subProcess/ratio_sub_%d",i);
  cout<<" hname= "<<hname<<endl;
  TH1D* h1=(TH1D*) f1->Get(hname);
  h1->Draw();
 }

 TCanvas *c4= new TCanvas("cscale");
 c4->Divide(2,2);
 f1->cd();
 scale->cd();
 for (Int_t i=0; i<3; i++) {
  c4->cd(i+1);
  char hname[100];
  sprintf(hname,"scale/ratio_scale_%d",i);
  cout<<" hname= "<<hname<<endl;
  TH1D* h2=(TH1D*) f1->Get(hname);
  h2->Draw();
 }


 f1->cd();
 subScale->cd();
 for (Int_t ip=0; ip<7; ip++) {
  char cname[100]; sprintf(cname,"cname%d",ip);
  TCanvas *c5= new TCanvas(cname,cname);
  c5->Divide(2,2);
  for (Int_t is=0; is<3; is++) {
   c5->cd(is+1);
   char hname[100];
   sprintf(hname,"subScale/ratio_subscale_%d_%d",ip,is);
   cout<<" hname= "<<hname<<endl;
   TH1D* h3=(TH1D*) f1->Get(hname);
   h3->Draw();
  }
 }

}
