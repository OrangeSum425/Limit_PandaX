#include "iostream"
#include "MakeNiceHisto.C"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1.h"
#include<time.h>

using namespace std;
 
  double md,dE,ratio_p;
  double Mbi,mA,abundance,eff,sigma_p;
// define the parameter for the gama baysian distribution graph;
  double Na,gamaq,gama_px;
  double gama_max,binw;
  int gama_baysian_bin;// the bin number of the baysian distribution of gama - the decay rate;
  gama_max=8E-26;//the boundary of the 0ubb decay rate region that Baysian method takes into acount[/year];
  binw=2E-29;// bin width in the gama baysian distribution;
  gama_baysian_bin= gama_max/binw;
  double gama_p[4000],gamax[4000];//for the PDF drawing;gama_p is the P desity;
// define the detector and background parameters; 
  md=200;//effective  mass [kg];
  Na=6.022E23;//Avogadro constant;
  mA=0.135;//mass per mole [kg/mol];
  abundance=0.9;//abundance of Xe136;
  eff=0.35;//detector efficiency;
  const double ubi=2447.7;//Bi241 Q value;
  const double uq=2458;//Xe136 Q value;
//  const double q_p=6E-27;
//  const double Emin=2470;//lower bound of ROI;
//  const double Emax=2570;//upper bound of ROI;
//  const double bin_number=Emax-Emin;
  const int hfinal_bin=7000000;//range of the half life histogram,*1E24 year;
  double lamda1,lamda2,lamda3;
  double lamda1_po,lamda2_po,lamda3_po;
  double content,content_po;
  double bincontent;
//  double local_value_1[100],local_value_2[100];
//  lamda1=BId*md*td*(Emax-Emin);
  double total_pro_h1,total_pro_h0;
  int r_times=10000;//MC time; 
  int gama_half,median;
  double gama_90,half_life,median_half_life,cutoff_bin,restsum;
  double rate_sum_1,rate_exp_1;
  double rate_sum_2,rate_exp_2;
  double local_rate;
/*
  double td=5;//time of experiment[year];
  double Rbi=0.428;//rate of Bi 241 [cts/year*kg];
  double BId=0.0102;//flat  background rate[cts/year*kev*kg]
  double s_de=2.1;//resolution-sigma[kev];
*/
  TCanvas *c1 = new TCanvas("c1", "graph", 200, 100, 700, 500);

  TH1F *hfinal = new TH1F("hfinal","",hfinal_bin,0,hfinal_bin*1E22);
//  TH1F *h1 = new TH1F("h1","",bin_number,Emin,Emax); 

//  double getmedian(){
  double getmedian(double Rbi,double BId,double td,double s_de, double ROI){
  double Emax = ROI + uq;
  double Emin = uq - ROI;
  const double bin_number=Emax-Emin;
  TH1F *h1 = new TH1F("h1","",bin_number,Emin,Emax);
  double local_value_1[bin_number],local_value_2[bin_number]; 
  lamda1=BId*md*td*(Emax-Emin);
  Mbi=Rbi*md*td;
  lamda2= Mbi;
  
  hfinal->Reset(); 
  TRandom2 r;

  rate_sum_2 = 0;
  rate_sum_1 = 0;
  rate_exp_1 = 1;
  rate_exp_2 = 1;
  for (int i=0;i<bin_number;i++){
      local_value_1[i]  = Mbi*exp(-((i+Emin+0.5-ubi)*(i+Emin+0.5-ubi))/(2*s_de*s_de))/(s_de*2.506)+BId*md*td;  
      local_value_2[i]  = (Na*abundance*eff*md*td/mA)*exp(-((i+Emin+0.5-uq)*(i+Emin+0.5-uq))/(2*(s_de*s_de)))/(sqrt(2*TMath::Pi())*s_de);
//      cout<<"i: "<<i<<"background rate: "<<local_value_1[i]<<endl; 
      rate_sum_1 += local_value_1[i]; 
      rate_sum_2 += local_value_2[i];
  }
//      cout<<"rate_sum_1"<<rate_sum_1<<endl;
  for (int k=0;k<r_times;k++){
    h1->Reset();
    lamda1_po = r.PoissonD(lamda1);
    lamda2_po = r.PoissonD(lamda2);  
//    cout<<lamda1_po<<endl;
//    cout<<lamda2_po<<endl;
//get the energy spectrum simulation
    for (int i=0;i<lamda1_po;i++){
      h1->Fill(r.Uniform(Emin,Emax));
    }
    for(int i =0;i<lamda2_po;i++){
      h1->Fill(r.Gaus(ubi,s_de));
    }
 
//   TCanvas *spectrum = new TCanvas("spectrum", "spectrum", 200, 100, 700, 500); 
//   h1->GetXaxis()->SetTitle("Energy[keV]"); 
//   h1->GetYaxis()->SetTitle("Event Number"); 
//   h1->GetXaxis()->CenterTitle(); 
//   h1->GetYaxis()->CenterTitle(); 
//   SetHistStyle(h1);

//    h1->Draw();

//get the P(gama|H1,E)
  gama_px=0;
  total_pro_h1=0;
  for (int j=0;j<gama_baysian_bin;j++){
    gama_p[j]=1;
    gama_px += binw;
    gamax[j] = gama_px;
//    if (j==199){ 
//      for (int i =0;i<100;i++){
//        local_value = myfunc_h1(i+2408.5);
//        cout<<i<<"local_value"<<local_value<<endl;
//      }
//    }  
    for (int i =0;i<bin_number;i++){
      bincontent = h1->GetBinContent(i+1);
      local_rate = local_value_2[i]*gama_px+local_value_1[i];
      gama_p[j]=gama_p[j]*factor(bincontent,local_rate);
//      gama_p[j]=gama_p[j]*(exp(-1*local_value));
    }
      gama_p[j]=gama_p[j]*exp(-1*(rate_sum_1+gama_px*rate_sum_2));
//    gama_p[j]=gama_p[j]*prior(gama_px);
    restsum = gama_p[j]*(gama_baysian_bin-j)*1000;// the rest smaller than 1/1000 of the tem_sum grantees a 0.01 accuracy of 90% cut;
    if ((restsum < total_pro_h1)&&(gama_px>5E-27)){
    cutoff_bin=j;
//    cout<<"cutoff_bin:"<<cutoff_bin<<endl;
    break;
    }else{
    total_pro_h1 +=gama_p[j];
    }
//   cout<<"gamma_p for:"<<j<<" is : "<<gama_p[j]<<endl;
  }
  
  
/*  for(int j=0;j<gama_baysian_bin;j++){
    total_pro_h1 += gama_p[j];
  }*/
  gama_90=0;
  for(int j=0;j<cutoff_bin;j++){
    gama_p[j]=  gama_p[j]/total_pro_h1 ;
//    cout<<"j:"<<j<<"gama_p:"<<gama_p[j]<<endl;
    gama_90 += gama_p[j];
    if(gama_90 < 0.9){
      gama_half=j; 
    }
  } 
//  TGraph *gr1 = new TGraph (gama_baysian_bin, gamax, gama_p);
//  c1->cd();
//  gr1->Draw();

// h1->Reset();
  cout<<gamax[gama_half]<<endl;
  half_life=0.693147180559945286/gamax[gama_half];
  hfinal->Fill(half_life);

  }

  median = 0;
  median_half_life = 0;
  for(int i=0;i<hfinal_bin;i++){
    if(median<(r_times/2)){
      median += (hfinal->GetBinContent(i+1));
    }else{
      median_half_life = i*1E22;
      break;
    }

  }

  cout<<" median: "<<median_half_life<<endl;

//  TCanvas *halflife_histo = new TCanvas("halflife_histo", "halflife_histo", 200, 100, 700, 500); 
//  halflife_histo->cd();
//  hfinal->Draw();
  return median_half_life;
  hfinal->Reset();
}
/*
double myfunc_h1(double x){
  double result_m = gama_px*(Na*abundance*eff*md*td/mA)*exp(-((x-2458)*(x-2458))/(2*(s_de*s_de)))/(sqrt(2*TMath::Pi())*s_de)+Mbi*exp(-((x-ubi)*(x-ubi))/(2*s_de*s_de))/(s_de*2.506)+BId*md*td;
;
  return result_m;
}
*/
/*
double prior(double x_p){
//  double result_p =exp(-((x_p-q_p)*(x_p-q_p))/(2*sigma_p*sigma_p))/(sigma_p*2.506);
  double result_p=1;
  if((x_p>4E-27)&&(x_p<5E-27)){
    result_p=ratio_p;
  }
  return result_p;
}
*/

double factor(double bincontent,double local){
  double result=1;
if(bincontent>0){
  while(bincontent>0){
    result = result*local/bincontent;
    bincontent--;
  }
}
  return result;
}

  


