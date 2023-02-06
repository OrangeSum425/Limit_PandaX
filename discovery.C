#include "iostream"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1.h"

using namespace std;

  double BId,md,td,dE,Rbi;
  double Mbi,s_de,mA,abundance,eff;
  double Na,gamaq,gama_px,Mq;
  double gama_max,binw,ergodic_binw;
  int h1_bin;
  gama_max=1E-25;
  binw=2E-28;
  ergodic_binw=1E-27;
  h1_bin= gama_max/binw;
  double gama_p[500];
  BId=0.01;
  md=741;
  s_de=2.1;
  td=5;
  Na=6.022E23;
  mA=0.1596;
  abundance=0.341;
  eff=0.813;
  Rbi=0;
//  gamaq=1E-25;
  Mq=gamaq*Na*abundance*eff*md*td/mA;
  Mbi=Rbi*md*td;
  const double ubi=2447.7;
  const double uq=2527;
  const double Emin=2521;
  const double Emax=2533;
  const int hfinal_bin=2000;
  double lamda1,lamda2,lamda3;
  double lamda1_po,lamda2_po,lamda3_po;
  double content,content_po;
  double bincontent,local_value;
  lamda1=BId*md*td*(Emax-Emin);
  lamda2= Mbi;
  lamda3= Mq;
  double total_pro_h1,total_pro_h0;
  int r_times=100;
  double PHE;
  int number_3sigma;
  double p_3sigma=0.0027;
  double n_3sigma = p_3sigma*hfinal_bin;

//  TCanvas *c1 = new TCanvas("c1", "graph", 200, 100, 700, 500);  
  TH1F *hfinal = new TH1F("hfinal","",hfinal_bin,0,1);
  TH1F *h1 = new TH1F("h1","energy spectrum simulation",100,Emin,Emax);
  TH1F *h1pre = new TH1F("h1pre","energy spectrum simulation",100,Emin,Emax);

void discovery(){
  int N_2;
    gamaq=8E-26;
  for (int j=0;j<h1_bin;j++){
    gamaq += ergodic_binw;
    Mq=gamaq*Na*abundance*eff*md*td/mA;
    lamda3= Mq;
    N_2 = 2*getnum();
//something wrong,,when P(H0|E)<0.0027 for N/2 experiments. 
// 4% error for N/2 
    if( N_2>r_times){
      cout<<"for t="<<td<<" T is "<<0.693147/gamaq<<endl;
      break;
    }
  }
}

int getnum(){

  TRandom2 r;

  for (int k=0;k<r_times;k++){
    lamda1_po = r.PoissonD(lamda1);
    lamda2_po = r.PoissonD(lamda2);
    lamda3_po = r.PoissonD(lamda3);

//get the energy spectrum simulation
    for (int i=0;i<lamda1_po;i++){
      h1->Fill(r.Uniform(Emin,Emax));
    }
    for(int i =0;i<lamda2_po;i++){
      h1->Fill(r.Gaus(ubi,s_de));
    }
    for(int i =0;i<lamda3_po;i++){
      h1->Fill(r.Gaus(uq,s_de));
    }
//an alternative way
/*
  for (int i=0;i<lamda1;i++){
    h1pre->Fill(r.Uniform(Emin,Emax));
  }
  for(int i =0;i<lamda2;i++){
    h1pre->Fill(r.Gaus(ubi,s_de));
  }
  for(int i =0;i<lamda3;i++){
    h1pre->Fill(r.Gaus(uq,s_de));
  }

  for(int i =0;i<100;i++){
    content = h1pre->GetBinContent(i+1);
    content_po= r.PoissonD(content);
    h1->Fill(2470+i,content_po);
  }
*/


//get the P(E|H1)
  gama_px=0;
  for (int j=0;j<h1_bin;j++){
    gama_p[j]=1;
    gama_px += binw;
    for (int i =0;i<12;i++){
      bincontent = h1->GetBinContent(i+1);
      local_value = myfunc_h1(i+2521.5);
      gama_p[j]=gama_p[j]*(exp(-1*local_value)*factor(bincontent,local_value));

    }
  }
  total_pro_h1=0;
  for(int j=0;j<h1_bin;j++){
    total_pro_h1 += gama_p[j];
  } 
  total_pro_h1= total_pro_h1*binw/gama_max;
//  cout<<total_pro_h1<<" "<<h1_bin<<endl;
  
//get the P(E|H0)
    total_pro_h0=1;
  for (int i =0;i<12;i++){
    bincontent = h1->GetBinContent(i+1);
    local_value = myfunc_h0(i+2521.5);
    total_pro_h0=total_pro_h0*(exp(-1*local_value)*factor(bincontent,local_value));
  }
//  cout<<total_pro_h0<<endl;

  h1->Reset();
  h1pre->Reset();
  PHE=total_pro_h0/(total_pro_h0+total_pro_h1);
//  cout<<PHE<<endl;
  hfinal->Fill(PHE);

  }

//get the number of P(h0|e)<0.0027 after r_times
//if more than N/2 then we find the gamaq, if not try another value.
  number_3sigma = 0;
  for(int i=0;i<n_3sigma;i++){
    number_3sigma += (hfinal->GetBinContent(i+1));
  }
  cout<<"N "<<r_times<<" "<<number_3sigma<<endl;
//  hfinal->Draw();
  hfinal->Reset();
  return number_3sigma;
}

double myfunc_h0(double x){
  double result_m = Mbi*exp(-((x-ubi)*(x-ubi))/(2*(s_de*s_de)))/(s_de*2.506)+BId*md*td;
;
  return result_m;
}

double myfunc_h1(double x){
  double result_m = gama_px*(Na*abundance*eff*md*td/mA)*exp(-((x-uq)*(x-uq))/(2*(s_de*s_de)))/(sqrt(2*TMath::Pi())*s_de)+Mbi*exp(-((x-ubi)*(x-ubi))/(2*(s_de*s_de)))/(s_de*2.506)+BId*md*td;
;
  return result_m;
}

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



