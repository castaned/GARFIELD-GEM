#include <iostream>
#include <fstream>
using namespace std;
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TText.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TCutG.h>
#include <TTree.h>

void analysis(char* gas, char* type)
{

  //========= Analysis of GARFIELD simulation output ==================================//
  // Macro to produce the plot  Gain vs Voltage for different gas mix combination and penning transfer
  // The macro works in the following way:
  // 1.- The directory with the GARFIELD root files OUTPUT_ROOT/   has to be located in the same directory this macro will be executed.
  // 2.- Change nSim to the number of voltage simulated points (i.e.  Int_t nsim=3;)
  // 3.- Below in the code change the vector with voltage point to simulated number (i.e.  Float voltage[n] = {v1,v2,....,vn})
  // 4.- Save and close the macro.
  // 5.- To run the macro from terminal execute :  root -l 'analysis.C("ar","voltage")'  (for arc02 and gain vs voltage)
  // 6.- The macro produces two TGraphs for effective and total gain
  //===================================================================================//

  TFile  *f[30]; 
  TTree  *t[30];
  TH1F   *h[30];

  Int_t nSim=4;  // Number of voltage points;
  
  Float_t voltage[4]={2679,2867,3055,3243};  // Define the values for voltage
  
  Float_t tgain[20];
  Float_t egain[20];
  
  for(int k=0;k<nSim;k++){
    
    
    if(type=="voltage"){
      char nameinput[100];
      double v_temp = voltage[k];
      if(gas=="ar") sprintf(nameinput,"OUTPUT_ROOT/output_%sco2_v%4d_pen%0.1f.root",gas,v_temp,0.4);
      if(gas=="ne") sprintf(nameinput,"OUTPUT_ROOT/output_%sco2_v%4d_pen%0.1f.root",gas,v_temp,0.48);
      f[k]  = new TFile(nameinput);
    }
    
    t[k] = (TTree*)(f[k]->Get("ntuple")); 
    h[k] = (TH1F*)(f[k]->Get("hElectrons"));
    
    // Declaration of leaves types 
    Float_t xe2,ye2,ze2; 
    
    t[k]->SetBranchAddress("xe2",&xe2);
    t[k]->SetBranchAddress("ye2",&ye2);
    t[k]->SetBranchAddress("ze2",&ze2);
    
    // determine number of events in the tree
    Int_t nentries = t[k]->GetEntries();
    
    int s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    
    for (Int_t j = 0; j < nentries;j++) 
      {
	t[k]->GetEntry(j);
	if (ye2 <= -0.4 && xe2 >= -0.122 && xe2 <= -0.054) s1 += 1;
	if (ye2 <= -0.4 && xe2 >= -0.034 && xe2 <=  0.034) s2 += 1;
	if (ye2 <= -0.4 && xe2 >=  0.054 && xe2 <=  0.102) s3 += 1;
	if (ye2 <= -0.4 && xe2 >=  0.122 && xe2 <=  0.190) s4 += 1;
      }
    
    //    cout<<" N = "<<h[k]->GetEntries()<<" , s1 = "<<s1<<", s2 = "<<s2<<", s3 = "<<s3<<", s4 = "<<s4<<endl;
    
    tgain[k] = nentries / h[k]->GetEntries(); 
    egain[k] = (s1 + s2 + s3 + s4) / h[k]->GetEntries(); 
  }
  
  TGraphErrors *gr_tgain = new TGraphErrors(nSim, voltage, tgain, 0, 0);
  TGraphErrors *gr_egain = new TGraphErrors(nSim, voltage, egain, 0, 0);
  
  TLegend *leg = new TLegend(0.22,0.72,0.56,0.90);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.043);

  char nameleg[100];
  char nameleg2[100];
  if(gas=="ar") sprintf(nameleg,"Gas mixture ArCO2 [70/30%]");
  if(gas=="ne") sprintf(nameleg,"Gas mixture NeCO2 [70/30%]");
  leg->AddEntry(gr_tgain,nameleg,"PL"); 
  if(gas=="ar") sprintf(nameleg2,"Penning Transfer=0.4");
  if(gas=="ne") sprintf(nameleg2,"Penning Transfer=0.48");
  leg->AddEntry((TObject*)0, nameleg2, "");


  gr_tgain->SetMarkerSize(1.2); 
  gr_egain->SetMarkerSize(1.2); 

  char namecanvas[100];
  if(gas=="ar") sprintf(namecanvas,"effgain_vs_voltage_%sco2_pen%0.1f.pdf",gas,0.4);
  if(gas=="ne") sprintf(namecanvas,"effgain_vs_voltage_%sco2_pen%0.1f.pdf",gas,0.48);

  char namecanvas2[100];
  if(gas=="ar") sprintf(namecanvas2,"totalgain_vs_voltage_%sco2_pen%0.1f.pdf",gas,0.4);
  if(gas=="ne") sprintf(namecanvas2,"totalgain_vs_voltage_%sco2_pen%0.1f.pdf",gas,0.48);


  TCanvas *c = new TCanvas("c","c",700,500);
  c->SetLogy(); 
  c->SetGridx();
  c->SetGridy();
  if(type="voltage") gr_egain->GetXaxis()->SetTitle("Drift Voltage [V]"); 
  gr_egain->GetYaxis()->SetTitle("Effective Gain");
  gr_egain->SetMinimum(1.0);
  gr_egain->SetMaximum(10000.0);
  gr_egain->Draw("APLsame");
  leg->Draw("same");
  c->SaveAs(namecanvas,"recreate");
  
  TCanvas *c1 = new TCanvas("c1","c1",700,500); 
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  if(type="voltage") gr_tgain->GetXaxis()->SetTitle("Drift Voltage [V]");
  gr_tgain->GetYaxis()->SetTitle("Total Gain");
  gr_tgain->SetMinimum(1.0);
  gr_tgain->SetMaximum(10000.0);
  gr_tgain->Draw("APLsame");
  leg->Draw("same");
  c1->SaveAs(namecanvas2,"recreate");

}
