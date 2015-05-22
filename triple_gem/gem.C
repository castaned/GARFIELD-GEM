using namespace std;
#include <iostream>
#include <fstream>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TNtuple.h>
#include "ComponentAnsys123.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ViewField.hh"
#include "ViewFEMesh.hh"
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) 
  {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
    
  // Dimensions of the GEM
  const double pitch = 0.014; // pitch of holes
  const double smear = pitch/2.;

  // Load the field map.
  ComponentAnsys123* fm = new ComponentAnsys123();

  const std::string efile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2867/ELIST.lis";
  const std::string nfile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2867/NLIST.lis";
  const std::string mfile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2867/MPLIST.lis";
  const std::string sfile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2867/PRNSOL.lis";

  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityZ();
  fm->SetMagneticField(0.,0.,0.);
  fm->PrintRange();
  const bool plotField = false;
  if (plotField) {
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(fm);
    fieldView->SetPlane(0., 0., -1., -0.007, 0., 0.);
    fieldView->SetArea(-.03, -0.4, 0.019, 0.5);
    fieldView->SetVoltageRange(-3400., 0.);
    TCanvas* cF = new TCanvas();
    fieldView->SetCanvas(cF);
    fieldView->PlotContour();
  }
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  //  gas->SetComposition("ar", 45., "co2", 15., "cf4", 40.);
gas->SetComposition("ar", 70., "co2", 30.);
  // Temperature and pressure
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->Initialise();  

  // Set the Penning transfer efficiency.
  const double rPenning = 0.57;
  const double lambdaPenning = 0.;
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  
  // Load the ion mobilities.
  gas->LoadIonMobility("/panfs/vol/y/yamaghr76/garfield/Data/IonMobility_Ar+_Ar.txt");
  
  // Associate the gas with the corresponding field map material. 
  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) 
  {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }
  fm->PrintMaterials();

  //Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-2*pitch, -0.5, -.048, 2*pitch, 0.3, 0.048);
  //  sensor->SetArea(-10.*pitch, -10.*pitch, -10., 10.*pitch, 10.*pitch, 10.);
  
  // avalanche 
  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  //aval->EnableAvalancheSizeLimit(1000);
  
  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);
  
  // //Canvas to display the drift lines and the field
  // TCanvas * c = new TCanvas("c", "c", 1000, 20, 700, 500);
  // ViewDrift* viewdrift = new ViewDrift();
  // //  viewdrift->SetArea(-0.01, -0.4, 0., 0.,0.3, 0.02);
  // viewdrift->SetArea(-2*pitch, -0.4, -.048, 2*pitch, 0.3, 0.048);
  // viewdrift->SetClusterMarkerSize(0.08);
  // viewdrift->SetCollisionMarkerSize(0.1);

  // // viewdrift->SetCanvas(c);
  // // new
  // ViewFEMesh* meshView = new ViewFEMesh();
  // meshView->SetComponent(fm);
  // meshView->SetArea(-2*pitch, -0.4, -.048, 2*pitch, 0.3, 0.04);
  // meshView->SetFillMesh(true);
  // //  meshView->SetCanvas(c);
  
  // aval->EnablePlotting(viewdrift);
  // drift->EnablePlotting(viewdrift);
  // // viewdrift->Plot();
  
  // open a root file to write results
  TFile* f = new TFile("test.root","RECREATE");
  TNtuple* ntuple  = new TNtuple("ntuple","","ne:xe1:ye1:ze1:te1:e1:xe2:ye2:ze2:te2:e2");
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons", 1000, 0, 1000);
  
  // intitial parameters of the electron
  const int nEvents = 2;
  for (int i = 0; i<nEvents; i++) 
  {    
    double x0 = -0.0025;
    double y0 = 0.025; 
    double z0 = 0.006;
    double t0 = 0.;  // time
    double e0 =1.0; //energy   0.5keV
       
    //make the avalanche for the defined electron
    aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    
    int ne = 0, ni = 0;
    aval->GetAvalancheSize(ne, ni);
    
    //fill histogram of the avalanche size
    hElectrons->Fill(ne);
    
    //track all the produced electrons 
    const int np = aval->GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    int status;
    
    for (int j = np; j--;) 
    {
	 //electrons end points
	 aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
	 
	 //fill the ntuple
	 ntuple->Fill(ne,xe1,ye1,ze1,te1,e1,xe2,ye2,ze2,te2,e2);		
    }
    cout<<"event : "<<i<<", total numebr of produced electrons = "<<ne<<endl;
    //    viewdrift->Plot();

    //    meshView->SetViewDrift(viewdrift);
    //    meshView->Plot();
    
  }
  //  viewdrift->Plot();
  //write the results in to the ntuple
  f->Write();
   
 
  app.Run(kTRUE);

  //end
}
