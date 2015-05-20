using namespace std;
#include <iostream>
#include <fstream>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
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
#include <TTree.h>
#include <TFile.h>
#include <TNtuple.h>

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"

using namespace Garfield;
using namespace std;

char Name[100]; 
ifstream in;
ofstream out;

int main(int argc, char * argv[]) 
{

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  const bool debug = true;

  // Load the field map.
  ComponentAnsys123* fm = new ComponentAnsys123();

  const std::string efile = "/panfs/vol/HEP/GEM/ansys_input/v3850/pitch560/ELIST.lis";
  const std::string nfile = "/panfs/vol/HEP/GEM/ansys_input/v3850/pitch560/NLIST.lis";
  const std::string mfile = "/panfs/vol/HEP/GEM/ansys_input/v3850/pitch560/MPLIST.lis";
  const std::string sfile = "/panfs/vol/HEP/GEM/ansys_input/v3850/pitch560/PRNSOL.lis";

  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  fm->PrintRange();
 
  // Dimensions of a single GEM
  const double kapton = 50.e-4;
  const double metal  = 5.e-4;
  const double outdia = 70.e-4;
  const double middia = 50.e-4;

  // Dimensions of the detector
  const double  singlecell  = 0.014;   // Distance between 2 GEM holes, in cm
  const double  nholes      = 4;   // Number of holes per strip pitch
  //const double  nholes      = 5;   
  //const double  nholes      = 6;   
  //const double  nholes      = 7;   
  //const double  nholes      = 8;   
  
  const double  pitch       = singlecell * nholes;   // strip pitch in cm

  const double  driftregion = 0.3 ;   // thickness of the drift region
  const double  transfer1   = 0.1 ;   // thickness of the transfer1 gap
  const double  transfer2   = 0.2 ;   // thickness of the transfer2 gap
  const double  induct      = 0.1 ;   // thickness of the induction gap
  double        total       = transfer1+transfer2+induct;  // thickness of the total Z regoin

  const double  stripwidth  = 0.036;  // width of the stripwidth in cm
  const double  thickness   = 0.0035; // thickness of the anode/stripwidth in cm
  
  const bool plotField = false;
  if (plotField) 
  {
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(fm);
    fieldView->PlotProfile(0., 0., 0.08, 0., 0., -0.08);
    fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
    fieldView->SetArea(-pitch/2 , -0.3, pitch/2 , 0.3);
    fieldView->SetVoltageRange(-160., 160.);
    TCanvas* cF = new TCanvas();
    fieldView->SetCanvas(cF);
    fieldView->PlotContour();
  }
  
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 45., "co2", 15., "cf4", 40.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->EnableDebugging();
  gas->Initialise();  
  gas->DisableDebugging();
  
  // Set the Penning transfer efficiency.
  const double rPenning = 0.1;
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

  // Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-5 * pitch, -5 * pitch, -5, 5 * pitch,  5 * pitch,  5);

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  
  const bool plotDrift = false;
  ViewDrift* driftView = new ViewDrift();
  
  if (plotDrift) 
  {
    driftView->SetArea(-3 * pitch, -3 * pitch, 4,  3 * pitch,  3 * pitch,  3);
    
    // Plot every 10 collisions (in microscopic tracking).
    aval->SetCollisionSteps(10); 
    aval->EnablePlotting(driftView);
    //drift->EnablePlotting(driftView);
  }

  // Histograms
  int nBinsGain = 100000;
  double gmin   = 0.;
  double gmax   = 100000.;
  
  // open a root file to write results
  TFile* f = new TFile("v3850_pitch560.root","RECREATE");

  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons", nBinsGain, gmin, gmax);
  TH1F* hIons      = new TH1F("hIons", "Number of ions          ", nBinsGain, gmin, gmax);
  TNtuple* ntuple  = new TNtuple("ntuple","ntuple","ne:xe1:ye1:ze1:xe2:ye2:ze2:e2");

  double sumElectronTotal     = 0.;
  double sumElectron1Strip    = 0.;
  double sumElectron2Strip    = 0.;
  double sumElectron3Strip    = 0.;
  double sumElectron4Strip    = 0.;
  double sumElectronHalfStrip = 0.;
  
  // open a data file to write results
  out.open("v3850_pitch560.dat");
 
  const int nEvents =1000;

  for (int i = nEvents; i--;) 
  {     
    //debug - nevents
    if (debug || i % 10 == 0) std::cout << i << "/" << nEvents << "\n";
    
    // Randomize the initial position.
    const double smear = pitch; 
    //double x0 = -smear + 2*RndmUniform() * smear;
    //double y0 = -smear + 2*RndmUniform() * smear;
    //double x0 = pitch/2.;
    //double y0 = sqrt(3)*singlecell/2;
  
    double x0 = RndmUniform()*smear;
    double y0 = RndmUniform()*sqrt(3)*singlecell;
    double z0 = 0.01; 
    double t0 = 0.;
    double e0 = 0.;

    aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    
    int ne = 0, ni = 0;
    aval->GetAvalancheSize(ne, ni);
    hElectrons->Fill(ne);
    hIons->Fill(ni); 
    
    const int np = aval->GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    
    if (ne>50)
    {
      
      for (int j = np; j--;) 
      {
	//electrons
	aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
	
	ntuple->Fill(ne,xe1,ye1,ze1,xe2,ye2,ze2,e2);
	
	sumElectronTotal += 1.;
	
	//usual number of electrons collected in the first half strips
	if (ze2 >= -total && ze2 <= -total + thickness 
	  && ( (xe2 >= 0. && xe2 <=  stripwidth/2.) || (xe2 >= pitch-stripwidth/2. && xe2 <= pitch) ) )
	  
	sumElectronHalfStrip += 1.;
	
	//number of electrons collected in different strips  
	 if (ze2 >= -total && ze2 <= -total + thickness 
	    && ( (xe2 >= -stripwidth/2. && xe2 <=  stripwidth/2.) || (xe2 >= pitch-stripwidth/2. && xe2 <= pitch+stripwidth/2.) ) )
	    
	 sumElectron1Strip += 1.;
	 
	 if (ze2 >= -total && ze2 <= -total + thickness 
	   && ( (xe2 >= -stripwidth/2.       && xe2 <= stripwidth/2.)         || (xe2 >= pitch-stripwidth/2.   && xe2 <= pitch+stripwidth/2.  ) 
	   ||   (xe2 >= -pitch-stripwidth/2. && xe2 <= -pitch+stripwidth/2.)  || (xe2 >= 2*pitch-stripwidth/2. && xe2 <= 2*pitch+stripwidth/2.) ) )
	   
	 sumElectron2Strip += 1.;
	 
	 if (ze2 >= -total && ze2 <= -total + thickness 
	   && ( (xe2 >= -stripwidth/2.         && xe2 <= stripwidth/2.)          || (xe2 >= pitch-stripwidth/2.   && xe2 <= pitch+stripwidth/2.   ) 
	   ||   (xe2 >= -pitch-stripwidth/2.   && xe2 <= -pitch+stripwidth/2.)   || (xe2 >= 2*pitch-stripwidth/2. && xe2 <= 2*pitch+stripwidth/2. ) 
	   ||   (xe2 >= -2*pitch-stripwidth/2. && xe2 <= -2*pitch+stripwidth/2.) || (xe2 >= 3*pitch-stripwidth/2. && xe2 <= 3*pitch+stripwidth/2. ) ) )
	   
	 sumElectron3Strip += 1.;
	 
	 if (ze2 >= -total && ze2 <= -total + thickness 
	   && ( (xe2 >= -stripwidth/2.         && xe2 <= stripwidth/2.)          || (xe2 >= pitch-stripwidth/2.   && xe2 <= pitch+stripwidth/2.   ) 
	   ||   (xe2 >= -pitch-stripwidth/2.   && xe2 <= -pitch+stripwidth/2.)   || (xe2 >= 2*pitch-stripwidth/2. && xe2 <= 2*pitch+stripwidth/2  ) 
	   ||   (xe2 >= -2*pitch-stripwidth/2. && xe2 <= -2*pitch+stripwidth/2.) || (xe2 >= 3*pitch-stripwidth/2. && xe2 <= 3*pitch+stripwidth/2. )  
	   ||   (xe2 >= -3*pitch-stripwidth/2. && xe2 <= -3*pitch+stripwidth/2.) || (xe2 >= 4*pitch-stripwidth/2. && xe2 <= 4*pitch+stripwidth/2. ) ) )
	   
	 sumElectron4Strip += 1.;
	 
      }
    } 
    
    const double AllElectron  = sumElectronTotal;
    const double SumHalfStrip = sumElectronHalfStrip;
    const double Sum1Strip    = sumElectron1Strip;
    const double Sum2Strip    = sumElectron2Strip;
    const double Sum3Strip    = sumElectron3Strip;
    const double Sum4Strip    = sumElectron4Strip;
    
    out<<ne<<"  "<<AllElectron<<"  "<<SumHalfStrip<<"  "<<Sum1Strip<<"  "<<Sum2Strip<<"  "<<Sum3Strip<<"  "<<Sum4Strip<<endl;
    //cout<<ne<<"  "<<AllElectron<<"  "<<SumHalfStrip<<"  "<<Sum1Strip<<"  "<<Sum2Strip<<"  "<<Sum3Strip<<"  "<<Sum4Strip<<endl;    
  }
  // close data and root files
  out.close();
  f->Write();

}
