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

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ComponentAnalyticField.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  const bool debug = true;

  // Load the field map.
  ComponentAnsys123* fm = new ComponentAnsys123();
  const std::string efile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2679/ELIST.lis";
  const std::string nfile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2679/NLIST.lis";
  const std::string mfile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2679/MPLIST.lis";
  const std::string sfile = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gain/v2679/PRNSOL.lis";
  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  fm->PrintRange();

  // Dimensions of the detector
  const double Pitch = 0.088; // Pitch of strips
  const double pitch = 0.014; // pitch of holes
  const double smear = pitch/2.;
  double singlecell  = 0.014; // Dimensions of the GEM

  // Magnetic field 
  const double MagX = 0.;
  const double MagY = 0.;
  const double MagZ = 0.;


  const bool plotField = false;
  if (plotField)
    {
      ViewField* fieldView = new ViewField();
      fieldView->SetComponent(fm);
      fieldView->PlotProfile(0., 0., 0.02, 0., 0., -0.02);
      fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
      fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
     fieldView->SetVoltageRange(-360., 360.);
     TCanvas* cF = new TCanvas();
     fieldView->SetCanvas(cF);
     fieldView->PlotContour();
    }
  
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 70., "co2", 30.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->EnableDebugging();
  gas->Initialise();  
  gas->DisableDebugging();

  // Set the Penning transfer efficiency.
  const double rPenning = 0.57;
  const double lambdaPenning = 0.;
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  gas->LoadIonMobility("IonMobility_Ar+_Ar.txt");
  
  // Associate the gas with the corresponding field map material. 
  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }
  fm->PrintMaterials();

  // Make a component with analytic electric field.
  ComponentAnalyticField* cmpAmp  = new ComponentAnalyticField();
  
  cmpAmp->AddPlaneY(0.30275, 1., "Driftplane");
  cmpAmp->AddPlaneY(-0.4    , 0., "striplane");
  
  //Next we construct the Strips for readout of te signal, also with labels
  double  Xstrip1 = -0.088, Xstrip2 = 0.0, Xstrip3 = 0.088, Xstrip4 = 0.176 ; //Store the center of strips
  cmpAmp->AddStripOnPlaneY('z', -0.4, -0.122 , -0.054, "Strip1");
  cmpAmp->AddStripOnPlaneY('z', -0.4, -0.034,   0.034, "Strip2");
  cmpAmp->AddStripOnPlaneY('z', -0.4,  0.054,   0.122, "Strip3");
  cmpAmp->AddStripOnPlaneY('z', -0.4,  0.142,   0.21, "Strip4");
 
  //calculate signal induced on the strip using ComponentAnalyticalField
  cmpAmp->AddReadout("Strip1");
  cmpAmp->AddReadout("Strip2");
  cmpAmp->AddReadout("Strip3");
  cmpAmp->AddReadout("Strip4");

  //Set constant magnetic field in [Tesla]
  fm->SetMagneticField(MagX, MagY, MagZ);
  cmpAmp->SetMagneticField(MagX, MagY, MagZ);

  // Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-10.*Pitch, -1., -10.*Pitch, 10.*Pitch, 1., 10.*Pitch);


  sensor->AddElectrode(cmpAmp, "Strip1"); 
  sensor->AddElectrode(cmpAmp, "Strip2"); 
  sensor->AddElectrode(cmpAmp, "Strip3"); 
  sensor->AddElectrode(cmpAmp, "Strip4"); 
  
  const double tStart = 0.;
  const double tStop  = 1000.;
  const int nSteps    = 1000;
  const double tStep  = (tStop - tStart) / nSteps;
  
  sensor->SetTimeWindow(tStart, tStep, nSteps);

  
  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  aval->EnableSignalCalculation(); 
  aval->SetTimeWindow(tStart,tStop);
  
  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

  sensor->ClearSignal();

  const bool plotDrift = true;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) 
    {
      driftView->SetArea(-5 * pitch, -1., -5 * pitch,
			 5 * pitch,  -1,  5 * pitch);
      
      // Plot every 10 collisions (in microscopic tracking).
      aval->SetCollisionSteps(10); 
      aval->EnablePlotting(driftView);
      drift->EnablePlotting(driftView);
    }
  
 // Histograms
  int nBinsGain = 100;
  double gmin =   0.;
  double gmax = 100.;
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons",
                              nBinsGain, gmin, gmax);

  const int nEvents = 1;
  for (int i = nEvents; i--;) { 
    if (debug || i % 10 == 0) std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position.
    const double smear = Pitch / 2.; 
    double x0 = -smear + RndmUniform() * smear;
    double y0 = 0.25;
    double z0 = -smear + RndmUniform() * smear;
    double t0 = 0.;
    double e0 = 0.5;
    aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval->GetAvalancheSize(ne, ni);
    hElectrons->Fill(ne);
  }


  // TCanvas* cD = new TCanvas();
 
  // if (plotDrift) {
  //   driftView->SetCanvas(cD);
  //   driftView->Plot();
  //  }

  const bool plotHistogram = true;
  if (plotHistogram) {
    TCanvas* cH = new TCanvas("cH", "Histograms", 800, 700);
    cH->Divide(2, 2);
    cH->cd(1);
    hElectrons->Draw();
  }
  
  app.Run(kTRUE);

}
