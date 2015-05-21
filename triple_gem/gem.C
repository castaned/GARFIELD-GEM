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


  // Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-10.*Pitch, -1., -10.*Pitch, 10.*Pitch, 1., 10.*Pitch);

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

 const bool plotDrift = true;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) 
  {
    driftView->SetArea(-2 * pitch, -2 * pitch, -0.02,
		       2 * pitch,  2 * pitch,  0.02);

    // Plot every 10 collisions (in microscopic tracking).
    aval->SetCollisionSteps(10); 
    aval->EnablePlotting(driftView);
    drift->EnablePlotting(driftView);
  }

 
 
  if (plotDrift) {
    driftView->SetCanvas(cD);
    driftView->Plot();
   }


  app.Run(kTRUE);

}
