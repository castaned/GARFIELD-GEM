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
#include "ComponentAnalyticField.hh"
using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) 
{

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
    
  // Dimensions of the GEM
  const double pitch = 0.014; // pitch of holes
  const double pitch_strip = 0.088; // pitch of strips
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


  // Make a component with analytic electric field.
  ComponentAnalyticField* cmpAmp  = new ComponentAnalyticField();

  cmpAmp->AddPlaneY(0.30275, 1., "Driftplane");
  cmpAmp->AddPlaneY(-0.4    , 0., "striplane");
  
  //Next we construct the Strips for readout of te signal, also with labels
  //  double  Xstrip1 = -0.088, Xstrip2 = 0.0, Xstrip3 = 0.088, Xstrip4 = 0.156 ; //Store the center of strips
  double  Xstrip1 = -0.088, Xstrip2 = 0.0, Xstrip3 = 0.088, Xstrip4 = 0.176 ; //Store the center of strips
  cmpAmp->AddStripOnPlaneY('z',  -0.4, -0.122, -0.054, "Strip1");
  cmpAmp->AddStripOnPlaneY('z',  -0.4, -0.034,  0.034, "Strip2");
  cmpAmp->AddStripOnPlaneY('z',  -0.4,  0.054,  0.122, "Strip3");
  cmpAmp->AddStripOnPlaneY('z',  -0.4,  0.142,  0.21,  "Strip4");
 
  //calculate signal induced on the strip using ComponentAnalyticalField
  cmpAmp->AddReadout("Strip1");
  cmpAmp->AddReadout("Strip2");
  cmpAmp->AddReadout("Strip3");
  cmpAmp->AddReadout("Strip4");
  
  //Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-2*pitch, -0.5, -.048, 2*pitch, 0.3, 0.048);
  //  sensor->SetArea(-10.*pitch, -10.*pitch, -10., 10.*pitch, 10.*pitch, 10.);

  sensor->AddElectrode(cmpAmp, "Strip1"); 
  sensor->AddElectrode(cmpAmp, "Strip2"); 
  sensor->AddElectrode(cmpAmp, "Strip3"); 
  sensor->AddElectrode(cmpAmp, "Strip4");

  //  const double tStart = 0.;
  // const double tStop  = 1000.;
  // const int nSteps    = 1000;

  //  const double tStop  = 1000.;
  //  const int nSteps    = 1000;
  //  const double tStep  = (tStop - tStart) / nSteps;
  
  //  sensor->SetTimeWindow(tStart, tStep, nSteps);

  
  // avalanche 
  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  aval->EnableSignalCalculation(); 
  //  aval->SetTimeWindow(tStart,tStop);
  //aval->EnableAvalancheSizeLimit(1000);
  
  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

  sensor->ClearSignal();
  
  const bool plotDrift = true;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) 
    {
      driftView->SetArea(-2 * pitch_strip, -0.4, -2 * pitch,
			 2 * pitch_strip,  0.3, 2*pitch);
      
      // Plot every 10 collisions (in microscopic tracking).
      aval->SetCollisionSteps(10); 
      aval->EnablePlotting(driftView);
      drift->EnablePlotting(driftView);
    }



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
  const int nEvents = 1;
  for (int i = 0; i<nEvents; i++) 
    {    
      double x0 = 0.0;
      double y0 = 0.25; 
      double z0 = 0.0;
      double t0 = 0.;  // time
      double e0 = 0.1; //energy   
      
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
      
	  cout<<" number of endpoint electron "<<j<<"  with status  "<<status<<endl;
      
      
	  //fill the ntuple
	  ntuple->Fill(ne,xe1,ye1,ze1,te1,e1,xe2,ye2,ze2,te2,e2);		
	}
      cout<<"event : "<<i<<", total numebr of produced electrons = "<<ne<<endl;
      //    viewdrift->Plot();

      //    meshView->SetViewDrift(viewdrift);
      //    meshView->Plot();
    
    }

  TCanvas* cD = new TCanvas();
  
  //  bool plot_geometry=false;
  bool plot_geometry=true;

  if(plot_geometry){

 // Dimensions of the GEM
  const double pitch = 0.014;
  const double kapton = 50.e-4;
  const double metal = 5.e-4;
  const double outdia = 70.e-4;
  const double middia = 50.e-4;
  


 // Build the geometry in Root.
  TGeoManager* geoman = new TGeoManager("world", "geometry");
  
  TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
  TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
  TGeoMaterial* matKapton = new TGeoMaterial("Kapton", 12, 6, 1.42);
  TGeoMedium* medKapton = new TGeoMedium("Kapton", 2, matKapton);
  TGeoMaterial* matCopper = new TGeoMaterial("Copper", 63, 29, 8.94);
  TGeoMedium* medCopper = new TGeoMedium("Copper", 3, matCopper);
  TGeoVolume* volTop = geoman->MakeBox("TOP", 
				       medVacuum, pitch, pitch, 0.02);
  volTop->SetVisibility(0);
  TGeoBBox* shpKapton = new TGeoBBox("K", pitch / 2., 
  				     pitch / 2., 
  				     kapton / 2.);
  TGeoPcon* shpHole = new TGeoPcon("H", 0., 360., 3);
  shpHole->DefineSection(0, -kapton / 2., 0., outdia / 2.);
  shpHole->DefineSection(1,           0., 0., middia / 2.);
  shpHole->DefineSection(2,  kapton / 2., 0., outdia / 2.);
  
  TGeoCompositeShape* shpGem = new TGeoCompositeShape("G", "K - H");
  TGeoVolume* volKapton = new TGeoVolume("Kapton", shpGem, medKapton);
  volKapton->SetLineColor(kGreen);
  volKapton->SetTransparency(50);
  
  TGeoBBox* shpMetal = new TGeoBBox("M", pitch / 2., 
  				    pitch / 2., 
                                           metal / 2.);
  TGeoTube* shpTube = new TGeoTube("T", 0., outdia / 2., metal / 2.);
  TGeoCompositeShape* shpElectrode = new TGeoCompositeShape("E", "M - T");
  TGeoVolume* volElectrode = new TGeoVolume("Electrode", 
  					    shpElectrode, medCopper);
  volElectrode->SetLineColor(kBlue);
  volElectrode->SetTransparency(50);
  
  TGeoVolumeAssembly* volGem = new TGeoVolumeAssembly("Gem");
  const double shift =  0.5 * (metal + kapton);

  volGem->AddNode(volKapton, 1);
  volGem->AddNode(volElectrode, 2, new TGeoTranslation(0., 0.,  shift));
  volGem->AddNode(volElectrode, 3, new TGeoTranslation(0., 0., -shift));

  double Pitch=0.08;
  
  TGeoBBox* DriftAnode = new TGeoBBox("Drift",2*Pitch, kapton / 2., 0.025);
  TGeoVolume* volAnode = new TGeoVolume("Anode", DriftAnode, medKapton);

  TGeoBBox* Readout = new TGeoBBox("Readout",2*Pitch, kapton / 2., 0.025);
  TGeoVolume* volReadout = new TGeoVolume("volReadout", Readout, medKapton);
  
  volAnode->SetLineColor(kMagenta);
  volReadout->SetLineColor(kMagenta);
  
  volAnode->SetTransparency(50);
  volReadout->SetTransparency(50);

  
  // volTop->AddNode(volAnode,22,new TGeoTranslation(0., 0.3,  0));
  // volTop->AddNode(volAnode,23,new TGeoTranslation(0., -0.4,  0));

  TGeoRotation r1;
  r1.SetAngles(0,90,0);

  Float_t pos_gem[3]={0.0,-0.1,-0.3};

  for(int m=0;m<3;m++){

    
    for(int j=0;j<6;j++){
      TGeoTranslation t1(j*pitch,pos_gem[m],0);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+1,ph1);
    }

    for(int j=1;j<6;j++){
      TGeoTranslation t1(-j*pitch,pos_gem[m],0);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+5,ph1);
    }


    for(int j=0;j<6;j++){
      TGeoTranslation t1(pitch*(1/2.)+pitch*(j*1.0),pos_gem[m],sqrt(3)*pitch/2.);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+11,ph1);
    }

    for(int j=0;j<6;j++){
      TGeoTranslation t1(pitch*(1/2.)+pitch*(j*1.0),pos_gem[m],-sqrt(3)*pitch/2.);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+17,ph1);
    }


    for(int j=0;j<6;j++){
      TGeoTranslation t1(-(pitch*(1/2.)+pitch*(j*1.0)),pos_gem[m],sqrt(3)*pitch/2.);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+23,ph1);
    }


    for(int j=0;j<6;j++){
      TGeoTranslation t1(-(pitch*(1/2.)+pitch*(j*1.0)),pos_gem[m],-sqrt(3)*pitch/2.);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+29,ph1);
    }

    for(int j=0;j<6;j++){
      TGeoTranslation t1(pitch*j,pos_gem[m],sqrt(3)*pitch);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+35,ph1);
    }

    for(int j=1;j<6;j++){
      TGeoTranslation t1(-pitch*j,pos_gem[m],sqrt(3)*pitch);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+40,ph1);
    }


    for(int j=0;j<6;j++){
      TGeoTranslation t1(pitch*j,pos_gem[m],-sqrt(3)*pitch);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+46,ph1);
    }

    for(int j=1;j<6;j++){
      TGeoTranslation t1(-pitch*j,pos_gem[m],-sqrt(3)*pitch);
      TGeoCombiTrans c1(t1,r1);
      TGeoHMatrix h1 = c1;
      TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
      volTop->AddNode(volGem,j+52,ph1);
    }

  }

 volTop->AddNode(volAnode,58,new TGeoTranslation(0., 0.3,  0.));
 volTop->AddNode(volReadout,59,new TGeoTranslation(0.,-0.4, 0.));


 TGeoBBox* strip = new TGeoBBox("strip",0.034,0.005,0.02);
 TGeoVolume* volstrip  = new TGeoVolume("strip",strip,medCopper);
 volstrip->SetLineColor(kRed);
 volstrip->SetTransparency(50);

 volTop->AddNode(volstrip,60,new TGeoTranslation(0,-0.395,0.));
 volTop->AddNode(volstrip,61,new TGeoTranslation(0.08,-0.395,0.));
 volTop->AddNode(volstrip,62,new TGeoTranslation(-0.08,-0.395,0.));


 //  volTop->AddNode(volstrip,9,new TGeoTranslation(2*0.08,-0.395,0.));
 //  volTop->AddNode(volstrip,10,new TGeoTranslation(2*(-0.08),-0.395,0.));





// TGeoBBox* DriftAnode = new TGeoBBox("TOP",2*pitch,metal/2.,2*pitch);
//   TGeoVolume* volAnode = new TGeoVolume("Anode", DriftAnode, medKapton);

//   TGeoBBox* Readout = new TGeoBBox("BOTTOM",2*pitch,metal/2.,2*pitch);
//   TGeoVolume* volReadout = new TGeoVolume("volReadout", DriftAnode, medKapton);


//   volTop->AddNode(volAnode,22,new TGeoTranslation(0., 0.3,  0));
//   volTop->AddNode(volAnode,23,new TGeoTranslation(0., -0.4,  0));
  
//   TGeoTranslation t1(0,0,0);
//   TGeoCombiTrans c1(t1,r1);
//   TGeoHMatrix h1 = c1;
//   TGeoHMatrix *ph1 = new TGeoHMatrix(h1);
//   volTop->AddNode(volGem, 1,ph1);

//   TGeoTranslation t2(-pitch, 0., 0.0);
//   TGeoCombiTrans c2(t2,r1);
//   TGeoHMatrix h2 = c2;
//   TGeoHMatrix *ph2 = new TGeoHMatrix(h2);
//   volTop->AddNode(volGem, 2, ph2);

//   TGeoTranslation t3(+pitch, 0., 0.);
//   TGeoCombiTrans c3(t3,r1);
//   TGeoHMatrix h3 = c3;
//   TGeoHMatrix *ph3 = new TGeoHMatrix(h3);
//   volTop->AddNode(volGem, 3, ph3);

//   TGeoTranslation t4(-pitch / 2.,0,sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c4(t4,r1);
//   TGeoHMatrix h4 = c4;
//   TGeoHMatrix *ph4 = new TGeoHMatrix(h4);
//   volTop->AddNode(volGem, 4, ph4);


//   TGeoTranslation t5(+pitch / 2.,0,sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c5(t5,r1);
//   TGeoHMatrix h5 = c5;
//   TGeoHMatrix *ph5 = new TGeoHMatrix(h5);
//   volTop->AddNode(volGem, 5, ph5);

//   TGeoTranslation t6(-pitch / 2.,0,-sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c6(t6,r1);
//   TGeoHMatrix h6 = c6;
//   TGeoHMatrix *ph6 = new TGeoHMatrix(h6);
//   volTop->AddNode(volGem, 6, ph6);

//   TGeoTranslation t7(+pitch / 2.,0,-sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c7(t7,r1);
//   TGeoHMatrix h7 = c7;
//   TGeoHMatrix *ph7 = new TGeoHMatrix(h7);
//   volTop->AddNode(volGem, 7, ph7);

//  TGeoTranslation t8(-2*pitch, 0., 0.0);
//  TGeoCombiTrans c8(t8,r1);
//  TGeoHMatrix h8 = c8;
//  TGeoHMatrix *ph8 = new TGeoHMatrix(h8);
//  volTop->AddNode(volGem, 24, ph8);

//  TGeoTranslation t9(2*pitch, 0., 0.0);
//  TGeoCombiTrans c9(t9,r1);
//  TGeoHMatrix h9 = c9;
//  TGeoHMatrix *ph9 = new TGeoHMatrix(h9);
//  volTop->AddNode(volGem, 25, ph9);

//  TGeoTranslation t10(-(3/2.)*pitch,0,-sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c10(t10,r1);
//  TGeoHMatrix h10 = c10;
//  TGeoHMatrix *ph10 = new TGeoHMatrix(h10);
//  volTop->AddNode(volGem, 26, ph10);

//  TGeoTranslation t11(+(3/2.)*pitch,0,-sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c11(t11,r1);
//  TGeoHMatrix h11 = c11;
//  TGeoHMatrix *ph11 = new TGeoHMatrix(h11);
//  volTop->AddNode(volGem, 27, ph11);


//  TGeoTranslation t12(-(3/2.)*pitch,0,sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c12(t12,r1);
//  TGeoHMatrix h12 = c12;
//  TGeoHMatrix *ph12 = new TGeoHMatrix(h12);
//  volTop->AddNode(volGem, 28, ph12);

//  TGeoTranslation t13(+(3/2.)*pitch,0,sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c13(t13,r1);
//  TGeoHMatrix h13 = c13;
//  TGeoHMatrix *ph13 = new TGeoHMatrix(h13);
//  volTop->AddNode(volGem, 29, ph13);


//  TGeoTranslation t14(0,0,sqrt(3)*pitch);
//  TGeoCombiTrans c14(t14,r1);
//  TGeoHMatrix h14 = c14;
//  TGeoHMatrix *ph14 = new TGeoHMatrix(h14);
//  volTop->AddNode(volGem, 30, ph14);

//  TGeoTranslation t15(0,0,-sqrt(3)*pitch);
//  TGeoCombiTrans c15(t15,r1);
//  TGeoHMatrix h15 = c15;
//  TGeoHMatrix *ph15 = new TGeoHMatrix(h15);
//  volTop->AddNode(volGem, 31, ph15);


//  TGeoTranslation t16(-pitch,0,sqrt(3)*pitch);
//  TGeoCombiTrans c16(t16,r1);
//  TGeoHMatrix h16 = c16;
//  TGeoHMatrix *ph16 = new TGeoHMatrix(h16);
//  volTop->AddNode(volGem, 32, ph16);

//  TGeoTranslation t17(-pitch,0,-sqrt(3)*pitch);
//  TGeoCombiTrans c17(t17,r1);
//  TGeoHMatrix h17 = c17;
//  TGeoHMatrix *ph17 = new TGeoHMatrix(h17);
//  volTop->AddNode(volGem, 33, ph17);

//  TGeoTranslation t18(+pitch,0,sqrt(3)*pitch);
//  TGeoCombiTrans c18(t18,r1);
//  TGeoHMatrix h18 = c18;
//  TGeoHMatrix *ph18 = new TGeoHMatrix(h18);
//  volTop->AddNode(volGem, 34, ph18);

//  TGeoTranslation t19(+pitch,0,-sqrt(3)*pitch);
//  TGeoCombiTrans c19(t19,r1);
//  TGeoHMatrix h19 = c19;
//  TGeoHMatrix *ph19 = new TGeoHMatrix(h19);
//  volTop->AddNode(volGem, 35, ph19);


//  ////////////////////////////////
 

//   // TGeoTranslation t3(+pitch, 0., 0.);
//   // TGeoCombiTrans c3(t3,r1);
//   // TGeoHMatrix h3 = c3;
//   // TGeoHMatrix *ph3 = new TGeoHMatrix(h3);
//   // volTop->AddNode(volGem, 3, ph3);


//   TGeoTranslation t12p(0,-0.1,0);
//   TGeoCombiTrans c12p(t12p,r1);
//   TGeoHMatrix h12p = c12p;
//   TGeoHMatrix *ph12p = new TGeoHMatrix(h12p);
//   volTop->AddNode(volGem, 8,ph12p);

//   TGeoTranslation t22(-pitch, -0.1, 0.0);
//   TGeoCombiTrans c22(t22,r1);
//   TGeoHMatrix h22 = c22;
//   TGeoHMatrix *ph22 = new TGeoHMatrix(h22);
//   volTop->AddNode(volGem, 9, ph22);

//   TGeoTranslation t32(+pitch, -0.1, 0.);
//   TGeoCombiTrans c32(t32,r1);
//   TGeoHMatrix h32 = c32;
//   TGeoHMatrix *ph32 = new TGeoHMatrix(h32);
//   volTop->AddNode(volGem, 10, ph32);

//   TGeoTranslation t42(-pitch / 2.,-0.1,sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c42(t42,r1);
//   TGeoHMatrix h42 = c42;
//   TGeoHMatrix *ph42 = new TGeoHMatrix(h42);
//   volTop->AddNode(volGem, 11, ph42);


//   TGeoTranslation t52(+pitch / 2.,-0.1,sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c52(t52,r1);
//   TGeoHMatrix h52 = c52;
//   TGeoHMatrix *ph52 = new TGeoHMatrix(h52);
//   volTop->AddNode(volGem, 12, ph52);

//   TGeoTranslation t62(-pitch / 2.,-0.1,-sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c62(t62,r1);
//   TGeoHMatrix h62 = c62;
//   TGeoHMatrix *ph62 = new TGeoHMatrix(h62);
//   volTop->AddNode(volGem, 13, ph62);

//   TGeoTranslation t72(+pitch / 2.,-0.1,-sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c72(t72,r1);
//   TGeoHMatrix h72 = c72;
//   TGeoHMatrix *ph72 = new TGeoHMatrix(h72);
//   volTop->AddNode(volGem, 14, ph72);



//   TGeoTranslation t82(-2*pitch, -0.1, 0.0);
//   TGeoCombiTrans c82(t82,r1);
//   TGeoHMatrix h82 = c82;
//   TGeoHMatrix *ph82 = new TGeoHMatrix(h82);
//   volTop->AddNode(volGem, 24, ph82);

//  TGeoTranslation t92(2*pitch, -0.1, 0.0);
//  TGeoCombiTrans c92(t92,r1);
//  TGeoHMatrix h92 = c92;
//  TGeoHMatrix *ph92 = new TGeoHMatrix(h92);
//  volTop->AddNode(volGem, 25, ph9);

//  TGeoTranslation t102(-(3/2.)*pitch,-0.1,-sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c102(t102,r1);
//  TGeoHMatrix h102 = c102;
//  TGeoHMatrix *ph102 = new TGeoHMatrix(h102);
//  volTop->AddNode(volGem, 26, ph102);

//  TGeoTranslation t112(+(3/2.)*pitch,-0.1,-sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c112(t112,r1);
//  TGeoHMatrix h112 = c112;
//  TGeoHMatrix *ph112 = new TGeoHMatrix(h112);
//  volTop->AddNode(volGem, 27, ph112);


//  TGeoTranslation t122(-(3/2.)*pitch,-0.1,sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c122(t122,r1);
//  TGeoHMatrix h122 = c122;
//  TGeoHMatrix *ph122 = new TGeoHMatrix(h122);
//  volTop->AddNode(volGem, 28, ph122);

//  TGeoTranslation t132(+(3/2.)*pitch,-0.1,sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c132(t132,r1);
//  TGeoHMatrix h132 = c132;
//  TGeoHMatrix *ph132 = new TGeoHMatrix(h132);
//  volTop->AddNode(volGem, 29, ph132);


//  TGeoTranslation t142(0,-0.1,sqrt(3)*pitch);
//  TGeoCombiTrans c142(t142,r1);
//  TGeoHMatrix h142 = c142;
//  TGeoHMatrix *ph142 = new TGeoHMatrix(h142);
//  volTop->AddNode(volGem, 30, ph142);

//  TGeoTranslation t152(0,-0.1,-sqrt(3)*pitch);
//  TGeoCombiTrans c152(t152,r1);
//  TGeoHMatrix h152 = c152;
//  TGeoHMatrix *ph152 = new TGeoHMatrix(h152);
//  volTop->AddNode(volGem, 31, ph152);


//  TGeoTranslation t162(-pitch,-0.1,sqrt(3)*pitch);
//  TGeoCombiTrans c162(t162,r1);
//  TGeoHMatrix h162 = c162;
//  TGeoHMatrix *ph162 = new TGeoHMatrix(h162);
//  volTop->AddNode(volGem, 32, ph162);

//  TGeoTranslation t172(-pitch,-0.1,-sqrt(3)*pitch);
//  TGeoCombiTrans c172(t172,r1);
//  TGeoHMatrix h172 = c172;
//  TGeoHMatrix *ph172 = new TGeoHMatrix(h172);
//  volTop->AddNode(volGem, 33, ph172);

//  TGeoTranslation t182(+pitch,-0.1,sqrt(3)*pitch);
//  TGeoCombiTrans c182(t182,r1);
//  TGeoHMatrix h182 = c182;
//  TGeoHMatrix *ph182 = new TGeoHMatrix(h182);
//  volTop->AddNode(volGem, 34, ph182);

//  TGeoTranslation t192(+pitch,-0.1,-sqrt(3)*pitch);
//  TGeoCombiTrans c192(t192,r1);
//  TGeoHMatrix h192 = c192;
//  TGeoHMatrix *ph192 = new TGeoHMatrix(h192);
//  volTop->AddNode(volGem, 35, ph192);



  

//   /////////////////////////////////////////////////////////////////
  
//   TGeoTranslation t13p(0,-0.3,0);
//   TGeoCombiTrans c13p(t13p,r1);
//   TGeoHMatrix h13p = c13p;
//   TGeoHMatrix *ph13p = new TGeoHMatrix(h13p);
//   volTop->AddNode(volGem, 15,ph13p);

//   TGeoTranslation t23(-pitch, -0.3, 0.0);
//   TGeoCombiTrans c23(t23,r1);
//   TGeoHMatrix h23 = c23;
//   TGeoHMatrix *ph23 = new TGeoHMatrix(h23);
//   volTop->AddNode(volGem, 16, ph23);

//   TGeoTranslation t33(+pitch, -0.3, 0.);
//   TGeoCombiTrans c33(t33,r1);
//   TGeoHMatrix h33 = c33;
//   TGeoHMatrix *ph33 = new TGeoHMatrix(h33);
//   volTop->AddNode(volGem, 17, ph33);

//   TGeoTranslation t43(-pitch / 2.,-0.3,sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c43(t43,r1);
//   TGeoHMatrix h43 = c43;
//   TGeoHMatrix *ph43 = new TGeoHMatrix(h43);
//   volTop->AddNode(volGem, 18, ph43);


//   TGeoTranslation t53(+pitch / 2.,-0.3,sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c53(t53,r1);
//   TGeoHMatrix h53 = c53;
//   TGeoHMatrix *ph53 = new TGeoHMatrix(h53);
//   volTop->AddNode(volGem, 19, ph53);

//   TGeoTranslation t63(-pitch / 2.,-0.3,-sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c63(t63,r1);
//   TGeoHMatrix h63 = c63;
//   TGeoHMatrix *ph63 = new TGeoHMatrix(h63);
//   volTop->AddNode(volGem, 20, ph63);

//   TGeoTranslation t73(+pitch / 2.,-0.3,-sqrt(3) * pitch / 2.);
//   TGeoCombiTrans c73(t73,r1);
//   TGeoHMatrix h73 = c73;
//   TGeoHMatrix *ph73 = new TGeoHMatrix(h73);
//   volTop->AddNode(volGem, 21, ph73);
  
//    TGeoTranslation t83(-2*pitch, -0.3, 0.0);
//   TGeoCombiTrans c83(t83,r1);
//   TGeoHMatrix h83 = c83;
//   TGeoHMatrix *ph83 = new TGeoHMatrix(h83);
//   volTop->AddNode(volGem, 24, ph83);

//  TGeoTranslation t93(2*pitch, -0.3, 0.0);
//  TGeoCombiTrans c93(t93,r1);
//  TGeoHMatrix h93 = c93;
//  TGeoHMatrix *ph93 = new TGeoHMatrix(h93);
//  volTop->AddNode(volGem, 25, ph93);

//  TGeoTranslation t103(-(3/2.)*pitch,-0.3,-sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c103(t103,r1);
//  TGeoHMatrix h103 = c103;
//  TGeoHMatrix *ph103 = new TGeoHMatrix(h103);
//  volTop->AddNode(volGem, 26, ph103);

//  TGeoTranslation t113(+(3/2.)*pitch,-0.3,-sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c113(t113,r1);
//  TGeoHMatrix h113 = c113;
//  TGeoHMatrix *ph113 = new TGeoHMatrix(h113);
//  volTop->AddNode(volGem, 27, ph113);


//  TGeoTranslation t123(-(3/2.)*pitch,-0.3,sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c123(t123,r1);
//  TGeoHMatrix h123 = c123;
//  TGeoHMatrix *ph123 = new TGeoHMatrix(h123);
//  volTop->AddNode(volGem, 28, ph123);

//  TGeoTranslation t133(+(3/2.)*pitch,-0.3,sqrt(3) * pitch / 2.);
//  TGeoCombiTrans c133(t133,r1);
//  TGeoHMatrix h133 = c133;
//  TGeoHMatrix *ph133 = new TGeoHMatrix(h133);
//  volTop->AddNode(volGem, 29, ph133);


//  TGeoTranslation t143(0,-0.3,sqrt(3)*pitch);
//  TGeoCombiTrans c143(t143,r1);
//  TGeoHMatrix h143 = c143;
//  TGeoHMatrix *ph143 = new TGeoHMatrix(h143);
//  volTop->AddNode(volGem, 30, ph143);

//  TGeoTranslation t153(0,-0.3,-sqrt(3)*pitch);
//  TGeoCombiTrans c153(t153,r1);
//  TGeoHMatrix h153 = c153;
//  TGeoHMatrix *ph153 = new TGeoHMatrix(h153);
//  volTop->AddNode(volGem, 31, ph153);


//  TGeoTranslation t163(-pitch,-0.3,sqrt(3)*pitch);
//  TGeoCombiTrans c163(t163,r1);
//  TGeoHMatrix h163 = c163;
//  TGeoHMatrix *ph163 = new TGeoHMatrix(h163);
//  volTop->AddNode(volGem, 32, ph163);

//  TGeoTranslation t173(-pitch,-0.3,-sqrt(3)*pitch);
//  TGeoCombiTrans c173(t173,r1);
//  TGeoHMatrix h173 = c173;
//  TGeoHMatrix *ph173 = new TGeoHMatrix(h173);
//  volTop->AddNode(volGem, 33, ph173);

//  TGeoTranslation t183(+pitch,-0.3,sqrt(3)*pitch);
//  TGeoCombiTrans c183(t183,r1);
//  TGeoHMatrix h183 = c183;
//  TGeoHMatrix *ph183 = new TGeoHMatrix(h183);
//  volTop->AddNode(volGem, 34, ph183);

//  TGeoTranslation t193(+pitch,-0.3,-sqrt(3)*pitch);
//  TGeoCombiTrans c193(t193,r1);
//  TGeoHMatrix h193 = c193;
//  TGeoHMatrix *ph193 = new TGeoHMatrix(h193);
//  volTop->AddNode(volGem, 35, ph193);



    geoman->SetVerboseLevel(0);
    geoman->SetTopVolume(volTop);
    geoman->CloseGeometry();
    geoman->CheckOverlaps(0.1e-4);
    geoman->SetNmeshPoints(100000);
    cD->cd();
    geoman->GetTopVolume()->Draw("ogl");
  }



  cD->cd();
  if (plotDrift) {
    driftView->SetCanvas(cD);
    driftView->Plot();
  }
  

  //  viewdrift->Plot();
  //write the results in to the ntuple
  f->cd();
  hElectrons->Write();
  ntuple->Write();
  f->Close();
   
 
  app.Run(kTRUE);

  //end
}
