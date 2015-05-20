using namespace std;
#include <mpi.h>

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
#include <TFile.h>

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include <TTree.h>
#include <TFile.h>
#include <TNtuple.h>
#include "RandomEngineMPIServer.hh"

#include <stdlib.h>
#include <stdio.h>

#include <vector>

#define DIETAG      1
#define WORKTAG     2
#define RESULT_TAG  3
#define RESULTE_TAG 4
#define RESULTI_TAG 5

using namespace Garfield;
using namespace std;

typedef struct EventResult_s {
  int eventID;

  int ne;
  int ni;
  double xe1;
  double ye1;
  double ze1;
  double xe2;
  double ye2;
  double ze2;
  double e2;
  

  double sumElectronTotal;
  double sumElectron1Strip;
  double sumElectron2Strip;
  double sumElectronHalfStrip;
  
} EventResult;

MPI_Datatype mpi_event_result_type;

/* Local functions */

static void master(int argc, char** argv, int size, int nEvents);
static void slave(int myrank, int size, int random_server);
static void random_number_server(void);
void create_mpi_event_result_type(void);

int
main(int argc, char **argv)
{
  int myrank, size;

  int nEvents = 0;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  nEvents = atoi(argv[1]);
  if(myrank==0) 
    cout << "NUMBER OF EVENTS: " << nEvents << endl;

  if(argc < 2) {
    cerr << "Command line error: Number of events required." << endl;
    MPI_Finalize();
    return 1;
  }

  if(size<3) {
    cerr << "Error: Requries at least 3 processes (1 master, 1 random_server, 1 or more slaves).\n";
    MPI_Finalize();
    return 1;
  }

  if(nEvents < size-2) {
    cerr << "Error: Number of events should not be less than number of slaves.\n";
    MPI_Finalize();
    return 1; 
  }
  
  // Create a MPI type for EventResult structure
  create_mpi_event_result_type();

  if (myrank == 0) {
    master(argc, argv, size, nEvents);
  } else if (myrank == size-1) {
    random_number_server();
  } else {
    slave(myrank, size, size-1);
  }

  /* Shut down MPI */

  MPI_Finalize();
  return 0;
}


void random_number_server(void) {
  RandomEngineMPIServer server;
  cout << "Random Server: starting the random number generator" << endl;
  server.Run();
  cout << "Random Server: random number generator killed" << endl;
}

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

// Field map
ComponentAnsys123* fm;
// Gas
MediumMagboltz* gas;
// Sensor
Sensor* sensor;
// Microsopic avalanche
AvalancheMicroscopic* aval;
// Drift
AvalancheMC* drift;

void master(int argc, char** argv, int size, int nEvents) {
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Histograms
  int nBinsGain = 100;
  double gmin =   0.;
  double gmax = 100.;
  
  TFile* f = new TFile("v3650_pitch560.root","RECREATE");
  
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons", nBinsGain, gmin, gmax);

  TNtuple* ntuple  = new TNtuple("ntuple","ntuple","ne:xe1:ye1:ze1:xe2:ye2:ze2:e2");

  vector<float> v_chrgE, v_chrgI;
  EventResult res;
  MPI_Status status, status2;
  
  ofstream out;
  out.open("progress.dat");

  int iEvent=0, iCompleted=0;
  for(; iEvent<size-2; iEvent++) 
  {
    int slave_id = iEvent + 1;
    cout << "######## SEND : " << iEvent << endl;
    MPI_Send(&iEvent, 1, MPI_INT, slave_id, WORKTAG, MPI_COMM_WORLD);
  }

  while(iCompleted<nEvents) 
  {
    
    MPI_Recv(&res, 1, mpi_event_result_type, 
      MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);

    out << "Event ID: " << res.eventID << ", ne= " <<  res.ne << endl;

    // Accumulate the results
    hElectrons->Fill(res.ne);
    ntuple->Fill(res.ne,res.xe1,res.ye1,res.ze1,res.xe2,res.ye2,res.ze2,res.e2);

    iCompleted++;

    // Send more work (if any left) to the same slave who just finished the work
    if(iEvent<nEvents) 
    {
      cout << "######## SEND : " << iEvent << endl;
      MPI_Send(&iEvent, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);    
      iEvent++;
    }
  }

  out.close();
  f->Write();
  
  int dummy=0;

  // Send DIETAG to all slaves
  for(int i=1; i<=size-2; i++)
    MPI_Send(&dummy, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);

  // Terminate the random number server
  MPI_Send(&dummy, 1, MPI_INT, size-1, RandomEngineMPI::KILL, MPI_COMM_WORLD );

}

void do_work(int eventID, EventResult& res, vector<float>& v_chrgE, vector<float>& v_chrgI);

void slave(int myrank, int size, int random_server) {

  const bool debug = true;

  // Load the field map.
  fm = new ComponentAnsys123();
  const std::string efile = "/panfs/vol/HEP/GEM/ansys_input/v4050/pitch840/ELIST.lis";
  const std::string nfile = "/panfs/vol/HEP/GEM/ansys_input/v4050/pitch840/NLIST.lis";
  const std::string mfile = "/panfs/vol/HEP/GEM/ansys_input/v4050/pitch840/MPLIST.lis";
  const std::string sfile = "/panfs/vol/HEP/GEM/ansys_input/v4050/pitch840/PRNSOL.lis";
  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  fm->PrintRange();

  // Setup the gas.
  gas = new MediumMagboltz();
  gas->SetComposition("ar", 45., "co2", 15., "cf4", 40);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->EnableDebugging();
  gas->Initialise();
  gas->DisableDebugging();
  // Set the Penning transfer efficiency.
  const double rPenning = 0.4;
  const double lambdaPenning = 0.;
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  gas->LoadIonMobility("/panfs/vol/y/yamaghr76/garfield/Data/IonMobility_Ar+_Ar.txt");

  
  // Associate the gas with the corresponding field map material.
  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }
  fm->PrintMaterials();

  // Create the sensor.
  sensor = new Sensor();
  sensor->AddComponent(fm);
		   
  sensor->SetArea(-5 * pitch, -5 * pitch, -5, 5 * pitch,  5 * pitch,  5);

  aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

  drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);

  
  EventResult res = {0};

  vector<float> v_chrgE, v_chrgI;
  int eventID;
  MPI_Status status;

  while(1) {
    /* Receive a message from the master */
    MPI_Recv(&eventID, 1, MPI_INT, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    /* Check the tag of the received message. */
    if (status.MPI_TAG == DIETAG) {
      break;
    } else if(status.MPI_TAG == WORKTAG) {
      /* Do the work */
      v_chrgE.clear();
      v_chrgI.clear();

      res = {0};
      do_work(eventID, res, v_chrgE, v_chrgI);
      
      /* Send the result back */
      MPI_Send(&res, 1, mpi_event_result_type, 0, RESULT_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_chrgE[0], v_chrgE.size(), MPI_FLOAT, 0, RESULTE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_chrgI[0], v_chrgI.size(), MPI_FLOAT, 0, RESULTI_TAG, MPI_COMM_WORLD);
    } else {
      cerr << "Slave: Unknow MPI tag" << endl;
    }
  }
}

void do_work(int eventID, EventResult& res, vector<float>& v_chrgE, vector<float>& v_chrgI) {
  // Randomize the initial position.
  res = {0};
  res.eventID = eventID;
  const double smear = pitch / 2.;
  double x0 = -smear + RndmUniform() * smear;
  double y0 = -smear + RndmUniform() * smear;
  double z0 = 0.025;
  double t0 = 0.;
  double e0 = 0.1;

  aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
  int ne = 0, ni = 0;
  
  double sumElectronTotal     = 0.;
  double sumElectron1Strip    = 0.;
  double sumElectron2Strip    = 0.;
  double sumElectronHalfStrip = 0.;
  
  
  aval->GetAvalancheSize(ne, ni);
  //hElectrons->Fill(ne);
  //hIons->Fill(ni);
  res.ne = ne;
  res.ni = ni;
  const int np = aval->GetNumberOfElectronEndpoints();
  double xe1, ye1, ze1, te1, e1;
  double xe2, ye2, ze2, te2, e2;
  double xi1, yi1, zi1, ti1;
  double xi2, yi2, zi2, ti2;
  int status;
  for (int j = np; j--;) 
  {
    aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
	
    sumElectronTotal += 1.;
    
    //usual number of electrons collected in the first half strips
    if (ze2 >= -total && ze2 <= -total + thickness && ( (xe2 >= 0. && xe2 <=  stripwidth/2.) || (xe2 >= pitch-stripwidth/2. && xe2 <= pitch) ) )
	  
    sumElectronHalfStrip += 1.;
	
    //number of electrons collected in different strips  
	 if (ze2 >= -total && ze2 <= -total + thickness 
	    && ( (xe2 >= -stripwidth/2. && xe2 <=  stripwidth/2.) || (xe2 >= pitch-stripwidth/2. && xe2 <= pitch+stripwidth/2.) ) )
	    
	 sumElectron1Strip += 1.;
	 
	 if (ze2 >= -total && ze2 <= -total + thickness 
	 && ( (xe2 >= -stripwidth/2.       && xe2 <= stripwidth/2.)         || (xe2 >= pitch-stripwidth/2.   && xe2 <= pitch+stripwidth/2.  ) 
	 ||   (xe2 >= -pitch-stripwidth/2. && xe2 <= -pitch+stripwidth/2.)  || (xe2 >= 2*pitch-stripwidth/2. && xe2 <= 2*pitch+stripwidth/2.) ) )
	   
	 sumElectron2Strip += 1.;

  }    
  
  res.xe1 = xe1;
  res.ye1 = ye1;
  res.ze1 = ze1;
  res.xe2 = xe2;
  res.ye2 = ye2;
  res.ze2 = ze2;
  res.sumElectronHalfStrip = sumElectronHalfStrip;
  res.sumElectron1Strip = sumElectron1Strip;
  res.sumElectron2Strip = sumElectron2Strip;
  


}      







       
void create_mpi_event_result_type(void) {
       
  EventResult res;
       
  const int nitems=14;
  int blocklengths[14] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  MPI_Datatype old_types[14] = {MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                            MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                            MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  MPI_Aint displ[14];
  MPI_Get_address(&res.eventID, &displ[0]);
  MPI_Get_address(&res.ne, &displ[1]);
  MPI_Get_address(&res.xe1, &displ[2]);
  MPI_Get_address(&res.ye1, &displ[3]);
  MPI_Get_address(&res.ze1, &displ[4]);
  MPI_Get_address(&res.xe2, &displ[5]);
  MPI_Get_address(&res.ye2, &displ[6]);
  MPI_Get_address(&res.ze2, &displ[7]);
  MPI_Get_address(&res.e2, &displ[8]);
  MPI_Get_address(&res.sumElectronTotal, &displ[9]);
  MPI_Get_address(&res.sumElectron1Strip, &displ[10]);
  MPI_Get_address(&res.sumElectron2Strip, &displ[11]);
  MPI_Get_address(&res.sumElectronHalfStrip, &displ[12]);
  MPI_Get_address(&res.sumElectron2Strip, &displ[13]);
  
 

  for(int i=13; i>=0; i--)
    displ[i] -= displ[0];

  MPI_Type_create_struct(nitems, blocklengths, displ, old_types, &mpi_event_result_type);
  MPI_Type_commit(&mpi_event_result_type);
}
