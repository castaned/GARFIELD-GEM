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
#include <TNtuple.h>
#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "RandomEngineMPIServer.hh"
#include "ComponentAnalyticField.hh"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <string>

#define DIETAG      1
#define WORKTAG     2
#define RESULT_TAG  3
#define RESULT_MORE_TAG 4

using namespace Garfield;
using namespace std;


void find_and_replace(string& source, string const& find, string const& replace)
{
  for(string::size_type i = 0; (i = source.find(find, i)) != string::npos;)
    {
      source.replace(i, find.length(), replace);
      i += replace.length();
    }
}

typedef struct EventResult_s {
  int eventID;
  
  int ne;
  int ni;
  
  // the ntuple consists of xe1,ye1,ze1,xe2,ye2,ze2,e2
  int size_ntuple;
  
  double start_time;
  double finish_time;
  
} EventResult;

MPI_Datatype mpi_event_result_type;

/* Local functions */
static void master(int argc, char** argv, int size, int nEvents, std::string ansys);
static void slave(int argc, char** argv, int myrank, int size, int random_server);
static void random_number_server(void);
void create_mpi_event_result_type(void);

int main(int argc, char **argv)
{
  int myrank, size;
  
  int nEvents = 0;
  
  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  nEvents = atoi(argv[1]);
  if(myrank==0) cout << "NUMBER OF EVENTS: " << nEvents << endl;

  if(argc<5) {
    cerr << "Command line error: Arguments required." << endl;
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
  
  std::string ansys = argv[4];
  
  if (myrank == 0) {
    master(argc, argv, size, nEvents, ansys);
  } else if (myrank == size-1) {
    random_number_server();
  } else {
    slave(argc, argv,myrank, size, size-1);
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

// Dimensions of the detector
const double Pitch = 0.088; // Pitch of strips
const double pitch = 0.014; // pitch of holes
const double smear = pitch/2.;
double singlecell  = 0.014;

// Magnetic field 
const double MagX = 0.;
const double MagY = 0.;
const double MagZ = 0.;

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

void master(int argc, char** argv, int size, int nEvents, std::string ansys) {
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Histograms
  int nBinsGain = 5000;
  double gmin   = 0.;
  double gmax   = 5000.;

  //  const std::string ansys = "/panfs/vol/HEP/GarfieldSim/PRO/ANSYS/gasgapscan/gap02pv2867/";
  
  std::string output_root_file  = ansys.substr(37)+"_out.root";
  
  find_and_replace(output_root_file,"/","_");
  
  char nameroot[200];
  strcpy(nameroot,output_root_file.c_str());
  
  TFile* f = new TFile(nameroot,"RECREATE");
  
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons", nBinsGain, gmin, gmax);
  
  TNtuple* ntuple  = new TNtuple("ntuple","ntuple","xe1:ye1:ze1:xe2:ye2:ze2:e2");
  
  vector<double> v_xe1, v_ye1, v_ze1, v_xe2, v_ye2, v_ze2, v_e2;

  EventResult res;
  MPI_Status status, status2;

  ofstream out;

  std::string out_log = ansys.substr(37)+"_out.dat";
  find_and_replace(out_log,"/","_");

  char namelog[200];
  strcpy(namelog,out_log.c_str());

  out.open(namelog);

  int iEvent=0, iCompleted=0;
  for(; iEvent<size-2; iEvent++) {
    int slave_id = iEvent + 1;
    cout << "######## SEND : " << iEvent << endl;
    MPI_Send(&iEvent, 1, MPI_INT, slave_id, WORKTAG, MPI_COMM_WORLD);
  }

  while(iCompleted<nEvents) {
    
    MPI_Recv(&res, 1, mpi_event_result_type, 
	     MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
    
    cout << "######## RECV : " << res.eventID << endl;
    
    v_xe1.resize(res.size_ntuple);
    v_ye1.resize(res.size_ntuple);
    v_ze1.resize(res.size_ntuple);
    v_xe2.resize(res.size_ntuple);
    v_ye2.resize(res.size_ntuple);
    v_ze2.resize(res.size_ntuple);
    v_e2.resize(res.size_ntuple);

    // receive xe1
    MPI_Recv(&v_xe1[0], res.size_ntuple, MPI_DOUBLE, status.MPI_SOURCE, 
	     RESULT_MORE_TAG, MPI_COMM_WORLD, &status2);
    // receive ye1
    MPI_Recv(&v_ye1[0], res.size_ntuple, MPI_DOUBLE, status.MPI_SOURCE, 
	     RESULT_MORE_TAG, MPI_COMM_WORLD, &status2);
    // receive ze1
    MPI_Recv(&v_ze1[0], res.size_ntuple, MPI_DOUBLE, status.MPI_SOURCE, 
	     RESULT_MORE_TAG, MPI_COMM_WORLD, &status2);
    // receive xe2
    MPI_Recv(&v_xe2[0], res.size_ntuple, MPI_DOUBLE, status.MPI_SOURCE, 
	     RESULT_MORE_TAG, MPI_COMM_WORLD, &status2);
    // receive ye2
    MPI_Recv(&v_ye2[0], res.size_ntuple, MPI_DOUBLE, status.MPI_SOURCE, 
	     RESULT_MORE_TAG, MPI_COMM_WORLD, &status2);
    // receive ze2
    MPI_Recv(&v_ze2[0], res.size_ntuple, MPI_DOUBLE, status.MPI_SOURCE, 
	     RESULT_MORE_TAG, MPI_COMM_WORLD, &status2);
    // receive e2
    MPI_Recv(&v_e2[0], res.size_ntuple, MPI_DOUBLE, status.MPI_SOURCE, 
	     RESULT_MORE_TAG, MPI_COMM_WORLD, &status2);

    out << (res.finish_time-res.start_time) << " "
        << "event=" << res.eventID << " "
        << "slave=" << status.MPI_SOURCE << " "
        << "ne=" <<  res.ne << endl;

    // Accumulate the results
    hElectrons->Fill(res.ne);

    for(int iTuple=0; iTuple<res.size_ntuple; iTuple++) {
      ntuple->Fill(v_xe1[iTuple],v_ye1[iTuple],v_ze1[iTuple],
		   v_xe2[iTuple],v_ye2[iTuple],v_ze2[iTuple],
		   v_e2[iTuple]);
    }

    iCompleted++;

    // Send more work (if any left) to the same slave who just finished the work
    if(iEvent<nEvents) {
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

void do_work(int eventID, EventResult& res, 
	     vector<double>& v_xe1, vector<double>& v_ye1, vector<double>& v_ze1,
	     vector<double>& v_xe2, vector<double>& v_ye2, vector<double>& v_ze2,
	     vector<double>& v_e2);

void slave(int argc, char** argv,int myrank, int size, int random_server) {

  // Load the field map.
  fm = new ComponentAnsys123();
  
  const std::string ansys_dir = argv[4];


  const std::string efile = ansys_dir+"/ELIST.lis";
  const std::string nfile = ansys_dir+"/NLIST.lis";
  const std::string mfile = ansys_dir+"/MPLIST.lis";
  const std::string sfile = ansys_dir+"/PRNSOL.lis";
  
  fm->Initialise(efile, nfile, mfile, sfile, "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityZ();
  fm->PrintRange();
  
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();


  const std::string noble_gas  = argv[2];
  
  
  // Set the Penning transfer efficiency.
  double rPenning = 0.0;
  rPenning =  atof(argv[3]);
  cout<<" penning transfer  "<<rPenning<<endl;

  cout<<" noble gas  "<<noble_gas<<endl;
  gas->SetComposition(noble_gas, 70., "co2", 30.);
  
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->Initialise();  

  const double lambdaPenning = 0.;
  if(noble_gas == "ar") gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  if(noble_gas == "ne") gas->EnablePenningTransfer(rPenning, lambdaPenning, "ne");
  
  // Load the ion mobilities.

  if(noble_gas == "ar")  gas->LoadIonMobility("IonMobility/IonMobility_Ar+_Ar.txt");
  if(noble_gas == "ne")  gas->LoadIonMobility("IonMobility/mob_Ne_Ne+.txt");
  
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
  
  //  cmpAmp->AddPlaneY(-0.30275, 1., "Driftplane");
  cmpAmp->AddPlaneY(0.30275, 1., "Driftplane");
  cmpAmp->AddPlaneY(-0.4    , 0., "striplane");
  
  //Next we construct the Strips for readout of te signal, also with labels
  //  double  Xstrip1 = -0.088, Xstrip2 = 0.0, Xstrip3 = 0.088, Xstrip4 = 0.156 ; //Store the center of strips
  double  Xstrip1 = -0.088, Xstrip2 = 0.0, Xstrip3 = 0.088, Xstrip4 = 0.176 ; //Store the center of strips
  cmpAmp->AddStripOnPlaneY('z', -0.4, -0.122 , -0.054, "Strip1");
  cmpAmp->AddStripOnPlaneY('z', -0.4, -0.034,   0.034, "Strip2");
  //  cmpAmp->AddStripOnPlaneY('z', -0.4,  0.054,   0.102, "Strip3");
  //  cmpAmp->AddStripOnPlaneY('z', -0.4,  0.122,   0.190, "Strip4");
  
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
  sensor = new Sensor();
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

  aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);
  aval->EnableSignalCalculation(); 
  aval->SetTimeWindow(tStart,tStop);
  //aval->EnableAvalancheSizeLimit(1000);
  
  drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);
  
  sensor->ClearSignal();
  
  EventResult res;
  memset(&res, 0, sizeof(EventResult_s));

  vector<double> v_xe1, v_ye1, v_ze1, v_xe2, v_ye2, v_ze2, v_e2;
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
      v_xe1.clear();
      v_ye1.clear();
      v_ze1.clear();
      v_xe2.clear();
      v_ye2.clear();
      v_ze2.clear();
      v_e2.clear();

      memset(&res, 0, sizeof(EventResult_s));
      do_work(eventID, res, v_xe1, v_ye1, v_ze1, v_xe2, v_ye2, v_ze2, v_e2);
      
      /* Send the result back */
      MPI_Send(&res, 1, mpi_event_result_type, 0, RESULT_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_xe1[0], v_xe1.size(), MPI_DOUBLE, 0, RESULT_MORE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ye1[0], v_ye1.size(), MPI_DOUBLE, 0, RESULT_MORE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ze1[0], v_ze1.size(), MPI_DOUBLE, 0, RESULT_MORE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_xe2[0], v_xe2.size(), MPI_DOUBLE, 0, RESULT_MORE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ye2[0], v_ye2.size(), MPI_DOUBLE, 0, RESULT_MORE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_ze2[0], v_ze2.size(), MPI_DOUBLE, 0, RESULT_MORE_TAG, MPI_COMM_WORLD);
      MPI_Send(&v_e2[0], v_e2.size(), MPI_DOUBLE, 0, RESULT_MORE_TAG, MPI_COMM_WORLD);
    } else {
      cerr << "Slave: Unknow MPI tag" << endl;
    }
  }
}

void do_work(int eventID, EventResult& res, 
	     vector<double>& v_xe1, vector<double>& v_ye1, vector<double>& v_ze1,
	     vector<double>& v_xe2, vector<double>& v_ye2, vector<double>& v_ze2,
	     vector<double>& v_e2) {

  // Randomize the initial position.
  memset(&res, 0, sizeof(EventResult_s));
  res.start_time = MPI_Wtime();
  res.eventID = eventID;

  double x0 = -Pitch + RndmUniform() * 2*Pitch;
  //double x0 = Pitch/2.;
  double y0 = 0.25; 
  double z0 = -smear + RndmUniform() * smear;
  //double z0 = pitch/2.;
  double t0 = 0.;
  double e0 = 0.5;  
  
  aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
  
  int ne = 0, ni = 0;
  aval->GetAvalancheSize(ne, ni);
   
  res.ne = ne;
  res.ni = ni;
  
  const int np = aval->GetNumberOfElectronEndpoints();
  double xe1, ye1, ze1, te1, e1;
  double xe2, ye2, ze2, te2, e2;
  int status;

  res.size_ntuple = aval->GetNumberOfElectronEndpoints();

  for (int j = np; j--;) {
    aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
    
    v_xe1.push_back(xe1);
    v_ye1.push_back(ye1);
    v_ze1.push_back(ze1);
    v_xe2.push_back(xe2);
    v_ye2.push_back(ye2);
    v_ze2.push_back(ze2);
    v_e2.push_back(e2);
  }
  
  res.finish_time = MPI_Wtime();
}

void create_mpi_event_result_type(void) {
  
  EventResult res;

  const int nitems=6;
  int blocklengths[nitems] = {1,1,1,1,1,1};
  MPI_Datatype old_types[nitems] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE};
  MPI_Aint displ[nitems];
  MPI_Get_address(&res.eventID, &displ[0]);
  MPI_Get_address(&res.ne, &displ[1]);
  MPI_Get_address(&res.ni, &displ[2]);
  MPI_Get_address(&res.size_ntuple, &displ[3]);
  MPI_Get_address(&res.start_time, &displ[4]);
  MPI_Get_address(&res.finish_time, &displ[5]);
  
  for(int i=nitems-1; i>=0; i--)
    displ[i] -= displ[0];

  MPI_Type_create_struct(nitems, blocklengths, displ, old_types, &mpi_event_result_type);
  MPI_Type_commit(&mpi_event_result_type);
}
