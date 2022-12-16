//
// Created by mdupont on 04/10/17.
//
#include <typeinfo>
#include <GateRunManager.hh>
#include "GateSiPM.hh"
#include "G4UIcmdWithADouble.hh"
#include "GateSiPMMessenger.hh"
#include "GateTools.hh"
#include "G4UnitsTable.hh"
#include "NumpyFile.hh"
#include "G4Run.hh"
#include "GateXMLDocument.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "GateRandomEngine.hh"
#include "Randomize.hh"
#include "GateClock.hh"

#include "GateObjectStore.hh"
#include "GateArrayParamsFinder.hh"
#include "GateSystemListManager.hh"
#include "GateVSystem.hh"
#include "stdint.h"
#include "GateVVolume.hh"
#include "GateBox.hh"
#include <math.h>
#include <stdio.h>
#include <random>
#include <gsl/gsl_integration.h>


uint64_t SiPMPulse::currentId = 0;


GateSiPM::GateSiPM(GatePulseProcessorChain* itsChain,
                   const G4String& itsName)
    : GateVPulseProcessor(itsChain,itsName),      
      m_histRes(0), m_startSignal(0), m_durationSignal(0),m_durationPulse(0),
      m_stepSignal(0), m_surface(NULL), m_deadTime(0.01), m_time_primary_value(0), m_t0(0.),
      m_computerSignal(NULL),m_SPTR(0)
{

  G4cout <<  itsName << "!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
  m_file_Pulses.open("pulses.npy");
  m_file_Pulses.register_variable("uid", &m_pulseID);
  m_file_Pulses.register_variable("time", &m_current_time_value);

  m_file_Pulses.register_variable("primary_time", &m_time_primary_value);
  m_file_Pulses.register_variable("value", &m_amplitude);
  m_file_Pulses.register_variable("type", &m_pulseType);
  m_file_Pulses.register_variable("generalDetId", &m_generalDetId);
  m_file_Pulses.register_variable("amplitude", &m_current_pulse_amplitude);
  m_file_Pulses.register_variable("runID", &m_runID);
//    m_file_Pulses.register_variable("eventID", &m_eventID);
  m_file_Pulses.register_variable("trackID", &m_trackID);
  m_file_Pulses.register_variable("nbCrosstalks", &m_nbCrosstalks);
  m_file_Pulses.register_variable("crosstalk_rand", &m_crosstalk_rand);
//  m_file_Pulses.register_variable("seed", &m_seed);
//  m_seed = 0;

  m_file.open("signal.npy");
  m_file.register_variable("time", &m_current_time_value);
  m_file.register_variable("value", &m_current_pulse_value);
  m_file.register_variable("generalDetId", &m_generalDetId);
  m_file.register_variable("runID", &m_runID);
  m_file.register_variable("uid", &m_pulseID);


  m_file_debug.open("debug.npy");
  m_file_debug.register_variable("runID", &m_runID);
  m_file_debug.register_variable("eventID", &m_eventID);
  m_file_debug.register_variable("trackID", &m_trackID);
  m_file_debug.register_variable("time", &m_time);
  m_file_debug.register_variable("generalDetId", &m_generalDetId);
  m_file_debug.register_variable("level1", &m_level1_debug);
  m_file_debug.register_variable("level2", &m_level2_debug);
  m_file_debug.register_variable("energy", &m_energy);
  m_file_debug.register_variable("PDGEncoding", &m_PDGEncoding);
  m_file_debug.register_variable("topVolume", &m_topVolume, 30);
  m_file_debug.register_variable("bottomVolume", &m_bottomVolume, 30);

  m_file_debug.writeHeader();

  if (nVerboseLevel >= 5)
    G4cout << "GateSiPM::GateSiPM" << G4endl;
  G4cout << "!!!!!GateSiPM::GateSiPM #" << "[" << GetObjectName() << "::GateSiPM]" << nVerboseLevel << G4endl;

  m_messenger = new GateSiPMMessenger(this);

}

void GateSiPM::ProcessOnePulse(const GatePulse *inputPulse, GatePulseList &outputPulseList)
{

}


GateSiPM::~GateSiPM()
{

    // intitialize the signal and add the dark noise
   initializeSignal();


    // intitialize afterpulings
   createPobsIntegrand();


  createSignal();
  delete m_messenger;

  //G4cout << "GateSiPM::~GateSiPM #" << "[" << GetObjectName() << "::~GateSiPM]" << nVerboseLevel << G4endl;

  auto system = GateSystemListManager::GetInstance()->GetSystem(0);
  auto level_max = system->GetTreeDepth();

  G4int *array_of_level = new G4int[level_max];

  for (unsigned  i = 1; i < level_max; ++i)
  {
    std::stringstream ss;
    ss << "level" << i;
    m_file.register_variable(ss.str(), &array_of_level[i]);

  }

  m_file.writeHeader();

  for(auto &&signal_in_one_pixel: m_signalTable) {
    auto signal = &signal_in_one_pixel.signal;

    for (auto current_data_signal_iterator = signal->begin();
         current_data_signal_iterator != signal->end(); current_data_signal_iterator++) {

      auto current_data_signal = *current_data_signal_iterator;
      auto current_time_value = current_data_signal.time;

      m_current_time_value = current_time_value + this->getStartSignal();
      m_current_pulse_value = current_data_signal.value;
      m_generalDetId = signal_in_one_pixel.generalDetId;
      //G4cout << current_data_signal.runID << G4endl;
      m_runID=current_data_signal.runID;


      for (unsigned i = 1; i < level_max; ++i)
      {
        array_of_level[i] = signal_in_one_pixel.outputVolumeID[i];
      }

      m_file.fill();
    }
  }

  m_file.close();
  m_file_debug.close();

  delete m_computerSignal;
}


/* permet d'initialiser  le SiPM. On créer un signal vide avec toutes les cellules allumés
à la nanoseconde zero, on donne l'id de la cellule etc ...*/
void GateSiPM::build_nested_volumeID(size_t level_current, GateOutputVolumeID outputVol)
{
  G4cout << "##[GateSiPM::build_nested_volumeID] " << "level_current = " << level_current << "\n";
  auto system = GateSystemListManager::GetInstance()->GetSystem(0);
  auto level_max = system->GetTreeDepth();
  auto nb_element_current_level = system->ComputeNofElementsAtLevel(level_current);

  //  G4cout << "##[GateSiPM::build_nested_volumeID] " << "nb_element_current_level = " << nb_element_current_level << "\n";
  // Le SiPM est au level le plus profond donc tant qu'on y est pas c'est quon parle pas
  // de lui
  if(level_current != level_max)
  {
    for (unsigned i_element_current_level = 0; i_element_current_level < nb_element_current_level; ++i_element_current_level)
    {
      outputVol[level_current] = i_element_current_level;
      this->build_nested_volumeID(level_current + 1, outputVol);
    }
    if (nb_element_current_level == 0)
      this->build_nested_volumeID(level_current + 1, outputVol);
  }
  else

    // sinon on créer le signal dans le pixel avec toutes les propriétés
  {
    //    G4cout << "##[GateSiPM::build_nested_volumeID] " << "outputVol = " << outputVol << "\n";

    // Generate SiPM signal table
    auto signal_in_one_pixel = SignalInOnePixel();
    signal_in_one_pixel.generalDetId = m_generalDetId;
    signal_in_one_pixel.outputVolumeID = outputVol;
    m_signalTable.push_back(signal_in_one_pixel);


    //G4cout << "ICI" << " " << microCells.number << G4endl;
  }
}



void GateSiPM::setVolume(const G4String &volume)
{
  m_volume = volume;

  GateObjectStore* m_store = GateObjectStore::GetInstance();
  if (m_store->FindCreator(volume))
  {
    m_volume = volume;

    //G4cout << "[GateSiPM::setVolume] ok == " <<  m_volume << "\n";

    GateVSystem* system = GateSystemListManager::GetInstance()->GetSystem(0);
    auto tree_depth_of_system = system->GetTreeDepth();
    m_numberOfComponentForLevel.resize(tree_depth_of_system);

    //G4cout << "##[GateSiPM::setVolume] " << "system->GetTreeDepth() = " << tree_depth_of_system << "\n";

    GateOutputVolumeID outputVol(tree_depth_of_system);

    m_generalDetId = 0;
    build_nested_volumeID(0, outputVol);

    //Find the array params
    m_arrayFinder = new GateArrayParamsFinder(m_store->FindCreator(m_volume),
                                              m_nbX, m_nbY, m_nbZ);
    // Find dimentions of the SiPM cells
    setCellsDimentions(m_store->FindCreator(m_volume));


    // How many levels higher than volumeName level ?
    m_numberOfHigherLevels = 0;
    auto creator_tmp = m_store->FindCreator(volume);
    while(creator_tmp->GetMotherList())
    {
      creator_tmp =  creator_tmp->GetMotherList()->GetCreator();
      m_numberOfHigherLevels ++;
    }

    //G4cout<<"Nof Higher level "<<m_numberOfHigherLevels;

    creator_tmp = m_store->FindCreator(volume);
    m_numberOfComponentForLevel[0] = creator_tmp->GetVolumeNumber();
    for (G4int i = 1; i < m_numberOfHigherLevels; i++)
    {

      creator_tmp = creator_tmp->GetMotherList()->GetCreator();
      //        G4cout << "GetVolumeNumber for " << i << " = " << creator_tmp->GetVolumeNumber();
      m_numberOfComponentForLevel[i] = creator_tmp->GetVolumeNumber();
    }


    //G4cout<<"[GateSiPM::setVolume] numberOfComponentForLevel ";
    for(auto &it: m_numberOfComponentForLevel  )
    {
      G4cout << it<<  " ";
    }
    G4cout << "\n";


    //    G4cout << "[GateSiPM::setVolume] numberTotalOfComponentInSystem = " << numberTotalOfComponentInSystem << "\n";
    //    m_signalTable.resize(numberTotalOfComponentInSystem);

  }
  else
  {
    G4cout << "Wrong Volume Name\n";
  }
}


void GateSiPM::FindLevelsParams()
{

}

void GateSiPM::saveAndPurge()
{
  G4cout << "GateSiPM::~GateSiPM #" << "[" << GetObjectName() << "::~GateSiPM]" << nVerboseLevel << G4endl;

  auto system = GateSystemListManager::GetInstance()->GetSystem(0);
  auto level_max = system->GetTreeDepth();

  G4int *array_of_level = new G4int[level_max];

  for (unsigned i = 1; i < level_max; ++i)
  {
    std::stringstream ss;
    ss << "level" << i;
    m_file.register_variable(ss.str(), &array_of_level[i]);
  }

  m_file.writeHeader();

  for(auto &&signal_in_one_pixel: m_signalTable) {
    auto signal = &signal_in_one_pixel.signal;

    for (auto current_data_signal_iterator = signal->begin();
         current_data_signal_iterator != signal->end(); current_data_signal_iterator++) {

      auto current_data_signal = *current_data_signal_iterator;
      auto current_time_value = current_data_signal.time;

      m_current_time_value = current_time_value + this->getStartSignal();
      m_current_pulse_value = current_data_signal.value;
      m_generalDetId = signal_in_one_pixel.generalDetId;

      for (unsigned i = 1; i < level_max; ++i)
      {
        array_of_level[i] = signal_in_one_pixel.outputVolumeID[i];
      }
      m_file.fill();
    }
  }
}

GatePulseList *GateSiPM::ProcessPulseList(const GatePulseList *inputPulseList)
{

  DescribeMyself(1);
  if(nVerboseLevel > 5)
    // DescribeMyself(1);

    if(nVerboseLevel > 3)
    {
      //G4cout << "[GateSiPM::ProcessPulseList]" << "TauRise: " << G4BestUnit(getTauRise(),"Time") << Gateendl;
      //G4cout << "[GateSiPM::ProcessPulseList]" << "TauFall: " << G4BestUnit(getTauFall(),"Time") << Gateendl;
    }



  //G4cout << "[GateSiPM::ProcessPulseList] GATE CURRENT TIME = " << this->GetCurrentTime() << " == " << G4BestUnit(this->GetCurrentTime(),"Time")<< G4endl;
  //G4cout << m_microCells << G4endl;

  size_t n_pulses = inputPulseList->size();

  G4cout << "[GateSiPM::ProcessPulseList]" << ": processing input list with " << n_pulses << " entries" << G4endl;

  if (nVerboseLevel>1)
    G4cout << "[GateSiPM::ProcessPulseList]" << ": processing input list with " << n_pulses << " entries\n";
if(n_pulses>0)
{
      G4int iterator =0;
      for(auto pulse_iterator = inputPulseList->begin(); pulse_iterator != inputPulseList->end(); pulse_iterator++)
      {
        iterator=iterator+1;
        GatePulse* pulse = *pulse_iterator;
        auto temps=pulse->GetTime();
        //G4cout << G4BestUnit(temps, "Time")<< G4endl;

        G4int m_depth =  (size_t)(pulse->GetVolumeID().GetCreatorDepth(m_volume));
        G4cout << "the depth of the interaction is: "<< m_depth  << "and the volumeid is "<< pulse->GetVolumeID() <<G4endl;
        if (m_depth!=-1)
        {
          const GateOutputVolumeID& blockID = pulse->GetOutputVolumeID();
          const GateVolumeID& volumeID = pulse->GetVolumeID();

          if (blockID.IsInvalid())
          {
            G4cout << "\t[GateReadout::ProcessPulseList]: out-of-block hit for \n"
                   <<  *pulse << Gateendl
                   << " -> pulse ignored\n\n";
            continue;
          }


          // G4cout << "\tdepth = " << m_depth << "\t Energy = " << G4BestUnit(pulse->GetEnergy(), "Energy") ;

          m_arrayFinder->FindInputPulseParams(pulse->GetVolumeID().GetCopyNo(m_depth), m_i, m_j, m_k);

          // that depends of the depth of application
          // of the dead time

          pulse->GetOutputVolumeID();

          m_generalDetId = pulse->GetVolumeID().GetCopyNo(m_depth);

          G4cout << "m_generalDetId = " << m_generalDetId << G4endl;

          auto one_pulse = SiPMPulse();
          one_pulse.generalDetId= (G4int) m_generalDetId; // pulse->GetVolumeID().GetCopyNo(m_depth);
          one_pulse.outputVolumeID=(*pulse).GetOutputVolumeID();
          one_pulse.microcell=wichMicroCell(pulse->GetLocalPos());
          one_pulse.time= pulse->GetTime();
          double sptr_noise = G4RandGauss::shoot(0, m_SPTR);
          one_pulse.time=one_pulse.time+sptr_noise;
          one_pulse.runID = pulse->GetRunID();
          one_pulse.trackID = pulse->GetTrackID();
          one_pulse.typePulse=primary;

          G4cout << "\tProcessPulseList:: pulse  = time = " << one_pulse.time
                 << " type = " <<  one_pulse.typePulse << " µcell = " << one_pulse.microcell
                 << "\n";





          auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), one_pulse);
          m_SiPMPulsesTable.insert(it, one_pulse);
          G4cout << "i added a pulse in processpulslist, pulse table is now of size: " << m_SiPMPulsesTable.size() << "\n";
          //m_pulsesTable.insert(m_pulsesTable.begin(), one_pulse);
          //G4cout << "in1" << " " << one_pulse.generalDetId << " " << (*it).generalDetId << " "  << G4endl;



          //pulses_in_pixel_concerned->push_back((one_pulse));
        }
        m_time = pulse->GetTime() / picosecond;



        //      G4cout << GateTools::Indent(1) << "[GateSiPM::ProcessPulseList]" << "pulse time: " << G4BestUnit(pulse_time ,"Time")  << Gateendl;
        //      G4cout << GateTools::Indent(1) << "[GateSiPM::ProcessPulseList]" << "compute pulse = " <<  compute_pulse_tmp << Gateendl;
        m_runID   = GateRunManager::GetRunManager()->GetCurrentRun()->GetRunID();
        m_runID   = pulse->GetRunID();
        m_eventID = pulse->GetEventID();
        m_energy = pulse->GetEnergy();
        m_trackID = pulse->GetTrackID();
        m_PDGEncoding = pulse->GetPDGEncoding();
        m_bottomVolume = m_topVolume = "";

        G4cout << "pulse->GetVolumeID().size() = " << pulse->GetVolumeID().size() << G4endl;


        if(pulse->GetVolumeID().GetTopVolume())
        {
          m_topVolume = pulse->GetVolumeID().GetTopVolume()->GetName();
    //          m_topVolume = "plop";
        }
    //
        if( pulse->GetVolumeID().GetBottomVolume())
        {
          m_bottomVolume = pulse->GetVolumeID().GetBottomVolume()->GetName();
    //          m_bottomVolume = "plop";
        }






        m_file_debug.fill();
        // G4cout << iterator << G4endl;
        // G4cout << temps << G4endl;
      }
    }

         return GateVPulseProcessor::ProcessPulseList(inputPulseList);
}



void GateSiPM::createSignal(){

  m_file_Pulses.writeHeader();
//G4cout << "/////////////////////////in createSignal pulse table size: " << m_SiPMPulsesTable.size() <<G4endl;
  //initialize SIPMs
  auto t_microCells_concerned = &((m_signalTable[0]).cellsList_ns);
  if (t_microCells_concerned->size()<=1)
    initialize();

  while (m_SiPMPulsesTable.size()>0)
    //for (auto pulse= pulses_in_pixel_concerned->begin(); pulse!=pulses_in_pixel_concerned->end(); pulse++ )
  {


    auto& pulse1= m_SiPMPulsesTable.front();
    m_pulse.microcell=pulse1.microcell;
    m_pulse.time= pulse1.time;
    m_pulse.runID = pulse1.runID;
    m_pulse.trackID = pulse1.trackID;
    m_pulse.typePulse=pulse1.typePulse;
    m_pulse.generalDetId=pulse1.generalDetId;
    m_pulse.outputVolumeID=pulse1.outputVolumeID;
    m_outputVolumeID=pulse1.outputVolumeID;
    m_generalDetId=pulse1.generalDetId;
    m_SiPMPulsesTable.erase(m_SiPMPulsesTable.begin());
    m_nbCrosstalks = 0;
    m_crosstalk_rand = -1;


    p_signal_concerned = &(m_signalTable[m_pulse.generalDetId].signal);

    auto t_microCells_concerned = &((m_signalTable[m_pulse.generalDetId]).cellsList_ns);


//      G4cout << "pulse.microcell = " << pulse.microcell << " pulse1.microcell = " << pulse1.microcell << G4endl;
//      G4cout << "pulse.typePulse = " << pulse.typePulse << G4endl;
    assert(m_pulse.microcell < t_microCells_concerned->size());

    auto last_time=(t_microCells_concerned->at(m_pulse.microcell)).time;

//      G4cout << "SET last_time = " << last_time << " pulse.time = " << pulse.time << "\n";

    //verfify if the microCells has begin to charge
    if ( m_pulse.time > last_time + m_deadTime) {


      //verifiy if triggering of the cell ( because of the recovery phase)
      G4double A = 1 - exp(-((m_pulse.time - (last_time + m_deadTime)) / m_tauRecovery_ns));

      if (m_pulse.typePulse != primary || G4UniformRand() < A )
      {

        if (m_pulse.typePulse != crosstalk)
        {
          generateCrosstalk(A);

        }
        if (true)
        {
          generateAfterpulses();
        }


        (t_microCells_concerned->at(m_pulse.microcell)).time = m_pulse.time;
        m_current_pulse_amplitude = writeSignal(last_time);


        if (m_pulse.typePulse == primary || m_pulse.typePulse == darkNoise )
          m_time_primary_value = m_pulse.time;

        m_pulseType = m_pulse.typePulse;
        m_runID = m_pulse.runID;
        m_trackID = m_pulse.trackID;
        m_amplitude=A;
        m_current_time_value = m_pulse.time;
        m_pulseID = m_pulse.m_id;
        m_file_Pulses.fill();
      }
    }
  }

  m_file_Pulses.close();
}


G4double GateSiPM::writeSignal (G4double last_time)
{

  G4double amplitude = 0. ;

  if(  m_pulse.time - (last_time + m_deadTime ) < 0. )
    amplitude = 0.;
  else
    amplitude= m_signalAmplitude*(1-exp(-(m_pulse.time - (last_time + m_deadTime ) )/m_tauRecovery_ns)) ;

  //amplitude = m_signalAmplitude;


//  std::normal_distribution<double> distribution(0,(m_signalAmplitudeSigma));
//  double noise = distribution(m_generator);
  double noise = G4RandGauss::shoot(0, m_signalAmplitudeSigma);
  
//    G4cout << "writeSignal:: amplitude = " << amplitude << " noise = " << noise << " m_signalAmplitudeSigma = " << m_signalAmplitudeSigma << " new amp = " << noise + amplitude << "\n";

  amplitude = noise + amplitude ;

  //find location to write the signal
  GateSiPMDataSignal dump;
  dump.time= m_pulse.time;
  auto it1 = std::lower_bound(p_signal_concerned->begin(), p_signal_concerned->end(), dump);
  dump.time= m_pulse.time +m_durationPulse;
  auto it2 = std::lower_bound(p_signal_concerned->begin(), p_signal_concerned->end(), dump);

  //iterate in signal and write
  for (auto current_data_signal_iterator = it1; current_data_signal_iterator != it2; current_data_signal_iterator++)
  {
    auto compute_pulse_tmp = m_computerSignal->compute(current_data_signal_iterator->time - m_pulse.time) *  amplitude;
    (*current_data_signal_iterator).value += compute_pulse_tmp;
    (*current_data_signal_iterator).runID = m_pulse.runID ;
    (*current_data_signal_iterator).trackID = m_pulse.trackID ;
  }

  return amplitude;
}





G4double GateSiPM::interp1(G4double x_new)
{
  G4double dx, dy, m, b;
  size_t x_max_idx = m_Pulse_time.size() - 1;

//  size_t idx = nearestNeighbourIndex(x, x_new);
  auto it = std::lower_bound(m_Pulse_time.begin(), m_Pulse_time.end(), x_new);
  std::size_t idx = std::distance(m_Pulse_time.begin(), it);

  if (m_Pulse_time[idx] > x_new)
  {
    dx = idx > 0 ? (m_Pulse_time[idx] - m_Pulse_time[idx - 1]) : (m_Pulse_time[idx + 1] - m_Pulse_time[idx]);
    dy = idx > 0 ? (m_Pulse_value[idx] - m_Pulse_value[idx - 1]) : (m_Pulse_value[idx + 1] - m_Pulse_value[idx]);
  }
  else
  {
    dx = idx < x_max_idx ? (m_Pulse_time[idx + 1] - m_Pulse_time[idx]) : (m_Pulse_time[idx] - m_Pulse_time[idx - 1]);
    dy = idx < x_max_idx ? (m_Pulse_value[idx + 1] - m_Pulse_value[idx]) : (m_Pulse_value[idx] - m_Pulse_value[idx - 1]);
  }
  m = dy / dx;
  b = m_Pulse_value[idx] - m_Pulse_time[idx] * m;
  return x_new * m + b;

}


double afterPulse_function(double t, void  *params){
  GateSiPM *sipm = (GateSiPM*) params;
  if (t- sipm->getT0() > 0) {
    G4double ret=sipm->getCap()*pow(t,sipm->geta())*exp(-t/sipm->getTauBuilk())*(1-exp(-(t-sipm->getT0())/sipm->getTauRecovery()));
    return double (ret);
  }
  else
    return 0.;
}

double afterCrosstalk_function(G4double t,void  *params){
  GateSiPM *sipm = (GateSiPM*) params;
  return sipm->getCct()*pow(t,sipm->getb())*exp(- t/sipm->getTauBuilk());
}


void GateSiPM::createPobsIntegrand(){
  G4double res1;
  G4double res2=0;
  G4double t=m_stepSignal;
  G4double resap, resct;
  double  error;

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000);

  gsl_function Ap;
  Ap.function = &afterPulse_function;
  Ap.params=this;

  gsl_function Act;
  Act.function = &afterCrosstalk_function;
  Act.params=this;

  SiPMProbs topush;
  topush.value=0.;
  topush.number=0.;
  m_afterPulse.push_back(topush);

  do
  {
    res1=res2;
    if (t>m_deadTime) {
      gsl_integration_qags(&Ap,m_deadTime,t,0,1e-8,100000, w,&resap,&error);
      gsl_integration_qags(&Act,m_histRes,t,0,1e-8,100000, w,&resct, &error);
    }
    else {
      resct=0;
      resap=0;
    }

    topush.value=resap;
    topush.number=t;
    m_afterPulse.push_back(topush);

    topush.value=resct;
    m_afterCrosstalk.push_back(topush);


    res2=resap+resct;
    t=t+m_stepSignal;
  }
  while (res2-res1 > 10e-8 || t <= (m_t0 + 2* m_stepSignal + m_deadTime ));

  gsl_integration_workspace_free (w);
  //for (auto i=m_afterPulse.begin(); i!=m_afterPulse.end(); i++)
  //G4cout << (*i).value << " " << (*i).number << " panzani" << G4endl;
}


void GateSiPM::generateAfterpulses ( )
{
//  G4cout << "generateAfterpulses started = " << G4endl;
  //unsigned int seed =
  //         static_cast<unsigned int>(*GateRandomEngine::GetInstance()->GetRandomEngine());
  //  srand(seed);

  SiPMProbs pos, posct; // a probaliblity
  G4double interpolated; // value resulting from interpolation
//  G4int* microCell_concernedPos; // position x and y of a microCell
  G4double theta;
  G4int x;
  G4int y;
  G4int i;

//  auto if_afterpulse_rand = std::bind(std::uniform_real_distribution<double>(0.,1.), m_generator);
//  auto theta_rand = std::bind(std::uniform_real_distribution<double>(0.,2 * M_PI), m_generator);



  //***** creating afterpulse
  //pos.value = rand()/(G4double)(RAND_MAX);  // fnum  ]0,1[
  pos.value = G4UniformRand();
  //G4cout << pos.value << " " << m_afterPulse.back().value << G4endl;
  while (pos.value <= m_afterPulse.back().value){

    auto itap= std::lower_bound(m_afterPulse.begin(), m_afterPulse.end(), pos);
    auto time=itap->number;
    SiPMPulse new_pulse;
    new_pulse.microcell=m_pulse.microcell;
    new_pulse.time= m_pulse.time + time ;
    new_pulse.runID = m_pulse.runID;
    new_pulse.trackID = m_pulse.trackID;
    new_pulse.typePulse=afterPulse;
    new_pulse.generalDetId=m_generalDetId;

//        for (auto&& p : m_SiPMPulsesTable)
//        {
//            G4cout << "pulse  = time = " << p.time << "\n";
//        }
//
//        G4cout << "generateAfterpulses::### INSERT !! " <<  new_pulse.time <<  " type = " << new_pulse.typePulse <<   "\n";

    auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), new_pulse);
    m_SiPMPulsesTable.insert(it, new_pulse);

//        for (auto&& p : m_SiPMPulsesTable)
//        {
//            G4cout << "generateAfterpulses:: pulse  = time = " << p.time << " type = " <<  p.typePulse << "\n";
//        }


//        G4cout << "### END INSERT !! " <<  new_pulse.time <<   "\n";

    //        pos.value= pos.value + rand()/(G4double)(RAND_MAX) ;
    pos.value= pos.value + G4UniformRand();
  }

//  G4cout << "generateAfterpulses creating aftercrosstalk = " << G4endl;

  //***** creating aftercrosstalk
  std::vector<std::pair<G4int, G4int> > fired_cells;
  //    posct.value = rand()/(G4double)(RAND_MAX);  // fnum  ]0,1[
  posct.value = G4UniformRand();
  while (posct.value <= m_afterCrosstalk.back().value){
//    G4cout << "generateAfterpulses creating aftercrosstalk  posct.value = " << posct.value << " m_afterCrosstalk.back().value = " <<  m_afterCrosstalk.back().value << G4endl;

    auto itap= std::lower_bound(m_afterCrosstalk.begin(), m_afterCrosstalk.end(), posct);
    auto time=itap->number;

    // calculate microCells to fire for crosstalk
    auto microCell_concernedPos=microCellPos_k_to_i_j(m_pulse.microcell);

    //        pos.value=rand()/(G4double)(RAND_MAX);
    pos.value = G4UniformRand();

    //interpolation for the pixel distance
    if (pos.value > m_CrosstalkDispertion.at(0).value)
    {
      auto bornInf = std::lower_bound(m_CrosstalkDispertion.begin(), m_CrosstalkDispertion.end(), pos);
      bornInf--;
      auto bornSup = std::upper_bound(m_CrosstalkDispertion.begin(), m_CrosstalkDispertion.end(), pos);
      interpolated = bornInf->number + (bornSup->number - bornInf->number)*(pos.value - bornInf->value)/(bornSup->value - bornInf->value);
      interpolated=interpolated+1;
    }
    else
    {
      interpolated = 2;
    }

    // find cell for crosstalk and empiling the new pulse
    //theta= rand()/(G4double)(RAND_MAX) * 2 * M_PI;
    theta =  G4RandFlat::shoot(0., 2*M_PI);
    x = interpolated * cos(theta) + microCell_concernedPos.first;
    y = interpolated * sin(theta) + microCell_concernedPos.second;
//    G4int* position= new G4int [2] {x,y};
    auto position = std::make_pair(x, y);
//    G4cout << " theta = " << theta << " interpolated = " << interpolated << " microCell_concernedPos.first = " << microCell_concernedPos.first << " m_pulse.microcell = " << m_pulse.microcell << "\n";
//    G4cout << " position = " << position.first << ", " << position.second << G4endl;
//    G4cout << " position = " << x << ", " << y << G4endl;

    // if the microcell exists (inside the sipm) then we add a pulse
    if ( (x >= 0) & (y >= 0) && (x < microCells.numberDim1) && (y < microCells.numberDim2) && !(std::find(fired_cells.begin(), fired_cells.end(), position) != fired_cells.end()) )
    {
      SiPMPulse new_pulse;
      new_pulse.microcell=microCellPos_i_j_to_k(std::make_pair(x,y));
      new_pulse.time= m_pulse.time + time;
      new_pulse.runID = m_pulse.runID;
      new_pulse.trackID = m_pulse.trackID;
      new_pulse.typePulse=afterCrosstalk;
      new_pulse.generalDetId=m_generalDetId;
      fired_cells.push_back(position);

      auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), new_pulse);
      m_SiPMPulsesTable.insert(it, new_pulse);

      //            posct.value= posct.value + rand()/(G4double)(RAND_MAX);
      posct.value= posct.value + G4UniformRand();
    }
  }

//  G4cout << "generateAfterpulses ended = " << G4endl;
}


void GateSiPM::generateCrosstalk (double A0)
{

//  G4cout << "generateCrosstalk started = " << G4endl;

  //unsigned int seed =
  //    static_cast<unsigned int>(*GateRandomEngine::GetInstance()->GetRandomEngine());
  //    srand(seed);

//  SiPMPulse new_pulse;
//  new_pulse.microcell=microCellPos_i_j_to_k(new G4int [2] {2,2});
//  new_pulse.time= m_pulse.time;
//  new_pulse.runID = m_pulse.runID;
//  new_pulse.trackID = m_pulse.trackID;
//  new_pulse.typePulse=crosstalk;
//  new_pulse.generalDetId=m_generalDetId;
//  //G4cout << " les microcells  " << new_pulse.microcell << " " << pulse.microcell << G4endl;
//  (m_SiPMPulsesTable).insert((m_SiPMPulsesTable).begin(), new_pulse);
//  return;




  SiPMProbs pos; // a probaliblity
  G4double interpolated; // value resulting from interpolation
//  G4int* microCell_concernedPos; // position x and y of a microCell
  G4double theta;
  G4int x;
  G4int y;
  G4int i;

//  auto crosstalk_rand = std::bind(std::uniform_real_distribution<double>(0.,1.), m_generator);
//  auto theta_rand = std::bind(std::uniform_real_distribution<double>(0.,2 * M_PI), m_generator);

  // calculate the number of crosstalks
  // pos.value = rand()/(G4double)(RAND_MAX);  // fnum  ]0,1[
  m_crosstalk_rand = pos.value = G4UniformRand();  // fnum  ]0,1[

//    double proba_with_amp = 1- (1 - pos.value)*m_amplitude/A0;
//    pos.value = proba_with_amp;


  auto itnbCrosstalk= std::upper_bound(m_CrosstalkTable.begin(), m_CrosstalkTable.end(), pos);
  auto nbCrosstalk=itnbCrosstalk->number;

  m_nbCrosstalks = nbCrosstalk;

  // calculate microCells to fire for crosstalk
  auto microCell_concernedPos=microCellPos_k_to_i_j(m_pulse.microcell);

  std::vector<std::pair<G4int, G4int> > fired_cells;
  i=0;
  while (i < nbCrosstalk)
  {
    //    pos.value=rand()/(G4double)(RAND_MAX);
    pos.value = G4UniformRand();  // fnum  ]0,1[

    //interpolation for the pixel distance

    if (pos.value > m_CrosstalkDispertion.at(0).value)
    {
      auto bornInf = std::lower_bound(m_CrosstalkDispertion.begin(), m_CrosstalkDispertion.end(), pos);
      bornInf--;
      auto bornSup = std::upper_bound(m_CrosstalkDispertion.begin(), m_CrosstalkDispertion.end(), pos);
      interpolated = bornInf->number + (bornSup->number - bornInf->number)*(pos.value - bornInf->value)/(bornSup->value - bornInf->value);
      interpolated=interpolated+1;
    }
    else
    {
      interpolated = 2;
    }
    // find cell for crosstalk and empiling the new pulse
    //        theta= rand()/(G4double)(RAND_MAX) * 2 * M_PI;
    theta = G4RandFlat::shoot(0.  ,2 * M_PI );
    x = interpolated * cos(theta) + microCell_concernedPos.first;
    y = interpolated * sin(theta) + microCell_concernedPos.second;
    auto position= std::make_pair(x,y);

    // if the microcell exists (inside the sipm) then we add a pulse
    if ( (x >= 0) && (y >= 0) && (x < microCells.numberDim1) && (y < microCells.numberDim2) && !(std::find(fired_cells.begin(), fired_cells.end(), position) != fired_cells.end()))
    {
      SiPMPulse new_pulse;
      new_pulse.microcell=microCellPos_i_j_to_k( std::make_pair(x,y));
      double sptr_noise = G4RandGauss::shoot(0, m_SPTR*sqrt(i+2));
      while ((m_SPTR*2.355*sqrt(i+1) + sptr_noise) < 0){
        sptr_noise = G4RandGauss::shoot(0, m_SPTR*sqrt(i+2));}
      new_pulse.time= m_pulse.time+abs(m_SPTR*2.355*sqrt(i+1)) + sptr_noise ;
      new_pulse.runID = m_pulse.runID;
      new_pulse.trackID = m_pulse.trackID;
      new_pulse.typePulse=crosstalk;
      new_pulse.generalDetId=m_generalDetId;
      //G4cout << " les microcells  " << new_pulse.microcell << " " << pulse.microcell << G4endl;
      (m_SiPMPulsesTable).insert((m_SiPMPulsesTable).begin(), new_pulse);

//            G4cout << "generateCrosstalk::### INSERT !! " <<  new_pulse.time <<  " type = " << new_pulse.typePulse <<   "\n";

      fired_cells.push_back(position);
      i++;
    }
  }

//  G4cout << "generateCrosstalk ended = " << G4endl;

}



void GateSiPM::initializeSignal()
{
    G4double durationSignal = this->getDurationSignal();
  GateSiPMCell tocreate; // micro_cells list ( for creation)

  //    unsigned int seed =static_cast<unsigned int>(*GateRandomEngine::GetInstance()->GetRandomEngine());
  //srand(seed);

  //std::default_random_engine generator(time(NULL));


//  auto whiteNoise_rand = std::bind(
//      std::normal_distribution<double> (0,m_whiteNoiseSigma), m_generator
//  );

//  auto deltatime_rand = std::bind(
//      std::uniform_real_distribution<double>(0,
//                                             this->getStepSignal()),
//      m_generator);

//  auto if_darknoise_rand = std::bind(std::uniform_real_distribution<double>(0.,1.), m_generator);

//  auto which_microcell_rand = std::bind(std::uniform_int_distribution<G4int>(0,microCells.number - 1), m_generator);


  //    std::std::uniform_int_distribution<int>(0, microCells.number);


  for(size_t i = 0; i < m_signalTable.size(); i++)
  {
    G4double prob;

    //Get the good sipm concerned
    p_signal_concerned = &((m_signalTable[i]).signal);

    m_signalTable[i].generalDetId=i;

    G4double startSignal = -m_durationPulse -3 * m_tauRecovery_ns;

    //creat microcells in SiPM
    auto t_microCells_concerned = &((m_signalTable[i]).cellsList_ns);
    if (t_microCells_concerned->size()<=1)
    {
      tocreate.time=startSignal;
      *t_microCells_concerned= *new std::vector<GateSiPMCell> (microCells.number,tocreate);
    }

    //auto signal = &signal_in_one_pixel.signal;

    // if first time of signal we start a bit before in case of darknoise


    //    signal_in_one_pixel.generalDetId = m_generalDetId;
    //    G4cout << "reserve because size = " << signal->size() << "\n";

//    SiPMPulse one_pulse;
//    //one_pulse.microcell= G4int((rand()/(G4double)(RAND_MAX))*microCells.number);
//    one_pulse.microcell =  G4RandFlat::shootInt( (long)0,  microCells.number - 1);
//    one_pulse.time= startSignal + 1;
//    one_pulse.runID = -1;
//    one_pulse.trackID = -1;
//    one_pulse.typePulse=darkNoise;
//    one_pulse.generalDetId=i;
//    auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), one_pulse);
//    m_SiPMPulsesTable.insert(it, one_pulse);

    if (p_signal_concerned->size() == 0)
    {


      uint64_t reserve = (uint64_t) (durationSignal-startSignal) / this->getStepSignal();
      //      if(nVerboseLevel > 1)
      //        G4cout << "[GateSiPM::ProcessPulseList]" << "Allocate signal, I reserve " << reserve << " elements" << Gateendl;
      p_signal_concerned->reserve(reserve);
      for (G4double current_time_value = startSignal; current_time_value < this->getDurationSignal(); current_time_value += this->getStepSignal()   )
      {
        prob = G4UniformRand();  // fnum  ]0,1[

        if (prob <= m_darkNoiseProb)
        {
//                    G4cout << prob << " " <<  m_darkNoiseProb << " " << m_darkNoise << G4endl;
          SiPMPulse one_pulse;
          //                    one_pulse.microcell= G4int((rand()/(G4double)(RAND_MAX))*microCells.number);
          one_pulse.microcell = G4RandFlat::shootInt( (long)0,  microCells.number - 1);
          //                    one_pulse.time= current_time_value;

          one_pulse.time = current_time_value + G4RandFlat::shoot(0., this->getStepSignal());

//                    G4cout << "dark noise, this->getStepSignal() = " << this->getStepSignal() << " current_time_value = " << current_time_value << " one_pulse.time = " << one_pulse.time << "\n";
          one_pulse.runID = -1;
          one_pulse.trackID = -1;
          one_pulse.typePulse=darkNoise;
          one_pulse.generalDetId=i;
          auto it = std::lower_bound(m_SiPMPulsesTable.begin(), m_SiPMPulsesTable.end(), one_pulse);
          if(m_SiPMPulsesTable.size() > 0) //iteraror can be deferenced
          {
            one_pulse.runID  = it->runID;
          }

          m_SiPMPulsesTable.insert(it, one_pulse);
        }
        // write white signal
        GateSiPMDataSignal d;
        d.time = current_time_value;
        d.value = G4RandGauss::shoot(0, m_whiteNoiseSigma);
//                d.value = 0.;
        d.runID = 0;
        p_signal_concerned->push_back(d);
      }
    }
  }

}

void GateSiPM::initialize()
{

//  m_generator.seed(time(NULL));
//  m_generator.seed(static_cast<unsigned int>(*GateRandomEngine::GetInstance()->GetRandomEngine()));
//  m_generator.seed(GateRandomEngine::GetInstance()->GetRandomEngine()->getSeed());
//  m_seed = GateRandomEngine::GetInstance()->GetRandomEngine()->getSeed();
//  m_engine = GateRandomEngine::GetInstance()->GetRandomEngine();

  //G4cout << m_signalAmplitudeSigma << " " << m_signalAmplitude << G4endl;
  //inititialize microCells dimentions
  setMircoCellsProperties();

  //initialize dar knoise probability
  m_darkNoiseProb =  (m_stepSignal/1e9) * m_darkNoise;





}

// Attaches a G4MaterialPropertiesTable to the optical surface.

void GateSiPM::DescribeMyself(size_t indent)
{
  G4cout << GateTools::Indent(indent) << "GateSiPM " << Gateendl;
  G4cout << GateTools::Indent(indent) << "TauRise: " << G4BestUnit(getTauRise(),"Time") << Gateendl;
  G4cout << GateTools::Indent(indent) << "TauFall: " << G4BestUnit(getTauFall(),"Time") << Gateendl;
  G4cout << GateTools::Indent(indent) << "Start signal: " << G4BestUnit(getStartSignal(),"Time") << Gateendl;
  G4cout << GateTools::Indent(indent) << "Duration signal: " << G4BestUnit(getDurationSignal(),"Time") << Gateendl;
  G4cout << GateTools::Indent(indent) << "Step signal: " << G4BestUnit(getStepSignal(),"Time") << Gateendl;
}

void  GateSiPM::setMircoCellsProperties()
{
  microCells.numberDim1= m_mircoCellsPitch[m_surface[0]]/ microCells.lengthDim1;
  microCells.numberDim2= m_mircoCellsPitch[m_surface[1]]/ microCells.lengthDim2;
  microCells.number= (G4int) microCells.numberDim1 * microCells.numberDim2;
}


void GateSiPM::setCellsDimentions(GateVVolume* volume)
{
  GateBox* box = static_cast<GateBox*>(volume);
  m_mircoCellsPitch=new G4double[3] {box->GetBoxXLength(), box->GetBoxYLength(), box->GetBoxZLength()};
}

void GateSiPM::setMircoCellsDim(G4double microCellsDim1,G4double microCellsDim2 )
{
  microCells.lengthDim1= microCellsDim1;
  //G4cout<< G4BestUnit(microCellsDim1,"Length") << G4endl;
  microCells.lengthDim2= microCellsDim2;
}

G4int GateSiPM::wichMicroCell(G4ThreeVector localPos) const
{
  G4double dump [3]={localPos.getX(), localPos.getY(),localPos.getZ()};
  G4int i=  (dump[m_surface[0]] + (m_mircoCellsPitch[m_surface[0]]/2)) / (microCells.lengthDim1);
  G4int j=  (dump[m_surface[1]] + (m_mircoCellsPitch[m_surface[1]] /2) ) / (microCells.lengthDim2);
  G4cout << "wichMicroCell::i = " << i << " j = " << j << " microCells.numberDim2 = " << microCells.numberDim2 << " \n";
  return i*microCells.numberDim2 + j;
}

G4int GateSiPM::microCellPos_i_j_to_k (const std::pair<G4int, G4int> &pos ){
  return pos.first*microCells.numberDim2 + pos.second;
}

std::pair<G4int, G4int> GateSiPM::microCellPos_k_to_i_j(G4int pos)
{
  G4int i= pos/microCells.numberDim2;
  G4int j=pos % microCells.numberDim2;

//  G4cout << " pos = " << pos << " i = " << i << " j = " << j << G4endl;

//  return  new G4int [2] {i,j};
  return std::make_pair(i, j);
}

G4double GateSiPM::getTauRise() const
{
  return m_tauRise_ns;
}

void GateSiPM::setTauRise(G4double tauRise)
{
  m_tauRise_ns = tauRise / nanosecond;
}

G4double GateSiPM::getTauFall() const
{
  return m_tauFall_ns;
}

void GateSiPM::setTauFall(G4double tauFall)
{
  m_tauFall_ns = tauFall / nanosecond;
}

G4double GateSiPM::getStartSignal() const
{
  return m_startSignal;
}

void GateSiPM::setStartSignal(G4double startSignal)
{
  m_startSignal = startSignal;
}

G4double GateSiPM::getDurationSignal() const
{
  return m_durationSignal;
}

void GateSiPM::setDurationSignal(G4double endSignal)
{
  m_durationSignal = endSignal;
}


G4double GateSiPM::getDurationPulse() const
{
  return m_durationPulse;
}

void GateSiPM::setDurationPulse(G4double duration)
{
  m_durationPulse = duration;
}


G4double GateSiPM::getSPTR() const
{
  return m_SPTR;
}

void GateSiPM::setSPTR(G4double SPTR)
{
  m_SPTR= SPTR;
}


void GateSiPM::setStepSignal(G4double stepSignal)
{
  m_stepSignal = stepSignal;
}


const G4String &GateSiPM::getVolume() const
{
  return m_volume;
}

const G4int* GateSiPM::getSurface() const
{
  return m_surface;
}

G4double* GateSiPM::getCellsDimentions() const
{
  return m_mircoCellsPitch;
}

G4double GateSiPM::getStepSignal() const
{
  return m_stepSignal;
}
G4double GateSiPM::getTauRecovery() const {
  return m_tauRecovery_ns;
}

void GateSiPM::setTauRecovery(G4double tauRecovery){
  m_tauRecovery_ns=tauRecovery;
}


void GateSiPM::setSurface(G4String surface){

  if ( surface=="ZY"){
    m_surface=new G4int [2] {2,1};
  }
  else if (surface=="XZ"){
    m_surface=new G4int [2] {0,2};
  }
  else if (surface=="ZX"){
    m_surface=new G4int [2] {2,0};
  }
  else if (surface=="YZ"){
    m_surface=new G4int [2] {1,2};
  }
  else if (surface=="XY"){
    m_surface=new G4int [2] {0,1};
  }
  else  {
    m_surface=new G4int [2] {1,0};
  }
}


void GateSiPM::setSiPMFromXml (G4String type){
  GateXMLDocument* doc = new GateXMLDocument("./SiPM.xml");
  if (doc->Ok())
  {
    doc->Enter();
    if (doc->Find("sipm",type))
    {
      /*G4String tauFall_ns = doc->GetProperty("tauFall_ns");
      this->setTauFall(G4UIcmdWithADouble::GetNewDoubleValue(tauFall_ns.c_str()));

      G4String tauRise_ns = doc->GetProperty("tauRise_ns");
      this->setTauRise(G4UIcmdWithADouble::GetNewDoubleValue(tauRise_ns.c_str()));

      G4String tauRecovery_ns = doc->GetProperty("tauRecovery_ns");
      this->setTauRecovery(G4UIcmdWithADouble::GetNewDoubleValue(tauRecovery_ns.c_str()));

      G4String microCellsDim1_micron = doc->GetProperty("microCellsDim1_micron");
      G4String microCellsDim2_micron = doc->GetProperty("microCellsDim2_micron");
      this->setMircoCellsDim(G4UIcmdWithADouble::GetNewDoubleValue(microCellsDim1_micron .c_str()), G4UIcmdWithADouble::GetNewDoubleValue(microCellsDim2_micron .c_str()));
*/
      doc->Enter();
      doc->First();
      if (doc->Find("propertiestable"))
      {   G4double value;
        G4double unit;
        doc->Enter();
        while (doc->Next())
        {
          if (doc->GetName() == "property")
          {
            G4String property = doc->GetProperty("name");
            G4String valuestr = doc->GetProperty("value");
            value = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
            G4String unitstr = "1 " + doc->GetProperty("unit");
            if (property == "tauFall"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());

              if(!m_computerSignal)
                m_computerSignal = new ComputeSignalFormula();
              static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauFall_ns = value*unit;

              //this->setTauFall(value*unit);
            }
            if (property == "histRes"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              m_histRes=(value*unit);
            }
            if (property == "tauRise"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              //this->setTauRise(value*unit);
              if(!m_computerSignal)
                m_computerSignal = new ComputeSignalFormula();
              static_cast<ComputeSignalFormula*>(m_computerSignal)->m_tauRise_ns = value*unit;
            }
            if (property == "durationPulse"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              this->setDurationPulse(value*unit);
            }
            if (property == "whiteNoiseSigma"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              m_whiteNoiseSigma=value*unit*10e5;
            }
            if (property == "tauRecovery"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              this->setTauRecovery(value*unit);
            }
            if (property == "SPTR"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              this->setSPTR(value*unit);
            }
            if (property == "signalDeconvolvedAmplitude"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              m_signalAmplitude=value*unit*10e5;
            }
            if (property == "signalDeconvolvedAmplitudeSigma"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              m_signalAmplitudeSigma=value*unit*10e5;
            }
            if (property == "darkNoise"){
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              m_darkNoise=value*unit*10e8;
            }
            if (property == "tauBuilk"){
              m_tauBuilk=(value*unit);
            }
            if (property == "Cap")
              m_Cap=value;
            if (property == "t0")
            {
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              m_t0=value*unit;
            }
            if (property == "deadTime")
            {
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              m_deadTime=value*unit;
            }
            if (property == "Cct")
              m_Cct=value;
            if (property == "a")
              m_a=value;
            if (property == "b")
              m_b=value;
          }
          else if (doc->GetName() == "propertyvector")
          {
            G4String property = doc->GetProperty("name");
            if (property == "CROSSTALK")
            {
              // read vector
              G4double accumulate=0;
              doc->Enter();
              G4int iterator=0;
              SiPMProbs topush;
              while (doc->Find("ve"))
              {
                G4String valuestr  = doc->GetProperty("value");
                G4double value     = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                accumulate+=value;
                topush.value=accumulate;
                topush.number=iterator;
                m_CrosstalkTable.push_back(topush);
                iterator ++;
              }
            }

            else if (property == "PULSE")
            {
              // read vector
              G4String unitstr = "1 " + doc->GetProperty("unit");
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              doc->Enter();

              m_computerSignal = new ComputeSignalVector();

              while (doc->Find("ve"))
              {
                G4String valuestr  = doc->GetProperty("time");
                G4double time = (G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str()))*unit;
                valuestr  = doc->GetProperty("value");
                G4double value = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                //m_Pulse_time.push_back(time);
                //m_Pulse_value.push_back(value);
                static_cast<ComputeSignalVector*>(m_computerSignal)->add_pair(time, value);
                G4cout << value << " " << time<< G4endl;
              }
            }
            else if (property == "DIMENTIONS")
            {
              G4String unitstr = "1 " + doc->GetProperty("unit");
              unit = G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(unitstr.c_str());
              doc->Enter();
              doc->Find("ve");
              G4String valuestr  = doc->GetProperty("value");
              G4double value1    = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
              valuestr  = doc->GetProperty("value");
              G4double value2 = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
              this->setMircoCellsDim(value1*unit,value2*unit);
            }
            else if (property == "CROSSTALK_DISPERTION")
            {
              // read vector
              G4double accumulate=0;
              doc->Enter();
              G4int iterator=0;
              SiPMProbs topush;
              while (doc->Find("ve"))
              {
                G4String valuestr  = doc->GetProperty("value");
                G4double value     = G4UIcmdWithADouble::GetNewDoubleValue(valuestr.c_str());
                accumulate+=value;
                topush.value=accumulate;
                topush.number=iterator;
                m_CrosstalkDispertion.push_back(topush);
                iterator ++;
              }
            }
            doc->Leave();
          }
        }
        doc->Leave();
      }
    }
  }
}


