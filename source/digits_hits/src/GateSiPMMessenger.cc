//
// Created by mdupont on 11/10/17.
//

#include "GateConfiguration.h"


#include "GateSiPMMessenger.hh"
#include "GateSiPM.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "NumpyFile.hh"
#include "stdint.h"


GateSiPMMessenger::GateSiPMMessenger(GateSiPM *itsSiPM)
: GatePulseProcessorMessenger(itsSiPM)
{

  G4String guidance;
  G4String cmdName;

  cmdName = GetDirectoryName() + "setTauRise";
  m_tauRiseCmd= new G4UIcmdWithADoubleAndUnit(cmdName,this);
  m_tauRiseCmd->SetGuidance("Set Time Rise (ns)");
  m_tauRiseCmd->SetUnitCategory("Time");

  cmdName = GetDirectoryName() + "type";
  m_type= new G4UIcmdWithAString(cmdName,this);
  m_type->SetGuidance("Set Type of SiPM");

  cmdName = GetDirectoryName() + "surface";
  m_surface= new G4UIcmdWithAString(cmdName,this);
  m_surface->SetGuidance("Set surface");

  cmdName = GetDirectoryName() + "setTauFall";
  m_tauFallCmd= new G4UIcmdWithADoubleAndUnit(cmdName,this);
  m_tauFallCmd->SetGuidance("Set Time Fall (ns)");
  m_tauFallCmd->SetUnitCategory("Time");

  cmdName = GetDirectoryName() + "setStartSignal";
  m_startSignalCmd= new G4UIcmdWithADoubleAndUnit(cmdName,this);
  m_startSignalCmd->SetUnitCategory("Time");

  cmdName = GetDirectoryName() + "setDurationSignal";
  m_durationSignalCmd= new G4UIcmdWithADoubleAndUnit(cmdName,this);
  m_durationSignalCmd->SetUnitCategory("Time");

  cmdName = GetDirectoryName() + "setStepSignal";
  m_stepSignalCmd= new G4UIcmdWithADoubleAndUnit(cmdName,this);
  m_stepSignalCmd->SetUnitCategory("Time");

  cmdName = GetDirectoryName() + "setVolume";
  m_volumeCmd = new G4UIcmdWithAString(cmdName,this);

  cmdName = GetDirectoryName() + "initialize";
  m_initialize = new G4UIcmdWithoutParameter(cmdName,this);
}

void GateSiPMMessenger::SetNewValue(G4UIcommand *icommand, G4String newValue)
{
  if (icommand== m_tauRiseCmd)
   GetSiPM()->setTauRise(m_tauRiseCmd->GetNewDoubleValue(newValue));
  else if (icommand== m_tauFallCmd) {
   GetSiPM()->setTauFall(m_tauFallCmd->GetNewDoubleValue(newValue));}
  else if (icommand== m_type) {
   GetSiPM()->setSiPMFromXml(newValue);}
  else if (icommand== m_surface) {
   GetSiPM()->setSurface(newValue);}
  else if (icommand== m_startSignalCmd)
    GetSiPM()->setStartSignal(m_startSignalCmd->GetNewDoubleValue(newValue));
  else if (icommand== m_durationSignalCmd)
    GetSiPM()->setDurationSignal(m_durationSignalCmd->GetNewDoubleValue(newValue));
  else if (icommand== m_stepSignalCmd)
    GetSiPM()->setStepSignal(m_stepSignalCmd->GetNewDoubleValue(newValue));
  else if (icommand == m_volumeCmd )
    GetSiPM()->setVolume(newValue);
  else if (icommand == m_initialize )
    GetSiPM()->initialize();
  else
    GatePulseProcessorMessenger::SetNewValue(icommand, newValue);
}

GateSiPMMessenger::~GateSiPMMessenger()
{
  delete m_tauRiseCmd;

}

