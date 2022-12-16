//
// Created by mdupont on 11/10/17.
//

#ifndef GATE_GATESIPMMESSENGER_HH
#define GATE_GATESIPMMESSENGER_HH



#include "GatePulseProcessorMessenger.hh"


class GateSiPM;


class GateSiPMMessenger : public GatePulseProcessorMessenger
{
public:
  GateSiPMMessenger(GateSiPM *itsSiPM);

  void SetNewValue(G4UIcommand *icommand, G4String string) override;

  inline GateSiPM* GetSiPM(){ return (GateSiPM*) GetPulseProcessor(); }

  virtual ~GateSiPMMessenger();


private:

  G4UIcmdWithADoubleAndUnit *m_tauRiseCmd;
  G4UIcmdWithADoubleAndUnit *m_tauFallCmd;
  G4UIcmdWithAString *m_type;
  G4UIcmdWithoutParameter *m_initialize;
  G4UIcmdWithADoubleAndUnit *m_startSignalCmd;
  G4UIcmdWithADoubleAndUnit *m_durationSignalCmd;
  G4UIcmdWithADoubleAndUnit *m_stepSignalCmd;
  G4UIcmdWithAString *m_volumeCmd;
  G4UIcmdWithAString *m_surface;
};


#endif //GATE_GATESIPMMESSENGER_HH
