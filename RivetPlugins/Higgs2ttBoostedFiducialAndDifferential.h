#ifndef _TRUTHRIVETTOOLS_HIGGS2TTBOOSTEDFIDUCIALANDDIFFERENTIAL_H_
#define _TRUTHRIVETTOOLS_HIGGS2TTBOOSTEDFIDUCIALANDDIFFERENTIAL_H_

#include "Rivet/Analysis.hh"

#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {
class Higgs2ttBoostedFiducialAndDifferential : public Analysis {
  public:
    //---ctors---
    Higgs2ttBoostedFiducialAndDifferential()
        : Analysis("Higgs2ttBoostedFiducialAndDifferential"), sumW_(0.) {};

    //---dtor---
    ~Higgs2ttBoostedFiducialAndDifferential() {};

    //---utils--
    void init();
    void analyze(const Event &event);
    void finalize();

  private:
    double sumW_;
    Histo1DPtr h_pt_h_;
    Histo1DPtr _h_sigma;
};

DECLARE_RIVET_PLUGIN(Higgs2ttBoostedFiducialAndDifferential);
}

#endif
