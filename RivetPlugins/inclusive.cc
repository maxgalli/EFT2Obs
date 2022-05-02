// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class inclusive : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(inclusive);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      book(_h["XXXX"], "myh1", 1, 0.0, 1.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      _h["XXXX"]->fill(0);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sf = crossSection() / picobarn / sumOfWeights();
      scale(_h["XXXX"], sf); // norm to cross section
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(inclusive);


}
