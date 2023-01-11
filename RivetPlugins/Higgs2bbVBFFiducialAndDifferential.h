#ifndef _TRUTHRIVETTOOLS_HIGGS2BBVBFFIDUCIALANDDIFFERENTIAL_H_
#define _TRUTHRIVETTOOLS_HIGGS2BBVBFFIDUCIALANDDIFFERENTIAL_H_

#include "Rivet/Analysis.hh"

#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Projections/FinalState.hh" 
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet
{
    class Higgs2bbVBFFiducialAndDifferential : public Analysis 
    {
    public:
    //---ctors---
        Higgs2bbVBFFiducialAndDifferential():
            Analysis("Higgs2bbVBFFiducialAndDifferential"),
            _sumW(0.)
            {};

    //---dtor---
        ~Higgs2bbVBFFiducialAndDifferential() {};

    //---utils--
        void init();
        void analyze(const Event& event);
        void finalize();

    private:
        double     _sumW;
        Histo1DPtr _h_pt_h;
        Histo1DPtr _h_mjj;
    };

    DECLARE_RIVET_PLUGIN(Higgs2bbVBFFiducialAndDifferential);
}

#endif
