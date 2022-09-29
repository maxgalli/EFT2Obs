#ifndef _TRUTHRIVETTOOLS_HIGGS2BBFIDUCIALANDDIFFERENTIAL_H_
#define _TRUTHRIVETTOOLS_HIGGS2BBFIDUCIALANDDIFFERENTIAL_H_

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
    class Higgs2bbFiducialAndDifferential : public Analysis 
    {
    public:
    //---ctors---
        Higgs2bbFiducialAndDifferential():
            Analysis("Higgs2bbFiducialAndDifferential"),
            _sumW(0.)
            {};

    //---dtor---
        ~Higgs2bbFiducialAndDifferential() {};

    //---utils--
        void init();
        void analyze(const Event& event);
        void finalize();

    private:
        double     _sumW;
        Histo1DPtr _h_pt_h;
    };

    DECLARE_RIVET_PLUGIN(Higgs2bbFiducialAndDifferential);
}

#endif
