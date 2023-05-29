#ifndef _TRUTHRIVETTOOLS_HIGGS2WWFIDUCIALANDDIFFERENTIAL_H_
#define _TRUTHRIVETTOOLS_HIGGS2WWFIDUCIALANDDIFFERENTIAL_H_
/*
 * This is based on both https://cds.cern.ch/record/2691268?ln=it 
 * and https://rivet.hepforge.org/analyses/CMS_2017_I1467451.html
 */

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
    class Higgs2WWFiducialAndDifferential : public Analysis 
    {
    public:
    //---ctors---
        Higgs2WWFiducialAndDifferential():
            Analysis("Higgs2WWFiducialAndDifferential"),
            sumW_(0.)
            {};

    //---dtor---
        ~Higgs2WWFiducialAndDifferential() {};

    //---utils--
        void init();
        void analyze(const Event& event);
        void finalize();

    private:
        double     sumW_;
        Histo1DPtr h_pt_h_;
        Histo1DPtr _h_njets;
        Histo1DPtr _h_sigma;
    };

    DECLARE_RIVET_PLUGIN(Higgs2WWFiducialAndDifferential);
}

#endif
