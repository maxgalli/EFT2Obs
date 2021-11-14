#ifndef _TRUTHRIVETTOOLS_HIGGS2GGFIDUCIALANDDIFFERENTIAL_H_
#define _TRUTHRIVETTOOLS_HIGGS2GGFIDUCIALANDDIFFERENTIAL_H_

#include "Rivet/Analysis.hh"

#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Projections/FinalState.hh" 
#include "Rivet/Projections/PromptFinalState.hh" 
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet
{
    class Higgs2GGFiducialAndDifferential : public Analysis 
    {
    public:
    //---ctors---
        Higgs2GGFiducialAndDifferential():
            Analysis("Higgs2GGFiducialAndDifferential"),
            sumW_(0.)
            {};

    //---dtor---
        ~Higgs2GGFiducialAndDifferential() {};

    //---utils--
        void init();
        void analyze(const Event& event);
        void finalize();

    private:
        double     sumW_;
        Histo1DPtr h_m_gg_;
        Histo1DPtr h_pt_gg_;
        Histo1DPtr h_eta_gg_;
        Histo1DPtr h_pt_g1_;
        Histo1DPtr h_pt_g2_;
        Histo1DPtr h_eta_g1_;
        Histo1DPtr h_eta_g2_;
        Histo1DPtr h_njets_eta2p5_;
        Histo1DPtr h_jet_pt_eta2p5_;
        Histo1DPtr h_njets_eta4p7_;
        Histo1DPtr h_jet_pt_eta4p7_;
    };

    DECLARE_RIVET_PLUGIN(Higgs2GGFiducialAndDifferential);
}

#endif
