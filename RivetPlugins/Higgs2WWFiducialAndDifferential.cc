#include "Higgs2WWFiducialAndDifferential.h"

namespace Rivet
{
    void Higgs2WWFiducialAndDifferential::init()
    {
        const double lepConeSize = 0.1;
        const double lepMaxEta = 2.5;
        const Cut lepton_cut = (Cuts::abseta < lepMaxEta);

        // Initialise and register projections
        FinalState fs((Cuts::etaIn(-2.5,2.5)));
        FinalState fsm((Cuts::etaIn(-5,5)));
        declare(fs, "FS");
        declare(fsm, "FSM");

        ChargedLeptons charged_leptons(fs);
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);

        PromptFinalState prompt_leptons(charged_leptons);
        prompt_leptons.acceptMuonDecays(true);
        prompt_leptons.acceptTauDecays(false);

        PromptFinalState prompt_photons(photons);
        prompt_photons.acceptMuonDecays(true);
        prompt_photons.acceptTauDecays(false);

        DressedLeptons dressed_leptons = DressedLeptons(prompt_photons, prompt_leptons, lepConeSize, lepton_cut, true);
        declare(dressed_leptons, "DressedLeptons");

        //---Jets
        FastJets fs_jets(fs, FastJets::ANTIKT, 0.4);
        declare(fs_jets, "JETS");

        MissingMomentum Met(fsm);
        declare(Met, "MET");
        
        //---Histograms
        book(h_pt_h_, "pt_h",{0,30,45,80,120,200,10000});
        book(_h_njets, "njets", {-0.5,0.5,1.5,2.5,3.5,100.5});
    }

    void Higgs2WWFiducialAndDifferential::analyze(const Event& event)
    {               
        const double weight = 1.0;
        
        // jets
        auto jets = apply<JetAlg>(event, "JETS").jetsByPt(Cuts::abseta < 4.7 && Cuts::pt > 30*GeV);

        Particles leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").particlesByPt(10.0*GeV);
        if (leptons.size() < 2) vetoEvent;
        if (leptons[0].pT() < 25*GeV || leptons[1].pT() < 13*GeV) vetoEvent;
        if (leptons[0].charge() == leptons[1].charge()) vetoEvent;
        if (leptons[0].abspid() == leptons[1].abspid()) vetoEvent;

        FourMomentum LL = (leptons[0].momentum() + leptons[1].momentum());
        if (LL.mass() < 12*GeV) vetoEvent;
        if (LL.pT() < 30*GeV) vetoEvent;

        FourMomentum EtMiss = applyProjection<MissingMomentum>(event,"MET").missingMomentum();
        FourMomentum P4H = LL + EtMiss;

        double dphi = deltaPhi(LL, EtMiss);

        double mT = sqrt(2*LL.pT()*EtMiss.pT()*(1-cos(dphi)));
        if (mT < 60*GeV) vetoEvent;

        //double mT_lead = sqrt(leptons[0].mass()*leptons[0].mass() + leptons[0].px()*leptons[0].px() + leptons[0].py()*leptons[0].py());
        double mT_trail = sqrt(leptons[1].mass()*leptons[1].mass() + leptons[1].px()*leptons[1].px() + leptons[1].py()*leptons[1].py());
        if (mT_trail < 30*GeV);

        h_pt_h_->fill(P4H.pT()/GeV, weight);
        _h_njets->fill(jets.size());

        sumW_ += event.weights()[0];
        
        return;
    }

    void Higgs2WWFiducialAndDifferential::finalize()
    {
        return;
    }
}
