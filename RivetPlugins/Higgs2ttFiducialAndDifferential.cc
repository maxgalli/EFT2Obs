#include "Higgs2ttFiducialAndDifferential.h"

namespace Rivet {
void Higgs2ttFiducialAndDifferential::init() {
    TauFinder tauleptonic(TauFinder::DecayMode::LEPTONIC);
    declare(tauleptonic, "TauLeptonic");
    TauFinder tauhadronic(TauFinder::DecayMode::HADRONIC);
    declare(tauhadronic, "TauHadronic");

    FinalState fs((Cuts::etaIn(-10,10)));
    FinalState fsm((Cuts::etaIn(-10,10)));
    declare(fs, "FS");
    declare(fsm, "FSM");
    MissingMomentum Met(fsm);
    declare(Met, "MET");

    declare(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");

    //---Histograms
    book(h_pt_h_, "pt_h", { 0, 45, 80, 120, 140, 170, 200, 350, 450, 10000 });
    book(_h_njets, "njets", {-0.5,0.5,1.5,2.5,3.5,100.5});
    book(_h_jet_pt, "pt_j0",{30,60,120,200,350,10000});
    book(_h_sigma, "h_sigma", 1, 0, 100000000);
}

void Higgs2ttFiducialAndDifferential::analyze(const Event &event) {
    const double weight = 1.0;
    sumW_ += event.weights()[0];
    _h_sigma->fill(sumW_);

    Particles tauhad =
        apply<TauFinder>(event, "TauHadronic").particlesByPt();
    Particles taulep =
        apply<TauFinder>(event, "TauLeptonic").particlesByPt();
    Particles tauhad_fh;
    Particles taulep_fh;

    for(Particle& t : tauhad) {
        if (t.hasParent(25)) tauhad_fh.push_back(t);
    }
    for(Particle& t : taulep) {
        if (t.hasParent(25)) taulep_fh.push_back(t);
    }

    if (tauhad_fh.size() == 0 && taulep_fh.size() == 0)
        vetoEvent;

    //std::cout << tauhad_fh.size();
    //std::cout << taulep_fh.size();

    FourMomentum EtMiss =
        applyProjection<MissingMomentum>(event, "MET").missingMomentum();
    FourMomentum P4H;
    
    auto all_jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::abseta < 4.7 && Cuts::pt > 30*GeV);

    if (tauhad_fh.size() == 0 && taulep_fh.size() != 0) // e mu
    {
        //std::cout << "Case 1\n";
        Particle e;
        Particle mu;
        if (taulep_fh.size() < 2)
            vetoEvent;
        Particle &par1 = taulep_fh[0];
        Particle &par2 = taulep_fh[1];
        if (par1.abspid() == 11 && par2.abspid() == 13) // par1 = e par2 = mu
        {
            e = par1;
            mu = par2;
        } else if (par1.abspid() == 13 &&
                   par2.abspid() == 11) // par1 = mu par2 = e
        {
            e = par2;
            mu = par1;
        } else {
            vetoEvent;
        }
        if (e.pT() < 15 * GeV)
            vetoEvent;
        if (mu.pT() < 24 * GeV)
            vetoEvent;
        if (e.abseta() > 2.4)
            vetoEvent;
        if (mu.abseta() > 2.4)
            vetoEvent;
        FourMomentum LL = e.momentum() + mu.momentum();
        double dphi = deltaPhi(LL, EtMiss);
        double mT = sqrt(2 * LL.pT() * EtMiss.pT() * (1 - cos(dphi)));
        if (mT > 60 * GeV)
            vetoEvent;
        P4H = LL;
        //std::cout << "Endo of case 1\n";
    } else if (tauhad_fh.size() != 0 && taulep_fh.size() == 0) // tau_h tau_h
    {
        //std::cout << "Case 2\n";
        if (tauhad_fh.size() < 2)
            vetoEvent;
        Particle &t1 = tauhad_fh[0];
        Particle &t2 = tauhad_fh[1];
        if (t1.pT() < 40 * GeV || t2.pT() < 40 * GeV)
            vetoEvent;
        if (t1.abseta() > 2.1 || t2.abseta() > 2.1)
            vetoEvent;
        const FastJets &jets = applyProjection<FastJets>(event, "Jets");
        if (jets.size() == 0)
            vetoEvent;
        P4H = t1.momentum() + t2.momentum();
        //std::cout << "Endo of case 2\n";
    } else // tau_h e or tau_h mu
    {
        //std::cout << "Case 3\n";
        Particle &lepton = taulep_fh[0];
        Particle &th = tauhad_fh[0];
        if (lepton.abseta() > 2.1)
            vetoEvent;
        if (lepton.abspid() == 11) {
            if (lepton.pT() < 25 * GeV)
                vetoEvent;
        } else if (lepton.abspid() == 13) {
            if (lepton.pT() < 20 * GeV)
                vetoEvent;
        }
        if (th.pT() < 30 * GeV)
            vetoEvent;
        if (th.abseta() > 2.3)
            vetoEvent;
        FourMomentum L = lepton.momentum();
        double dphi = deltaPhi(L, EtMiss);
        double mT = sqrt(2 * L.pT() * EtMiss.pT() * (1 - cos(dphi)));
        if (mT > 50 * GeV)
            vetoEvent;
        P4H = L;
        //std::cout << "Endo of case 3\n";
    }

    h_pt_h_->fill(P4H.pT() / GeV, weight);
    _h_njets->fill(all_jets.size());
    if(all_jets.size() > 0)
      _h_jet_pt->fill(all_jets[0].pt());

    sumW_ += event.weights()[0];

    return;
}

void Higgs2ttFiducialAndDifferential::finalize() { return; }
}
