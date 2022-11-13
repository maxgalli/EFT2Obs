/*
 * Based on what is described at page 3 of the paper
 * Still missing:
 * HT > 200 GeV
 */

#include "Higgs2ttBoostedFiducialAndDifferential.h"

namespace Rivet {
void Higgs2ttBoostedFiducialAndDifferential::init() {
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
    
    declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");
    declare(FastJets(fs, FastJets::ANTIKT, 0.8), "JetsAK8");

    //---Histograms
    book(h_pt_h_, "pt_h", { 0, 450, 600, 10000 });
}

void Higgs2ttBoostedFiducialAndDifferential::analyze(const Event &event) {
    const double weight = 1.0;

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

    // no tau candidates with pt below 30 GeV
    bool no_taus_above_30 = true;
    for(Particle& t : tauhad) {
        if (t.pT() > 30 * GeV) no_taus_above_30 = false;
    }
    for(Particle& t : taulep) {
        if (t.pT() > 30 * GeV) no_taus_above_30 = false;
    }
    if (no_taus_above_30 == true) vetoEvent;

    //std::cout << tauhad_fh.size();
    //std::cout << taulep_fh.size();
    
    const FastJets& fjAK4 = applyProjection<FastJets>(event, "JetsAK4");
    const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(30*GeV, 10000.0*GeV));
    const FastJets& fjAK8 = applyProjection<FastJets>(event, "JetsAK8");
    const Jets& jetsAK8 = fjAK8.jetsByPt(Cuts::ptIn(30*GeV, 10000.0*GeV));
    Jets jets;
    for (const Jet& j : jetsAK4) {
        jets.push_back(j);    
    }
    for (const Jet& j : jetsAK8) {
        jets.push_back(j);    
    }
    
    // no jets from b quarks in the event
    for (const Jet& j : jets) {
        if (j.bTagged()) vetoEvent;
    }
    std::cout << "Size of jets " << jets.size() << std::endl;

    FourMomentum EtMiss =
        applyProjection<MissingMomentum>(event, "MET").missingMomentum();
    FourMomentum P4H;

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
        // opposite charge
        if (e.charge() == mu.charge()) vetoEvent;
        // deltaR < 0.8
        if (deltaR(e, mu) > 0.8) vetoEvent;
        FourMomentum LL = e.momentum() + mu.momentum();
        double dphi = deltaPhi(LL, EtMiss);
        double mT = sqrt(2 * LL.pT() * EtMiss.pT() * (1 - cos(dphi)));
        if (mT < 80 * GeV)
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
        // deltaR < 0.8
        if (deltaR(t1, t2) > 0.8) vetoEvent;
        if (jets.size() == 0)
            vetoEvent;
        P4H = t1.momentum() + t2.momentum();
        //std::cout << "Endo of case 2\n";
    } else // tau_h e or tau_h mu
    {
        //std::cout << "Case 3\n";
        Particle &lepton = taulep_fh[0];
        Particle &th = tauhad_fh[0];
        // deltaR < 0.8
        if (deltaR(th, lepton) > 0.8) vetoEvent;
        FourMomentum L = lepton.momentum();
        double dphi = deltaPhi(L, EtMiss);
        double mT = sqrt(2 * L.pT() * EtMiss.pT() * (1 - cos(dphi)));
        if (mT > 80 * GeV)
            vetoEvent;
        P4H = L + th.momentum();
        //std::cout << "Endo of case 3\n";
    }

    // invariant mass of di-tau above 250 GeV
    //if (P4H.mass() < 250 * GeV) vetoEvent;
    std::cout << "We got something" << std::endl;

    h_pt_h_->fill(P4H.pT() / GeV, weight);

    sumW_ += event.weights()[0];

    return;
}

void Higgs2ttBoostedFiducialAndDifferential::finalize() { return; }
}
