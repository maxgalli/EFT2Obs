#include "Higgs2GGFiducialAndDifferential.h"

namespace Rivet
{
    void Higgs2GGFiducialAndDifferential::init()
    {
        //---All particle final state
        FinalState fs;       
        declare(fs, "FS");       

        //---Photons
        IdentifiedFinalState fs_photons;
        fs_photons.acceptId(PID::PHOTON);
        declare(fs_photons, "FS_PHOTONS");

        //---Jets
        FastJets fs_jets(fs, FastJets::ANTIKT, 0.4);
        declare(fs_jets, "JETS");

        //---Histograms
        book(h_m_gg_, "mass_gg", 100, 120, 130);
        book(h_pt_gg_, "pt_gg",{0,5,10,15,20,25,30,35,45,60,80,100,120,140,170,200,250,350,450,10000});
        book(h_eta_gg_, "eta_gg", {-5, -4, -3, -2,-1.5, -1, 0, 1,1.5, 2, 3, 4, 5});
        book(h_pt_g1_, "pt_g1",{0,10,20,30,40,50,70,90,100,150,200,250,300,400, 500});
        book(h_pt_g2_, "pt_g2", {0,10,20,30,40,50,70,90,100,150,200,250,300,400, 500});
        book(h_eta_g1_, "eta_g1", {-5, -4, -3, -2,-1.5, -1, 0, 1,1.5, 2, 3, 4, 5});
        book(h_eta_g2_, "eta_g2", {-5, -4, -3, -2,-1.5, -1, 0, 1,1.5, 2, 3, 4, 5});
        book(h_njets_eta2p5_, "n_jets_eta2p5", {-0.5,0.5,1.5,2.5,3.5,100.5});
        book(h_jet_pt_eta2p5_, "h_jet_pt_eta2p5", {30,40,50,70,90,100,150,200,250,300,400, 500});
        book(h_njets_eta4p7_, "n_jets_eta4p7", {-0.5,0.5,1.5,2.5,3.5,100.5});
        book(h_jet_pt_eta4p7_, "h_jet_pt_eta4p7",{30,40,50,70,90,100,150,200,250,300,400, 500});
        book(h_deltaphijj_, "h_deltaphijj", {0.0, 0.5, 0.9, 1.3, 1.7, 2.5, M_PI});
    }

    void Higgs2GGFiducialAndDifferential::analyze(const Event& event)
    {               
        auto photons = apply<IdentifiedFinalState>(event, "FS_PHOTONS").particlesByPt();       

        if (photons.size() < 2) 
            vetoEvent;       

        //---Isolate photons with ET_sum in cone
        Particles isolated_photons;
        Particles fs = apply<FinalState>(event, "FS").particlesByPt();
        for(auto& photon : photons) 
        {
            FourMomentum mom_in_cone;
            for(auto& particle : fs) 
            {
                if(deltaR(photon, particle) < 0.4 && photon.p() != particle.p())
                    mom_in_cone += particle.momentum();
            }
            if(mom_in_cone.pt() < 10*GeV)
                isolated_photons.push_back(photon);
        }

        if(isolated_photons.size() < 2) 
            vetoEvent;

        //---kinematic photon selection
        FourMomentum mom_PP = isolated_photons[0].momentum() + isolated_photons[1].momentum();
        if(isolated_photons[0].pT() < mom_PP.mass() / 3. || isolated_photons[1].pT() < mom_PP.mass() / 4.) 
            vetoEvent;

        h_m_gg_->fill(mom_PP.mass());
        h_pt_gg_->fill(mom_PP.pt());
        h_eta_gg_->fill(mom_PP.eta());
        h_pt_g1_->fill(isolated_photons[0].pT());
        h_pt_g2_->fill(isolated_photons[1].pT());
        h_eta_g1_->fill(isolated_photons[0].eta());
        h_eta_g2_->fill(isolated_photons[1].eta());

        //---jets
        auto jets_eta2p5 = apply<JetAlg>(event, "JETS").jetsByPt(Cuts::abseta < 2.5 && Cuts::pt > 30*GeV);
        auto jets_eta4p7 = apply<JetAlg>(event, "JETS").jetsByPt(Cuts::abseta < 4.7 && Cuts::pt > 30*GeV);
        
        h_njets_eta2p5_->fill(jets_eta2p5.size());
        if(jets_eta2p5.size()>0)
            h_jet_pt_eta2p5_->fill(jets_eta2p5[0].pt());
        h_njets_eta4p7_->fill(jets_eta4p7.size());
        if(jets_eta4p7.size()>0)
            h_jet_pt_eta4p7_->fill(jets_eta4p7[0].pt());

        for (size_t i = 0; i < jets_eta4p7.size()-1; ++i) {
            for (size_t j = i+1; j<jets_eta4p7.size(); ++j) {
                const Jet& jet1 = jets_eta4p7[i];
                const Jet& jet2 = jets_eta4p7[j];
                const double deltaPhiJJ = deltaPhi(jet1.momentum(), jet2.momentum());
                h_deltaphijj_->fill(deltaPhiJJ);
            }
        }

        sumW_ += event.weights()[0];

        return;
    }

    void Higgs2GGFiducialAndDifferential::finalize()
    {
        return;
    }
}
