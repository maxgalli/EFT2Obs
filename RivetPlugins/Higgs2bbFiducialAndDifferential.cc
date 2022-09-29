#include "Higgs2bbFiducialAndDifferential.h"

/* 
 * Based on https://arxiv.org/pdf/1709.05543.pdf
*/


namespace Rivet
{
    void Higgs2bbFiducialAndDifferential::init()
    {
        // Initialize the projections
        const FinalState fs;
        declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");
        declare(FastJets(fs, FastJets::ANTIKT, 0.8), "JetsAK8");
        
        //---Histograms
        book(_h_pt_h, "pt_h",{0,200,300,450,650,10000});
    }

    void Higgs2bbFiducialAndDifferential::analyze(const Event& event)
    {               
        const double weight = 1.0;
        
        const FastJets& fjAK4 = applyProjection<FastJets>(event, "JetsAK4");
        const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(30*GeV, 10000.0*GeV));
        Jets jetsAK4_good;
        Jets jetsAK8_good;
        for (const Jet& j : jetsAK4) {
            if (!j.bTagged()) continue; 
            if (j.neutralEnergy()/j.totalEnergy() > 0.9) continue;
            if (j.hadronicEnergy() == 0) continue;
            
            Particles constituents = j.constituents();
            bool refuse = false;
            for (const Particle& c : constituents) {
                auto id = c.abspid();
                    if (id == 11 && c.pT() > 10*GeV && c.abseta() < 2.5) refuse = true;
                    if (id == 13 && c.pT() > 10*GeV && c.abseta() < 2.4) refuse = true;
                    if (id == 17 && c.pT() > 18*GeV && c.abseta() < 2.3) refuse = true;
            }
            if (refuse == true) continue;
            jetsAK4_good.push_back(j);
        }
      
        const FastJets& fjAK8 = applyProjection<FastJets>(event, "JetsAK8");
        const Jets& jetsAK8 = fjAK8.jetsByPt(Cuts::ptIn(450*GeV, 10000.0*GeV));
        //const Jets& jetsAK8 = fjAK8.jetsByPt();
        for (const Jet& j : jetsAK8) {
            if (!j.bTagged()) continue; 
            if (j.neutralEnergy()/j.totalEnergy() > 0.9) continue;
            if (j.hadronicEnergy() == 0) continue;
            Particles constituents = j.constituents();
            bool refuse = false;
            for (const Particle& c : constituents) {
                auto id = c.abspid();
                    if (id == 11 && c.pT() > 10*GeV && c.abseta() < 2.5) refuse = true;
                    if (id == 13 && c.pT() > 10*GeV && c.abseta() < 2.4) refuse = true;
                    if (id == 17 && c.pT() > 18*GeV && c.abseta() < 2.3) refuse = true;
                    if (id == 22 && c.pT() > 175*GeV) refuse = true;
            }
            if (refuse == true) continue;
            jetsAK8_good.push_back(j);
        }

        Jets final_candidates;
        for (const Jet& j8 : jetsAK8_good) {
            bool too_far = false;
            for (const Jet& j4 : jetsAK4_good) {
                auto dp = deltaPhi(j4, j8);
                //std::cout << dp << std::endl;
                if (dp > 1.57075) too_far = true;
            }
            if (!too_far) final_candidates.push_back(j8);
        }
        //std::cout << jetsAK8_good.size();
        std::cout << final_candidates.size() << std::endl;

        if (final_candidates.size() > 0){
            Jet final_candidate = final_candidates[0];
            //std::cout << final_candidate.pT() << std::endl;
            _h_pt_h->fill(final_candidate.pT() / GeV, weight);
            _sumW += event.weights()[0];
        }
       
        return;
    }

    void Higgs2bbFiducialAndDifferential::finalize()
    {
        return;
    }
}
