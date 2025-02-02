// -*- C++ -*-
// Based on ATLAS_2015_I1394865.cc
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Math/MathUtils.hh"

namespace Rivet {


  /// Inclusive 4-lepton lineshape
  class Higgs2ZZFiducialAndDifferential : public Analysis {
  public:

    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(Higgs2ZZFiducialAndDifferential);


    void init() {
      sumW_ = 0.;
      // Everything super loose
      FinalState fs(Cuts::abseta < 10.0);

      IdentifiedFinalState photon(fs, PID::PHOTON);
      IdentifiedFinalState bare_EL(fs, {PID::ELECTRON, -PID::ELECTRON});
      IdentifiedFinalState bare_MU(fs, {PID::MUON, -PID::MUON});

      // Selection 1: ZZ-> llll selection
      Cut etaranges_el = Cuts::abseta < 2.5 && Cuts::pT > 7*GeV;
      Cut etaranges_mu = Cuts::abseta < 2.7 && Cuts::pT > 5*GeV;

      DressedLeptons electron_sel4l(photon, bare_EL, 0.3, etaranges_el);
      declare(electron_sel4l, "ELECTRON_sel4l");
      DressedLeptons muon_sel4l(photon, bare_MU, 0.3, etaranges_mu);
      declare(muon_sel4l, "MUON_sel4l");

      //---Jets
      FastJets fs_jets(fs, FastJets::ANTIKT, 0.4);
      declare(fs_jets, "JETS");

      // Both ZZ on-shell histos
      book(_h_ZZ_pTZZ, "pt_h", {0,10,20,30,45,60,80,120,200,10000});
      book(_h_ZZ_pTZZ_coarser, "mom_hc", {0,30,60,120,10000});
      book(_h_jet_pt, "pt_j0",{30,55,95,200,500});
      book(_h_njets, "njets", {-0.5,0.5,1.5,2.5,3.5,100.5});
      //book(_h_rapidity, "h_rapidity", {0.0,0.15,0.3,0.45,0.6,0.75,0.9,1.2,1.6,2.5});
      book(_h_deta, "deta_jj", {0.0,1.6,3.0,1000});
      book(_h_deltaphijj, "deltaphijj", {-M_PI, -M_PI/2, 0, M_PI/2, M_PI});
      //book(_h_deltaphijj, "deltaphijj", {0, M_PI/2, M_PI, 3*M_PI/2, 2*M_PI});
      book(_h_sigma, "h_sigma", 1, 0, 100000000);
      book(_h_sigma_post, "h_sigma_post", 1, 0, 100000000);
    }


    /// Do the analysis
    void analyze(const Event& e) {
       sumW_ += e.weights()[0];
      _h_sigma->fill(sumW_);
     
      // jets
      auto jets = apply<JetAlg>(e, "JETS").jetsByPt(Cuts::abseta < 4.7 && Cuts::pt > 30*GeV);

      ////////////////////////////////////////////////////////////////////
      // Preselection of leptons for ZZ-> llll final state
      ////////////////////////////////////////////////////////////////////

      Particles leptons_sel4l;
      const vector<DressedLepton>& mu_sel4l = apply<DressedLeptons>(e, "MUON_sel4l").dressedLeptons();
      const vector<DressedLepton>& el_sel4l = apply<DressedLeptons>(e, "ELECTRON_sel4l").dressedLeptons();
      const vector<DressedLepton> leptonsFS_sel4l = mu_sel4l + el_sel4l;
      // leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), mu_sel4l.begin(), mu_sel4l.end() );
      // leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), el_sel4l.begin(), el_sel4l.end() );

      // mu: pT > 6 GeV, eta < 2.7; ele: pT > 7 GeV, eta < 2.5
      // Completely pointless since the selection in the constructor is stricter
      for (const DressedLepton& l : leptonsFS_sel4l) {
        if (l.abspid() == PID::ELECTRON) leptons_sel4l.push_back(l);  // REDUNDANT: if (l.pT() > 7*GeV && l.abseta() < 2.5)
        else if (l.abspid() == PID::MUON) leptons_sel4l.push_back(l); // REDUNDANT: if (l.pT() > 6*GeV && l.abseta() < 2.7)
      }

      //////////////////////////////////////////////////////////////////
      // Exactly two opposite charged leptons
      //////////////////////////////////////////////////////////////////

      // Calculate total 'flavour' charge
      double totalcharge = 0;
      for (const Particle& l : leptons_sel4l)  totalcharge += l.pid();

      // Analyze 4 lepton events
      // Get only events with two OS pairs
      if (leptons_sel4l.size() != 4 || totalcharge != 0) vetoEvent;

      // Identify Z states from 4 lepton pairs
      Zstate Z1, Z2, Z1_alt, Z2_alt;
      if ( !identifyZstates(Z1, Z2, Z1_alt, Z2_alt, leptons_sel4l) )  vetoEvent;

      const double mZ1 = Z1.mom().mass();
      const double mZ2 = Z2.mom().mass();
      const double mZ1_alt = Z1_alt.mom().mass();
      const double mZ2_alt = Z2_alt.mom().mass();
      const double pTZ1 = Z1.mom().pT();
      const double pTZ2 = Z2.mom().pT();
      const double mZZ = (Z1.mom() + Z2.mom()).mass();
      const double pTZZ = (Z1.mom() + Z2.mom()).pT();
      //std::cout << pTZZ << std::endl;
      //std::cout << mZ1 << " " << mZ2 << " " << mZ1_alt << " " << mZ2_alt << std::endl;

      // Event selections
      // pT(Z) > 2 GeV
      bool pass = pTZ1 > 2*GeV && pTZ2 > 2*GeV;
      if (!pass) vetoEvent;

      // Lepton kinematics: pT > 20, 10, 7 (5 if muon) GeV
      int n1 = 0, n2 = 0, n3 = 0;
      for (Particle& l : leptons_sel4l) {
        if (l.pT() > 20*GeV) ++n1;
        if (l.pT() > 10*GeV) ++n2;
        if (l.pT() > 7*GeV && l.abspid() == PID::ELECTRON) ++n3;
        if (l.pT() > 5*GeV && l.abspid() == PID::MUON) ++n3;
      }
      pass = pass && n1>=1 && n2>=2 && n3>=3;
      if (!pass) vetoEvent;

      // Dilepton mass: 40 < mZ1 < 120 GeV, 12 < mZ2 < 120 GeV
      pass = pass && mZ1 > 40*GeV && mZ1 < 120*GeV;
      pass = pass && mZ2 > 12*GeV && mZ2 < 120*GeV;
      if (!pass) vetoEvent;

      // Lepton separation: deltaR(l, l') > 0.02 
      for (size_t i = 0; i < leptons_sel4l.size(); ++i) {
        for (size_t j = i + 1; j < leptons_sel4l.size(); ++j) {
          const Particle& l1 = leptons_sel4l[i];
          const Particle& l2 = leptons_sel4l[j];
          //pass = pass && deltaR(l1, l2) > (l1.abspid() == l2.abspid() ? 0.1 : 0.2);
          pass = pass && deltaR(l1, l2) > 0.02;
          if (!pass) vetoEvent;
        }
      }

      // J/Psi veto: m(l+l-) > 4 GeV
      pass = pass && mZ1 > 4*GeV && mZ2 > 4*GeV && mZ1_alt > 4*GeV && mZ2_alt > 4*GeV;
      if (!pass) vetoEvent;

      // 105 < m4l < 160 GeV
      pass = pass && mZZ > 105*GeV && mZZ < 160*GeV;
      if (!pass) vetoEvent;

      // Fill histograms
      //std::cout << pTZZ << std::endl;
      _h_ZZ_pTZZ->fill(pTZZ);
      _h_ZZ_pTZZ_coarser->fill(pTZZ);
      _h_njets->fill(jets.size());
      if(jets.size() > 0)
        _h_jet_pt->fill(jets[0].pt());
      if(jets.size() > 1)
          _h_deta->fill(fabs(deltaEta(jets[0], jets[1])));
      if(jets.size() > 1)
          //std::cout << jets[0].phi() << std::endl;
          //std::cout << jets[1].phi() << std::endl;
          //std::cout << deltaPhi(jets[0].ptvec(), jets[1].ptvec()) << std::endl;
          //std::cout << std::endl;
          _h_deltaphijj->fill(deltaPhiCustom(jets[0].phi(), jets[1].phi()));
      /*for (size_t i = 0; i < jets.size()-1; ++i) {
          for (size_t j = i+1; j<jets.size(); ++j) {
              const Jet& jet1 = jets[i];
              const Jet& jet2 = jets[j];
              const double deltaPhiJJ = deltaPhi(jet1.momentum(), jet2.momentum());
              _h_deltaphijj->fill(deltaPhiJJ);
          }
      }*/
      _h_sigma_post->fill(sumW_);
   
    }


    /// Finalize
    void finalize() {
      //const double norm = crossSection()/sumOfWeights()/femtobarn/TeV;
      //scale(_h_ZZ_pTZZ, norm);
    }


    double deltaPhiCustom(double phi1, double phi2) {
      const double x = mapAngleMPiToPi(phi1 - phi2);
      return x;
    }


    /// Generic Z candidate
    struct Zstate : public ParticlePair {
      Zstate() { }
      Zstate(ParticlePair _particlepair) : ParticlePair(_particlepair) { }
      FourMomentum mom() const { return first.momentum() + second.momentum(); }
      operator FourMomentum() const { return mom(); }
      static bool cmppT(const Zstate& lx, const Zstate& rx) { return lx.mom().pT() < rx.mom().pT(); }
    };


    /// @brief 4l to ZZ assignment algorithm
    ///
    /// ZZ->4l pairing
    /// - At least two same flavour opposite sign (SFOS) lepton pairs
    /// - Ambiguities in pairing are resolved following the procedure
    ///   1. the leading Z (Z1) is choosen as the SFOS with dilepton mass closer to Z mass
    ///   2. the subleading Z (Z2) is choosen as the remaining SFOS dilepton pair
    ///
    /// Z1, Z2: the selected pairing
    /// Z1_alt, Z2_alt: the alternative pairing (the same as Z1, Z2 in 2e2m case)
    bool identifyZstates(Zstate& Z1, Zstate& Z2, Zstate& Z1_alt, Zstate& Z2_alt, const Particles& leptons_sel4l) {
      const double ZMASS = 91.1876*GeV;
      bool findZZ = false;

      Particles part_pos_el, part_neg_el, part_pos_mu, part_neg_mu;
      for (const Particle& l : leptons_sel4l) {
        if (l.abspid() == PID::ELECTRON) {
          if (l.pid() < 0) part_neg_el.push_back(l);
          if (l.pid() > 0) part_pos_el.push_back(l);
        }
        else if (l.abspid() == PID::MUON) {
          if (l.pid() < 0) part_neg_mu.push_back(l);
          if (l.pid() > 0) part_pos_mu.push_back(l);
        }
      }

      // eeee/mmmm channel
      if ((part_neg_el.size() == 2 && part_pos_el.size() == 2) || (part_neg_mu.size() == 2 && part_pos_mu.size() == 2)) {
        findZZ = true;

        Zstate Zcand_1, Zcand_2, Zcand_3, Zcand_4;
        Zstate Zcand_1_tmp, Zcand_2_tmp, Zcand_3_tmp, Zcand_4_tmp;
        if (part_neg_el.size() == 2) { // eeee
          Zcand_1_tmp = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[0] ) );
          Zcand_2_tmp = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[1] ) );
          Zcand_3_tmp = Zstate( ParticlePair( part_neg_el[1],  part_pos_el[0] ) );
          Zcand_4_tmp = Zstate( ParticlePair( part_neg_el[1],  part_pos_el[1] ) );
        }
        else { // mmmm
          Zcand_1_tmp = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[0] ) );
          Zcand_2_tmp = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[1] ) );
          Zcand_3_tmp = Zstate( ParticlePair( part_neg_mu[1],  part_pos_mu[0] ) );
          Zcand_4_tmp = Zstate( ParticlePair( part_neg_mu[1],  part_pos_mu[1] ) );
        }

        // We can have the following pairs: (Z1 + Z4) || (Z2 + Z3)
        // Firstly, reorder withing each quadruplet to have
        //  - fabs(mZ1 - ZMASS) < fabs(mZ4 - ZMASS)
        //  - fabs(mZ2 - ZMASS) < fabs(mZ3 - ZMASS)
        if (fabs(Zcand_1_tmp.mom().mass() - ZMASS) < fabs(Zcand_4_tmp.mom().mass() - ZMASS)) {
          Zcand_1 = Zcand_1_tmp;
          Zcand_4 = Zcand_4_tmp;
        } else {
          Zcand_1 = Zcand_4_tmp;
          Zcand_4 = Zcand_1_tmp;
        }
        if (fabs(Zcand_2_tmp.mom().mass() - ZMASS) < fabs(Zcand_3_tmp.mom().mass() - ZMASS)) {
          Zcand_2 = Zcand_2_tmp;
          Zcand_3 = Zcand_3_tmp;
        } else {
          Zcand_2 = Zcand_3_tmp;
          Zcand_3 = Zcand_2_tmp;
        }

        // We can have the following pairs: (Z1 + Z4) || (Z2 + Z3)
        // Secondly, select the leading and subleading Z following
        //   1. the leading Z (Z1) is choosen as the SFOS with dilepton mass closet to Z mass
        //   2. the subleading Z (Z2) is choosen as the remaining SFOS dilepton pair
        if (fabs(Zcand_1.mom().mass() - ZMASS) < fabs(Zcand_2.mom().mass() - ZMASS)) {
          Z1 = Zcand_1;
          Z2 = Zcand_4;
          Z1_alt = Zcand_2;
          Z2_alt = Zcand_3;
        } else {
          Z1 = Zcand_2;
          Z2 = Zcand_3;
          Z1_alt = Zcand_1;
          Z2_alt = Zcand_4;
        }
      } // end of eeee/mmmm channel
      else if (part_neg_el.size() == 1 && part_pos_el.size() == 1 && part_neg_mu.size() == 1 && part_pos_mu.size() == 1) { // 2e2m channel
        findZZ = true;

        Zstate Zcand_1, Zcand_2;

        Zcand_1 = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[0] ) );
        Zcand_2 = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[0] ) );

        if (fabs(Zcand_1.mom().mass() - ZMASS) < fabs(Zcand_2.mom().mass() - ZMASS)) {
          Z1 = Zcand_1;
          Z2 = Zcand_2;
        } else {
          Z1 = Zcand_2;
          Z2 = Zcand_1;
        }
        Z1_alt = Z1;
        Z2_alt = Z2;
      }

      return findZZ;
    }


  private:

    double     sumW_;
    Histo1DPtr _h_ZZ_pTZZ;
    Histo1DPtr _h_ZZ_pTZZ_coarser;
    Histo1DPtr _h_jet_pt;
    Histo1DPtr _h_njets;
    //Histo1DPtr _h_rapidity;
    Histo1DPtr _h_deta;
    Histo1DPtr _h_deltaphijj;
    Histo1DPtr _h_sigma;
    Histo1DPtr _h_sigma_post;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(Higgs2ZZFiducialAndDifferential);

}
