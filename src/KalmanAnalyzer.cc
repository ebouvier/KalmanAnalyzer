// -*- C++ -*-
//
// Package:    KalmanAnalyzer
// Class:      KalmanAnalyzer
// 
/**\class KalmanAnalyzer KalmanAnalyzer.cc UserCode/KalmanAnalyzer/src/KalmanAnalyzer.cc

Description: Reconstruct Bpm -> D0 + mupm (+ nu) -> K Pi + mu (+ nu)

Implementation: Using a simple Kalman vertex fitter for the D0
*/


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "UserCode/KalmanAnalyzer/interface/TopTriggerEfficiencyProvider.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
// #include "RecoVertex/KinematicFitPrimitives/interface/"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TH1D.h"
#include "TLorentzVector.h"

//
// class declaration
//

class KalmanAnalyzer : public edm::EDAnalyzer {
  public:
    explicit KalmanAnalyzer(const edm::ParameterSet&);
    ~KalmanAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    // ----------member data ---------------------------
    
    // evts properties
    unsigned int nEvts;
    unsigned int nEvts2;
    double nEvts3;
    bool hasGoodMuons;
    bool hasGoodElectrons;
    bool hasGoodLeptons;
    bool hasGoodJets;
    bool isGoodEvt;

    double weight;

    TH1D* h_nJets;
    TH1D* h_CSV;

    TH1D* h_genD0_n;
    TH1D* h_genD0_p;
    TH1D* h_genD0_pT;
    TH1D* h_genD0_eta;
    TH1D* h_genD0_phi;

    // simple Kalman Vertex Fitter
    TH1D* h_B_cuts;

    TH1D* h_candGenRecoD0dR;
    TH1D* h_candGenRecoTrdR;

    TH1D* h_D0Cand_Chi2NDOF;

    TH1D* h_D0Cand_MassChi2Inf1;
    TH1D* h_D0Cand_MassChi2Inf1p5;
    TH1D* h_D0Cand_MassChi2Inf2;
    TH1D* h_D0Cand_MassChi2Inf2p5;
    TH1D* h_D0Cand_MassChi2Inf3;
    TH1D* h_D0Cand_MassChi2Inf3p5;
    TH1D* h_D0Cand_MassChi2Inf4;
    TH1D* h_D0Cand_MassChi2Inf4p5;
    TH1D* h_D0Cand_MassChi2Inf5;
    TH1D* h_D0Cand_MassChi2Inf5p5;
    TH1D* h_D0Cand_MassChi2Inf6;
    TH1D* h_D0Cand_MassChi2Inf6p5;
    TH1D* h_D0Cand_MassChi2Inf7;
    TH1D* h_D0Cand_MassChi2Inf7p5;
    TH1D* h_D0Cand_MassChi2Inf8;

    TH1D* h_D0Cand_L;
    TH1D* h_D0Cand_SigmaL;
    TH1D* h_D0Cand_LOverSigmaL;

    TH1D* h_D0Cand_pT;

    TH1D* h_D0_p;
    TH1D* h_D0_pT;
    TH1D* h_D0_eta;
    TH1D* h_D0_phi;
    TH1D* h_D0_L;
    TH1D* h_D0_SigmaL;
    TH1D* h_D0_dRJet;
    TH1D* h_D0_Mass;
    TH1D* h_D0_paired_p;
    TH1D* h_D0_paired_pT;
    TH1D* h_D0_paired_eta;
    TH1D* h_D0_paired_phi;
    TH1D* h_D0_paired_L;
    TH1D* h_D0_paired_SigmaL;
    TH1D* h_D0_paired_dRJet;
    TH1D* h_D0_paired_Mass;
    TH1D* h_D0_switched_p;
    TH1D* h_D0_switched_pT;
    TH1D* h_D0_switched_eta;
    TH1D* h_D0_switched_phi;
    TH1D* h_D0_switched_L;
    TH1D* h_D0_switched_SigmaL;
    TH1D* h_D0_switched_dRJet;
    TH1D* h_D0_switched_Mass;
    TH1D* h_D0_unmatched_p;
    TH1D* h_D0_unmatched_pT;
    TH1D* h_D0_unmatched_eta;
    TH1D* h_D0_unmatched_phi;
    TH1D* h_D0_unmatched_L;
    TH1D* h_D0_unmatched_SigmaL;
    TH1D* h_D0_unmatched_dRJet;
    TH1D* h_D0_unmatched_Mass;

    TH1D* h_BCand_DeltaRD0Mu;
    TH1D* h_B_DeltaRD0Mu;
    TH1D* h_B_pD0pMu;
    TH1D* h_B_pTD0pTMu;
    TH1D* h_B_D0Mass;
    TH1D* h_B_Mass;
    TH1D* h_B_p;
    TH1D* h_B_pT;
    TH1D* h_B_eta;
    TH1D* h_B_phi;
    TH1D* h_B_paired_DeltaRD0Mu;
    TH1D* h_B_paired_pD0pMu;
    TH1D* h_B_paired_pTD0pTMu;
    TH1D* h_B_paired_D0Mass;
    TH1D* h_B_paired_Mass;
    TH1D* h_B_paired_p;
    TH1D* h_B_paired_pT;
    TH1D* h_B_paired_eta;
    TH1D* h_B_paired_phi;
    TH1D* h_B_switched_DeltaRD0Mu;
    TH1D* h_B_switched_pD0pMu;
    TH1D* h_B_switched_pTD0pTMu;
    TH1D* h_B_switched_D0Mass;
    TH1D* h_B_switched_Mass;
    TH1D* h_B_switched_p;
    TH1D* h_B_switched_pT;
    TH1D* h_B_switched_eta;
    TH1D* h_B_switched_phi;
    TH1D* h_B_unmatched_DeltaRD0Mu;
    TH1D* h_B_unmatched_pD0pMu;
    TH1D* h_B_unmatched_pTD0pTMu;
    TH1D* h_B_unmatched_D0Mass;
    TH1D* h_B_unmatched_Mass;
    TH1D* h_B_unmatched_p;
    TH1D* h_B_unmatched_pT;
    TH1D* h_B_unmatched_eta;
    TH1D* h_B_unmatched_phi;

    //constrained Kalman Vertex Fitter
    TH1D* h_D0consCand_L;
    TH1D* h_D0consCand_SigmaL;
    TH1D* h_D0consCand_LOverSigmaL;

    TH1D* h_D0consCand_pT;

    TH1D* h_D0cons_Chi2NDOF;
    TH1D* h_D0cons_L;
    TH1D* h_D0cons_SigmaL;
    TH1D* h_D0cons_p;
    TH1D* h_D0cons_pT;
    TH1D* h_D0cons_eta;
    TH1D* h_D0cons_phi;
    TH1D* h_D0cons_dRJet;
    TH1D* h_D0cons_Mass;
    TH1D* h_D0cons_paired_Chi2NDOF;
    TH1D* h_D0cons_paired_L;
    TH1D* h_D0cons_paired_SigmaL;
    TH1D* h_D0cons_paired_p;
    TH1D* h_D0cons_paired_pT;
    TH1D* h_D0cons_paired_eta;
    TH1D* h_D0cons_paired_phi;
    TH1D* h_D0cons_paired_dRJet;
    TH1D* h_D0cons_paired_Mass;
    TH1D* h_D0cons_switched_Chi2NDOF;
    TH1D* h_D0cons_switched_L;
    TH1D* h_D0cons_switched_SigmaL;
    TH1D* h_D0cons_switched_p;
    TH1D* h_D0cons_switched_pT;
    TH1D* h_D0cons_switched_eta;
    TH1D* h_D0cons_switched_phi;
    TH1D* h_D0cons_switched_dRJet;
    TH1D* h_D0cons_switched_Mass;
    TH1D* h_D0cons_unmatched_Chi2NDOF;
    TH1D* h_D0cons_unmatched_L;
    TH1D* h_D0cons_unmatched_SigmaL;
    TH1D* h_D0cons_unmatched_p;
    TH1D* h_D0cons_unmatched_pT;
    TH1D* h_D0cons_unmatched_eta;
    TH1D* h_D0cons_unmatched_phi;
    TH1D* h_D0cons_unmatched_dRJet;
    TH1D* h_D0cons_unmatched_Mass;

    TH1D* h_BCand_DeltaRD0consMu;
    TH1D* h_B_DeltaRD0consMu;
    TH1D* h_B_pD0conspMu;
    TH1D* h_B_pTD0conspTMu;
    TH1D* h_B_D0consMass;
    TH1D* h_B_consMass;
    TH1D* h_B_consp;
    TH1D* h_B_conspT;
    TH1D* h_B_conseta;
    TH1D* h_B_consphi;
    TH1D* h_B_paired_DeltaRD0consMu;
    TH1D* h_B_paired_pD0conspMu;
    TH1D* h_B_paired_pTD0conspTMu;
    TH1D* h_B_paired_D0consMass;
    TH1D* h_B_paired_consMass;
    TH1D* h_B_paired_consp;
    TH1D* h_B_paired_conspT;
    TH1D* h_B_paired_conseta;
    TH1D* h_B_paired_consphi;
    TH1D* h_B_switched_DeltaRD0consMu;
    TH1D* h_B_switched_pD0conspMu;
    TH1D* h_B_switched_pTD0conspTMu;
    TH1D* h_B_switched_D0consMass;
    TH1D* h_B_switched_consMass;
    TH1D* h_B_switched_consp;
    TH1D* h_B_switched_conspT;
    TH1D* h_B_switched_conseta;
    TH1D* h_B_switched_consphi;
    TH1D* h_B_unmatched_DeltaRD0consMu;
    TH1D* h_B_unmatched_pD0conspMu;
    TH1D* h_B_unmatched_pTD0conspTMu;
    TH1D* h_B_unmatched_D0consMass;
    TH1D* h_B_unmatched_consMass;
    TH1D* h_B_unmatched_consp;
    TH1D* h_B_unmatched_conspT;
    TH1D* h_B_unmatched_conseta;
    TH1D* h_B_unmatched_consphi;

    // simple invariant combination
    TH1D* h_candGenRecoD0combidR;
    TH1D* h_candGenRecoTrcombidR;
    
    TH1D* h_D0combiCand_pT;

    TH1D* h_dRTr1Tr2combi;
    TH1D* h_dEtaTr1Tr2combi;
    TH1D* h_dPhiTr1Tr2combi;
    TH1D* h_pTrCombi;
    TH1D* h_pTTrCombi;
    TH1D* h_etaTrCombi;
    TH1D* h_phiTrCombi;
    TH1D* h_pTr1pTr2combi;
    TH1D* h_pTTr1pTTr2combi;
    TH1D* h_dRTr1Tr2combi_paired;
    TH1D* h_dEtaTr1Tr2combi_paired;
    TH1D* h_dPhiTr1Tr2combi_paired;
    TH1D* h_pTrCombi_paired;
    TH1D* h_pTTrCombi_paired;
    TH1D* h_etaTrCombi_paired;
    TH1D* h_phiTrCombi_paired;
    TH1D* h_pTr1pTr2combi_paired;
    TH1D* h_pTTr1pTTr2combi_paired;
    TH1D* h_dRTr1Tr2combi_switched;
    TH1D* h_dEtaTr1Tr2combi_switched;
    TH1D* h_dPhiTr1Tr2combi_switched;
    TH1D* h_pTrCombi_switched;
    TH1D* h_pTTrCombi_switched;
    TH1D* h_etaTrCombi_switched;
    TH1D* h_phiTrCombi_switched;
    TH1D* h_pTr1pTr2combi_switched;
    TH1D* h_pTTr1pTTr2combi_switched;
    TH1D* h_dRTr1Tr2combi_unmatched;
    TH1D* h_dEtaTr1Tr2combi_unmatched;
    TH1D* h_dPhiTr1Tr2combi_unmatched;
    TH1D* h_pTrCombi_unmatched;
    TH1D* h_pTTrCombi_unmatched;
    TH1D* h_etaTrCombi_unmatched;
    TH1D* h_phiTrCombi_unmatched;
    TH1D* h_pTr1pTr2combi_unmatched;
    TH1D* h_pTTr1pTTr2combi_unmatched;

    TH1D* h_D0combi_p;
    TH1D* h_D0combi_pT;
    TH1D* h_D0combi_eta;
    TH1D* h_D0combi_phi;
    TH1D* h_D0combi_dRJet;
    TH1D* h_D0combi_Mass;
    TH1D* h_D0combi_paired_p;
    TH1D* h_D0combi_paired_pT;
    TH1D* h_D0combi_paired_eta;
    TH1D* h_D0combi_paired_phi;
    TH1D* h_D0combi_paired_dRJet;
    TH1D* h_D0combi_paired_Mass;
    TH1D* h_D0combi_switched_p;
    TH1D* h_D0combi_switched_pT;
    TH1D* h_D0combi_switched_eta;
    TH1D* h_D0combi_switched_phi;
    TH1D* h_D0combi_switched_dRJet;
    TH1D* h_D0combi_switched_Mass;
    TH1D* h_D0combi_unmatched_p;
    TH1D* h_D0combi_unmatched_pT;
    TH1D* h_D0combi_unmatched_eta;
    TH1D* h_D0combi_unmatched_phi;
    TH1D* h_D0combi_unmatched_dRJet;
    TH1D* h_D0combi_unmatched_Mass;

    TH1D* h_DeltaR_trCand_Mu;
    TH1D* h_DeltaR_trCand_El;
    TH1D* h_BCand_DeltaRD0combiMu;
    TH1D* h_B_DeltaRD0combiMu;
    TH1D* h_B_pD0combipMu;
    TH1D* h_B_pTD0combipTMu;
    TH1D* h_B_D0combiMass;
    TH1D* h_B_combiMass;
    TH1D* h_B_combip;
    TH1D* h_B_combipT;
    TH1D* h_B_combieta;
    TH1D* h_B_combiphi;
    TH1D* h_B_paired_DeltaRD0combiMu;
    TH1D* h_B_paired_pD0combipMu;
    TH1D* h_B_paired_pTD0combipTMu;
    TH1D* h_B_paired_D0combiMass;
    TH1D* h_B_paired_combiMass;
    TH1D* h_B_paired_combip;
    TH1D* h_B_paired_combipT;
    TH1D* h_B_paired_combieta;
    TH1D* h_B_paired_combiphi;
    TH1D* h_B_switched_DeltaRD0combiMu;
    TH1D* h_B_switched_pD0combipMu;
    TH1D* h_B_switched_pTD0combipTMu;
    TH1D* h_B_switched_D0combiMass;
    TH1D* h_B_switched_combiMass;
    TH1D* h_B_switched_combip;
    TH1D* h_B_switched_combipT;
    TH1D* h_B_switched_combieta;
    TH1D* h_B_switched_combiphi;
    TH1D* h_B_unmatched_DeltaRD0combiMu;
    TH1D* h_B_unmatched_pD0combipMu;
    TH1D* h_B_unmatched_pTD0combipTMu;
    TH1D* h_B_unmatched_D0combiMass;
    TH1D* h_B_unmatched_combiMass;
    TH1D* h_B_unmatched_combip;
    TH1D* h_B_unmatched_combipT;
    TH1D* h_B_unmatched_combieta;
    TH1D* h_B_unmatched_combiphi;

    // optimized invariant combination
    TH1D* h_D0optcombiCand_pT;

    TH1D* h_D0optcombi_p;
    TH1D* h_D0optcombi_pT;
    TH1D* h_D0optcombi_eta;
    TH1D* h_D0optcombi_phi;
    TH1D* h_D0optcombi_dRJet;
    TH1D* h_D0optcombi_Mass;
    TH1D* h_D0optcombi_paired_p;
    TH1D* h_D0optcombi_paired_pT;
    TH1D* h_D0optcombi_paired_eta;
    TH1D* h_D0optcombi_paired_phi;
    TH1D* h_D0optcombi_paired_dRJet;
    TH1D* h_D0optcombi_paired_Mass;
    TH1D* h_D0optcombi_switched_p;
    TH1D* h_D0optcombi_switched_pT;
    TH1D* h_D0optcombi_switched_eta;
    TH1D* h_D0optcombi_switched_phi;
    TH1D* h_D0optcombi_switched_dRJet;
    TH1D* h_D0optcombi_switched_Mass;
    TH1D* h_D0optcombi_unmatched_p;
    TH1D* h_D0optcombi_unmatched_pT;
    TH1D* h_D0optcombi_unmatched_eta;
    TH1D* h_D0optcombi_unmatched_phi;
    TH1D* h_D0optcombi_unmatched_dRJet;
    TH1D* h_D0optcombi_unmatched_Mass;

    TH1D* h_BCand_DeltaRD0optcombiMu;
    TH1D* h_B_DeltaRD0optcombiMu;
    TH1D* h_B_pD0optcombipMu;
    TH1D* h_B_pTD0optcombipTMu;
    TH1D* h_B_D0optcombiMass;
    TH1D* h_B_optcombiMass;
    TH1D* h_B_optcombip;
    TH1D* h_B_optcombipT;
    TH1D* h_B_optcombieta;
    TH1D* h_B_optcombiphi;
    TH1D* h_B_paired_DeltaRD0optcombiMu;
    TH1D* h_B_paired_pD0optcombipMu;
    TH1D* h_B_paired_pTD0optcombipTMu;
    TH1D* h_B_paired_D0optcombiMass;
    TH1D* h_B_paired_optcombiMass;
    TH1D* h_B_paired_optcombip;
    TH1D* h_B_paired_optcombipT;
    TH1D* h_B_paired_optcombieta;
    TH1D* h_B_paired_optcombiphi;
    TH1D* h_B_switched_DeltaRD0optcombiMu;
    TH1D* h_B_switched_pD0optcombipMu;
    TH1D* h_B_switched_pTD0optcombipTMu;
    TH1D* h_B_switched_D0optcombiMass;
    TH1D* h_B_switched_optcombiMass;
    TH1D* h_B_switched_optcombip;
    TH1D* h_B_switched_optcombipT;
    TH1D* h_B_switched_optcombieta;
    TH1D* h_B_switched_optcombiphi;
    TH1D* h_B_unmatched_DeltaRD0optcombiMu;
    TH1D* h_B_unmatched_pD0optcombipMu;
    TH1D* h_B_unmatched_pTD0optcombipTMu;
    TH1D* h_B_unmatched_D0optcombiMass;
    TH1D* h_B_unmatched_optcombiMass;
    TH1D* h_B_unmatched_optcombip;
    TH1D* h_B_unmatched_optcombipT;
    TH1D* h_B_unmatched_optcombieta;
    TH1D* h_B_unmatched_optcombiphi;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
KalmanAnalyzer::KalmanAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  nEvts = 0;
  nEvts2 = 0;
  nEvts3 = 0.;

}


KalmanAnalyzer::~KalmanAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "Initial number of events                    = " << nEvts << std::endl;
  std::cout << "Initial number of semileptonic ttbar events = " << nEvts2 << " (" << nEvts3 << ")" << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
KalmanAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  nEvts++;
  hasGoodMuons = false;
  hasGoodElectrons = false;
  hasGoodLeptons = false;
  hasGoodJets = false;
  isGoodEvt = false;
  weight = 1.;

  //-------------------------------------------
  // Apply the recommended ttbar semilept sel
  //-------------------------------------------
  
  // Muon selection
  edm::Handle<std::vector<pat::Muon>>  muonHandle;
  edm::InputTag tagMuon("selectedPatMuonsPFlow");
  iEvent.getByLabel(tagMuon, muonHandle);
  const std::vector<pat::Muon>& muon = *(muonHandle.product());
	std::vector<const pat::Muon*> myMuons;  
  int ngoodmuon = 0;
  int nvetomuon = 0;

  edm::Handle<std::vector<reco::Vertex>> pvHandle;
  iEvent.getByLabel("goodOfflinePrimaryVertices", pvHandle);
 
	for (unsigned int iMuon = 0; iMuon < muon.size(); ++iMuon) {
	  if (muon[iMuon].pt() < 26.) continue;
	  if (fabs(muon[iMuon].eta()) > 2.1) continue;
    if (muon[iMuon].normChi2() >= 10) continue;
    if (muon[iMuon].track()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
    if (muon[iMuon].globalTrack()->hitPattern().numberOfValidMuonHits() <= 0) continue;
    if (muon[iMuon].numberOfMatchedStations() <= 1) continue;
    if (muon[iMuon].dB() >= 0.2) continue;
    if (fabs(muon[iMuon].muonBestTrack()->dz(pvHandle->at(0).position())) >= 0.5) continue;
    if (muon[iMuon].innerTrack()->hitPattern().numberOfValidPixelHits() <= 0) continue; 
	  myMuons.push_back(&muon[iMuon]);
	  ++ngoodmuon;
	}  
  
  for (unsigned int iMuon = 0; iMuon < muon.size(); ++iMuon) {
	  if (muon[iMuon].pt() < 10.) continue;
	  if (fabs(muon[iMuon].eta()) > 2.5) continue;
    ++nvetomuon;
  }
      
  // Electron selection
  edm::Handle<std::vector<pat::Electron>>  electronHandle;
  edm::InputTag tagElectron("selectedPatElectronsPFlow");
  iEvent.getByLabel(tagElectron, electronHandle);
  const std::vector<pat::Electron>& electron = *(electronHandle.product());
	std::vector<const pat::Electron*> myElectrons;  
  int ngoodelectron = 0;
  int nvetoelectron = 0;

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  edm::Handle<reco::BeamSpot> hBeamspot;
  iEvent.getByLabel("offlineBeamSpot", hBeamspot);
  const reco::BeamSpot &beamSpot = *hBeamspot;
  edm::Handle<reco::VertexCollection> hVtx;
  iEvent.getByLabel("goodOfflinePrimaryVertices", hVtx);
  edm::Handle<double> hRhoIso;
  iEvent.getByLabel(edm::InputTag("kt6PFJets", "rho", "RECO"), hRhoIso); 
  double rhoIso = *hRhoIso;
  ElectronEffectiveArea::ElectronEffectiveAreaTarget EAtarget;
  EAtarget  =ElectronEffectiveArea::kEleEANoCorr; // for MC

  for (unsigned int iElectron = 0; iElectron < electron.size(); ++iElectron) {
    if (electron[iElectron].pt() < 30.) continue;
    if (fabs(electron[iElectron].eta()) > 2.5) continue;
    if (!EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, electron[iElectron], hConversions, beamSpot, hVtx, electron[iElectron].chargedHadronIso(), electron[iElectron].photonIso(), electron[iElectron].neutralHadronIso(), rhoIso, EAtarget)) continue;
    if (fabs(electron[iElectron].superCluster()->eta()) >= 1.4442 && fabs(electron[iElectron].superCluster()->eta()) < 1.5660) continue;
	  myElectrons.push_back(&electron[iElectron]);    
    ++ngoodelectron;
  }

  for (unsigned int iElectron = 0; iElectron < electron.size(); ++iElectron) {
    if (electron[iElectron].pt() < 20.) continue;
    if (fabs(electron[iElectron].eta()) > 2.5) continue;
    if (fabs(electron[iElectron].superCluster()->eta()) >= 1.4442 && fabs(electron[iElectron].superCluster()->eta()) < 1.5660) continue;
    ++nvetoelectron;
  }

  // Lepton selection
  if (ngoodmuon == 1 && nvetoelectron == 0) hasGoodMuons = true;
  if (ngoodelectron == 1 && nvetomuon == 0) hasGoodElectrons = true;
  if (hasGoodMuons || hasGoodElectrons) hasGoodLeptons =true;

  // Jets selection
  edm::Handle<std::vector<pat::Jet>> jetHandle;
  edm::InputTag tagJet("selectedPatJetsPFlow","","PAT");
  iEvent.getByLabel(tagJet, jetHandle);
  const std::vector<pat::Jet>& jet = *(jetHandle.product());
	std::vector<const pat::Jet*> myJets;  
  int n55jet = 0;
  int n45jet = 0;
  int n35jet = 0;
  int n20jet = 0;

  for (unsigned int iJet = 0; iJet < jet.size(); ++iJet) {
    if (jet[iJet].pt() < 20.) continue;
    if (fabs(jet[iJet].eta()) > 2.5) continue;
    if (jet[iJet].pt() > 55.) ++n55jet;
    if (jet[iJet].pt() > 45.) ++n45jet;
    if (jet[iJet].pt() > 35.) ++n35jet;
	  myJets.push_back(&jet[iJet]);    
    ++n20jet;
  }

  if (n55jet > 0 && n45jet > 1 && n35jet > 2 && n20jet > 3) hasGoodJets = true;

  if (hasGoodLeptons && hasGoodJets) isGoodEvt = true;

  if (isGoodEvt) {
    nEvts2++;

    // Compute weight
    edm::InputTag vtxTag("offlinePrimaryVertices","","RECO");
    edm::Handle<std::vector<reco::Vertex> > Hvertex;
    iEvent.getByLabel(vtxTag,Hvertex);
    const std::vector<reco::Vertex>& vertex = *(Hvertex.product());
    int nvtx = 0;
    for (unsigned int iVtx=0; iVtx<vertex.size(); ++iVtx ) {
      if (vertex[iVtx].isFake()) continue;
      if (vertex[iVtx].ndof() <= 4) continue;  
      ++nvtx;
    }
    double lumitab[4]= {888.7,4446,7021 ,7221}; //when running on all the stat
    TopTriggerEfficiencyProvider *weight_provider = new TopTriggerEfficiencyProvider(true,lumitab);
    if (hasGoodMuons)
      weight *= weight_provider->get_weight(myMuons[0]->pt(), myMuons[0]->eta(), myJets[3]->pt(), myJets[3]->eta(), nvtx, n20jet,true, TopTriggerEfficiencyProvider::NOMINAL);
    if (hasGoodElectrons)
      weight = weight_provider->get_weight(myElectrons[0]->pt(), myElectrons[0]->eta(), myJets[3]->pt(), myJets[3]->eta(), nvtx, n20jet,false, TopTriggerEfficiencyProvider::NOMINAL);

    double puWeight[41] = {0.343966,0.421778,0.436096,0.244907,0.159864,0.301344,0.241472,0.274829,0.44007,0.758224,1.17234,1.57047,1.65648,1.48344,1.25881,1.09602,1.02284,1.01614,1.05619,1.11854,1.17075,1.1998,1.20717,1.1982,1.17317,1.13123,1.0737,1.00772,0.928615,0.834017,0.723992,0.605741,0.485348,0.371787,0.270933,0.187491,0.124365,0.0791913,0.0484192,0.0288752,0.0127159};
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float m_nTrueInteractions = 0.;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) { 
        m_nTrueInteractions = PVI->getTrueNumInteractions();
        continue;
      }
    }
    int npu = (int)m_nTrueInteractions;
    if ( npu > 40) npu = 40;
    weight *= puWeight[npu];

    nEvts3 = nEvts3 + weight;

    //-------------------------------------------
    // Access the vertex
    //-------------------------------------------

    edm::Handle<reco::VertexCollection>  vtxHandle;
    edm::InputTag tagVtx("offlinePrimaryVertices");
    iEvent.getByLabel(tagVtx, vtxHandle);
    const reco::VertexCollection vtx = *(vtxHandle.product());

    if (vtx.size() == 0) {
      std::cout << " WARNING : no PV for this event ... " << std::endl;
      return;
    }

    //-------------------------------------------
    // Access the Gen D0 -> K Pi
    //-------------------------------------------
  
    edm::Handle<std::vector<reco::GenParticle> > mcHandle;
    edm::InputTag tagGen("genParticles","","SIM");
    iEvent.getByLabel(tagGen,mcHandle);
    const reco::GenParticleCollection& mcparts = *(mcHandle.product());

    // create a vector of pointers to D0/JPSI MC particles
    std::vector<const reco::GenParticle*> genPartD0;

    for (reco::GenParticleCollection::const_iterator it = mcparts.begin(); it != mcparts.end(); ++it)  {
      if (abs((*it).pdgId()) != 421) continue;
      if ((*it).numberOfDaughters() != 2) continue;
      if ((abs((*it).daughter(0)->pdgId()) != 321 || abs((*it).daughter(1)->pdgId()) != 211) 
          && (abs((*it).daughter(0)->pdgId()) != 211 || abs((*it).daughter(1)->pdgId()) != 321)) 
        continue;
      genPartD0.push_back(&(*it));
      h_genD0_p->Fill((*it).p(),weight);
      h_genD0_pT->Fill((*it).pt(),weight);
      h_genD0_eta->Fill((*it).eta(),weight);
      h_genD0_phi->Fill((*it).phi(),weight);
    }  
    h_genD0_n->Fill(genPartD0.size(),weight);

    //--------------------------------------------------
    // Access the PF candidates for non-isolated mu/e
    //--------------------------------------------------

    edm::Handle<reco::PFCandidateCollection>  pfHandle;
    edm::InputTag tagPF("particleFlow","","RECO");
    iEvent.getByLabel(tagPF, pfHandle);
    reco::PFCandidateCollection pfs = *pfHandle;
    if (!pfHandle.isValid()) {
      std::cout << "=> pfHandle is not valid..." << std::endl;
      return;
    }

    // Select good PF muons :

    std::vector<const reco::PFCandidate*> myPFmu;
    std::vector<const reco::PFCandidate*> myPFel;
    for (unsigned int i = 0; i < pfs.size(); ++i) {

      if (pfs[i].pt() < 4.) continue;

      if (abs(pfs[i].pdgId()) == 13) { 
        if (fabs(pfs[i].eta()) < 2.4) myPFmu.push_back(&pfs[i]);
      }
      else {
        if (abs(pfs[i].pdgId()) == 11) {
          if (fabs(pfs[i].eta()) < 2.6) myPFel.push_back(&pfs[i]);
        }
        else continue;
      }
    }

    //-------------------------------------------
    // Access the jets 
    //-------------------------------------------

    //  double d0mass = 1.86484;
    //  double bmass = 5.2796;

    ParticleMass gMassD0 = 1.86483;
    //float      gSigmaD0 = 0.00014;

    //ParticleMass gMassW = 80.399;
    //float        gSigmaW = 0.023;
    //float        gResoW = 10.;

    ParticleMass gMassMu  = 0.105658367;
    //float      gSigmaMu = 0.000000004; 

    ParticleMass gMassK  = 0.493677;
    float        gSigmaK = 0.000001;

    ParticleMass gMassPi  = 0.13957018;
    float        gSigmaPi = 0.00000001;

    // Track setup

    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);  

    edm::InputTag jetTag("selectedPatJetsPFlow","","PAT");

    edm::Handle<pat::JetCollection> jetsHandle;
    iEvent.getByLabel(jetTag, jetsHandle);
    pat::JetCollection jets = *jetsHandle;

    float maxcsv = -1.;
    int maxind = -1;
    float second_max = -1.0;
    int maxind2 = -1;
    int iSelJet = 0;

    for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
      double btag = (*it).bDiscriminator("combinedSecondaryVertexBJetTags");
      if ( (*it).pt() >= 30. ) {
        if(btag >= maxcsv) {
          second_max = maxcsv;
          maxcsv = btag;
          maxind2 = maxind; 
          maxind = iSelJet;
        }
        else if(btag > second_max) {
          second_max = btag;
          maxind2 = iSelJet;
        }
      }
      iSelJet++;
    }
    h_nJets->Fill(iSelJet,weight);

    iSelJet = -1;

    for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
      iSelJet++;

      if (iSelJet == maxind || iSelJet == maxind2){
        h_CSV->Fill((*it).bDiscriminator("combinedSecondaryVertexBJetTags"),weight);

        TLorentzVector p_Jet, p_D0cons;
        p_Jet.SetPtEtaPhiM((*it).pt(), (*it).eta(), (*it).phi(), (*it).mass());
        p_D0cons.SetPtEtaPhiM(0., 0., 0., 0.);
        double D0cons_chi2NDOF = 200.;
        double chargeK = 0;
        double D0cons_ctau[2] = {0., 0.};
        bool isD0consPaired = false;
        bool isD0consSwitched = false;
        bool isD0consUnmatched = false;

        double pt_trCand_D0combi[3]  = {0., 0., 0.};
        double eta_trCand_D0combi[3] = {0., 0., 0.};
        double phi_trCand_D0combi[3] = {0., 0., 0.};
        double id_trCand_D0combi[3] = {0., 0., 0.};

        reco::TrackRefVector jetTracks = (*it).associatedTracks();

        for (reco::track_iterator iter1 = jetTracks.begin(); iter1 != jetTracks.end(); ++iter1) {

          const reco::Track& Track1 = **iter1;

          if ((**iter1).pt() < 4.) continue;
          if (!Track1.quality(reco::Track::highPurity)) continue;

          // store the 3 traces with highest pT, excluding muons and electrons
          bool trCandIsMu = false;
          for (unsigned int iMuCand = 0; iMuCand < myPFmu.size(); iMuCand++) {
            TLorentzVector p_MuCand, p_trCand;
            p_MuCand.SetPtEtaPhiM(myPFmu[iMuCand]->pt(), myPFmu[iMuCand]->eta(), myPFmu[iMuCand]->phi(), gMassMu);
            p_trCand.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassPi);
            h_DeltaR_trCand_Mu->Fill(p_trCand.DeltaR(p_MuCand),weight);
            if (p_trCand.DeltaR(p_MuCand) < 0.0005) {
              trCandIsMu = true;
              break;
            }
          }
          bool trCandIsEl = false;
          for (unsigned int iElCand = 0; iElCand < myPFel.size(); iElCand++) {
            TLorentzVector p_ElCand, p_trCand;
            p_ElCand.SetPtEtaPhiM(myPFel[iElCand]->pt(), myPFel[iElCand]->eta(), myPFel[iElCand]->phi(), 0.);
            p_trCand.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassPi);
            h_DeltaR_trCand_El->Fill(p_trCand.DeltaR(p_ElCand),weight);
            if (p_trCand.DeltaR(p_ElCand) < 0.005) {
              trCandIsEl = true;
              break;
            }
          }
          if (!trCandIsMu && !trCandIsEl) {
            if ((**iter1).pt() >= pt_trCand_D0combi[0]) {
              pt_trCand_D0combi[2]  = pt_trCand_D0combi[1];
              eta_trCand_D0combi[2] = eta_trCand_D0combi[1];
              phi_trCand_D0combi[2] = phi_trCand_D0combi[1];
              id_trCand_D0combi[2]  = id_trCand_D0combi[1];
              pt_trCand_D0combi[1]  = pt_trCand_D0combi[0];
              eta_trCand_D0combi[1] = eta_trCand_D0combi[0];
              phi_trCand_D0combi[1] = phi_trCand_D0combi[0];
              id_trCand_D0combi[1]  = id_trCand_D0combi[0];
              pt_trCand_D0combi[0]  = (**iter1).pt();
              eta_trCand_D0combi[0] = (**iter1).eta();
              phi_trCand_D0combi[0] = (**iter1).phi();
              id_trCand_D0combi[0]  = (**iter1).charge();
            }
            else {
              if ((**iter1).pt() >= pt_trCand_D0combi[1]) {
                pt_trCand_D0combi[2]  = pt_trCand_D0combi[1];
                eta_trCand_D0combi[2] = eta_trCand_D0combi[1];
                phi_trCand_D0combi[2] = phi_trCand_D0combi[1];
                id_trCand_D0combi[2]  = id_trCand_D0combi[1];
                pt_trCand_D0combi[1]  = (**iter1).pt();
                eta_trCand_D0combi[1] = (**iter1).eta();
                phi_trCand_D0combi[1] = (**iter1).phi();
                id_trCand_D0combi[1]  = (**iter1).charge();
              }
              else  
                if ((**iter1).pt() >= pt_trCand_D0combi[2]) {
                  pt_trCand_D0combi[2]  = (**iter1).pt();
                  eta_trCand_D0combi[2] = (**iter1).eta();
                  phi_trCand_D0combi[2] = (**iter1).phi();
                  id_trCand_D0combi[2]  = (**iter1).charge();
                }
            }
          }

          reco::TransientTrack tr1 = (*theB).build((**iter1));

          for (reco::track_iterator iter2 = jetTracks.begin(); iter2 != jetTracks.end(); ++iter2) {
            const reco::Track& Track2 = **iter2;

            if (iter2 == iter1) continue;
            if ((**iter2).pt() < 4.) continue;
            if (!Track2.quality(reco::Track::highPurity)) continue;

            int iBCut = 0;
            h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
            h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"2 tracks with p_{T} > 4 GeV/c");

            reco::TransientTrack tr2 = (*theB).build((**iter2));

            //============================
            // simple Kalman Vertex Fitter
            //============================

            //~~~~~~~~~~~~~~~~~~~~~~~~~~
            // reconstruct D^0 -> K Pi
            //~~~~~~~~~~~~~~~~~~~~~~~~~~

            // Select OS tracks
            if ((**iter1).charge()*(**iter2).charge() > 0) continue;
            h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
            h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... of opposite signs");

            // Compute the mass
            TLorentzVector p_tr1_D0, p_tr2_D0, p_D0;
            p_tr1_D0.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassK);
            p_tr2_D0.SetPtEtaPhiM((**iter2).pt(), (**iter2).eta(), (**iter2).phi(), gMassPi);

            bool tr1IsMu = false;
            bool tr2IsMu = false;
            for (unsigned int iMuCand = 0; iMuCand < myPFmu.size(); iMuCand++) {
              TLorentzVector p_MuCand;
              p_MuCand.SetPtEtaPhiM(myPFmu[iMuCand]->pt(), myPFmu[iMuCand]->eta(), myPFmu[iMuCand]->phi(), gMassMu);
              if (p_tr1_D0.DeltaR(p_MuCand) < 0.0005) tr1IsMu = true;
              if (p_tr2_D0.DeltaR(p_MuCand) < 0.0005) tr2IsMu = true;
              if (tr1IsMu && tr2IsMu) break;
            }
            if (tr1IsMu || tr2IsMu) continue;
            h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
            h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... not identified as #mu");
            bool tr1IsEl = false;
            bool tr2IsEl = false;
            for (unsigned int iElCand = 0; iElCand < myPFel.size(); iElCand++) {
              TLorentzVector p_ElCand;
              p_ElCand.SetPtEtaPhiM(myPFel[iElCand]->pt(), myPFel[iElCand]->eta(), myPFel[iElCand]->phi(), 0.);
              if (p_tr1_D0.DeltaR(p_ElCand) < 0.005) tr1IsEl = true;
              if (p_tr2_D0.DeltaR(p_ElCand) < 0.005) tr2IsEl = true;
              if (tr1IsEl && tr2IsEl) break;
            }
            if (tr1IsEl || tr2IsEl) continue;
            h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
            h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... not identified as e");

            // min DR between the kaon and the muon
            /*
               if (p_tr1_D0.DeltaR(p_tr2_D0) > 0.2) continue;
               h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
               h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... within #DeltaR < 0.2");
               */

            p_D0 = p_tr1_D0 + p_tr2_D0;

            bool isD0paired = false;
            bool isD0switched = false;
            bool isD0unmatched = false;
            int iGenD0 = -1;
            double GenRecoD0dR = 200.;
            double GenRecoTr1dR = 200.;
            double GenRecoTr2dR = 200.;

            for (unsigned int iGen = 0; iGen < genPartD0.size(); iGen++) {
              TLorentzVector p_genD0;
              p_genD0.SetPtEtaPhiM(genPartD0[iGen]->pt(), genPartD0[iGen]->eta(), genPartD0[iGen]->phi(), genPartD0[iGen]->mass());
              double GenRecoD0dRtmp = p_D0.DeltaR(p_genD0);
              if (GenRecoD0dRtmp < GenRecoD0dR) {
                GenRecoD0dR = GenRecoD0dRtmp;
                iGenD0 = iGen;
              }
            }

            if (fabs(GenRecoD0dR - 200.) > 1e-10) h_candGenRecoD0dR->Fill(GenRecoD0dR,weight);
            if (GenRecoD0dR > 0.5) 
              isD0unmatched = true;
            else {
              TLorentzVector p_genK, p_genPi;
              if (abs(genPartD0[iGenD0]->daughter(0)->pdgId()) == 321) {
                p_genK.SetPtEtaPhiM(genPartD0[iGenD0]->daughter(0)->pt(), genPartD0[iGenD0]->daughter(0)->eta(), genPartD0[iGenD0]->daughter(0)->phi(), genPartD0[iGenD0]->daughter(0)->mass());
                p_genPi.SetPtEtaPhiM(genPartD0[iGenD0]->daughter(1)->pt(), genPartD0[iGenD0]->daughter(1)->eta(), genPartD0[iGenD0]->daughter(1)->phi(), genPartD0[iGenD0]->daughter(1)->mass());
              }
              else {
                p_genK.SetPtEtaPhiM(genPartD0[iGenD0]->daughter(1)->pt(), genPartD0[iGenD0]->daughter(1)->eta(), genPartD0[iGenD0]->daughter(1)->phi(), genPartD0[iGenD0]->daughter(1)->mass());
                p_genPi.SetPtEtaPhiM(genPartD0[iGenD0]->daughter(0)->pt(), genPartD0[iGenD0]->daughter(0)->eta(), genPartD0[iGenD0]->daughter(0)->phi(), genPartD0[iGenD0]->daughter(0)->mass());
              }
              double GenKRecoTr1dR = p_tr1_D0.DeltaR(p_genK);
              double GenPiRecoTr1dR = p_tr1_D0.DeltaR(p_genPi);
              double GenKRecoTr2dR = p_tr2_D0.DeltaR(p_genK);
              double GenPiRecoTr2dR = p_tr2_D0.DeltaR(p_genPi);
              if (GenKRecoTr1dR < GenPiRecoTr1dR) GenRecoTr1dR = GenKRecoTr1dR;
              else GenRecoTr1dR = GenPiRecoTr1dR;
              if (GenPiRecoTr2dR < GenKRecoTr2dR) GenRecoTr2dR = GenPiRecoTr2dR;
              else GenRecoTr2dR = GenKRecoTr2dR;
              h_candGenRecoTrdR->Fill(GenRecoTr1dR,weight);
              h_candGenRecoTrdR->Fill(GenRecoTr2dR,weight);
              if (GenRecoTr1dR > 0.005 || GenRecoTr2dR > 0.005) 
                isD0unmatched = true;
              else {
                if (GenKRecoTr1dR < GenPiRecoTr1dR && GenPiRecoTr2dR < GenKRecoTr2dR) isD0paired = true;
                else isD0switched = true;
              }
            }

            //Creating a KinematicParticleFactory
            KinematicParticleFactoryFromTransientTrack pFactory;

            //initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
            float chi_D0 = 0.;
            float ndf_D0 = 0.;

            //making particles
            std::vector<RefCountedKinematicParticle> D0Particles;
            D0Particles.push_back(pFactory.particle (tr1,gMassK,chi_D0,ndf_D0,gSigmaK));
            D0Particles.push_back(pFactory.particle (tr2,gMassPi,chi_D0,ndf_D0,gSigmaPi));

            /* Example of a simple vertex fit, without other constraints
             * The reconstructed decay tree is a result of the kinematic fit
             * The KinematicParticleVertexFitter fits the final state particles to their vertex and
             * reconstructs the decayed state
             */

            // creating the vertex fitter
            KinematicParticleVertexFitter D0fitter;
            RefCountedKinematicTree D0vertexFitTree = D0fitter.fit(D0Particles);

            if (D0vertexFitTree->isValid()) {

              //accessing the tree components, move pointer to top
              D0vertexFitTree->movePointerToTheTop();

              //We are now at the top of the decay tree getting the d0 reconstructed KinematicPartlcle
              RefCountedKinematicParticle D0 = D0vertexFitTree->currentParticle();
              RefCountedKinematicVertex D0_vertex = D0vertexFitTree->currentDecayVertex();

              if (D0_vertex->vertexIsValid()) {

                h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
                h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... with a valid D^{0} vertex");

                h_D0Cand_Chi2NDOF->Fill(D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom(),weight);

                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 1.) 
                  h_D0Cand_MassChi2Inf1->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 1.5) 
                  h_D0Cand_MassChi2Inf1p5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 2.)
                  h_D0Cand_MassChi2Inf2->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 2.5)
                  h_D0Cand_MassChi2Inf2p5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 3.)
                  h_D0Cand_MassChi2Inf3->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 3.5)
                  h_D0Cand_MassChi2Inf3p5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 4.)
                  h_D0Cand_MassChi2Inf4->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 4.5)
                  h_D0Cand_MassChi2Inf4p5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 5.)
                  h_D0Cand_MassChi2Inf5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 5.5)
                  h_D0Cand_MassChi2Inf5p5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 6.)
                  h_D0Cand_MassChi2Inf6->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 6.5)
                  h_D0Cand_MassChi2Inf6p5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 7.)
                  h_D0Cand_MassChi2Inf7->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 7.5)
                  h_D0Cand_MassChi2Inf7p5->Fill(p_D0.M(),weight);
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 8.)
                  h_D0Cand_MassChi2Inf8->Fill(p_D0.M(),weight);

                // cut on chi2/NDOF
                if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 4.) {
                  h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
                  h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... with #chi^{2}/NDOF < 4");

                  // Distance to PV :
                  GlobalPoint D0_svPos    = D0_vertex->position();
                  GlobalError D0_svPosErr = D0_vertex->error();

                  double sigmax_vtx_D0vtx = sqrt(pow(vtx[0].xError(), 2.) + pow(D0_svPosErr.cxx(), 2.));
                  double sigmay_vtx_D0vtx = sqrt(pow(vtx[0].yError(), 2.) + pow(D0_svPosErr.cyy(), 2.));
                  double sigmaz_vtx_D0vtx = sqrt(pow(vtx[0].zError(), 2.) + pow(D0_svPosErr.czz(), 2.));

                  double D0_interx = pow((p_D0.Px()/p_D0.M())/sigmax_vtx_D0vtx, 2.);
                  double D0_intery = pow((p_D0.Py()/p_D0.M())/sigmay_vtx_D0vtx, 2.);
                  double D0_interz = pow((p_D0.Pz()/p_D0.M())/sigmaz_vtx_D0vtx, 2.);

                  double D0_sigmaL3D = pow(D0_interx + D0_intery + D0_interz, -0.5);

                  double D0_part1 = (p_D0.Px()/p_D0.M())*pow(D0_sigmaL3D/sigmax_vtx_D0vtx,2.)*( D0_svPos.x() - vtx[0].x());
                  double D0_part2 = (p_D0.Py()/p_D0.M())*pow(D0_sigmaL3D/sigmay_vtx_D0vtx,2.)*( D0_svPos.y() - vtx[0].y());
                  double D0_part3 = (p_D0.Pz()/p_D0.M())*pow(D0_sigmaL3D/sigmaz_vtx_D0vtx,2.)*( D0_svPos.z() - vtx[0].z());

                  double D0_L3D = fabs(D0_part1 + D0_part2 + D0_part3);

                  double D0_L3DoverSigmaL3D = D0_L3D/D0_sigmaL3D;

                  h_D0Cand_L->Fill(D0_L3D,weight);
                  h_D0Cand_SigmaL->Fill(D0_sigmaL3D,weight);
                  h_D0Cand_LOverSigmaL->Fill(D0_L3DoverSigmaL3D,weight);

                  // cut on L/SigmaL
                  //if (D0_L3DoverSigmaL3D > 50.) {
                  if (D0_L3DoverSigmaL3D > 100.) {
                    h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
                    h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... with c#tau/#sigma(c#tau) > 100");
                    //h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... with c#tau/#sigma(c#tau) > 50");

                    h_D0Cand_pT->Fill(p_D0.Pt(),weight);

                    // cut on pT
                    if (p_D0.Pt() > 12.) {
                      h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
                      h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... with p_{T} > 12 GeV/c");

                      h_D0_L->Fill(D0_L3D,weight);
                      h_D0_SigmaL->Fill(D0_sigmaL3D,weight);
                      h_D0_dRJet->Fill(p_D0.DeltaR(p_Jet),weight);
                      h_D0_Mass->Fill(p_D0.M(),weight);
                      h_D0_p->Fill(p_D0.P(),weight);
                      h_D0_pT->Fill(p_D0.Pt(),weight);
                      h_D0_eta->Fill(p_D0.Eta(),weight);
                      h_D0_phi->Fill(p_D0.Phi(),weight);
                      if (isD0paired) {
                        h_D0_paired_L->Fill(D0_L3D,weight);
                        h_D0_paired_SigmaL->Fill(D0_sigmaL3D,weight);
                        h_D0_paired_dRJet->Fill(p_D0.DeltaR(p_Jet),weight);
                        h_D0_paired_Mass->Fill(p_D0.M(),weight);
                        h_D0_paired_p->Fill(p_D0.P(),weight);
                        h_D0_paired_pT->Fill(p_D0.Pt(),weight);
                        h_D0_paired_eta->Fill(p_D0.Eta(),weight);
                        h_D0_paired_phi->Fill(p_D0.Phi(),weight);
                      }
                      if (isD0switched) {
                        h_D0_switched_L->Fill(D0_L3D,weight);
                        h_D0_switched_SigmaL->Fill(D0_sigmaL3D,weight);
                        h_D0_switched_dRJet->Fill(p_D0.DeltaR(p_Jet),weight);
                        h_D0_switched_Mass->Fill(p_D0.M(),weight);
                        h_D0_switched_p->Fill(p_D0.P(),weight);
                        h_D0_switched_pT->Fill(p_D0.Pt(),weight);
                        h_D0_switched_eta->Fill(p_D0.Eta(),weight);
                        h_D0_switched_phi->Fill(p_D0.Phi(),weight);
                      }
                      if (isD0unmatched) {
                        h_D0_unmatched_L->Fill(D0_L3D,weight);
                        h_D0_unmatched_SigmaL->Fill(D0_sigmaL3D,weight);
                        h_D0_unmatched_dRJet->Fill(p_D0.DeltaR(p_Jet),weight);
                        h_D0_unmatched_Mass->Fill(p_D0.M(),weight);
                        h_D0_unmatched_p->Fill(p_D0.P(),weight);
                        h_D0_unmatched_pT->Fill(p_D0.Pt(),weight);
                        h_D0_unmatched_eta->Fill(p_D0.Eta(),weight);
                        h_D0_unmatched_phi->Fill(p_D0.Phi(),weight);
                      }

                      //~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      // associate D^0 to a PF muon
                      //~~~~~~~~~~~~~~~~~~~~~~~~~~
                      TLorentzVector p_Mu, p_WGaus, p_WDirac;
                      int iMu = -1;
                      double deltaRD0Mu = 2000.;

                      // find closest muon with opposite charged wrt the kaon
                      for (unsigned int iMuCand = 0; iMuCand < myPFmu.size(); iMuCand++) {
                        if (myPFmu[iMuCand]->pdgId()*(**iter1).charge() > 0) continue;
                        TLorentzVector p_MuCand;
                        p_MuCand.SetPtEtaPhiM(myPFmu[iMuCand]->pt(), myPFmu[iMuCand]->eta(), myPFmu[iMuCand]->phi(), gMassMu);
                        double tmp_deltaRD0Mu = p_D0.DeltaR(p_MuCand);
                        if (tmp_deltaRD0Mu < deltaRD0Mu) {
                          deltaRD0Mu = tmp_deltaRD0Mu;
                          iMu = iMuCand;
                        }
                      }

                      if (iMu >= 0) {
                        h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
                        h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... at least a non isolated #mu");
                        h_BCand_DeltaRD0Mu->Fill(deltaRD0Mu,weight);

                        // keep going if closest muon is close enough
                        if (deltaRD0Mu < 0.4) {
                          h_B_cuts->Fill((double)iBCut,weight); ++iBCut;
                          h_B_cuts->GetXaxis()->SetBinLabel(iBCut,"... within #DeltaR < 0.4");
                          h_B_DeltaRD0Mu->Fill(deltaRD0Mu,weight);

                          h_B_D0Mass->Fill(p_D0.M(),weight);

                          p_Mu.SetPtEtaPhiM(myPFmu[iMu]->pt(), myPFmu[iMu]->eta(), myPFmu[iMu]->phi(), gMassMu);
                          h_B_pD0pMu->Fill(p_D0.P() / myPFmu[iMu]->p(),weight);
                          h_B_pTD0pTMu->Fill(p_D0.Pt() / myPFmu[iMu]->pt(),weight);

                          TLorentzVector p_B = p_Mu + p_D0;
                          h_B_Mass->Fill(p_B.M(),weight);
                          h_B_p->Fill(p_B.P(),weight);
                          h_B_pT->Fill(p_B.Pt(),weight);
                          h_B_eta->Fill(p_B.Eta(),weight);
                          h_B_phi->Fill(p_B.Phi(),weight);

                          if (isD0paired) {
                            h_B_paired_DeltaRD0Mu->Fill(deltaRD0Mu,weight);
                            h_B_paired_D0Mass->Fill(p_D0.M(),weight);
                            h_B_paired_pD0pMu->Fill(p_D0.P() / myPFmu[iMu]->p(),weight);
                            h_B_paired_pTD0pTMu->Fill(p_D0.Pt() / myPFmu[iMu]->pt(),weight);
                            h_B_paired_Mass->Fill(p_B.M(),weight);
                            h_B_paired_p->Fill(p_B.P(),weight);
                            h_B_paired_pT->Fill(p_B.Pt(),weight);
                            h_B_paired_eta->Fill(p_B.Eta(),weight);
                            h_B_paired_phi->Fill(p_B.Phi(),weight);
                          }
                          if (isD0switched) {
                            h_B_switched_DeltaRD0Mu->Fill(deltaRD0Mu,weight);
                            h_B_switched_D0Mass->Fill(p_D0.M(),weight);
                            h_B_switched_pD0pMu->Fill(p_D0.P() / myPFmu[iMu]->p(),weight);
                            h_B_switched_pTD0pTMu->Fill(p_D0.Pt() / myPFmu[iMu]->pt(),weight);
                            h_B_switched_Mass->Fill(p_B.M(),weight);
                            h_B_switched_p->Fill(p_B.P(),weight);
                            h_B_switched_pT->Fill(p_B.Pt(),weight);
                            h_B_switched_eta->Fill(p_B.Eta(),weight);
                            h_B_switched_phi->Fill(p_B.Phi(),weight);
                          }
                          if (isD0unmatched) {
                            h_B_unmatched_DeltaRD0Mu->Fill(deltaRD0Mu,weight);
                            h_B_unmatched_D0Mass->Fill(p_D0.M(),weight);
                            h_B_unmatched_pD0pMu->Fill(p_D0.P() / myPFmu[iMu]->p(),weight);
                            h_B_unmatched_pTD0pTMu->Fill(p_D0.Pt() / myPFmu[iMu]->pt(),weight);
                            h_B_unmatched_Mass->Fill(p_B.M(),weight);
                            h_B_unmatched_p->Fill(p_B.P(),weight);
                            h_B_unmatched_pT->Fill(p_B.Pt(),weight);
                            h_B_unmatched_eta->Fill(p_B.Eta(),weight);
                            h_B_unmatched_phi->Fill(p_B.Phi(),weight);
                          }

                        } // non iso mu close enough
                      } // non iso mu
                    } // pT cut
                  } // L3D/sigma cut
                } // chi2/NDOF cut
              } // D0_vertex is valid
            } // Fit is valid

            //=================================
            // constrained Kalman Vertex Fitter
            //=================================

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // reconstruct D^0 to K Pi
            //~~~~~~~~~~~~~~~~~~~~~~~~~~
            //creating the two track mass constraint
            MultiTrackKinematicConstraint* D0cons = new TwoTrackMassKinematicConstraint(gMassD0);

            //creating the fitter
            KinematicConstrainedVertexFitter D0consFitter;

            //obtaining the resulting tree
            RefCountedKinematicTree D0vertexConsFitTree = D0consFitter.fit(D0Particles, D0cons);

            if (D0vertexConsFitTree->isValid()) {

              D0vertexConsFitTree->movePointerToTheTop();
              RefCountedKinematicParticle D0cons = D0vertexConsFitTree->currentParticle();
              RefCountedKinematicVertex D0cons_vertex = D0vertexConsFitTree->currentDecayVertex();

              if (D0cons_vertex->vertexIsValid()) {

                // Distance to PV :
                GlobalPoint D0cons_svPos    = D0cons_vertex->position();
                GlobalError D0cons_svPosErr = D0cons_vertex->error();

                double sigmax_vtx_D0consvtx = sqrt(pow(vtx[0].xError(), 2.) + pow(D0cons_svPosErr.cxx(), 2.));
                double sigmay_vtx_D0consvtx = sqrt(pow(vtx[0].yError(), 2.) + pow(D0cons_svPosErr.cyy(), 2.));
                double sigmaz_vtx_D0consvtx = sqrt(pow(vtx[0].zError(), 2.) + pow(D0cons_svPosErr.czz(), 2.));

                double D0cons_interx = pow((p_D0.Px()/p_D0.M())/sigmax_vtx_D0consvtx, 2.);
                double D0cons_intery = pow((p_D0.Py()/p_D0.M())/sigmay_vtx_D0consvtx, 2.);
                double D0cons_interz = pow((p_D0.Pz()/p_D0.M())/sigmaz_vtx_D0consvtx, 2.);

                double D0cons_sigmaL3D = pow(D0cons_interx + D0cons_intery + D0cons_interz, -0.5);

                double D0cons_part1 = (p_D0.Px()/p_D0.M())*pow(D0cons_sigmaL3D/sigmax_vtx_D0consvtx,2.)*( D0cons_svPos.x() - vtx[0].x());
                double D0cons_part2 = (p_D0.Py()/p_D0.M())*pow(D0cons_sigmaL3D/sigmay_vtx_D0consvtx,2.)*( D0cons_svPos.y() - vtx[0].y());
                double D0cons_part3 = (p_D0.Pz()/p_D0.M())*pow(D0cons_sigmaL3D/sigmaz_vtx_D0consvtx,2.)*( D0cons_svPos.z() - vtx[0].z());

                double D0cons_L3D = fabs(D0cons_part1 + D0cons_part2 + D0cons_part3);

                double D0cons_L3DoverSigmaL3D = D0cons_L3D/D0cons_sigmaL3D;

                h_D0consCand_L->Fill(D0cons_L3D,weight);
                h_D0consCand_SigmaL->Fill(D0cons_sigmaL3D,weight);
                h_D0consCand_LOverSigmaL->Fill(D0cons_L3DoverSigmaL3D,weight);

                // cut on L/SigmaL
                //if (D0cons_L3DoverSigmaL3D > 50.) {
                if (D0cons_L3DoverSigmaL3D > 100.) {

                  h_D0consCand_pT->Fill(p_D0.Pt(),weight);

                  // cut on pT
                  if (p_D0.Pt() > 12.) {

                    if (D0cons_vertex->chiSquared()/(double)D0cons_vertex->degreesOfFreedom() < D0cons_chi2NDOF) {
                      D0cons_chi2NDOF = D0cons_vertex->chiSquared()/(double)D0cons_vertex->degreesOfFreedom();
                      p_D0cons.SetPtEtaPhiM(p_D0.Pt(), p_D0.Eta(), p_D0.Phi(), p_D0.M());
                      D0cons_ctau[0] = D0cons_L3D;
                      D0cons_ctau[1] = D0cons_sigmaL3D;
                      chargeK = (**iter1).charge();
                      isD0consPaired = isD0paired;
                      isD0consSwitched = isD0switched;
                      isD0consUnmatched = isD0unmatched;
                    }

                  }  // pT cut
                } // L3D/sigma cut
              } // vertex is valid
            } // fit is valid

          } // 2nd jet's track loop
        } // 1st jet's track loop

        if (fabs(p_D0cons.Pt()) > 0.) {
          h_D0cons_Chi2NDOF->Fill(D0cons_chi2NDOF,weight);
          h_D0cons_L->Fill(D0cons_ctau[0],weight); 
          h_D0cons_SigmaL->Fill(D0cons_ctau[1],weight); 
          h_D0cons_dRJet->Fill(p_D0cons.DeltaR(p_Jet),weight);
          h_D0cons_Mass->Fill(p_D0cons.M(),weight);
          h_D0cons_p->Fill(p_D0cons.P(),weight);
          h_D0cons_pT->Fill(p_D0cons.Pt(),weight);
          h_D0cons_eta->Fill(p_D0cons.Eta(),weight);
          h_D0cons_phi->Fill(p_D0cons.Phi(),weight);
          if (isD0consPaired) {
            h_D0cons_paired_Chi2NDOF->Fill(D0cons_chi2NDOF,weight);
            h_D0cons_paired_L->Fill(D0cons_ctau[0],weight); 
            h_D0cons_paired_SigmaL->Fill(D0cons_ctau[1],weight); 
            h_D0cons_paired_dRJet->Fill(p_D0cons.DeltaR(p_Jet),weight);
            h_D0cons_paired_Mass->Fill(p_D0cons.M(),weight);
            h_D0cons_paired_p->Fill(p_D0cons.P(),weight);
            h_D0cons_paired_pT->Fill(p_D0cons.Pt(),weight);
            h_D0cons_paired_eta->Fill(p_D0cons.Eta(),weight);
            h_D0cons_paired_phi->Fill(p_D0cons.Phi(),weight);
          }
          if (isD0consSwitched) {
            h_D0cons_switched_Chi2NDOF->Fill(D0cons_chi2NDOF,weight);
            h_D0cons_switched_L->Fill(D0cons_ctau[0],weight); 
            h_D0cons_switched_SigmaL->Fill(D0cons_ctau[1],weight); 
            h_D0cons_switched_dRJet->Fill(p_D0cons.DeltaR(p_Jet),weight);
            h_D0cons_switched_Mass->Fill(p_D0cons.M(),weight);
            h_D0cons_switched_p->Fill(p_D0cons.P(),weight);
            h_D0cons_switched_pT->Fill(p_D0cons.Pt(),weight);
            h_D0cons_switched_eta->Fill(p_D0cons.Eta(),weight);
            h_D0cons_switched_phi->Fill(p_D0cons.Phi(),weight);
          }
          if (isD0consUnmatched) {
            h_D0cons_unmatched_Chi2NDOF->Fill(D0cons_chi2NDOF,weight);
            h_D0cons_unmatched_L->Fill(D0cons_ctau[0],weight); 
            h_D0cons_unmatched_SigmaL->Fill(D0cons_ctau[1],weight); 
            h_D0cons_unmatched_dRJet->Fill(p_D0cons.DeltaR(p_Jet),weight);
            h_D0cons_unmatched_Mass->Fill(p_D0cons.M(),weight);
            h_D0cons_unmatched_p->Fill(p_D0cons.P(),weight);
            h_D0cons_unmatched_pT->Fill(p_D0cons.Pt(),weight);
            h_D0cons_unmatched_eta->Fill(p_D0cons.Eta(),weight);
            h_D0cons_unmatched_phi->Fill(p_D0cons.Phi(),weight);
          }

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~
          // associate D^0 to a PF muon
          //~~~~~~~~~~~~~~~~~~~~~~~~~~
          TLorentzVector p_MuCons;
          int iMuCons = -1;
          double deltaRD0consMu = 2000.;

          // find closest muon with opposite charged wrt the kaon
          for (unsigned int iMuConsCand = 0; iMuConsCand < myPFmu.size(); iMuConsCand++) {
            if (myPFmu[iMuConsCand]->pdgId()*chargeK > 0) continue;
            TLorentzVector p_MuConsCand;
            p_MuConsCand.SetPtEtaPhiM(myPFmu[iMuConsCand]->pt(), myPFmu[iMuConsCand]->eta(), myPFmu[iMuConsCand]->phi(), gMassMu);
            double tmp_deltaRD0consMu = p_D0cons.DeltaR(p_MuConsCand);
            if (tmp_deltaRD0consMu < deltaRD0consMu) {
              deltaRD0consMu = tmp_deltaRD0consMu;
              iMuCons = iMuConsCand;
            }
          }

          if (iMuCons >= 0) {
            h_BCand_DeltaRD0consMu->Fill(deltaRD0consMu,weight);

            // keep going if closest muon is close enough
            if (deltaRD0consMu < 0.4) {
              h_B_DeltaRD0consMu->Fill(deltaRD0consMu,weight);

              h_B_D0consMass->Fill(p_D0cons.M(),weight);

              p_MuCons.SetPtEtaPhiM(myPFmu[iMuCons]->pt(), myPFmu[iMuCons]->eta(), myPFmu[iMuCons]->phi(), gMassMu);
              h_B_pD0conspMu->Fill(p_D0cons.P() / myPFmu[iMuCons]->p(),weight);
              h_B_pTD0conspTMu->Fill(p_D0cons.Pt() / myPFmu[iMuCons]->pt(),weight);

              TLorentzVector p_BCons = p_MuCons + p_D0cons;
              h_B_consMass->Fill(p_BCons.M(),weight);
              h_B_consp->Fill(p_BCons.P(),weight);
              h_B_conspT->Fill(p_BCons.Pt(),weight);
              h_B_conseta->Fill(p_BCons.Eta(),weight);
              h_B_consphi->Fill(p_BCons.Phi(),weight);
              if (isD0consPaired) {
                h_B_paired_DeltaRD0consMu->Fill(deltaRD0consMu,weight);
                h_B_paired_D0consMass->Fill(p_D0cons.M(),weight);
                h_B_paired_pD0conspMu->Fill(p_D0cons.P() / myPFmu[iMuCons]->p(),weight);
                h_B_paired_pTD0conspTMu->Fill(p_D0cons.Pt() / myPFmu[iMuCons]->pt(),weight);
                h_B_paired_consMass->Fill(p_BCons.M(),weight);
                h_B_paired_consp->Fill(p_BCons.P(),weight);
                h_B_paired_conspT->Fill(p_BCons.Pt(),weight);
                h_B_paired_conseta->Fill(p_BCons.Eta(),weight);
                h_B_paired_consphi->Fill(p_BCons.Phi(),weight);
              }
              if (isD0consSwitched) {
                h_B_switched_DeltaRD0consMu->Fill(deltaRD0consMu,weight);
                h_B_switched_D0consMass->Fill(p_D0cons.M(),weight);
                h_B_switched_pD0conspMu->Fill(p_D0cons.P() / myPFmu[iMuCons]->p(),weight);
                h_B_switched_pTD0conspTMu->Fill(p_D0cons.Pt() / myPFmu[iMuCons]->pt(),weight);
                h_B_switched_consMass->Fill(p_BCons.M(),weight);
                h_B_switched_consp->Fill(p_BCons.P(),weight);
                h_B_switched_conspT->Fill(p_BCons.Pt(),weight);
                h_B_switched_conseta->Fill(p_BCons.Eta(),weight);
                h_B_switched_consphi->Fill(p_BCons.Phi(),weight);
              }
              if (isD0consUnmatched) {
                h_B_unmatched_DeltaRD0consMu->Fill(deltaRD0consMu,weight);
                h_B_unmatched_D0consMass->Fill(p_D0cons.M(),weight);
                h_B_unmatched_pD0conspMu->Fill(p_D0cons.P() / myPFmu[iMuCons]->p(),weight);
                h_B_unmatched_pTD0conspTMu->Fill(p_D0cons.Pt() / myPFmu[iMuCons]->pt(),weight);
                h_B_unmatched_consMass->Fill(p_BCons.M(),weight);
                h_B_unmatched_consp->Fill(p_BCons.P(),weight);
                h_B_unmatched_conspT->Fill(p_BCons.Pt(),weight);
                h_B_unmatched_conseta->Fill(p_BCons.Eta(),weight);
                h_B_unmatched_consphi->Fill(p_BCons.Phi(),weight);
              }
            } // non iso mu close enough
          } // non iso mu
        }

        //=================================
        // simple invariant combination
        //=================================

        int p1[6] = {0, 0, 1, 1, 2, 2};
        int p2[6] = {1, 2, 2, 0, 0, 1};

        TLorentzVector p_track1_D0combi, p_track2_D0combi, p_D0combi, p_D0optcombi;
        p_D0optcombi.SetPtEtaPhiM(0., 0., 0., 200.);
        int tk2charge = 0;
        bool isD0optcombiPaired = false;
        bool isD0optcombiSwitched = false;
        bool isD0optcombiUnmatched = false;

        if (fabs(pt_trCand_D0combi[0]) > 1e-10 && fabs(pt_trCand_D0combi[1]) > 1e-10 && fabs(pt_trCand_D0combi[2]) > 1e-10) {
          for (unsigned int iD0combi = 0; iD0combi < 6; iD0combi++) {

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // reconstruct D^0 to K Pi
            //~~~~~~~~~~~~~~~~~~~~~~~~~~
            int tk1 = p1[iD0combi];
            int tk2 = p2[iD0combi];

            // Opposite sign
            if (id_trCand_D0combi[tk1]*id_trCand_D0combi[tk2] > 0) continue;

            p_track1_D0combi.SetPtEtaPhiM(pt_trCand_D0combi[tk1], eta_trCand_D0combi[tk1], phi_trCand_D0combi[tk1], gMassPi);
            p_track2_D0combi.SetPtEtaPhiM(pt_trCand_D0combi[tk2], eta_trCand_D0combi[tk2], phi_trCand_D0combi[tk2], gMassK);
            p_D0combi = p_track1_D0combi + p_track2_D0combi;

            bool isD0combiPaired = false;
            bool isD0combiSwitched = false;
            bool isD0combiUnmatched = false;
            int iGenD0combi = -1;
            double GenRecoD0combidR = 200.;
            double GenRecoTr1combidR = 200.;
            double GenRecoTr2combidR = 200.;

            for (unsigned int iGen = 0; iGen < genPartD0.size(); iGen++) {
              TLorentzVector p_genD0;
              p_genD0.SetPtEtaPhiM(genPartD0[iGen]->pt(), genPartD0[iGen]->eta(), genPartD0[iGen]->phi(), genPartD0[iGen]->mass());
              double GenRecoD0combidRtmp = p_D0combi.DeltaR(p_genD0);
              if (GenRecoD0combidRtmp < GenRecoD0combidR) {
                GenRecoD0combidR = GenRecoD0combidRtmp;
                iGenD0combi = iGen;
              }
            }

            if (fabs(GenRecoD0combidR - 200.) > 1e-10) h_candGenRecoD0combidR->Fill(GenRecoD0combidR,weight);
            if (GenRecoD0combidR > 0.5) 
              isD0combiUnmatched = true;
            else {
              TLorentzVector p_genK, p_genPi;
              if (abs(genPartD0[iGenD0combi]->daughter(0)->pdgId()) == 321) {
                p_genK.SetPtEtaPhiM(genPartD0[iGenD0combi]->daughter(0)->pt(), genPartD0[iGenD0combi]->daughter(0)->eta(), genPartD0[iGenD0combi]->daughter(0)->phi(), genPartD0[iGenD0combi]->daughter(0)->mass());
                p_genPi.SetPtEtaPhiM(genPartD0[iGenD0combi]->daughter(1)->pt(), genPartD0[iGenD0combi]->daughter(1)->eta(), genPartD0[iGenD0combi]->daughter(1)->phi(), genPartD0[iGenD0combi]->daughter(1)->mass());
              }
              else {
                p_genK.SetPtEtaPhiM(genPartD0[iGenD0combi]->daughter(1)->pt(), genPartD0[iGenD0combi]->daughter(1)->eta(), genPartD0[iGenD0combi]->daughter(1)->phi(), genPartD0[iGenD0combi]->daughter(1)->mass());
                p_genPi.SetPtEtaPhiM(genPartD0[iGenD0combi]->daughter(0)->pt(), genPartD0[iGenD0combi]->daughter(0)->eta(), genPartD0[iGenD0combi]->daughter(0)->phi(), genPartD0[iGenD0combi]->daughter(0)->mass());
              }
              double GenKRecoTr1combidR = p_track1_D0combi.DeltaR(p_genK);
              double GenPiRecoTr1combidR = p_track1_D0combi.DeltaR(p_genPi);
              double GenKRecoTr2combidR = p_track2_D0combi.DeltaR(p_genK);
              double GenPiRecoTr2combidR = p_track2_D0combi.DeltaR(p_genPi);
              if (GenKRecoTr1combidR < GenPiRecoTr1combidR) GenRecoTr1combidR = GenKRecoTr1combidR;
              else GenRecoTr1combidR = GenPiRecoTr1combidR;
              if (GenPiRecoTr2combidR < GenKRecoTr2combidR) GenRecoTr2combidR = GenPiRecoTr2combidR;
              else GenRecoTr2combidR = GenKRecoTr2combidR;
              h_candGenRecoTrcombidR->Fill(GenRecoTr1combidR,weight);
              h_candGenRecoTrcombidR->Fill(GenRecoTr2combidR,weight);
              if (GenRecoTr1combidR > 0.005 || GenRecoTr2combidR > 0.005) 
                isD0combiUnmatched = true;
              else {
                if (GenKRecoTr1combidR < GenPiRecoTr1combidR && GenPiRecoTr2combidR < GenKRecoTr2combidR) isD0combiPaired = true;
                else isD0combiSwitched = true;
              }
            }

            if (fabs(p_D0combi.M() - gMassD0) < fabs(p_D0optcombi.M() - gMassD0)) {
              p_D0optcombi.SetPtEtaPhiM(p_D0combi.Pt(), p_D0combi.Eta(), p_D0combi.Phi(), p_D0combi.M());
              tk2charge = id_trCand_D0combi[tk2];
              isD0optcombiPaired = isD0combiPaired;
              isD0optcombiSwitched = isD0combiSwitched;
              isD0optcombiUnmatched = isD0combiUnmatched;
            }

            h_D0combiCand_pT->Fill(p_D0combi.Pt(),weight);

            // cut on pT
            if (p_D0combi.Pt() < 15.) continue;

            h_dRTr1Tr2combi->Fill(p_track1_D0combi.DeltaR(p_track2_D0combi),weight);
            h_dEtaTr1Tr2combi->Fill(fabs(p_track1_D0combi.Eta()-p_track2_D0combi.Eta()),weight);
            h_dPhiTr1Tr2combi->Fill(p_track1_D0combi.DeltaPhi(p_track2_D0combi),weight);
            h_pTrCombi->Fill(p_track1_D0combi.P(),weight);
            h_pTrCombi->Fill(p_track2_D0combi.P(),weight);
            h_pTTrCombi->Fill(p_track1_D0combi.Pt(),weight);
            h_pTTrCombi->Fill(p_track2_D0combi.Pt(),weight);
            h_etaTrCombi->Fill(p_track1_D0combi.Eta(),weight);
            h_etaTrCombi->Fill(p_track2_D0combi.Eta(),weight);
            h_phiTrCombi->Fill(p_track1_D0combi.Phi(),weight);
            h_phiTrCombi->Fill(p_track2_D0combi.Phi(),weight);
            h_pTr1pTr2combi->Fill(p_track1_D0combi.P()/p_track2_D0combi.P(),weight);
            h_pTTr1pTTr2combi->Fill(p_track1_D0combi.Pt()/p_track2_D0combi.Pt(),weight);
            h_D0combi_Mass->Fill(p_D0combi.M(),weight);
            h_D0combi_dRJet->Fill(p_D0combi.DeltaR(p_Jet),weight);
            h_D0combi_p->Fill(p_D0combi.P(),weight);
            h_D0combi_pT->Fill(p_D0combi.Pt(),weight);
            h_D0combi_eta->Fill(p_D0combi.Eta(),weight);
            h_D0combi_phi->Fill(p_D0combi.Phi(),weight);
            if (isD0combiPaired) {
              h_dRTr1Tr2combi_paired->Fill(p_track1_D0combi.DeltaR(p_track2_D0combi),weight);
              h_dEtaTr1Tr2combi_paired->Fill(fabs(p_track1_D0combi.Eta()-p_track2_D0combi.Eta()),weight);
              h_dPhiTr1Tr2combi_paired->Fill(p_track1_D0combi.DeltaPhi(p_track2_D0combi),weight);
              h_pTrCombi_paired->Fill(p_track1_D0combi.P(),weight);
              h_pTrCombi_paired->Fill(p_track2_D0combi.P(),weight);
              h_pTTrCombi_paired->Fill(p_track1_D0combi.Pt(),weight);
              h_pTTrCombi_paired->Fill(p_track2_D0combi.Pt(),weight);
              h_etaTrCombi_paired->Fill(p_track1_D0combi.Eta(),weight);
              h_etaTrCombi_paired->Fill(p_track2_D0combi.Eta(),weight);
              h_phiTrCombi_paired->Fill(p_track1_D0combi.Phi(),weight);
              h_phiTrCombi_paired->Fill(p_track2_D0combi.Phi(),weight);
              h_pTr1pTr2combi_paired->Fill(p_track1_D0combi.P()/p_track2_D0combi.P(),weight);
              h_pTTr1pTTr2combi_paired->Fill(p_track1_D0combi.Pt()/p_track2_D0combi.Pt(),weight);
              h_D0combi_paired_Mass->Fill(p_D0combi.M(),weight);
              h_D0combi_paired_dRJet->Fill(p_D0combi.DeltaR(p_Jet),weight);
              h_D0combi_paired_p->Fill(p_D0combi.P(),weight);
              h_D0combi_paired_pT->Fill(p_D0combi.Pt(),weight);
              h_D0combi_paired_eta->Fill(p_D0combi.Eta(),weight);
              h_D0combi_paired_phi->Fill(p_D0combi.Phi(),weight);
            }
            if (isD0combiSwitched) {
              h_dRTr1Tr2combi_switched->Fill(p_track1_D0combi.DeltaR(p_track2_D0combi),weight);
              h_dEtaTr1Tr2combi_switched->Fill(fabs(p_track1_D0combi.Eta()-p_track2_D0combi.Eta()),weight);
              h_dPhiTr1Tr2combi_switched->Fill(p_track1_D0combi.DeltaPhi(p_track2_D0combi),weight);
              h_pTrCombi_switched->Fill(p_track1_D0combi.P(),weight);
              h_pTrCombi_switched->Fill(p_track2_D0combi.P(),weight);
              h_pTTrCombi_switched->Fill(p_track1_D0combi.Pt(),weight);
              h_pTTrCombi_switched->Fill(p_track2_D0combi.Pt(),weight);
              h_etaTrCombi_switched->Fill(p_track1_D0combi.Eta(),weight);
              h_etaTrCombi_switched->Fill(p_track2_D0combi.Eta(),weight);
              h_phiTrCombi_switched->Fill(p_track1_D0combi.Phi(),weight);
              h_phiTrCombi_switched->Fill(p_track2_D0combi.Phi(),weight);
              h_pTr1pTr2combi_switched->Fill(p_track1_D0combi.P()/p_track2_D0combi.P(),weight);
              h_pTTr1pTTr2combi_switched->Fill(p_track1_D0combi.Pt()/p_track2_D0combi.Pt(),weight);
              h_D0combi_switched_Mass->Fill(p_D0combi.M(),weight);
              h_D0combi_switched_dRJet->Fill(p_D0combi.DeltaR(p_Jet),weight);
              h_D0combi_switched_p->Fill(p_D0combi.P(),weight);
              h_D0combi_switched_pT->Fill(p_D0combi.Pt(),weight);
              h_D0combi_switched_eta->Fill(p_D0combi.Eta(),weight);
              h_D0combi_switched_phi->Fill(p_D0combi.Phi(),weight);
            }
            if (isD0combiUnmatched) {
              h_dRTr1Tr2combi_unmatched->Fill(p_track1_D0combi.DeltaR(p_track2_D0combi),weight);
              h_dEtaTr1Tr2combi_unmatched->Fill(fabs(p_track1_D0combi.Eta()-p_track2_D0combi.Eta()),weight);
              h_dPhiTr1Tr2combi_unmatched->Fill(p_track1_D0combi.DeltaPhi(p_track2_D0combi),weight);
              h_pTrCombi_unmatched->Fill(p_track1_D0combi.P(),weight);
              h_pTrCombi_unmatched->Fill(p_track2_D0combi.P(),weight);
              h_pTTrCombi_unmatched->Fill(p_track1_D0combi.Pt(),weight);
              h_pTTrCombi_unmatched->Fill(p_track2_D0combi.Pt(),weight);
              h_etaTrCombi_unmatched->Fill(p_track1_D0combi.Eta(),weight);
              h_etaTrCombi_unmatched->Fill(p_track2_D0combi.Eta(),weight);
              h_phiTrCombi_unmatched->Fill(p_track1_D0combi.Phi(),weight);
              h_phiTrCombi_unmatched->Fill(p_track2_D0combi.Phi(),weight);
              h_pTr1pTr2combi_unmatched->Fill(p_track1_D0combi.P()/p_track2_D0combi.P(),weight);
              h_pTTr1pTTr2combi_unmatched->Fill(p_track1_D0combi.Pt()/p_track2_D0combi.Pt(),weight);
              h_D0combi_unmatched_Mass->Fill(p_D0combi.M(),weight);
              h_D0combi_unmatched_dRJet->Fill(p_D0combi.DeltaR(p_Jet),weight);
              h_D0combi_unmatched_p->Fill(p_D0combi.P(),weight);
              h_D0combi_unmatched_pT->Fill(p_D0combi.Pt(),weight);
              h_D0combi_unmatched_eta->Fill(p_D0combi.Eta(),weight);
              h_D0combi_unmatched_phi->Fill(p_D0combi.Phi(),weight);
            }

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // associate D^0 to a PF muon
            //~~~~~~~~~~~~~~~~~~~~~~~~~~
            TLorentzVector p_Mu;
            int iMu = -1;
            double deltaRD0combiMu = 2000.;

            // find closest muon with opposite charged wrt the kaon
            for (unsigned int iMuCand = 0; iMuCand < myPFmu.size(); iMuCand++) {
              if (myPFmu[iMuCand]->pdgId()*id_trCand_D0combi[tk2] > 0) continue;
              TLorentzVector p_MuCand;
              p_MuCand.SetPtEtaPhiM(myPFmu[iMuCand]->pt(), myPFmu[iMuCand]->eta(), myPFmu[iMuCand]->phi(), gMassMu);
              double tmp_deltaRD0combiMu = p_D0combi.DeltaR(p_MuCand);
              if (tmp_deltaRD0combiMu < deltaRD0combiMu) {
                deltaRD0combiMu = tmp_deltaRD0combiMu;
                iMu = iMuCand;
              }
            }

            if (iMu < 0) continue;
            h_BCand_DeltaRD0combiMu->Fill(deltaRD0combiMu,weight);

            // keep going if closest muon is close enough
            if (deltaRD0combiMu > 0.4) continue;
            h_B_DeltaRD0combiMu->Fill(deltaRD0combiMu,weight);

            h_B_D0combiMass->Fill(p_D0combi.M(),weight);

            p_Mu.SetPtEtaPhiM(myPFmu[iMu]->pt(), myPFmu[iMu]->eta(), myPFmu[iMu]->phi(), gMassMu);
            h_B_pD0combipMu->Fill(p_D0combi.P() / myPFmu[iMu]->p(),weight);
            h_B_pTD0combipTMu->Fill(p_D0combi.Pt() / myPFmu[iMu]->pt(),weight);

            TLorentzVector p_Bcombi = p_Mu + p_D0combi;
            h_B_combiMass->Fill(p_Bcombi.M(),weight);
            h_B_combip->Fill(p_Bcombi.P(),weight);
            h_B_combipT->Fill(p_Bcombi.Pt(),weight);
            h_B_combieta->Fill(p_Bcombi.Eta(),weight);
            h_B_combiphi->Fill(p_Bcombi.Phi(),weight);
            if (isD0combiPaired) {
              h_B_paired_DeltaRD0combiMu->Fill(deltaRD0combiMu,weight);
              h_B_paired_D0combiMass->Fill(p_D0combi.M(),weight);
              h_B_paired_pD0combipMu->Fill(p_D0combi.P() / myPFmu[iMu]->p(),weight);
              h_B_paired_pTD0combipTMu->Fill(p_D0combi.Pt() / myPFmu[iMu]->pt(),weight);
              h_B_paired_combiMass->Fill(p_Bcombi.M(),weight);
              h_B_paired_combip->Fill(p_Bcombi.P(),weight);
              h_B_paired_combipT->Fill(p_Bcombi.Pt(),weight);
              h_B_paired_combieta->Fill(p_Bcombi.Eta(),weight);
              h_B_paired_combiphi->Fill(p_Bcombi.Phi(),weight);
            }
            if (isD0combiSwitched) {
              h_B_switched_DeltaRD0combiMu->Fill(deltaRD0combiMu,weight);
              h_B_switched_D0combiMass->Fill(p_D0combi.M(),weight);
              h_B_switched_pD0combipMu->Fill(p_D0combi.P() / myPFmu[iMu]->p(),weight);
              h_B_switched_pTD0combipTMu->Fill(p_D0combi.Pt() / myPFmu[iMu]->pt(),weight);
              h_B_switched_combiMass->Fill(p_Bcombi.M(),weight);
              h_B_switched_combip->Fill(p_Bcombi.P(),weight);
              h_B_switched_combipT->Fill(p_Bcombi.Pt(),weight);
              h_B_switched_combieta->Fill(p_Bcombi.Eta(),weight);
              h_B_switched_combiphi->Fill(p_Bcombi.Phi(),weight);
            }
            if (isD0combiUnmatched) {
              h_B_unmatched_DeltaRD0combiMu->Fill(deltaRD0combiMu,weight);
              h_B_unmatched_D0combiMass->Fill(p_D0combi.M(),weight);
              h_B_unmatched_pD0combipMu->Fill(p_D0combi.P() / myPFmu[iMu]->p(),weight);
              h_B_unmatched_pTD0combipTMu->Fill(p_D0combi.Pt() / myPFmu[iMu]->pt(),weight);
              h_B_unmatched_combiMass->Fill(p_Bcombi.M(),weight);
              h_B_unmatched_combip->Fill(p_Bcombi.P(),weight);
              h_B_unmatched_combipT->Fill(p_Bcombi.Pt(),weight);
              h_B_unmatched_combieta->Fill(p_Bcombi.Eta(),weight);
              h_B_unmatched_combiphi->Fill(p_Bcombi.Phi(),weight);
            }
          }
        }

        //=================================
        // optimized invariant combination
        //=================================

        if (fabs(p_D0optcombi.Pt()) > 1e-10) {

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~
          // reconstruct D^0 to K Pi
          //~~~~~~~~~~~~~~~~~~~~~~~~~~

          h_D0optcombiCand_pT->Fill(p_D0optcombi.Pt(),weight);

          // cut on pT
          if (p_D0optcombi.Pt() < 15.) continue;

          h_D0optcombi_Mass->Fill(p_D0optcombi.M(),weight);
          h_D0optcombi_dRJet->Fill(p_D0optcombi.DeltaR(p_Jet),weight);
          h_D0optcombi_p->Fill(p_D0optcombi.P(),weight);
          h_D0optcombi_pT->Fill(p_D0optcombi.Pt(),weight);
          h_D0optcombi_eta->Fill(p_D0optcombi.Eta(),weight);
          h_D0optcombi_phi->Fill(p_D0optcombi.Phi(),weight);
          if (isD0optcombiPaired) {
            h_D0optcombi_paired_Mass->Fill(p_D0optcombi.M(),weight);
            h_D0optcombi_paired_dRJet->Fill(p_D0optcombi.DeltaR(p_Jet),weight);
            h_D0optcombi_paired_p->Fill(p_D0optcombi.P(),weight);
            h_D0optcombi_paired_pT->Fill(p_D0optcombi.Pt(),weight);
            h_D0optcombi_paired_eta->Fill(p_D0optcombi.Eta(),weight);
            h_D0optcombi_paired_phi->Fill(p_D0optcombi.Phi(),weight);
          }
          if (isD0optcombiSwitched) {
            h_D0optcombi_switched_Mass->Fill(p_D0optcombi.M(),weight);
            h_D0optcombi_switched_dRJet->Fill(p_D0optcombi.DeltaR(p_Jet),weight);
            h_D0optcombi_switched_p->Fill(p_D0optcombi.P(),weight);
            h_D0optcombi_switched_pT->Fill(p_D0optcombi.Pt(),weight);
            h_D0optcombi_switched_eta->Fill(p_D0optcombi.Eta(),weight);
            h_D0optcombi_switched_phi->Fill(p_D0optcombi.Phi(),weight);
          }
          if (isD0optcombiUnmatched) {
            h_D0optcombi_unmatched_Mass->Fill(p_D0optcombi.M(),weight);
            h_D0optcombi_unmatched_dRJet->Fill(p_D0optcombi.DeltaR(p_Jet),weight);
            h_D0optcombi_unmatched_p->Fill(p_D0optcombi.P(),weight);
            h_D0optcombi_unmatched_pT->Fill(p_D0optcombi.Pt(),weight);
            h_D0optcombi_unmatched_eta->Fill(p_D0optcombi.Eta(),weight);
            h_D0optcombi_unmatched_phi->Fill(p_D0optcombi.Phi(),weight);
          }

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~
          // associate D^0 to a PF muon
          //~~~~~~~~~~~~~~~~~~~~~~~~~~
          TLorentzVector p_Mu;
          int iMu = -1;
          double deltaRD0optcombiMu = 2000.;

          // find closest muon with opposite charged wrt the kaon
          for (unsigned int iMuCand = 0; iMuCand < myPFmu.size(); iMuCand++) {
            if (myPFmu[iMuCand]->pdgId()*tk2charge > 0) continue;  
            TLorentzVector p_MuCand;
            p_MuCand.SetPtEtaPhiM(myPFmu[iMuCand]->pt(), myPFmu[iMuCand]->eta(), myPFmu[iMuCand]->phi(), gMassMu);
            double tmp_deltaRD0optcombiMu = p_D0optcombi.DeltaR(p_MuCand);
            if (tmp_deltaRD0optcombiMu < deltaRD0optcombiMu) {
              deltaRD0optcombiMu = tmp_deltaRD0optcombiMu;
              iMu = iMuCand;
            }
          }

          if (iMu < 0) continue;
          h_BCand_DeltaRD0optcombiMu->Fill(deltaRD0optcombiMu,weight);

          // keep going if closest muon is close enough
          if (deltaRD0optcombiMu > 0.4) continue;
          h_B_DeltaRD0optcombiMu->Fill(deltaRD0optcombiMu,weight);

          h_B_D0optcombiMass->Fill(p_D0optcombi.M(),weight);
          
          p_Mu.SetPtEtaPhiM(myPFmu[iMu]->pt(), myPFmu[iMu]->eta(), myPFmu[iMu]->phi(), gMassMu);
          h_B_pD0optcombipMu->Fill(p_D0optcombi.P() / myPFmu[iMu]->p(),weight);
          h_B_pTD0optcombipTMu->Fill(p_D0optcombi.Pt() / myPFmu[iMu]->pt(),weight);

          TLorentzVector p_Boptcombi = p_Mu + p_D0optcombi;
          h_B_optcombiMass->Fill(p_Boptcombi.M(),weight);
          h_B_optcombip->Fill(p_Boptcombi.P(),weight);
          h_B_optcombipT->Fill(p_Boptcombi.Pt(),weight);
          h_B_optcombieta->Fill(p_Boptcombi.Eta(),weight);
          h_B_optcombiphi->Fill(p_Boptcombi.Phi(),weight);
          if (isD0optcombiPaired) {
            h_B_paired_DeltaRD0optcombiMu->Fill(deltaRD0optcombiMu,weight);
            h_B_paired_D0optcombiMass->Fill(p_D0optcombi.M(),weight);
            h_B_paired_pD0optcombipMu->Fill(p_D0optcombi.P() / myPFmu[iMu]->p(),weight);
            h_B_paired_pTD0optcombipTMu->Fill(p_D0optcombi.Pt() / myPFmu[iMu]->pt(),weight);
            h_B_paired_optcombiMass->Fill(p_Boptcombi.M(),weight);
            h_B_paired_optcombip->Fill(p_Boptcombi.P(),weight);
            h_B_paired_optcombipT->Fill(p_Boptcombi.Pt(),weight);
            h_B_paired_optcombieta->Fill(p_Boptcombi.Eta(),weight);
            h_B_paired_optcombiphi->Fill(p_Boptcombi.Phi(),weight);
          }
          if (isD0optcombiSwitched) {
            h_B_switched_DeltaRD0optcombiMu->Fill(deltaRD0optcombiMu,weight);
            h_B_switched_D0optcombiMass->Fill(p_D0optcombi.M(),weight);
            h_B_switched_pD0optcombipMu->Fill(p_D0optcombi.P() / myPFmu[iMu]->p(),weight);
            h_B_switched_pTD0optcombipTMu->Fill(p_D0optcombi.Pt() / myPFmu[iMu]->pt(),weight);
            h_B_switched_optcombiMass->Fill(p_Boptcombi.M(),weight);
            h_B_switched_optcombip->Fill(p_Boptcombi.P(),weight);
            h_B_switched_optcombipT->Fill(p_Boptcombi.Pt(),weight);
            h_B_switched_optcombieta->Fill(p_Boptcombi.Eta(),weight);
            h_B_switched_optcombiphi->Fill(p_Boptcombi.Phi(),weight);
          }
          if (isD0optcombiUnmatched) {
            h_B_unmatched_DeltaRD0optcombiMu->Fill(deltaRD0optcombiMu,weight);
            h_B_unmatched_D0optcombiMass->Fill(p_D0optcombi.M(),weight);
            h_B_unmatched_pD0optcombipMu->Fill(p_D0optcombi.P() / myPFmu[iMu]->p(),weight);
            h_B_unmatched_pTD0optcombipTMu->Fill(p_D0optcombi.Pt() / myPFmu[iMu]->pt(),weight);
            h_B_unmatched_optcombiMass->Fill(p_Boptcombi.M(),weight);
            h_B_unmatched_optcombip->Fill(p_Boptcombi.P(),weight);
            h_B_unmatched_optcombipT->Fill(p_Boptcombi.Pt(),weight);
            h_B_unmatched_optcombieta->Fill(p_Boptcombi.Eta(),weight);
            h_B_unmatched_optcombiphi->Fill(p_Boptcombi.Phi(),weight);
          }
        }
      }
      //iSelJet++;
    } // jet loop
  }
}


// ------------ method called once each job just before starting event loop  ------------
  void 
KalmanAnalyzer::beginJob()
{

  //  std::cout << "Creating histos..." << std::endl;

  edm::Service<TFileService> fs;

  // evts properties
  h_nJets = fs->make<TH1D>("h_nJets","h_nJets", 20, 0., 20.);
  h_CSV = fs->make<TH1D>("h_CSV","h_CSV", 100, 0., 1.);
  h_DeltaR_trCand_Mu = fs->make<TH1D>("h_DeltaR_trCand_Mu","h_DeltaR_trCand_Mu",1000,0.,0.5);
  h_DeltaR_trCand_El = fs->make<TH1D>("h_DeltaR_trCand_El","h_DeltaR_trCand_El",1000,0.,0.5);

  h_genD0_n    = fs->make<TH1D>("h_genD0_n","h_genD0_n",4,0.,4.);
  h_genD0_p    = fs->make<TH1D>("h_genD0_p","h_genD0_p",1000,0.,500.);
  h_genD0_pT   = fs->make<TH1D>("h_genD0_pT","h_genD0_pT",1000,0.,500.);
  h_genD0_eta  = fs->make<TH1D>("h_genD0_eta","h_genD0_eta",60,-3.,3.);
  h_genD0_phi  = fs->make<TH1D>("h_genD0_phi","h_genD0_phi",64,-3.2,3.2);

  // simple Kalman vertex fitter
  h_B_cuts = fs->make<TH1D>("h_B_cuts","h_B_cuts",30,0.,30.);
  h_B_cuts->SetOption("bar");
  h_B_cuts->SetBarWidth(0.75);
  h_B_cuts->SetBarOffset(0.125);

  h_candGenRecoD0dR   = fs->make<TH1D>("h_candGenRecoD0dR","h_candGenRecoD0dR",200,0.,1.);
  h_candGenRecoTrdR   = fs->make<TH1D>("h_candGenRecoTrdR","h_candGenRecoTrdR",200,0.,1.);

  h_D0Cand_Chi2NDOF       = fs->make<TH1D>("h_D0Cand_Chi2NDOF","h_D0Cand_Chi2NDOF",110,0.,11.);
  h_D0Cand_MassChi2Inf1   = fs->make<TH1D>("h_D0Cand_MassChi2Inf1","h_D0Cand_MassChi2Inf1",1000,0.,10.);
  h_D0Cand_MassChi2Inf1p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf1p5","h_D0Cand_MassChi2Inf1p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf2   = fs->make<TH1D>("h_D0Cand_MassChi2Inf2","h_D0Cand_MassChi2Inf2",1000,0.,10.);
  h_D0Cand_MassChi2Inf2p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf2p5","h_D0Cand_MassChi2Inf2p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf3   = fs->make<TH1D>("h_D0Cand_MassChi2Inf3","h_D0Cand_MassChi2Inf3",1000,0.,10.);
  h_D0Cand_MassChi2Inf3p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf3p5","h_D0Cand_MassChi2Inf3p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf4   = fs->make<TH1D>("h_D0Cand_MassChi2Inf4","h_D0Cand_MassChi2Inf4",1000,0.,10.);
  h_D0Cand_MassChi2Inf4p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf4p5","h_D0Cand_MassChi2Inf4p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf5   = fs->make<TH1D>("h_D0Cand_MassChi2Inf5","h_D0Cand_MassChi2Inf5",1000,0.,10.);
  h_D0Cand_MassChi2Inf5p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf5p5","h_D0Cand_MassChi2Inf5p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf6   = fs->make<TH1D>("h_D0Cand_MassChi2Inf6","h_D0Cand_MassChi2Inf6",1000,0.,10.);
  h_D0Cand_MassChi2Inf6p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf6p5","h_D0Cand_MassChi2Inf6p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf7   = fs->make<TH1D>("h_D0Cand_MassChi2Inf7","h_D0Cand_MassChi2Inf7",1000,0.,10.);
  h_D0Cand_MassChi2Inf7p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf7p5","h_D0Cand_MassChi2Inf7p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf8   = fs->make<TH1D>("h_D0Cand_MassChi2Inf8","h_D0Cand_MassChi2Inf8",1000,0.,10.);

  h_D0Cand_L           = fs->make<TH1D>("h_D0Cand_L","h_D0Cand_L",1000,0.,1.);
  h_D0Cand_SigmaL      = fs->make<TH1D>("h_D0Cand_SigmaL","h_D0Cand_SigmaL",5000,0.,0.005);
  h_D0Cand_LOverSigmaL = fs->make<TH1D>("h_D0Cand_LOverSigmaL","h_D0Cand_LOverSigmaL",21000,0.,7000.);

  h_D0Cand_pT      = fs->make<TH1D>("h_D0Cand_pT","h_D0Cand_pT",1000,0.,500.);

  h_D0_Mass    = fs->make<TH1D>("h_D0_Mass","h_D0_Mass",1000,0.,10.);
  h_D0_L       = fs->make<TH1D>("h_D0_L","h_D0_L",1000,0.,1.);
  h_D0_SigmaL  = fs->make<TH1D>("h_D0_SigmaL","h_D0_SigmaL",5000,0.,0.005);
  h_D0_dRJet   = fs->make<TH1D>("h_D0_dRJet","h_D0_dRJet",200,0.,1.);
  h_D0_p       = fs->make<TH1D>("h_D0_p","h_D0_p",1000,0.,500.);
  h_D0_pT      = fs->make<TH1D>("h_D0_pT","h_D0_pT",1000,0.,500.);
  h_D0_eta     = fs->make<TH1D>("h_D0_eta","h_D0_eta",60,-3.,3.);
  h_D0_phi     = fs->make<TH1D>("h_D0_phi","h_D0_phi",64,-3.2,3.2);
  h_D0_paired_Mass    = fs->make<TH1D>("h_D0_paired_Mass","h_D0_paired_Mass",1000,0.,10.);
  h_D0_paired_L       = fs->make<TH1D>("h_D0_paired_L","h_D0_paired_L",1000,0.,1.);
  h_D0_paired_SigmaL  = fs->make<TH1D>("h_D0_paired_SigmaL","h_D0_paired_SigmaL",5000,0.,0.005);
  h_D0_paired_dRJet   = fs->make<TH1D>("h_D0_paired_dRJet","h_D0_paired_dRJet",200,0.,1.);
  h_D0_paired_p       = fs->make<TH1D>("h_D0_paired_p","h_D0_paired_p",1000,0.,500.);
  h_D0_paired_pT      = fs->make<TH1D>("h_D0_paired_pT","h_D0_paired_pT",1000,0.,500.);
  h_D0_paired_eta     = fs->make<TH1D>("h_D0_paired_eta","h_D0_paired_eta",60,-3.,3.);
  h_D0_paired_phi     = fs->make<TH1D>("h_D0_paired_phi","h_D0_paired_phi",64,-3.2,3.2);
  h_D0_switched_Mass    = fs->make<TH1D>("h_D0_switched_Mass","h_D0_switched_Mass",1000,0.,10.);
  h_D0_switched_L       = fs->make<TH1D>("h_D0_switched_L","h_D0_switched_L",1000,0.,1.);
  h_D0_switched_SigmaL  = fs->make<TH1D>("h_D0_switched_SigmaL","h_D0_switched_SigmaL",5000,0.,0.005);
  h_D0_switched_dRJet   = fs->make<TH1D>("h_D0_switched_dRJet","h_D0_switched_dRJet",200,0.,1.);
  h_D0_switched_p       = fs->make<TH1D>("h_D0_switched_p","h_D0_switched_p",1000,0.,500.);
  h_D0_switched_pT      = fs->make<TH1D>("h_D0_switched_pT","h_D0_switched_pT",1000,0.,500.);
  h_D0_switched_eta     = fs->make<TH1D>("h_D0_switched_eta","h_D0_switched_eta",60,-3.,3.);
  h_D0_switched_phi     = fs->make<TH1D>("h_D0_switched_phi","h_D0_switched_phi",64,-3.2,3.2);
  h_D0_unmatched_Mass    = fs->make<TH1D>("h_D0_unmatched_Mass","h_D0_unmatched_Mass",1000,0.,10.);
  h_D0_unmatched_L       = fs->make<TH1D>("h_D0_unmatched_L","h_D0_unmatched_L",1000,0.,1.);
  h_D0_unmatched_SigmaL  = fs->make<TH1D>("h_D0_unmatched_SigmaL","h_D0_unmatched_SigmaL",5000,0.,0.005);
  h_D0_unmatched_dRJet   = fs->make<TH1D>("h_D0_unmatched_dRJet","h_D0_unmatched_dRJet",200,0.,1.);
  h_D0_unmatched_p       = fs->make<TH1D>("h_D0_unmatched_p","h_D0_unmatched_p",1000,0.,500.);
  h_D0_unmatched_pT      = fs->make<TH1D>("h_D0_unmatched_pT","h_D0_unmatched_pT",1000,0.,500.);
  h_D0_unmatched_eta     = fs->make<TH1D>("h_D0_unmatched_eta","h_D0_unmatched_eta",60,-3.,3.);
  h_D0_unmatched_phi     = fs->make<TH1D>("h_D0_unmatched_phi","h_D0_unmatched_phi",64,-3.2,3.2);

  h_BCand_DeltaRD0Mu = fs->make<TH1D>("h_BCand_DeltaRD0Mu","h_BCand_DeltaRD0Mu",100,0.,5.);
  h_B_D0Mass     = fs->make<TH1D>("h_B_D0Mass","h_B_D0Mass",1000,0.,10.);
  h_B_DeltaRD0Mu = fs->make<TH1D>("h_B_DeltaRD0Mu","h_B_DeltaRD0Mu",80,0.,0.4);
  h_B_pD0pMu     = fs->make<TH1D>("h_B_pD0pMu","h_B_pD0pMu",500,0.,10.);
  h_B_pTD0pTMu   = fs->make<TH1D>("h_B_pTD0pTMu","h_B_pTD0pTMu",500,0.,10.);
  h_B_Mass       = fs->make<TH1D>("h_B_Mass","h_B_Mass",1000,0.,10.);
  h_B_p          = fs->make<TH1D>("h_B_p","h_B_p",1000,0.,500.);
  h_B_pT         = fs->make<TH1D>("h_B_pT","h_B_pT",1000,0.,500.);
  h_B_eta        = fs->make<TH1D>("h_B_eta","h_B_eta",60,-3.,3.);
  h_B_phi        = fs->make<TH1D>("h_B_phi","h_B_phi",64,-3.2,3.2);
  h_B_paired_D0Mass     = fs->make<TH1D>("h_B_paired_D0Mass","h_B_paired_D0Mass",1000,0.,10.);
  h_B_paired_DeltaRD0Mu = fs->make<TH1D>("h_B_paired_DeltaRD0Mu","h_B_paired_DeltaRD0Mu",80,0.,0.4);
  h_B_paired_pD0pMu     = fs->make<TH1D>("h_B_paired_pD0pMu","h_B_paired_pD0pMu",500,0.,10.);
  h_B_paired_pTD0pTMu   = fs->make<TH1D>("h_B_paired_pTD0pTMu","h_B_paired_pTD0pTMu",500,0.,10.);
  h_B_paired_Mass       = fs->make<TH1D>("h_B_paired_Mass","h_B_paired_Mass",1000,0.,10.);
  h_B_paired_p          = fs->make<TH1D>("h_B_paired_p","h_B_paired_p",1000,0.,500.);
  h_B_paired_pT         = fs->make<TH1D>("h_B_paired_pT","h_B_paired_pT",1000,0.,500.);
  h_B_paired_eta        = fs->make<TH1D>("h_B_paired_eta","h_B_paired_eta",60,-3.,3.);
  h_B_paired_phi        = fs->make<TH1D>("h_B_paired_phi","h_B_paired_phi",64,-3.2,3.2);
  h_B_switched_D0Mass     = fs->make<TH1D>("h_B_switched_D0Mass","h_B_switched_D0Mass",1000,0.,10.);
  h_B_switched_DeltaRD0Mu = fs->make<TH1D>("h_B_switched_DeltaRD0Mu","h_B_switched_DeltaRD0Mu",80,0.,0.4);
  h_B_switched_pD0pMu     = fs->make<TH1D>("h_B_switched_pD0pMu","h_B_switched_pD0pMu",500,0.,10.);
  h_B_switched_pTD0pTMu   = fs->make<TH1D>("h_B_switched_pTD0pTMu","h_B_switched_pTD0pTMu",500,0.,10.);
  h_B_switched_Mass       = fs->make<TH1D>("h_B_switched_Mass","h_B_switched_Mass",1000,0.,10.);
  h_B_switched_p          = fs->make<TH1D>("h_B_switched_p","h_B_switched_p",1000,0.,500.);
  h_B_switched_pT         = fs->make<TH1D>("h_B_switched_pT","h_B_switched_pT",1000,0.,500.);
  h_B_switched_eta        = fs->make<TH1D>("h_B_switched_eta","h_B_switched_eta",60,-3.,3.);
  h_B_switched_phi        = fs->make<TH1D>("h_B_switched_phi","h_B_switched_phi",64,-3.2,3.2);
  h_B_unmatched_D0Mass     = fs->make<TH1D>("h_B_unmatched_D0Mass","h_B_unmatched_D0Mass",1000,0.,10.);
  h_B_unmatched_DeltaRD0Mu = fs->make<TH1D>("h_B_unmatched_DeltaRD0Mu","h_B_unmatched_DeltaRD0Mu",80,0.,0.4);
  h_B_unmatched_pD0pMu     = fs->make<TH1D>("h_B_unmatched_pD0pMu","h_B_unmatched_pD0pMu",500,0.,10.);
  h_B_unmatched_pTD0pTMu   = fs->make<TH1D>("h_B_unmatched_pTD0pTMu","h_B_unmatched_pTD0pTMu",500,0.,10.);
  h_B_unmatched_Mass       = fs->make<TH1D>("h_B_unmatched_Mass","h_B_unmatched_Mass",1000,0.,10.);
  h_B_unmatched_p          = fs->make<TH1D>("h_B_unmatched_p","h_B_unmatched_p",1000,0.,500.);
  h_B_unmatched_pT         = fs->make<TH1D>("h_B_unmatched_pT","h_B_unmatched_pT",1000,0.,500.);
  h_B_unmatched_eta        = fs->make<TH1D>("h_B_unmatched_eta","h_B_unmatched_eta",60,-3.,3.);
  h_B_unmatched_phi        = fs->make<TH1D>("h_B_unmatched_phi","h_B_unmatched_phi",64,-3.2,3.2);

  // constrained Kalman Vertex fitter
  h_D0consCand_L           = fs->make<TH1D>("h_D0consCand_L","h_D0consCand_L",1000,0.,1.);
  h_D0consCand_SigmaL      = fs->make<TH1D>("h_D0consCand_SigmaL","h_D0consCand_SigmaL",5000,0.,0.005);
  h_D0consCand_LOverSigmaL = fs->make<TH1D>("h_D0consCand_LOverSigmaL","h_D0consCand_LOverSigmaL",21000,0.,7000.);

  h_D0consCand_pT      = fs->make<TH1D>("h_D0consCand_pT","h_D0consCand_pT",1000,0.,500.);

  h_D0cons_Chi2NDOF = fs->make<TH1D>("h_D0cons_Chi2NDOF","h_D0cons_Chi2NDOF",110,0.,11.);
  h_D0cons_L       = fs->make<TH1D>("h_D0cons_L","h_D0cons_L",1000,0.,1.);
  h_D0cons_SigmaL  = fs->make<TH1D>("h_D0cons_SigmaL","h_D0cons_SigmaL",5000,0.,0.005);
  h_D0cons_Mass    = fs->make<TH1D>("h_D0cons_Mass","h_D0cons_Mass",1000,0.,10.);
  h_D0cons_dRJet   = fs->make<TH1D>("h_D0cons_dRJet","h_D0cons_dRJet",200,0.,1.);
  h_D0cons_p       = fs->make<TH1D>("h_D0cons_p","h_D0cons_p",1000,0.,500.);
  h_D0cons_pT      = fs->make<TH1D>("h_D0cons_pT","h_D0cons_pT",1000,0.,500.);
  h_D0cons_eta     = fs->make<TH1D>("h_D0cons_eta","h_D0cons_eta",60,-3.,3.);
  h_D0cons_phi     = fs->make<TH1D>("h_D0cons_phi","h_D0cons_phi",64,-3.2,3.2);
  h_D0cons_paired_Chi2NDOF = fs->make<TH1D>("h_D0cons_paired_Chi2NDOF","h_D0cons_paired_Chi2NDOF",110,0.,11.);
  h_D0cons_paired_L       = fs->make<TH1D>("h_D0cons_paired_L","h_D0cons_paired_L",1000,0.,1.);
  h_D0cons_paired_SigmaL  = fs->make<TH1D>("h_D0cons_paired_SigmaL","h_D0cons_paired_SigmaL",5000,0.,0.005);
  h_D0cons_paired_Mass    = fs->make<TH1D>("h_D0cons_paired_Mass","h_D0cons_paired_Mass",1000,0.,10.);
  h_D0cons_paired_dRJet   = fs->make<TH1D>("h_D0cons_paired_dRJet","h_D0cons_paired_dRJet",200,0.,1.);
  h_D0cons_paired_p       = fs->make<TH1D>("h_D0cons_paired_p","h_D0cons_paired_p",1000,0.,500.);
  h_D0cons_paired_pT      = fs->make<TH1D>("h_D0cons_paired_pT","h_D0cons_paired_pT",1000,0.,500.);
  h_D0cons_paired_eta     = fs->make<TH1D>("h_D0cons_paired_eta","h_D0cons_paired_eta",60,-3.,3.);
  h_D0cons_paired_phi     = fs->make<TH1D>("h_D0cons_paired_phi","h_D0cons_paired_phi",64,-3.2,3.2);
  h_D0cons_switched_Chi2NDOF = fs->make<TH1D>("h_D0cons_switched_Chi2NDOF","h_D0cons_switched_Chi2NDOF",110,0.,11.);
  h_D0cons_switched_L       = fs->make<TH1D>("h_D0cons_switched_L","h_D0cons_switched_L",1000,0.,1.);
  h_D0cons_switched_SigmaL  = fs->make<TH1D>("h_D0cons_switched_SigmaL","h_D0cons_switched_SigmaL",5000,0.,0.005);
  h_D0cons_switched_Mass    = fs->make<TH1D>("h_D0cons_switched_Mass","h_D0cons_switched_Mass",1000,0.,10.);
  h_D0cons_switched_dRJet   = fs->make<TH1D>("h_D0cons_switched_dRJet","h_D0cons_switched_dRJet",200,0.,1.);
  h_D0cons_switched_p       = fs->make<TH1D>("h_D0cons_switched_p","h_D0cons_switched_p",1000,0.,500.);
  h_D0cons_switched_pT      = fs->make<TH1D>("h_D0cons_switched_pT","h_D0cons_switched_pT",1000,0.,500.);
  h_D0cons_switched_eta     = fs->make<TH1D>("h_D0cons_switched_eta","h_D0cons_switched_eta",60,-3.,3.);
  h_D0cons_switched_phi     = fs->make<TH1D>("h_D0cons_switched_phi","h_D0cons_switched_phi",64,-3.2,3.2);
  h_D0cons_unmatched_Chi2NDOF = fs->make<TH1D>("h_D0cons_unmatched_Chi2NDOF","h_D0cons_unmatched_Chi2NDOF",110,0.,11.);
  h_D0cons_unmatched_L       = fs->make<TH1D>("h_D0cons_unmatched_L","h_D0cons_unmatched_L",1000,0.,1.);
  h_D0cons_unmatched_SigmaL  = fs->make<TH1D>("h_D0cons_unmatched_SigmaL","h_D0cons_unmatched_SigmaL",5000,0.,0.005);
  h_D0cons_unmatched_Mass    = fs->make<TH1D>("h_D0cons_unmatched_Mass","h_D0cons_unmatched_Mass",1000,0.,10.);
  h_D0cons_unmatched_dRJet   = fs->make<TH1D>("h_D0cons_unmatched_dRJet","h_D0cons_unmatched_dRJet",200,0.,1.);
  h_D0cons_unmatched_p       = fs->make<TH1D>("h_D0cons_unmatched_p","h_D0cons_unmatched_p",1000,0.,500.);
  h_D0cons_unmatched_pT      = fs->make<TH1D>("h_D0cons_unmatched_pT","h_D0cons_unmatched_pT",1000,0.,500.);
  h_D0cons_unmatched_eta     = fs->make<TH1D>("h_D0cons_unmatched_eta","h_D0cons_unmatched_eta",60,-3.,3.);
  h_D0cons_unmatched_phi     = fs->make<TH1D>("h_D0cons_unmatched_phi","h_D0cons_unmatched_phi",64,-3.2,3.2);

  h_BCand_DeltaRD0consMu = fs->make<TH1D>("h_BCand_DeltaRD0consMu","h_BCand_DeltaRD0consMu",100,0.,5.);
  h_B_D0consMass     = fs->make<TH1D>("h_B_D0consMass","h_B_D0consMass",1000,0.,10.);
  h_B_DeltaRD0consMu = fs->make<TH1D>("h_B_DeltaRD0consMu","h_B_DeltaRD0consMu",80,0.,0.4);
  h_B_pD0conspMu     = fs->make<TH1D>("h_B_pD0conspMu","h_B_pD0conspMu",500,0.,10.);
  h_B_pTD0conspTMu   = fs->make<TH1D>("h_B_pTD0conspTMu","h_B_pTD0conspTMu",500,0.,10.);
  h_B_consMass       = fs->make<TH1D>("h_B_consMass","h_B_consMass",1000,0.,10.);
  h_B_consp          = fs->make<TH1D>("h_B_consp","h_B_consp",1000,0.,500.);
  h_B_conspT         = fs->make<TH1D>("h_B_conspT","h_B_conspT",1000,0.,500.);
  h_B_conseta        = fs->make<TH1D>("h_B_conseta","h_B_conseta",60,-3.,3.);
  h_B_consphi        = fs->make<TH1D>("h_B_consphi","h_B_consphi",64,-3.2,3.2);
  h_B_paired_D0consMass     = fs->make<TH1D>("h_B_paired_D0consMass","h_B_paired_D0consMass",1000,0.,10.);
  h_B_paired_DeltaRD0consMu = fs->make<TH1D>("h_B_paired_DeltaRD0consMu","h_B_paired_DeltaRD0consMu",80,0.,0.4);
  h_B_paired_pD0conspMu     = fs->make<TH1D>("h_B_paired_pD0conspMu","h_B_paired_pD0conspMu",500,0.,10.);
  h_B_paired_pTD0conspTMu   = fs->make<TH1D>("h_B_paired_pTD0conspTMu","h_B_paired_pTD0conspTMu",500,0.,10.);
  h_B_paired_consMass       = fs->make<TH1D>("h_B_paired_consMass","h_B_paired_consMass",1000,0.,10.);
  h_B_paired_consp          = fs->make<TH1D>("h_B_paired_consp","h_B_paired_consp",1000,0.,500.);
  h_B_paired_conspT         = fs->make<TH1D>("h_B_paired_conspT","h_B_paired_conspT",1000,0.,500.);
  h_B_paired_conseta        = fs->make<TH1D>("h_B_paired_conseta","h_B_paired_conseta",60,-3.,3.);
  h_B_paired_consphi        = fs->make<TH1D>("h_B_paired_consphi","h_B_paired_consphi",64,-3.2,3.2);
  h_B_switched_D0consMass     = fs->make<TH1D>("h_B_switched_D0consMass","h_B_switched_D0consMass",1000,0.,10.);
  h_B_switched_DeltaRD0consMu = fs->make<TH1D>("h_B_switched_DeltaRD0consMu","h_B_switched_DeltaRD0consMu",80,0.,0.4);
  h_B_switched_pD0conspMu     = fs->make<TH1D>("h_B_switched_pD0conspMu","h_B_switched_pD0conspMu",500,0.,10.);
  h_B_switched_pTD0conspTMu   = fs->make<TH1D>("h_B_switched_pTD0conspTMu","h_B_switched_pTD0conspTMu",500,0.,10.);
  h_B_switched_consMass       = fs->make<TH1D>("h_B_switched_consMass","h_B_switched_consMass",1000,0.,10.);
  h_B_switched_consp          = fs->make<TH1D>("h_B_switched_consp","h_B_switched_consp",1000,0.,500.);
  h_B_switched_conspT         = fs->make<TH1D>("h_B_switched_conspT","h_B_switched_conspT",1000,0.,500.);
  h_B_switched_conseta        = fs->make<TH1D>("h_B_switched_conseta","h_B_switched_conseta",60,-3.,3.);
  h_B_switched_consphi        = fs->make<TH1D>("h_B_switched_consphi","h_B_switched_consphi",64,-3.2,3.2);
  h_B_unmatched_D0consMass     = fs->make<TH1D>("h_B_unmatched_D0consMass","h_B_unmatched_D0consMass",1000,0.,10.);
  h_B_unmatched_DeltaRD0consMu = fs->make<TH1D>("h_B_unmatched_DeltaRD0consMu","h_B_unmatched_DeltaRD0consMu",80,0.,0.4);
  h_B_unmatched_pD0conspMu     = fs->make<TH1D>("h_B_unmatched_pD0conspMu","h_B_unmatched_pD0conspMu",500,0.,10.);
  h_B_unmatched_pTD0conspTMu   = fs->make<TH1D>("h_B_unmatched_pTD0conspTMu","h_B_unmatched_pTD0conspTMu",500,0.,10.);
  h_B_unmatched_consMass       = fs->make<TH1D>("h_B_unmatched_consMass","h_B_unmatched_consMass",1000,0.,10.);
  h_B_unmatched_consp          = fs->make<TH1D>("h_B_unmatched_consp","h_B_unmatched_consp",1000,0.,500.);
  h_B_unmatched_conspT         = fs->make<TH1D>("h_B_unmatched_conspT","h_B_unmatched_conspT",1000,0.,500.);
  h_B_unmatched_conseta        = fs->make<TH1D>("h_B_unmatched_conseta","h_B_unmatched_conseta",60,-3.,3.);
  h_B_unmatched_consphi        = fs->make<TH1D>("h_B_unmatched_consphi","h_B_unmatched_consphi",64,-3.2,3.2);

  // simple invariant combination
  h_candGenRecoD0combidR   = fs->make<TH1D>("h_candGenRecoD0combidR","h_candGenRecoD0combidR",200,0.,1.);
  h_candGenRecoTrcombidR   = fs->make<TH1D>("h_candGenRecoTrcombidR","h_candGenRecoTrcombidR",200,0.,1.);

  h_D0combiCand_pT = fs->make<TH1D>("h_D0combiCand_pT","h_D0combiCand_pT",1000,0.,500.);

  h_dRTr1Tr2combi = fs->make<TH1D>("h_dRTr1Tr2combi","h_dRTr1Tr2combi",200,0.,1.);
  h_dEtaTr1Tr2combi = fs->make<TH1D>("h_dEtaTr1Tr2combi","h_dEtaTr1Tr2combi",200,0.,1.);
  h_dPhiTr1Tr2combi = fs->make<TH1D>("h_dPhiTr1Tr2combi","h_dPhiTr1Tr2combi",200,0.,1.);
  h_pTrCombi      = fs->make<TH1D>("h_pTrCombi","h_pTrCombi",1000,0.,500.);
  h_pTTrCombi     = fs->make<TH1D>("h_pTTrCombi","h_pTTrCombi",1000,0.,500.);
  h_etaTrCombi    = fs->make<TH1D>("h_etaTrCombi","h_etaTrCombi",60,-3.,3.);
  h_phiTrCombi    = fs->make<TH1D>("h_phiTrCombi","h_phiTrCombi",64,-3.2,3.2);
  h_pTr1pTr2combi   = fs->make<TH1D>("h_pTr1pTr2combi","h_pTr1pTr2combi",500,0.,10.);
  h_pTTr1pTTr2combi = fs->make<TH1D>("h_pTTr1pTTr2combi","h_pTTr1pTTr2combi",500,0.,10.);
  h_dRTr1Tr2combi_paired = fs->make<TH1D>("h_dRTr1Tr2combi_paired","h_dRTr1Tr2combi_paired",200,0.,1.);
  h_dEtaTr1Tr2combi_paired = fs->make<TH1D>("h_dEtaTr1Tr2combi_paired","h_dEtaTr1Tr2combi_paired",200,0.,1.);
  h_dPhiTr1Tr2combi_paired = fs->make<TH1D>("h_dPhiTr1Tr2combi_paired","h_dPhiTr1Tr2combi_paired",200,0.,1.);
  h_pTrCombi_paired      = fs->make<TH1D>("h_pTrCombi_paired","h_pTrCombi_paired",1000,0.,500.);
  h_pTTrCombi_paired     = fs->make<TH1D>("h_pTTrCombi_paired","h_pTTrCombi_paired",1000,0.,500.);
  h_etaTrCombi_paired    = fs->make<TH1D>("h_etaTrCombi_paired","h_etaTrCombi_paired",60,-3.,3.);
  h_phiTrCombi_paired    = fs->make<TH1D>("h_phiTrCombi_paired","h_phiTrCombi_paired",64,-3.2,3.2);
  h_pTr1pTr2combi_paired   = fs->make<TH1D>("h_pTr1pTr2combi_paired","h_pTr1pTr2combi_paired",500,0.,10.);
  h_pTTr1pTTr2combi_paired = fs->make<TH1D>("h_pTTr1pTTr2combi_paired","h_pTTr1pTTr2combi_paired",500,0.,10.);
  h_dRTr1Tr2combi_switched = fs->make<TH1D>("h_dRTr1Tr2combi_switched","h_dRTr1Tr2combi_switched",200,0.,1.);
  h_dEtaTr1Tr2combi_switched = fs->make<TH1D>("h_dEtaTr1Tr2combi_switched","h_dEtaTr1Tr2combi_switched",200,0.,1.);
  h_dPhiTr1Tr2combi_switched = fs->make<TH1D>("h_dPhiTr1Tr2combi_switched","h_dPhiTr1Tr2combi_switched",200,0.,1.);
  h_pTrCombi_switched      = fs->make<TH1D>("h_pTrCombi_switched","h_pTrCombi_switched",1000,0.,500.);
  h_pTTrCombi_switched     = fs->make<TH1D>("h_pTTrCombi_switched","h_pTTrCombi_switched",1000,0.,500.);
  h_etaTrCombi_switched    = fs->make<TH1D>("h_etaTrCombi_switched","h_etaTrCombi_switched",60,-3.,3.);
  h_phiTrCombi_switched    = fs->make<TH1D>("h_phiTrCombi_switched","h_phiTrCombi_switched",64,-3.2,3.2);
  h_pTr1pTr2combi_switched   = fs->make<TH1D>("h_pTr1pTr2combi_switched","h_pTr1pTr2combi_switched",500,0.,10.);
  h_pTTr1pTTr2combi_switched = fs->make<TH1D>("h_pTTr1pTTr2combi_switched","h_pTTr1pTTr2combi_switched",500,0.,10.);
  h_dRTr1Tr2combi_unmatched = fs->make<TH1D>("h_dRTr1Tr2combi_unmatched","h_dRTr1Tr2combi_unmatched",200,0.,1.);
  h_dEtaTr1Tr2combi_unmatched = fs->make<TH1D>("h_dEtaTr1Tr2combi_unmatched","h_dEtaTr1Tr2combi_unmatched",200,0.,1.);
  h_dPhiTr1Tr2combi_unmatched = fs->make<TH1D>("h_dPhiTr1Tr2combi_unmatched","h_dPhiTr1Tr2combi_unmatched",200,0.,1.);
  h_pTrCombi_unmatched      = fs->make<TH1D>("h_pTrCombi_unmatched","h_pTrCombi_unmatched",1000,0.,500.);
  h_pTTrCombi_unmatched     = fs->make<TH1D>("h_pTTrCombi_unmatched","h_pTTrCombi_unmatched",1000,0.,500.);
  h_etaTrCombi_unmatched    = fs->make<TH1D>("h_etaTrCombi_unmatched","h_etaTrCombi_unmatched",60,-3.,3.);
  h_phiTrCombi_unmatched    = fs->make<TH1D>("h_phiTrCombi_unmatched","h_phiTrCombi_unmatched",64,-3.2,3.2);
  h_pTr1pTr2combi_unmatched   = fs->make<TH1D>("h_pTr1pTr2combi_unmatched","h_pTr1pTr2combi_unmatched",500,0.,10.);
  h_pTTr1pTTr2combi_unmatched = fs->make<TH1D>("h_pTTr1pTTr2combi_unmatched","h_pTTr1pTTr2combi_unmatched",500,0.,10.);

  h_D0combi_Mass  = fs->make<TH1D>("h_D0combi_Mass","h_D0combi_Mass",1000,0.,10.);
  h_D0combi_dRJet = fs->make<TH1D>("h_D0combi_dRJet","h_D0combi_dRJet",200,0.,1.);
  h_D0combi_p     = fs->make<TH1D>("h_D0combi_p","h_D0combi_p",1000,0.,500.);
  h_D0combi_pT    = fs->make<TH1D>("h_D0combi_pT","h_D0combi_pT",1000,0.,500.);
  h_D0combi_eta   = fs->make<TH1D>("h_D0combi_eta","h_D0combi_eta",60,-3.,3.);
  h_D0combi_phi   = fs->make<TH1D>("h_D0combi_phi","h_D0combi_phi",64,-3.2,3.2);
  h_D0combi_paired_Mass  = fs->make<TH1D>("h_D0combi_paired_Mass","h_D0combi_paired_Mass",1000,0.,10.);
  h_D0combi_paired_dRJet = fs->make<TH1D>("h_D0combi_paired_dRJet","h_D0combi_paired_dRJet",200,0.,1.);
  h_D0combi_paired_p     = fs->make<TH1D>("h_D0combi_paired_p","h_D0combi_paired_p",1000,0.,500.);
  h_D0combi_paired_pT    = fs->make<TH1D>("h_D0combi_paired_pT","h_D0combi_paired_pT",1000,0.,500.);
  h_D0combi_paired_eta   = fs->make<TH1D>("h_D0combi_paired_eta","h_D0combi_paired_eta",60,-3.,3.);
  h_D0combi_paired_phi   = fs->make<TH1D>("h_D0combi_paired_phi","h_D0combi_paired_phi",64,-3.2,3.2);
  h_D0combi_switched_Mass  = fs->make<TH1D>("h_D0combi_switched_Mass","h_D0combi_switched_Mass",1000,0.,10.);
  h_D0combi_switched_dRJet = fs->make<TH1D>("h_D0combi_switched_dRJet","h_D0combi_switched_dRJet",200,0.,1.);
  h_D0combi_switched_p     = fs->make<TH1D>("h_D0combi_switched_p","h_D0combi_switched_p",1000,0.,500.);
  h_D0combi_switched_pT    = fs->make<TH1D>("h_D0combi_switched_pT","h_D0combi_switched_pT",1000,0.,500.);
  h_D0combi_switched_eta   = fs->make<TH1D>("h_D0combi_switched_eta","h_D0combi_switched_eta",60,-3.,3.);
  h_D0combi_switched_phi   = fs->make<TH1D>("h_D0combi_switched_phi","h_D0combi_switched_phi",64,-3.2,3.2);
  h_D0combi_unmatched_Mass  = fs->make<TH1D>("h_D0combi_unmatched_Mass","h_D0combi_unmatched_Mass",1000,0.,10.);
  h_D0combi_unmatched_dRJet = fs->make<TH1D>("h_D0combi_unmatched_dRJet","h_D0combi_unmatched_dRJet",200,0.,1.);
  h_D0combi_unmatched_p     = fs->make<TH1D>("h_D0combi_unmatched_p","h_D0combi_unmatched_p",1000,0.,500.);
  h_D0combi_unmatched_pT    = fs->make<TH1D>("h_D0combi_unmatched_pT","h_D0combi_unmatched_pT",1000,0.,500.);
  h_D0combi_unmatched_eta   = fs->make<TH1D>("h_D0combi_unmatched_eta","h_D0combi_unmatched_eta",60,-3.,3.);
  h_D0combi_unmatched_phi   = fs->make<TH1D>("h_D0combi_unmatched_phi","h_D0combi_unmatched_phi",64,-3.2,3.2);

  h_BCand_DeltaRD0combiMu = fs->make<TH1D>("h_BCand_DeltaRD0combiMu","h_BCand_DeltaRD0combiMu",100,0.,5.);
  h_B_D0combiMass     = fs->make<TH1D>("h_B_D0combiMass","h_B_D0combiMass",1000,0.,10.);
  h_B_DeltaRD0combiMu = fs->make<TH1D>("h_B_DeltaRD0combiMu","h_B_DeltaRD0combiMu",80,0.,0.4);
  h_B_pD0combipMu     = fs->make<TH1D>("h_B_pD0combipMu","h_B_pD0combipMu",500,0.,10.);
  h_B_pTD0combipTMu   = fs->make<TH1D>("h_B_pTD0combipTMu","h_B_pTD0combipTMu",500,0.,10.);
  h_B_combiMass       = fs->make<TH1D>("h_B_combiMass","h_B_combiMass",1000,0.,10.);
  h_B_combip          = fs->make<TH1D>("h_B_combip","h_B_combip",1000,0.,500.);
  h_B_combipT         = fs->make<TH1D>("h_B_combipT","h_B_combipT",1000,0.,500.);
  h_B_combieta        = fs->make<TH1D>("h_B_combieta","h_B_combieta",60,-3.,3.);
  h_B_combiphi        = fs->make<TH1D>("h_B_combiphi","h_B_combiphi",64,-3.2,3.2);
  h_B_paired_D0combiMass     = fs->make<TH1D>("h_B_paired_D0combiMass","h_B_paired_D0combiMass",1000,0.,10.);
  h_B_paired_DeltaRD0combiMu = fs->make<TH1D>("h_B_paired_DeltaRD0combiMu","h_B_paired_DeltaRD0combiMu",80,0.,0.4);
  h_B_paired_pD0combipMu     = fs->make<TH1D>("h_B_paired_pD0combipMu","h_B_paired_pD0combipMu",500,0.,10.);
  h_B_paired_pTD0combipTMu   = fs->make<TH1D>("h_B_paired_pTD0combipTMu","h_B_paired_pTD0combipTMu",500,0.,10.);
  h_B_paired_combiMass       = fs->make<TH1D>("h_B_paired_combiMass","h_B_paired_combiMass",1000,0.,10.);
  h_B_paired_combip          = fs->make<TH1D>("h_B_paired_combip","h_B_paired_combip",1000,0.,500.);
  h_B_paired_combipT         = fs->make<TH1D>("h_B_paired_combipT","h_B_paired_combipT",1000,0.,500.);
  h_B_paired_combieta        = fs->make<TH1D>("h_B_paired_combieta","h_B_paired_combieta",60,-3.,3.);
  h_B_paired_combiphi        = fs->make<TH1D>("h_B_paired_combiphi","h_B_paired_combiphi",64,-3.2,3.2);
  h_B_switched_D0combiMass     = fs->make<TH1D>("h_B_switched_D0combiMass","h_B_switched_D0combiMass",1000,0.,10.);
  h_B_switched_DeltaRD0combiMu = fs->make<TH1D>("h_B_switched_DeltaRD0combiMu","h_B_switched_DeltaRD0combiMu",80,0.,0.4);
  h_B_switched_pD0combipMu     = fs->make<TH1D>("h_B_switched_pD0combipMu","h_B_switched_pD0combipMu",500,0.,10.);
  h_B_switched_pTD0combipTMu   = fs->make<TH1D>("h_B_switched_pTD0combipTMu","h_B_switched_pTD0combipTMu",500,0.,10.);
  h_B_switched_combiMass       = fs->make<TH1D>("h_B_switched_combiMass","h_B_switched_combiMass",1000,0.,10.);
  h_B_switched_combip          = fs->make<TH1D>("h_B_switched_combip","h_B_switched_combip",1000,0.,500.);
  h_B_switched_combipT         = fs->make<TH1D>("h_B_switched_combipT","h_B_switched_combipT",1000,0.,500.);
  h_B_switched_combieta        = fs->make<TH1D>("h_B_switched_combieta","h_B_switched_combieta",60,-3.,3.);
  h_B_switched_combiphi        = fs->make<TH1D>("h_B_switched_combiphi","h_B_switched_combiphi",64,-3.2,3.2);
  h_B_unmatched_D0combiMass     = fs->make<TH1D>("h_B_unmatched_D0combiMass","h_B_unmatched_D0combiMass",1000,0.,10.);
  h_B_unmatched_DeltaRD0combiMu = fs->make<TH1D>("h_B_unmatched_DeltaRD0combiMu","h_B_unmatched_DeltaRD0combiMu",80,0.,0.4);
  h_B_unmatched_pD0combipMu     = fs->make<TH1D>("h_B_unmatched_pD0combipMu","h_B_unmatched_pD0combipMu",500,0.,10.);
  h_B_unmatched_pTD0combipTMu   = fs->make<TH1D>("h_B_unmatched_pTD0combipTMu","h_B_unmatched_pTD0combipTMu",500,0.,10.);
  h_B_unmatched_combiMass       = fs->make<TH1D>("h_B_unmatched_combiMass","h_B_unmatched_combiMass",1000,0.,10.);
  h_B_unmatched_combip          = fs->make<TH1D>("h_B_unmatched_combip","h_B_unmatched_combip",1000,0.,500.);
  h_B_unmatched_combipT         = fs->make<TH1D>("h_B_unmatched_combipT","h_B_unmatched_combipT",1000,0.,500.);
  h_B_unmatched_combieta        = fs->make<TH1D>("h_B_unmatched_combieta","h_B_unmatched_combieta",60,-3.,3.);
  h_B_unmatched_combiphi        = fs->make<TH1D>("h_B_unmatched_combiphi","h_B_unmatched_combiphi",64,-3.2,3.2);

  // optimized invariant combination
  h_D0optcombiCand_pT = fs->make<TH1D>("h_D0optcombiCand_pT","h_D0optcombiCand_pT",1000,0.,500.);

  h_D0optcombi_Mass  = fs->make<TH1D>("h_D0optcombi_Mass","h_D0optcombi_Mass",1000,0.,10.);
  h_D0optcombi_dRJet = fs->make<TH1D>("h_D0optcombi_dRJet","h_D0optcombi_dRJet",200,0.,1.);
  h_D0optcombi_p     = fs->make<TH1D>("h_D0optcombi_p","h_D0optcombi_p",1000,0.,500.);
  h_D0optcombi_pT    = fs->make<TH1D>("h_D0optcombi_pT","h_D0optcombi_pT",1000,0.,500.);
  h_D0optcombi_eta   = fs->make<TH1D>("h_D0optcombi_eta","h_D0optcombi_eta",60,-3.,3.);
  h_D0optcombi_phi   = fs->make<TH1D>("h_D0optcombi_phi","h_D0optcombi_phi",64,-3.2,3.2);
  h_D0optcombi_paired_Mass  = fs->make<TH1D>("h_D0optcombi_paired_Mass","h_D0optcombi_paired_Mass",1000,0.,10.);
  h_D0optcombi_paired_dRJet = fs->make<TH1D>("h_D0optcombi_paired_dRJet","h_D0optcombi_paired_dRJet",200,0.,1.);
  h_D0optcombi_paired_p     = fs->make<TH1D>("h_D0optcombi_paired_p","h_D0optcombi_paired_paired_p",1000,0.,500.);
  h_D0optcombi_paired_pT    = fs->make<TH1D>("h_D0optcombi_paired_pT","h_D0optcombi_paired_pT",1000,0.,500.);
  h_D0optcombi_paired_eta   = fs->make<TH1D>("h_D0optcombi_paired_eta","h_D0optcombi_paired_eta",60,-3.,3.);
  h_D0optcombi_paired_phi   = fs->make<TH1D>("h_D0optcombi_paired_phi","h_D0optcombi_paired_phi",64,-3.2,3.2);
  h_D0optcombi_switched_Mass  = fs->make<TH1D>("h_D0optcombi_switched_Mass","h_D0optcombi_switched_Mass",1000,0.,10.);
  h_D0optcombi_switched_dRJet = fs->make<TH1D>("h_D0optcombi_switched_dRJet","h_D0optcombi_switched_dRJet",200,0.,1.);
  h_D0optcombi_switched_p     = fs->make<TH1D>("h_D0optcombi_switched_p","h_D0optcombi_switched_switched_p",1000,0.,500.);
  h_D0optcombi_switched_pT    = fs->make<TH1D>("h_D0optcombi_switched_pT","h_D0optcombi_switched_pT",1000,0.,500.);
  h_D0optcombi_switched_eta   = fs->make<TH1D>("h_D0optcombi_switched_eta","h_D0optcombi_switched_eta",60,-3.,3.);
  h_D0optcombi_switched_phi   = fs->make<TH1D>("h_D0optcombi_switched_phi","h_D0optcombi_switched_phi",64,-3.2,3.2);
  h_D0optcombi_unmatched_Mass  = fs->make<TH1D>("h_D0optcombi_unmatched_Mass","h_D0optcombi_unmatched_Mass",1000,0.,10.);
  h_D0optcombi_unmatched_dRJet = fs->make<TH1D>("h_D0optcombi_unmatched_dRJet","h_D0optcombi_unmatched_dRJet",200,0.,1.);
  h_D0optcombi_unmatched_p     = fs->make<TH1D>("h_D0optcombi_unmatched_p","h_D0optcombi_unmatched_unmatched_p",1000,0.,500.);
  h_D0optcombi_unmatched_pT    = fs->make<TH1D>("h_D0optcombi_unmatched_pT","h_D0optcombi_unmatched_pT",1000,0.,500.);
  h_D0optcombi_unmatched_eta   = fs->make<TH1D>("h_D0optcombi_unmatched_eta","h_D0optcombi_unmatched_eta",60,-3.,3.);
  h_D0optcombi_unmatched_phi   = fs->make<TH1D>("h_D0optcombi_unmatched_phi","h_D0optcombi_unmatched_phi",64,-3.2,3.2);

  h_BCand_DeltaRD0optcombiMu = fs->make<TH1D>("h_BCand_DeltaRD0optcombiMu","h_BCand_DeltaRD0optcombiMu",100,0.,5.);
  h_B_D0optcombiMass     = fs->make<TH1D>("h_B_D0optcombiMass","h_B_D0optcombiMass",1000,0.,10.);
  h_B_DeltaRD0optcombiMu = fs->make<TH1D>("h_B_DeltaRD0optcombiMu","h_B_DeltaRD0optcombiMu",80,0.,0.4);
  h_B_pD0optcombipMu     = fs->make<TH1D>("h_B_pD0optcombipMu","h_B_pD0optcombipMu",500,0.,10.);
  h_B_pTD0optcombipTMu   = fs->make<TH1D>("h_B_pTD0optcombipTMu","h_B_pTD0optcombipTMu",500,0.,10.);
  h_B_optcombiMass       = fs->make<TH1D>("h_B_optcombiMass","h_B_optcombiMass",1000,0.,10.);
  h_B_optcombip          = fs->make<TH1D>("h_B_optcombip","h_B_optcombip",1000,0.,500.);
  h_B_optcombipT         = fs->make<TH1D>("h_B_optcombipT","h_B_optcombipT",1000,0.,500.);
  h_B_optcombieta        = fs->make<TH1D>("h_B_optcombieta","h_B_optcombieta",60,-3.,3.);
  h_B_optcombiphi        = fs->make<TH1D>("h_B_optcombiphi","h_B_optcombiphi",64,-3.2,3.2);
  h_B_paired_D0optcombiMass     = fs->make<TH1D>("h_B_paired_D0optcombiMass","h_B_paired_D0optcombiMass",1000,0.,10.);
  h_B_paired_DeltaRD0optcombiMu = fs->make<TH1D>("h_B_paired_DeltaRD0optcombiMu","h_B_paired_DeltaRD0optcombiMu",80,0.,0.4);
  h_B_paired_pD0optcombipMu     = fs->make<TH1D>("h_B_paired_pD0optcombipMu","h_B_paired_pD0optcombipMu",500,0.,10.);
  h_B_paired_pTD0optcombipTMu   = fs->make<TH1D>("h_B_paired_pTD0optcombipTMu","h_B_paired_pTD0optcombipTMu",500,0.,10.);
  h_B_paired_optcombiMass       = fs->make<TH1D>("h_B_paired_optcombiMass","h_B_paired_optcombiMass",1000,0.,10.);
  h_B_paired_optcombip          = fs->make<TH1D>("h_B_paired_optcombip","h_B_paired_optcombip",1000,0.,500.);
  h_B_paired_optcombipT         = fs->make<TH1D>("h_B_paired_optcombipT","h_B_paired_optcombipT",1000,0.,500.);
  h_B_paired_optcombieta        = fs->make<TH1D>("h_B_paired_optcombieta","h_B_paired_optcombieta",60,-3.,3.);
  h_B_paired_optcombiphi        = fs->make<TH1D>("h_B_paired_optcombiphi","h_B_paired_optcombiphi",64,-3.2,3.2);
  h_B_switched_D0optcombiMass     = fs->make<TH1D>("h_B_switched_D0optcombiMass","h_B_switched_D0optcombiMass",1000,0.,10.);
  h_B_switched_DeltaRD0optcombiMu = fs->make<TH1D>("h_B_switched_DeltaRD0optcombiMu","h_B_switched_DeltaRD0optcombiMu",80,0.,0.4);
  h_B_switched_pD0optcombipMu     = fs->make<TH1D>("h_B_switched_pD0optcombipMu","h_B_switched_pD0optcombipMu",500,0.,10.);
  h_B_switched_pTD0optcombipTMu   = fs->make<TH1D>("h_B_switched_pTD0optcombipTMu","h_B_switched_pTD0optcombipTMu",500,0.,10.);
  h_B_switched_optcombiMass       = fs->make<TH1D>("h_B_switched_optcombiMass","h_B_switched_optcombiMass",1000,0.,10.);
  h_B_switched_optcombip          = fs->make<TH1D>("h_B_switched_optcombip","h_B_switched_optcombip",1000,0.,500.);
  h_B_switched_optcombipT         = fs->make<TH1D>("h_B_switched_optcombipT","h_B_switched_optcombipT",1000,0.,500.);
  h_B_switched_optcombieta        = fs->make<TH1D>("h_B_switched_optcombieta","h_B_switched_optcombieta",60,-3.,3.);
  h_B_switched_optcombiphi        = fs->make<TH1D>("h_B_switched_optcombiphi","h_B_switched_optcombiphi",64,-3.2,3.2);
  h_B_unmatched_D0optcombiMass     = fs->make<TH1D>("h_B_unmatched_D0optcombiMass","h_B_unmatched_D0optcombiMass",1000,0.,10.);
  h_B_unmatched_DeltaRD0optcombiMu = fs->make<TH1D>("h_B_unmatched_DeltaRD0optcombiMu","h_B_unmatched_DeltaRD0optcombiMu",80,0.,0.4);
  h_B_unmatched_pD0optcombipMu     = fs->make<TH1D>("h_B_unmatched_pD0optcombipMu","h_B_unmatched_pD0optcombipMu",500,0.,10.);
  h_B_unmatched_pTD0optcombipTMu   = fs->make<TH1D>("h_B_unmatched_pTD0optcombipTMu","h_B_unmatched_pTD0optcombipTMu",500,0.,10.);
  h_B_unmatched_optcombiMass       = fs->make<TH1D>("h_B_unmatched_optcombiMass","h_B_unmatched_optcombiMass",1000,0.,10.);
  h_B_unmatched_optcombip          = fs->make<TH1D>("h_B_unmatched_optcombip","h_B_unmatched_optcombip",1000,0.,500.);
  h_B_unmatched_optcombipT         = fs->make<TH1D>("h_B_unmatched_optcombipT","h_B_unmatched_optcombipT",1000,0.,500.);
  h_B_unmatched_optcombieta        = fs->make<TH1D>("h_B_unmatched_optcombieta","h_B_unmatched_optcombieta",60,-3.,3.);
  h_B_unmatched_optcombiphi        = fs->make<TH1D>("h_B_unmatched_optcombiphi","h_B_unmatched_optcombiphi",64,-3.2,3.2);
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
KalmanAnalyzer::endJob() 
{

  //  std::cout << "Closing histos..." << std::endl;

}

// ------------ method called when starting to processes a run  ------------
  void 
KalmanAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
  void 
KalmanAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
  void 
KalmanAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
  void 
KalmanAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
KalmanAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KalmanAnalyzer);
