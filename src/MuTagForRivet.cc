// -*- C++ -*-
//
// Package:    MuTagForRivet
// Class:      MuTagForRivet
// 
/**\class MuTagForRivet MuTagForRivet.cc UserCode/KalmanAnalyzer/src/MuTagForRivet.cc

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

// #include "DataFormats/PatCandidates/interface/PFParticle.h"

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
#include "TH2D.h"
#include "TTree.h"
#include "TLorentzVector.h"

//
// class declaration
//

class MuTagForRivet : public edm::EDAnalyzer {
  public:
    explicit MuTagForRivet(const edm::ParameterSet&);
    ~MuTagForRivet();

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
    bool _isCSVbased;

    // evts properties
    unsigned int nEvts;
    unsigned int nEvts2;
    double nEvts3;
    bool hasGoodMuons;
    bool hasGoodElectrons;
    bool hasGoodLeptons;
    bool hasGoodJets;
    bool isGoodEvt;

    TH1D* _h_idPFel;
    TH2D* _h_id_pT_PFel;
    TH2D* _h_dR_dpT_PFel;
    TH1D* _h_idPFmu;
    TH2D* _h_id_pT_PFmu;
    TH2D* _h_dR_dpT_PFmu;
    TH1D* _h_id1stPFmu;
    TH2D* _h_id_pT_1stPFmu;
    TH2D* _h_dR_dpT_1stPFmu;
    TH2D* _h_origin_id_1stPFmu_Inf100;
    TH2D* _h_eta_dpz_1stPFmu_Inf100;
    TH2D* _h_origin_id_1stPFmu_Sup100;
    TH2D* _h_eta_dpz_1stPFmu_Sup100;
    TH2D* _h_genId_pT_NonSelJets;

    TH1D* _h_nVtx;

    TH1D* _h_CSVSelJets;
    TH1D* _h_pTSelJets;
    TH1D* _h_etaSelJets;
    TH1D* _h_genIdSelJets;
    TH1D* _h_etach;
    TH1D* _h_pTch;

    TH1D* _h_Nch;
    TH1D* _h_sump;
    TLorentzVector _sumpvec;
    TH1D* _h_sumpvec;
    TH1D* _h_sum1p;
    TH1D* _h_sum2p;
    TH1D* _h_sum3p;
    TH1D* _h_mass3;
    TH1D* _h_R1;
    TH2D* _h_R1_Nch;
    TH1D* _h_R3;
    TH2D* _h_R3_Nch;
    TH1D* _h_sum1p_nomu;
    TH1D* _h_sum2p_nomu;
    TH1D* _h_sum3p_nomu;
    TH1D* _h_mass3_nomu;
    TH1D* _h_R1_nomu;
    TH2D* _h_R1_Nch_nomu;
    TH1D* _h_R3_nomu;
    TH2D* _h_R3_Nch_nomu;
    TH1D* _h_D0Mass;
    TH1D* _h_D0p;
    TH1D* _h_D0pT;
    TH1D* _h_D0eta;
    TH1D* _h_BMomentum_unbiased;
    TH1D* _h_D0MassClean;
    TH1D* _h_D0pClean;
    TH1D* _h_D0pTClean;
    TH1D* _h_D0etaClean;
    TH1D* _h_BMomentum;
    TH1D* _h_mup;
    TH1D* _h_BMomentumClean;

    TTree* _t_bjets;
    double weight;
    double _CSV;
    double _vecP[4];
    int _Nch;
    double _sump;
    double _sumpt;
    double _tr1[4];
    double _tr2[4];
    double _tr3[4];
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
MuTagForRivet::MuTagForRivet(const edm::ParameterSet& iConfig) :
  _isCSVbased(iConfig.getUntrackedParameter<bool>("isCSVbased", false))
{
  // now do what ever initialization is needed
  nEvts = 0;
  nEvts2 = 0;
  nEvts3 = 0.;

}


MuTagForRivet::~MuTagForRivet()
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
MuTagForRivet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    double lumitab[4]= {888.7,4446,7021 ,7221}; // when running on all the stat
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

    if (_isCSVbased) {
      double vtxWeight[45] = {23./16.24817, 157./96.10187, 608./389.7191, 1653./1092.023, 3518./2428.246, 6185./4495.628, 9848./7287.974, 13907./10632.4, 17801./14183.71, 21165./17578., 24344./20656.52, 26132./22919.96, 26440./24470.79, 26610./25254.18, 25605./25303.6, 23974./24804.25, 21937./23480.46, 19587./21797.11, 17451./19784.64, 14764./17591.69, 12440./15221.65, 10121./12837., 8188./10599.46, 6449./6746.303, 4873./8565.719, 3690./5220.74, 2741./3919.616, 1977./2855.812, 1451./2064.941, 990./1440.364, 656./1021.42, 479./686.967, 315./451.8202, 210./290.9677, 140./184.5276, 71./118.4836, 72./73.5648, 41./44.23643, 23./26.39862, 13./16.03296, 7./9.068202, 10./6.257565, 4./3.187427, 1./1.576181, 1./1.232338};
      if (vtx.size() < 46)
        weight *= vtxWeight[vtx.size()];
      else 
        weight *= vtxWeight[44];
    }
    else {
      double vtxWeight[45] = {23./16.24817, 157./96.10187, 608./389.7191, 1653./1092.023, 3518./2428.246, 6185./4495.628, 9848./7287.974, 13907./10632.4, 17801./14183.71, 21165./17578., 24344./20656.52, 26132./22919.96, 26440./24470.79, 26610./25254.18, 25605./25303.6, 23974./24804.25, 21937./23480.46, 19587./21797.11, 17451./19784.64, 14764./17591.69, 12440./15221.65, 10121./12837., 8188./10599.46, 6449./8565.719, 4873./6746.303, 3690./5220.74, 2741./3919.616, 1977./2855.812, 1451./2064.941, 990./1440.364, 656./1021.42, 479./686.967, 315./451.8202, 210./290.9677, 140./184.5276, 71./118.4836, 72./73.5648, 41./44.23643, 23./26.39862, 13./16.03296, 7./9.068202, 10./6.257565, 4./3.187427, 1./1.576181, 1./1.232338};
      if (vtx.size() < 46)
        weight *= vtxWeight[vtx.size()-1];
      else 
        weight *= vtxWeight[44];
    }

    nEvts3 = nEvts3 + weight;
    _h_nVtx->Fill((double)vtx.size(), weight);

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

    // Select good PF muons and electrons :

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

    // Match PF e/mu to genParticles
    edm::Handle<std::vector<reco::GenParticle> > mcHandle;
    edm::InputTag tagGen("genParticles", "", "SIM");
    iEvent.getByLabel(tagGen, mcHandle);
    const reco::GenParticleCollection& mcparts = *(mcHandle.product());
    for (unsigned int iPFel = 0; iPFel < myPFel.size(); iPFel++) {
      TLorentzVector p_PFel;
      p_PFel.SetPtEtaPhiM(myPFel[iPFel]->pt(), myPFel[iPFel]->eta(), myPFel[iPFel]->phi(), 0.);
      double minDeltaRPFelGen = 200.;
      double minDeltaPtPFelGen = 200.;
      int idGen = 0;
      double pTGen = 0.;
      for (reco::GenParticleCollection::const_iterator it = mcparts.begin(); it != mcparts.end(); it++) {
        if (fabs((*it).pt()) < 1e-10) continue;
        TLorentzVector p_Gen;
        p_Gen.SetPtEtaPhiM((*it).pt(), (*it).eta(), (*it).phi(), 0.);
        if (p_Gen.DeltaR(p_PFel) < minDeltaRPFelGen && fabs(p_Gen.Pt()-p_PFel.Pt()) < 0.15*p_PFel.Pt()) {
          minDeltaRPFelGen = p_Gen.DeltaR(p_PFel);
          minDeltaPtPFelGen = fabs(p_Gen.Pt()-p_PFel.Pt());
          idGen = abs((*it).pdgId());
          pTGen = (*it).pt();
        }
      }
      if (abs(idGen) > 0) {
        _h_idPFel->Fill((double)idGen, weight);
        _h_id_pT_PFel->Fill((double)idGen, pTGen, weight);
        _h_dR_dpT_PFel->Fill(minDeltaRPFelGen, minDeltaPtPFelGen, weight);
      }
    }
    for (unsigned int iPFmu = 0; iPFmu < myPFmu.size(); iPFmu++) {
      TLorentzVector p_PFmu;
      p_PFmu.SetPtEtaPhiM(myPFmu[iPFmu]->pt(), myPFmu[iPFmu]->eta(), myPFmu[iPFmu]->phi(), 0.);
      double minDeltaRPFmuGen = 200.;
      double minDeltaPtPFmuGen = 200.;
      int idGen = 0;
      double pTGen = 0.;
      for (reco::GenParticleCollection::const_iterator it = mcparts.begin(); it != mcparts.end(); it++) {
        if (fabs((*it).pt()) < 1e-10) continue;
        TLorentzVector p_Gen;
        p_Gen.SetPtEtaPhiM((*it).pt(), (*it).eta(), (*it).phi(), 0.);
        if (p_Gen.DeltaR(p_PFmu) < minDeltaRPFmuGen && fabs(p_Gen.Pt()-p_PFmu.Pt()) < 0.15*p_PFmu.Pt()) {
          minDeltaRPFmuGen = p_Gen.DeltaR(p_PFmu);
          minDeltaPtPFmuGen = fabs(p_Gen.Pt()-p_PFmu.Pt());
          idGen = abs((*it).pdgId());
          pTGen = (*it).pt();
        }
      }
      if (abs(idGen) > 0) {
        _h_idPFmu->Fill((double)idGen, weight);
        _h_id_pT_PFmu->Fill((double)idGen, pTGen, weight);
        _h_dR_dpT_PFmu->Fill(minDeltaRPFmuGen, minDeltaPtPFmuGen, weight);
      }
    }

    //-------------------------------------------
    // Access the jets 
    //-------------------------------------------

    // double d0mass = 1.86484;
    // double bmass = 5.2796;

    ParticleMass gMassD0 = 1.86483;
    // float        gSigmaD0 = 0.00014;

    // ParticleMass gMassW = 80.399;
    // float        gSigmaW = 0.023;
    // float        gResoW = 10.;

    ParticleMass gMassMu  = 0.105658367;
    // float        gSigmaMu = 0.000000004; 

    ParticleMass gMassK  = 0.493677;
    // float        gSigmaK = 0.000001;

    ParticleMass gMassPi  = 0.13957018;
    // float        gSigmaPi = 0.00000001;

    // Track setup

    edm::ESHandle<TransientTrackBuilder> theB;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);  

    edm::InputTag jetTag("selectedPatJetsPFlow","","PAT");

    edm::Handle<pat::JetCollection> jetsHandle;
    iEvent.getByLabel(jetTag, jetsHandle);
    pat::JetCollection jets = *jetsHandle;

    /*
       if ( (*it).pt() >= 30. ) {
       */    
    _CSV = -1; 
    _vecP[0] = -100.; _vecP[1] = -100.; _vecP[2] = -100.; _vecP[3] = -100.;
    _Nch = -1; 
    _sump = -1.; 
    _sumpt = -1.; 
    _sumpvec.SetPtEtaPhiM(0.,0.,0.,0.);
    _tr1[0] = -100.; _tr1[1] = -100.; _tr1[2] = -100.; _tr1[3] = -1000.;
    _tr2[0] = -100.; _tr2[1] = -100.; _tr2[2] = -100.; _tr2[3] = -1000.;
    _tr3[0] = -100.; _tr3[1] = -100.; _tr3[2] = -100.; _tr3[3] = -1000.;


    for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {

      // select jet with pT > 20 GeV/c and a muon inside
      bool hasMuInside = false;
      reco::TrackRefVector jetTracks = (*it).associatedTracks();
      for (reco::track_iterator iter1 = jetTracks.begin(); iter1 != jetTracks.end(); ++iter1) {

        const reco::Track& Track1 = **iter1;

        if ((**iter1).pt() < 4.) continue;
        if (!Track1.quality(reco::Track::highPurity)) continue;

        // look for muons 
        for (unsigned int iMuCand = 0; iMuCand < myPFmu.size(); iMuCand++) {
          TLorentzVector p_MuCand, p_trCand;
          p_MuCand.SetPtEtaPhiM(myPFmu[iMuCand]->pt(), myPFmu[iMuCand]->eta(), myPFmu[iMuCand]->phi(), gMassMu);
          p_trCand.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassPi);
          if (p_trCand.DeltaR(p_MuCand) < 0.0005) {
            hasMuInside = true;
            break;
          }
        }
      }
      if (!hasMuInside) continue;
      if ((*it).genParticle() &&  (*it).pt() < 30.)
        _h_genId_pT_NonSelJets->Fill((double)abs(((*it).genParticle())->pdgId()), (*it).pt(), weight); 
      if ((*it).pt() < 20.) continue;

      _CSV = (*it).bDiscriminator("combinedSecondaryVertexBJetTags");
      _vecP[0] = (*it).pt();
      _vecP[1] = (*it).eta();
      _vecP[2] = (*it).phi();
      _vecP[3] = (*it).mass();
      _h_CSVSelJets->Fill((*it).bDiscriminator("combinedSecondaryVertexBJetTags"), weight);
      _h_pTSelJets->Fill((*it).pt(), weight);
      _h_etaSelJets->Fill((*it).eta(), weight);
      if ((*it).genParticle())
        _h_genIdSelJets->Fill((double)abs(((*it).genParticle())->pdgId()), weight); 

      TLorentzVector p_Jet;
      p_Jet.SetPtEtaPhiM((*it).pt(), (*it).eta(), (*it).phi(), (*it).mass());

      double pt_trCand_nomu[3]  = {0., 0., 0.};
      double eta_trCand_nomu[3] = {0., 0., 0.};
      double phi_trCand_nomu[3] = {0., 0., 0.};
      double id_trCand_nomu[3] = {0., 0., 0.};
      TLorentzVector p_trCand_nomu[3];

      double pt_trCand[3]  = {0., 0., 0.};
      double eta_trCand[3] = {0., 0., 0.};
      double phi_trCand[3] = {0., 0., 0.};
      double id_trCand[3] = {0., 0., 0.};
      TLorentzVector p_trCand[3];

      _Nch = 0; 
      _sump = 0.; 
      _sumpt = 0.; 
      std::vector<const reco::PFCandidate*> myPFmuInSelJet;
      for (reco::track_iterator iter1 = jetTracks.begin(); iter1 != jetTracks.end(); ++iter1) {

        const reco::Track& Track1 = **iter1;

        if ((**iter1).pt() < 4.) continue;
        if (!Track1.quality(reco::Track::highPurity)) continue;

        // look for muons and electrons
        bool trCandIsMu = false;
        int iMuInSelJet = 0;
        for (unsigned int iMuCand = 0; iMuCand < myPFmu.size(); iMuCand++) {
          TLorentzVector p_MuCand, p_trCand;
          p_MuCand.SetPtEtaPhiM(myPFmu[iMuCand]->pt(), myPFmu[iMuCand]->eta(), myPFmu[iMuCand]->phi(), gMassMu);
          p_trCand.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassPi);
          if (p_trCand.DeltaR(p_MuCand) < 0.0005) {
            trCandIsMu = true;
            iMuInSelJet = iMuCand;
            break;
          }
        }
        if (trCandIsMu) myPFmuInSelJet.push_back(myPFmu[iMuInSelJet]);
        bool trCandIsEl = false;
        for (unsigned int iElCand = 0; iElCand < myPFel.size(); iElCand++) {
          TLorentzVector p_ElCand, p_trCand;
          p_ElCand.SetPtEtaPhiM(myPFel[iElCand]->pt(), myPFel[iElCand]->eta(), myPFel[iElCand]->phi(), 0.);
          p_trCand.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassPi);
          if (p_trCand.DeltaR(p_ElCand) < 0.005) {
            trCandIsEl = true;
            break;
          }
        }
        // store the 3 traces with highest pT
        if ((**iter1).pt() >= pt_trCand[0]) {
          pt_trCand[2]  = pt_trCand[1];
          eta_trCand[2] = eta_trCand[1];
          phi_trCand[2] = phi_trCand[1];
          id_trCand[2]  = id_trCand[1];
          pt_trCand[1]  = pt_trCand[0];
          eta_trCand[1] = eta_trCand[0];
          phi_trCand[1] = phi_trCand[0];
          id_trCand[1]  = id_trCand[0];
          pt_trCand[0]  = (**iter1).pt();
          eta_trCand[0] = (**iter1).eta();
          phi_trCand[0] = (**iter1).phi();
          id_trCand[0]  = (**iter1).charge();
          _tr3[3] = _tr2[3];
          _tr2[3] = _tr1[3];
          if (trCandIsEl) 
            _tr1[3] = (**iter1).charge()*11.;
          if (trCandIsMu)
            _tr1[3] = (**iter1).charge()*13.;
          if (!trCandIsMu && !trCandIsEl)
            _tr1[3] = (**iter1).charge()*211.;
        }
        else {
          if ((**iter1).pt() >= pt_trCand[1]) {
            pt_trCand[2]  = pt_trCand[1];
            eta_trCand[2] = eta_trCand[1];
            phi_trCand[2] = phi_trCand[1];
            id_trCand[2]  = id_trCand[1];
            pt_trCand[1]  = (**iter1).pt();
            eta_trCand[1] = (**iter1).eta();
            phi_trCand[1] = (**iter1).phi();
            id_trCand[1]  = (**iter1).charge();
            _tr3[3] = _tr3[3];
            if (trCandIsEl) 
              _tr2[3] = (**iter1).charge()*11.;
            if (trCandIsMu)
              _tr2[3] = (**iter1).charge()*13.;
            if (!trCandIsMu && !trCandIsEl)
              _tr2[3] = (**iter1).charge()*211.;
          }
          else  
            if ((**iter1).pt() >= pt_trCand[2]) {
              pt_trCand[2]  = (**iter1).pt();
              eta_trCand[2] = (**iter1).eta();
              phi_trCand[2] = (**iter1).phi();
              id_trCand[2]  = (**iter1).charge();
              if (trCandIsEl) 
                _tr3[3] = (**iter1).charge()*11.;
              if (trCandIsMu)
                _tr3[3] = (**iter1).charge()*13.;
              if (!trCandIsMu && !trCandIsEl)
                _tr3[3] = (**iter1).charge()*211.;
            }
        }
        if (fabs(_tr1[3]) < 999) {
          _tr1[0] = pt_trCand[0];
          _tr1[1] = eta_trCand[0];
          _tr1[2] = phi_trCand[0];
          if (fabs(_tr2[3]) < 999) {
            _tr2[0] = pt_trCand[1];
            _tr2[1] = eta_trCand[1];
            _tr2[2] = phi_trCand[1];
            if (fabs(_tr3[3]) < 999) {
              _tr3[0] = pt_trCand[2];
              _tr3[1] = eta_trCand[2];
              _tr3[2] = phi_trCand[2];
            }
          }
        }
        // store the 3 traces with highest pT, excluding muons
        if (!trCandIsMu) {
          if ((**iter1).pt() >= pt_trCand_nomu[0]) {
            pt_trCand_nomu[2]  = pt_trCand_nomu[1];
            eta_trCand_nomu[2] = eta_trCand_nomu[1];
            phi_trCand_nomu[2] = phi_trCand_nomu[1];
            id_trCand_nomu[2]  = id_trCand_nomu[1];
            pt_trCand_nomu[1]  = pt_trCand_nomu[0];
            eta_trCand_nomu[1] = eta_trCand_nomu[0];
            phi_trCand_nomu[1] = phi_trCand_nomu[0];
            id_trCand_nomu[1]  = id_trCand_nomu[0];
            pt_trCand_nomu[0]  = (**iter1).pt();
            eta_trCand_nomu[0] = (**iter1).eta();
            phi_trCand_nomu[0] = (**iter1).phi();
            id_trCand_nomu[0]  = (**iter1).charge();
          }
          else {
            if ((**iter1).pt() >= pt_trCand_nomu[1]) {
              pt_trCand_nomu[2]  = pt_trCand_nomu[1];
              eta_trCand_nomu[2] = eta_trCand_nomu[1];
              phi_trCand_nomu[2] = phi_trCand_nomu[1];
              id_trCand_nomu[2]  = id_trCand_nomu[1];
              pt_trCand_nomu[1]  = (**iter1).pt();
              eta_trCand_nomu[1] = (**iter1).eta();
              phi_trCand_nomu[1] = (**iter1).phi();
              id_trCand_nomu[1]  = (**iter1).charge();
            }
            else  
              if ((**iter1).pt() >= pt_trCand_nomu[2]) {
                pt_trCand_nomu[2]  = (**iter1).pt();
                eta_trCand_nomu[2] = (**iter1).eta();
                phi_trCand_nomu[2] = (**iter1).phi();
                id_trCand_nomu[2]  = (**iter1).charge();
              }
          }
        }

        TLorentzVector p_tr1;
        p_tr1.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassPi);
        if (fabs(pt_trCand[0]) > 1e-10) {
          p_trCand[0].SetPtEtaPhiM(pt_trCand[0], eta_trCand[0], phi_trCand[0], gMassPi);
          if (fabs(pt_trCand[1]) > 1e-10 && fabs(pt_trCand[2]) > 1e-10) {
            p_trCand[1].SetPtEtaPhiM(pt_trCand[1], eta_trCand[1], phi_trCand[1], gMassPi);
            p_trCand[2].SetPtEtaPhiM(pt_trCand[2], eta_trCand[2], phi_trCand[2], gMassPi);
          }
          else {
            p_trCand[1].SetPtEtaPhiM(0., 0., 0., 0.);
            p_trCand[2].SetPtEtaPhiM(0., 0., 0., 0.);
          }
        }
        else p_trCand[0].SetPtEtaPhiM(0., 0., 0., 0.);
        if (fabs(pt_trCand_nomu[0]) > 1e-10) {
          p_trCand_nomu[0].SetPtEtaPhiM(pt_trCand_nomu[0], eta_trCand_nomu[0], phi_trCand_nomu[0], gMassPi);
          if (fabs(pt_trCand_nomu[1]) > 1e-10 && fabs(pt_trCand_nomu[2]) > 1e-10) {
            p_trCand_nomu[1].SetPtEtaPhiM(pt_trCand_nomu[1], eta_trCand_nomu[1], phi_trCand_nomu[1], gMassPi);
            p_trCand_nomu[2].SetPtEtaPhiM(pt_trCand_nomu[2], eta_trCand_nomu[2], phi_trCand_nomu[2], gMassPi);
          }
          else {
            p_trCand_nomu[1].SetPtEtaPhiM(0., 0., 0., 0.);
            p_trCand_nomu[2].SetPtEtaPhiM(0., 0., 0., 0.);
          }
        }
        else p_trCand_nomu[0].SetPtEtaPhiM(0., 0., 0., 0.);

        _Nch++;
        _sump = _sump + p_tr1.P();
        _sumpt = _sumpt + p_tr1.Pt();
        _sumpvec = _sumpvec + p_tr1;
        _h_etach->Fill(p_tr1.Eta(), weight);
        _h_pTch->Fill(p_tr1.Pt(), weight);

      } // 1st jet's track loop

      _h_Nch->Fill((double)_Nch, weight);
      _h_sump->Fill(_sump, weight);
      _h_sumpvec->Fill(_sumpvec.P(), weight);
      if (p_trCand[0].M() > 1e-10) {
        _h_sum1p->Fill(p_trCand[0].P(), weight);
        _h_R1->Fill(p_trCand[0].P() / _sump, weight);
        _h_R1_Nch->Fill(p_trCand[0].P() / _sump, (double)_Nch, weight);
        if (p_trCand[1].M() > 1e-10) {
          _h_sum2p->Fill(p_trCand[0].P() + p_trCand[1].P(), weight);
          if (p_trCand[2].M() > 1e-10) {
            _h_sum3p->Fill(p_trCand[0].P() + p_trCand[1].P() + p_trCand[2].P(), weight);
            _h_mass3->Fill((p_trCand[0] + p_trCand[1] + p_trCand[2]).M(), weight);
            _h_R3->Fill((p_trCand[0].P() + p_trCand[1].P() + p_trCand[2].P()) / _sump, weight);
            _h_R3_Nch->Fill((p_trCand[0].P() + p_trCand[1].P() + p_trCand[2].P()) / _sump, (double)_Nch, weight);
          }
        }
      }
      if (p_trCand_nomu[0].M() > 1e-10) {
        _h_sum1p_nomu->Fill(p_trCand_nomu[0].P(), weight);
        _h_R1_nomu->Fill(p_trCand_nomu[0].P() / _sump, weight);
        _h_R1_Nch_nomu->Fill(p_trCand_nomu[0].P() / _sump, (double)_Nch, weight);
        if (p_trCand_nomu[1].M() > 1e-10) {
          _h_sum2p_nomu->Fill(p_trCand_nomu[0].P() + p_trCand_nomu[1].P(), weight);
          if (p_trCand_nomu[2].M() > 1e-10) {
            _h_sum3p_nomu->Fill(p_trCand_nomu[0].P() + p_trCand_nomu[1].P() + p_trCand_nomu[2].P(), weight);
            _h_mass3_nomu->Fill((p_trCand_nomu[0] + p_trCand_nomu[1] + p_trCand_nomu[2]).M(), weight);
            _h_R3_nomu->Fill((p_trCand_nomu[0].P() + p_trCand_nomu[1].P() + p_trCand_nomu[2].P()) / _sump, weight);
            _h_R3_Nch_nomu->Fill((p_trCand_nomu[0].P() + p_trCand_nomu[1].P() + p_trCand_nomu[2].P()) / _sump, (double)_Nch, weight);
          }
        }
      }

      //=================================
      // simple invariant combination
      //=================================

      int p1[6] = {0, 0, 1, 1, 2, 2};
      int p2[6] = {1, 2, 2, 0, 0, 1};

      TLorentzVector p_track1_D0combi, p_track2_D0combi, p_D0combi, p_D0optcombi;
      p_D0optcombi.SetPtEtaPhiM(0., 0., 0., 200.);
      int tk2charge = 0;

      if (fabs(pt_trCand_nomu[0]) > 1e-10 && fabs(pt_trCand_nomu[1]) > 1e-10 && fabs(pt_trCand_nomu[2]) > 1e-10) {
        for (unsigned int iD0combi = 0; iD0combi < 6; iD0combi++) {

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~
          // reconstruct D^0 to K Pi
          //~~~~~~~~~~~~~~~~~~~~~~~~~~
          int tk1 = p1[iD0combi];
          int tk2 = p2[iD0combi];

          // Opposite sign
          if (id_trCand_nomu[tk1]*id_trCand_nomu[tk2] > 0) continue;

          p_track1_D0combi.SetPtEtaPhiM(pt_trCand_nomu[tk1], eta_trCand_nomu[tk1], phi_trCand_nomu[tk1], gMassPi);
          p_track2_D0combi.SetPtEtaPhiM(pt_trCand_nomu[tk2], eta_trCand_nomu[tk2], phi_trCand_nomu[tk2], gMassK);
          p_D0combi = p_track1_D0combi + p_track2_D0combi;

          int iMaxMuInSelJet = -1;
          double maxMuInSelJet = -1;
          for (unsigned int iMuCand = 0; iMuCand < myPFmuInSelJet.size(); iMuCand++) {
            if (myPFmuInSelJet[iMuCand]->pdgId()*id_trCand_nomu[tk2] > 0) continue;
            if (myPFmuInSelJet[iMuCand]->pt() > maxMuInSelJet) {
              iMaxMuInSelJet = iMuCand;
              maxMuInSelJet = myPFmuInSelJet[iMaxMuInSelJet]->pt();
            }
          }

          if (fabs(p_D0combi.M() - gMassD0) < fabs(p_D0optcombi.M() - gMassD0)) {
            p_D0optcombi.SetPtEtaPhiM(p_D0combi.Pt(), p_D0combi.Eta(), p_D0combi.Phi(), p_D0combi.M());
            // tk2charge = id_trCand_nomu[tk2];
            tk2charge = iMaxMuInSelJet;
          }

          // cut on pT
          if (p_D0combi.Pt() < 15.) continue;

          _h_D0Mass->Fill(p_D0combi.M(), weight);
          _h_D0p->Fill(p_D0combi.P(), weight);
          _h_D0pT->Fill(p_D0combi.Pt(), weight);
          _h_D0eta->Fill(p_D0combi.Eta(), weight);

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~
          // associate D^0 to a PF muon
          //~~~~~~~~~~~~~~~~~~~~~~~~~~
          TLorentzVector p_Mu;

          if (iMaxMuInSelJet < 0) continue;

          p_Mu.SetPtEtaPhiM(myPFmuInSelJet[iMaxMuInSelJet]->pt(), myPFmuInSelJet[iMaxMuInSelJet]->eta(), myPFmuInSelJet[iMaxMuInSelJet]->phi(), gMassMu);
          // Just for MC : look at the genId of the selected PF mu
          double minDeltaRPFmuGen = 200.;
          double minDeltaPtPFmuGen = 200.;
          int idGen = 0;
          double pTGen = 0.;
          double dpz = 0.;
          double origin = 0.; // will b 1 if mu comes from t->W
          for (reco::GenParticleCollection::const_iterator it = mcparts.begin(); it != mcparts.end(); it++) {
            if (fabs((*it).pt()) < 1e-10) continue;
            TLorentzVector p_Gen;
            p_Gen.SetPtEtaPhiM((*it).pt(), (*it).eta(), (*it).phi(), 0.);
            if (p_Gen.DeltaR(p_Mu) < minDeltaRPFmuGen && fabs(p_Gen.Pt()-p_Mu.Pt()) < 0.15*p_Mu.Pt()) {
              minDeltaRPFmuGen = p_Gen.DeltaR(p_Mu);
              minDeltaPtPFmuGen = fabs(p_Gen.Pt()-p_Mu.Pt());
              idGen = abs((*it).pdgId());
              pTGen = p_Gen.Pt();
              dpz = fabs((p_Gen.Pz()-p_Mu.Pz())/p_Gen.Pz());
              origin = 0.;
              reco::GenParticle mutmp = *it;
              reco::GenParticle *mothertmp = NULL;
              for (int i=0; i<100; i++) {
                if (mutmp.mother(i) != 0 && fabs(mutmp.mother(i)->pdgId())==24) {
                  mothertmp = (reco::GenParticle*) mutmp.mother(i);
                  for (int i=0; i<100; i++) {
                    if (mothertmp->mother(i) != 0 && fabs(mothertmp->mother(i)->pdgId()) == 6) {
                      origin = 1.;
                      break;
                    }
                  }
                  break;
                } 
              }
            }
          }
          if (abs(idGen) > 0) {
            _h_id1stPFmu->Fill((double)idGen, weight);
            _h_id_pT_1stPFmu->Fill((double)idGen, pTGen, weight);
            _h_dR_dpT_1stPFmu->Fill(minDeltaRPFmuGen, minDeltaPtPFmuGen, weight);
            if (pTGen < 100) {
              _h_origin_id_1stPFmu_Inf100->Fill(origin, (double)idGen, weight);
              _h_eta_dpz_1stPFmu_Inf100->Fill(p_Mu.Eta(), dpz, weight);
            }
            else {
              _h_origin_id_1stPFmu_Sup100->Fill(origin, (double)idGen, weight);
              _h_eta_dpz_1stPFmu_Sup100->Fill(p_Mu.Eta(), dpz, weight);
            }
          }

          TLorentzVector p_Bcombi = p_Mu + p_D0combi;
          _h_BMomentum_unbiased->Fill(p_Bcombi.P(), weight);
        }
      }

      //=================================
      // optimized invariant combination
      //=================================

      if (fabs(p_D0optcombi.Pt()) > 1e-10) {

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // reconstruct D^0 to K Pi
        //~~~~~~~~~~~~~~~~~~~~~~~~~~

        // cut on pT
        if (p_D0optcombi.Pt() < 15.) continue;

        _h_D0MassClean->Fill(p_D0optcombi.M(), weight);
        _h_D0pClean->Fill(p_D0optcombi.P(), weight);
        _h_D0pTClean->Fill(p_D0optcombi.Pt(), weight);
        _h_D0etaClean->Fill(p_D0optcombi.Eta(), weight);

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // associate D^0 to a PF muon
        //~~~~~~~~~~~~~~~~~~~~~~~~~~

        TLorentzVector p_Mu;

        if (tk2charge < 0) continue;

        p_Mu.SetPtEtaPhiM(myPFmuInSelJet[tk2charge]->pt(), myPFmuInSelJet[tk2charge]->eta(), myPFmuInSelJet[tk2charge]->phi(), gMassMu);
        TLorentzVector p_Boptcombi = p_Mu + p_D0optcombi;
        _h_BMomentum->Fill(p_Boptcombi.P(), weight);
        _h_mup->Fill(p_Mu.P(), weight);
        if (p_D0optcombi.M() > 1.7 && p_D0optcombi.M() < 2.) 
          _h_BMomentumClean->Fill(p_Boptcombi.P(), weight);
      }
    } // jet loop
    _t_bjets->Fill();
  }
  }


  // ------------ method called once each job just before starting event loop  ------------
  void 
    MuTagForRivet::beginJob()
    {

      // std::cout << "Creating histos..." << std::endl;

      edm::Service<TFileService> fs;
      _h_idPFel = fs->make<TH1D>("pdgIdPFel", "pdgIdPFel", 500, 0., 500.);
      _h_id_pT_PFel = fs->make<TH2D>("pdgId-pT-PFel", "pdgId-pT-PFel", 500, 0., 500., 150, 0, 300);
      _h_dR_dpT_PFel = fs->make<TH2D>("DeltaR-DeltapT-PFel-Gen", "DeltaR-DeltapT-PFel-Gen", 100, 0., 0.5, 100, 0., 5.);
      _h_idPFmu = fs->make<TH1D>("pdgIdPFmu", "pdgIdPFmu", 500, 0., 500.);
      _h_id_pT_PFmu = fs->make<TH2D>("pdgId-pT-PFmu", "pdgId-pT-PFmu", 500, 0., 500., 150, 0, 300);
      _h_dR_dpT_PFmu = fs->make<TH2D>("DeltaR-DeltapT-PFmu-Gen", "DeltaR-DeltapT-PFmu-Gen", 100, 0., 0.5, 100, 0., 5.);
      _h_id1stPFmu = fs->make<TH1D>("pdgId1stPFmu", "pdgId1stPFmu", 500, 0., 500.);
      _h_id_pT_1stPFmu = fs->make<TH2D>("pdgId-pT-1stPFmu", "pdgId-pT-1stPFmu", 500, 0., 500., 150, 0, 300);
      _h_dR_dpT_1stPFmu = fs->make<TH2D>("DeltaR-DeltapT-1stPFmu-Gen", "DeltaR-DeltapT-1stPFmu-Gen", 100, 0., 0.5, 100, 0., 5.);
      _h_origin_id_1stPFmu_Inf100 = fs->make<TH2D>("Origin-pdgId-1stPFmu-Under100GeV", "Origin-pdgId-1stPFmu-Under100GeV", 2, 0., 2., 500, 0., 500.); 
      _h_eta_dpz_1stPFmu_Inf100 = fs->make<TH2D>("Eta-RelDeltaPz-1stPFmu-Under100GeV", "Eta-RelDeltaPz-1stPFmu-Under100GeV", 60, -3., 3., 100, 0., 1.); 
      _h_origin_id_1stPFmu_Sup100 = fs->make<TH2D>("Origin-pdgId-1stPFmu-Above100GeV", "Origin-pdgId-1stPFmu-Above100GeV", 2, 0., 2., 500, 0., 500.); 
      _h_eta_dpz_1stPFmu_Sup100 = fs->make<TH2D>("Eta-RelDeltaPz-1stPFmu-Above100GeV", "Eta-RelDeltaPz-1stPFmu-Above100GeV", 60, -3., 3., 100, 0., 1.); 
      _h_genId_pT_NonSelJets = fs->make<TH2D>("GenID-pT-below30GeVjets", "GenID-pT-below30GeVjets", 22, 0., 22., 30, 0., 30.);

      _h_nVtx = fs->make<TH1D>("NPrimaryVtx", "NPrimaryVtx", 50, 0., 50.); 

      _h_CSVSelJets = fs->make<TH1D>("CSV-b-jets", "CSV-b-jets", 100, 0., 1.);
      _h_genIdSelJets = fs->make<TH1D>("GenID-b-jets", "GenID-b-jets", 22, 0., 22.);
      _h_pTSelJets = fs->make<TH1D>("TransverseMomentum-b-jets", "TransverseMomentum-b-jets", 100, 0., 500.);
      _h_etaSelJets = fs->make<TH1D>("Eta-b-jets", "Eta-b-jets", 60, -3., 3.);
      _h_etach = fs->make<TH1D>("Etach-b-jets", "Etach-b-jets", 60, -3., 3.);
      _h_pTch = fs->make<TH1D>("TransverseMomentumch-b-jets", "TransverseMomentumch-b-jets", 100, 0., 100.);

      _h_Nch = fs->make<TH1D>("Nch-b-jets", "Nch-b-jets", 45, 0, 45);
      _h_sump = fs->make<TH1D>("Sump-b-jets", "Sump-b-jets", 200, 0, 1000);
      _h_sumpvec = fs->make<TH1D>("VectorialSump-b-jets", "VectorialSump-b-jets", 300, 0, 1500);
      _h_sum1p = fs->make<TH1D>("Highestp-b-jets", "Highestp-b-jets", 150, 0, 300);
      _h_sum2p = fs->make<TH1D>("Sum2p-b-jets", "Sum2p-b-jets", 150, 0, 300);
      _h_sum3p = fs->make<TH1D>("Sum3p-b-jets", "Sum3p-b-jets", 150, 0, 300);
      _h_mass3 = fs->make<TH1D>("Mass3-b-jets", "Mass3-b-jets", 400, 0., 10.);
      _h_R1 = fs->make<TH1D>("R1-b-jets", "R1-b-jets", 51, 0, 1.02);
      _h_R1_Nch = fs->make<TH2D>("R1-Nch-b-jets", "R1-Nch-b-jets", 51, 0, 1.02, 45, 0, 45);
      _h_R3 = fs->make<TH1D>("R3-b-jets", "R3-b-jets", 51, 0, 1.02);
      _h_R3_Nch = fs->make<TH2D>("R3-Nch-b-jets", "R3-Nch-b-jets", 51, 0, 1.02, 45, 0, 45);
      _h_sum1p_nomu = fs->make<TH1D>("Highestp-nomu-b-jets", "Highestp-nomu-b-jets", 150, 0, 300);
      _h_sum2p_nomu = fs->make<TH1D>("Sum2p-nomu-b-jets", "Sum2p-nomu-b-jets", 150, 0, 300);
      _h_sum3p_nomu = fs->make<TH1D>("Sum3p-nomu-b-jets", "Sum3p-nomu-b-jets", 150, 0, 300);
      _h_mass3_nomu = fs->make<TH1D>("Mass3-nomu-b-jets", "Mass3-nomu-b-jets", 400, 0., 10.);
      _h_R1_nomu = fs->make<TH1D>("R1-nomu-b-jets", "R1-nomu-b-jets", 51, 0, 1.02);
      _h_R1_Nch_nomu = fs->make<TH2D>("R1-Nch-nomu-b-jets", "R1-Nch-nomu-b-jets", 51, 0, 1.02, 45, 0, 45);
      _h_R3_nomu = fs->make<TH1D>("R3-nomu-b-jets", "R3-nomu-b-jets", 51, 0, 1.02);
      _h_R3_Nch_nomu = fs->make<TH2D>("R3-Nch-nomu-b-jets", "R3-Nch-nomu-b-jets", 51, 0, 1.02, 45, 0, 45);
      _h_D0Mass = fs->make<TH1D>("D0Mass-b-jets", "D0Mass-b-jets", 400, 0, 8);
      _h_D0p = fs->make<TH1D>("D0p-b-jets", "D0p-b-jets", 150, 0, 300);
      _h_D0pT = fs->make<TH1D>("D0pT-b-jets", "D0pT-b-jets", 100, 0, 400);
      _h_D0eta = fs->make<TH1D>("D0eta-b-jets", "D0eta-b-jets", 60, -3, 3);
      _h_BMomentum_unbiased = fs->make<TH1D>("BMomentum-nobias-b-jets", "BMomentum-nobias-b-jets", 100, 0, 400);
      _h_D0MassClean = fs->make<TH1D>("D0MassClean-b-jets", "D0MassClean-b-jets", 400, 0, 8);
      _h_D0pClean = fs->make<TH1D>("D0pClean-b-jets", "D0pClean-b-jets", 150, 0, 300);
      _h_D0pTClean = fs->make<TH1D>("D0pTClean-b-jets", "D0pTClean-b-jets", 100, 0, 400);
      _h_D0etaClean = fs->make<TH1D>("D0etaClean-b-jets", "D0etaClean-b-jets", 60, -3, 3);
      _h_BMomentum = fs->make<TH1D>("BMomentum-b-jets", "BMomentum-b-jets", 100, 0, 400);
      _h_mup = fs->make<TH1D>("Muonp-b-jets", "Muonp-b-jets", 150, 0, 300);
      _h_BMomentumClean = fs->make<TH1D>("BMomentum-D0cut-b-jets", "BMomentum-D0cut-b-jets", 100, 0, 400);

      _t_bjets = fs->make<TTree>("b-jets", "b-jets", 1);
      _t_bjets->Branch("Weight", &weight, "Weight/D");
      _t_bjets->Branch("CSV", &_CSV, "CSV/D");
      _t_bjets->Branch("vecP", _vecP, "vecP[4]/D");
      _t_bjets->Branch("Nch", &_Nch, "Nch/I");
      _t_bjets->Branch("Sump", &_sump, "Sump/D");
      _t_bjets->Branch("Tr1", _tr1, "Tr1[4]/D");
      _t_bjets->Branch("Tr2", _tr2, "Tr2[4]/D");
      _t_bjets->Branch("Tr3", _tr3, "Tr3[4]/D");
    }

  // ------------ method called once each job just after ending the event loop  ------------
  void 
    MuTagForRivet::endJob() 
    {

      // std::cout << "Closing histos..." << std::endl;
    }

  // ------------ method called when starting to processes a run  ------------
  void 
    MuTagForRivet::beginRun(edm::Run const&, edm::EventSetup const&)
    {
    }

  // ------------ method called when ending the processing of a run  ------------
  void 
    MuTagForRivet::endRun(edm::Run const&, edm::EventSetup const&)
    {
    }

  // ------------ method called when starting to processes a luminosity block  ------------
  void 
    MuTagForRivet::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }

  // ------------ method called when ending the processing of a luminosity block  ------------
  void 
    MuTagForRivet::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void
    MuTagForRivet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      // The following says we do not know what parameters are allowed so do no validation
      // Please change this to state exactly what you do use, even if it is no parameters
      edm::ParameterSetDescription desc;
      desc.setUnknown();
      descriptions.addDefault(desc);
    }

  // define this as a plug-in
  DEFINE_FWK_MODULE(MuTagForRivet);

