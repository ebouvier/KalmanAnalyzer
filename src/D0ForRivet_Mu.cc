// -*- C++ -*-
//
// Package:    D0ForRivet_Mu
// Class:      D0ForRivet_Mu
// 
/**\class D0ForRivet_Mu D0ForRivet_Mu.cc UserCode/KalmanAnalyzer/src/D0ForRivet_Mu.cc

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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

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

class D0ForRivet_Mu : public edm::EDAnalyzer {
  public:
    explicit D0ForRivet_Mu(const edm::ParameterSet&);
    ~D0ForRivet_Mu();

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
    bool hasGoodMuons;
    bool hasGoodElectrons;
    bool hasGoodLeptons;
    bool hasGoodJets;
    bool isGoodEvt;
    bool passTrigger;

    int _Nch[2];
    TH1D* _h_Nch[2];
    double _sump[2];
    TH1D* _h_sump[2];
    TLorentzVector _sumpvec[2];
    TH1D* _h_sumpvec[2];
    TH1D* _h_sum1p[2];
    TH1D* _h_sum3p[2];
    TH1D* _h_R3[2];
    TH1D* _h_D0Mass[2];
    TH1D* _h_D0MassClean[2];
    TH1D* _h_BMomentum[2];
    TH1D* _h_D0MassBlow[2];
    TH1D* _h_D0MassCleanBlow[2];

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
D0ForRivet_Mu::D0ForRivet_Mu(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  nEvts = 0;
  nEvts2 = 0;

  _Nch[0] = 0; _Nch[1] = 0;
  _sump[0] = 0.; _sump[1] = 0.;
  _sumpvec[0].SetPtEtaPhiM(0.,0.,0.,0.);

}


D0ForRivet_Mu::~D0ForRivet_Mu()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "Initial number of events                    = " << nEvts << std::endl;
  std::cout << "Initial number of semileptonic ttbar events = " << nEvts2 << std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
D0ForRivet_Mu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  nEvts++;
  _Nch[0] = 0; _Nch[1] = 0;
  _sump[0] = 0.; _sump[1] = 0.;
  _sumpvec[0].SetPtEtaPhiM(0.,0.,0.,0.);
  hasGoodMuons = false;
  hasGoodElectrons = false;
  hasGoodLeptons = false;
  hasGoodJets = false;
  isGoodEvt = false;
  passTrigger = false;

  //-------------------------------------------
  // Apply the recommended ttbar semilept sel
  //-------------------------------------------
  
  // Muon selection
  edm::Handle<std::vector<pat::Muon>>  muonHandle;
  edm::InputTag tagMuon("selectedPatMuonsPFlow");
  iEvent.getByLabel(tagMuon, muonHandle);
  const std::vector<pat::Muon>& muon = *(muonHandle.product());
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
  EAtarget  =ElectronEffectiveArea::kEleEAData2012; // for DATA

  for (unsigned int iElectron = 0; iElectron < electron.size(); ++iElectron) {
    if (electron[iElectron].pt() < 30.) continue;
    if (fabs(electron[iElectron].eta()) > 2.5) continue;
    if (!EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, electron[iElectron], hConversions, beamSpot, hVtx, electron[iElectron].chargedHadronIso(), electron[iElectron].photonIso(), electron[iElectron].neutralHadronIso(), rhoIso, EAtarget)) continue;
    if (fabs(electron[iElectron].superCluster()->eta()) >= 1.4442 && fabs(electron[iElectron].superCluster()->eta()) < 1.5660) continue;
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
  if (hasGoodMuons) hasGoodLeptons =true;

  // Jets selection
  edm::Handle<std::vector<pat::Jet>> jetHandle;
  edm::InputTag tagJet("selectedPatJetsPFlow","","PAT");
  iEvent.getByLabel(tagJet, jetHandle);
  const std::vector<pat::Jet>& jet = *(jetHandle.product());
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
    ++n20jet;
  }

  if (n55jet > 0 && n45jet > 1 && n35jet > 2 && n20jet > 3) hasGoodJets = true;

  if (hasGoodLeptons && hasGoodJets) isGoodEvt = true;

  // Trigger selection
  edm::Handle<edm::TriggerResults> triggerResults ;
  edm::InputTag tag("TriggerResults", "", "HLT");
  iEvent.getByLabel(tag,triggerResults);

  if (triggerResults.isValid()) {
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);

    for(int i = 0 ; i < static_cast<int>(triggerResults->size()); ++i) {
      if (triggerResults->accept(i)!=0) {
        std::string triggerName = triggerNames.triggerName(i);
        if (triggerName == "HLTriggerFinalPath") continue; // This one is pretty useless...
        if (triggerName[0] == 'A') continue;     // Remove AlCa HLT paths

        TString ThistrigName= triggerName.c_str();

        int run = static_cast<int>((iEvent.eventAuxiliary()).run());

        if ( run <  193834                 && ThistrigName.Contains("HLT_IsoMu17_eta2p1_TriCentralPFJet30")          )  passTrigger = true;
        if ( run >= 193834 && run < 194270 && ThistrigName.Contains("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30")      )  passTrigger = true;
        if ( run >= 194270 && run < 199698 && ThistrigName.Contains("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20"))  passTrigger = true;
        if ( run >= 199698                 && ThistrigName.Contains("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25"))  passTrigger = true;
      }
    }
  }
  
  if (isGoodEvt && passTrigger) {
    nEvts2++;
  
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
    //float        gSigmaK = 0.000001;

    ParticleMass gMassPi  = 0.13957018;
    //float        gSigmaPi = 0.00000001;

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
      //double btag = (*it).bDiscriminator("combinedSecondaryVertexBJetTags");  // save 2 jets of highest CSV among those with pT > 30 GeV/c
      double btag = (*it).pt(); //save 2 jets of highest pT among those with pT > 30 GeV/c
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

    iSelJet = -1;

    for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
      iSelJet++;

      if (iSelJet == maxind || iSelJet == maxind2){
        int indJet = -1;
        if (iSelJet == maxind) indJet = 0;
        if (iSelJet == maxind2) indJet = 1;

        TLorentzVector p_Jet;
        p_Jet.SetPtEtaPhiM((*it).pt(), (*it).eta(), (*it).phi(), (*it).mass());

        double pt_trCand_D0combi[3]  = {0., 0., 0.};
        double eta_trCand_D0combi[3] = {0., 0., 0.};
        double phi_trCand_D0combi[3] = {0., 0., 0.};
        double id_trCand_D0combi[3] = {0., 0., 0.};

        TLorentzVector p_trCand[3];

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

          TLorentzVector p_tr1;
          p_tr1.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassPi);
          if (fabs(pt_trCand_D0combi[0]) > 1e-10) {
            p_trCand[0].SetPtEtaPhiM(pt_trCand_D0combi[0], eta_trCand_D0combi[0], phi_trCand_D0combi[0], gMassPi);
            if (fabs(pt_trCand_D0combi[1]) > 1e-10 && fabs(pt_trCand_D0combi[2]) > 1e-10) {
              p_trCand[1].SetPtEtaPhiM(pt_trCand_D0combi[1], eta_trCand_D0combi[1], phi_trCand_D0combi[1], gMassPi);
              p_trCand[2].SetPtEtaPhiM(pt_trCand_D0combi[2], eta_trCand_D0combi[2], phi_trCand_D0combi[2], gMassPi);
            }
            else {
              p_trCand[1].SetPtEtaPhiM(0., 0., 0., 0.);
              p_trCand[2].SetPtEtaPhiM(0., 0., 0., 0.);
            }
          }
          else p_trCand[0].SetPtEtaPhiM(0., 0., 0., 0.);

          _Nch[indJet]++;
          _sump[indJet] = _sump[indJet] + p_tr1.P();
          _sumpvec[indJet] = _sumpvec[indJet] + p_tr1;
          
        } // 1st jet's track loop

        _h_Nch[indJet]->Fill((double)_Nch[indJet]);
        _h_sump[indJet]->Fill(_sump[indJet]);
        _h_sumpvec[indJet]->Fill(_sumpvec[indJet].P());
        if (p_trCand[0].M() > 1e-10) {
          _h_sum1p[indJet]->Fill(p_trCand[0].P());
          if (p_trCand[1].M() > 1e-10 && p_trCand[2].M() > 1e-10) {
            _h_sum3p[indJet]->Fill(p_trCand[0].P() + p_trCand[1].P() + p_trCand[2].P());
            _h_R3[indJet]->Fill((p_trCand[0].P() + p_trCand[1].P() + p_trCand[2].P()) / _sump[indJet]);
          }
        }

        //=================================
        // simple invariant combination
        //=================================

        int p1[6] = {0, 0, 1, 1, 2, 2};
        int p2[6] = {1, 2, 2, 0, 0, 1};

        TLorentzVector p_track1_D0combi, p_track2_D0combi, p_D0combi, p_D0optcombi;
        p_D0optcombi.SetPtEtaPhiM(0., 0., 0., 200.);
        //int tk2charge = 0;

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

            if (fabs(p_D0combi.M() - gMassD0) < fabs(p_D0optcombi.M() - gMassD0)) {
              p_D0optcombi.SetPtEtaPhiM(p_D0combi.Pt(), p_D0combi.Eta(), p_D0combi.Phi(), p_D0combi.M());
              //tk2charge = id_trCand_D0combi[tk2];
            }


            // cut on pT
            if (p_D0combi.Pt() < 15.) continue;

            _h_D0Mass[indJet]->Fill(p_D0combi.M());
            _h_D0MassBlow[indJet]->Fill(p_D0combi.M());

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

            // keep going if closest muon is close enough
            if (deltaRD0combiMu > 0.4) continue;

            p_Mu.SetPtEtaPhiM(myPFmu[iMu]->pt(), myPFmu[iMu]->eta(), myPFmu[iMu]->phi(), gMassMu);
            TLorentzVector p_Bcombi = p_Mu + p_D0combi;
            _h_BMomentum[indJet]->Fill(p_Bcombi.P());            
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

          _h_D0MassClean[indJet]->Fill(p_D0optcombi.M());
          _h_D0MassCleanBlow[indJet]->Fill(p_D0optcombi.M());

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~
          // associate D^0 to a PF muon
          //~~~~~~~~~~~~~~~~~~~~~~~~~~
          /*
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

          // keep going if closest muon is close enough
          if (deltaRD0optcombiMu > 0.4) continue;

          p_Mu.SetPtEtaPhiM(myPFmu[iMu]->pt(), myPFmu[iMu]->eta(), myPFmu[iMu]->phi(), gMassMu);
          TLorentzVector p_Boptcombi = p_Mu + p_D0optcombi;
          */
        }
      }
//      iSelJet++;
    } // jet loop
  }
}


// ------------ method called once each job just before starting event loop  ------------
  void 
D0ForRivet_Mu::beginJob()
{

  //  std::cout << "Creating histos..." << std::endl;

  edm::Service<TFileService> fs;

  _h_Nch[0] = fs->make<TH1D>("Nch-b-jet1", "Nch-b-jet1", 45, 0, 45);
  _h_Nch[1] = fs->make<TH1D>("Nch-b-jet2", "Nch-b-jet2", 45, 0, 45);
  _h_sump[0] = fs->make<TH1D>("Sump-b-jet1", "Sump-b-jet1", 200, 0, 1000);
  _h_sump[1] = fs->make<TH1D>("Sump-b-jet2", "Sump-b-jet2", 200, 0, 1000);
  _h_sumpvec[0] = fs->make<TH1D>("VectorialSump-b-jet1", "VectorialSump-b-jet1", 300, 0, 1500);
  _h_sumpvec[1] = fs->make<TH1D>("VectorialSump-b-jet2", "VectorialSump-b-jet2", 300, 0, 1500);
  _h_sum1p[0] = fs->make<TH1D>("Highestp-b-jet1", "Highestp-b-jet1", 150, 0, 300);
  _h_sum1p[1] = fs->make<TH1D>("Highestp-b-jet2", "Highestp-b-jet2", 150, 0, 300);
  _h_sum3p[0] = fs->make<TH1D>("Sum3p-b-jet1", "Sum3p-b-jet1", 150, 0, 300);
  _h_sum3p[1] = fs->make<TH1D>("Sum3p-b-jet2", "Sum3p-b-jet2", 150, 0, 300);
  _h_R3[0] = fs->make<TH1D>("R3-b-jet1", "R3-b-jet1", 50, 0, 1);
  _h_R3[1] = fs->make<TH1D>("R3-b-jet2", "R3-b-jet2", 50, 0, 1);
  _h_D0Mass[0] = fs->make<TH1D>("D0Mass-b-jet1", "D0Mass-b-jet1", 160, 0, 8);
  _h_D0Mass[1] = fs->make<TH1D>("D0Mass-b-jet2", "D0Mass-b-jet2", 160, 0, 8);
  _h_D0MassClean[0] = fs->make<TH1D>("D0MassClean-b-jet1", "D0MassClean-b-jet1", 160, 0, 8);
  _h_D0MassClean[1] = fs->make<TH1D>("D0MassClean-b-jet2", "D0MassClean-b-jet2", 160, 0, 8);
  _h_BMomentum[0] = fs->make<TH1D>("BMomentum-b-jet1", "BMomentum-b-jet1", 100, 0, 400);
  _h_BMomentum[1] = fs->make<TH1D>("BMomentum-b-jet2", "BMomentum-b-jet2", 100, 0, 400);
  _h_D0MassBlow[0] = fs->make<TH1D>("D0Mass-b-jet1-blow", "D0Mass-b-jet1-blow", 30, 1.7, 2.);
  _h_D0MassBlow[1] = fs->make<TH1D>("D0Mass-b-jet2-blow", "D0Mass-b-jet2-blow", 30, 1.7, 2.);
  _h_D0MassCleanBlow[0] = fs->make<TH1D>("D0MassClean-b-jet1-blow", "D0MassClean-b-jet1-blow", 30, 1.7, 2.);
  _h_D0MassCleanBlow[1] = fs->make<TH1D>("D0MassClean-b-jet2-blow", "D0MassClean-b-jet2-blow", 30, 1.7, 2.);

}

// ------------ method called once each job just after ending the event loop  ------------
  void 
D0ForRivet_Mu::endJob() 
{

  //  std::cout << "Closing histos..." << std::endl;
  /* Uncomment for Integral() = 1.
  _h_Nch[0]->Scale(1./_h_Nch[0]->Integral());
  _h_Nch[1]->Scale(1./_h_Nch[1]->Integral());
  _h_sump[0]->Scale(1./_h_sump[0]->Integral());
  _h_sump[1]->Scale(1./_h_sump[1]->Integral());
  _h_sumpvec[0]->Scale(1./_h_sumpvec[0]->Integral());
  _h_sumpvec[1]->Scale(1./_h_sumpvec[1]->Integral());
  _h_sum1p[0]->Scale(1./_h_sum1p[0]->Integral());
  _h_sum1p[1]->Scale(1./_h_sum1p[1]->Integral());
  _h_sum3p[0]->Scale(1./_h_sum3p[0]->Integral());
  _h_sum3p[1]->Scale(1./_h_sum3p[1]->Integral());
  _h_R3[0]->Scale(1./_h_R3[0]->Integral());
  _h_R3[1]->Scale(1./_h_R3[1]->Integral());
  _h_D0Mass[0]->Scale(1./_h_D0Mass[0]->Integral());
  _h_D0Mass[1]->Scale(1./_h_D0Mass[1]->Integral());
  _h_D0MassClean[0]->Scale(1./_h_D0MassClean[0]->Integral());
  _h_D0MassClean[1]->Scale(1./_h_D0MassClean[1]->Integral());
  _h_BMomentum[0]->Scale(1./_h_BMomentum[0]->Integral());
  _h_BMomentum[1]->Scale(1./_h_BMomentum[1]->Integral());
  _h_D0MassBlow[0]->Scale(1./_h_D0MassBlow[0]->Integral());
  _h_D0MassBlow[1]->Scale(1./_h_D0MassBlow[1]->Integral());
  _h_D0MassCleanBlow[0]->Scale(1./_h_D0MassCleanBlow[0]->Integral());
  _h_D0MassCleanBlow[1]->Scale(1./_h_D0MassCleanBlow[1]->Integral());
  */
}

// ------------ method called when starting to processes a run  ------------
  void 
D0ForRivet_Mu::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
  void 
D0ForRivet_Mu::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
  void 
D0ForRivet_Mu::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
  void 
D0ForRivet_Mu::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
D0ForRivet_Mu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(D0ForRivet_Mu);
