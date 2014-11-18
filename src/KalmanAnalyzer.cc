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

    //  TFile* tfile_;
    TH1D* h_nJets;
    TH1D* h_CSV;

    TH1D* h_B0_cuts;

    TH1D* h_D0Cand_Chi2NDOF;

    TH1D* h_D0Cand_MassChi2Inf1;
    TH1D* h_D0Cand_MassChi2Inf1p5;
    TH1D* h_D0Cand_MassChi2Inf2;
    TH1D* h_D0Cand_MassChi2Inf2p5;
    TH1D* h_D0Cand_MassChi2Inf3;
    TH1D* h_D0Cand_MassChi2Inf3p5;
    TH1D* h_D0Cand_MassChi2Inf4;

    TH1D* h_D0Cand_L;
    TH1D* h_D0Cand_SigmaL;
    TH1D* h_D0Cand_LOverSigmaL;

    TH1D* h_D0_L;
    TH1D* h_D0_SigmaL;
    TH1D* h_D0_Mass;

    TH1D* h_BCand_DeltaRD0Mu;
    TH1D* h_B_D0Mass;
    TH1D* h_B_Mass;

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

}


KalmanAnalyzer::~KalmanAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
KalmanAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  //-------------------------------------------
  // Access the vertex
  //-------------------------------------------

  edm::Handle<reco::VertexCollection>  vtxHandle;
  edm::InputTag tagVtx("offlinePrimaryVertices");
  iEvent.getByLabel(tagVtx, vtxHandle);
  const reco::VertexCollection vtx = *(vtxHandle.product());

  if ( vtx.size() == 0 ) {
    std::cout << " WARNING : no PV for this event ... " << std::endl;
    return;
  }

  //--------------------------------------------------
  // Access the PF candidates for non-isolated muons
  //--------------------------------------------------

  edm::Handle<reco::PFCandidateCollection>  pfHandle;
  edm::InputTag tagPF("particleFlow","","RECO");
  iEvent.getByLabel(tagPF, pfHandle);
  reco::PFCandidateCollection pfs = *pfHandle;
  if ( !pfHandle.isValid() ) {
    std::cout << "=> pfHandle is not valid..." << std::endl;
    return;
  }

  // Select good PF muons :

  std::vector<const reco::PFCandidate*> myPFparts;
  for ( unsigned int i = 0; i < pfs.size(); ++i ) {

    if ( abs(pfs[i].pdgId()) != 13  ) continue;
    if ( pfs[i].pt()         <   4. ) continue;

    myPFparts.push_back(&pfs[i]);
  }

  //-------------------------------------------
  // Access the jets 
  //-------------------------------------------

  //  double d0mass = 1.86484;
  //  double bmass = 5.2796;
  ParticleMass gMassMu  = 0.105658367;
  //float        gSigmaMu = 0.000000004; 

  ParticleMass gMassK  = 0.493677;
  float        gSigmaK = 0.000001;

  ParticleMass gMassPi  = 0.13957018;
  float        gSigmaPi = 0.00000001;

  // Track setup

  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);  

  edm::InputTag jetTag("selectedPatJetsPFlow","","PAT");

  edm::Handle<pat::JetCollection> jetHandle;
  iEvent.getByLabel(jetTag, jetHandle);
  pat::JetCollection jets = *jetHandle;

  float maxcsv = -1.;
  unsigned int maxind = -1;
  float second_max = -1.0;
  unsigned int maxind2 = -1;
  unsigned int iSelJet = 0;

  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    double btag = (*it).bDiscriminator("combinedSecondaryVertexBJetTags");
    if ( (*it).pt() >= 30. ) {
      if(btag >= maxcsv) {
        second_max = maxcsv;
        maxcsv = btag;
        maxind = iSelJet;
        maxind2 = maxind;
      }
      else if(btag > second_max) {
        second_max = btag;
        maxind2 = iSelJet;
      }
    }
    iSelJet++;
  }
  h_nJets->Fill(iSelJet);

  iSelJet = 0;

  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {

    if (iSelJet == maxind || iSelJet == maxind2){
      h_CSV->Fill((*it).bDiscriminator("combinedSecondaryVertexBJetTags"));

      reco::TrackRefVector jetTracks = (*it).associatedTracks();

      for (reco::track_iterator iter1 = jetTracks.begin(); iter1 != jetTracks.end(); ++iter1) {

        const reco::Track& Track1 = **iter1;

        if ( (**iter1).pt() < 4.                 ) continue;
        if ( !Track1.quality(reco::Track::highPurity)) continue;

        reco::TransientTrack tr1 = (*theB).build((**iter1));

        for (reco::track_iterator iter2 = jetTracks.begin(); iter2 != jetTracks.end(); ++iter2) {
          const reco::Track& Track2 = **iter2;

          if ( iter2 == iter1      ) continue;
          if ( (**iter2).pt() < 4. ) continue;
          if ( !Track2.quality(reco::Track::highPurity)) continue;

          int iB0Cut = 0;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"2 tracks with p_{T} > 4 GeV/c");

          reco::TransientTrack tr2 = (*theB).build((**iter2));

          //~~~~~~~~~~~~~~~~~~~~~~~~~~
          // reconstruct D^0 -> K Pi
          //~~~~~~~~~~~~~~~~~~~~~~~~~~

          // Select OS tracks
          if ( (**iter1).charge()*(**iter2).charge() > 0 ) continue;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... of opposite signs");

          // Compute the mass
          TLorentzVector p_tr1_D0, p_tr2_D0, p_D0;
          p_tr1_D0.SetPtEtaPhiM((**iter1).pt(), (**iter1).eta(), (**iter1).phi(), gMassK);
          p_tr2_D0.SetPtEtaPhiM((**iter2).pt(), (**iter2).eta(), (**iter2).phi(), gMassPi);

          if (p_tr1_D0.DeltaR(p_tr2_D0) > 0.2) continue;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... within #DeltaR < 0.2");

          p_D0 = p_tr1_D0 + p_tr2_D0;

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

          if (!D0vertexFitTree->isValid()) continue;

          //accessing the tree components, move pointer to top
          D0vertexFitTree->movePointerToTheTop();

          //We are now at the top of the decay tree getting the d0 reconstructed KinematicPartlcle
          RefCountedKinematicParticle D0 = D0vertexFitTree->currentParticle();
          RefCountedKinematicVertex D0_vertex = D0vertexFitTree->currentDecayVertex();

          if ( !D0_vertex->vertexIsValid()) continue;

          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... with a valid D^{0} vertex");

          h_D0Cand_Chi2NDOF->Fill(D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom());

          if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 1.) 
            h_D0Cand_MassChi2Inf1->Fill(p_D0.M());
          if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 1.5) 
            h_D0Cand_MassChi2Inf1p5->Fill(p_D0.M());
          if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 2.)
            h_D0Cand_MassChi2Inf2->Fill(p_D0.M());
          if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 2.5)
            h_D0Cand_MassChi2Inf2p5->Fill(p_D0.M());
          if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 3.)
            h_D0Cand_MassChi2Inf3->Fill(p_D0.M());
          if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 3.5)
            h_D0Cand_MassChi2Inf3p5->Fill(p_D0.M());
          if (D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() < 4.)
            h_D0Cand_MassChi2Inf4->Fill(p_D0.M());

          // cut on chi2/NDOF
          if ( D0_vertex->chiSquared()/(double)D0_vertex->degreesOfFreedom() > 3.5) continue;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... with #chi^{2}/NDOF < 3.5");

          // Distance to PV :
          GlobalPoint D0_svPos    = D0_vertex->position();
          GlobalError D0_svPosErr = D0_vertex->error();

          double sigmax_vtx_D0vtx = sqrt( vtx[0].xError()*vtx[0].xError() + D0_svPosErr.cxx()*D0_svPosErr.cxx() );
          double sigmay_vtx_D0vtx = sqrt( vtx[0].yError()*vtx[0].yError() + D0_svPosErr.cyy()*D0_svPosErr.cyy() );
          double sigmaz_vtx_D0vtx = sqrt( vtx[0].zError()*vtx[0].zError() + D0_svPosErr.czz()*D0_svPosErr.czz() );

          double D0_interx = pow( (p_D0.Px()/p_D0.M())/sigmax_vtx_D0vtx, 2.);
          double D0_intery = pow( (p_D0.Py()/p_D0.M())/sigmay_vtx_D0vtx, 2.);
          double D0_interz = pow( (p_D0.Pz()/p_D0.M())/sigmaz_vtx_D0vtx, 2.);

          double D0_sigmaL3D = pow( D0_interx + D0_intery + D0_interz , -0.5);

          double D0_part1 = (p_D0.Px()/p_D0.M())*pow(D0_sigmaL3D/sigmax_vtx_D0vtx,2.)*( D0_svPos.x() - vtx[0].x());
          double D0_part2 = (p_D0.Py()/p_D0.M())*pow(D0_sigmaL3D/sigmay_vtx_D0vtx,2.)*( D0_svPos.y() - vtx[0].y());
          double D0_part3 = (p_D0.Pz()/p_D0.M())*pow(D0_sigmaL3D/sigmaz_vtx_D0vtx,2.)*( D0_svPos.z() - vtx[0].z());

          double D0_L3D = fabs(D0_part1 + D0_part2 + D0_part3);

          double D0_L3DoverSigmaL3D = D0_L3D/D0_sigmaL3D;

          h_D0Cand_L->Fill(D0_L3D);
          h_D0Cand_SigmaL->Fill(D0_sigmaL3D);
          h_D0Cand_LOverSigmaL->Fill(D0_L3DoverSigmaL3D);

          // cut on L/SigmaL
          if ( D0_L3DoverSigmaL3D < 50. ) continue;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... with c#tau/#sigma(c#tau) > 50");

          h_D0_L->Fill(D0_L3D);
          h_D0_SigmaL->Fill(D0_sigmaL3D);
          h_D0_Mass->Fill(p_D0.M());

          // 3 sigma window around the peak signal
          if (p_D0.M() < 1.8089 || p_D0.M() > 1.9229) continue;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... 1.8089 < M(D^{0}) < 1.9229 GeV/c^{2}");

          //~~~~~~~~~~~~~~~~~~~~~~~~~~~
          // associate D^0 to a PF muon
          //~~~~~~~~~~~~~~~~~~~~~~~~~~
          TLorentzVector p_Mu;
          int iMu = -1;
          double deltaRD0Mu = 2000.;

          // find closest muon
          for (unsigned int iMuCand = 0; iMuCand < myPFparts.size(); iMuCand++) {
            TLorentzVector p_MuCand;
            p_MuCand.SetPtEtaPhiM(myPFparts[iMuCand]->pt(), myPFparts[iMuCand]->eta(), myPFparts[iMuCand]->phi(), gMassMu);
            double tmp_deltaRD0Mu = p_D0.DeltaR(p_MuCand);
            if (tmp_deltaRD0Mu < deltaRD0Mu) {
              deltaRD0Mu = tmp_deltaRD0Mu;
              iMu = iMuCand;
            }
          }

          if (iMu < 0) continue;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... at least a non isolated #mu");
          h_BCand_DeltaRD0Mu->Fill(deltaRD0Mu);

          // keep going if closest muon is close enough
          if (deltaRD0Mu > 0.4) continue;
          h_B0_cuts->Fill((double)iB0Cut); ++iB0Cut;
          h_B0_cuts->GetXaxis()->SetBinLabel(iB0Cut,"... within #DeltaR < 0.4");

          h_B_D0Mass->Fill(p_D0.M());
          
          p_Mu.SetPtEtaPhiM(myPFparts[iMu]->pt(), myPFparts[iMu]->eta(), myPFparts[iMu]->phi(), gMassMu);
          TLorentzVector p_B = p_Mu + p_D0;
          h_B_Mass->Fill(p_B.M());

          //reco::GsfTrackRef tr_Mu = myPFparts[iMu]->gsfTrackRef();
          //*tr_Mu is a "const reco::GsfTrack" 

        } // 2nd jet's track loop
      } // 1st jet's track loop
    }
    iSelJet++;
  } // jet loop

}


// ------------ method called once each job just before starting event loop  ------------
  void 
KalmanAnalyzer::beginJob()
{

  //  std::cout << "Creating histos..." << std::endl;

  edm::Service<TFileService> fs;

  h_nJets = fs->make<TH1D>("h_nJets","h_nJets", 20, 0., 20.);
  h_CSV = fs->make<TH1D>("h_CSV","h_CSV", 100, 0., 1.);

  h_B0_cuts = fs->make<TH1D>("h_B0_cuts","h_B0_cuts",30,0.,30.);
  h_B0_cuts->SetOption("bar");
  h_B0_cuts->SetBarWidth(0.75);
  h_B0_cuts->SetBarOffset(0.125);

  h_D0Cand_Chi2NDOF       = fs->make<TH1D>("h_D0Cand_Chi2NDOF","h_D0Cand_Chi2NDOF",110,0.,11.);
  h_D0Cand_MassChi2Inf1   = fs->make<TH1D>("h_D0Cand_MassChi2Inf1","h_D0Cand_MassChi2Inf1",1000,0.,10.);
  h_D0Cand_MassChi2Inf1p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf1p5","h_D0Cand_MassChi2Inf1p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf2   = fs->make<TH1D>("h_D0Cand_MassChi2Inf2","h_D0Cand_MassChi2Inf2",1000,0.,10.);
  h_D0Cand_MassChi2Inf2p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf2p5","h_D0Cand_MassChi2Inf2p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf3   = fs->make<TH1D>("h_D0Cand_MassChi2Inf3","h_D0Cand_MassChi2Inf3",1000,0.,10.);
  h_D0Cand_MassChi2Inf3p5 = fs->make<TH1D>("h_D0Cand_MassChi2Inf3p5","h_D0Cand_MassChi2Inf3p5",1000,0.,10.);
  h_D0Cand_MassChi2Inf4   = fs->make<TH1D>("h_D0Cand_MassChi2Inf4","h_D0Cand_MassChi2Inf4",1000,0.,10.);

  h_D0Cand_L           = fs->make<TH1D>("h_D0Cand_L","h_D0Cand_L",1000,0.,1.);
  h_D0Cand_SigmaL      = fs->make<TH1D>("h_D0Cand_SigmaL","h_D0Cand_SigmaL",5000,0.,0.005);
  h_D0Cand_LOverSigmaL = fs->make<TH1D>("h_D0Cand_LOverSigmaL","h_D0Cand_LOverSigmaL",21000,0.,7000.);

  h_D0_Mass    = fs->make<TH1D>("h_D0_Mass","h_D0_Mass",1000,0.,10.);
  h_D0_L       = fs->make<TH1D>("h_D0_L","h_D0_L",1000,0.,1.);
  h_D0_SigmaL  = fs->make<TH1D>("h_D0_SigmaL","h_D0_SigmaL",5000,0.,0.005);

  h_BCand_DeltaRD0Mu = fs->make<TH1D>("h_BCand_DeltaRD0Mu","h_BCand_DeltaRD0Mu",100,0.,5.);
  h_B_D0Mass    = fs->make<TH1D>("h_B_D0Mass","h_B_D0Mass",114,1.8089,1.9229);
  h_B_Mass    = fs->make<TH1D>("h_B_Mass","h_B_Mass",1000,0.,10.);
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
