// -*- C++ -*-
//
// Package:    KalmanAnalyzer
// Class:      KalmanAnalyzer
// 
/**\class KalmanAnalyzer KalmanAnalyzer.cc UserCode/KalmanAnalyzer/src/KalmanAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
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

  TH1D* HCHI2;

  TH1D* HMASS1;
  TH1D* HMASS2;
  TH1D* HMASS3;
  TH1D* HMASS4;

  TH1D* CAND_L;
  TH1D* CAND_L2;
  TH1D* CAND_SIGMAL;
  TH1D* CAND_LOVERSIG;

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
  event.getByLabel(tagPF, pfHandle);
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
  
  std::cout << " Number of PF muons = " << myPFparts.size() << std::endl;

  //-------------------------------------------
  // Access the jets 
  //-------------------------------------------
  
  //  double d0mass = 1.86484;

  ParticleMass kaon_mass  = 0.493677;
  float        kaon_sigma = 0.000001;

  ParticleMass pion_mass  = 0.13957018;
  float        pion_sigma = 0.00000001;

  // Track setup

  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);  

  edm::InputTag jetTag("selectedPatJetsPFlow","","PAT");

  edm::Handle<pat::JetCollection> jetHandle;
  iEvent.getByLabel(jetTag, jetHandle);
  pat::JetCollection jets = *jetHandle;
  
  //  std::cout << "Number of jets = " << jets.size() << std::endl;
  
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    reco::TrackRefVector jetTracks = (*it).associatedTracks();
    //    std::cout << " ->Number of Tracks = " << jetTracks.size() << std::endl;

    double btag = (*it).bDiscriminator("combinedSecondaryVertexBJetTags");

    // Btag medium requirement (medium = 0.679, tight = 0.898
    if ( (*it).pt() < 40.   ) continue;
    if ( btag       < 0.679 ) continue;

    int ngood  = 0;
    int ngood2 = 0;

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
      
	reco::TransientTrack tr2 = (*theB).build((**iter2));
	
	// Select OS tracks
	if ( (**iter1).charge()*(**iter2).charge() > 0 ) continue;

	// Compute the mass
	double px = (**iter1).px() + (**iter2).px();
	double py = (**iter1).py() + (**iter2).py();
	double pz = (**iter1).pz() + (**iter2).pz();
	double p  = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
	double e1 = sqrt(pow((**iter1).p(),2) + pow(kaon_mass,2));
	double e2 = sqrt(pow((**iter2).p(),2) + pow(pion_mass,2));
	double e  = e1 + e2;
	double m  = pow(e,2)-pow(p,2);
	if ( m>0 ) m = sqrt(m);
	else       m = 0.;

	// To select Jpsi
	//	if ( m < 3. || m > 3.4 ) continue;

	//Creating a KinematicParticleFactory
	KinematicParticleFactoryFromTransientTrack pFactory;
	
	//initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
	float chi = 0.;
	float ndf = 0.;
	
	//making particles
	std::vector<RefCountedKinematicParticle> d0Particles;
	d0Particles.push_back(pFactory.particle (tr1,kaon_mass,chi,ndf,kaon_sigma));
	d0Particles.push_back(pFactory.particle (tr2,pion_mass,chi,ndf,pion_sigma));
	
	/* Example of a simple vertex fit, without other constraints
	 * The reconstructed decay tree is a result of the kinematic fit
	 * The KinematicParticleVertexFitter fits the final state particles to their vertex and
	 * reconstructs the decayed state
	 */
	
	// creating the vertex fitter
	KinematicParticleVertexFitter fitter;
	// reconstructing a J/Psi decay
	RefCountedKinematicTree vertexFitTree = fitter.fit(d0Particles);
	
	if (!vertexFitTree->isValid()) continue;
	//	std::cout << "Fit is valid !!" << std::endl;

	//accessing the tree components, move pointer to top
	vertexFitTree->movePointerToTheTop();
	
	//We are now at the top of the decay tree getting the d0 reconstructed KinematicPartlcle
	RefCountedKinematicParticle d0 = vertexFitTree->currentParticle();
	//	  AlgebraicVector7 par0 = d0->currentState().kinematicParameters().vector();
	
	RefCountedKinematicVertex d0_vertex = vertexFitTree->currentDecayVertex();

	if ( !d0_vertex->vertexIsValid()) continue;
	//	std::cout << "Vertex fit is valid: " 
	//		  << d0_vertex->chiSquared() << " / " 
	//		  << d0_vertex->degreesOfFreedom() << std::endl;

	++ngood;
	if (d0_vertex->chiSquared()/(double)d0_vertex->degreesOfFreedom() < 1. ) ++ngood2;

	HCHI2->Fill(d0_vertex->chiSquared()/(double)d0_vertex->degreesOfFreedom());
	HMASS1->Fill(m);
	
	if (d0_vertex->chiSquared()/(double)d0_vertex->degreesOfFreedom() < 4. )
	  HMASS2->Fill(m);
	if (d0_vertex->chiSquared()/(double)d0_vertex->degreesOfFreedom() < 2. )
	  HMASS3->Fill(m);
	
	if ( d0_vertex->chiSquared()/(double)d0_vertex->degreesOfFreedom() > 2. ) continue;

	// Distance to PV :
	
	GlobalPoint svPos    = d0_vertex->position();
        GlobalError svPosErr = d0_vertex->error();

	double sigmax = sqrt( vtx[0].xError()*vtx[0].xError() + svPosErr.cxx()*svPosErr.cxx() );
        double sigmay = sqrt( vtx[0].yError()*vtx[0].yError() + svPosErr.cyy()*svPosErr.cyy() );
        double sigmaz = sqrt( vtx[0].zError()*vtx[0].zError() + svPosErr.czz()*svPosErr.czz() );

	double interx = pow( (px/m)/sigmax, 2.);
        double intery = pow( (py/m)/sigmay, 2.);
        double interz = pow( (pz/m)/sigmaz, 2.);
	
	double sigmaL3D = pow( interx + intery + interz , -0.5);

	double part1 = (px/m)*pow(sigmaL3D/sigmax,2.)*( svPos.x() - vtx[0].x());
        double part2 = (py/m)*pow(sigmaL3D/sigmay,2.)*( svPos.y() - vtx[0].y());
        double part3 = (pz/m)*pow(sigmaL3D/sigmaz,2.)*( svPos.z() - vtx[0].z());

	double L3D = fabs(part1 + part2 + part3);

	double L3DoverSigmaL3D = L3D/sigmaL3D;

	if ( L3DoverSigmaL3D < 50. ) continue;
	
	CAND_L->Fill(L3D);
	CAND_L2->Fill(L3D);
	CAND_SIGMAL->Fill(sigmaL3D);
	CAND_LOVERSIG->Fill(L3DoverSigmaL3D);
	
	HMASS4->Fill(m);

	
	//        std::cout << "L3D         = " << L3D      << std::endl;
	//        std::cout << "sigmaL3D    = " << sigmaL3D << std::endl;
	//        std::cout << "(L/sigma)3D = " << L3DoverSigmaL3D << std::endl;    


      } // 2nd jet's track loop

      //      std::cout << (**iter).pt() << " , " << (**iter).charge() << std::endl;
      
      
    } // 1st jet's track loop
    
    //    std::cout << "Number of good fits = " << ngood << std::endl;
    //    std::cout << "Number of good fits with chi2 cut = " << ngood2 << std::endl;

  } // jet loop

}


// ------------ method called once each job just before starting event loop  ------------
void 
KalmanAnalyzer::beginJob()
{

  //  std::cout << "Creating histos..." << std::endl;

  edm::Service<TFileService> fs;
  HCHI2  = fs->make<TH1D>("HCHI2","HCHI2",110,0.,11.);
  HMASS1 = fs->make<TH1D>("HMASS1","HMASS1",1000,0.,10.);
  HMASS2 = fs->make<TH1D>("HMASS2","HMASS2",1000,0.,10.);
  HMASS3 = fs->make<TH1D>("HMASS3","HMASS3",1000,0.,10.);
  HMASS4 = fs->make<TH1D>("HMASS4","HMASS4",1000,0.,10.);
  
  CAND_L        = fs->make<TH1D>("CAND_L","CAND_L",1000,0.,1.);
  CAND_L2       = fs->make<TH1D>("CAND_L2","CAND_L2",1000,0.,10.);
  CAND_SIGMAL   = fs->make<TH1D>("CAND_SIGMAL","CAND_SIGMAL",5000,0.,0.005);
  CAND_LOVERSIG = fs->make<TH1D>("CAND_LOVERSIG","CAND_LOVERSIG",21000,0.,7000.);

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
