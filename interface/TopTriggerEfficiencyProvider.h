#ifndef TOPTRIGGERWEIGHTPROVIDER_H
#define TOPTRIGGERWEIGHTPROVIDER_H

/**                                                                                                                                                                        
 * \class TopTriggerWeightProvider                                                                                                                                   
 *                                                                                                                                                                         
 * \brief Abstract: Class to read the trigger efficiency weights and
 * apply them on MC. Note that this class does not check whether you
 * are running on Data or MC, this is the responsibility of the
 * analyst.  In the class the various datataking periods have been
 * taken into account. One should note however that also the cut
 * selection changes during datataking, and this class assumes the
 * analyst is running the correct cutset for each run period. The
 * provided weights also assume that the analyst is comparing with
 * data with the correct trigger bit, depending on the run period
 *                                                                                                                                                                         
 * \author Kelly Beernaert, Ghent University, v1 08/11/2012                                                                                                                
 */


class TopTriggerEfficiencyProvider {
 public:
  
  TopTriggerEfficiencyProvider(bool verbose = true, double *lumis = 0);
  ~TopTriggerEfficiencyProvider() {};
  
  enum Runs { RunBEGIN=0, RunA=0, RunB, RunC, RunD, RunEND};
  enum JES { NOMINAL=0, DOWN=-1, UP=1};

  void setLumi(Runs r, double lumi) {luminosity[r]=lumi;}

  double get_weight(double lep_pt, double lep_eta, 
		    double jet_pt, double jet_eta, 
		    int npvertices, int njets, 
		    bool LepIsMuon, JES jes=NOMINAL) const;
  
 protected:
  double GetLeptonWeight(bool isMuon, double pt, double eta, int npvertices) const;
  double GetJetWeight(double pt, double eta, int njets, JES jes) const;
  double TurnOn(double x, double par[4]) const;
  double VFunction(int x, const double par[2]) const;
  void warn(const char*) const;

  double luminosity[RunEND];
  const bool verbose;
  
  const static unsigned
    nMuEtabins, nMuPtbins,
    nElEtabins, nElPtbins,
    nJetEtabins,nJetPtbins,
    nJetMin,nJetNbins;

  const static double 
    lumiDefaults[RunEND],
    jetPtMin, muPtMin, elPtMin,
    jetPtMax, muPtMax, elPtMax,
    jetEtaMax,muEtaMax,elEtaMax,
    JetEtaEdges[],JetPtEdges[],
    ElEtaEdges[],  ElPtEdges[],
    MuEtaEdges[],  MuPtEdges[],
    Mu[RunEND][6][9][2],
    El[RunEND][11][8][2],
    Jet[RunEND][4][5][11],
    JetDown[RunEND][4][5][11],
    JetUp[RunEND][4][5][11];
};

#endif
