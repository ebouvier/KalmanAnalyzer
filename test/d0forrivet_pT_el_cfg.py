import FWCore.ParameterSet.Config as cms

process = cms.Process("pTbasedSelection")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# To build Transient Tracks
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'START53_V19PR::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_100_1_sBt.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_101_1_ocf.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_102_1_hgt.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_103_1_0RL.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_104_1_OrG.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_105_1_nkj.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_106_1_ZXF.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_107_1_hRN.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_108_1_fOU.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_109_1_RYV.root',
'/store/user/verdier/SingleElectron/SingleElectron_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_10_1_ucY.root',
    )
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('D0ForRivet_El.root')
)

process.pTbasedSelection = cms.EDAnalyzer('D0ForRivet_El',
        isCSVbased = cms.untracked.bool(False)
)


process.p = cms.Path(process.pTbasedSelection)
