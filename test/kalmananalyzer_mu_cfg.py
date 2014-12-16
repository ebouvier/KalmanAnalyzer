import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# To build Transient Tracks
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'START53_V19PR::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_100_1_7tp.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_101_2_HIX.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_103_3_wXC.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_104_1_9x0.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_105_2_zw2.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_108_1_5TX.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_109_2_Umz.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_10_2_7e2.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_110_1_Vps.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_112_1_6kI.root',
'/store/user/verdier/SingleMu/SingleMu_Run2012B-22Jan2013_02Oct13-v3/7ed5d64fb39097b01209acde5484d3b2/patTuple_113_1_MeL.root',
    )
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('kalmanAnalyzed_Mu.root')
)

process.ana = cms.EDAnalyzer('KalmanAnalyzer_Mu'
)


process.p = cms.Path(process.ana)
