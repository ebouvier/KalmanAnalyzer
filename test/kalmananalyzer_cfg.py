import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# To build Transient Tracks
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'START53_V19PR::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/gridgroup/cms/verdier/patTuple_871_1_WR1.root'
        
        # TTbar semi-lept madgraph with Jpsi filter:
        #'/store/user/verdier/TTbar173_16Feb13-v1/pv-patprod-TTbar173_18Mar13-v1/92b2e837d814b85b002ea027156a7dca/patTuple_1_1_SdN.root',
        #'/store/user/verdier/TTbar173_16Feb13-v1/pv-patprod-TTbar173_18Mar13-v1/92b2e837d814b85b002ea027156a7dca/patTuple_2_1_1mY.root',
        #'/store/user/verdier/TTbar173_16Feb13-v1/pv-patprod-TTbar173_18Mar13-v1/92b2e837d814b85b002ea027156a7dca/patTuple_3_1_TQW.root',
        #'/store/user/verdier/TTbar173_16Feb13-v1/pv-patprod-TTbar173_18Mar13-v1/92b2e837d814b85b002ea027156a7dca/patTuple_4_1_l2J.root',
        #'/store/user/verdier/TTbar173_16Feb13-v1/pv-patprod-TTbar173_18Mar13-v1/92b2e837d814b85b002ea027156a7dca/patTuple_5_1_63z.root'

        # TTbar semi-lept madgraph without Jpsi filter:
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1000_1_lIh.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1001_1_61N.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1002_1_jif.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1003_1_YSz.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1004_1_3xc.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1005_1_ZeB.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1006_1_VLU.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1007_1_3Wn.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1008_1_XOo.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1009_1_LsP.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_100_1_GBz.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1010_1_H41.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1011_1_yzA.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1012_1_Btf.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1013_1_o0W.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1014_1_fSM.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1015_1_C4Y.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1016_1_OfX.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1017_1_IPC.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1018_1_PGR.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1019_1_zxE.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_101_1_fL0.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1020_1_jLx.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1021_1_D7w.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1022_1_ezF.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1023_1_GnN.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1024_1_wy5.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1025_1_9Us.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1026_1_tH5.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1027_1_V64.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1028_1_Gia.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1029_1_vTe.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_102_1_IXB.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1030_1_Om0.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1031_1_Yrt.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1032_1_zks.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1033_1_AFy.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1034_1_fzd.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1035_1_DiN.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1036_1_qPz.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1037_1_BTQ.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1038_1_icd.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1039_1_aZZ.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_103_1_qnl.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1040_1_tS2.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1041_1_lEr.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1042_1_Zg5.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1043_1_zR8.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1044_1_Zfi.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1045_1_QPa.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1046_1_sCn.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1047_1_I6k.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1048_1_NtG.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1049_1_kPo.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_104_1_iXu.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1050_1_M7x.root'

    )
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('kalmanAnalyzed.root')
)

process.ana = cms.EDAnalyzer('KalmanAnalyzer'
)


process.p = cms.Path(process.ana)
