import FWCore.ParameterSet.Config as cms

process = cms.Process("muTagBasedSelection")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# To build Transient Tracks
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'START53_V19PR::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

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
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1050_1_M7x.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1051_1_cXy.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1052_1_r6P.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1053_1_OmC.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1054_1_0nl.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1055_1_bxe.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1056_1_suL.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1057_1_h2s.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1058_1_v9v.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1059_1_kAo.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_105_1_wRt.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1060_1_M8i.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1061_1_Pqf.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1062_1_lC1.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1063_1_8Aj.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1064_1_rV8.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1065_1_le0.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1066_1_Thq.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1067_1_D9G.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1068_1_Q9B.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1069_1_CmP.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_106_1_1oz.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1070_1_ooP.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1071_1_fat.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1072_1_yxO.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1073_1_rzz.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1074_1_V47.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1075_1_dv4.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1076_1_2as.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1077_1_0sP.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1078_1_Mdc.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1079_1_3WS.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_107_1_eQD.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1080_1_Zk3.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1081_1_2pJ.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1082_1_DWo.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1083_1_LTG.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1084_1_qgx.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1085_1_Cmb.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1086_1_EXA.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1087_1_glL.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1088_1_uZn.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1089_1_s7D.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_108_1_gOB.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1090_1_8qV.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1091_1_tVm.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1092_1_ZM9.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1093_1_x0P.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1094_1_dc6.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1095_1_zET.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1096_1_cnE.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1097_1_u2e.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1098_1_zAH.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1099_1_J5a.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_109_1_Mgg.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_10_1_YlN.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1100_1_Upw.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1101_1_9ZC.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1102_1_oWf.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1103_1_m4y.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1104_1_x9m.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1105_1_jDL.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1106_1_G0b.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1107_1_mum.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1108_1_N7D.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1109_1_UZU.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_110_1_mXt.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1110_1_Tmh.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1111_1_Dm3.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1112_1_eNl.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1113_1_Ryk.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1114_1_OtX.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1115_1_oHa.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1116_1_AJZ.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1117_1_AOi.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1118_1_oad.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1119_1_B4n.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_111_1_fKQ.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1120_1_NeD.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1121_1_T0g.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1122_1_ra1.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1123_1_diN.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1124_1_nNo.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1125_1_7jW.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1126_1_95a.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1127_1_JYH.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1128_1_zeJ.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1129_1_cim.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_112_2_6HY.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1130_1_YlY.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1131_1_LlC.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1132_1_VNd.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1133_1_2ky.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1134_1_O3C.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1135_1_k0c.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1136_1_SXU.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1137_1_i1B.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1138_1_p0R.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1139_1_qix.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_113_1_xoP.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1140_1_k4y.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1141_1_q3M.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1142_1_KXN.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1143_1_3Yb.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1144_1_77m.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1145_1_tfN.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1146_1_qdf.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1147_1_Ubr.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1148_1_fk3.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1149_1_Wcs.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_114_1_M8z.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1150_1_XEr.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1151_1_Fuh.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1152_1_U18.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1153_1_n6A.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1154_1_w6x.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1155_1_Zc8.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1156_1_RTt.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1157_1_lVu.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1158_1_8OZ.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1159_1_Nci.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_115_1_TkA.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1160_1_pwC.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1161_1_P8A.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1162_1_YKa.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1163_1_rvz.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1164_1_nem.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1165_1_JCS.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1166_1_hNp.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1167_1_bY8.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1168_1_f0F.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1169_1_3oc.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_116_1_Fuv.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1170_1_vqp.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1171_1_wRc.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1172_1_Dm8.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1173_1_9qj.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1174_1_NaS.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1175_1_Mr5.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1176_1_TcX.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1177_1_uGn.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1178_1_9Zs.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1179_1_rve.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_117_1_2iP.root',
'/store/user/verdier/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTbar-semilept_05Oct13-v3/80a352c723b543c465f0aa98c4bbef52/patTuple_1180_1_5us.root'
    )
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('MuTagForRivet.root')
)

process.muTagBasedSelection = cms.EDAnalyzer('MuTagForRivet',
        # selection = cms.untracked.int32(1) # Top PAG jets selection and pT(tr) > 4 GeV/c
        selection = cms.untracked.int32(2) # 4 jets with pT > 30 GeV/c and pT(tr) > 0.5 GeV/c
)


process.p = cms.Path(process.muTagBasedSelection)
