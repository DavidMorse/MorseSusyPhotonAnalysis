import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')
# add a list of strings for events to process
options.register ('eventsToProcess',
                                  '',
                                  VarParsing.multiplicity.list,
                                  VarParsing.varType.string,
                                  "Events to process")
options.parseArguments()

process = cms.Process("PickEvent")
process.source = cms.Source ("PoolSource",
          fileNames = cms.untracked.vstring (
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/194/199/ACAFEE65-80A0-E111-A69F-002481E0D90C.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/776/3455DA19-82B3-E111-8965-0030486780B8.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/194/912/9098C039-33A8-E111-BC73-5404A640A639.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/930/0E6DF538-7DB4-E111-B3ED-001D09F25109.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/194/050/92F841B0-2E9E-E111-9FA4-5404A63886CE.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/656/749AF91D-ADB2-E111-AF87-003048D2BE12.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/194/050/7632C65D-2A9E-E111-9B95-BCAEC53296F7.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/378/F4BACCE4-E4AE-E111-BBFD-E0CB4E4408C4.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/113/58BA6DA3-BCAA-E111-B1E6-00237DDC5C24.root',
#    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/194/151/926D4DDC-0CA0-E111-9197-0030486780AC.root',
##    '/store/data/Run2012A/Photon/RECO/PromptReco-v1/000/190/895/020B0954-1485-E111-AE98-BCAEC518FF80.root',
##    '/store/data/Run2012A/Photon/RECO/PromptReco-v1/000/190/906/60AD0A95-BC85-E111-9E7B-003048F118D4.root',
##    '/store/data/Run2012A/Photon/RECO/PromptReco-v1/000/191/247/10D76EB4-4888-E111-9921-5404A63886C4.root',

    
    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/776/6AD71B1B-82B3-E111-90E7-003048D2C174.root',
    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/930/1C00583A-7DB4-E111-A6E8-BCAEC532971B.root',
    #'/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/194/050/76D4165E-2D9E-E111-8450-BCAEC518FF74.root',
    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/656/7C053F6A-8AB2-E111-B8E9-001D09F29146.root',
    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/378/F4F488B7-E4AE-E111-A135-BCAEC518FF68.root',
    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/195/113/56CCAC9F-BCAA-E111-B262-003048F024C2.root',
    '/store/data/Run2012B/DoublePhoton/RECO/PromptReco-v1/000/194/151/B46EEDE1-0CA0-E111-B005-001D09F276CF.root',
    '/store/data/Run2012A/Photon/RECO/23May2012-v2/0000/AA1EABA1-C5A5-E111-8574-0025B31E3CBE.root',
    '/store/data/Run2012A/Photon/RECO/23May2012-v2/0000/C4EA840B-C4A5-E111-9B4E-9C8E991A143E.root',

    
    ),
          eventsToProcess = cms.untracked.VEventRange (
    #options.eventsToProcess
    
    '190895:885:873605097',
    '190906:302:309419123',
    '191247:417:587935144',
    '195776:272:225003419',
    '195930:487:432208932',
    #'194050:617:584876809',
    '195656:172:156861545',
    '195378:164:179871821',
    '195113:477:566362664',
#    '194108:213:184780031',
    '194151:76:74934725',
    )                               
)

process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string ('DiPhotonMet100SUSY_JetReq_2012B_Jun28.root')
)

process.end = cms.EndPath(process.Out)
