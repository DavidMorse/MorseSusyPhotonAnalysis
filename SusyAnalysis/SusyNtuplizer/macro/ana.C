// Jet energy correction is possible at ntuple level.
// $ cd ../jec/JetMETObjects
// $ make
// This will create a shared library in jec/lib
// which is included below as libJetMETObjects.so
//
// Come back to this directory and do
// $ make
// $ root -b -q -l ana.C
// will produce hist_"physics"_"ds".root

//void ana(TString ds="", TString physics="", TString fileName="") {
//void ana(TString ds="Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelVetoOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_PUjetID_NewDoGG_Sep16_allTrigs_bEta2p0_TRASH", TString physics="Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos") {
//void ana(TString ds="Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_EleVeto_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_May16", TString physics="Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos") {
//void ana(TString ds="Photon_cms533v1_Prompt_MatrixMethod_eventCounts_May6", TString physics="Data2012C_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos") {
//void ana(TString ds="MS_1100_MG_720_MN_375", TString physics="SignalMCAcceptance_Bino_2012IDloose_15GevCHfakeCut_sihih012_pixelCutNoEleVeto_2Jetreq_Apr17") {

//void ana(TString ds="DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12-PU_S7_START52_V9-v1_June24",TString physics="Photon"){

//void ana(TString ds="ZZTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1_May6",TString physics="cms533v1"){
//void ana(TString ds="WZTo3LNu_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1_May6",TString physics="cms533v1"){

//void ana(TString ds="GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31", TString physics="Higgs_cms533v1"){
//void ana(TString ds="VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31", TString physics="Higgs_cms533v1"){
//void ana(TString ds="WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31", TString physics="Higgs_cms533v1"){
//void ana(TString ds="TTH_HToGG_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Jul31", TString physics="Higgs_cms533v1"){

//void ana(TString ds="GluGluToHToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_noMegVeto", TString physics="Higgs_cms533v1"){
//void ana(TString ds="VBF_HToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1", TString physics="Higgs_cms533v1"){
//void ana(TString ds="WH_ZH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1", TString physics="Higgs_cms533v1"){
//void ana(TString ds="TTH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1", TString physics="Higgs_cms533v1"){

  
//void ana(TString ds="handMade_DiPho_Born_and_Box_Pt20_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Nov7",TString physics="Photon"){
//void ana(TString ds="handMade_DiPho_Born_and_Box_Pt20to350_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Nov2",TString physics="Photon"){
//void ana(TString ds="handMade_DiPho_Born_and_Box_Pt350_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug8",TString physics="Photon"){
//void ana(TString ds="handMade_GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22",TString physics="Photon"){
//void ana(TString ds="handMade_GJet_Pt40to200_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_PhoIsoFakeDef",TString physics="Photon"){
//void ana(TString ds="handMade_GJet_Pt-200to350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20",TString physics="Photon"){
//void ana(TString ds="handMade_GJet_Pt-350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20",TString physics="Photon"){

//void ana(TString ds="handMade_QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_PhoIsoFakeDef",TString physics="Photon"){
//void ana(TString ds="handMade_QCD_Pt-40to100_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_PhoIsoFakeDef",TString physics="Photon"){
//void ana(TString ds="handMade_QCD_Pt-100to200_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20_PhoIsoFakeDef",TString physics="Photon"){
//void ana(TString ds="handMade_QCD_Pt-200to350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20_PhoIsoFakeDef",TString physics="Photon"){
//void ana(TString ds="handMade_QCD_Pt-350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug26_PhoIsoFakeDef",TString physics="Photon"){

//void ana(TString ds="ZJetToEE_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_TRASH",TString physics="Photon"){
//void ana(TString ds="handMade_Z_ZJet_Pt-20-inf_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22",TString physics="Photon"){



//void ana(TString ds="SMS_TChiWH_WincHgg_2J_TRASH",TString physics="Photon"){
//void ana(TString ds="SMS_TChiZH_ZincHgg_2J",TString physics="Photon"){
//void ana(TString ds="SMS_TChiHH_2b2g_2J",TString physics="Photon"){
//void ana(TString ds="SMS_TChiHH_2Z2g_2J_InvMass",TString physics="Photon"){
//void ana(TString ds="SMS_TChiHH_2W2g_2J_InvMass",TString physics="Photon"){
void ana(TString ds="SMS_TChiHH_2tau2g_2J_InvMass",TString physics="Photon"){


//void ana(TString ds="cms533_WinoNLSP_chargino130_bino1_5_10_15_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino150_bino1_5_10_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino175_bino1_5_10_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino200_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino225_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino250_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino275_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino300_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino325_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino350_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino375_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino400_bino1_5_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino425_bino1_5_10_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino450_bino1_5_10_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino475_bino1_5_10_hw_aaw_full_Sep11",TString physics="Photon"){
//void ana(TString ds="cms533_WinoNLSP_chargino500_bino1_5_10_hw_aaw_full_Sep11",TString physics="Photon"){


//void ana(TString ds="Filter_Ana_HLT_JSON_Two40-25GeVbarrelPhotons_ALL-Higgsino-selectedEvents_mgmuVeto",TString physics="Data2012A_13July_06Aug_recover_Data2012B_13July_Data2012C_Prompt_Runs198022-198903_Runs198941-203742_Data2012D_Prompt_Runs207920-209151_Runs203777-207905"){

//void ana(TString ds="ZGToLLG",TString physics="Photon"){

  TStopwatch ts;
  ts.Start();
  
  gSystem->Load("libSusyEvent.so");

  // Look ../jec/JetMETObjects/README
  gSystem->Load("../jec/lib/libJetMETObjects.so");

  // Printing utility for ntuple variables
  //gROOT->LoadMacro("SusyEventPrinter.cc+");

  // Main analysis code
  gROOT->LoadMacro("SusyEventAnalyzer.cc+");

  // chain of inputs
  TChain* chain = new TChain("susyTree");
  
  std::cout<<"Grabbing Files..."<<std::endl;
   
  //chain->Add(fileName.Data());

  //--------SignalMC-------------
  //2012 - 8TeV
  
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/cms533/Bino/cms533v1_v1/tree_1700_1720_375.root");
 
  //New ID higgs SM files
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/cms533/Higgs/GluGlu/susyEvents_GluGluH_SM_ALL.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/cms533/Higgs/VBF/susyEvents_VBFH_SM_ALL.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/cms533/Higgs/WZH/susyEvents_WZH_SM_ALL.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/cms533/Higgs/TTH/susyEvents_TTH_SM_ALL.root");


  //126 GeV Higgs SM files
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Higgs_cms533v1_GluGluToHToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Higgs_cms533v1_VBF_HToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Higgs_cms533v1_WH_ZH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Higgs_cms533v1_TTH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  //eos files
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/SMhiggs126/Higgs_cms533v1_GluGluToHToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/SMhiggs126/Higgs_cms533v1_VBF_HToGG_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/SMhiggs126/Higgs_cms533v1_WH_ZH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/SMhiggs126/Higgs_cms533v1_TTH_HToGG_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");


  //------Other MC------
  //lots more of these
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12-PU_S7_START52_V9-v1/susyEvents_100_1_clc.root");

  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/cms533/ZZTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/SignalMC/cms533/WZTo3LNu_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1.root");

  //--------cms533v1--------------------
 

  //Apr19
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_50GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_Apr19-SelectedEvents.root");
  //seem to be missing some....add old one
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelCutOnFakes_19499pb-SelectedEvents.root");

  //output of that Apr23
  //first 50GeV fake cut - using this for now
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_50GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_Apr23-SelectedEvents.root");
  //now 15GeV fake cut
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Data2012_Filter_HLTJsonNVertexTwo40-25GeVBarrelPhos_Photon_cms533v1_ReRecoPlusPrompt_15GevCHfakeCut_PixelCutOnFakes_NewDiEMPtBins_19499pb_NewJetMatching_2JetReq_Apr23-SelectedEvents.root");




  //Handmade samples for closure study
 
  //DiPho box and born 20-inf 1080*50000 + 10M +10M = 74M events
  
  ////chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_DiPho_Born_and_Box_Pt20_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_1.root");
  ////chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_DiPho_Born_and_Box_Pt20_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_2.root");
  /*
    chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_DiPho_Born_and_Box_Pt20to350_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug8-SelectedEvents.root");
    chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_DiPho_Born_and_Box_Pt20to350_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug8-SelectedEvents_1.root");
  */
  //DiPho box and born 350-inf 712 jobs @ 10k/job = 7.12M events
  /*
    chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_DiPho_Born_and_Box_Pt350_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug8-SelectedEvents.root");
    chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_DiPho_Born_and_Box_Pt350_doubleEMEnriched_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug8-SelectedEvents_1.root");
  */
  
  //GJet 20-40 2278 jobs @ 3M/job --> 6.834B events
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22.root");
  //GJet 40-200  300M events  + 43 jobs @ 1M/job --> 343M events
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_GJet_Pt40to200_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22-SelectedEvents.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_PhoIsoFakeDef.root");
  //GJet 200-350
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_GJet_Pt-200to350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20.root");
  //GJet 350-inf  300M events  + 43 jobs @ 1M/job --> 343M events
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_GJet_Pt-350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug8.root");

  //QCD 30-40 3888 jobs @3.1M -->12.0528B events
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July12.root");
  //QCD 40-100 3981 jobs @ 1M/job + 8871 jobs @1.5M/job + 244 jobs @1.5M/job  --> 3.981B + 13.305B + .366B  = 17.6535B events
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-40to100_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22-SelectedEvents.root");
  //QCD 100-200 10014 jobs @ 50k/job + 10009 jobs @ 400k/job --> 4.5043B events
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-100to200_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20-SelectedEvents.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-100to200_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20.root");
  //QCD 200-350	6498 jobs @ 200000/job -->
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-200to350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug20.root");
  //QCD 350-inf 4997 jobs @ 60k/job + 10366 jobs @ 200k/job --> 299.82 M + 2.0732 B = 2.37302 B events
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-350_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_Aug19.root");

  //ZtoEE
  /*
  chain->Add("/data/ndpc3/c/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_ZJetToEE_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  chain->Add("/data/ndpc3/c/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_ZJetToEE_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_1.root");
  */

  //FastSim Z+ZJet
  //chain->Add("/data/ndpc2/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_Z_ZJet_Pt-20-inf_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22.root");

  //----------------------------------------------------------------------------------
  
  //aaw EWKino
  
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_cms533_WinoNLSP_chargino130_bino1_5_10_15_hw_aaw_full_Sep11.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_cms533_WinoNLSP_chargino150_bino1_5_10_hw_aaw_full_Sep11.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_cms533_WinoNLSP_chargino175_bino1_5_10_hw_aaw_full_Sep11.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino200_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino225_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino250_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino275_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino300_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino325_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino350_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino375_bino1_5_hw_aaw_full.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino400_bino1_5_hw_aaw_full.root");//fix this
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_cms533_WinoNLSP_chargino425_bino1_5_10_hw_aaw_full_Sep11.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_cms533_WinoNLSP_chargino450_bino1_5_10_hw_aaw_full_Sep11.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_cms533_WinoNLSP_chargino475_bino1_5_10_hw_aaw_full_Sep11.root");
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_cms533_WinoNLSP_chargino500_bino1_5_10_hw_aaw_full_Sep11.root");


  //chain->Add("~/work/RA3/cms533v0/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/susyEvents.root");

  //fix this file
  // chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/2013/AODSIM/electroHiggs/hw/aaw/ntuples/cms533_WinoNLSP_chargino400_bino1_5_hw_aaw_full.root");

  //SMS_TChiWH_WincHgg_2J
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_SMS_TChiWH_WincHgg_2J_Redo.root");
  //SMS_TChiZH_ZincHgg_2J
  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_SMS_TChiZH_ZincHgg_2J.root");
  //SMS_TChiHH_2b2g_2J
  //chain->Add("/data/ndpc3/c/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_SMS_TChiHH_2b2g_2J.root");
  //SMS_TChiHH_2Z2g_2J
  //chain->Add("/data/ndpc3/c/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_SMS_TChiHH_2Z2g_2J.root");
  //SMS_TChiHH_2W2g_2J
  //chain->Add("/data/ndpc3/c/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_SMS_TChiHH_2W2g_2J.root");
 
  //eos files
  //SMS_TChiWH_WincHgg_2J
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/EWKinoWZHgg/Photon_SMS_TChiWH_WincHgg_2J_Redo.root");
  //SMS_TChiZH_ZincHgg_2J
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/EWKinoWZHgg/Photon_SMS_TChiZH_ZincHgg_2J.root");
  //SMS_TChiHH_2b2g_2J
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/EWKinoWZHgg/Photon_SMS_TChiHH_2b2g_2J.root");
  //SMS_TChiHH_2Z2g_2J
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/EWKinoWZHgg/Photon_SMS_TChiHH_2Z2g_2J.root");
  //SMS_TChiHH_2W2g_2J
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/EWKinoWZHgg/Photon_SMS_TChiHH_2W2g_2J.root");
  //SMS_TChiHH_2tau2g_2J
  chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/EWKinoWZHgg/Photon_SMS_TChiHH_2tau2g_2J.root");


  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/handmade/Photon_handMade_QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22.root");
  
  /*
  chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-100_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22.root");
  chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-100_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_1.root");
  chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_QCD_Pt-100_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22_2.root");
  */

  //chain->Add("/data/ndpc3/b/dmorse/RA3/RA3Ntuples/FilteredNtuples/cms533v1/Photon_handMade_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_July22.root");

  //chain->Add("~/work/RA3/cms533v0/CMSSW_5_3_3/src/susyEvents.root");

  //newly re-skimmed ntuples
  //chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/Data/Data2012A_13July_06Aug_recover_Data2012B_13July_Data2012C_Prompt_Runs198022-198903_Runs198941-203742_Data2012D_Prompt_Runs207920-209151_Runs203777-207905_Filter_Ana_HLT_JSON_Two40-25GeVbarrelPhotons_ALL-Higgsino-selectedEvents.root");


//lots more like this
//chain->Add("root://eoscms//eos/cms/store/user/dmorse/RA3Ntuples/MC/cms533v1/ZGToLLG/susyEvents_1000_1_En0.root");

  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);


  // configuration parameters
  // any values given here will replace the default values
  sea->SetDataset(physics+"_"+ds);        // dataset name
  sea->SetPrintInterval(2.5e5);             // print frequency
  sea->SetPrintLevel(0);                  // print level for event contents
  sea->SetOutputEventNumbers(false);      // print run and event numbers
  sea->SetUseTrigger(true);
  sea->SetUseJSON(false);
  sea->DoRhoCorrection(true);
  sea->DoNvertexCorrection(false);
  sea->SetDR03Rho25Corr(0.081,0.022,0.011);//Ecal,Hcal.Track
  sea->SetPFisoRho25Corr(0.031,0.013,0.078);//chargedHadronIso,neutralHadronIso,photonIso
  sea->SetFilter(false);                  // filter events passing final cuts
  sea->isFastSim(true);    //only matters for HiggsAna()
  sea->isFullSim(false);   //only matters for HiggsAna()
  sea->SetProcessNEvents(-1);     // number of events to be processed
  
  // HLT trigger path names
  sea->AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  sea->AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
  sea->AddHltName("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
  sea->AddHltName("HLT_Photon36_R9Id85_Photon22_R9Id85_v");
  

  //fully combined json
  //sea->IncludeAJson("/afs/cern.ch/user/d/dmorse/scratch0/RA3/JsonFiles/2012/Cert_190456-190781_8TeV_13JulReReco_190782-190949_8TeV_May23ReReco_190950-196531_13JulReReco_196532-208686_8TeV_PromptReco_Collisions12_JSON_19499pb.txt");
  //sea->IncludeAJson("/afs/cern.ch/user/d/dmorse/scratch0/RA3/JsonFiles/2012/Cert_2012A.txt");

  //combined first 13JulyReReco:0-190781, then 23MayReReco:190782-190949, then 13JulyReReco:190950-196531, then prompt:196532-999999
  //ignoring 06AugReReco JSON because it covers same run range as 23MayReReco but drops some lumi sections in 190949 - why?

  //sea->Loop();
  //sea->DR03();
  //sea->Pileup();
  //sea->Filter();
  //sea->PhotonId();
  sea->HiggsAna();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
