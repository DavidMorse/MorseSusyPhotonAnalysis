// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.h
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.h,v 1.3 2011/10/27 13:16:01 dmorse Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <map>

#include "../../../TMVA-v4.1.2/TMVA/Tools.h"
#include "../../../TMVA-v4.1.2/TMVA/Reader.h"
#include "../../../TMVA-v4.1.2/TMVA/MethodCuts.h"

#include "../src/SusyEvent.h"

class SusyEventAnalyzer {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  susy::Event     *event;

  // List of branches
  TBranch        *b_Event;

  SusyEventAnalyzer(TTree *tree=0);
  virtual ~SusyEventAnalyzer();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();                          // event loop for main analysis
  virtual void     DR03();
  virtual void     Pileup();
  virtual void     Filter();
  virtual void     PhotonId();
  virtual void     HiggsAna();
  TMVA::Reader *reader;


  // utility functions
  void CategorizeEvents(susy::Photon* pho1, susy::Photon* pho2, float Rho25, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg, bool useTrigger);
  void CategorizeEventsPhoIso(susy::Photon* pho1, susy::Photon* pho2, float Rho25, bool &gogg, bool &goee, bool &goeg, bool &goff, bool &gogammafake, bool &gogf, bool &gofg, bool useTrigger);
  bool isSameObject(TVector3& p1, TLorentzVector& p2, float dR_Cut);
  bool isSameObject(TLorentzVector& p1, TVector3& p2, float dR_Cut);
  bool isSameObject(TLorentzVector& p1, TLorentzVector& p2, float dR_Cut);
  bool isSameObject(TVector3& p1, TVector3& p2, float dR_Cut);
  float getDR(TVector3& p1, TLorentzVector& p2);
  float getDR(TLorentzVector& p1, TLorentzVector& p2);
  float getDphi(float p1, float p2);
  float d0correction(TVector3& beamSpot, susy::Track& track) const;
  void IncludeAJson(std::string jsonfile);  // Call to pull in a json file 
  bool isInJson(Int_t run,Int_t lumi);      // JSON based good run list cut...
  bool PassTrigger(TString v); // return true if path v is fired
  bool PassTriggers(); // return true if any of names in hltNames are fired
  float GetDiEmPt(susy::Photon* Pho1, susy::Photon* Pho2);//calculates and returns diEmPt
  float GetDiEmPt(TLorentzVector Pho1, TLorentzVector Pho2);
  float GetDiJetPt(susy::PFJet* Jet1, susy::PFJet* Jet2);//calculates and returns diJetPt
  float GetTriEmPt(susy::Photon* Pho1, susy::Photon* Pho2, susy::Photon* Pho3);//calculates and returns TriEmPt
  float InvariantMass(TLorentzVector P1, TLorentzVector P2);//calculates and returns Invariant Mass  
  float InvariantMass(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3);//calculates and returns Invariant Mass of three object system
  float TransverseMass(TLorentzVector part, TVector2 met);//calculates and returns Transverse Mass
  float TransverseMassSquare3body(TLorentzVector lep1, TLorentzVector lep2, TVector2 met);//calculates and returns Transverse Mass Squared for 3 object system
  float InvariantMass(Float_t ggInvMass, TLorentzVector part, TVector2 met);//calculates and returns Invariant Mass
  float MT2(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3);//calculates and returns mt2 for 3 body system, e.g. (e+)(e-)(gamma)
  bool tooClosePhi(TVector3& p1, TVector3& p2, float phi_Cut);//checks whether gg or ff are too close to each other in phi
  bool tooClosePhi(TVector3& p1, TVector2& p2, float phi_Cut);//checks whether gg or ff are too close to MET in phi
  std::pair<float,float> GetMetReweight(float diEMPT,std::string type,std::vector< std::pair<float,float> > binEE[4],std::vector< std::pair<float,float> > binFF[4],std::vector< std::pair<float,float> > binFG[4],std::vector< std::pair<float,float> > binGF[4],std::vector< std::pair<float,float> > binGammaFake[4],std::vector< std::pair<float,float> > binEEsidebandLowJet[4],std::vector< std::pair<float,float> > binEEsidebandHighJet[4],int numJets);//gets reweighting for met plots from diempt ratio
  //pair<float,float> GetMetReweight(float diEMPT,string type,vector< pair<float,float> > binEE[4],vector< pair<float,float> > binFF[4],vector< pair<float,float> > binEEsidebandLowJet[4],vector< pair<float,float> > binEEsidebandHighJet[4],int numJets);
  float GetRazrMr(susy::Photon* Pho1, susy::Photon* Pho2);
  float GetRazrR2(susy::Photon* Pho1, susy::Photon* Pho2, susy::MET* Met);
  void MatchPhosToJets(susy::Photon* pOne, susy::Photon* pTwo, std::vector<susy::PFJet*> jets, susy::PFJet* &jet1, susy::PFJet* &jet2, bool &hasdijetpt, float dR);
  float GetAlphaT(TLorentzVector pOne,TLorentzVector pTwo);
  float GetPhotonLessHt(float Ht, TLorentzVector pOne, TLorentzVector pTwo);
  std::vector<float> GetMetAndInvMassProbs(float invmass, float met, TH1F* SMmet, TH1F* SMinvmass, TH1F* BGmet, TH1F* BGinvmass, TH1F* NewHiggsMet);
  // parameter configuration functions
  void Initialize();         // global variables needed to be initialized just once
  void InitializePerEvent(); // global variables needed to be initialized per event
  void SetDataset(TString& v) {          ds = v; }
  void SetPrintInterval(int v) {         printInterval = v; }
  void SetPrintLevel(int v) {            printLevel = v; }
  void SetProcessNEvents(int v) {        processNEvents = v; }
  void SetUseTrigger(bool v) {           useTrigger = v; }
  void SetUseJSON(bool v) {              useJSON = v; }
  void AddHltName(TString v) {           hltNames.push_back(v); }
  void SetFilter(bool v) {               enableFilter = v; }
  void SetFilteredFileName(TString v) {  filtered_file_name = v; }
  void SetOutputEventNumbers(bool v) {   outputEventNumbers = v; }
  void DoRhoCorrection(bool v) {         doRhoCorrection = v; }
  void DoNvertexCorrection(bool v) {     doNVertexCorrection = v; }
  void SetDR03Rho25Corr(float ecal, float hcal, float track){ PUCorr_ECAL=ecal;PUCorr_HCAL=hcal;PUCorr_TRACK = track; }
  void SetPFisoRho25Corr(float ch, float nh, float ph){ PUCorr_chargedHadron=ch;PUCorr_neutralHadron=nh;PUCorr_photon=ph; }
  void isFastSim(bool isFast){FastSim=isFast;}
  void isFullSim(bool isFull){FullSim=isFull;}
  float FastSimSmear(susy::Photon* pho, TRandom* rand);
  float FullSimSmear(susy::Photon* pho, TRandom* rand);
  susy::PFJet* JECup  (susy::PFJet* jet);
  susy::PFJet* JECdown(susy::PFJet* jet);
  TVector2 CalcMet(std::vector<susy::PFJet*> jets, susy::Photon* pho1, susy::Photon* pho2, std::vector<susy::Electron*> eles, std::vector<susy::Muon*> mus);
  TVector2 CalcMetFromPFandJets(std::vector<susy::PFJet*> jets, TVector2 pfMet);
  float GetPhoScaleFactor(susy::Photon* pho);
  float GetEleScaleFactor(susy::Electron* ele);
  float GetMuScaleFactor(susy::Muon* mu);


 private:

  TString ds;               // dataset name to be used for output hitfile name

  // printLevel
  // 0 : default - no printout
  // 1 : print functional step in every event
  // 2 : print values in collections
  int printLevel;           // print frequency

  int printInterval;        // print level for event content: defined in Event.h
  int processNEvents;       // number of events to be processed
  bool useTrigger;          // flag for using trigger bit selection.
  bool useJSON;             // flag for using JSON selection
  std::vector<TString> hltNames;          // HLT trigger path names
  bool enableFilter;        // filter events of interest
  TString filtered_file_name; // filtered output file name
  bool outputEventNumbers;  // print run and event numbers for gg, eg, ee, ff to txt 
  bool doRhoCorrection;
  bool doNVertexCorrection;
  float PUCorr_ECAL;
  float PUCorr_HCAL;
  float PUCorr_TRACK;
  float PUCorr_chargedHadron;
  float PUCorr_neutralHadron;
  float PUCorr_photon;
  bool FastSim;
  bool FullSim;

  typedef std::map<int,std::map<int,bool> > RunLumiFlagHolder;  //define map that holds json list
  RunLumiFlagHolder goodrunlumilist;  // instantiate it

};

#endif

#ifdef SusyEventAnalyzer_cxx
SusyEventAnalyzer::SusyEventAnalyzer(TTree *tree)
{
  if (tree == 0) {
    std::cout << "Error!!! There is no file containing a tree." << std::endl;
  }
  Init(tree);
  Initialize();
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t SusyEventAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry,0);
}
Long64_t SusyEventAnalyzer::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
  }
  return centry;
}

void SusyEventAnalyzer::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  //fChain->SetMakeClass(1);

  event = new susy::Event;

  fChain->SetBranchAddress("susyEvent", &event, &b_Event);
  //fChain->SetBranchStatus("l1Map*",0);
  //fChain->SetBranchStatus("superClusters*",0);
  //fChain->SetBranchStatus("clusters*",0);
  //fChain->SetBranchStatus("electrons*",0);
  //fChain->SetBranchStatus("muons*",0);
  //fChain->SetBranchStatus("caloJets*",0);
  //fChain->SetBranchStatus("jptJets*",0);
  //fChain->SetBranchStatus("generalTracks*",0);
}

void SusyEventAnalyzer::Initialize() {

  ds = "test";
  printLevel = 0;
  printInterval = 1000;
  processNEvents = -1;
  useTrigger = false;
  useJSON = false;
  enableFilter = false;
  filtered_file_name = "filtered.root";

 reader = new TMVA::Reader( "!Color:!Silent" );
  /*
  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;
  
  // --- Cut optimisation
  Use["Cuts"]            = 1;
  Use["CutsD"]           = 1;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  // 
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 1;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 1;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 1;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 1; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 1; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
  Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine 
  Use["SVM"]             = 1;
  // 
  // --- Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  // 
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 1;
  // ---------------------------------------------------------------
  Use["Plugin"]          = 0;
  Use["Category"]        = 0;
  Use["SVM_Gauss"]       = 0;
  Use["SVM_Poly"]        = 0;
  Use["SVM_Lin"]         = 0;

  */
  /* Float_t* varsumEt, vardPhi, varMETdPhiLead, varMETdPhiTrail, varAlphaT, varPhotonLessHT, varDiEMPt, varInvarMass, varMet, varMR;

  reader->AddVariable( "sumEt", varsumEt);
  reader->AddVariable( "dPhi", vardPhi);
  reader->AddVariable( "METdPhiLead", varMETdPhiLead);
  reader->AddVariable( "METdPhiTrail", varMETdPhiTrail);
  reader->AddVariable( "AlphaT", varAlphaT);
  reader->AddVariable( "PhotonLessHT", varPhotonLessHT);
  reader->AddVariable( "DiEMPt", varDiEMPt);
  reader->AddVariable( "InvarMass", varInvarMass);
  reader->AddVariable( "Met", varMet);
  reader->AddVariable( "MR", varMR);
  */
  //TString dir    = "../../../TMVA-v4.1.2/TMVA/test/weights/";
  //TString prefix = "TMVAClassification";
  /*
  // Book method(s)
  for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
    if (it->second) {
      TString methodName = TString(it->first) + TString(" method");
      TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
      reader->BookMVA( methodName, weightfile ); 
    }
  }
  */
  /*
  // Book output histograms
  UInt_t nbin = 100;
  TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
  TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
  TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
  TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
  TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);
  
  if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
  if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
  if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
  if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
  if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
  if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
  if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
  if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
  if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
  if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
  if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
  if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
  if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
  if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
  if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
  if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
  if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
  if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
  if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
  if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
  if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
  if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
  if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
  if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
  if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
  if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
  if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
  if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
  if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
  if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

  // PDEFoam also returns per-event error, fill in histogram, and also fill significance
  if (Use["PDEFoam"]) {
    histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
    histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
    histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }
  
  // Book example histogram for probability (the other methods are done similarly)
  TH1F *probHistFi(0), *rarityHistFi(0);
  if (Use["Fisher"]) {
     probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
     rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }
  */
  // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
  //TFile *input(0);



}

void SusyEventAnalyzer::IncludeAJson(std::string jsonfile) {


// Fairly primitive brute force json parser -- opens the json file named in the argument
// and adds that to the goodrunlumilist map.  Overlapping jsons are merged inclusively.

 char thing;

 ifstream jsonInput;

 std::cout << "Sucking in Json file: " << jsonfile << " which includes: " << std::endl;

 jsonInput.open(jsonfile.c_str());
 
 if (!jsonInput.good()) {
   std::cout << "Problem reading Json file...  Didn't suck anything in... " << std::endl;
   return;
 }
 
 jsonInput.get(thing);
 
 while (jsonInput.good()) {
   if (thing=='{') {  // start of list
     while (thing != '}') {
       int runnum;
       if (thing == '"') {
         std::string srunnum;
         jsonInput.get(thing); // get stuff inside ""

         while (thing != '"') {
           srunnum+=thing; // get stuff inside ""
           jsonInput.get(thing);

	   }
         sscanf(srunnum.c_str(),"%i",&runnum);
         std::cout << " runnum: " << runnum << std::endl;
         bool newrun=true;
         
       } // inside ""
       if (thing == '[') {
          jsonInput.get(thing); // get stuff inside []
	 while (thing != ']') {
           if (thing == '[') {
             jsonInput.get(thing); // get stuff inside series []

             std::string lumiseries;
             int firstlumi,lastlumi;
             while (thing !=']') {
               lumiseries+=thing;
                jsonInput.get(thing); // get stuff inside series []
             }
             sscanf(lumiseries.c_str(),"%i,%i",&firstlumi,&lastlumi);
             std::cout << "  lumis  " << firstlumi << " to " << lastlumi << std::endl;

	     // At this point have runnum, first lumi, last lumi -- so can fill map here...
	     for (int l=firstlumi;l<=lastlumi;l++) {
               goodrunlumilist[runnum][l]=true;
	     }

           } // inside actual series []
             jsonInput.get(thing); // get stuff inside []
         }
       } // inside []
         jsonInput.get(thing); // get another char looking for "

     } 
   } // inside {}
    jsonInput.get(thing); // get another char looking for {

 } // EOF 

 jsonInput.close();

}


bool SusyEventAnalyzer::isInJson(Int_t run,Int_t lumi) {

//#ifdef MC
//  return 1;
//#endif

if (goodrunlumilist[run][lumi]) return true;

return false;

}


bool SusyEventAnalyzer::isSameObject(TVector3& p1, TVector3& p2, float dR_Cut) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}

bool SusyEventAnalyzer::isSameObject(TVector3& p1, TLorentzVector& p2, float dR_Cut) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}
bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TVector3& p2, float dR_Cut) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}
bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TLorentzVector& p2, float dR_Cut) {
  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < dR_Cut) return true;
  //if(dR < 0.5) return true;
  return false;
}


float SusyEventAnalyzer::getDR(TVector3& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}
float SusyEventAnalyzer::getDR(TLorentzVector& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  return dR;
}

float SusyEventAnalyzer::getDphi(float p1, float p2)
{
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1 - p2));
  return dPhi;
}


float SusyEventAnalyzer::d0correction(TVector3& beamSpot, susy::Track& track) const {
  float d0 = track.d0() - beamSpot.X()*std::sin(track.phi()) + beamSpot.Y()*std::cos(track.phi());
  return d0;
}

bool SusyEventAnalyzer::PassTrigger(TString path) {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event->hltMap.begin(); it != event->hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second)) && it->second.first==1 ) {
      //std::cout<<"Path "<<it->first<<" passed!  Prescale: "<<it->second.first<<endl;
      pass = true;
      break;
    }
  }
  return pass;
}


bool SusyEventAnalyzer::PassTriggers() {
  bool pass = false;
  for(std::vector<TString>::iterator it = hltNames.begin(); it != hltNames.end(); it++) {
    if(PassTrigger(*it)) {
      pass = true;
      break;
    }
  }
  return pass;
}

bool SusyEventAnalyzer::tooClosePhi(TVector3& p1, TVector3& p2, float phi_Cut) {
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi()));
  if(dPhi<phi_Cut){
    //std::cout<<"p1Phi= "<<p1.Phi()<<"  p2Phi= "<<p2.Phi()<<"  dPhi= "<<dPhi<<endl; 
    return true;
  }
  return false;
}

bool SusyEventAnalyzer::tooClosePhi(TVector3& p1, TVector2& p2, float phi_Cut) {
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi()));
  if(dPhi<phi_Cut){
    //std::cout<<"p1Phi= "<<p1.Phi()<<"  p2Phi= "<<p2.Phi()<<"  dPhi= "<<dPhi<<endl; 
    return true;
  }
  return false;
}

float SusyEventAnalyzer::InvariantMass(TLorentzVector P1, TLorentzVector P2){
  //TLorentzVector P4=Pho1->momentum + Pho2->momentum;
  float InvarMass = (P1 + P2).M();
  return InvarMass;
}

float SusyEventAnalyzer::MT2(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3){
  float mt2 = (P1 + P2 + P3).Mt2();
  return mt2;
}

float SusyEventAnalyzer::InvariantMass(TLorentzVector P1, TLorentzVector P2, TLorentzVector P3){
  //TLorentzVector P4=Pho1->momentum + Pho2->momentum;
  float InvarMass = (P1 + P2 + P3).M();
  return InvarMass;
}

float SusyEventAnalyzer::InvariantMass(Float_t ggInvMass, TLorentzVector part, TVector2 met){

  float dPhi = std::fabs(TVector2::Phi_mpi_pi(part.Phi()-met.Phi()));
  float Mt2 = ggInvMass*ggInvMass + 2*part.Pt()*met.Mod()*(1-cos(dPhi));
  float InvMass=sqrt(Mt2);
  return InvMass;

}


float SusyEventAnalyzer::TransverseMass(TLorentzVector part, TVector2 met){
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(part.Phi()-met.Phi()));
  float Mt2 = 2*part.Pt()*met.Mod()*(1-cos(dPhi));
  float Mt=sqrt(Mt2);
  return Mt;
}


float SusyEventAnalyzer::TransverseMassSquare3body(TLorentzVector lep1, TLorentzVector lep2, TVector2 met){
  TLorentzVector temp = lep1+lep2;
  float dPhi = std::fabs(TVector2::Phi_mpi_pi(temp.Phi()-met.Phi()));
  float Mt2 = 2*temp.Pt()*met.Mod()*(1-cos(dPhi));
  //float Mt=sqrt(Mt2);
  return Mt2;
}

float SusyEventAnalyzer::GetDiEmPt(susy::Photon* Pho1, susy::Photon* Pho2){
  /*float PtX = Pho1->momentum.Pt()*cos(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*cos(Pho2->caloPosition.Phi());
    float PtY = Pho1->momentum.Pt()*sin(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*sin(Pho2->caloPosition.Phi());
    float diEmPt = sqrt(PtX*PtX + PtY*PtY);*/
  float diEmPt = (Pho1->momentum + Pho2->momentum).Pt();
  /*cout<<"PhoOnePT: "<<Pho1->momentum.Pt()<<endl;
    cout<<"PhoTwoPT: "<<Pho2->momentum.Pt()<<endl;
    cout<<"PtX: "<<PtX<<endl;
    cout<<"PtY: "<<PtY<<endl;
    cout<<"diEmPt: "<<diEmPt<<endl;*/
  return diEmPt;
}

float SusyEventAnalyzer::GetDiEmPt(TLorentzVector Pho1, TLorentzVector Pho2){
  /*float PtX = Pho1->momentum.Pt()*cos(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*cos(Pho2->caloPosition.Phi());
    float PtY = Pho1->momentum.Pt()*sin(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*sin(Pho2->caloPosition.Phi());
    float diEmPt = sqrt(PtX*PtX + PtY*PtY);*/
  float diEmPt = (Pho1 + Pho2).Pt();
  /*cout<<"PhoOnePT: "<<Pho1->momentum.Pt()<<endl;
    cout<<"PhoTwoPT: "<<Pho2->momentum.Pt()<<endl;
    cout<<"PtX: "<<PtX<<endl;
    cout<<"PtY: "<<PtY<<endl;
    cout<<"diEmPt: "<<diEmPt<<endl;*/
  return diEmPt;
}

float SusyEventAnalyzer::GetDiJetPt(susy::PFJet* Jet1, susy::PFJet* Jet2){
  /*float PtX = Jet1->momentum.Pt()*cos(Jet1->momentum.Phi()) + Jet2->momentum.Pt()*cos(Jet2->momentum.Phi());
    float PtY = Jet1->momentum.Pt()*sin(Jet1->momentum.Phi()) + Jet2->momentum.Pt()*sin(Jet2->momentum.Phi());
    float diJetPt = sqrt(PtX*PtX + PtY*PtY);*/
  float diJetPt = (Jet1->momentum + Jet2->momentum).Pt();
  /*cout<<"PhoOnePT: "<<Pho1->momentum.Pt()<<endl;
    cout<<"PhoTwoPT: "<<Pho2->momentum.Pt()<<endl;
    cout<<"PtX: "<<PtX<<endl;
    cout<<"PtY: "<<PtY<<endl;
    cout<<"diEmPt: "<<diEmPt<<endl;*/
  return diJetPt;
}

float SusyEventAnalyzer::GetTriEmPt(susy::Photon* Pho1, susy::Photon* Pho2,susy::Photon* Pho3){
  /*float PtX = Pho1->momentum.Pt()*cos(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*cos(Pho2->caloPosition.Phi()) + Pho3->momentum.Pt()*cos(Pho3->caloPosition.Phi());
    float PtY = Pho1->momentum.Pt()*sin(Pho1->caloPosition.Phi()) + Pho2->momentum.Pt()*sin(Pho2->caloPosition.Phi()) + Pho3->momentum.Pt()*sin(Pho3->caloPosition.Phi());*/
  float TriEmPt = (Pho1->momentum + Pho2->momentum + Pho3->momentum).Pt();
  /*cout<<"PhoOnePT: "<<Pho1->momentum.Pt()<<endl;
    cout<<"PhoTwoPT: "<<Pho2->momentum.Pt()<<endl;
    cout<<"PtX: "<<PtX<<endl;
    cout<<"PtY: "<<PtY<<endl;
    cout<<"diEmPt: "<<diEmPt<<endl;*/
  return TriEmPt;
}

void SusyEventAnalyzer::InitializePerEvent() {

}


float SusyEventAnalyzer::GetRazrMr(susy::Photon* Pho1, susy::Photon* Pho2){
  float E1 = Pho1->momentum.E();
  float E2 = Pho2->momentum.E();
  float pz1 = Pho1->momentum.Pz();
  float pz2 = Pho2->momentum.Pz();
  float MR = sqrt((E1+E2)*(E1+E2) - (pz1+pz2)*(pz1+pz2));
  return MR;
}
float SusyEventAnalyzer::GetRazrR2(susy::Photon* Pho1, susy::Photon* Pho2, susy::MET* Met){
  float Mr = GetRazrMr(Pho1,Pho2);
  float met = Met->met();
  float pt1 = Pho1->momentum.Pt();
  float pt2 = Pho2->momentum.Pt();
  float Pt1plusPt2 = pt1+pt2;
  TLorentzVector p = Pho1->momentum + Pho2->momentum;
  float MxdotPx = Met->metX()*p.Px();
  float MydotPy = Met->metY()*p.Py();
  float MdotP = MxdotPx + MydotPy;
  float MT = sqrt((met*Pt1plusPt2 - MdotP)/2);
  float R = MT/Mr;
  float R2 = R*R;
  return R2;
}


float SusyEventAnalyzer::GetAlphaT(TLorentzVector pOne, TLorentzVector pTwo)
{
  float minEt = (( pOne.Et()<pTwo.Et())?pOne.Et():pTwo.Et() );
  float Et2 = (pOne.Et()+pTwo.Et())*(pOne.Et()+pTwo.Et());
  float Px2 = (pOne.Px()+pTwo.Px())*(pOne.Px()+pTwo.Px());
  float Py2 = (pOne.Py()+pTwo.Py())*(pOne.Py()+pTwo.Py());
  float Mt = sqrt(Et2 - Px2 - Py2);
  float alphaT = minEt/Mt;
  return alphaT;
}

float SusyEventAnalyzer::GetPhotonLessHt(float Ht, TLorentzVector pOne, TLorentzVector pTwo)
{
  float pLessHt = Ht - pOne.Et() - pTwo.Et();
  return pLessHt;
}

template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}

std::vector<float> SusyEventAnalyzer::GetMetAndInvMassProbs(float invmass, float met, TH1F* SMmet, TH1F* SMinvmass, TH1F* BGmet, TH1F* BGinvmass, TH1F* NewHiggsMet){
  int metBin = SMmet->FindBin(met);
  int invmassBin = SMinvmass->FindBin(invmass);
  float tempSM = SMmet->GetBinContent(metBin);if(tempSM==0)tempSM=1e-6;
  float tempBG = BGmet->GetBinContent(metBin);if(tempBG==0)tempBG=1e-6;
  float tempNH = NewHiggsMet->GetBinContent(metBin);if(tempNH==0)tempNH=1e-6;
  float SM = tempSM*SMinvmass->GetBinContent(invmassBin);
  float BG = tempBG*BGinvmass->GetBinContent(invmassBin);
  float NewHiggs = tempNH*SMinvmass->GetBinContent(invmassBin);
  std::vector<float> vals;
  vals.push_back(SM);
  vals.push_back(BG);
  vals.push_back(NewHiggs);
  return vals;
}

float SusyEventAnalyzer::FastSimSmear(susy::Photon* pho, TRandom *rand){
  float eta = fabs(pho->caloPosition.Eta());
  float r9 = pho->r9;
  float sigma=0.;
  float smearedPt=0.;
  if(eta<1){
    if(r9>0.94) sigma = 0.009213;//0.0092 in AN
    else sigma = 0.010046;//0.01 in AN
  }
  else{
    if(r9>0.94) sigma = 0.016318;//0.0163 in AN
    else sigma = 0.017477;//0.0175 in AN
  }
  smearedPt=rand->Gaus(1,sigma);
  return smearedPt;
}

float SusyEventAnalyzer::FullSimSmear(susy::Photon* pho, TRandom *rand){
  float eta = fabs(pho->caloPosition.Eta());
  float r9 = pho->r9;
  float sigma=0.;
  float smearedPt=0.;
  if(eta<1){
    if(r9>0.94) sigma = 0.0113;//0.0113 in AN
    else sigma = 0.0109;//0.0109 in AN
  }
  else{
    if(r9>0.94) sigma = 0.0171;//0.0171 in AN
    else sigma = 0.0203;//0.0203 in AN
  }
  smearedPt=rand->Gaus(1,sigma);
  return smearedPt;
}

susy::PFJet* SusyEventAnalyzer::JECup(susy::PFJet* jet){
  //cout<<"JECup function"<<endl;
  susy::PFJet* newJet = new susy::PFJet();
  double x=jet->momentum.Px(), y=jet->momentum.Py(), z=jet->momentum.Pz(), E=jet->momentum.E();
  //TLorentzVector p = TLorentzVector(x,y,z,E);
  newJet->momentum=TLorentzVector(x,y,z,E);
  TLorentzVector newP = newJet->momentum;
  //cout<<"pT in JECup:"<<newP.Pt()<<endl;
  /*
    float eta = jet->momentum.Eta();
    if(eta<0.5)newP*=1.115;
    else if(eta<1.1)newP*=1.114;
    else if(eta<1.7)newP*=1.161;
    else if(eta<2.3)newP*=1.228;
    else if(eta<5)newP*=1.4888;
  */
  newP*=1.02;
  //cout<<"pT in JECup:"<<newP.Pt()<<endl;
  newJet->momentum = newP;
  //cout<<"pT in JECup:"<<newJet->momentum.Pt()<<endl;
  return newJet;
}
susy::PFJet* SusyEventAnalyzer::JECdown(susy::PFJet* jet){
  //cout<<"JECdown function"<<endl;
  susy::PFJet* newJet = new susy::PFJet();
  double x=jet->momentum.Px(), y=jet->momentum.Py(), z=jet->momentum.Pz(), E=jet->momentum.E();
  //TLorentzVector p = TLorentzVector(x,y,z,E);
  newJet->momentum=TLorentzVector(x,y,z,E);
  TLorentzVector newP = newJet->momentum;
  //cout<<"pT in JECdown:"<<newP.Pt()<<endl;
  /*
    float eta = jet->momentum.Eta();
    if(eta<0.5)newP*=0.990;
    else if(eta<1.1)newP*=1.001;
    else if(eta<1.7)newP*=1.032;
    else if(eta<2.3)newP*=1.042;
    else if(eta<5)newP*=1.089;
  */
  newP*=0.98;
  //cout<<"pT in JECdown:"<<newP.Pt()<<endl;
  newJet->momentum = newP;
  //cout<<"pT in JECdown:"<<newJet->momentum.Pt()<<endl;
  return newJet;
}

//susy::Muon* SusyEventAnalyzer::MuEffDown(susy::Muon* mu){
//
//}

TVector2 SusyEventAnalyzer::CalcMet(std::vector<susy::PFJet*> jets, susy::Photon* pho1, susy::Photon* pho2, std::vector<susy::Electron*> eles, std::vector<susy::Muon*> mus){
  TLorentzVector p;
  TVector2 met;
  for(std::vector<susy::PFJet*>::iterator jet_it = jets.begin(); jet_it != jets.end(); jet_it++){
    p += (*jet_it)->momentum;
  }
  
  p += pho1->momentum;
  p += pho2->momentum;
  
  for(std::vector<susy::Electron*>::iterator ele_it = eles.begin(); ele_it != eles.end(); ele_it++){
    p += (*ele_it)->momentum;
  }
  for(std::vector<susy::Muon*>::iterator mu_it = mus.begin(); mu_it != mus.end(); mu_it++){
    p += (*mu_it)->momentum;
  }
  met.Set(p.Px(),p.Py());
  return met;
}
TVector2 SusyEventAnalyzer::CalcMetFromPFandJets(std::vector<susy::PFJet*> jets, TVector2 pfMet){
  TLorentzVector p;
  TVector2 met;
  for(std::vector<susy::PFJet*>::iterator jet_it = jets.begin(); jet_it != jets.end(); jet_it++){
    p += (*jet_it)->momentum;
  }
  met.Set(p.Px(),p.Py());
  if(printLevel>0)cout<<"PFmet px:"<<pfMet.Px()<<" py:"<<pfMet.Py()<<endl;
  if(printLevel>0)cout<<"  met px:"<<met.Px()<<" py:"<<met.Py()<<endl;
  met=pfMet-met;
  if(printLevel>0)cout<<"addition"<<endl;
  if(printLevel>0)cout<<"  met px:"<<met.Px()<<" py:"<<met.Py()<<endl;
  return met;
}

float SusyEventAnalyzer::GetPhoScaleFactor(susy::Photon* pho){
  float pT = pho->momentum.Pt();
  float eta = pho->caloPosition.Eta();

  if(pT<0){cout<<"--------- Photon pT less than 1 !!!!!!! -----------"<<endl;return -1.;}
  if(pT<20)return .99;
  else if(pT<30) return .99;
  else  return .99;


}

float SusyEventAnalyzer::GetEleScaleFactor(susy::Electron* ele){
  float pT = ele->momentum.Pt();
  float eta = fabs(ele->momentum.Eta());

  if(pT<15){cout<<"--------- Electron pT less than 15 !!!!!!! -----------"<<endl;return -1.;}
  if(fabs(eta)>2.4){cout<<"--------- Electron |eta| greater than 2.4 !!!!!!! -----------"<<endl;return -1.;}
  if(pT>10.00 && pT<15.00 && eta>0.00 && eta<0.80)return 0.914;
  else if(pT>15.00 && pT<20.00 && eta>0.00 && eta<0.80)return 1.001;
  else if(pT>20.00 && pT<30.00 && eta>0.00 && eta<0.80)return 1.023;
  else if(pT>30.00 && pT<40.00 && eta>0.00 && eta<0.80)return 1.010;
  else if(pT>40.00 && pT<50.00 && eta>0.00 && eta<0.80)return 1.011;
  else if(pT>50.00 /*&& pT<200.00*/ && eta>0.00 && eta<0.80)return 1.010;
  else if(pT>10.00 && pT<15.00 && eta>0.80 && eta<1.44)return 0.899;
  else if(pT>15.00 && pT<20.00 && eta>0.80 && eta<1.44)return 0.984;
  else if(pT>20.00 && pT<30.00 && eta>0.80 && eta<1.44)return 0.996;
  else if(pT>30.00 && pT<40.00 && eta>0.80 && eta<1.44)return 0.997;
  else if(pT>40.00 && pT<50.00 && eta>0.80 && eta<1.44)return 0.995;
  else if(pT>50.00 /*&& pT<200.00*/ && eta>0.80 && eta<1.44)return 1.000;
  else if(pT>10.00 && pT<15.00 && eta>1.44 && eta<1.56)return 1.126;
  else if(pT>15.00 && pT<20.00 && eta>1.44 && eta<1.56)return 0.964;
  else if(pT>20.00 && pT<30.00 && eta>1.44 && eta<1.56)return 1.070;
  else if(pT>30.00 && pT<40.00 && eta>1.44 && eta<1.56)return 1.007;
  else if(pT>40.00 && pT<50.00 && eta>1.44 && eta<1.56)return 0.991;
  else if(pT>50.00 /*&& pT<200.00*/ && eta>1.44 && eta<1.56)return 0.997;
  else if(pT>10.00 && pT<15.00 && eta>1.56 && eta<2.00)return 0.877;
  else if(pT>15.00 && pT<20.00 && eta>1.56 && eta<2.00)return 1.016;
  else if(pT>20.00 && pT<30.00 && eta>1.56 && eta<2.00)return 0.993;
  else if(pT>30.00 && pT<40.00 && eta>1.56 && eta<2.00)return 0.999;
  else if(pT>40.00 && pT<50.00 && eta>1.56 && eta<2.00)return 1.007;
  else if(pT>50.00 /*&& pT<200.00*/ && eta>1.56 && eta<2.00)return 1.007;
  else if(pT>10.00 && pT<15.00 && eta>2.00 && eta<2.50)return 1.044;
  else if(pT>15.00 && pT<20.00 && eta>2.00 && eta<2.50)return 1.044;
  else if(pT>20.00 && pT<30.00 && eta>2.00 && eta<2.50)return 1.042;
  else if(pT>30.00 && pT<40.00 && eta>2.00 && eta<2.50)return 1.032;
  else if(pT>40.00 && pT<50.00 && eta>2.00 && eta<2.50)return 1.024;
  else if(pT>50.00 /*&& pT<200.00*/ && eta>2.00 && eta<2.50)return 1.013;
  else {cout<<"---------Can't find electron scale factor!!! ------"<<endl;return -1;}
}

float SusyEventAnalyzer::GetMuScaleFactor(susy::Muon* mu){
  float pT = mu->momentum.Pt();
  float eta = mu->momentum.Eta();

  if(pT<15){cout<<"--------- Muon pT less than 15 !!!!!!! -----------"<<endl;return -1.;}
  if(fabs(eta)>2.4){cout<<"--------- Muon |eta| greater than 2.4 !!!!!!! -----------"<<endl;return -1.;}
  //cout<<"muon eta: "<<eta<<endl;
  if(eta<-2.1)return 1.061544;
  else if(eta<-1.6)return 1.00709;
  else if(eta<-1.2)return 1.003202;
  else if(eta<-.9)return 1.001565;
  else if(eta<-.6)return 1.001009;
  else if(eta<-.3)return 0.9992277;
  else if(eta<-.2)return 0.9973085;
  else if(eta<.2)return 0.9985061;
  else if(eta<.3)return 1.000557;
  else if(eta<.6)return 1.001074;
  else if(eta<.9)return 1.003233;
  else if(eta<1.2)return 1.00432;
  else if(eta<1.6)return 1.001738;
  else if(eta<2.1)return 1.006988;
  else if(eta>=2.1)return 1.068847;
  else {cout<<"---------Can't find muon scale factor!!! ------"<<endl;return -1;}

}

#endif // #ifdef SusyEventAnalyzer_cxx
