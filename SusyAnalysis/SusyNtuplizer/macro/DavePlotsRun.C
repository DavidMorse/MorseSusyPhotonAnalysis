void DavePlotsRun()
{

  //load source code and precompiled shared libraries
  gSystem->Load("/afs/cern.ch/user/d/dmorse/scratch0/RA3/2012/cms533v1/CMSSW_5_3_3/src/PhysicsTools/lib/libTagAndProbe.so");
  gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/lib/libRooFit.so");
  gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/lib/libRooFitCore.so");
  //gSystem->Load("libRooFit.so");
  //gSystem->Load("libRooFitCore.so");


  gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/include -I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/root/5.32.03-cms9/include");
  gROOT->LoadMacro("DavePlots.C++");

  DavePlots MakePlots;

  TStopwatch ts;
  ts.Start();
  MakePlots.Loop();
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
