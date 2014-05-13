void HiggsTagging()
{

  gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/lib/libRooFit.so");
  gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/lib/libRooFitCore.so");

  gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/include -I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/root/5.32.03-cms9/include");
  
  //gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.34.01-cms8/lib/libRooFit.so");
  //gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.34.01-cms8/lib/libRooFitCore.so");
  //gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.34.01-cms8/include -I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/root/5.34.01-cms8/include");
  
  gROOT->LoadMacro("HiggsTools.C++");

  HiggsTools MakePlots;

  TStopwatch ts;
  ts.Start();
  MakePlots.Loop();
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}






