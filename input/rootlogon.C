// This rootlogon.C is copied and modified from namespace_ejungwoo/input/rootlogon.C
{
  cout << "Loading namespace ejungwoo" << endl;

  if (KEBIPATHISSET) {
    gROOT -> LoadMacro("PWD/KBGlobal.hh");
    gROOT -> LoadMacro("PWD/KBParameterContainer.hh");
    gROOT -> LoadMacro("PWD/KBParameterContainer.cc");
  }
  gSystem -> Setenv("NSEJWINPUTPATH","PWD/input");
  gROOT -> LoadMacro("PWD/ejungwoo.h");
  gROOT -> ProcessLine("using namespace ejungwoo;"); // for automatic using namespace call
}
