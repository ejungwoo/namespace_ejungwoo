{
  cout << "Loading namespace ejungwoo" << endl;

  if (KEBIPATHISSET)
    gROOT -> LoadMacro("PWD/input/KBParameterContainer.h");
  gSystem -> Setenv("NSEJWINPUTPATH","PWD/input");
  gROOT -> LoadMacro("PWD/ejungwoo.h");
  //gROOT -> ProcessLine("using namespace ejungwoo;"); // for automatic using namespace call
}
