{
  cout << "Loading namespace ejungwoo" << endl;

  gSystem -> Setenv("EJUNGWOOINPUTPATH","PWD/input");
  gROOT -> LoadMacro("PWD/ejungwoo.h");
  //gROOT -> ProcessLine("using namespace ejungwoo;"); // for automatic using namespace call
}
