#ifndef KBPARAMETERCONTAINER_HH
#define KBPARAMETERCONTAINER_HH

#include "TObjArray.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TFormula.h"

#include <fstream>
#include <iostream>
#include <string>
//#include <sstream>
#include <stdlib.h>
#include <vector>

using namespace std;

/**
 * List of parameters <[parameter name], [parameter values]>
 * parameter type features Bool_t, Int_t, Double_t, TString.
 *
 * Structure of parameter file should be list of : [name] [type_initial] [value]
 * Each elemets are divided by space.
 * Comments are used by #.
 *
 * @param name Name of parameter file with no space
 * @param type_initial
 *   - b for Bool_t
 *   - i for Int_t
 *   - d for Double_t
 *   - s for TString
 * @param value Value of parameter. TString value do not accept space.
 *
 * ex)\n
 * \#example parameter file\n
 * worldSize    d  1000   # [mm]\n
 * nTbs         i  512\n
 * specialFile  s  /path/to/specialFile.dat  \#special\n 
 *
 * ============================================================================
 *
 * With fDebugMode true, KBParameterContainer will not terminate at attempt of
 * getting non-existing paramter, but print message and create empty parameter.
 *
*/

class KBParameterContainer : public TObjArray
{
  public:
    KBParameterContainer(Bool_t debug = false);
    KBParameterContainer(const char *parName, Bool_t debug = false); ///< Constructor with input parameter file name
    virtual ~KBParameterContainer() {}

    void SaveAs(const char *filename, Option_t *option = "") const;

    void SetDebugMode(Bool_t val = true);

    virtual void Print(Option_t *option ="") const;

    /**
     * Add parameter by given fileName.
     * If fileName does not include path, file will be searched in path/to/KEBI/input.
     *
     * fileName will also be registered as parameter. 
     * If parNameForFile is not set, parameter name will be set as 
     * INPUT_PARAMETER_FILE[fNumInputFiles]
    */
    virtual Int_t AddFile(TString fileName, TString parNameForFile = "");
    virtual Int_t AddPar(KBParameterContainer *parc, TString parNameForFile = "");
    Int_t GetNumInputFiles(); ///< Get number of input parameter files

    Bool_t SetPar(string line);
    Bool_t SetPar(TString name, Bool_t val, Bool_t overwrite = false);       ///< Set Bool_t   type parameter with given name
    Bool_t SetPar(TString name, Int_t val, Bool_t overwrite = false);        ///< Set Int_t    type parameter with given name
    Bool_t SetPar(TString name, Double_t val, Bool_t overwrite = false);     ///< Set Double_t type parameter with given name
    Bool_t SetPar(TString name, TString val, Bool_t overwrite = false);      ///< Set TString  type parameter with given name
    Bool_t SetPar(TString name, const char* val, Bool_t overwrite = false);  ///< Set TString  type parameter with given name
    Bool_t SetParColor(TString name, TString valColor, Bool_t overwrite = false);

    Int_t    GetParN     (TString name);                ///< Get number of parameters in array of given name.
    Bool_t   GetParBool  (TString name, Int_t idx=-1);  ///< Get Bool_t   type parameter by given name.
    Int_t    GetParInt   (TString name, Int_t idx=-1);  ///< Get Int_t    type parameter by given name.
    Double_t GetParDouble(TString name, Int_t idx=-1);  ///< Get Double_t type parameter by given name.
    TString  GetParString(TString name, Int_t idx=-1);  ///< Get TString  type parameter by given name.

    vector<Bool_t>   GetParVBool  (TString name);
    vector<Int_t>    GetParVInt   (TString name);
    vector<Double_t> GetParVDouble(TString name);
    vector<TString>  GetParVString(TString name);

    TVector3 GetParV3(TString name);

    Double_t GetParX(TString name) { return GetParDouble(name,0); }
    Double_t GetParY(TString name) { return GetParDouble(name,1); }
    Double_t GetParZ(TString name) { return GetParDouble(name,2); }

    Int_t GetParWidth(TString name)   { return GetParInt(name); }
    Int_t GetParColor(TString name)   { return GetParInt(name); }
    Double_t GetParSize(TString name) { return GetParDouble(name); }

    Bool_t CheckPar(TString name);

    void ReplaceEnvironmentVariable(TString &val);

    bool IsEmpty() const;

  private:
    Bool_t fDebugMode = false;
    Int_t fNumInputFiles = 0;
};

KBParameterContainer::KBParameterContainer(Bool_t debug)
:TObjArray(), fDebugMode(debug)
{
  fName = "ParameterContainer";
}

KBParameterContainer::KBParameterContainer(const char *parName, Bool_t debug)
:KBParameterContainer(debug)
{
  AddFile(TString(parName));
}

void KBParameterContainer::SetDebugMode(Bool_t val) { fDebugMode = val; }

void KBParameterContainer::SaveAs(const char *fileName, Option_t *) const
{
  TString fielName0(fileName);
  if (!(fielName0.EndsWith(".conf") || fielName0.EndsWith(".par")))
    Print(Form("%s.conf",fileName));
  else
    Print(fileName);
}

void KBParameterContainer::ReplaceEnvironmentVariable(TString &val)
{
  if (val[0] == '$') {
    TString env = val;
    Ssiz_t nenv = env.First("/");
    env.Resize(nenv);
    env.Remove(0, 1);
    TString path = getenv(env);
    val.Replace(0, nenv, path);
  }
}

Int_t KBParameterContainer::AddFile(TString fileName, TString parNameForFile)
{
  TString fileNameFull;

  if (fileName[0] != '/' && fileName[0] != '$' && fileName != '~')
    fileNameFull = TString(gSystem -> Getenv("KEBIPATH")) + "/input/" + fileName;
  else
    fileNameFull = fileName;

  ReplaceEnvironmentVariable(fileNameFull);

  if (TString(gSystem -> Which(".", fileNameFull.Data())).IsNull())
    fileNameFull = TString(gSystem -> Getenv("PWD")) + "/" + fileName;

  if (TString(gSystem -> Which(".", fileNameFull.Data())).IsNull()) {
    cout << "Parameter file " << fileNameFull << " does not exist!" << endl;
    return 0;
  }
  cout << "Adding parameter file " << fileNameFull << endl;

  if (parNameForFile.IsNull())
    parNameForFile = Form("INPUT_FILE_%d", fNumInputFiles);
  else if (parNameForFile.Index("INPUT_FILE")<0)
    parNameForFile = Form("INPUT_FILE_%d", fNumInputFiles);

  fNumInputFiles++;
  SetPar(parNameForFile, fileNameFull);

  Int_t countParameters = 0;

  ifstream file(fileNameFull);
  string line;

  while (getline(file, line)) {
    if (SetPar(line))
      countParameters++;
  }

  if (countParameters == 0) {
    this -> Remove(FindObject(parNameForFile));
    fNumInputFiles--;
  }

  return countParameters;
}

Int_t KBParameterContainer::AddPar(KBParameterContainer *parc, TString parNameForFile)
{
  cout << "Adding parameter container " << parc -> GetName() << endl;

  if (parNameForFile.IsNull())
    parNameForFile = Form("INPUT_PAR_CONTAINER%d", fNumInputFiles);
  fNumInputFiles++;
  //SetPar(parNameForFile, ""); //@todo

  Int_t countParameters = 0;
  Int_t countSameParameters = 0;

  TIter iterator(parc);
  TObject *obj;
  while ((obj = dynamic_cast<TObject*>(iterator())))
  {
    TString name = obj -> GetName();

    TObject *found = FindObject(name);
    if (found != nullptr) {
      if (name.Index("INPUT_FILE")==0)
        ((TNamed *) obj) -> SetName(name+"_");
      else {
        cout << "Parameter with name " << name << " already exist!" << endl;
        ++countSameParameters ;
        continue;
      }
    }

    Add(obj);
    ++countParameters;
  }

  if (countParameters == 0) {
    //this -> Remove(FindObject(parNameForFile));
    fNumInputFiles--;
  }

  return countParameters;
}

Int_t KBParameterContainer::GetNumInputFiles() { return fNumInputFiles; }

void KBParameterContainer::Print(Option_t *option) const
{
  if (TString(option)=="raw") {
    TObjArray::Print();
    return;
  }

  bool printout = true;
  ofstream fileOut;

  TString fileName(option);
  if (fileName.EndsWith(".conf") || fileName.EndsWith(".par"))
    printout = false;

  if (printout) {
    if (fDebugMode) cout << "[" << fName << "]" << "(debug mode)" << " List of parameters :" << endl;
    else cout << "[" << fName << "]" << " List of parameters :" << endl;
  }
  if (!printout) {
    cout << "Writting " << fileName << " as parameter file." << endl;
    fileOut.open(fileName);
    fileOut << "# " << fileName << endl;
    fileOut << "# created from method KBParameterContainer" << endl;
    fileOut << endl;
  }

  TString arrayTitleIsSet;
  Int_t countDownArrayIdx = 0;
  TIter iterator(this);
  TObject *obj;
  while ((obj = dynamic_cast<TObject*>(iterator())))
  {
    bool skipPrint = false;

    TString className = obj -> ClassName();
    TString key = obj -> GetName();
    TString parType;
    TString valueString;
    TString comment;
    Bool_t parameterIsNew = false;
    if (key.Index("NEWPAR")==0) { key.Remove(0,6); parameterIsNew = true; }

    if (className == "TNamed") {
      TNamed *par = (TNamed *) obj;
      valueString = par -> GetTitle();
      if (valueString.Index("AXIS_PARAMETER_")==0) {
        valueString.ReplaceAll("AXIS_PARAMETER_","");
        parType = "axis";
      }
      else if (valueString.Index("INPUT_FILE_")==0) {
        valueString.ReplaceAll("INPUT_FILE_","");
        parType = "file";
      }
      else
        parType = "string";
    }
    else if (className == "TParameter<int>") {
      TParameter<Int_t> *par = (TParameter<Int_t> *) obj;
      valueString += par -> GetVal();
      parType = "int";
      if (key.Index("NUM_VALUES_")==0) {
        countDownArrayIdx = valueString.Atoi();
        arrayTitleIsSet = "[]";
        skipPrint = true;
      }
      if (key.Index("VECTOR3_")==0) {
        countDownArrayIdx = 3;
        arrayTitleIsSet = "v3";
        skipPrint = true;
      }
    }
    else if (className == "TParameter<double>") {
      TParameter<Double_t> *par = (TParameter<Double_t> *) obj;
      valueString += par -> GetVal();
      parType = "double";
    }
    else if (className == "TParameter<bool>") {
      TParameter<Bool_t> *par = (TParameter<Bool_t> *) obj;
      valueString = ( (par->GetVal()==true) ? "true" : "false");
      parType = "bool";
    }
    else
      continue;

    if (parameterIsNew)
      comment += "YOU MUST MODIFY THIS PARAMETER VALUE";

    if (!comment.IsNull())
      comment = TString("#  ") + comment;

    if (!skipPrint) {
      ostringstream ssLine;
      bool thisIsNewLine = false;
      if (countDownArrayIdx!=0) {
        --countDownArrayIdx;
        key = TString(key(0,key.First("[")));
        if (!arrayTitleIsSet.IsNull()) {
          if (arrayTitleIsSet=="[]")
            arrayTitleIsSet = parType + arrayTitleIsSet;
          thisIsNewLine = true;
          ssLine << left << setw(28) << key << setw(10) << arrayTitleIsSet << valueString;
          arrayTitleIsSet = "";
          if (countDownArrayIdx==0)
            ssLine << comment << endl;
        }
        else if (countDownArrayIdx==0) ssLine << ",  " << valueString << comment << endl;
        else ssLine << ",  " << valueString;
      }
      else {
        thisIsNewLine = true;
        ssLine << left << setw(28) << key << setw(10) << parType << valueString << comment << endl;
      }

      if (printout)
        if (thisIsNewLine) cout << ssLine.str();
        else cout << ssLine.str();
      else
        fileOut << ssLine.str();
    }
  }
}


Bool_t KBParameterContainer::SetPar(string line)
{
  if (line.find("#") == 0)
    return false;

  TString parName;
  TString parType;

  istringstream ss(line);
  ss >> parName >> parType;
  parType.ToLower();

  Bool_t overwrite = false;
  if (parType == "o" || parType == "overwrite") {
    overwrite = true;
    ss >> parType;
    parType.ToLower();
  }

  if (parType == "f" || parType == "file") {
    TString val;
    ss >> val;
    ReplaceEnvironmentVariable(val);
    AddFile(val, parName);
  }
  else if (parType == "v3" || parType == "vector3" || parType == "tvector3" || parType == "kbvector3") {
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("3")+2, valuesInString.Sizeof()-valuesInString.First("3")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (numValues != 3)
      return false;
    SetPar(TString("VECTOR3_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal)
      SetPar(parName+"["+iVal+"]",TString(((TObjString *) valueTokens->At(iVal))->GetString()).Atof(),overwrite);
  }
  else if (parType == "b" || parType == "bool" || parType == "bool_t") {
    TString valueBoolean;
    ss >> valueBoolean;
    valueBoolean.ToLower();
    Bool_t val = false;
    if (valueBoolean == "true" || valueBoolean == "1" || valueBoolean == "ktrue")
      val = true;
    SetPar(parName, val, overwrite);
  }
  else if (parType.Index("b[")==0 || parType.Index("bool[")>=0 || parType.Index("bool_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (arrayLength!=0 && numValues > arrayLength) numValues = arrayLength;
    if (arrayLength>1 && numValues==1) { numValues = arrayLength; arrayLength = -1; }
    SetPar(TString("NUM_VALUES_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal) {
      int idx = ((arrayLength<0) ? 0 : iVal);
      TString valueBoolean = TString(((TObjString *) valueTokens->At(idx))->GetString());
      Bool_t val = false;
      if (valueBoolean == "true" || valueBoolean == "1" || valueBoolean == "ktrue")
        val = true;
      SetPar(parName+"["+iVal+"]",val,overwrite);
    }
  }
  else if (parType == "i" || parType == "int" || parType == "int_t" || parType == "w" || parType == "width" || parType == "width_t") {
    Int_t val;
    ss >> val;
    SetPar(parName, val, overwrite);
  }
  else if (parType.Index("i[")==0 || parType.Index("int[")>=0 || parType.Index("int_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (arrayLength!=0 && numValues > arrayLength) numValues = arrayLength;
    if (arrayLength>1 && numValues==1) { numValues = arrayLength; arrayLength = -1; }
    SetPar(TString("NUM_VALUES_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal) {
      int idx = ((arrayLength<0) ? 0 : iVal);
      SetPar(parName+"["+iVal+"]",TString(((TObjString *) valueTokens->At(idx))->GetString()).Atoi(),overwrite);
    }
  }
  else if (parType == "d" || parType == "double" || parType == "double_t" || parType == "size" || parType == "size_t") {
    TString valFormula;
    ss >> valFormula;
    Double_t val = TFormula("formula",valFormula).Eval(0);
    SetPar(parName, val, overwrite);
  }
  else if (parType.Index("d[")==0 || parType.Index("double[")>=0 || parType.Index("double_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (arrayLength!=0 && numValues > arrayLength) numValues = arrayLength;
    if (arrayLength>1 && numValues==1) { numValues = arrayLength; arrayLength = -1; }
    SetPar(TString("NUM_VALUES_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal) {
      int idx = ((arrayLength<0) ? 0 : iVal);
      TString valFormula = TString(((TObjString *) valueTokens->At(idx))->GetString());
      Double_t val = TFormula("formula",valFormula).Eval(0);
      SetPar(parName+"["+iVal+"]",val,overwrite);
    }
  }
  else if (parType == "s" || parType == "string" || parType == "tstring") {
    TString val;
    ss >> val;
    ReplaceEnvironmentVariable(val);
    SetPar(parName, val, overwrite);
  }
  else if (parType.Index("s[")==0 || parType.Index("string[")>=0 || parType.Index("tstring[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (arrayLength!=0 && numValues > arrayLength) numValues = arrayLength;
    if (arrayLength>1 && numValues==1) { numValues = arrayLength; arrayLength = -1; }
    SetPar(TString("NUM_VALUES_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal) {
      int idx = ((arrayLength<0) ? 0 : iVal);
      TString val(((TObjString *) valueTokens->At(idx))->GetString());
      ReplaceEnvironmentVariable(val);
      SetPar(parName+"["+iVal+"]",val,overwrite);
    }
  }
  else if (parType == "a" || parType == "axis" || parType == "kbvector3::axis") {
    TString val;
    ss >> val;
    if (val.Index("AXIS_PARAMETER_")<0)
      val = TString("AXIS_PARAMETER_") + val;
    SetPar(parName, val, overwrite);
  }
  else if (parType.Index("a[")==0 || parType.Index("axis[")>=0 || parType.Index("kbvector3::axis[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (arrayLength!=0 && numValues > arrayLength) numValues = arrayLength;
    if (arrayLength>1 && numValues==1) { numValues = arrayLength; arrayLength = -1; }
    SetPar(TString("NUM_VALUES_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal) {
      int idx = ((arrayLength<0) ? 0 : iVal);
      TString val(((TObjString *) valueTokens->At(idx))->GetString());
      if (val.Index("AXIS_PARAMETER_")<0)
        val = TString("AXIS_PARAMETER_") + val;
      SetPar(parName+"["+iVal+"]",val,overwrite);
    }
  }
  else if (parType == "c" || parType == "color" || parType == "color_t") {
    TString val;
    ss >> val;
    SetParColor(parName, val, overwrite);
  }
  else if (parType.Index("c[")==0 || parType.Index("color[")>=0 || parType.Index("color_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (arrayLength!=0 && numValues > arrayLength) numValues = arrayLength;
    if (arrayLength>1 && numValues==1) { numValues = arrayLength; arrayLength = -1; }
    SetPar(TString("NUM_VALUES_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal) {
      int idx = ((arrayLength<0) ? 0 : iVal);
      TString val(((TObjString *) valueTokens->At(idx))->GetString());
      SetParColor(parName+"["+iVal+"]",val,overwrite);
    }
  }
  else
    return false;

  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, Bool_t val, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      cout << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Add(new TParameter<Bool_t>(name, val));
  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, Int_t val, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      cout << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Add(new TParameter<Int_t>(name, val));
  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, Double_t val, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      cout << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Add(new TParameter<Double_t>(name, val));
  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, TString val, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      cout << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Add(new TNamed(name, val));
  return true;
}

Bool_t KBParameterContainer::SetParColor(TString name, TString valColor, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      cout << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Int_t val = 0;
  if (valColor.IsDec()) {
    val = valColor.Atoi();
  }
  else if (valColor.Index("k")==0) {
    valColor.ReplaceAll("kWhite"  ,"0");
    valColor.ReplaceAll("kBlack"  ,"1");
    valColor.ReplaceAll("kGray"   ,"920");
    valColor.ReplaceAll("kRed"    ,"632");
    valColor.ReplaceAll("kGreen"  ,"416");
    valColor.ReplaceAll("kBlue"   ,"600");
    valColor.ReplaceAll("kYellow" ,"400");
    valColor.ReplaceAll("kMagenta","616");
    valColor.ReplaceAll("kCyan"   ,"432");
    valColor.ReplaceAll("kOrange" ,"800");
    valColor.ReplaceAll("kSpring" ,"820");
    valColor.ReplaceAll("kTeal"   ,"840");
    valColor.ReplaceAll("kAzure"  ,"860");
    valColor.ReplaceAll("kViolet" ,"880");
    valColor.ReplaceAll("kPink"   ,"900");

    if (valColor.Index("+")>0) {
      auto colorCombi = valColor.Tokenize("+");
      int val1 = (((TObjString *) colorCombi->At(0)) -> GetString()).Atoi();
      int val2 = 0;
      if (colorCombi->GetEntriesFast() > 1)
        val2 = (((TObjString *) colorCombi->At(1)) -> GetString()).Atoi();
      val = val1 + val2;
    }
    else if (valColor.Index("-")>0) {
      auto colorCombi = valColor.Tokenize("-");
      int val1 = (((TObjString *) colorCombi->At(0)) -> GetString()).Atoi();
      int val2 = 0;
      if (colorCombi->GetEntriesFast() > 1)
        val2 = (((TObjString *) colorCombi->At(1)) -> GetString()).Atoi();
      val = val1 - val2;
    }
    else
      val = valColor.Atoi();
  }
  else
    return false;

  Add(new TParameter<Int_t>(name, val));
  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, const char* val, Bool_t overwrite)
{
  return SetPar(name, TString(val), overwrite);
}

Int_t KBParameterContainer::GetParN(TString name)
{
  name = TString("NUM_VALUES_") + name;

  TObject *obj = FindObject(name);
  if (obj == nullptr)
    return 0;

  TString className = obj -> ClassName();
  if (className != "TParameter<int>") {
    cout << name << " parameter type is " << className << ", not int!" << endl;
  }

  return ((TParameter<Int_t> *) obj) -> GetVal();
}

Bool_t KBParameterContainer::GetParBool(TString name, Int_t idx)
{
  if (idx>=0) return GetParBool(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    cout << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, false);
    return false;
  }

  TString className = obj -> ClassName();
  if (className != "TParameter<bool>") {
    cout << name << " parameter type is " << className << ", not bool!" << endl;
  }

  return ((TParameter<Bool_t> *) obj) -> GetVal();
}

Int_t KBParameterContainer::GetParInt(TString name, Int_t idx)
{
  if (idx>=0) return GetParInt(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    cout << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, -999);
    return -999;
  }

  TString className = obj -> ClassName();
  if (className != "TParameter<int>") {
    cout << name << " parameter type is " << className << ", not int!" << endl;
  }

  return ((TParameter<Int_t> *) obj) -> GetVal();
}

Double_t KBParameterContainer::GetParDouble(TString name, Int_t idx)
{
  if (idx>=0) return GetParDouble(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    cout << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, -999.999);
    return -999.999;
  }

  TString className = obj -> ClassName();
  if (className != "TParameter<double>") {
    cout << name << " parameter type is " << className << ", not double!" << endl;
  }

  return ((TParameter<Double_t> *) obj) -> GetVal();
}

TString KBParameterContainer::GetParString(TString name, Int_t idx)
{
  if (idx>=0) return GetParString(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    cout << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, "DOES_NOT_EXIST");
    return "DOES_NOT_EXIST";
  }

  TString className = obj -> ClassName();
  if (className != "TNamed") {
    cout << name << " parameter type is " << className << ", not string!" << endl;
  }

  return ((TNamed *) obj) -> GetTitle();
}

vector<Bool_t> KBParameterContainer::GetParVBool(TString name)
{
  vector<Bool_t> array;
  auto npar = GetParN(name);
  for (auto i=0; i<npar; ++i)
    array.push_back(GetParBool(name,i));
  return array;
}

vector<Int_t> KBParameterContainer::GetParVInt(TString name)
{
  vector<Int_t> array;
  auto npar = GetParN(name);
  for (auto i=0; i<npar; ++i)
    array.push_back(GetParInt(name,i));
  return array;
}

vector<Double_t> KBParameterContainer::GetParVDouble(TString name)
{
  vector<Double_t> array;
  auto npar = GetParN(name);
  for (auto i=0; i<npar; ++i)
    array.push_back(GetParDouble(name,i));
  return array;
}

vector<TString> KBParameterContainer::GetParVString(TString name)
{
  vector<TString> array;
  auto npar = GetParN(name);
  for (auto i=0; i<npar; ++i)
    array.push_back(GetParString(name,i));
  return array;
}

TVector3 KBParameterContainer::GetParV3(TString name)
{
  TString xname = name + "[0]";
  TString yname = name + "[1]";
  TString zname = name + "[2]";

  TObject *xobj = FindObject(xname);
  if (xobj == nullptr) {
    cout << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, -999.999);
    return TVector3(-999,-999,-999);
  }

  TString className = xobj -> ClassName();
  if (className != "TParameter<double>") {
    cout << name << " parameter type is " << className << ", not (v3)double!" << endl;
  }

  auto x = ((TParameter<Double_t> *) FindObject(xname)) -> GetVal();
  auto y = ((TParameter<Double_t> *) FindObject(yname)) -> GetVal();
  auto z = ((TParameter<Double_t> *) FindObject(zname)) -> GetVal();

  return TVector3(x,y,z);
}

Bool_t KBParameterContainer::CheckPar(TString name)
{
  if (FindObject(name) != nullptr) return true;
  if (FindObject(TString("VECTOR3_")+name) != nullptr) return true;
  if (FindObject(TString("NUM_VALUES_")+name) != nullptr) return true;
  return false;
}

bool KBParameterContainer::IsEmpty() const
{
  if (GetEntriesFast()>0)
    return false;
  return true;
}

#endif
