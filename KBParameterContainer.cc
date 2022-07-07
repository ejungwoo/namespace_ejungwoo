#include "TROOT.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TApplication.h"
#include "KBParameterContainer.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "TFormula.h"
#include "TObjString.h"

using namespace std;

ClassImp(KBParameterContainer)

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

void KBParameterContainer::ReplaceEnvironmentVariable(TString &valInput)
{
  if (valInput[0] == '$') {
    TString env = valInput;
    Ssiz_t nenv = env.First("/");
    env.Resize(nenv);
    env.Remove(0, 1);
    TString path = getenv(env);
    valInput.Replace(0, nenv, path);
  }
}

TString KBParameterContainer::ReplaceVariables(TString &valInput)
{
  ReplaceEnvironmentVariable(valInput);

  if (valInput.Index("@")>=0)
  for (auto ita : {0,1})
  {
    TIter iterator(this);
    TObject *obj;
    while ((obj = dynamic_cast<TObject*>(iterator())))
    {
      if (valInput.Index("@")<0) break;

      TString className = obj -> ClassName();
      TString parName = obj -> GetName();
      TString parValue;

      if (className == "TNamed") {
        parValue = ((TNamed *) obj) -> GetTitle();
        if (parValue.Index("AXIS_PARAMETER_")==0) parValue.ReplaceAll("AXIS_PARAMETER_","");
        else if (parName.Index("INPUT_FILE_")==0) parName.ReplaceAll("INPUT_FILE_","");
      }
      else if (className == "TParameter<double>") parValue += ((TParameter<double> *) obj) -> GetVal();
      else if (className == "TParameter<bool>")   parValue = ( (((TParameter<bool> *) obj) -> GetVal()==true) ? "true" : "false");
      else if (className == "TParameter<int>") {
        parValue += ((TParameter<int> *) obj) -> GetVal();
        if (parName.Index("NUM_VALUES_")==0) continue;
        if (parName.Index("VECTOR3_")==0) {
          parName.ReplaceAll("VECTOR3_","");
          auto numV = parValue.Atoi();
          parValue = "";
          for (auto iVal=0; iVal<numV; ++iVal) {
            parValue += GetParDouble(parName,iVal);
            parValue += " ";
          }
        }
      }
      else
        continue;

           if (ita==0) parName = Form("@%s@",parName.Data());
      else if (ita==1) parName = Form("@%s" ,parName.Data());
      valInput.ReplaceAll(parName,parValue);
    }
  }

  return valInput;
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
    kr_error(0) << "Parameter file " << fileNameFull << " does not exist!" << endl;
    return 0;
  }
  kr_info(0) << "Adding parameter file " << fileNameFull << endl;

  if (parNameForFile.IsNull())                   parNameForFile = Form("INPUT_FILE_%d", fNumInputFiles);
  else if (parNameForFile.Index("INPUT_FILE")<0) parNameForFile = Form("INPUT_FILE_%d", fNumInputFiles);

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
  kr_info(0) << "Adding parameter container " << parc -> GetName() << endl;

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
        kr_error(0) << "Parameter with name " << name << " already exist!" << endl;
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
    if (fDebugMode) kr_info(0) << "[" << fName << "]" << "(debug mode)" << " List of parameters :" << endl;
    else kr_info(0) << "[" << fName << "]" << " List of parameters :" << endl;
  }
  if (!printout) {
    kr_info(0) << "Writting " << fileName << " as parameter file." << endl;
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
      else if (key.Index("INPUT_FILE_")==0) {
        //valueString.ReplaceAll("INPUT_FILE_","");
        parType = "file";
      }
      else
        parType = "string";
    }
    else if (className == "TParameter<int>") {
      TParameter<int> *par = (TParameter<int> *) obj;
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
      TParameter<double> *par = (TParameter<double> *) obj;
      valueString += par -> GetVal();
      parType = "double";
    }
    else if (className == "TParameter<bool>") {
      TParameter<bool> *par = (TParameter<bool> *) obj;
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
        if (thisIsNewLine) kr_info(0) << ssLine.str();
        else kb_out << ssLine.str();
      else
        fileOut << ssLine.str();
    }
  }
}


Bool_t KBParameterContainer::SetPar(std::string line)
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
    val = ReplaceVariables(val);
    AddFile(val, parName);
  }
  else if (parType == "v3" || parType == "vector3" || parType == "tvector3" || parType == "kbvector3") {
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("3")+2, valuesInString.Sizeof()-valuesInString.First("3")));
    valuesInString = ReplaceVariables(valuesInString);
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
    TString valSingle;
    ss >> valSingle;
    valSingle = ReplaceVariables(valSingle);
    valSingle.ToLower();
    Bool_t val = false;
    if (valSingle == "true" || valSingle == "1" || valSingle == "ktrue")
      val = true;
    SetPar(parName, val, overwrite);
  }
  else if (parType.Index("b[")==0 || parType.Index("bool[")>=0 || parType.Index("bool_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    valuesInString = ReplaceVariables(valuesInString);
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
      TString valSingle = TString(((TObjString *) valueTokens->At(idx))->GetString());
      Bool_t val = false;
      if (valSingle == "true" || valSingle == "1" || valSingle == "ktrue")
        val = true;
      SetPar(parName+"["+iVal+"]",val,overwrite);
    }
  }
  else if (parType == "i" || parType == "int" || parType == "int_t" || parType == "w" || parType == "width" || parType == "width_t") {
    TString valSingle;
    ss >> valSingle;
    valSingle = ReplaceVariables(valSingle);
    if (valSingle.Index("k")==0)
      SetParColor(parName, valSingle, overwrite);
    else {
      Int_t valInt = Int_t(TFormula("formula",valSingle).Eval(0));
      SetPar(parName, valInt, overwrite);
    }
  }
  else if (parType.Index("i[")==0 || parType.Index("int[")>=0 || parType.Index("int_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    if (valuesInString.First("#")>=0)
      valuesInString = TString(valuesInString(0, valuesInString.First("#")));
    valuesInString = ReplaceVariables(valuesInString);
    valuesInString.ReplaceAll("{",""); valuesInString.ReplaceAll("}",""); valuesInString.ReplaceAll(","," ");
    auto valueTokens = valuesInString.Tokenize(" ");
    Int_t numValues = valueTokens -> GetEntriesFast();
    if (arrayLength!=0 && numValues > arrayLength) numValues = arrayLength;
    if (arrayLength>1 && numValues==1) { numValues = arrayLength; arrayLength = -1; }
    SetPar(TString("NUM_VALUES_")+parName,numValues,overwrite);
    for (auto iVal=0; iVal<numValues; ++iVal) {
      int idx = ((arrayLength<0) ? 0 : iVal);
      TString valString = TString(((TObjString *) valueTokens->At(idx))->GetString());
      if (valString.Index("k")==0)
        SetParColor(parName+"["+iVal+"]",valString,overwrite);
      else {
        Int_t valInt = TFormula("formula",valString).Eval(0);
        SetPar(parName+"["+iVal+"]",valInt,overwrite);
      }
    }
  }
  else if (parType == "d" || parType == "double" || parType == "double_t" || parType == "size" || parType == "size_t") {
    TString valSingle;
    ss >> valSingle;
    valSingle = ReplaceVariables(valSingle);
    Double_t val = TFormula("formula",valSingle).Eval(0);
    SetPar(parName, val, overwrite);
  }
  else if (parType.Index("d[")==0 || parType.Index("double[")>=0 || parType.Index("double_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    valuesInString = ReplaceVariables(valuesInString);
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
      TString valSingle = TString(((TObjString *) valueTokens->At(idx))->GetString());
      if (valSingle.Index("k")==0)
        SetParColor(parName+"["+iVal+"]",valSingle,overwrite);
      else {
        Double_t val = TFormula("formula",valSingle).Eval(0);
        SetPar(parName+"["+iVal+"]",val,overwrite);
      }
    }
  }
  else if (parType == "s" || parType == "string" || parType == "tstring") {
    TString valSingle;
    ss >> valSingle;
    valSingle = ReplaceVariables(valSingle);
    SetPar(parName, valSingle, overwrite);
  }
  else if (parType.Index("s[")==0 || parType.Index("string[")>=0 || parType.Index("tstring[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    valuesInString = ReplaceVariables(valuesInString);
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
      SetPar(parName+"["+iVal+"]",val,overwrite);
    }
  }
  else if (parType == "a" || parType == "axis" || parType == "kbvector3::axis") {
    TString valSingle;
    ss >> valSingle;
    valSingle = ReplaceVariables(valSingle);
    if (valSingle.Index("AXIS_PARAMETER_")<0)
      valSingle = TString("AXIS_PARAMETER_") + valSingle;
    SetPar(parName, valSingle, overwrite);
  }
  else if (parType.Index("a[")==0 || parType.Index("axis[")>=0 || parType.Index("kbvector3::axis[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    valuesInString = ReplaceVariables(valuesInString);
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
    TString valSingle;
    ss >> valSingle;
    valSingle = ReplaceVariables(valSingle);
    SetParColor(parName, valSingle, overwrite);
  }
  else if (parType.Index("c[")==0 || parType.Index("color[")>=0 || parType.Index("color_t[")==0) {
    Int_t arrayLength = TString(parType(parType.First("[")+1,parType.First("]")-parType.First("[")-1)).Atoi();
    TString valuesInString = line;
    valuesInString = TString(valuesInString(valuesInString.First("]")+2, valuesInString.Sizeof()-valuesInString.First("]")));
    valuesInString = ReplaceVariables(valuesInString);
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
      kr_error(0) << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Add(new TParameter<bool>(name, val));
  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, Int_t val, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      kr_error(0) << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Add(new TParameter<int>(name, val));
  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, Double_t val, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      kr_error(0) << "Parameter with name " << name << " already exist!" << endl;
      return false;
    }
  }

  Add(new TParameter<double>(name, val));
  return true;
}

Bool_t KBParameterContainer::SetPar(TString name, TString val, Bool_t overwrite)
{
  if (FindObject(name) != nullptr) {
    if (overwrite)
      this -> Remove(FindObject(name));
    else {
      kr_error(0) << "Parameter with name " << name << " already exist!" << endl;
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
      kr_error(0) << "Parameter with name " << name << " already exist!" << endl;
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

    /*
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
     */

    val = Int_t(TFormula("formula",valColor).Eval(0));
  }
  else
    return false;

  Add(new TParameter<int>(name, val));
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
  if (fAllowConversionFromString&&className=="TNamed") {
    TString title = ((TNamed *) obj) -> GetTitle();
    return title.Atoi();
  }
  if (className != "TParameter<int>") {
    kr_error(0) << name << " parameter type is " << className << ", not int!" << endl;
  }

  return ((TParameter<int> *) obj) -> GetVal();
}

Bool_t KBParameterContainer::GetParBool(TString name, Int_t idx)
{
  if (idx>=0) return GetParBool(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    kr_error(0) << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, false);
    return false;
  }

  TString className = obj -> ClassName();
  if (fAllowConversionFromString&&className=="TNamed") {
    TString title = ((TNamed *) obj) -> GetTitle();
    title.ToLower();
    if (title.Index("true")>=0) return true;
    else if (title.Index("false")>=0) return false;
    else if (title.IsDec()) return Bool_t(title.Atoi());
    else return false;
  }
  if (className != "TParameter<bool>") {
    kr_error(0) << name << " parameter type is " << className << ", not bool!" << endl;
  }

  return ((TParameter<bool> *) obj) -> GetVal();
}

Int_t KBParameterContainer::GetParInt(TString name, Int_t idx)
{
  if (idx>=0) return GetParInt(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    kr_error(0) << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, -999);
    return -999;
  }

  TString className = obj -> ClassName();
  if (fAllowConversionFromString&&className=="TNamed") {
    TString title = ((TNamed *) obj) -> GetTitle();
    return title.Atoi();
  }
  if (className != "TParameter<int>") {
    kr_error(0) << name << " parameter type is " << className << ", not !" << endl;
  }

  return ((TParameter<int> *) obj) -> GetVal();
}

Double_t KBParameterContainer::GetParDouble(TString name, Int_t idx)
{
  if (idx>=0) return GetParDouble(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    kr_error(0) << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, -999.999);
    return -999.999;
  }

  TString className = obj -> ClassName();
  if (fAllowConversionFromString&&className=="TNamed") {
    TString title = ((TNamed *) obj) -> GetTitle();
    return title.Atof();
  }
  if (className != "TParameter<double>") {
    kr_error(0) << name << " parameter type is " << className << ", not double!" << endl;
  }

  return ((TParameter<double> *) obj) -> GetVal();
}

TString KBParameterContainer::GetParString(TString name, Int_t idx)
{
  if (idx>=0) return GetParString(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    kr_error(0) << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, "DOES_NOT_EXIST");
    return "DOES_NOT_EXIST";
  }

  TString className = obj -> ClassName();
  if (className != "TNamed") {
    kr_error(0) << name << " parameter type is " << className << ", not string!" << endl;
  }

  return ((TNamed *) obj) -> GetTitle();
}

#ifdef KEBI
kbaxis KBParameterContainer::GetParAxis(TString name, Int_t idx)
{
  if (idx>=0) return GetParAxis(name+"["+idx+"]");

  TObject *obj = FindObject(name);
  if (obj == nullptr) {
    kr_error(0) << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, TString("AXIS_PARAMETER_DOES_NOT_EXIST"));
    return KBVector3::kNon;
  }

  TString className = obj -> ClassName();
  if (className != "TNamed") {
    kr_error(0) << name << " parameter type is " << className << ", not axis(string)!" << endl;
  }

  TString value = ((TNamed *) obj) -> GetTitle();
  if (value.Index("AXIS_PARAMETER_")==0)
    value.ReplaceAll("AXIS_PARAMETER_","");
  return KBVector3::GetAxis(value);
}
#endif

vector<bool> KBParameterContainer::GetParVBool(TString name)
{
  vector<bool> array;
  auto npar = GetParN(name);
  for (auto i=0; i<npar; ++i)
    array.push_back(GetParBool(name,i));
  return array;
}

vector<int> KBParameterContainer::GetParVInt(TString name)
{
  vector<int> array;
  auto npar = GetParN(name);
  for (auto i=0; i<npar; ++i)
    array.push_back(GetParInt(name,i));
  return array;
}

vector<double> KBParameterContainer::GetParVDouble(TString name)
{
  vector<double> array;
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

#ifdef KEBI
vector<kbaxis> KBParameterContainer::GetParVAxis(TString name)
{
  vector<kbaxis> array;
  auto npar = GetParN(name);
  for (auto i=0; i<npar; ++i)
    array.push_back(GetParAxis(name,i));
  return array;
}
#endif



TVector3 KBParameterContainer::GetParV3(TString name)
{
  TString xname = name + "[0]";
  TString yname = name + "[1]";
  TString zname = name + "[2]";

  TObject *xobj = FindObject(xname);
  if (xobj == nullptr) {
    kr_error(0) << "parameter with name " << name << " does not exist!" << endl;
    if (!fDebugMode) gApplication -> Terminate();
    SetPar(TString("NEWPAR")+name, -999.999);
    return TVector3(-999,-999,-999);
  }

  TString className = xobj -> ClassName();
  if (className != "TParameter<double>") {
    kr_error(0) << name << " parameter type is " << className << ", not (v3)double!" << endl;
  }

  auto x = ((TParameter<double> *) FindObject(xname)) -> GetVal();
  auto y = ((TParameter<double> *) FindObject(yname)) -> GetVal();
  auto z = ((TParameter<double> *) FindObject(zname)) -> GetVal();

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
