#ifndef KBPARAMETERCONTAINER
#define KBPARAMETERCONTAINER

#include "TObjArray.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TFormula.h"
#include <vector>
using namespace std;
#ifdef KEBI
#include "KBVector3.hh"
#include "KBGlobal.hh"
#endif
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>

#ifdef KEBI
typedef KBVector3::Axis kbaxis;
#else
typedef TVector3 kbaxis;
#endif

/**
 * List of parameters with various types
 *
 * Structure of parameter file should be list of : [name] [type] [value]
 *
 * [name] is the name of parameter and [value] is the value of parameter. [type] is one of below :
 *
 * @param f   (file) for parameter file of KBParameterContainer 
 * @param b   (bool, Bool_t) : Bool_t, Features 1, 0, true, false, kTrue, kFalse
 * @param i   (int, Int_t) Int_t
 * @param d   (double, Double_t) Double_t
 * @param s   (string, TString) TString (do not accept space)
 * @param a   (axis, KBVector3::Axis) KBVector3::Axis
 * @param v3  (vector3, TVector3, KBVector3) KBVector3::Axis
 * @param c   (color, Color_t) Color_t saved as Int_t. Features TColor keyword such as kRed+1.
 * @param w   (width, Width_t) Width_t saved as Int_t.
 * @param s   (size, Size_t) Size_t saved as Double_t.
 *
 * Comments starts with #.
 *
 * ex)\n
 * \#example parameter file\n
 * worldSize    d  1000
 * specialFile  s  /path/to/specialFile.dat
 *
 * Generally, the parameters are stored with variable type.
 * However, parameters can be picked up from the string parameter with corresponding return type of the getter, once the method SetAllowConversionFromString() is called.
 * In this case, One should be carefull that the parameter value is safe to be converted to the chosen return type.
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

    void SetAllowConversionFromString(bool val=true) { fAllowConversionFromString = val; }

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
#ifdef KEBI
    kbaxis   GetParAxis  (TString name, Int_t idx=-1);  ///< Get KBVector3::Axis type parameter by given name.
#endif

    vector<bool>   GetParVBool  (TString name);
    vector<int>    GetParVInt   (TString name);
    vector<double> GetParVDouble(TString name);
    vector<TString>  GetParVString(TString name);
#ifdef KEBI
    vector<kbaxis>   GetParVAxis  (TString name);
#endif

    TVector3 GetParV3(TString name);
    Double_t GetParX(TString name) { return GetParDouble(name,0); }
    Double_t GetParY(TString name) { return GetParDouble(name,1); }
    Double_t GetParZ(TString name) { return GetParDouble(name,2); }

    Int_t    GetParWidth(TString name, int idx=-1) { return GetParInt(name,idx); }
    Int_t    GetParColor(TString name, int idx=-1) { return GetParInt(name,idx); }
    Double_t GetParSize (TString name, int idx=-1) { return GetParDouble(name,idx); }

    vector<int>    GetParVWidth(TString name) { return GetParVInt(name); }
    vector<int>    GetParVColor(TString name) { return GetParVInt(name); }
    vector<double> GetParVSize (TString name) { return GetParVDouble(name); }

    Bool_t CheckPar(TString name);

    void ReplaceEnvironmentVariable(TString &val);
    TString ReplaceVariables(TString &val);

    bool IsEmpty() const;

  private:
    Bool_t fDebugMode = false;
    Bool_t fAllowConversionFromString = false;
    Int_t fNumInputFiles = 0;

  ClassDef(KBParameterContainer, 1)
};

#endif
