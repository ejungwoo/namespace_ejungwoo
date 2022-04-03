#ifndef ejungwo

#include "TClonesArray.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TObject.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

namespace ejungwoo
{
  class binning; ///< 1-dimensional binning class
  class binning2; ///< 2-dimensional binning class
  typedef KBParameterContainer parContainer; ///< parameter container (for configuration) class.
  class TParCanvas {
    public:
      int textFont=133, textSize=0, menuSize=25, padSize[2], marginInn[4], marginTop[4], marginDiv[2];
      int axisNdivision[3]={0}, titleAlign[4], titleColor[4], numContours=200, sidePad[2]={0};
      bool removeInnerPadAxis[2], removeInnerMainTitle, removeInnerZaxis, titleRotate[4], setStats=false, setDivIndexYX=false;
      double tickSize[3], labelSize[3], labelOffset[3], titleSize[4], titleOffset[4];
  };
  class TParAttribute {
    public:
      vector<int> lineStyle, lineWidth, lineColor, markerStyle, markerColor, fillColor, fillStyle, textColor, allColor;
      vector<double> markerSize;
      int textFont=133, textSize=25, textAlign=22;
  };

  enum verboseLevel : int { kQuiet=0, kNormal=1, kAll=2 };
  verboseLevel fVerboseLevel = kNormal;
  void verbose(int level);

  TString fNameVersion;
  int fHistIndex = 0;
  int fCanvasIndex = 0;
  TObjArray *fCanvasArray = nullptr;
  TClonesArray *fParameterArray = nullptr;
  const char *fNameAttributeConf = "att1";
  const char *fNameCanvasConf = "single";
  TParAttribute fParAtt;
  TParCanvas fParCvs;
  void setNextHistIndex(int val) { fHistIndex = val; }
  parContainer *conf(const char *nameConf);
  bool findConf(TVirtualPad *vpad, TString &nameConfCvs, int &nx, int &ny, int &cidx);
  TParCanvas getParCanvas(TString name);
  TParAttribute getParAttribute(TString name);
  TParCanvas getParCanvas(parContainer *par);
  TParAttribute getParAttribute(parContainer *par);
  TObject *att(TObject *obj, int idx=0, const char *nameConf="");

  void        findxy(TVirtualPad *pad, double &x1, double &y1);
  void        findxy(TVirtualPad *pad, double &x1, double &y1, double &dx, double &dy, double &xUnit, double &yUnit);
  //TCanvas*    canvas(const char *nameCvs, int dx, int dy, TArrayD poss);
  TCanvas*    canvas(const char *nameCvs, int nx, int ny, const char *nameConf="");
  TCanvas*    canvas(const char *nameCvs, int nn, const char *nameConf="");
  TCanvas*    canvas(const char *nameCvs="", const char *nameConf="") { return canvas(nameCvs, 1, 1, nameConf); }
  int         icanvas() { return fCanvasIndex; }
  TPad*       innerpad(TVirtualPad *vpad, double xr1, double yr1, double xr2, double yr2);
  TH1*        make(TH1* hist, TVirtualPad *vpad, TString drawOption="");
  TH1*        draw(TH1* hist, TVirtualPad *vpad=nullptr, TString drawOption="", TString legendContent="");
  TGraph*     draw(TGraph *graph, TVirtualPad *vpad=nullptr, TString drawOption="", TString legendContent="");
  TF1*        draw(TF1 *func, TVirtualPad *vpad=nullptr, TString drawOption="", TString legendContent="");
  TPad*       redraw(TPad *pad);
  TLegend*    getlg(TVirtualPad *vpad);
  TGaxis*     drawz(TH1* hist, TVirtualPad *vpad, const char *titlez="", bool logaxis=false);
  void        drawef(TVirtualPad *vpad);
  TLegend*    make(TLegend* lg,TVirtualPad *vp,double x1=-1,double y1=-1,double dx=0,double dy=0,double mgo=-1,int tf=-1,double tz=-1,int ta=-1,int fs=-1,int fc=-1,int bs=-1);
  TLegend*    draw(TLegend* lg,TVirtualPad *vp,double x1=-1,double y1=-1,double dx=0,double dy=0,double mgo=-1,int tf=-1,double tz=-1,int ta=-1,int fs=-1,int fc=-1,int bs=-1);
  TPad*       newpad(TVirtualPad *vpad, TString name, TString title, double x1, double y1, double dx, double dy);
  TPaveText*  newpt(TString content,TVirtualPad *vpad,double x1=-1,double y1=-1,double dx=0,double dy=0,int tf=133,double tz=20,int ta=12,int fs=0,int fc=0,int bs=0);
  TLatex*     newtt(TString content, double x, double y, int tf=133, double tz=20, int ta=12, int tc=kBlack);
  TH1*        project(TTree *tree, TString formula, TCut cut="", TString name="", TString title="", int nx=0, double x1=0, double x2=0, int ny=0, double y1=0, double y2=0);
  TH1*        project(TTree *tree, TString formula, int nx, double x1=0, double x2=0, int ny=0, double y1=0, double y2=0, TCut cut="", TString name="", TString title="")
              { return project(tree, formula, cut, name, title, nx, x1, x2, ny, y1, y2); }

  TGraphErrors *tograph(TH1D *,TString opt="ex;ey",int iatt=-1,TString nameConf="");
  double get_x(TGraph *graph,int i);
  double get_y(TGraph *graph,int i);
  TString tok(TString line, TString token, int i);
  TString tok(TObjArray *line, int i);
  TF1 *fitg(TH1 *hist, double distSigma=2.5, Option_t *opt="RQ0");
  TGraph *extractPoints(TGraph *graph_hand_pointed, double x1, double x2, double y1, double y2, TString saveName="");

  TString get_value(TString array, TString skey, bool exact=false, TString token=";");
  TString get_value(TObjArray *elements, TString skey, bool exact=false);
  bool exist_key   (TObjArray *elements, TString skey);
  bool exist_keyv  (TObjArray *elements, TString skey);
  int get_valuei   (TObjArray *elements, TString skey);
  double get_valued(TObjArray *elements, TString skey);

  void colorwheel();
  void markers();
  TString toString(double value, int ndigit=-1);
  TDatabasePDG *particleDB();
  TParticlePDG *particle(TString name);
  TParticlePDG *particle(int pdg);
  void version(TString nameVersion) { fNameVersion = nameVersion; }
  void saveRoot(TObject *obj, TString nameFile="", TString nameVersion="", bool savePrimitives=false, bool simplifyNames=false);
  void savePDF(TCanvas *cvs, TString nameVersion="");
  void savePNG(TCanvas *cvs, TString nameVersion="");
  void saveEPS(TCanvas *cvs, TString nameVersion="");
  void saveRoot(TString nameVersion="", bool savePrimitives=false, bool simplifyNames=false);
  void savePDF(TString nameVersion="");
  void savePNG(TString nameVersion="");
  void saveEPS(TString nameVersion="");
  void saveAll(TString nameVersion="");
  void write(TObject *obj);
};

class ejungwoo::binning
{
  private:
    int fNbins = 0; ///< number of bins
    double fMin = 0; ///< lower bound
    double fMax = 0; ///< upper bound
    double fWidth = 0; ///< binning space width
    double fValue = 0; ///< value will be set after iteration using next(), back(), nextb()
    int fBindex = 0; ///< index for iteration
    const char *fTitle = "";
    const char *fSelection = "";
    TString fExpression;
    TString fFormula;
    vector<TString> fParameters;

    enum bntype : int { kEqualSpacing=0, kContinuousSpacing=1, kDescreteSpacing=2 };
    bntype fSpacing = kEqualSpacing;
    TArrayD fBinArray;

  public:
    binning() : fNbins(0), fMin(0), fMax(0) {}
    binning(int nb, double mi, double mx, const char *tl="", const char *fml="", const char *sel="") : fNbins(nb), fMin(mi), fMax(mx), fTitle(tl), fSelection(sel) { setFormula(fml); init(); }
    binning(int nb, TArrayD arr, const char *tl="", const char *fml="", const char *sel="") : fTitle(tl), fSelection(sel) { setFormula(fml); initArray(nb,arr); }
    binning(TArrayD arr, const char *ttl="", const char *fml="", const char *sel="") : binning(0, arr, ttl, fml, sel) {}
    binning(TTree *tree, const char *branchName) : binning() { make(tree, branchName); }
    binning(TH1 *hist, int axis_123=1) : binning() { make(hist, axis_123); }
    ~binning() {}

    void init();
    void initArray(int nbins, TArrayD array);
    void make(TH1 *hist, int axis_123=1);
    void make(TTree *tree, const char* branchName);
    void set(double nbins, double min, double max, const char *ttl="", const char *fml="", const char *sel="");
    void set(int nbins, TArrayD array, const char *ttl="", const char *fml="", const char *sel="");
    void set(TArrayD array, const char *ttl="", const char *fml="", const char *sel="");
    void setN         (double nbins)     { fNbins = nbins; fWidth = (fMax-fMin)/fNbins; }
    void setMin       (double min)       { fMin = min; fWidth = (fMax-fMin)/fNbins; }
    void setMax       (double max)       { fMax = max; fWidth = (fMax-fMin)/fNbins; }
    void setW         (double w)         { fWidth = w; fNbins = int((fMax-fMin)/fWidth); }
    void setTitle     (const char *ttl)  { fTitle = ttl; }
    void setSelection (const char *sel)  { fSelection = sel; }
    void setExpression(const char *fml)  { setFormula(fml); }
    void setFormula   (const char *fml)  {
      fFormula = fml;
      fExpression = fml;
      fParameters.clear();
      int posf = 0;
      while (1) {
        auto posi = fFormula.Index("[",1,posf,TString::kExact); if (posi<0) break;
        auto pose = fFormula.Index("=",1,posi,TString::kExact);
             posf = fFormula.Index("]",1,posi,TString::kExact);
        if (pose>=0&&pose<posf) fParameters.push_back(TString(fFormula(posi+1,pose-posi-1)));
        else fParameters.push_back(TString(fFormula(posi+1,posf-posi-1)));
      }
    }

    bntype getBinningType()       const { return fSpacing; }
    bool   esp()                  const { if (fSpacing==kEqualSpacing) return true; return false; }
    bool   csp()                  const { if (fSpacing==kContinuousSpacing) return true; return false; }
    bool   dsp()                  const { if (fSpacing==kDescreteSpacing) return true; return false; }
    bool   isNull()               const { if (fNbins<1||(fMin==0&&fMax==0)) return true; return false; }
    int    getN()                 const { return fNbins; }
    double getMin()               const { return fMin; }
    double getMax()               const { return fMax; }
    double getW(int bin=0)        const;
    const char *getTitle()        const { return fTitle; }
    TString getFormula()          const { return fFormula; }
    const char *getExpression()   const { return fExpression.Data(); }
    const char *getSelection()    const { return fSelection; }

    TArrayD getBinArray()         const { return fBinArray; }
    bool   null()                 const { return isNull(); }
    int    n()                    const { return getN(); }
    double min()                  const { return getMin(); }
    double max()                  const { return getMax(); }
    double width(int bin=-1)      const { return getW(bin); }
    const char *title()           const { return getTitle(); }
    const char *expr()            const { return getExpression(); }
    TString formula()             const { return getFormula(); }
    const char *selection()       const { return getSelection(); }

    void setPar(const char *skey, double value) {
      int posi=-1,pose,posf=0;
      while (1) {
        posi = fFormula.Index(Form("[%s",skey),posi+1,TString::kExact); if (posi<0) break;
        posf = fFormula.Index("]",1,posi,TString::kExact);
        fFormula.Replace(posi,posf-posi+1,Form("[%s=%f]",skey,value));
      }
      fExpression = fFormula;
      posf = 0;
      for (auto par : fParameters) {
        posi = fExpression.Index(Form("[%s",par.Data()),posi+1,TString::kExact);
        posf = fExpression.Index("]",1,posi,TString::kExact);
        pose = fExpression.Index("=",1,posi,TString::kExact);
        fExpression.Replace(posi,posf-posi+1,TString(fExpression(pose+1,posf-pose-1)));
      }
    }

    void setPar(const char *skey, int value) {
      int posi=-1,pose,posf=0;
      while (1) {
        posi = fFormula.Index(Form("[%s",skey),posi+1,TString::kExact); if (posi<0) break;
        posf = fFormula.Index("]",1,posi,TString::kExact);
        fFormula.Replace(posi,posf-posi+1,Form("[%s=%d]",skey,value));
      }
      fExpression = fFormula;
      posf = 0;
      for (auto par : fParameters) {
        posi = fExpression.Index(Form("[%s",par.Data()),posi+1,TString::kExact);
        posf = fExpression.Index("]",1,posi,TString::kExact);
        pose = fExpression.Index("=",1,posi,TString::kExact);
        fExpression.Replace(posi,posf-posi+1,TString(fExpression(pose+1,posf-pose-1)));
      }
    }

    void setPar(int ikey, double value) { return setPar(Form("%d",ikey),value); }
    void setPar(int ikey, int value)    { return setPar(Form("%d",ikey),value); }

    //iterator
    void   reset();
    bool   next();
    void   end();
    bool   prev();
    bool   isf()    const { return ((fBindex==1)?true:false); }
    bool   isl()    const { return ((fBindex==fNbins)?true:false); }
    int    ii()     const { return fBindex-1; }
    int    bi()     const { return fBindex; }
    int    hi()     const { if (dsp()) return (2*fBindex-1); return fBindex; }
    double val()    const { return fValue; }
    double low()    const { return lowEdge(hi()); }
    double high()   const { return highEdge(hi()); }
    double center() const { return val(); }
    void   toi(int idx) { reset(); for (auto i=0; i<idx; ++i) next(); }

    double lowEdge(int bin)       const;
    double highEdge(int bin)      const;
    int    findBin(double val)    const;
    double getFullWidth()         const { return (fMax - fMin); }
    double getFullCenter()        const { return (fMax + fMin)/2.; }
    double getCenter(int bin)     const { return (fMin + (bin-.5)*fWidth); }
    bool   isInside(double value) const { if(value>fMin && value<fMax) return true; return false; }
    double xByRatio(double ratio) const { return ((1.-ratio)*fMin + ratio*fMax); }
    double wByRatio(double ratio) const { return ratio*(fMax-fMin); }
    double fullw()                const { return getFullWidth(); }
    double fullc()                const { return getFullCenter(); }
    double center(int bin)        const { return getCenter(bin); }
    bool   inside(double value)   const { return isInside(value); }
    double xratio(double ratio)   const { return xByRatio(ratio); }
    double wratio(double ratio)   const { return wByRatio(ratio); }

    const char *xxmm(int bin=0, int ndigit=-1) const {
      const char *namemm;
      if (bin==0)     namemm = Form("%s_%s", toString(fMin,ndigit).Data(),         toString(fMax,ndigit).Data()          );
      else if (bin<0) namemm = Form("%s_%s", toString(fMin,ndigit).Data(),         toString(highEdge(-bin),ndigit).Data());
      else            namemm = Form("%s_%s", toString(lowEdge(bin),ndigit).Data(), toString(highEdge(bin),ndigit).Data() );
      return namemm;
    }

    const char *nnmm(int bin=0, int ndigit=-1) const {
      const char *namemm;
      if (bin==0)     namemm = Form("%s_%d_%s_%s",getExpression(),fNbins, toString(fMin,ndigit).Data(),         toString(fMax,ndigit).Data()          );
      else if (bin<0) namemm = Form("%s_%d_%s_%s",getExpression(),fNbins, toString(fMin,ndigit).Data(),         toString(highEdge(-bin),ndigit).Data());
      else            namemm = Form("%s_%d_%s_%s",getExpression(),fNbins, toString(lowEdge(bin),ndigit).Data(), toString(highEdge(bin),ndigit).Data() );
      return namemm;
    }

    const char *ttmm(int bin=0, int ndigit=-1) const {
      const char *titlemm;
      if (bin==0)     titlemm = Form("%s = %s - %s",fTitle,toString(fMin,ndigit).Data(),         toString(fMax,ndigit).Data());
      else if (bin<0) titlemm = Form("%s = %s - %s",fTitle,toString(fMin,ndigit).Data(),         toString(highEdge(-bin),ndigit).Data());
      else            titlemm = Form("%s = %s - %s",fTitle,toString(lowEdge(bin),ndigit).Data(), toString(highEdge(bin),ndigit).Data());
      return titlemm;
    }

    TCut cut(int bin=0, int ndigit=-1) const {
      TCut cut0;
      if (bin==0)     cut0 = TCut(xxmm(bin,ndigit),Form("%s>=%f&&%s<%f",getExpression(),fMin,getExpression(),fMax));
      else if (bin<0) cut0 = TCut(xxmm(bin,ndigit),Form("%s>=%f&&%s<%f",getExpression(),fMin,getExpression(),highEdge(-bin)));
      else            cut0 = TCut(xxmm(bin,ndigit),Form("%s>=%f&&%s<%f",getExpression(),lowEdge(bin),getExpression(),highEdge(bin)));
      return cut0;
    }

    TLine *newly(double yval, int lc=1, int ls=1) const { auto line = new TLine(fMin,yval,fMax,yval); line -> SetLineColor(lc); line -> SetLineStyle(ls); return line; }
    TLine *newlx(double xval, int lc=1, int ls=1) const { auto line = new TLine(xval,fMin,xval,fMax); line -> SetLineColor(lc); line -> SetLineStyle(ls); return line; }

    TH1D *newHist(TString nameHist="", TString title="");
    TH1D *project(TTree *tree, TCut selection="", TString name="", TString title="");

    TString print(bool pout=1) const;
    TString printi(bool pout=1) const;
    binning copy() const;
    void operator=(const binning binn);
    ejungwoo::binning2 mult(const binning *binn);
    ejungwoo::binning2 operator*(const binning binn);
};

void ejungwoo::binning::init() {
  fSpacing = kEqualSpacing;
  if (fWidth>0&&fNbins<1) setW(fWidth);
  else if (fNbins>0&&fWidth<1) setN(fNbins);
}

void ejungwoo::binning::initArray(int nbins, TArrayD array) {
  fWidth = -1;
  fNbins = nbins;
  if (fNbins==0) fNbins = array.GetSize()-1;
  fSpacing = (fNbins==(array.GetSize()-1))?kContinuousSpacing:kDescreteSpacing;
  int nValues = (fSpacing==kContinuousSpacing)?nbins+1:2*nbins;
  fBinArray = array;
  fMin = fBinArray[0];
  fMax = fBinArray[nValues-1];
}

void ejungwoo::binning::make(TH1 *hist, int axis_123) {
  TAxis *axis;
  if (axis_123==2) axis = hist -> GetYaxis();
  else if (axis_123==3) axis = hist -> GetZaxis();
  else axis = hist -> GetXaxis();

  fTitle = axis -> GetTitle();
  TArrayD xbins(axis->GetXbins()->GetSize(),axis->GetXbins()->GetArray());
  if (xbins.GetSize()>0) {
    initArray(xbins.GetSize()-1,xbins);
  } else {
    fNbins = axis -> GetNbins();
    fMin = axis -> GetBinLowEdge(1);
    fMax = axis -> GetBinUpEdge(fNbins);
    fWidth = -1;
    init();
  }
}

void ejungwoo::binning::make(TTree *tree, const char* branchName) {
  if (fNbins<1) fNbins = 100;
  if (tree==nullptr) return;
  auto min = tree -> GetMinimum(branchName);
  auto max = tree -> GetMaximum(branchName);
  auto dmm = (max - min) / 100.;
  fMax = max + dmm;
  fMin = min;
  fTitle = branchName;
  fFormula = branchName;
  fExpression = branchName;
  init();
}

void ejungwoo::binning::set(double nbins, double min, double max, const char *ttl, const char *fml, const char *selection) {
  fNbins = nbins;
  fMin = min;
  fMax = max;
  fWidth = (fMax-fMin)/fNbins;
  init();
  if (strcmp(ttl,"")!=0) fTitle = ttl;
  if (strcmp(fml,"")!=0) { setFormula(fml); }
  if (strcmp(selection,"")!=0) fSelection = selection;
}

void ejungwoo::binning::set(int nbins, TArrayD array, const char *ttl, const char *fml, const char *selection) {
  initArray(nbins,array);
  if (strcmp(ttl,"")!=0) fTitle = ttl;
  if (strcmp(fml,"")!=0) { setFormula(fml); }
  if (strcmp(selection,"")!=0) fSelection = selection;
}

void ejungwoo::binning::set(TArrayD array, const char *ttl, const char *fml, const char *selection) {
  initArray(0,array);
  if (strcmp(ttl,"")!=0) fTitle = ttl;
  if (strcmp(fml,"")!=0) { setFormula(fml); }
  if (strcmp(selection,"")!=0) fSelection = selection;
}

double ejungwoo::binning::getW(int bin) const {
  if (bin==0) return (fMax - fMin);
  if (bin<0) bin = fBindex;
       if (esp()) { return fWidth; }
  else if (csp()) { return (fBinArray[bin] - fBinArray[bin-1]); }
  else if (dsp()) { return (fBinArray[bin*2-1] - fBinArray[bin*2-2]); }
  return fWidth;
}

void ejungwoo::binning::reset() { fBindex = 0; }
bool ejungwoo::binning::next()  {
  if (fBindex>fNbins-1) return false;
  fBindex++;
  if (esp()) fValue = fMin + (fBindex-1)*fWidth + .5*fWidth;
  else if (csp()) fValue = (.5*(fBinArray[fBindex]     + fBinArray[fBindex-1]));
  else if (dsp()) fValue = (.5*(fBinArray[fBindex*2-1] + fBinArray[fBindex*2-2]));
  return true;
}
void ejungwoo::binning::end()   { fBindex = fNbins+1; }
bool ejungwoo::binning::prev()  {
  if (fBindex<2) return false;
  fBindex--;
  if (esp()) fValue = fMin + (fBindex-1)*fWidth + .5*fWidth;
  else if (csp()) fValue = (.5*(fBinArray[fBindex]     + fBinArray[fBindex-1]));
  else if (dsp()) fValue = (.5*(fBinArray[fBindex*2-1] + fBinArray[fBindex*2-2]));
  return true;
}

double ejungwoo::binning::lowEdge(int bin) const {
  if (bin==1) return fMin;
  double edgeValue = -999;
       if (esp()) edgeValue = fMin+(bin-1)*(fMax-fMin)/fNbins;
  else if (csp()) edgeValue = fBinArray[bin-1];
  else if (dsp()) edgeValue = fBinArray[bin-1];
  return edgeValue;
}
double ejungwoo::binning::highEdge(int bin) const {
  if (bin==-1) return fMax;
  double edgeValue = -999;
       if (esp()) edgeValue = fMin+(bin)*(fMax-fMin)/fNbins;
  else if (csp()) edgeValue = fBinArray[bin];
  else if (dsp()) edgeValue = fBinArray[bin];
  return edgeValue;
}
int ejungwoo::binning::findBin(double val) const {
  int binf = 0;
       if (esp()) binf = int((val-fMin)/fWidth);
  else if (csp()) binf = 0;// TODO //binf = int((val-fMin)/fWidth);
  else if (dsp()) binf = 0;// TODO //binf = int((val-fMin)/fWidth);
  return binf = int((val-fMin)/fWidth);
}

TH1D *ejungwoo::binning::newHist(TString nameHist, TString title) {
  if (nameHist.IsNull())
    nameHist = Form("hist_%d",fHistIndex++);
  auto titlexy = Form("%s;%s;",title.Data(),fTitle);
  TH1D *hist;
  if (esp()) hist = new TH1D(nameHist,titlexy,fNbins,fMin,fMax);
  else if (csp()) hist = new TH1D(nameHist,titlexy,fNbins,fBinArray.GetArray());
  else if (dsp()) hist = new TH1D(nameHist,titlexy,2*fNbins-1,fBinArray.GetArray());
  return hist;
}

TH1D *ejungwoo::binning::project(TTree *tree, TCut selection, TString name, TString title) {
  if (name.IsNull()) {
    name = Form("hist_%d",fHistIndex++);
    //name.ReplaceAll("/","DV");
    //name.ReplaceAll("(","L");
    //name.ReplaceAll(")","R");
  }
  auto hist = newHist(name,title);
  auto nn = tree -> Project(name,getExpression(),selection);
  cout << name << "   " << getExpression() << "   " << selection << "   " << nn << endl;
  return hist;
}

TString ejungwoo::binning::printi(bool pout) const {
  TString line = Form("%s(%d)=%s~%s(%s)",getExpression(),ii(),toString(low()).Data(),toString(high()).Data(),toString(val()).Data());
  if (pout)
    cout << line << endl;
  return line;
}

TString ejungwoo::binning::print(bool pout) const {
  TString minString = ejungwoo::toString(fMin);
  TString maxString = ejungwoo::toString(fMax);
  TString wString = ejungwoo::toString(fWidth);
  TString line = Form("[%s,%s] %d,%s,%s",getExpression(),fTitle,fNbins,minString.Data(),maxString.Data());
  if (esp()) {}
  else if (csp()) {
    line = line + " : ";
    for (auto i=0; i<fNbins+1; ++i)
      line = line + fBinArray[i] + ", ";
  }
  else if (dsp()) {
    line = line + " : ";
    for (auto i=0; i<fNbins; ++i)
      line = line + "(" + fBinArray[2*i] + "," + fBinArray[2*i+1] + "), ";
  }
  if (pout)
    cout << line << endl;
  return line;
}

ejungwoo::binning ejungwoo::binning::copy() const {
  binning bncopy;
  if (esp()) bncopy = binning(fNbins, fMin, fMax, fTitle, getFormula(), fSelection);
  else       bncopy = binning(fNbins, fBinArray, fTitle, getFormula(), fSelection);
  return bncopy;
}

void ejungwoo::binning::operator=(const binning binn) {
  if (binn.esp())
    set(binn.getN(),binn.getMin(),binn.getMax(),binn.getTitle(),binn.getExpression(),binn.getSelection());
  else
    set(binn.getN(),binn.getBinArray(),binn.getTitle(),binn.getExpression(),binn.getSelection());
}

class ejungwoo::binning2
{
  private:
    binning fbnX;
    binning fbnY;

    const char *justTitle(const char *val) const {
      if (strcmp(val,"")==0)
        return "";
      auto elements = TString(val).Tokenize(";");
      const char *title = ((TObjString *) elements->At(0))->GetString();
      return title;
    }

  public :
    binning2(binning bnX=binning(), binning bnY=binning()) : fbnX(bnX), fbnY(bnY) {}
    binning2(int nx, double x1, double x2, int ny, double y1, double y2) : fbnX(binning(nx,x1,x2)), fbnY(binning(ny,y1,y2)) {}
    binning2(double x1, double x2, double y1, double y2) : fbnX(binning(100,x1,x2)), fbnY(binning(100,y1,y2)) {}
    binning2(TH2D *hist) { fbnX = binning(hist,1);  fbnY = binning(hist,2); }
    ~binning2() {}

    TH2D *newHist(TString nameHist="", TString titleHist="") {
      if (nameHist.IsNull())
        nameHist = Form("hist_%d",fHistIndex++);
      auto titlexy = Form("%s;%s;%s;",titleHist.Data(),justTitle(fbnX.getTitle()),justTitle(fbnY.getTitle()));
      TH2D *hist = nullptr;
           if ( fbnX.esp() && fbnY.esp()) hist = new TH2D(nameHist,titlexy, fbnX.getN(),fbnX.getMin(),fbnX.getMax(), fbnY.getN(),fbnY.getMin(),fbnY.getMax());
      else if ( fbnX.esp() && fbnY.csp()) hist = new TH2D(nameHist,titlexy, fbnX.getN(),fbnX.getMin(),fbnX.getMax(), fbnY.getN(),fbnY.getBinArray().GetArray());
      else if ( fbnX.csp() && fbnY.esp()) hist = new TH2D(nameHist,titlexy, fbnX.getN(),fbnX.getBinArray().GetArray(), fbnY.getN(),fbnY.getMin(),fbnY.getMax());
      else if ( fbnX.csp() && fbnY.csp()) hist = new TH2D(nameHist,titlexy, fbnX.getN(),fbnX.getBinArray().GetArray(), fbnY.getN(),fbnY.getBinArray().GetArray());
      else if ( fbnX.esp() && fbnY.dsp()) hist = new TH2D(nameHist,titlexy, fbnX.getN(),fbnX.getMin(),fbnX.getMax(), 2*fbnY.getN()-1,fbnY.getBinArray().GetArray());
      else if ( fbnX.dsp() && fbnY.esp()) hist = new TH2D(nameHist,titlexy, 2*fbnX.getN()-1,fbnX.getBinArray().GetArray(), fbnY.getN(),fbnY.getMin(),fbnY.getMax());
      else if ( fbnX.dsp() && fbnY.dsp()) hist = new TH2D(nameHist,titlexy, 2*fbnX.getN()-1,fbnX.getBinArray().GetArray(), 2*fbnY.getN()-1,fbnY.getBinArray().GetArray());
      return hist;
    }

    const char *getExpression() const { return Form("%s:%s",fbnY.getExpression(),fbnX.getExpression()); }
    const char *getFormula()    const { return Form("%s:%s",fbnY.getExpression(),fbnX.getExpression()); }
    const char *getSelection()  const { return fbnY.getSelection(); }
    const char *expr()          const { return getExpression(); }
    const char *formula()       const { return getFormula(); }
    const char *name()          const { return namexy(); }
    const char *namexy()        const { return Form("%s%s",fbnX.getExpression(),fbnY.getExpression()); }
    const char *namexynn()      const { return Form("%s%d%s%d",fbnX.getExpression(),fbnX.getN(),fbnY.getExpression(),fbnY.getN()); }
    const char *nnmm()          const { return Form("%s__%s",fbnX.nnmm(),fbnY.nnmm()); }

    binning bx() const { return fbnX; }
    binning by() const { return fbnY; }

    TH2D *project(TTree *tree, TCut selection="", TString name="", TString title="") {
      if (name.IsNull())
        name = getExpression();
      auto hist = newHist(name,title);
      tree -> Project(name,getExpression(),selection);
      return hist;
    }

    double content(TH2D *hist) { return hist -> GetBinContent(fbnX.hi(),fbnY.hi()); }

    //iterator
    void reset() { fbnX.reset(); fbnY.reset(); fbnY.next(); }
    bool next()  {
      if (fbnX.next()) return true;
      if (!fbnY.next()) return false;
      fbnX.reset();
      return next();
    }

    int    ii(bool xy)     const { return ((xy)?fbnY.ii():fbnX.ii());         }
    int    bi(bool xy)     const { return ((xy)?fbnY.bi():fbnX.bi());         }
    int    hi(bool xy)     const { return ((xy)?fbnY.hi():fbnX.hi());         }
    double val(bool xy)    const { return ((xy)?fbnY.val():fbnX.val());       }
    double low(bool xy)    const { return ((xy)?fbnY.low():fbnX.low());       }
    double high(bool xy)   const { return ((xy)?fbnY.high():fbnX.high());     }
    double center(bool xy) const { return ((xy)?fbnY.center():fbnX.center()); }

    bool   inside(double xVal, double yVal) const { return (fbnX.isInside(xVal)&&fbnY.isInside(yVal)); }

    int icvs(int iX, int iY) { return (fbnY.getN()-iY-1)*fbnX.getN()+iX+1; }

    TObjArray *rangeBoxGrid() const;
    TGraph *rangeBox(int binX=-1, int binY=-1) const;

    TString print(bool pout=1) const {
      TString line1 = fbnX.print(pout);
      TString line2 = fbnY.print(pout);
      return (line1 + ", " + line2);
    }

    void operator=(const binning2 binn2) {
      fbnX = binn2.bx();
      fbnY = binn2.by();
    }
};

ejungwoo::binning2 ejungwoo::binning::operator*(const binning binn) {
  return binning2(copy(), binn);
}
ejungwoo::binning2 ejungwoo::binning::mult(const binning *binn) {
  return binning2(copy(), binn->copy());
}

TObjArray *ejungwoo::binning2::rangeBoxGrid() const {
  auto lineArray = new TObjArray(100);
  auto bnx = bx();
  auto bny = by();
  bnx.reset();
  while (bnx.next()) {
                    lineArray -> Add(new TLine(bnx.low(), bny.min(),bnx.low(), bny.max()));
    if (bnx.dsp())  lineArray -> Add(new TLine(bnx.high(),bny.min(),bnx.high(),bny.max()));
  } if (!bnx.dsp()) lineArray -> Add(new TLine(bnx.high(),bny.min(),bnx.high(),bny.max()));
  bny.reset();
  while (bny.next()) {
                    lineArray -> Add(new TLine(bnx.min(),bny.low(), bnx.max(),bny.low()));
    if (bnx.dsp())  lineArray -> Add(new TLine(bnx.min(),bny.high(),bnx.max(),bny.high()));
  } if (!bnx.dsp()) lineArray -> Add(new TLine(bnx.min(),bny.high(),bnx.max(),bny.high()));
  return lineArray;
}

TGraph *ejungwoo::binning2::rangeBox(int binX, int binY) const {
  auto graph = new TGraph();
  auto x1 = fbnX.lowEdge(binX);
  auto x2 = fbnX.highEdge(binX);
  auto y1 = fbnY.lowEdge(binY);
  auto y2 = fbnY.highEdge(binY);
  if (binX < 0) {
    x1 = fbnX.getMin();
    x2 = fbnX.getMax();
  }
  if (binY < 0)  {
    y1 = fbnY.getMin();
    y2 = fbnY.getMax();
  }
  graph -> SetPoint(graph->GetN(),x1,y1);
  graph -> SetPoint(graph->GetN(),x1,y2);
  graph -> SetPoint(graph->GetN(),x2,y2);
  graph -> SetPoint(graph->GetN(),x2,y1);
  graph -> SetPoint(graph->GetN(),x1,y1);
  return graph;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void ejungwoo::verbose(int level) {
  fVerboseLevel = (verboseLevel)level;
  cout << "verbose is kQuiet(0) " << endl;
  cout << "verbose is kNormal(1)" << endl;
  cout << "verbose is kAll(2)   " << endl;
}

double ejungwoo::get_x(TGraph *graph, int i) { double x, y; graph->GetPoint(i,x,y); return x; }
double ejungwoo::get_y(TGraph *graph, int i) { double x, y; graph->GetPoint(i,x,y); return y; }

/// Get configuration parameter container from [nameConf].conf file.
KBParameterContainer *ejungwoo::conf(const char *nameConf)
{
  KBParameterContainer *par = nullptr;
  if (fParameterArray==nullptr) fParameterArray = new TClonesArray("KBParameterContainer",10);
  else par = (KBParameterContainer *) fParameterArray -> FindObject(nameConf);
  if (par==nullptr) {
    par = (KBParameterContainer *) fParameterArray -> ConstructedAt(fParameterArray -> GetEntriesFast());
    const char *inputFull = Form("%s/%s.conf",gSystem -> Getenv("NSEJWINPUTPATH"),nameConf);
    const char *localConf = Form("conf_ej/%s.conf",nameConf);
    if (gSystem -> Which(".",inputFull)!=nullptr)
      par -> AddFile(inputFull);
    else if (gSystem -> Which(".",localConf)!=nullptr)
      par -> AddFile(localConf);
    else
      par -> AddFile(Form("%s.conf",nameConf));
    par -> SetName(nameConf);
  }
  return par;
}

bool ejungwoo::findConf(TVirtualPad *vpad, TString &nameConfCvs, int &nx, int &ny, int &cidx)
{
  if (vpad!=nullptr) {
    TString titleCvs = vpad -> GetTitle();
    if (titleCvs.Index("____  ")<0)
      return false;
    int indexConf = 1;
    if (titleCvs.Index("____")==0)
      indexConf = 0;
    //cout << "========== " << titleCvs << endl;
    titleCvs.ReplaceAll("  ____  ",";;");
    titleCvs.ReplaceAll("____  ",";;");
    auto tokens1 = titleCvs.Tokenize(";;");
    TString titleConf = ((TObjString *) tokens1->At(indexConf))->GetString();
    //cout << "========== " << titleConf << endl;
    auto tokens2 = titleConf.Tokenize(".");
    if (tokens2->GetEntries()>3) {
      nameConfCvs = ((TObjString *) tokens2->At(0))->GetString();
      nx = ((TObjString *) tokens2->At(1))->GetString().Atoi();
      ny = ((TObjString *) tokens2->At(2))->GetString().Atoi();
      cidx = ((TObjString *) tokens2->At(3))->GetString().Atoi();
      return true;
    }
  }
  return false;
}

ejungwoo::TParCanvas ejungwoo::getParCanvas(parContainer *par)
{
  if (par -> CheckPar("set_stats")) fParCvs.setStats = par -> GetParBool("set_stats");
  if (par -> CheckPar("set_div_idx_yx")) fParCvs.setDivIndexYX = par -> GetParBool("set_div_idx_yx");
  fParCvs.textFont = par -> GetParInt("font");
  if (par->CheckPar("text_size")) fParCvs.textSize = par -> GetParInt("text_size");
  else fParCvs.textSize = 0;
  fParCvs.menuSize = par -> GetParInt("menu_size");
  fParCvs.removeInnerZaxis = par -> GetParBool("remove_inn_zaxis");
  fParCvs.removeInnerMainTitle = par -> GetParBool("remove_inn_main_title");
  if (par -> GetParN("remove_inn_pad_axis")>0) {
    fParCvs.removeInnerPadAxis[0] = par -> GetParBool("remove_inn_pad_axis",0);
    fParCvs.removeInnerPadAxis[1] = par -> GetParBool("remove_inn_pad_axis",1);
  } else {
    fParCvs.removeInnerPadAxis[0] = par -> GetParBool("remove_inn_pad_axis");
    fParCvs.removeInnerPadAxis[1] = par -> GetParBool("remove_inn_pad_axis");
  }
  for (int ixy : {0,1}) {
    fParCvs.padSize[ixy] = par -> GetParInt("pad_size",ixy);
    fParCvs.marginDiv[ixy] = par -> GetParInt("margin_div",ixy);
  }
  for (int ixyz : {0,1,2}) {
    if (par->CheckPar("axis_ndivision"))
      fParCvs.axisNdivision[ixyz] = par -> GetParInt("axis_ndivision",ixyz);
    else
      fParCvs.axisNdivision[ixyz] = 510;
    fParCvs.labelSize[ixyz] = par -> GetParDouble("label_size",ixyz);
    fParCvs.labelOffset[ixyz] = par -> GetParDouble("label_offset",ixyz);
    fParCvs.tickSize[ixyz] = par -> GetParDouble("tick_size",ixyz);
  }
  for (int ixyzm : {0,1,2,3}) {
    fParCvs.titleAlign[ixyzm] = par -> GetParInt("title_align",ixyzm);
    if (fParCvs.titleAlign[ixyzm]<0) {
      fParCvs.titleAlign[ixyzm] = -fParCvs.titleAlign[ixyzm];
      fParCvs.titleRotate[ixyzm] = true;
    } else
      fParCvs.titleRotate[ixyzm] = false;
    fParCvs.titleSize[ixyzm] = par -> GetParDouble("title_size",ixyzm);
    fParCvs.titleOffset[ixyzm] = par -> GetParDouble("title_offset",ixyzm);
    fParCvs.titleColor[ixyzm] = par -> GetParInt("title_color",ixyzm);
  }
  for (int img : {0,1,2,3}) {
    fParCvs.marginTop[img] = par -> GetParInt("margin_top",img);
    fParCvs.marginInn[img] = par -> GetParInt("margin_inn",img);
  }
  if (par -> CheckPar("num_contours"))
    fParCvs.numContours = par -> GetParInt("num_contours");
  else
    fParCvs.numContours = 100;
  if (par -> CheckPar("side_pad")) {
    fParCvs.sidePad[0] = par -> GetParInt("side_pad",0);
    fParCvs.sidePad[1] = par -> GetParInt("side_pad",1);
  } else {
    fParCvs.sidePad[0] = 0;
    fParCvs.sidePad[1] = 0;
  }

  return fParCvs;
}

ejungwoo::TParAttribute ejungwoo::getParAttribute(parContainer *par)
{
  fParAtt.lineStyle = par -> GetParVInt("line_style");
  fParAtt.lineWidth = par -> GetParVInt("line_width");
  fParAtt.markerSize = par -> GetParVDouble("marker_size");
  fParAtt.markerStyle = par -> GetParVInt("marker_style");
  fParAtt.fillStyle = par -> GetParVInt("fill_style");
  if (par -> CheckPar("text_font")) fParAtt.textFont = par -> GetParInt("text_font");
  if (par -> CheckPar("text_size")) fParAtt.textSize = par -> GetParInt("text_size");
  if (par -> CheckPar("text_align")) fParAtt.textAlign = par -> GetParInt("text_align");
  if (par -> CheckPar("all_color")) {
    fParAtt.lineColor   = par -> GetParVInt("all_color");
    fParAtt.markerColor = par -> GetParVInt("all_color");
    fParAtt.fillColor   = par -> GetParVInt("all_color");
    fParAtt.textColor   = par -> GetParVInt("all_color");
  } else {
    fParAtt.lineColor = par -> GetParVInt("line_color");
    fParAtt.markerColor = par -> GetParVInt("marker_color");
    fParAtt.fillColor = par -> GetParVInt("fill_color");
    if (par -> CheckPar("text_color")) fParAtt.textColor = par -> GetParVInt("text_color");
  }
  return fParAtt;
}

ejungwoo::TParCanvas    ejungwoo::getParCanvas   (TString name) { return getParCanvas   (conf(name)); }
ejungwoo::TParAttribute ejungwoo::getParAttribute(TString name) { return getParAttribute(conf(name)); }

/// Make new pave text (TPaveText)
/// @param content  text
/// @param vpad  canvas
/// @param x1  x1 (left position) in ratio (0.~1.)
/// @param y1  y1 (bottom position) in ratio (0.~1.)
/// @param dx  dx in ratio (0.~1.)
/// @param dy  dy in ratio (0.~1.)
/// @param tf  text font
/// @param tz  text size
/// @param ta  text align
/// @param fc  font color
/// @param fs  font size
/// @param bs  border size
TPaveText* ejungwoo::newpt(TString content, TVirtualPad *vpad, double x1, double y1, double dx, double dy, int tf, double tz, int ta, int fs, int fc, int bs)
{
  TPaveText *paveText = nullptr;
  if (content.Index(";cord")>=0) {
    content.ReplaceAll(";cord","");
    double x2 = x1 + dx;
    double y2 = y1 + dy;
    if (content.Index(";xyxy")) {
      content.ReplaceAll(";xyxy","");
      x2 = dx;
      y2 = dy;
    }
    paveText = new TPaveText(x1,y1,x2,y2,"NB");
  } else {
    double xu, yu;
    findxy(vpad, x1, y1, dx, dy, xu, yu);
    paveText = new TPaveText(x1,y1,x1+dx,y1+dy,"NDCNB");
  }
  paveText -> AddText(content);
  paveText -> SetTextFont(tf);
  paveText -> SetTextSize(tz);
  paveText -> SetTextAlign(ta);
  paveText -> SetFillColor(fc);
  paveText -> SetFillStyle(fs);
  paveText -> SetBorderSize(bs);
  return paveText;
}

TPad *ejungwoo::newpad(TVirtualPad *vpad, TString name, TString title, double x1, double y1, double dx, double dy)
{
  double xu, yu;
  findxy(vpad, x1, y1, dx, dy, xu, yu);
  cout << x1 << " " << y1 << " " << x1+dx << " " << y1+dy << endl;
  auto pad2 = new TPad(name,title,x1,y1,x1+dx,y1+dy);
  return pad2;
}

TLatex* ejungwoo::newtt(TString content, double x, double y, int tf, double tz, int ta, int tc) {
  auto tt = new TLatex(x,y,content);
  tt -> SetTextFont(tf);
  tt -> SetTextSize(tz);
  tt -> SetTextAlign(ta);
  tt -> SetTextColor(tc);
  return tt;
}

TCanvas *ejungwoo::canvas(const char *nameCvs, int nn, const char *nameConf)
{
  int nx=10, ny=10;
  if (nn==1) return canvas(nameCvs,nameConf);
  else if (nn<=2) { nx = 2; ny = 1; }
  else if (nn<=4) { nx = 2; ny = 2; }
  else if (nn<=6) { nx = 3; ny = 2; }
  else if (nn<=8) { nx = 4; ny = 2; }
  else if (nn<=12) { nx = 4; ny = 3; }
  else if (nn<=15) { nx = 5; ny = 3; }
  else if (nn<=20) { nx = 5; ny = 4; }
  else if (nn<=24) { nx = 6; ny = 4; }
  else if (nn<=30) { nx = 6; ny = 5; }
  else if (nn<=35) { nx = 7; ny = 5; }
  else if (nn<=40) { nx = 8; ny = 5; }
  else if (nn<=50) { nx = 10; ny = 5; }
  else if (nn<=60) { nx = 10; ny = 6; }
  else if (nn<=70) { nx = 10; ny = 7; }
  else if (nn<=80) { nx = 10; ny = 8; }
  else if (nn<=90) { nx = 10; ny = 9; }
  return canvas(nameCvs,nx,ny,nameConf);
}

/// Create canvas with nameCvs and division number nx and ny.
/// Set nx and ny is 1 for single pad canvas.
/// Configuration of the canvas(+ histogram and legened) is set by [nameConf].conf file.
TCanvas *ejungwoo::canvas(const char *nameCvs, int nx, int ny, const char *nameConf)
{
  if (strcmp(nameCvs,"")==0)
    nameCvs = Form("canvas_%d",fCanvasIndex);
  fCanvasIndex++;

  if (TString(nameConf).IsNull()) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  getParCanvas(par);

  bool sidePadExist = true;
  if (fParCvs.sidePad[0]==0)
    sidePadExist = false;

  // size
  int innSize[2] =  {
    fParCvs.padSize[0]+fParCvs.marginInn[0]+fParCvs.marginInn[1],
    fParCvs.padSize[1]+fParCvs.marginInn[2]+fParCvs.marginInn[3]};
  int topSize[2] =  {
    fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[0]*nx+(nx-1+sidePadExist*1)*fParCvs.marginDiv[0]+sidePadExist*(fParCvs.sidePad[0]),
    fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1]*ny+(ny-1)*fParCvs.marginDiv[1]};
  int topSize0[2] = {
    fParCvs.marginTop[0]+fParCvs.marginTop[1]+fParCvs.padSize[0],
    fParCvs.marginTop[2]+fParCvs.marginTop[3]+fParCvs.padSize[1]};

  int fullSize[2] = {topSize[0], topSize[1]+fParCvs.menuSize};
  //if (nx==1&&ny==1) { fullSize[0] = topSize0[0]; fullSize[1] = topSize0[1]+fParCvs.menuSize; }

  // margin
  double marginRatioDiv[2] = {double(fParCvs.marginDiv[0])/topSize[0], double(fParCvs.marginDiv[1])/topSize[1]};
  double marginRatioTop[4]; for (auto img : {0,1,2,3}) marginRatioTop[img] = double(fParCvs.marginTop[img])/topSize[int(img<2?0:1)];
  //if (nx==1&&ny==1)         for (auto img : {0,1,2,3}) marginRatioTop[img] = double(fParCvs.marginTop[img])/topSize0[int(img<2?0:1)];
  double marginRatioPad[4]; for (auto img : {0,1,2,3}) marginRatioPad[img] = double(fParCvs.marginInn[img])/innSize[int(img<2?0:1)]; 
  double marginRatioEdg[4]; for (auto img : {0,1,2,3}) marginRatioEdg[img] = double(fParCvs.marginTop[img]+fParCvs.marginInn[img])/(fParCvs.marginTop[img]+innSize[int(img<2?0:1)]);
  double marginRatioEd2[4] = {
    double(fParCvs.marginTop[0]+fParCvs.marginInn[0])/(fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[0]),
    double(fParCvs.marginTop[1]+fParCvs.marginInn[1])/(fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[0]),
    double(fParCvs.marginTop[2]+fParCvs.marginInn[2])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1]),
    double(fParCvs.marginTop[3]+fParCvs.marginInn[3])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1]),
  };

  // dxInRatio dyInRatio pad
  double padRatioInn[2] = {double(innSize[0])/topSize[0], double(innSize[1])/topSize[1]};
  double padRatioInnEdge1[4]; for (auto img : {0,1,2,3}) padRatioInnEdge1[img] = double(fParCvs.marginTop[img]+innSize[int(img<2?0:1)])/topSize[int(img<2?0:1)];
  double padRatioInnEdge2[2] = {double(fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[0])/topSize[0], double(fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1])/topSize[1]};
  double padRatioSidePad[4] = {double(fParCvs.sidePad[0])/topSize[0], double(fParCvs.sidePad[1])/topSize[1]};

  //const char *titleCvs = Form("%s.%d.%d.%d.   %s",nameConf,nx,ny,0,nameCvs);
  const char *titleCvs = Form("%s  ____  %s.%d.%d.%d.",nameCvs,nameConf,nx,ny,0);
  auto cvs = new TCanvas(nameCvs, titleCvs, fCanvasIndex*10, 25+fCanvasIndex*10, fullSize[0], fullSize[1]);
  cvs -> SetMargin(marginRatioTop[0],marginRatioTop[1],marginRatioTop[2],marginRatioTop[3]);

  if (fCanvasArray==nullptr) fCanvasArray = new TObjArray(50);
  fCanvasArray -> Add(cvs);

  if (nx>1||ny>1||TString(nameConf).Index("nn")>=0||sidePadExist)
  //if (TString(nameConf).Index("nn")>=0)
  {
    int countPad = 0;

    vector<int> idxX;
    vector<int> idxY;
    if (fParCvs.setDivIndexYX) {
      for (auto ix=0; ix<nx; ++ix) {
        for (auto iy=0; iy<ny; ++iy) {
          idxX.push_back(ix);
          idxY.push_back(iy);
        }
      }
    } else {
      for (auto iy=0; iy<ny; ++iy) {
        for (auto ix=0; ix<nx; ++ix) {
          idxX.push_back(ix);
          idxY.push_back(iy);
        }
      }
    }

    int numIdexes = idxX.size();
    for (int idx=0; idx<numIdexes; ++idx)
    {
      auto ix = idxX[idx];
      auto iy = idxY[idx];

      bool isTop = (iy==0);
      bool isBottom = (iy==ny-1);
      bool isLeft = (ix==0);
      bool isRight = ((fParCvs.sidePad[0]==0)&&(ix==nx-1));

      cvs -> cd();
      countPad++;

      double y2 = 1 - (iy>0)*(padRatioInnEdge1[3] + (iy-1)*padRatioInn[1]) - iy*marginRatioDiv[1];
      double y1 = y2 - padRatioInn[1];
      if (isTop&&isBottom) y1 = y2 - padRatioInnEdge2[1];
      else if (isTop)      y1 = y2 - padRatioInnEdge1[3];
      else if (isBottom)   y1 = y2 - padRatioInnEdge1[2];
      if (y1<0) y1 = 0;

      double x1 = (ix>0)*(padRatioInnEdge1[0] + (ix-1)*padRatioInn[0]) + ix*marginRatioDiv[0];
      double x2 = x1 + padRatioInn[0];
      if (isLeft&&isRight) x2 = x1 + padRatioInnEdge2[0];
      else if (isLeft)     x2 = x1 + padRatioInnEdge1[0];
      else if (isRight)    x2 = x1 + padRatioInnEdge1[1];
      if (x1<0) x1 = 0;

      //cout << par->GetName() << endl;
      //if (TString(par->GetName())=="paper_nn6") {
      //  if (idx<2) y1 = 0.40;
      //  else y2 = 0.40;
      //}

      const char *namePad = Form("%s_%d",cvs->GetName(),countPad);
      const char *titlePad = Form("____  %s.%d.%d.%d",nameConf,nx,ny,countPad);
      auto padi = new TPad(namePad,titlePad,x1,y1,x2,y2);
      padi -> SetNumber(countPad);
      padi -> SetMargin(marginRatioPad[0],marginRatioPad[1],marginRatioPad[2],marginRatioPad[3]);
      padi -> Draw();

      if (isLeft&&isRight) {
        padi -> SetLeftMargin (marginRatioEd2[0]);
        padi -> SetRightMargin(marginRatioEd2[1]);
      }
      else if (isLeft)    padi -> SetLeftMargin (marginRatioEdg[0]);
      else if (isRight) padi -> SetRightMargin(marginRatioEdg[1]);

      if (isTop&&isBottom) {
        padi -> SetBottomMargin(marginRatioEd2[2]);
        padi -> SetTopMargin   (marginRatioEd2[3]);
      }
      else if (isBottom) padi -> SetBottomMargin(marginRatioEdg[2]);
      else if (isTop)    padi -> SetTopMargin   (marginRatioEdg[3]);
    }

    if (sidePadExist) {
      cvs -> cd();
      countPad=99;

      const char *namePad = Form("%s_%d",cvs->GetName(),countPad);
      const char *titlePad = Form("____  %s.%d.%d.%d",nameConf,countPad,countPad,countPad);
      double x1 = (padRatioInnEdge1[0] + (nx-1)*padRatioInn[0]) + nx*marginRatioDiv[0];
      double x2 = x1 + padRatioSidePad[0];
      double y2 = 1 - (0>0)*(padRatioInnEdge1[3] + (0-1)*padRatioInn[1]) - 0*marginRatioDiv[1];
      double y1 = y2 - padRatioInnEdge2[1] - padRatioInn[1];;

      auto padi = new TPad(namePad,titlePad,x1,y1,x2,y2);
      cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
      padi -> SetNumber(countPad);
      padi -> SetLeftMargin (0);
      padi -> SetRightMargin (double(fParCvs.marginTop[1]+fParCvs.marginInn[1])/(fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[1]));
      padi -> SetBottomMargin(double(fParCvs.marginTop[2]+fParCvs.marginInn[2])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+2*innSize[1]));
      padi -> SetTopMargin   (double(fParCvs.marginTop[3]+fParCvs.marginInn[3])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+2*innSize[1]));

      padi -> Draw();
    }
  }

  return cvs;
}

TPad *ejungwoo::innerpad(TVirtualPad *vpad, double xr1, double yr1, double xr2, double yr2)
{
  auto pad = (TPad *) vpad;
  int padH = pad -> GetWh();
  int padW = pad -> GetWw();
  double ml = pad -> GetLeftMargin();
  double mr = pad -> GetRightMargin();
  double mb = pad -> GetBottomMargin();
  double mt = pad -> GetTopMargin();
  double xMin = ((0. + ml) * padW) + 1;
  double xMax = ((1. - mr) * padW) - 1;
  double yMin = ((0. + mb) * padH) + 1;
  double yMax = ((1. - mt) * padH) - 1;
  double x1 = ((1.-xr1)*xMin + xr1*(xMax)) / padW;
  double x2 = ((1.-xr2)*xMin + xr2*(xMax)) / padW;
  double y1 = ((1.-yr1)*yMin + yr1*(yMax)) / padH;
  double y2 = ((1.-yr2)*yMin + yr2*(yMax)) / padH;
  auto pad2 = new TPad(Form("%s_innerpad",pad->GetName()),"",x1,y1,x2,y2);
  return pad2;
}

TH1 *ejungwoo::project(TTree *tree, TString formula, TCut cut, TString name, TString title, int nx, double x1, double x2, int ny, double y1, double y2)
{
  TString bnamex;
  TString bnamey;
  if(formula.Index(":")>=0) {
    TObjArray *elements = formula.Tokenize(":");
    TString bnamex=((TObjString *)elements->At(1))->GetString();
    TString bnamey=((TObjString *)elements->At(0))->GetString();
  }
  else
    bnamex = formula;

  if(nx<1) nx=200;
  if(formula.Index(":")>=0 && ny<1) ny=200;

  if (x1==0&&x2==0) {
    x1 = tree -> GetMinimum(bnamex);
    x2 = tree -> GetMaximum(bnamex);
    auto x3 = (x2-x1)/10.;
    x2 = x2+x3;
    x1 = x1-x3;
  }
  if (ny>0&&y1==0&&y2==0) {
    y1 = tree->GetMinimum(bnamey);
    y2 = tree->GetMaximum(bnamey);
    auto y3 = (y2-y1)/10.;
    y2 = y2+y3;
    y1 = y1-y3;
  }

  TString title0;
  if (name.IsNull())
    name = Form("hist_%d",fHistIndex++);
  //else title0 = name;

  TH1 *histProjected;
  if(ny>0) {
    if(title.IsNull()) title = title0+";"+bnamex+";"+bnamey;
    else if (title.Index(";")<0)
      title = title+";"+bnamex+";"+bnamey;
    histProjected = new TH2D(name,title,nx,x1,x2,ny,y1,y2);
  }
  else {
    if(title.IsNull()) title = title0+";"+formula;
    histProjected = new TH1D(name,title,nx,x1,x2);
    if (y1!=0||y2!=0) {
      histProjected -> SetMinimum(y1);
      histProjected -> SetMinimum(y2);
    }
  }

  Long64_t entries = tree->Project(name,formula,cut);

  return histProjected;
}

TPad *ejungwoo::redraw(TPad *pad)
{
  auto list = pad -> GetListOfPrimitives();
  TIter next(list);

  int topIsFound = 0;
  TH1 *histTop1 = nullptr;
  TH2 *histTop2 = nullptr;
  TObjArray objectList;
  double xMin=DBL_MAX, xMax=-DBL_MAX, yMin=DBL_MAX, yMax=-DBL_MAX;

  while (auto prim = next())
  {
    TString nameObj = prim -> GetName();
    TString className = prim -> ClassName();

    if (prim->InheritsFrom(TFrame::Class()))
      continue;

    else if (prim->InheritsFrom(TH2::Class())) {
      auto hist = (TH2 *) prim;
      if (topIsFound==0) {
        histTop2 = hist;
        topIsFound = 2;
      } else {
        double x1,x2,y1,y2;
        x1 = hist -> GetXaxis() -> GetBinLowEdge(1);
        x2 = hist -> GetXaxis() -> GetBinUpEdge(hist->GetXaxis()->GetNbins());
        y1 = hist -> GetYaxis() -> GetBinLowEdge(1);
        y2 = hist -> GetYaxis() -> GetBinUpEdge(hist->GetYaxis()->GetNbins());
        if (x1<xMin) xMin = x1;
        if (x2>xMax) xMax = x2;
        if (y1<yMin) yMin = y1;
        if (y2>yMax) yMax = y2;
      }
    }

    else if (prim->InheritsFrom(TH1::Class())) {
      auto hist = (TH1 *) prim;
      if (topIsFound==0) {
        histTop1 = hist;
        topIsFound = 1;
      } else {
        double x1,x2,y1,y2;
        x1 = hist -> GetXaxis() -> GetBinLowEdge(1);
        x2 = hist -> GetXaxis() -> GetBinUpEdge(hist->GetXaxis()->GetNbins());
        y1 = hist -> GetMinimum();
        y2 = hist -> GetMaximum();
        if (x1<xMin) xMin = x1;
        if (x2>xMax) xMax = x2;
        if (y1<yMin) yMin = y1;
        if (y2>yMax) yMax = y2;
      }
    }

    else if (prim->InheritsFrom(TGraph::Class())) {
      auto graph = (TGraph *) prim;
      double x1,x2,y1,y2;
      graph -> ComputeRange(x1,y1,x2,y2);
      if (x1<xMin) xMin = x1;
      if (x2>xMax) xMax = x2;
      if (y1<yMin) yMin = y1;
      if (y2>yMax) yMax = y2;
      if (topIsFound==0) {
        topIsFound = 3;
      }
    }

    else if (prim->InheritsFrom(TF1::Class())) {
      auto f1 = (TF1 *) prim;
      double x1,x2,y1,y2;
      f1 -> GetRange(x1, x2);
      y1 = f1 -> GetMinimum();
      y2 = f1 -> GetMaximum();
      if (x1<xMin) xMin = x1;
      if (x2>xMax) xMax = x2;
      if (y1<yMin) yMin = y1;
      if (y2>yMax) yMax = y2;
      if (topIsFound==0) {
        topIsFound = 4;
      }
    }
    else
      continue;

    objectList.Add(prim);
  }

  if (topIsFound==1) {
    histTop1 -> SetMaximum(yMax*1.1);
    histTop1 -> Print();
    histTop1 -> Draw();
    pad -> Modified();
    pad -> Update();
  }

  if (topIsFound==2) {
    auto yMin0 = histTop2 -> GetYaxis() -> GetBinLowEdge(1);
    histTop2 -> GetYaxis() -> SetRangeUser(yMin0,yMax*1.1);
    pad -> Modified();
    pad -> Update();
  }

  //objectList.ls();

  return pad;
}

TGraphErrors *ejungwoo::tograph(TH1D *hist, TString opt, int iatt, TString nameConf)
{
  binning bnx(hist);
  auto graph = new TGraphErrors();
  bool ignore0 = false;
  bool reversex = false;
  bool reversey = false;
  bool copyex = true;
  bool copyey = true;
  double xlow = -DBL_MAX;
  double xhigh = DBL_MAX;
  double ylow = -DBL_MAX;
  double yhigh = DBL_MAX;
  vector<double> addatxarray;

  auto optArray = opt.Tokenize(";");
  auto numArray = optArray->GetEntries();
  for (auto i=0; i<numArray; ++i) {
    TString option = tok(optArray,i);
    if (option=="0") ignore0 = true;
    if (option=="-x") reversex = true;
    if (option=="-y") reversey = true;
    if (option=="ex") copyex = true;
    if (option=="ey") copyey = true;
    if (option=="!ex") copyex = false;
    if (option=="!ey") copyey = false;
    if (option.Index("x<")==0) { option.ReplaceAll("x<",""); xhigh=option.Atof(); }
    if (option.Index("x>")==0) { option.ReplaceAll("x>",""); xlow =option.Atof(); }
    if (option.Index("y<")==0) { option.ReplaceAll("y<",""); yhigh=option.Atof(); }
    if (option.Index("y>")==0) { option.ReplaceAll("y>",""); ylow =option.Atof(); }
    if (option.Index("+x=")==0) { option.ReplaceAll("+x=",""); addatxarray.push_back(option.Atof()); }
  }

  while (bnx.next()) {
    if (bnx.val()<xlow ) continue;
    if (bnx.val()>xhigh) continue;
    if (bnx.val()<xlow ) continue;
    if (bnx.val()>xhigh) continue;

    auto content = hist -> GetBinContent(bnx.bi());
    if (ignore0 && content==0) {}
    else {
      double xvalue = bnx.val();
      double yvalue = content;
      if (reversex) xvalue = -xvalue;
      if (reversey) yvalue = -yvalue;
      graph -> SetPoint(graph->GetN(), xvalue, yvalue);
      double xerror = bnx.width()/2.;
      double yerror = hist -> GetBinError(bnx.bi());
      if (!copyex) xerror = 0;
      if (!copyey) yerror = 0;
      graph -> SetPointError(graph->GetN()-1, xerror, yerror);
    }
  }

  auto npoints = graph -> GetN();
  if (npoints>=2) {
    for (auto addatx : addatxarray) {
      int ibin1,ibin2;
      double dxmin = DBL_MAX;
      for (auto ii=0; ii<npoints; ++ii) {
        auto x0 = get_x(graph,ii);
        if (abs(x0-addatx)<dxmin) {
          ibin2 = ibin1;
          ibin1 = ii;
          dxmin = abs(x0-addatx);
        }
      }
      if (ibin1==npoints-1) {
        ibin2 = ibin1-1;
        ibin1 = ibin2+1;
      } else if (ibin1==0) {
        ibin1 = 0;
        ibin2 = 1;
      }

      auto yatx = graph -> Eval(addatx);
      double xvalue1 = get_x(graph,ibin1);
      double xerror1 = graph -> GetErrorX(ibin1);
      double yerror1 = graph -> GetErrorY(ibin1);
      double xvalue2 = get_x(graph,ibin2);
      double xerror2 = graph -> GetErrorX(ibin2);
      double yerror2 = graph -> GetErrorY(ibin2);
      double xerror = ((addatx-xvalue1)*(xerror2-xerror1)/(xvalue2-xvalue1))+xerror1;
      double yerror = ((addatx-xvalue1)*(yerror2-yerror1)/(xvalue2-xvalue1))+yerror1;
      if (!copyex) xerror = 0;
      if (!copyey) yerror = 0;
      graph -> SetPoint(graph->GetN(),addatx,yatx);
      graph -> SetPointError(graph->GetN()-1, xerror, yerror);
    }
  }

  if (iatt>=0) {
    if (nameConf.IsNull())
      nameConf = opt;
    ejungwoo::att(graph,iatt,nameConf);
  }

  return graph;
}

TString ejungwoo::tok(TString line, TString token, int i)
{
  auto array = line.Tokenize(token);
  if (array -> GetEntries()==0)
    return "";
  return tok(line.Tokenize(token),i);
}

TString ejungwoo::tok(TObjArray *line, int i) {
  return ((TObjString *) line->At(i))->GetString();
}

void ejungwoo::savePDF(TCanvas *cvs, TString nameVersion) {
  cvs -> Modified();
  cvs -> Update();
  if (nameVersion.IsNull()) nameVersion = fNameVersion;
  if (nameVersion.IsNull()) nameVersion = ".";
  else gSystem -> mkdir(nameVersion);
  TString pathToData = nameVersion + "/pdf/";
  gSystem -> mkdir(pathToData);
  cvs -> SaveAs(pathToData+cvs->GetName()+".pdf");
}

void ejungwoo::saveEPS(TCanvas *cvs, TString nameVersion) {
  cvs -> Modified();
  cvs -> Update();
  if (nameVersion.IsNull()) nameVersion = fNameVersion;
  if (nameVersion.IsNull()) nameVersion = ".";
  else gSystem -> mkdir(nameVersion);
  TString pathToData = nameVersion + "/eps/";
  gSystem -> mkdir(pathToData);
  cvs -> SaveAs(pathToData+cvs->GetName()+".eps");
}

void ejungwoo::savePNG(TCanvas *cvs, TString nameVersion) {
  cvs -> Modified();
  cvs -> Update();
  if (nameVersion.IsNull()) nameVersion = fNameVersion;
  if (nameVersion.IsNull()) nameVersion = ".";
  else gSystem -> mkdir(nameVersion);
  TString pathToData = nameVersion + "/figures/";
  gSystem -> mkdir(pathToData);
  cvs -> SaveAs(pathToData+cvs->GetName()+".png");
}

void ejungwoo::savePDF(TString nameVersion) {
  if (nameVersion.IsNull()) nameVersion = fNameVersion;
  if (nameVersion.IsNull()) nameVersion = ".";
  else gSystem -> mkdir(nameVersion);
  TString pathToData = nameVersion + "/pdf/";
  gSystem -> mkdir(pathToData);
  auto numCanvases = fCanvasArray -> GetEntriesFast();
  for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
    auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
    cvs -> Modified();
    cvs -> Update();
    cvs -> SaveAs(pathToData+cvs->GetName()+".pdf");
  }
}

void ejungwoo::savePNG(TString nameVersion) {
  if (nameVersion.IsNull()) nameVersion = fNameVersion;
  if (nameVersion.IsNull()) nameVersion = ".";
  else gSystem -> mkdir(nameVersion);
  TString pathToData = nameVersion + "/figures/";
  gSystem -> mkdir(pathToData);
  auto numCanvases = fCanvasArray -> GetEntriesFast();
  for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
    auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
    cvs -> Modified();
    cvs -> Update();
    cvs -> SaveAs(pathToData+cvs->GetName()+".png");
  }
}

void ejungwoo::saveEPS(TString nameVersion) {
  if (nameVersion.IsNull()) nameVersion = fNameVersion;
  if (nameVersion.IsNull()) nameVersion = ".";
  else gSystem -> mkdir(nameVersion);
  TString pathToData = nameVersion + "/eps/";
  gSystem -> mkdir(pathToData);
  auto numCanvases = fCanvasArray -> GetEntriesFast();
  for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
    auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
    cvs -> Modified();
    cvs -> Update();
    cvs -> SaveAs(pathToData+cvs->GetName()+".eps");
  }
}

void ejungwoo::saveRoot(TObject *obj, TString nameFile, TString nameVersion, bool savePrimitives, bool simplifyNames)
{
  if (nameVersion.IsNull()) nameVersion = fNameVersion;
  if (nameVersion.IsNull()) nameVersion = ".";
  else gSystem -> mkdir(nameVersion);
  TString pathToData = nameVersion + "/rooto/";
  gSystem -> mkdir(pathToData);

  TFile *fileOut = nullptr;
  if (!nameFile.IsNull()) {
    if (fVerboseLevel>=verboseLevel::kNormal)
      cout_info << "Creating " << nameFile << endl;
    fileOut = new TFile(nameFile, "recreate");
  }

  auto saveCanvasRoot = [fileOut,savePrimitives,pathToData,simplifyNames](TPad *cvs)
  {
    cvs -> Modified();
    cvs -> Update();
    TFile *fileCvsOut = fileOut;
    if (fileOut==nullptr) {
      TString nameCvsFile = pathToData+cvs->GetName()+".root";
      if (fVerboseLevel>=verboseLevel::kNormal)
        cout_info << "Creating " << nameCvsFile << endl;
      fileCvsOut = new TFile(nameCvsFile,"recreate");
    }
    if (simplifyNames) cvs -> Write("cvs");
    else cvs -> Write();
    if (savePrimitives) {
      TObjArray nameArray;
      TObjArray nameArray2;
      TObjArray laterArray;
      auto list = cvs -> GetListOfPrimitives();
      TIter next(list);
      while (auto prim = next()) {
        TString nameObj = prim -> GetName();
        TString className = prim -> ClassName();
        if (nameObj.IsNull()) nameObj = className;

        if (simplifyNames) {
               if (prim->InheritsFrom(TH1::Class()))     nameObj = "hist";
          else if (prim->InheritsFrom(TPad::Class()))    nameObj = "pad";
          else if (prim->InheritsFrom(TF1::Class()))     nameObj = "func";
          else if (prim->InheritsFrom(TGraph::Class()))  nameObj = "graph";
        }
             if (nameObj=="TLine") nameObj = "line";
        else if (nameObj=="TText") nameObj = "text";
        else if (nameObj=="TLatex") nameObj = "text";
        else if (nameObj=="TFrame") nameObj = "frame";
        else if (nameObj=="TGraph") nameObj = "graph";
        else if (nameObj=="TMarker") nameObj = "marker";
        else if (nameObj=="TPaveText") nameObj = "pavet";
        else if (nameObj=="TGraphErrors") nameObj = "graphe";
        else if (nameObj=="TPave") { 
          if (prim->InheritsFrom(TLegend::Class())) nameObj = "legend";
          else nameObj = "pave";
        }
        TString nameObjA = nameObj;
        int idxSameName = 0;
        while (1) {
          if (idxSameName!=0)
            nameObj = nameObjA + "_" + idxSameName;
          if (nameArray.FindObject(nameObj)!=nullptr) {
            idxSameName++;
          } else
            break;
        }
        auto named = new TNamed(nameObj.Data(),"");
        nameArray.Add(named);
        if (className.Sizeof()>9) {
          nameArray2.Add(named);
          laterArray.Add(prim);
        }
        else
          prim -> Write(nameObj);
      }
      auto numLater = laterArray.GetEntriesFast();
      for (auto iLater=0; iLater<numLater; ++iLater) {
        auto prim = laterArray.At(iLater);
        prim -> Write(nameArray2.At(iLater) -> GetName());
      }
    }
    if (fileOut==nullptr)
      fileCvsOut -> Close();
  };

  if (obj==nullptr) {
    int numCanvases = 0;
    if (fCanvasArray!=nullptr)
      numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      saveCanvasRoot(cvs);
    }
  }
  else if (obj->InheritsFrom(TPad::Class()))
    saveCanvasRoot((TPad *) obj);
  else {
    if (nameFile.IsNull())
      nameFile = pathToData+obj->GetName()+".root";
    if (fVerboseLevel>=verboseLevel::kNormal)
      cout_info << "Creating " << nameFile << endl;
    TFile *fileObjOut = fileOut;
    if (fileOut==nullptr) fileObjOut = new TFile(nameFile,"recreate");
    fileObjOut -> cd();

    TString nameWrite = obj -> GetName();
    if (simplifyNames) {
           if (obj->InheritsFrom(TH1::Class()))     nameWrite = "hist";
      else if (obj->InheritsFrom(TPad::Class()))    nameWrite = "pad";
      else if (obj->InheritsFrom(TF1::Class()))     nameWrite = "func";
      else if (obj->InheritsFrom(TGraph::Class()))  nameWrite = "graph";
    }
    obj -> Write(nameWrite);
    if (fileOut==nullptr) fileObjOut -> Close();
  }
  if (fileOut!=nullptr) fileOut -> Close();
}

void ejungwoo::saveRoot(TString nameVersion, bool savePrimitives, bool simplifyNames) { saveRoot((TObject*)nullptr, "", nameVersion, savePrimitives, simplifyNames); }

/// save all canvases created so far in [nameVersion]/figures/[cvs-name].png and [nameVersion]/pdf/[cvs-name].pdf files
void ejungwoo::saveAll(TString nameVersion) {
  savePNG(nameVersion);
  saveEPS(nameVersion);
  savePDF(nameVersion);
  saveRoot(nameVersion,1,1);
}

void ejungwoo::write(TObject *obj)
{
  gSystem -> mkdir(Form("rooto"));
  TString fileName = Form("rooto/%s.root",obj -> GetName());
  auto ofile = new TFile(fileName,"recreate");
  if (fVerboseLevel>=verboseLevel::kNormal)
    cout_info << "Writting file " << fileName << "!" << endl;
  obj -> Write();
}

/// Make hist fit to the pad. If idx>0, pad where legend is drawn is selected by vpad->cd(idx). 
/// Draw histogram before make(TH1* hist, TVirtualPad *vpad) to apply main title attribute.
TH1 *ejungwoo::make(TH1 *hist, TVirtualPad *vpad, TString drawOption)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd();

  int nx = 1;
  int ny = 1;
  int cidx = 0;
  TString nameConf = "";
  if (!findConf(padi,nameConf,nx,ny,cidx)) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  getParCanvas(par);

  // size
  int innSize[2] = {fParCvs.padSize[0]+fParCvs.marginInn[0]+fParCvs.marginInn[1], fParCvs.padSize[1]+fParCvs.marginInn[2]+fParCvs.marginInn[3]};
  int topSize[2] = {fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[0]*nx+2*nx*fParCvs.marginDiv[0], fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1]*ny+2*ny*fParCvs.marginDiv[1]};
  int topSize0[2] = {fParCvs.marginTop[0]+fParCvs.marginTop[1]+fParCvs.padSize[0], fParCvs.marginTop[2]+fParCvs.marginTop[3]+fParCvs.padSize[1]};
  int fullSize[2] = {topSize[0], topSize[1]+fParCvs.menuSize};
  if (nx==1&&ny==1) {
    fullSize[0] = topSize0[0];
    fullSize[1] = topSize0[1]+fParCvs.menuSize;
  }

  // margin
  double marginRatioDiv[2] = {double(fParCvs.marginDiv[0])/topSize[0], double(fParCvs.marginDiv[1])/topSize[1]};
  double marginRatioTop[4]; for (auto img : {0,1,2,3}) marginRatioTop[img] = double(fParCvs.marginTop[img])/topSize[int(img<2?0:1)];
  if (nx==1&&ny==1)         for (auto img : {0,1,2,3}) marginRatioTop[img] = double(fParCvs.marginTop[img])/topSize0[int(img<2?0:1)];
  double marginRatioPad[4]; for (auto img : {0,1,2,3}) marginRatioPad[img] = double(fParCvs.marginInn[img])/innSize[int(img<2?0:1)]; 
  double marginRatioEdg[4]; for (auto img : {0,1,2,3}) marginRatioEdg[img] = double(fParCvs.marginTop[img]+fParCvs.marginInn[img])/(fParCvs.marginTop[img]+innSize[int(img<2?0:1)]);
  double marginRatioCur[4] = {padi->GetLeftMargin(), padi->GetRightMargin(), padi->GetBottomMargin(), padi->GetTopMargin()};

  int countPad = 0, ix=0, iy=0;
  if (!(nx==1&&ny==1)&&cidx>0)
  {
    for (iy=0; iy<ny; ++iy) {
      bool breaky = false;
      for (ix=0; ix<nx; ++ix) {
        countPad++;
        if (countPad==cidx) {
          breaky = true;
          break;
        }
      }
      if (breaky)
        break;
    }
  }

  bool isTop = (iy==0);
  bool isBottom = (iy==ny-1);
  bool isLeft = (ix==0);
  bool isRight = ((fParCvs.sidePad[0]==0)&&(ix==nx-1));

  double scaleOffset = double(fullSize[1])/(topSize0[1]+fParCvs.menuSize);
  double titleOffset[4]; for (auto ixyzm : {0,1,2,3}) titleOffset[ixyzm] = fParCvs.titleOffset[ixyzm]*scaleOffset;
  
  hist -> SetStats(fParCvs.setStats);

  for (auto ixyz : {0,1,2}) {
    TAxis *axis;
         if (ixyz==0) axis = (TAxis *) hist -> GetXaxis();
    else if (ixyz==1) axis = (TAxis *) hist -> GetYaxis();
    else if (ixyz==2) axis = (TAxis *) hist -> GetZaxis();
    if (fParCvs.titleAlign[ixyz]==2) axis -> CenterTitle();
    else axis -> CenterTitle(false);
    axis -> RotateTitle(fParCvs.titleRotate[ixyz]);
    axis -> SetTitleOffset(titleOffset[ixyz]);
    axis -> SetTitleSize(fParCvs.titleSize[ixyz]);
    axis -> SetTitleColor(fParCvs.titleColor[ixyz]);
    axis -> SetTitleFont(fParCvs.textFont);
    axis -> SetLabelFont(fParCvs.textFont);
    axis -> SetLabelSize(fParCvs.labelSize[ixyz]);
    axis -> SetLabelOffset(fParCvs.labelOffset[ixyz]);
    axis -> SetTickSize(fParCvs.tickSize[ixyz]);
    int ndivision = fParCvs.axisNdivision[ixyz];
    if (ndivision!=0)
      axis -> SetNdivisions(ndivision);
  }

  hist -> SetContour(fParCvs.numContours);

  if (!(nx==1&&ny==1)&&cidx>0)
  {
    if (isLeft&&isRight) {}
    else  {
      if (ix!=0 && fParCvs.removeInnerPadAxis[1]) {
        hist -> GetYaxis() -> SetLabelOffset(100);
        hist -> GetYaxis() -> SetTitleOffset(100);
      }
      if (isRight&&fParCvs.removeInnerZaxis) {
        TString drawOptionStr = drawOption;
        if (drawOptionStr.Index("z")>=0)
          drawOptionStr.ReplaceAll("z","");
        drawOption = drawOptionStr;
      }
    }
    if (isTop&&isBottom) {}
    else {
      if (iy!=ny-1&&fParCvs.removeInnerPadAxis[0]) {
        hist -> GetXaxis() -> SetLabelOffset(100);
        hist -> GetXaxis() -> SetTitleOffset(100);
      }
      if (iy!=0&&fParCvs.removeInnerMainTitle) {
        hist -> SetTitle("");
        hist -> SetTitle("");
      }
    }
  }

  padi -> Modified();
  padi -> Update();
  if (padi!=nullptr) {
    auto list_primitive = padi -> GetListOfPrimitives();
    auto mainTitle = (TPaveText*) list_primitive -> FindObject("title");
    if (mainTitle!=nullptr)
    {
      mainTitle -> SetTextFont(fParCvs.textFont);
      mainTitle -> SetTextAlign(fParCvs.titleAlign[3]);
      mainTitle -> SetTextSizePixels(fParCvs.titleSize[3]);
      mainTitle -> SetTextColor(fParCvs.titleColor[3]);
      mainTitle -> SetX1NDC(marginRatioCur[0]);
      mainTitle -> SetX2NDC(1.-marginRatioCur[1]);
      mainTitle -> SetY1NDC(1.-marginRatioCur[3]);
      mainTitle -> SetY2NDC(1.);
    }
  }
  padi -> Modified();
  padi -> Update();

  return hist;
}

/// Make z axis
TGaxis *ejungwoo::drawz(TH1* hist, TVirtualPad *vpad, const char *titlez, bool logaxis)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd();

  int nx = 1;
  int ny = 1;
  int cidx = 0;
  TString nameConf = "";
  if (!findConf(padi,nameConf,nx,ny,cidx)) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  getParCanvas(par);

  int countPad = 0, ix=0, iy=0;
  if (!(nx==1&&ny==1)&&cidx>0)
  {
    for (iy=0; iy<ny; ++iy) {
      bool breaky = false;
      for (ix=0; ix<nx; ++ix) {
        countPad++;
        if (countPad==cidx) {
          breaky = true;
          break;
        }
      }
      if (breaky)
        break;
    }
  }

  auto zaxis = new TGaxis(padi->GetUxmax(), hist->GetYaxis()->GetXmin(), padi->GetUxmax(), hist->GetYaxis()->GetXmax(),
      hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax(), fParCvs.axisNdivision[2] ,(logaxis?"+LG":"+L"));

  //cout
  //  << padi->GetUxmax() << " "
  //  << hist->GetYaxis()->GetXmin() << " "
  //  << padi->GetUxmax() << " "
  //  << hist->GetYaxis()->GetXmax() << " "
  //  << hist->GetYaxis()->GetXmin() << " "
  //  << hist->GetYaxis()->GetXmax() << " "
  //  << fParCvs.axisNdivision[2] << endl;

  bool isTop = (iy==0);
  bool isBottom = (iy==ny-1);
  bool isLeft = (ix==0);
  bool isRight = ((fParCvs.sidePad[0]==0)&&(ix==nx-1));

  zaxis -> SetTitle(titlez);
  if (fParCvs.titleAlign[2]==2) zaxis -> CenterTitle();
  else zaxis -> CenterTitle(false);
  zaxis -> SetBit(TAxis::kRotateTitle);
  zaxis -> SetTitleOffset(hist->GetYaxis()->GetTitleOffset());
  //zaxis -> SetMaxDigits(hist->GetYaxis()->GetMaxDigits()); // for higher version of root
  zaxis -> SetTitleSize(fParCvs.titleSize[2]);
  zaxis -> SetTitleColor(fParCvs.titleColor[2]);
  zaxis -> SetLabelSize(fParCvs.labelSize[2]);
  zaxis -> SetLabelOffset(fParCvs.labelOffset[2]);
  zaxis -> SetTickSize(fParCvs.tickSize[2]);
  //zaxis -> SetNdivisions(fParCvs.axisNdivision[2]);
  //zaxis -> SetNdivisions(hist->GetYaxis()->GetNdivisions());
  zaxis -> SetTitleFont(fParCvs.textFont);
  zaxis -> SetLabelFont(fParCvs.textFont);

  if (!(nx==1&&ny==1)&&cidx>0) {
    if (isLeft&&isRight) {}
    else if (ix!=nx-1&&fParCvs.removeInnerZaxis) {
      zaxis -> SetTitleOffset(100);
      zaxis -> SetLabelOffset(100);
    }
  }

  zaxis -> SetTickSize(0);
  padi -> cd();
  zaxis -> Draw();

  return zaxis;
}

TLegend* ejungwoo::getlg(TVirtualPad *vpad)
{
  auto legend = (TLegend *) vpad -> GetPrimitive("elegend");
  if (legend==nullptr) {
    legend = new TLegend();
    legend -> SetName("elegend");
    vpad -> GetListOfPrimitives() -> Add(legend);
  }
  //if (legend->GetListOfPrimitives()!=nullptr) legend -> GetListOfPrimitives() -> ls();
  return legend;
}

TH1 *ejungwoo::draw(TH1 *obj, TVirtualPad *vpad, TString drawOption, TString legendContent)
{
  if (vpad==nullptr)
    vpad = canvas();
  vpad -> cd();
  if (vpad->GetPrimitive("TFrame")!=nullptr) {
    if (drawOption.Index("same")<0) drawOption += "same";
  }
  obj -> Draw(drawOption);
  make(obj,vpad,drawOption);
  if (!legendContent.IsNull()) {
    auto legend = getlg(vpad);
    auto legendTokens = legendContent.Tokenize(";");
    auto legendTitle = ((TObjString *) legendTokens->At(0))->GetString();
    auto legendOption = (legendTokens->GetEntries()<2)?"":((TObjString *) legendTokens->At(1))->GetString();
    legend -> AddEntry(obj,legendTitle,legendOption);
  }
  return obj;
}

TGraph *ejungwoo::draw(TGraph *obj, TVirtualPad *vpad, TString drawOption, TString legendContent)
{
  if (vpad==nullptr)
    vpad = canvas();
  vpad -> cd();
  if (vpad->GetPrimitive("TFrame")!=nullptr) {
    if (drawOption.Index("same")<0) drawOption += "same";
  } else {
    drawOption.ReplaceAll("same","");
    if (drawOption.Index("a")<0) drawOption += "a";
    if (drawOption.Index("p")<0&&drawOption.Index("*")<0) drawOption += "*";
  }
  obj -> Draw(drawOption);
  if (!legendContent.IsNull()) {
    auto legend = getlg(vpad);
    auto legendTokens = legendContent.Tokenize(";");
    auto legendTitle = ((TObjString *) legendTokens->At(0))->GetString();
    auto legendOption = (legendTokens->GetEntries()<2)?"":((TObjString *) legendTokens->At(1))->GetString();
    legend -> AddEntry(obj,legendTitle,legendOption);
  }
  return obj;
}

TF1 *ejungwoo::draw(TF1 *obj, TVirtualPad *vpad, TString drawOption, TString legendContent)
{
  if (vpad==nullptr)
    vpad = canvas();
  vpad -> cd();
  if (vpad->GetPrimitive("TFrame")!=nullptr) {
    if (drawOption.Index("same")<0) drawOption += "same";
  } else {
    drawOption.ReplaceAll("same","");
    //if (drawOption.Index("a")<0) drawOption += "a";
    if (drawOption.Index("l")<0) drawOption += "l";
  }
  obj -> Draw(drawOption);
  if (!legendContent.IsNull()) {
    auto legend = getlg(vpad);
    auto legendTokens = legendContent.Tokenize(";");
    auto legendTitle = ((TObjString *) legendTokens->At(0))->GetString();
    auto legendOption = (legendTokens->GetEntries()<2)?"":((TObjString *) legendTokens->At(1))->GetString();
    legend -> AddEntry(obj,legendTitle,legendOption);
  }
  return obj;
}

void ejungwoo::findxy(TVirtualPad *pad, double &x1, double &y1)
{
  auto x1Frame = 0. + pad -> GetLeftMargin();
  auto x2Frame = 1. - pad -> GetRightMargin();
  auto y1Frame = 0. + pad -> GetBottomMargin();
  auto y2Frame = 1. - pad -> GetTopMargin();
  double x1New = (1.-x1)*x1Frame + x1*(x2Frame);
  double y1New = (1.-y1)*y1Frame + y1*(y2Frame);
  x1 = x1New;
  y1 = y1New;
}

void ejungwoo::findxy(TVirtualPad *pad, double &x1, double &y1, double &dx, double &dy, double &xUnit, double &yUnit)
{
  auto x1Frame = 0. + pad -> GetLeftMargin();
  auto x2Frame = 1. - pad -> GetRightMargin();
  auto y1Frame = 0. + pad -> GetBottomMargin();
  auto y2Frame = 1. - pad -> GetTopMargin();
  xUnit = x2Frame - x1Frame;
  yUnit = y2Frame - y1Frame;
  auto dxNew = dx * xUnit;
  auto dyNew = dy * yUnit;

  double x1New, x2New, y1New, y2New;

  if (x1<0&&x1>-1) {
    x2New = (1.-(-x1))*x2Frame + (-x1)*(x1Frame);
    x1New = x2New - dxNew;
  } else {
    x1New = (1.-x1)*x1Frame + x1*(x2Frame);
    if (x1<0) x1New = x2Frame-dxNew;
    x2New = x1New + dxNew;
  }

  if (y1<0&&y1>-1) {
    y2New = (1.-(-y1))*y2Frame + (-y1)*(y1Frame);
    y1New = y2New - dyNew;
  } else {
    y1New = (1.-y1)*y1Frame + y1*(y2Frame);
    if (y1<0) y1New = y2Frame-dyNew;
    y2New = y1New + dyNew;
  }

  x1 = x1New;
  y1 = y1New;
  dx = dxNew;
  dy = dyNew;
}

/// Make legend fit to the pad.
/// @param legend
/// @param vpad
/// @param x1InRatio  Set x1 (legend left, x-position) by the ratio inside the histogram frame (0~1). If x1InRatio<0 (by default), the legend position is set to stick to the right most frame.
/// @param y1InRatio  Set y1 (legend bottom, y-position) by the ratio inside the histogram frame (0~1). If y1InRatio<0 (by default), the legend position is set to stick to the top most frame.
/// @param dxInRatio  Set dy (legend width) by the ratio inside the histogram frame (0~1).  If dxInRatio<=0 (by default), the legend width is calcuated from the text length.
/// @param dyInRatio  Set dy (legend height) by the ratio inside the histogram frame (0~1). If dyInRatio<=0 (by default), the legend height is calcuated from the number of legend entries.
/// @param marginObj  Set ratio of the space occupied by object in the legend (compared to the descriptions). By default, it is >0.25
TLegend *ejungwoo::make(TLegend *legend, TVirtualPad *vpad, double x1InRatio=-1, double y1InRatio=-1, double dxInRatio=0, double dyInRatio=0, double marginObj=-1, int tf=-1, double tz=-1, int ta=-1, int fs=-1, int fc=-1, int bs=-1)
{
  auto pad = (TPad *) vpad;

  double fLegendTextScale = 0.028;
  double fWidthPerLengthLegend = .8;
  double fWidthDefaultLegend = fWidthPerLengthLegend*2.5;
  double fLegendHeightPerEntry = 2.5;
  double xOffset = 0;
  double yOffset = 0;

  int dummy;
  TString nameConf = "";
  if (!findConf(pad,nameConf,dummy,dummy,dummy)) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  getParCanvas(par);

  double dxTextMax = 0.;
  bool nullLegend = true;
  auto dxNorm = TLatex(0,0,"0").GetXsize();
  TIter next_entry(legend->GetListOfPrimitives());
  while (TLegendEntry *le=(TLegendEntry*)next_entry()) {
    TString legendString = le->GetLabel();
    auto dxText = TLatex(0,0,legendString).GetXsize()/dxNorm;
    if (dxTextMax<dxText) dxTextMax = dxText;
    if (le->GetObject()!=nullptr) nullLegend = false;
  }

  if (nullLegend) {
    fWidthDefaultLegend = 0;
  }

  auto numColumns = legend -> GetNColumns();
  auto numRows = legend -> GetNRows();
  if (numColumns>1) {
    numRows = ceil(double(numRows) / numColumns);
  }

  if (dxInRatio<=0) dxInRatio = (fWidthDefaultLegend + dxTextMax * fWidthPerLengthLegend) * fLegendTextScale * numColumns;
  if (dxInRatio>1) dxInRatio = 1;
  if (dyInRatio<=0) dyInRatio = (numRows * fLegendHeightPerEntry) * fLegendTextScale;

  double x1Legend = x1InRatio;
  double y1Legend = y1InRatio;
  double dxLegend = dxInRatio;
  double dyLegend = dyInRatio;
  double xUnit, yUnit;
  findxy(pad, x1Legend, y1Legend, dxLegend, dyLegend, xUnit, yUnit);
  double x2Legend = x1Legend + dxLegend;
  double y2Legend = y1Legend + dyLegend;
  //cout << x1Legend << " " << x2Legend << " " << y1Legend << " " << y2Legend << endl;

  double marginLegend = marginObj;
  if (nullLegend)
    marginLegend = 0;
  else if (marginLegend<0||marginLegend>1) {
    marginLegend = 1.-((dxTextMax*fWidthPerLengthLegend*fLegendTextScale)/(dxLegend/xUnit));
    //if (marginLegend<0.25) marginLegend = 0.25;
    if (marginLegend>.5) marginLegend = .5;
    if (marginLegend<.1) marginLegend = 0.1;
  }

  if (fs<0) fs=0;
  if (fc<0) fc=kWhite;
  if (bs<0) bs=0;

  legend -> SetX1(x1Legend + xOffset*xUnit);
  legend -> SetX2(x2Legend + xOffset*xUnit);
  legend -> SetY1(y1Legend + yOffset*yUnit);
  legend -> SetY2(y2Legend + yOffset*yUnit);

  //legend -> SetX1NDC(x1Legend + xOffset*xUnit);
  //legend -> SetX2NDC(x2Legend + xOffset*xUnit);
  //legend -> SetY1NDC(y1Legend + yOffset*yUnit);
  //legend -> SetY2NDC(y2Legend + yOffset*yUnit);

  legend -> SetFillStyle(fs);
  legend -> SetFillColor(fc);
  legend -> SetBorderSize(bs);

  if (tf>0) legend -> SetTextFont(tf);
  else legend -> SetTextFont(fParCvs.textFont);

  if (tz>0) legend -> SetTextSize(tz);
  else if (fParCvs.textSize>0)
    legend -> SetTextSize(fParCvs.textSize);

  if (ta>0) legend -> SetTextAlign(ta);
  legend -> SetMargin(marginLegend);

  return legend;
}

/// Make legend fit to the pad and draw.
TLegend* ejungwoo::draw(TLegend* legend, TVirtualPad *vpad, double x1InRatio, double y1InRatio, double dxInRatio, double dyInRatio, double marginObj, int tf, double tz, int ta, int fs, int fc, int bs)
{
  vpad -> cd();
  make(legend,vpad,x1InRatio,y1InRatio,dxInRatio,dyInRatio,marginObj,tf,tz,ta,fs,fc,bs);
  legend -> Draw();
  return legend;
}

void ejungwoo::drawef(TVirtualPad *vpad) {
  vpad -> cd();
  vpad -> SetFrameLineColor(0);
  auto histef = vpad -> DrawFrame(0,10,0,10);
  //auto histef = new TH1D(vpad->GetName(),"",100,0,10);
  histef -> SetStats(0);
  histef -> GetXaxis() -> SetLabelOffset(100);
  histef -> GetYaxis() -> SetLabelOffset(100);
  histef -> GetXaxis() -> SetNdivisions(0);
  histef -> GetYaxis() -> SetNdivisions(0);
  histef -> Draw("a");
}

/**
 * Set attribute of (line, marker and filling) to [histogram, graph, function, marker and line] from given nameConf
 *
 * ll : line
 * mm : marker
 * ff : fill
 * tt : fill
 *
 * xy : x is one of l(line), m(marker), f(fill), t(text). y is one of s(style), z(size), c(color).
 */
TObject *ejungwoo::att(TObject *obj, int idx, const char *nameConf)
{
  TString sNameConf = nameConf;
  int olz=idx, ols=idx, olc=idx, omz=idx, oms=idx, omc=idx, ofz=idx, ofs=idx, ofc=idx, otz=idx, otf=idx, otc=idx;
  double xmz=1, xtz=1;

  if (sNameConf.Index(";")>=0) {
    auto elements = sNameConf.Tokenize(";");
    if (sNameConf.Index(";")!=0) {
      sNameConf = ((TObjString *) elements->At(0))->GetString();
      elements -> RemoveAt(0);
      elements -> Compress();
    }
    else {
      sNameConf = "";
    }

    if (!exist_key(elements,"ll")) {
      if(exist_key(elements,"lz")) olz=idx; else if(exist_keyv(elements,"lz")) olz=get_valuei(elements,"lz"); else olz=-1; 
      if(exist_key(elements,"ls")) olz=idx; else if(exist_keyv(elements,"ls")) ols=get_valuei(elements,"ls"); else ols=-1; 
      if(exist_key(elements,"lc")) olz=idx; else if(exist_keyv(elements,"lc")) olc=get_valuei(elements,"lc"); else olc=-1; 
    }
    if (!exist_key(elements,"mm")) {
      if(exist_key(elements,"mz")) omz=idx; else if(exist_keyv(elements,"mz")) omz=get_valuei(elements,"mz"); else omz=-1; 
      if(exist_key(elements,"ms")) omz=idx; else if(exist_keyv(elements,"ms")) oms=get_valuei(elements,"ms"); else oms=-1; 
      if(exist_key(elements,"mc")) omz=idx; else if(exist_keyv(elements,"mc")) omc=get_valuei(elements,"mc"); else omc=-1; 
    }
    if (!exist_key(elements,"ff")) {
      if(exist_key(elements,"fz")) ofz=idx; else if(exist_keyv(elements,"fz")) ofz=get_valuei(elements,"fz"); else ofz=-1; 
      if(exist_key(elements,"fs")) ofz=idx; else if(exist_keyv(elements,"fs")) ofs=get_valuei(elements,"fs"); else ofs=-1; 
      if(exist_key(elements,"fc")) ofz=idx; else if(exist_keyv(elements,"fc")) ofc=get_valuei(elements,"fc"); else ofc=-1; 
    }
    if (!exist_key(elements,"tm")) {
      if(exist_key(elements,"tz")) otz=idx; else if(exist_keyv(elements,"tz")) otz=get_valuei(elements,"tz"); else otz=-1; 
      if(exist_key(elements,"ts")) otz=idx; else if(exist_keyv(elements,"ts")) otf=get_valuei(elements,"ts"); else otf=-1; 
      if(exist_key(elements,"tf")) otz=idx; else if(exist_keyv(elements,"tf")) otf=get_valuei(elements,"tf"); else otf=-1; 
      if(exist_key(elements,"tc")) otz=idx; else if(exist_keyv(elements,"tc")) otc=get_valuei(elements,"tc"); else otc=-1; 
    }
    if (exist_keyv(elements,"xmz")) xmz = get_valued(elements,"xmz");
    if (exist_keyv(elements,"xtz")) xtz = get_valued(elements,"xtz");
    if (exist_key(elements,"cc")) { olc=idx; omc=idx; ofc=idx; otc=idx; }
  }

  if (sNameConf.IsNull()) sNameConf = fNameAttributeConf;
  auto par = conf(sNameConf);
  getParAttribute(par);

  if (obj->InheritsFrom(TH1::Class())) {
    auto hist = (TH1 *) obj;
    if (olz>=0) hist -> SetLineWidth(fParAtt.lineWidth[olz]);
    if (ols>=0) hist -> SetLineStyle(fParAtt.lineStyle[ols]);
    if (olc>=0) hist -> SetLineColor(fParAtt.lineColor[olc]);
    if (omz>=0) hist -> SetMarkerSize(fParAtt.markerSize[omz]*xmz);
    if (oms>=0) hist -> SetMarkerStyle(fParAtt.markerStyle[oms]);
    if (omc>=0) hist -> SetMarkerColor(fParAtt.markerColor[omc]);
    if (ofs>=0) hist -> SetFillStyle(fParAtt.fillStyle[ofs]);
    if (ofc>=0) hist -> SetFillColor(fParAtt.fillColor[ofc]);
  }

  else if (obj->InheritsFrom(TGraph::Class())||obj->InheritsFrom(TGraph2D::Class())) {
    auto graph = (TGraph *) obj;
    if (olz>=0) graph -> SetLineWidth(fParAtt.lineWidth[olz]);
    if (ols>=0) graph -> SetLineStyle(fParAtt.lineStyle[ols]);
    if (olc>=0) graph -> SetLineColor(fParAtt.lineColor[olc]);
    if (omz>=0) graph -> SetMarkerSize(fParAtt.markerSize[omz]*xmz);
    if (oms>=0) graph -> SetMarkerStyle(fParAtt.markerStyle[oms]);
    if (omc>=0) graph -> SetMarkerColor(fParAtt.markerColor[omc]);
    if (ofs>=0) graph -> SetFillStyle(fParAtt.fillStyle[ofs]);
    if (ofc>=0) graph -> SetFillColor(fParAtt.fillColor[ofc]);
  }

  else if (obj->InheritsFrom(TF1::Class())) {
    auto f1 = (TF1 *) obj;
    if (olz>=0) f1 -> SetLineWidth(fParAtt.lineWidth[olz]);
    if (ols>=0) f1 -> SetLineStyle(fParAtt.lineStyle[ols]);
    if (olc>=0) f1 -> SetLineColor(fParAtt.lineColor[olc]);
    if (ofs>=0) f1 -> SetFillStyle(fParAtt.fillStyle[ofs]);
    if (ofc>=0) f1 -> SetFillColor(fParAtt.fillColor[ofc]);
  }

  else if (obj->InheritsFrom(TMarker::Class())) {
    auto marker = (TMarker *) obj;
    if (omz>=0) marker -> SetMarkerSize(fParAtt.markerSize[omz]*xmz);
    if (oms>=0) marker -> SetMarkerStyle(fParAtt.markerStyle[oms]);
    if (omc>=0) marker -> SetMarkerColor(fParAtt.markerColor[omc]);
  }

  else if (obj->InheritsFrom(TLine::Class())) {
    auto line = (TLine *) obj;
    if (ols>=0) line -> SetLineWidth(fParAtt.lineWidth[ols]);
    if (olz>=0) line -> SetLineStyle(fParAtt.lineStyle[olz]);
    if (olc>=0) line -> SetLineColor(fParAtt.lineColor[olc]);
  }

  else if (obj->InheritsFrom(TText::Class())) {
    auto text = (TText *) obj;
    text -> SetTextAlign(fParAtt.textAlign);
    if (otz>=0) text -> SetTextSize(fParAtt.textSize*xtz);
    if (otf>=0) text -> SetTextFont(fParAtt.textFont);
    if (otc>=0) if (otc<fParAtt.textColor.size()) text -> SetTextColor(fParAtt.textColor[otc]);
  }

  return obj;
}

/// Remove trailing zeros and return as a string.
TString ejungwoo::toString(double value, int ndigit)
{
  TString vstring = Form("%f",value);
  if (ndigit==0) vstring = Form("%.0f",value);
  if (ndigit==1) vstring = Form("%.1f",value);
  if (ndigit==2) vstring = Form("%.2f",value);
  if (ndigit==3) vstring = Form("%.3f",value);
  if (ndigit==4) vstring = Form("%.4f",value);
  if (ndigit==5) vstring = Form("%.5f",value);

  auto posdot = vstring.Index(".");
  auto poslast = vstring.Sizeof()-2;
  auto lstring = TString(vstring(poslast));
  if (posdot>=0&&posdot<poslast) {
    while (lstring=="0") {
      vstring.Remove(poslast);
      poslast = vstring.Sizeof()-2;
      lstring = TString(vstring(poslast));
    }
    if (lstring==".")
      vstring.Remove(poslast);
  }

  return vstring;
}

TString ejungwoo::get_value(TString array, TString skey, bool exact=false, TString token) {
  auto elements = array.Tokenize(token);
  return get_value(elements,skey,exact);
}

TString ejungwoo::get_value(TObjArray *array, TString skey, bool exact)
{
  TIter next(array);
  while (auto elemento = (TObjString*) next()) {
    TString element = elemento -> GetString();
    if (exact) {
      if (element==skey)
        return element;
    }
    else {
      if (element.Index(skey)==0) {
        int poseq = element.Index("=");
        if (int(skey.Sizeof()-1)!=int(poseq)) 
          continue;
        else if (poseq>=0) {
          auto fulle = element.Sizeof();
          TString value = element(poseq+1,fulle-poseq-2);
          return value;
        }
        return element;
      }
    }
  }
  return "";
}

bool ejungwoo::exist_key(TObjArray *elements, TString skey) {
  auto value = get_value(elements,skey,true);
  if (value.IsNull())
    return false;
  return true;
}
bool ejungwoo::exist_keyv(TObjArray *elements, TString skey) {
  auto value = get_value(elements,skey,false);
  if (value.IsNull())
    return false;
  return true;
}
double ejungwoo::get_valued(TObjArray *elements, TString skey) {
  auto value = get_value(elements,skey,false);
  if (value.IsNull())
    return -999;
  return value.Atof();
}
int ejungwoo::get_valuei(TObjArray *elements, TString skey) {
  auto value = get_value(elements,skey,false);
  if (value.IsNull())
    return -999;
  return value.Atoi();
}

void ejungwoo::colorwheel() {
  TColorWheel *w = new TColorWheel();
  w -> SetCanvas(new TCanvas("cw","",800,825));
  w -> Draw();
}

void ejungwoo::markers() {
  auto cvs = new TCanvas("markers","",600,400);
  cvs->SetMargin(0.02,0.02,0.02,0.02);
  auto hist = new TH2D("hist","",100,0.2,10.8,100,0.1,5.6);
  hist -> SetStats(0);
  hist -> SetStats(0);
  hist -> GetXaxis() -> SetLabelOffset(100);
  hist -> GetYaxis() -> SetLabelOffset(100);
  hist -> GetXaxis() -> SetNdivisions(0);
  hist -> GetYaxis() -> SetNdivisions(0);
  hist -> Draw();
  int i = 0;
  for (auto y=5; y>=1; --y) {
    for (auto x=1; x<=10; ++x) {
      if (i==0) { i++; continue; }
      auto m = new TMarker(x,y,i);
      m -> SetMarkerSize(3.5);
      auto t = new TText(x,y-0.42,Form("%d",i));
      t -> SetTextSize(0.035);
      t -> SetTextAlign(22);
      if ((i>=20&&i<=29)||i==33||i==34) t -> SetTextColor(kBlue);
      if ((i>=24&&i<=28)||i==30||i==32) t -> SetTextColor(kRed);
      if (i>=9&&i<=19) {
        t -> SetTextColor(kGray);
        m -> SetMarkerColor(kGray);
      }
      m -> Draw();
      t -> Draw();
      i++;
    }
  }
}

TF1 *ejungwoo::fitg(TH1 *hist, double distSigma, Option_t *opt)
{
  if (distSigma<=0)
    distSigma = 2.5;

  double binmax = hist -> GetMaximumBin();
  double max = hist -> GetBinContent(binmax);
  double xmax = hist -> GetXaxis() -> GetBinCenter(binmax);
  double xerr = hist -> GetStdDev();
  double frange1, frange2, dfrangeCut;
  auto nbins = hist -> GetXaxis() -> GetNbins();
  auto hrange1 = hist -> GetXaxis() -> GetBinLowEdge(1);
  auto hrange2 = hist -> GetXaxis() -> GetBinUpEdge(nbins);
  bool outofrange;

  auto fit = new TF1(Form("%s_fitg",hist -> GetName()),"gaus(0)",xmax-xerr*distSigma,xmax+xerr*distSigma);
  fit -> SetParameters(max,xmax,xerr);

  outofrange = false;
  fit -> GetRange(frange1, frange2);
  if (frange1 < hrange1) {
    if (frange2 < hrange1) outofrange = true;
    else if (frange2 - hrange1 < dfrangeCut) outofrange = true;
    else frange1 = hrange1;
  } else if (frange2 > hrange2) {
    if (frange1 > hrange2) outofrange = true;
    else if (hrange2 - frange1 < dfrangeCut) outofrange = true;
    else frange2 = hrange2;
  } fit -> SetRange(frange1, frange2);

  if (outofrange) {
    fit -> SetParameters(0,0,0);
    return fit;
  }

  hist -> Fit(fit,opt);
  xmax=fit -> GetParameter(1);
  xerr=fit -> GetParameter(2);
  fit -> SetRange(xmax-distSigma*xerr,xmax+distSigma*xerr);

  outofrange = false;
  fit -> GetRange(frange1, frange2);
  if (frange1 < hrange1) {
    if (frange2 < hrange1) outofrange = true;
    else if (frange2 - hrange1 < dfrangeCut) outofrange = true;
    else frange1 = hrange1;
  } else if (frange2 > hrange2) {
    if (frange1 > hrange2) outofrange = true;
    else if (hrange2 - frange1 < dfrangeCut) outofrange = true;
    else frange2 = hrange2;
  } fit -> SetRange(frange1, frange2);

  if (outofrange) {
    fit  ->  SetParameters(0,0,0);
    return fit;
  }

  int fitcount = 1;
  double xerr2 = 0.;
  while (fitcount<20&&TMath::Abs(xerr-xerr2)/xerr>0.2) {
    xerr2=xerr;
    fit -> SetRange(xmax-distSigma*xerr2,xmax+distSigma*xerr2);
    fit -> GetRange(frange1, frange2);

    nbins = hist -> GetXaxis()->GetNbins();
    hrange1 = hist -> GetXaxis()->GetBinLowEdge(1);
    hrange2 = hist -> GetXaxis()->GetBinUpEdge(nbins);

    outofrange = false;
    fit -> GetRange(frange1, frange2);
    if (frange1 < hrange1) {
      if (frange2 < hrange1) outofrange = true;
      else if (frange2 - hrange1 < dfrangeCut) outofrange = true;
      else frange1 = hrange1;
    } else if (frange2 > hrange2) {
      if (frange1 > hrange2) outofrange = true;
      else if (hrange2 - frange1 < dfrangeCut) outofrange = true;
      else frange2 = hrange2;
    } fit -> SetRange(frange1, frange2);
    if (outofrange)
      break;

    hist -> Fit(fit,opt);
    xmax=fit -> GetParameter(1);
    xerr=fit -> GetParameter(2);
    ++fitcount;
  }

  return fit;
}

/// The first two points of graph_hand_pointed should be the two opposite corner of histogram (x1,y1) and (x2,y2);
/// reconstructed  points will be written out if saveName is given
TGraph *ejungwoo::extractPoints(TGraph *graph_hand_pointed, double x1, double x2, double y1, double y2, TString saveName)
{
  ofstream outFile;
  if (!saveName.IsNull()) {
    if (saveName==".")
      saveName = Form("%s.dat",graph_hand_pointed -> GetName());
    if (fVerboseLevel>=verboseLevel::kNormal)
      cout_info << "writting reconstructed points to " << saveName << endl;
    outFile.open(saveName);
  }
  auto graphReco = new TGraph();
  int numPoints = graph_hand_pointed -> GetN();
  double x0, y0, x1LC, x2LC, y1LC, y2LC;
  graph_hand_pointed -> GetPoint(0,x1LC,y1LC);
  graph_hand_pointed -> GetPoint(1,x2LC,y2LC);
  for (auto iPoint=2; iPoint<numPoints; ++iPoint) {
    graph_hand_pointed -> GetPoint(iPoint,x0,y0);
    double dx1 = x0 - x1LC;
    double dx2 = x2LC - x0;
    double xReco = (dx1*x2+dx2*x1) / (dx1+dx2);
    double dy1 = y0 - y1LC;
    double dy2 = y2LC - y0;
    double yReco = (dy1*y2+dy2*y1) / (dy1+dy2);
    graphReco -> GetPoint(graphReco->GetN(),xReco,yReco);
    if (!saveName.IsNull())
      outFile << xReco << " " << yReco << endl;
  }
  return graphReco;
}

TDatabasePDG *ejungwoo::particleDB()
{
  TDatabasePDG *db = TDatabasePDG::Instance();
  if (db->GetParticle("deuteron")==nullptr) {
    db -> AddParticle("deuteron","2H"  , 1.87561 ,1,0, 3,"Ion",1000010020);
    db -> AddParticle("triton"  ,"3H"  , 2.80892 ,1,0, 3,"Ion",1000010030);
    db -> AddParticle("he3"     ,"3He" , 2.80839 ,1,0, 6,"Ion",1000020030);
    db -> AddParticle("he4"     ,"4He" , 3.72738 ,1,0, 6,"Ion",1000020040);
    db -> AddParticle("li6"     ,"6Li" , 5.6     ,1,0, 9,"Ion",1000030060);
    db -> AddParticle("li7"     ,"7Li" , 6.5     ,1,0, 9,"Ion",1000030070);
    db -> AddParticle("be7"     ,"7Be" , 6.5     ,1,0,12,"Ion",1000040070);
    db -> AddParticle("be9"     ,"9Be" , 8.4     ,1,0,12,"Ion",1000040090);
    db -> AddParticle("be10"    ,"10Be", 9.3     ,1,0,12,"Ion",1000040100);
    db -> AddParticle("bo10"    ,"10Bo", 9.3     ,1,0,15,"Ion",1000050100);
    db -> AddParticle("bo11"    ,"11Bo",10.2     ,1,0,15,"Ion",1000050110);
    db -> AddParticle("c11"     ,"11C" ,10.2     ,1,0,18,"Ion",1000060110);
    db -> AddParticle("c12"     ,"12C" ,11.17793 ,1,0,18,"Ion",1000060120);
    db -> AddParticle("c13"     ,"13C" ,12.11255 ,1,0,18,"Ion",1000060130);
    db -> AddParticle("c14"     ,"14C" ,13.04394 ,1,0,18,"Ion",1000060140);
    db -> AddParticle("n13"     ,"13N" ,12.1     ,1,0,21,"Ion",1000070130);
    db -> AddParticle("n14"     ,"14N" ,13.0     ,1,0,21,"Ion",1000070140);
    db -> AddParticle("n15"     ,"15N" ,14.0     ,1,0,21,"Ion",1000070150);
    db -> AddParticle("o16"     ,"16O" ,14.89917 ,1,0,24,"Ion",1000080160);
    db -> AddParticle("o17"     ,"17O" ,15.83459 ,1,0,24,"Ion",1000080170);
    db -> AddParticle("o18"     ,"18O" ,16.76611 ,1,0,24,"Ion",1000080180);
  }
  return db;
}
TParticlePDG *ejungwoo::particle(TString name) { return (particleDB() -> GetParticle(name)); }
TParticlePDG *ejungwoo::particle(int pdg)      { return (particleDB() -> GetParticle(pdg)); }

#endif
