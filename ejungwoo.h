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

  class parCanvas {
    public:
      int font, menuSize=25, padSize[2], marginInn[4], marginTop[4], marginDiv[2] ,axisNdivision[3], titleAlign[4], titleColor[4], numContours=200;
      bool removeInnerPadAxis[2], removeInnerMainTitle, removeInnerZaxis, titleRotate[4], setStats=false, setDivIndexYX=false;
      double tickSize[3], labelSize[3], labelOffset[3], titleSize[4], titleOffset[4];
  };
  class parAttribute {
    public:
      vector<int> lineStyle, lineWidth, lineColor, markerStyle, markerColor, fillColor, fillStyle, textColor;
      vector<double> markerSize;
      int textFont=133, textSize=25, textAlign=22;
  };

  enum verboseLevel : int { kQuiet=0, kNormal=1, kAll=2 };
  verboseLevel fVerboseLevel = kNormal;
  void verbose(int level) {
    fVerboseLevel = (verboseLevel)level;
    cout << "verbose is kQuiet(0) " << endl;
    cout << "verbose is kNormal(1)" << endl;
    cout << "verbose is kAll(2)   " << endl;
  }

  int fCanvasIndex = 0;
  TObjArray *fCanvasArray = nullptr;
  TClonesArray *fParameterArray = nullptr;
  const char *fNameAttributeConf = "att1";
  const char *fNameCanvasConf = "single";
  parAttribute fParAtt;
  parCanvas fParCvs;

  parContainer *conf(const char *nameConf);
  bool findConf(TVirtualPad *vpad, TString &nameConfCvs, int &nx, int &ny, int &cidx);
  void setCanvasPar(parContainer *par);
  void setAttributePar(parContainer *par);
  TObject *att(TObject *obj, int idx=0, const char *nameConf="");

  void        findxy(TVirtualPad *pad, double &x1, double &y1);
  void        findxy(TVirtualPad *pad, double &x1, double &y1, double &dx, double &dy, double &xUnit, double &yUnit);
  TCanvas*    canvas(const char *nameCvs, int dx, int dy, TArrayD poss);
  TCanvas*    canvas(const char *nameCvs, int nx, int ny, const char *nameConf="");
  TCanvas*    canvas(const char *nameCvs="", const char *nameConf="") { return canvas(nameCvs, 1, 1, nameConf); }
  TPad*       innerpad(TVirtualPad *vpad, double xr1, double yr1, double xr2, double yr2);
  TH1*        make (TH1* hist, TVirtualPad *vpad, const char *drawOption="");
  TH1*        draw (TH1* hist, TVirtualPad *vpad, const char *drawOption="");
  TGaxis*     drawz(TH1* hist, TVirtualPad *vpad, const char *titlez="");
  TLegend*    make(TLegend* legend, TVirtualPad *vpad, double x1ir=-1, double y1ir=-1, double dxir=0, double dyir=0, double marginObj=-1, int fillStyle=0, int fillColor=kWhite, int borderSize=0);
  TLegend*    draw(TLegend* legend, TVirtualPad *vpad, double x1ir=-1, double y1ir=-1, double dxir=0, double dyir=0, double marginObj=-1, int fillStyle=0, int fillColor=kWhite, int borderSize=0);
  TPaveText*  newpt(TString content, TVirtualPad *vpad, double x1ir=-1, double y1ir=-1, double dxir=0, double dyir=0, int tf=133, int ts=25, int ta=12, int fc=0, int fs=0, int bs=0);
  TH1*        project(TTree *tree, TString formula, TCut cut="", TString name="", TString title="", int nx=0, int x1=0, int x2=0, int ny=0, int y1=0, int y2=0);
  TPad*       redraw(TPad *pad);

  void saveRoot(TObject *obj, TString nameFile="", TString nameVersion="", bool savePrimitives=false, bool simplifyNames=false);
  void saveRoot(TString nameVersion="", bool savePrimitives=false, bool simplifyNames=false);
  void savePDF(TString nameVersion="");
  void savePNG(TString nameVersion="");
  void saveAll(TString nameVersion="");
  void write(TObject *obj);

  TString tok(TString line, TString token, int i);
  TString tok(TObjArray *line, int i);
  TF1 *fitg(TH1 *hist, double distSigma=2.5, Option_t *opt="RQ0");
  TGraph *extractPoints(TGraph *graph_hand_pointed, double x1, double x2, double y1, double y2, TString saveName="");

  void colorwheel();
  void markers();
  TString toString(double value);
  TDatabasePDG *particleDB();
  TParticlePDG *particle(TString name);
  TParticlePDG *particle(int pdg);
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
    const char *fExpression = "";
    const char *fSelection = "";

    enum bntype : int { kEqualSpacing=0, kContinuousSpacing=1, kDescreteSpacing=2 };
    bntype fSpacing = kEqualSpacing;
    TArrayD fBinArray;

  public:
    binning() : fNbins(0), fMin(0), fMax(0) {}
    binning(int nbins, double min, double max, const char *ttl="", const char *expression="", const char *selection="") : fNbins(nbins), fMin(min), fMax(max), fTitle(ttl), fExpression(expression), fSelection(selection) { init(); }
    binning(int nbins, TArrayD array, const char *ttl="", const char *expression="", const char *selection="") : fTitle(ttl), fExpression(expression), fSelection(selection) { initArray(nbins,array); }
    binning(           TArrayD array, const char *ttl="", const char *expression="", const char *selection="") : binning(0, array, ttl, expression, selection) {}
    binning(TTree *tree, const char *branchName) : binning() { make(tree, branchName); }
    binning(TH1 *hist, int axis_123=1) : binning() { make(hist, axis_123); }
    ~binning() {}

    void init();
    void initArray(int nbins, TArrayD array);
    void make(TH1 *hist, int axis_123=1);
    void make(TTree *tree, const char* branchName);
    void set(double nbins, double min, double max, const char *ttl="", const char *expression="", const char *selection="");
    void set(int nbins, TArrayD array, const char *ttl="", const char *expression="", const char *selection="");
    void set(TArrayD array, const char *ttl="", const char *expression="", const char *selection="");
    void setN         (double nbins)             { fNbins = nbins; fWidth = (fMax-fMin)/fNbins; }
    void setMin       (double min)               { fMin = min; fWidth = (fMax-fMin)/fNbins; }
    void setMax       (double max)               { fMax = max; fWidth = (fMax-fMin)/fNbins; }
    void setW         (double w)                 { fWidth = w; fNbins = int((fMax-fMin)/fWidth); }
    void setTitle     (const char *ttl)          { fTitle = ttl; }
    void setExpression(const char *expression)   { fExpression = expression; }
    void setSelection (const char *selection)    { fSelection = selection; }

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
    const char *getExpression()   const { return fExpression; }
    const char *getSelection()    const { return fSelection; }
    TArrayD getBinArray()         const { return fBinArray; }
    bool   null()                 const { return isNull(); }
    int    n()                    const { return getN(); }
    double min()                  const { return getMin(); }
    double max()                  const { return getMax(); }
    double width(int bin=-1)      const { return getW(bin); }
    const char *title()           const { return getTitle(); }
    const char *expr()            const { return getExpression(); }
    const char *selection()       const { return getSelection(); }

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

    double lowEdge(int bin)       const;
    double highEdge(int bin)      const;
    int    findBin(double val)    const;
    double getFullWidth()         const { return (fMax - fMin); }
    double getFullCenter()        const { return (fMax + fMin)/2.; }
    double getCenter(int bin)     const { return (fMin + (bin-.5)*fWidth); }
    bool   isInside(double value) const { if(value>fMin && value<fMax) return true; return false; }
    double xByRatio(double ratio) const { return ((1.-ratio)*fMin + ratio*fMax); }
    double fullw()                const { return getFullWidth(); }
    double fullc()                const { return getFullCenter(); }
    double center(int bin)        const { return getCenter(bin); }
    bool   inside(double value)   const { return isInside(value); }
    double xratio(double ratio)   const { return xByRatio(ratio); }

    TCut cut() const { return TCut(Form("%s>=%f&&%s<%f",fExpression,fMin,fExpression,fMax)); }
    const char *nnmm() { return Form("%s_%d_%.2f_%.2f",fExpression,fNbins,fMin,fMax); }
    TH1D *newHist(const char *name, const char *title="");
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
  else                  axis = hist -> GetXaxis();
  fNbins = axis -> GetNbins();
  fMin = axis -> GetBinLowEdge(1);
  fMax = axis -> GetBinUpEdge(fNbins);
  init();
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
  fExpression = branchName;
  init();
}

void ejungwoo::binning::set(double nbins, double min, double max, const char *ttl, const char *expression, const char *selection) {
  fNbins = nbins;
  fMin = min;
  fMax = max;
  fWidth = (fMax-fMin)/fNbins;
  init();
  if (strcmp(ttl,"")!=0) fTitle = ttl;
  if (strcmp(expression,"")!=0) fExpression = expression;
  if (strcmp(selection,"")!=0) fSelection = selection;
}

void ejungwoo::binning::set(int nbins, TArrayD array, const char *ttl, const char *expression, const char *selection) {
  initArray(nbins,array);
  if (strcmp(ttl,"")!=0) fTitle = ttl;
  if (strcmp(expression,"")!=0) fExpression = expression;
  if (strcmp(selection,"")!=0) fSelection = selection;
}

void ejungwoo::binning::set(TArrayD array, const char *ttl, const char *expression, const char *selection) {
  initArray(0,array);
  if (strcmp(ttl,"")!=0) fTitle = ttl;
  if (strcmp(expression,"")!=0) fExpression = expression;
  if (strcmp(selection,"")!=0) fSelection = selection;
}

double ejungwoo::binning::getW(int bin) const {
  if (bin==0) return (fMax - fMin);
  if (bin<0) bin = fBindex;
       if (esp()) return fWidth;
  else if (csp()) return (fBinArray[bin] - fBinArray[bin-1]);
  else if (dsp()) return (fBinArray[bin*2-1] - fBinArray[bin*2-2]);
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

TH1D *ejungwoo::binning::newHist(const char *name, const char *title) {
  auto titlexy = Form("%s;%s;",title,fTitle);
  TH1D *hist;
  if (esp()) hist = new TH1D(name,titlexy,fNbins,fMin,fMax);
  else if (csp()) hist = new TH1D(name,titlexy,fNbins,fBinArray.GetArray());
  else if (dsp()) hist = new TH1D(name,titlexy,2*fNbins-1,fBinArray.GetArray());
  return hist;
}

TH1D *ejungwoo::binning::project(TTree *tree, TCut selection, TString name, TString title) {
  if (name.IsNull())
    name = getExpression();
  auto hist = newHist(name,title);
  tree -> Project(name,fExpression,selection);
  return hist;
}

TString ejungwoo::binning::printi(bool pout) const {
  TString line = Form("%s(%d)=%s~%s",fExpression,ii(),toString(min()).Data(),toString(max()).Data());
  if (pout)
    cout << line << endl;
  return line;
}

TString ejungwoo::binning::print(bool pout) const {
  TString minString = ejungwoo::toString(fMin);
  TString maxString = ejungwoo::toString(fMax);
  TString wString = ejungwoo::toString(fWidth);
  TString line = Form("[%s,%s] %d,%s,%s",fExpression,fTitle,fNbins,minString.Data(),maxString.Data());
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
  if (esp()) bncopy = binning(fNbins, fMin, fMax, fTitle, fExpression, fSelection);
  else       bncopy = binning(fNbins, fBinArray, fTitle, fExpression, fSelection);
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
      auto tokens = TString(val).Tokenize(";");
      const char *title = ((TObjString *) tokens->At(0))->GetString();
      return title;
    }

  public :
    binning2(binning bnX, binning bnY) : fbnX(bnX), fbnY(bnY) {}
    binning2(int nx, double x1, double x2, int ny, double y1, double y2) : fbnX(binning(nx,x1,x2)), fbnY(binning(ny,y1,y2)) {}
    binning2(double x1, double x2, double y1, double y2) : fbnX(binning(100,x1,x2)), fbnY(binning(100,y1,y2)) {}
    ~binning2() {}

    TH2D *newHist(const char *nameHist, const char *titleHist="") {
      auto titlexy = Form("%s;%s;%s;",titleHist,justTitle(fbnX.getTitle()),justTitle(fbnY.getTitle()));
      TH2D *hist = nullptr;
           if ( fbnX.esp()&& fbnY.esp()) hist = new TH2D(nameHist,titlexy, fbnX.getN(),fbnX.getMin(),fbnX.getMax(), fbnY.getN(),fbnY.getMin(),fbnY.getMax());
      else if ( fbnX.esp()&&!fbnY.esp()) hist = new TH2D(nameHist,titlexy, fbnX.getN(),fbnX.getMin(),fbnX.getMax(), 2*fbnY.getN()-1,fbnY.getBinArray().GetArray());
      else if (!fbnX.esp()&& fbnY.esp()) hist = new TH2D(nameHist,titlexy, 2*fbnX.getN()-1,fbnX.getBinArray().GetArray(), fbnY.getN(),fbnY.getMin(),fbnY.getMax());
      else                               hist = new TH2D(nameHist,titlexy, 2*fbnX.getN()-1,fbnX.getBinArray().GetArray(), 2*fbnY.getN()-1,fbnY.getBinArray().GetArray());
      return hist;
    }

    const char *getExpression() const { return Form("%s:%s",fbnY.getExpression(),fbnX.getExpression()); }
    const char *getSelection()  const { return fbnY.getSelection(); }
    const char *expr()          const { return getExpression(); }
    const char *namexy()        const { return Form("%s%s",fbnX.getExpression(),fbnY.getExpression()); }
    const char *namexynn()      const { return Form("%s%d%s%d",fbnX.getExpression(),fbnX.getN(),fbnY.getExpression(),fbnY.getN()); }
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

    int icvs(int iX, int iY) { return (fbnY.getN()-iY-1)*fbnX.getN()+iX+1; }

    TObjArray *rangeBoxGrid();
    TGraph *rangeBox(int binX=-1, int binY=-1);

    TString print(bool pout=1) const {
      TString line1 = fbnX.print(pout);
      TString line2 = fbnY.print(pout);
      return (line1 + ", " + line2);
    }
};

ejungwoo::binning2 ejungwoo::binning::operator*(const binning binn) {
  return binning2(copy(), binn);
}
ejungwoo::binning2 ejungwoo::binning::mult(const binning *binn) {
  return binning2(copy(), binn->copy());
}

TObjArray *ejungwoo::binning2::rangeBoxGrid() {
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

TGraph *ejungwoo::binning2::rangeBox(int binX, int binY) {
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

void ejungwoo::setCanvasPar(parContainer *par)
{
  if (par -> CheckPar("set_stats")) fParCvs.setStats = par -> GetParBool("set_stats");
  if (par -> CheckPar("set_div_idx_yx")) fParCvs.setDivIndexYX = par -> GetParBool("set_div_idx_yx");
  fParCvs.font = par -> GetParInt("font");
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
    fParCvs.axisNdivision[ixyz] = par -> GetParInt("axis_ndivision",ixyz);
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
  if (par -> CheckPar("num_contours")) fParCvs.numContours = par -> GetParInt("num_contours");
}

void ejungwoo::setAttributePar(parContainer *par)
{
  fParAtt.lineStyle = par -> GetParVInt("line_style");
  fParAtt.lineWidth = par -> GetParVInt("line_width");
  fParAtt.lineColor = par -> GetParVInt("line_color");
  fParAtt.markerSize = par -> GetParVDouble("marker_size");
  fParAtt.markerStyle = par -> GetParVInt("marker_style");
  fParAtt.markerColor = par -> GetParVInt("marker_color");
  fParAtt.fillStyle = par -> GetParVInt("fill_style");
  fParAtt.fillColor = par -> GetParVInt("fill_color");
  if (par -> CheckPar("text_font")) fParAtt.textFont = par -> GetParInt("text_font");
  if (par -> CheckPar("text_size")) fParAtt.textSize = par -> GetParInt("text_size");
  if (par -> CheckPar("text_align")) fParAtt.textAlign = par -> GetParInt("text_align");
  if (par -> CheckPar("text_color")) fParAtt.textColor = par -> GetParVInt("text_color");
}

/// Make new pave text (TPaveText)
/// @param content  text
/// @param vpad  canvas
/// @param x1  x1 (left position) in ratio (0.~1.)
/// @param y1  y1 (bottom position) in ratio (0.~1.)
/// @param dx  dx in ratio (0.~1.)
/// @param dy  dy in ratio (0.~1.)
/// @param tf  text font
/// @param ts  text size
/// @param ta  text align
/// @param fc  font color
/// @param fs  font size
/// @param bs  border size
TPaveText* ejungwoo::newpt(TString content, TVirtualPad *vpad, double x1, double y1, double dx, double dy, int tf, int ts, int ta, int fc, int fs, int bs)
{
  TPaveText *paveText = nullptr;
  if (content.Index(";cord")) {
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
  paveText -> SetTextSize(ts);
  paveText -> SetTextAlign(ta);
  paveText -> SetFillColor(fc);
  paveText -> SetFillStyle(fs);
  paveText -> SetBorderSize(bs);
  return paveText;
}

TCanvas* ejungwoo::canvas(const char *nameCvs, int dx, int dy, TArrayD poss)
{
  if (strcmp(nameCvs,"")==0)
    nameCvs = Form("canvas_%d",fCanvasIndex);
  fCanvasIndex++;

  auto cvs = new TCanvas(nameCvs, nameCvs, dx, dy);

  int countPad = 0;
  int numPads = poss.GetSize()/4;
  for (auto iPad=0; iPad<numPads; ++iPad)
  {
    countPad++;
    auto x1 = poss.At(numPads*iPad+0);
    auto y1 = poss.At(numPads*iPad+1);
    auto x2 = poss.At(numPads*iPad+2);
    auto y2 = poss.At(numPads*iPad+3);

    const char *namePad = Form("%s_%d",cvs->GetName(),countPad);
    auto padi = new TPad(namePad,"",x1,y1,x2,y2);
    padi -> SetNumber(countPad);
    padi -> Draw();
  }

  return cvs;
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
  setCanvasPar(par);

  // size
  int innSize[2] = {fParCvs.padSize[0]+fParCvs.marginInn[0]+fParCvs.marginInn[1], fParCvs.padSize[1]+fParCvs.marginInn[2]+fParCvs.marginInn[3]};
  int topSize[2] = {fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[0]*nx+(nx-1)*fParCvs.marginDiv[0], fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1]*ny+(ny-1)*fParCvs.marginDiv[1]};
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

  //const char *titleCvs = Form("%s.%d.%d.%d.   %s",nameConf,nx,ny,0,nameCvs);
  const char *titleCvs = Form("%s  ____  %s.%d.%d.%d.",nameCvs,nameConf,nx,ny,0);
  auto cvs = new TCanvas(nameCvs, titleCvs, fCanvasIndex*10, 25+fCanvasIndex*10, fullSize[0], fullSize[1]);
  cvs -> SetMargin(marginRatioTop[0],marginRatioTop[1],marginRatioTop[2],marginRatioTop[3]);

  if (fCanvasArray==nullptr) fCanvasArray = new TObjArray(50);
  fCanvasArray -> Add(cvs);

  if (nx>1||ny>1||TString(nameConf).Index("nn")>=0)
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

      cvs -> cd();
      countPad++;

      double y2 = 1 - (iy>0)*(padRatioInnEdge1[3] + (iy-1)*padRatioInn[1]) - iy*marginRatioDiv[1];
      double y1 = y2 - padRatioInn[1];
      if (iy==0&&iy==ny-1) y1 = y2 - padRatioInnEdge2[1];
      else if (iy==0)      y1 = y2 - padRatioInnEdge1[3];
      else if (iy==ny-1)   y1 = y2 - padRatioInnEdge1[2];
      if (y1<0) y1 = 0;

      double x1 = (ix>0)*(padRatioInnEdge1[0] + (ix-1)*padRatioInn[0]) + ix*marginRatioDiv[0];
      double x2 = x1 + padRatioInn[0];
      if (ix==0&&ix==nx-1) x2 = x1 + padRatioInnEdge2[0];
      else if (ix==0)      x2 = x1 + padRatioInnEdge1[0];
      else if (ix==nx-1)   x2 = x1 + padRatioInnEdge1[1];
      if (x1<0) x1 = 0;

      const char *namePad = Form("%s_%d",cvs->GetName(),countPad);
      const char *titlePad = Form("____  %s.%d.%d.%d",nameConf,nx,ny,countPad);
      auto padi = new TPad(namePad,titlePad,x1,y1,x2,y2);
      padi -> SetNumber(countPad);
      padi -> SetMargin(marginRatioPad[0],marginRatioPad[1],marginRatioPad[2],marginRatioPad[3]);
      padi -> Draw();

      if (ix==0&&ix==nx-1) {
        padi -> SetLeftMargin (marginRatioEd2[0]);
        padi -> SetRightMargin(marginRatioEd2[1]);
      }
      else if (ix==0)    padi -> SetLeftMargin (marginRatioEdg[0]);
      else if (ix==nx-1) padi -> SetRightMargin(marginRatioEdg[1]);

      if (iy==0&&iy==ny-1) {
        padi -> SetBottomMargin(marginRatioEd2[2]);
        padi -> SetTopMargin   (marginRatioEd2[3]);
      }
      else if (iy==ny-1) padi -> SetBottomMargin(marginRatioEdg[2]);
      else if (iy==0)    padi -> SetTopMargin   (marginRatioEdg[3]);
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

TH1 *ejungwoo::project(TTree *tree, TString formula, TCut cut, TString name, TString title, int nx, int x1, int x2, int ny, int y1, int y2)
{
  TString bnamex;
  TString bnamey;
  if(formula.Index(":")>=0) {
    TObjArray *tokens = formula.Tokenize(":");
    TString bnamex=((TObjString *)tokens->At(1))->GetString();
    TString bnamey=((TObjString *)tokens->At(0))->GetString();
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

  TH1 *histProjected;
  if(ny>0) {
    if(title.IsNull()) title = name+";"+bnamex+";"+bnamey;
    else if (title.Index(";")<0)
      title = title+";"+bnamex+";"+bnamey;
    histProjected = new TH2D(name,title,nx,x1,x2,ny,y1,y2);
  }
  else {
    if(title.IsNull()) title = name+";"+formula;
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

void ejungwoo::savePDF(TString nameVersion) {
  if (!nameVersion.IsNull()) {
    gSystem -> mkdir(nameVersion);
    TString pathToData = nameVersion + "/pdf/";
    gSystem -> mkdir(pathToData);
    auto numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      cvs -> SaveAs(pathToData+cvs->GetName()+".pdf");
    }
  }
  else {
    gSystem -> mkdir(Form("pdf"));
    auto numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      cvs -> SaveAs(Form("pdf/%s.pdf",cvs->GetName()));
    }
  }
}

void ejungwoo::savePNG(TString nameVersion) {
  if (!nameVersion.IsNull()) {
    gSystem -> mkdir(nameVersion);
    TString pathToData = nameVersion + "/figures/";
    gSystem -> mkdir(pathToData);
    auto numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      cvs -> SaveAs(pathToData+cvs->GetName()+".png");
    }
  }
  else {
    gSystem -> mkdir(Form("figures"));
    auto numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      cvs -> SaveAs(Form("figures/%s.png",cvs->GetName()));
    }
  }
}

void ejungwoo::saveRoot(TObject *obj, TString nameFile, TString nameVersion, bool savePrimitives, bool simplifyNames)
{
  TString pathToData = "rooto/";
  if (!nameVersion.IsNull()) {
    gSystem -> mkdir(nameVersion);
    pathToData = nameVersion + "/rooto/";
  }
  gSystem -> mkdir(pathToData);

  TFile *fileOut = nullptr;
  if (!nameFile.IsNull()) {
    if (fVerboseLevel>=verboseLevel::kNormal)
      cout_info << "Creating " << nameFile << endl;
    fileOut = new TFile(nameFile, "recreate");
  }

  auto saveCanvasRoot = [fileOut,savePrimitives,pathToData,simplifyNames](TCanvas *cvs)
  {
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
    auto numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      saveCanvasRoot(cvs);
    }
  }
  else if (obj->InheritsFrom(TCanvas::Class()))
    saveCanvasRoot((TCanvas *) obj);
  else {
    if (nameFile.IsNull())
      nameFile = pathToData+obj->GetName()+".root";
    if (fVerboseLevel>=verboseLevel::kNormal)
      cout_info << "Creating " << nameFile << endl;
    TFile *fileObjOut = fileOut;
    if (fileOut==nullptr) fileObjOut = new TFile(nameFile,"recreate");
    fileObjOut -> cd();
    obj -> Write();
    if (fileOut==nullptr) fileObjOut -> Close();
  }
  if (fileOut!=nullptr) fileOut -> Close();
}

void ejungwoo::saveRoot(TString nameVersion, bool savePrimitives, bool simplifyNames) { saveRoot((TObject*)nullptr, "", nameVersion, savePrimitives, simplifyNames); }

/// save all canvases created so far in [nameVersion]/figures/[cvs-name].png and [nameVersion]/pdf/[cvs-name].pdf files
void ejungwoo::saveAll(TString nameVersion) {
  savePNG(nameVersion);
  savePDF(nameVersion);
  saveRoot(nameVersion);
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
TH1 *ejungwoo::make(TH1 *hist, TVirtualPad *vpad, const char *drawOption)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd();

  int nx = 1;
  int ny = 1;
  int cidx = 0;
  TString nameConf = "";
  if (!findConf(padi,nameConf,nx,ny,cidx)) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  setCanvasPar(par);

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
    axis -> SetTitleFont(fParCvs.font);
    axis -> SetLabelFont(fParCvs.font);
    axis -> SetLabelSize(fParCvs.labelSize[ixyz]);
    axis -> SetLabelOffset(fParCvs.labelOffset[ixyz]);
    axis -> SetTickSize(fParCvs.tickSize[ixyz]);
    axis -> SetNdivisions(fParCvs.axisNdivision[ixyz]);
  }

  hist -> SetContour(fParCvs.numContours);

  if (!(nx==1&&ny==1)&&cidx>0)
  {
    if (ix==0&&ix==nx-1) {}
    else  {
      if (ix!=0 && fParCvs.removeInnerPadAxis[1]) {
        hist -> GetYaxis() -> SetLabelOffset(100);
        hist -> GetYaxis() -> SetTitleOffset(100);
      }
      if (ix==nx-1&&fParCvs.removeInnerZaxis) {
        TString drawOptionStr = drawOption;
        if (drawOptionStr.Index("z")>=0)
          drawOptionStr.ReplaceAll("z","");
        drawOption = drawOptionStr;
      }
    }
    if (iy==0&&iy==ny-1) {}
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
      mainTitle -> SetTextFont(fParCvs.font);
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
TGaxis *ejungwoo::drawz(TH1* hist, TVirtualPad *vpad, const char *titlez)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd();

  int nx = 1;
  int ny = 1;
  int cidx = 0;
  TString nameConf = "";
  if (!findConf(padi,nameConf,nx,ny,cidx)) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  setCanvasPar(par);

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

  auto zaxis = new TGaxis(padi->GetUxmax(),hist->GetYaxis()->GetXmin(),padi->GetUxmax(), hist->GetYaxis()->GetXmax(), hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax(), fParCvs.axisNdivision[2] ,"+L");
  zaxis -> SetTitle(titlez);
  if (fParCvs.titleAlign[2]==2) zaxis -> CenterTitle();
  else zaxis -> CenterTitle(false);
  zaxis -> SetBit(TAxis::kRotateTitle);
  zaxis -> SetTitleOffset(hist->GetZaxis()->GetTitleOffset());
  zaxis -> SetTitleSize(fParCvs.titleSize[2]);
  zaxis -> SetTitleColor(fParCvs.titleColor[2]);
  zaxis -> SetLabelSize(fParCvs.labelSize[2]);
  zaxis -> SetLabelOffset(fParCvs.labelOffset[2]);
  zaxis -> SetTickSize(fParCvs.tickSize[2]);
  //zaxis -> SetNdivisions(fParCvs.axisNdivision[2]);
  zaxis -> SetNdivisions(hist->GetYaxis()->GetNdivisions());
  zaxis -> SetTitleFont(fParCvs.font);
  zaxis -> SetLabelFont(fParCvs.font);

  if (!(nx==1&&ny==1)&&cidx>0) {
    if (ix==0&&ix==nx-1) {}
    else if (ix!=nx-1&&fParCvs.removeInnerZaxis) {
      zaxis -> SetTitleOffset(100);
      zaxis -> SetLabelOffset(100);
    }
  }

  zaxis -> Draw();

  return zaxis;
}

/// Make hist fit to the pad and draw with drawOption.
TH1 *ejungwoo::draw(TH1 *hist, TVirtualPad *vpad, const char *drawOption)
{
  vpad -> cd();
  hist -> Draw(drawOption);
  make(hist,vpad,drawOption);
  return hist;
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
  double x1New = (1.-x1)*x1Frame + x1*(x2Frame);
  double y1New = (1.-y1)*y1Frame + y1*(y2Frame);
  if (x1<0) x1New = x2Frame-dxNew;
  if (y1<0) y1New = y2Frame-dyNew;
  double x2New = x1New + dxNew;
  double y2New = y1New + dyNew;

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
TLegend *ejungwoo::make(TLegend *legend, TVirtualPad *vpad, double x1InRatio, double y1InRatio, double dxInRatio, double dyInRatio, double marginObj, int fillStyle, int fillColor, int borderSize)
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
  setCanvasPar(par);

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

  if (dxInRatio<=0) dxInRatio = (fWidthDefaultLegend + dxTextMax * fWidthPerLengthLegend) * fLegendTextScale;
  if (dxInRatio>1) dxInRatio = 1;
  if (dyInRatio<=0) dyInRatio = (legend -> GetNRows() * fLegendHeightPerEntry) * fLegendTextScale;

  double x1Legend = x1InRatio;
  double y1Legend = y1InRatio;
  double dxLegend = dxInRatio;
  double dyLegend = dyInRatio;
  double xUnit, yUnit;
  findxy(pad, x1Legend, y1Legend, dxLegend, dyLegend, xUnit, yUnit);
  double x2Legend = x1Legend + dxLegend;
  double y2Legend = y1Legend + dyLegend;

  double marginLegend = marginObj;
  if (nullLegend)
    marginLegend = 0;
  else if (marginLegend<0) {
    marginLegend = 1.-((dxTextMax*fWidthPerLengthLegend*fLegendTextScale)/(dxLegend/xUnit));
    if (marginLegend<0.25) marginLegend = 0.25;
  }

  legend -> SetX1(x1Legend + xOffset*xUnit);
  legend -> SetX2(x2Legend + xOffset*xUnit);
  legend -> SetY1(y1Legend + yOffset*yUnit);
  legend -> SetY2(y2Legend + yOffset*yUnit);
  legend -> SetFillStyle(fillStyle);
  legend -> SetFillColor(fillColor);
  legend -> SetBorderSize(borderSize);
  legend -> SetTextFont(fParCvs.font);
  legend -> SetMargin(marginLegend);

  return legend;
}

/// Make legend fit to the pad and draw.
TLegend *ejungwoo::draw(TLegend* legend, TVirtualPad *vpad, double x1InRatio, double y1InRatio, double dxInRatio, double dyInRatio, double marginObj, int fillStyle, int fillColor, int borderSize)
{
  vpad -> cd();
  make(legend,vpad,x1InRatio,y1InRatio,dxInRatio,dyInRatio,marginObj,fillStyle,fillColor,borderSize);
  legend -> Draw();
  return legend;
}

/// Set attribute of (line, marker and filling) to [histogram, graph, function, marker and line] from given nameConf
TObject *ejungwoo::att(TObject *obj, int idx, const char *nameConf)
{
  if (TString(nameConf).IsNull()) nameConf = fNameAttributeConf;
  auto par = conf(nameConf);
  setAttributePar(par);

  if (obj->InheritsFrom(TH1::Class())) {
    auto hist = (TH1 *) obj;
    hist -> SetLineStyle(fParAtt.lineStyle[idx]);
    hist -> SetLineWidth(fParAtt.lineWidth[idx]);
    hist -> SetLineColor(fParAtt.lineColor[idx]);
    hist -> SetMarkerSize(fParAtt.markerSize[idx]);
    hist -> SetMarkerStyle(fParAtt.markerStyle[idx]);
    hist -> SetMarkerColor(fParAtt.markerColor[idx]);
    hist -> SetFillStyle(fParAtt.fillStyle[idx]);
    hist -> SetFillColor(fParAtt.fillColor[idx]);
  }

  else if (obj->InheritsFrom(TGraph::Class())) {
    auto graph = (TGraph *) obj;
    graph -> SetLineStyle(fParAtt.lineStyle[idx]);
    graph -> SetLineWidth(fParAtt.lineWidth[idx]);
    graph -> SetLineColor(fParAtt.lineColor[idx]);
    graph -> SetMarkerSize(fParAtt.markerSize[idx]);
    graph -> SetMarkerStyle(fParAtt.markerStyle[idx]);
    graph -> SetMarkerColor(fParAtt.markerColor[idx]);
    graph -> SetFillStyle(fParAtt.fillStyle[idx]);
    graph -> SetFillColor(fParAtt.fillColor[idx]);
  }

  else if (obj->InheritsFrom(TF1::Class())) {
    auto f1 = (TF1 *) obj;
    f1 -> SetLineStyle(fParAtt.lineStyle[idx]);
    f1 -> SetLineWidth(fParAtt.lineWidth[idx]);
    f1 -> SetLineColor(fParAtt.lineColor[idx]);
    f1 -> SetMarkerSize(fParAtt.markerSize[idx]);
    f1 -> SetMarkerStyle(fParAtt.markerStyle[idx]);
    f1 -> SetMarkerColor(fParAtt.markerColor[idx]);
    f1 -> SetFillStyle(fParAtt.fillStyle[idx]);
    f1 -> SetFillColor(fParAtt.fillColor[idx]);
  }

  else if (obj->InheritsFrom(TMarker::Class())) {
    auto marker = (TMarker *) obj;
    marker -> SetMarkerSize(fParAtt.markerSize[idx]);
    marker -> SetMarkerStyle(fParAtt.markerStyle[idx]);
    marker -> SetMarkerColor(fParAtt.markerColor[idx]);
  }

  else if (obj->InheritsFrom(TLine::Class())) {
    auto line = (TLine *) obj;
    line -> SetLineStyle(fParAtt.lineStyle[idx]);
    line -> SetLineWidth(fParAtt.lineWidth[idx]);
    line -> SetLineColor(fParAtt.lineColor[idx]);
  }

  else if (obj->InheritsFrom(TText::Class())) {
    auto text = (TText *) obj;
    text -> SetTextFont(fParAtt.textFont);
    text -> SetTextSize(fParAtt.textSize);
    text -> SetTextAlign(fParAtt.textAlign);
    if (idx<fParAtt.textColor.size())
      text -> SetTextColor(fParAtt.textColor[idx]);
  }

  return obj;
}

/// Remove trailing zeros and return as a string.
TString ejungwoo::toString(double value)
{
  TString vstring = Form("%f",value);
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
    /*
    db -> AddParticle("deuteron","d"   , 1.87561 ,1,0, 3,"Ion",1000010020);
    db -> AddParticle("triton"  ,"t"   , 2.80892 ,1,0, 3,"Ion",1000010030);
    db -> AddParticle("He3"     ,"he3" , 2.80839 ,1,0, 6,"Ion",1000020030);
    db -> AddParticle("He4"     ,"he4" , 3.72738 ,1,0, 6,"Ion",1000020040);
    db -> AddParticle("Li6"     ,"li6" , 5.6     ,1,0, 9,"Ion",1000030060);
    db -> AddParticle("Li7"     ,"li7" , 6.5     ,1,0, 9,"Ion",1000030070);
    db -> AddParticle("Be7"     ,"be7" , 6.5     ,1,0,12,"Ion",1000040070);
    db -> AddParticle("Be9"     ,"be9" , 8.4     ,1,0,12,"Ion",1000040090);
    db -> AddParticle("Be10"    ,"be10", 9.3     ,1,0,12,"Ion",1000040100);
    db -> AddParticle("Bo10"    ,"bo10", 9.3     ,1,0,15,"Ion",1000050100);
    db -> AddParticle("Bo11"    ,"bo11",10.2     ,1,0,15,"Ion",1000050110);
    db -> AddParticle("C11"     ,"c11" ,10.2     ,1,0,18,"Ion",1000060110);
    db -> AddParticle("C12"     ,"c12" ,11.17793 ,1,0,18,"Ion",1000060120);
    db -> AddParticle("C13"     ,"c13" ,12.11255 ,1,0,18,"Ion",1000060130);
    db -> AddParticle("C14"     ,"c14" ,13.04394 ,1,0,18,"Ion",1000060140);
    db -> AddParticle("N13"     ,"n13" ,12.1     ,1,0,21,"Ion",1000070130);
    db -> AddParticle("N14"     ,"n14" ,13.0     ,1,0,21,"Ion",1000070140);
    db -> AddParticle("N15"     ,"n15" ,14.0     ,1,0,21,"Ion",1000070150);
    db -> AddParticle("O16"     ,"o16" ,14.89917 ,1,0,24,"Ion",1000080160);
    db -> AddParticle("O17"     ,"o17" ,15.83459 ,1,0,24,"Ion",1000080170);
    db -> AddParticle("O18"     ,"o18" ,16.76611 ,1,0,24,"Ion",1000080180);
    */
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
