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
#include "TTimeStamp.h"
#include "TCutG.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

typedef ROOT::TSeqI qi;

namespace ejungwoo
{
  class efile;
  class ecut;
  class epoint;
  class edata;
  class binning; ///< 1-dimensional binning class
  class binning2; ///< 2-dimensional binning class
  typedef KBParameterContainer parContainer; ///< parameter container (for configuration) class.
  class TParCanvas {
    public:
      int canvasType=1;
      int textFont=133, textSize=0, menuSize=25, padSize[2], marginInn[4], marginTop[4], marginDiv[2];
      int axisNdivision[3]={0}, titleAlign[4], titleColor[4], numContours=200, sidePad[2]={0};
      bool removeInnerPadAxis[2], removeInnerMainTitle, removeInnerZaxis, titleRotate[4], setStats=false, setDivIndexYX=false;
      double tickSize[3], labelSize[3], labelOffset[3], titleSize[4], titleOffset[4];
  };
  class TParCanvas2 {
    public:
      bool setStats;
      int numPads, menuSize, numContours;
      int sizeCvsX, sizeCvsY, sizeWindowX, sizeWindowY;
      int marginCvsL, marginCvsR, marginCvsB, marginCvsT;
      double marginRatioCvsL, marginRatioCvsR, marginRatioCvsB, marginRatioCvsT;
      vector<int> posPadX, posPadY;
      vector<int> sizeFrameX,sizeFrameY;
      vector<int> sizePadX, sizePadY;
      vector<int> marginPadL, marginPadR, marginPadB, marginPadT;
      vector<double> marginRatioPadL, marginRatioPadR, marginRatioPadB, marginRatioPadT;
      vector<int> legendFont, legendSize;
      vector<int> titleFont[3], titleSize[3], titleColor[3], titleAlign[3];
      vector<int> labelFont[3], labelSize[3], labelColor[3], labelAlign[3];
      vector<double> titleOffset[3], labelOffset[3];
      vector<double> axisNdivision[3], axisTickSize[3];
      vector<int> mtitleFont, mtitleSize, mtitleAlign, mtitleColor;
      void init() {
        posPadX.clear(); posPadY.clear();
        sizeFrameX.clear(); sizeFrameY.clear();
        sizePadX.clear(); sizePadY.clear();
        marginPadL.clear(); marginPadR.clear(); marginPadB.clear(); marginPadT.clear();
        marginRatioPadL.clear(); marginRatioPadR.clear(); marginRatioPadB.clear(); marginRatioPadT.clear();
        legendFont.clear(); legendSize.clear();
        for (auto ixyz : {0,1,2}) {
          axisNdivision[ixyz].clear(); axisTickSize[ixyz].clear();
          titleFont[ixyz].clear(); titleSize[ixyz].clear(); titleOffset[ixyz].clear(); titleColor[ixyz].clear(); titleAlign[ixyz].clear();
          labelFont[ixyz].clear(); labelSize[ixyz].clear(); labelOffset[ixyz].clear(); labelColor[ixyz].clear(); labelAlign[ixyz].clear();
        }
        mtitleFont.clear(); mtitleSize.clear(); mtitleAlign.clear(); mtitleColor.clear();
        for (auto iPad=0; iPad<numPads; iPad++) {
          posPadX.push_back(0); posPadY.push_back(0);
          sizeFrameX.push_back(0); sizeFrameY.push_back(0);
          sizePadX.push_back(0); sizePadY.push_back(0);
          marginPadL.push_back(0); marginPadR.push_back(0); marginPadB.push_back(0); marginPadT.push_back(0);
          marginRatioPadL.push_back(0); marginRatioPadR.push_back(0); marginRatioPadB.push_back(0); marginRatioPadT.push_back(0);
          legendFont.push_back(0); legendSize.push_back(0);
          for (auto ixyz : {0,1,2}) {
            axisNdivision[ixyz].push_back(0); axisTickSize[ixyz].push_back(0);
            titleFont[ixyz].push_back(0); titleSize[ixyz].push_back(0); titleOffset[ixyz].push_back(0); titleColor[ixyz].push_back(0); titleAlign[ixyz].push_back(0);
            labelFont[ixyz].push_back(0); labelSize[ixyz].push_back(0); labelOffset[ixyz].push_back(0); labelColor[ixyz].push_back(0); labelAlign[ixyz].push_back(0);
          }
          mtitleFont.push_back(0); mtitleSize.push_back(0); mtitleAlign.push_back(0); mtitleColor.push_back(0);
        }
      }
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

  int fHistIndex = 0;
  int fCanvasIndex = 0;
  TObjArray *fCanvasArray = nullptr;
  TClonesArray *fParameterArray = nullptr;
  const char *fNameAttributeConf = "att_color";
  const char *fNameCanvasConf = "single";
  const char *fNameCanvasConf2 = "cvs2";
  TParAttribute fParAtt;
  TParCanvas fParCvs;
  TParCanvas2 fParCvs2;
  void setNextHistIndex(int val) { fHistIndex = val; }
  TString histName() { return Form("hist_%d",fHistIndex++); }
  parContainer *conf(const char *nameConf);
  bool findConf(TVirtualPad *vpad, TString &nameConfCvs, int &nx, int &ny, int &cidx);
  TParCanvas getParCanvas(TString name);
  TParCanvas getParCanvas(parContainer *par);
  TParCanvas2 getParCanvas2(TString name);
  TParCanvas2 getParCanvas2(parContainer *par);
  TParAttribute getParAttribute(TString name);
  TParAttribute getParAttribute(parContainer *par);
  TObject *att(TObject *obj, int idx=0, TString nameConf="");

  void        findxy(TVirtualPad *pad, double &x1, double &y1);
  void        findxy(TVirtualPad *pad, double &x1, double &y1, double &dx, double &dy, double &xUnit, double &yUnit);
  TCanvas*    canvas(const char *nameCvs, int nx, int ny, const char *nameConf="");
  TCanvas*    canvas(const char *nameCvs, int nn, const char *nameConf="");
  TCanvas*    canvas(const char *nameCvs="", const char *nameConf="") { return canvas(nameCvs, 1, 1, nameConf); }
  TCanvas*    canvas2(const char *nameCvs="", const char *nameConf="");
  int         icanvas() { return fCanvasIndex; }
  TPad*       innerpad(TVirtualPad *vpad, double xr1, double yr1, double xr2, double yr2);
  TH1*        make(TH1* hist, TVirtualPad *vpad, TString drawOption="");
  TH1*        make2(TH1* hist, TVirtualPad *vpad, TString drawOption="");
  TH1*        draw(TH1* hist, TVirtualPad *vpad=nullptr, TString drawOption="", int iatt=-1, TString drawConf="");
  TGraph*     draw(TGraph *graph, TVirtualPad *vpad=nullptr, TString drawOption="", int iatt=-1, TString drawConf="");
  TF1*        draw(TF1 *func, TVirtualPad *vpad=nullptr, TString drawOption="", int iatt=-1, TString drawConf="");
  TPad*       redraw(TPad *pad);
  TLegend*    getlg(TVirtualPad *vpad);
  TGaxis*     drawz(TH1* hist, TVirtualPad *vpad, const char *titlez="", bool logaxis=false);
  TGaxis*     drawz2(TH1* hist, TVirtualPad *vpad, const char *titlez="", bool logaxis=false);
  void        drawef(TPad *vpad);
  TLegend*    make(TLegend* lg,TVirtualPad *vp,double x1=-1,double y1=-1,double dx=0,double dy=0,double mgo=-1,int tf=-1,double tz=-1,int ta=-1,int fs=-1,int fc=-1,int bs=-1);
  TLegend*    make2(TLegend* lg,TVirtualPad *vp,double x1=-1,double y1=-1,double dx=0,double dy=0,double mgo=-1,int tf=-1,double tz=-1,int ta=-1,int fs=-1,int fc=-1,int bs=-1);
  TLegend*    draw(TLegend* lg,TVirtualPad *vp,double x1=-1,double y1=-1,double dx=0,double dy=0,double mgo=-1,int tf=-1,double tz=-1,int ta=-1,int fs=-1,int fc=-1,int bs=-1);
  TPad*       newpad(TVirtualPad *vpad, TString name, TString title, double x1, double y1, double dx, double dy);
  TPaveText*  newpt(TString content,TVirtualPad *vpad,double x1=-1,double y1=-1,double dx=0,double dy=0,int tf=133,double tz=20,int ta=12,int fs=0,int fc=0,int bs=0);
  TLatex*     newtt(TString content, double x, double y, int tf=133, double tz=20, int ta=12, int tc=kBlack);
  TH1*        project(TTree *tree, TString formula, TCut cut="", TString name="", TString title="", int nx=0, double x1=0, double x2=0, int ny=0, double y1=0, double y2=0);
  TH1*        project(TTree *tree, TString formula, int nx, double x1=0, double x2=0, int ny=0, double y1=0, double y2=0, TCut cut="", TString name="", TString title="")
              { return project(tree, formula, cut, name, title, nx, x1, x2, ny, y1, y2); }

  TGraphErrors *zbgraph(TH2D *hist, double z_boundary, bool include0=false);
  TGraphErrors *tograph(TH1D *,TString opt=";ex;ey",int iatt=-1,TString nameConf="");
  TGraphErrors *tograph(TF1 *,TString opt=";ex;ey",int iatt=-1,TString nameConf="");
  double get_x(TGraph *graph,int i);
  double get_y(TGraph *graph,int i);
  TString tok(TString line, TString token, int i);
  TString tok(TObjArray *line, int i);
  TF1 *fitg(TH1 *hist, double distSigma=2.5, Option_t *opt="RQ0");
  TF1 *fitg(TGraph *graph, double distSigma=2.5, Option_t *opt="RQ0");
  TGraph *extractPoints(TGraph *graph_hand_pointed, double x1, double x2, double y1, double y2, TString saveName="");

  TCutG *cutg(TString f, TString cut_name, TString x="x", TString y="y"); ///< set TCutG from file name
  TCutG *cutg(TFile  *f, TString cut_name, TString x="x", TString y="y"); ///< set TCutG from file
  TH1 *cutg(TH1 *h, TCutG *cut); ///< recreate histogram satisfying graphical cut

  TString get_value(TString array, TString skey, bool exact=false, TString token=";");
  TString get_value(TObjArray *elements, TString skey, bool exact=false);
  bool exist_key   (TObjArray *elements, TString skey);
  bool exist_keyv  (TObjArray *elements, TString skey);
  int get_valuei   (TObjArray *elements, TString skey);
  double get_valued(TObjArray *elements, TString skey);

  TString readline(std::ifstream& infile); ///< read line and return as TString of stream infile

  void colorwheel();
  void markers();
  TString toString(double value, int ndigit=-1);
  TDatabasePDG *particleDB();
  TParticlePDG *particle(TString name);
  TParticlePDG *particle(int pdg);
  TString makeNameVersion(TString nameVersion="", bool createDirectory=true);
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

class ejungwoo::efile
{
  private:
    std::ifstream fFile;
    std::string fLine;
    istringstream fStream;

  public:
    efile(const char *fileName) { fFile.open(fileName); }

    bool good() { return (!fFile.fail()); }
    bool bad() { return (fFile.fail()); }
    bool next(int iit=1) {
      if (std::getline(fFile, fLine)) {
        fStream = istringstream(fLine);
        if (iit==1)
          return true;
        else
          return next(iit-1);
      }
      return false;
    }
    TString line() { return TString(fLine); }
    TString ws(int iit=1) { TString values; fStream >> values; if (iit==1) return values; return ws(iit-1); }
        int wi(int iit=1) {     int valuei; fStream >> valuei; /*cout << valuei << endl;*/ if (iit==1) return valuei; return wi(iit-1); }
     double wd(int iit=1) {  double valued; fStream >> valued; /*cout << valued << endl;*/ if (iit==1) return valued; return wd(iit-1); }
};

class ejungwoo::ecut
{
  private:
    vector<TString> fName;
    vector<double> fMin;
    vector<double> fMax;

  public:
    ecut(TString name, double min1, double max1, double min2=0, double max2=0, double min3=0, double max3=0, double min4=0, double max4=0) {
      addcut(name,min1,max1);
      if (!(min2==0&&max2==0)) addcut(name,min2,max2);
      if (!(min3==0&&max3==0)) addcut(name,min3,max3);
      if (!(min4==0&&max4==0)) addcut(name,min4,max4);
    }
    ecut(double min1, double max1, double min2=0, double max2=0, double min3=0, double max3=0, double min4=0, double max4=0)
      : ecut("",min1,max1,min2,max2,min3,max3,min4,max4) {}

    TCut getCut() {
      TString scut;
      int n = fName.size();
      for (int i=0 ; i<n; ++i) {
        if (i!=0) scut += "&&";
             if (fMax[i]== DBL_MAX) scut += Form("(%s<%s)",fName[i].Data(),toString(fMin[i]).Data());
        else if (fMin[i]==-DBL_MAX) scut += Form("(%s>%s)",fName[i].Data(),toString(fMax[i]).Data());
        else                        scut += Form("(%s>%s&&%s<%s)",fName[i].Data(),toString(fMax[i]).Data(),fName[i].Data(),toString(fMin[i]).Data());
      }
      return TCut(scut);
    }

    void addcut(TString name, double min, double max) { fName.push_back(name); fMin.push_back(min); fMax.push_back(max); }
    void addll (TString name, double min) { addcut(name,  min, DBL_MAX); }
    void addhl (TString name, double max) { addcut(name, -DBL_MAX, max); }

    bool in(TString name, double value) {
      int n = fName.size();
      for (int i=0 ; i<n; ++i) {
        if (name==fName[i]) {
          if (value>=fMin[i]||value<fMax[i])
            return true;
        }
      }
      return false;
    }

    bool in(double value) {
      int n = fName.size();
      for (int i=0 ; i<n; ++i) {
        if (value>fMin[i]||value<fMax[i])
          return true;
      }
      return false;
    }
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
    bool   next(int ioffset=0);
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
bool ejungwoo::binning::next(int ioffset)  {
  if (fBindex-ioffset>fNbins-1) return false;
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
    nameHist = histName();
  auto titlexy = Form("%s;%s;",title.Data(),fTitle);
  TH1D *hist;
  if (esp()) hist = new TH1D(nameHist,titlexy,fNbins,fMin,fMax);
  else if (csp()) hist = new TH1D(nameHist,titlexy,fNbins,fBinArray.GetArray());
  else if (dsp()) hist = new TH1D(nameHist,titlexy,2*fNbins-1,fBinArray.GetArray());
  return hist;
}

TH1D *ejungwoo::binning::project(TTree *tree, TCut selection, TString name, TString title) {
  if (name.IsNull()) {
    name = histName();
    //name.ReplaceAll("/","DV");
    //name.ReplaceAll("(","L");
    //name.ReplaceAll(")","R");
  }
  auto hist = newHist(name,title);
  auto nn = tree -> Project(name,getExpression(),selection);
  cout << name << "  |  " << getExpression() << "  |  " << selection << endl;
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
    binning2(TH1 *hist, double yscale = 1.) {
      if (hist->InheritsFrom(TH1::Class())) {
        fbnX = binning(hist,1);
        fbnY.setN(100);
        fbnY.setMin(hist->GetMinimum());
        fbnY.setMax(hist->GetMaximum()*yscale);
      }
      else if (hist->InheritsFrom(TH2::Class())) {
        fbnX = binning(hist,1);
        fbnY = binning(hist,2);
      }
    }
    ~binning2() {}

    TH2D *newHist(TString nameHist="", TString titleHist="") {
      if (nameHist.IsNull())
        nameHist = histName();
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
      name.ReplaceAll(":","_vs_");
      name.ReplaceAll("/","_o_");
      name.ReplaceAll("(","_");
      name.ReplaceAll(")","_");
      auto hist = newHist(name,title);
      cout << name << "  |  " << getExpression() << "  |  " << selection << endl;
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

class ejungwoo::epoint
{
  private:
    double fValue = 0;
    double fError = 0;
    double fXYZ[3][3];
    double fErrors[5];

    void init() {
      for (auto i=0; i<3; ++i) { fXYZ[i][0] = 0; fXYZ[i][1] = 0; fXYZ[i][2] = 0; }
      for (auto i=0; i<5; ++i) fErrors[i] = 0;
    }
    void makee() { fError = 0; for (auto i=0; i<5; ++i) fError += fErrors[i]*fErrors[i]; fError = sqrt(fError); }

  public:
    epoint() { init(); }
    epoint(double xx,  double val,  double e1) { init(); setx(xx); setv(val,e1,0); }
    epoint(double xx, double x1, double x2, double vv, double e1, double e2, double e3=0, double e4=0, double e5=0) { init(); setxve(xx,x1,x2,vv,e1,e2,e3,e4,e5); }

    void setx(double xx, double x1=0, double x2=0) { fXYZ[0][0] = xx; fXYZ[0][1] = x1; fXYZ[0][2] = x2; }
    void sety(double yy, double y1=0, double y2=0) { fXYZ[1][0] = yy; fXYZ[1][1] = y1; fXYZ[1][2] = y2; }
    void setz(double zz, double z1=0, double z2=0) { fXYZ[2][0] = zz; fXYZ[2][1] = z1; fXYZ[2][2] = z2; }
    void setv(double vv, double e1=0, double e2=0) { fValue     = vv; fErrors[0] = e1; fErrors[1] = e2; makee(); }
    void sete(int ie, int ee) { fErrors[ie] = ee; makee(); }

    void setxve(double xx, double x1, double x2, double vv, double e1, double e2, double e3=0, double e4=0, double e5=0) {
      setx(xx,x1,x2);
      setv(vv);
      sete(0,e1);
      sete(1,e2);
      sete(2,e3);
      sete(3,e4);
      sete(4,e5);
    }

    double x(int i=0) { return fXYZ[0][i]; }
    double y(int i=0) { return fXYZ[1][i]; }
    double z(int i=0) { return fXYZ[2][i]; }
    double wx() { return (fXYZ[0][2] - fXYZ[0][1]); }
    double wy() { return (fXYZ[1][2] - fXYZ[1][1]); }
    double wz() { return (fXYZ[2][2] - fXYZ[2][1]); }
    double value() { return fValue; }
    double error(int i=-1) { if (i>=0) return fErrors[i]; return fError; }
    double errorv(int i=-1) { if (i>=0) return fValue*0.01*fErrors[i]; return fValue*0.01*fError; }

    TString fill(TString format, TString sformat="%f") {
      format.ReplaceAll("[v]",    Form(sformat,fValue));
      format.ReplaceAll("[e]",    Form(sformat,fError));
      format.ReplaceAll("[x]",    Form(sformat,fXYZ[0][0]));
      format.ReplaceAll("[x1]",   Form(sformat,fXYZ[0][1]));
      format.ReplaceAll("[x2]",   Form(sformat,fXYZ[0][2]));
      format.ReplaceAll("[dx]",   Form(sformat,fXYZ[0][2]-fXYZ[0][1]));
      format.ReplaceAll("[dx/2]", Form(sformat,(fXYZ[0][2]-fXYZ[0][1])/2.));
      format.ReplaceAll("[y]",    Form(sformat,fXYZ[1][0]));
      format.ReplaceAll("[y1]",   Form(sformat,fXYZ[1][1]));
      format.ReplaceAll("[y2]",   Form(sformat,fXYZ[1][2]));
      format.ReplaceAll("[dy]",   Form(sformat,fXYZ[1][2]-fXYZ[1][1]));
      format.ReplaceAll("[dy/2]", Form(sformat,(fXYZ[1][2]-fXYZ[1][1])/2.));
      format.ReplaceAll("[z]",    Form(sformat,fXYZ[2][0]));
      format.ReplaceAll("[z1]",   Form(sformat,fXYZ[2][1]));
      format.ReplaceAll("[z2]",   Form(sformat,fXYZ[2][2]));
      format.ReplaceAll("[dz]",   Form(sformat,fXYZ[2][2]-fXYZ[2][1]));
      format.ReplaceAll("[dz/2]", Form(sformat,(fXYZ[2][2]-fXYZ[2][1])/2.));
      format.ReplaceAll("[e1]",   Form(sformat,fErrors[0]));
      format.ReplaceAll("[e2]",   Form(sformat,fErrors[1]));
      format.ReplaceAll("[e3]",   Form(sformat,fErrors[2]));
      format.ReplaceAll("[e4]",   Form(sformat,fErrors[3]));
      format.ReplaceAll("[e5]",   Form(sformat,fErrors[4]));
      format.ReplaceAll("[ev1]",  Form(sformat,fValue*0.01*fErrors[0]));
      format.ReplaceAll("[ev2]",  Form(sformat,fValue*0.01*fErrors[1]));
      format.ReplaceAll("[ev3]",  Form(sformat,fValue*0.01*fErrors[2]));
      format.ReplaceAll("[ev4]",  Form(sformat,fValue*0.01*fErrors[3]));
      format.ReplaceAll("[ev5]",  Form(sformat,fValue*0.01*fErrors[4]));
      return format;
    }
};

class ejungwoo::edata
{
  private:
    TString fName;
    ejungwoo::binning fBinning;
    vector<epoint> fPoints;
    TH1D *fDataHist = nullptr;
    TGraphErrors *fDataGraph = nullptr;

  public:
    edata(TString name="") { fName = name; }
    edata(TH1D *hist) { initData(hist); }
    void init() { fPoints.clear(); fDataHist = nullptr;}

    void addPoint(epoint point1) { fPoints.push_back(point1); }//fPoints.back().fill(TString()); }
    void setPoint(int ipoint, epoint point1) { fPoints[ipoint] = point1; }
    void setHist(TH1D *hist) { fDataHist = (TH1D *) hist -> Clone(Form("%s_c",hist->GetName())); }

    void initData(TH1D *hist);
    void initError(TH1D *hist);
    edata makeError(TH1D *hist);

    void make(); ///< make fDataHist, fDataGraph
    TH1D *getHist() { return fDataHist; }
    TGraphErrors *getGraph() { return fDataGraph; }
    TGraphErrors *getErrorGraph(int ierror);

    void print(TString format="[x] #pm [dx/2] && [v] && [e1] \% && [e2] \% \\\\", TString sformat = "%f");
    void write(TString nameFile, TString format="[x] #pm [dx/2] && [v] && [e1] \% && [e2] \% \\\\");
    void read(TString nameFile);
};

void ejungwoo::edata::read(TString nameFile)
{
  init();
  TFile *file = new TFile(nameFile,"read");
  TTree *tree = (TTree *) file -> Get("data1");
  double bx, bx1, bx2, bvalue, berror, berror1, berror2, berror3, berror4, berror5;
  tree -> Branch("x",&bx);
  tree -> Branch("x1",&bx1);
  tree -> Branch("x2",&bx2);
  tree -> Branch("value",&bvalue);
  tree -> Branch("error",&berror);
  tree -> Branch("error1",&berror1);
  tree -> Branch("error2",&berror2);
  tree -> Branch("error3",&berror3);
  tree -> Branch("error4",&berror4);
  tree -> Branch("error5",&berror5);
  auto nPoints = tree -> GetEntries();
  for (auto iPoint=0; iPoint<nPoints; ++iPoint) {
    tree -> GetEntry(iPoint);
    epoint point1(bx, bx1, bx2, bvalue, berror1, berror2, berror3, berror4, berror5);
    fPoints.push_back(point1);
  }
  make();
}

void ejungwoo::edata::write(TString nameFile, TString format)
{
  auto nPoints = fPoints.size();
  if (nPoints==0) return;
  if (fDataHist==nullptr) make();

  if (nameFile.Index(".")<0) nameFile = nameFile + ".root";
  if (nameFile.Index("/")<0) {
    auto nameVersion = makeNameVersion();
    nameFile = nameVersion + "/data/" + nameFile;
  }

  if (nameFile.EndsWith(".root")) {
    TFile *fileData = new TFile(nameFile,"recreate");
    TTree *tree = new TTree("data1","ejungwoo::edata");
    double bx, bx1, bx2, bvalue, berror, berror1, berror2, berror3, berror4, berror5;
    tree -> SetBranchAddress("x",&bx);
    tree -> SetBranchAddress("x1",&bx1);
    tree -> SetBranchAddress("x2",&bx2);
    tree -> SetBranchAddress("value",&bvalue);
    tree -> SetBranchAddress("error",&berror);
    tree -> SetBranchAddress("error1",&berror1);
    tree -> SetBranchAddress("error2",&berror2);
    tree -> SetBranchAddress("error3",&berror3);
    tree -> SetBranchAddress("error4",&berror4);
    tree -> SetBranchAddress("error5",&berror5);
    fBinning.reset();
    while (fBinning.next()) {
      bx = fPoints[fBinning.ii()].x();
      bx1 = fPoints[fBinning.ii()].x(1);
      bx2 = fPoints[fBinning.ii()].x(2);
      bvalue = fPoints[fBinning.ii()].value();
      berror = fPoints[fBinning.ii()].error();
      berror1 = fPoints[fBinning.ii()].error(1);
      berror2 = fPoints[fBinning.ii()].error(2);
      berror3 = fPoints[fBinning.ii()].error(3);
      berror4 = fPoints[fBinning.ii()].error(4);
      berror5 = fPoints[fBinning.ii()].error(5);
      tree -> Fill();
    }
    fileData -> cd();
    if (fDataHist!=nullptr) fDataHist -> Write();
    tree -> Write();
  }
  else {
    ofstream fileData(nameFile.Data());
    fileData << format.Data() << endl;
    fBinning.reset();
    while (fBinning.next()) {
      auto line = fPoints[fBinning.ii()].fill(format);
      fileData << line << endl;
    }
  }
}

void ejungwoo::edata::print(TString format, TString sformat) {
  fBinning.reset();
  while (fBinning.next()) {
    auto line = fPoints[fBinning.ii()].fill(format,sformat);
    cout << line << endl;
  }
}

TGraphErrors *ejungwoo::edata::getErrorGraph(int ierror) {
  auto graph = new TGraphErrors();
  fBinning.reset();
  while (fBinning.next()) {
    cout << fBinning.val() << " " << fPoints[fBinning.hi()].value() << " " << fPoints[fBinning.hi()].error() << endl;
    graph -> SetPoint(fBinning.ii(),fBinning.val(),fPoints[fBinning.hi()].error(ierror));
  }
  return graph;
}

void ejungwoo::edata::make()
{
  auto nPoints = fPoints.size();
  bool equalW = true;
  double wx1 = fPoints[0].wx();
  for (auto iPoint=1; iPoint<nPoints; ++iPoint) {
    auto wx2 = fPoints[iPoint].wx();
    if (wx1!=wx2) {
      equalW = false;
      break;
    }
  }
  if (equalW) {
    auto x1 = fPoints[0].x(1);
    auto x2 = fPoints[nPoints-1].x(2);
    if (x1==0&&x2==0) {
      auto dx = (fPoints[nPoints-1].x()-fPoints[0].x())/(nPoints-1)/2.;
      x1 = fPoints[0].x() - dx;
      x2 = fPoints[nPoints-1].x() + dx;
    }
    fBinning = binning(nPoints,x1,x2);
    fDataHist = fBinning.newHist();
    fDataGraph = new TGraphErrors();
    fBinning.reset();
    while (fBinning.next()) {
      fDataHist -> SetBinContent(fBinning.hi(),fPoints[fBinning.hi()].value());
      fDataHist -> SetBinError  (fBinning.hi(),fPoints[fBinning.hi()].error());
      fDataGraph -> SetPoint(fBinning.ii(),fBinning.val(),fPoints[fBinning.hi()].value());
      fDataGraph -> SetPointError(fBinning.ii(),fBinning.width()/2.,fPoints[fBinning.hi()].error());
      cout << "make " << fBinning.val() << " " << fPoints[fBinning.hi()].value() << " " << fPoints[fBinning.hi()].error() << endl;
    }
  }
  else {// TODO
  }
}

void ejungwoo::edata::initData(TH1D *hist) {
  init();
  //fName = hist -> GetName();
  setHist(hist);
  fBinning = binning(hist);
  fBinning.reset();
  while (fBinning.next()) {
    double xx = fBinning.val();
    double x1 = fBinning.low();
    double x2 = fBinning.high();
    double vv = hist -> GetBinContent(fBinning.hi());
    double e1 = 100 * vv * hist -> GetBinError(fBinning.hi());
    epoint point1(xx,x1,x2,vv,e1,0);
    fPoints.push_back(point1);
  }
}

void ejungwoo::edata::initError(TH1D *hist) {
  init();
  setHist(hist);
  fBinning = binning(hist);
  fBinning.reset();
  while (fBinning.next()) {
    double xx = fBinning.val();
    double x1 = fBinning.low();
    double x2 = fBinning.high();
    double vv = 1;
    double e1 = hist -> GetBinContent(fBinning.hi());
    epoint point1(xx,x1,x2,vv,e1,0);
    fPoints.push_back(point1);
  }
}

ejungwoo::edata ejungwoo::edata::makeError(TH1D *hist)
{
  ejungwoo::edata data1;
  fBinning = binning(hist);
  fBinning.reset();
  while (fBinning.next()) {
    double xx = fBinning.val();
    double x1 = fBinning.low();
    double x2 = fBinning.high();
    double vv = hist -> GetBinContent(fBinning.hi());
    double e1 = 100 * vv * hist -> GetBinError(fBinning.hi());
    double e2 = fDataHist -> Interpolate(xx);
    epoint point1(xx,x1,x2,vv,e1,e2);
    double ee = point1.errorv();
    hist -> SetBinError(fBinning.hi(),ee);
    data1.addPoint(point1);
  }
  data1.setHist(hist);
  return data1;
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
    //if (titleCvs.Index("____  ")<0)
    if (titleCvs.Index("____  ",6,0,TString::kExact)<0)
      return false;
    int indexConf = 1;
    //if (titleCvs.Index("____")==0)
    if (titleCvs.Index("____",4,0,TString::kExact)==0)
      indexConf = 0;
    titleCvs.ReplaceAll("  ____  ",";;");
    titleCvs.ReplaceAll("____  ",";;");
    auto tokens1 = titleCvs.Tokenize(";;");
    TString titleConf = ((TObjString *) tokens1->At(indexConf))->GetString();
    auto tokens2 = titleConf.Tokenize(".");
    if (tokens2->GetEntries()==3) {
      nameConfCvs = ((TObjString *) tokens2->At(0))->GetString();
      nx = ((TObjString *) tokens2->At(1))->GetString().Atoi();
      ny = ((TObjString *) tokens2->At(1))->GetString().Atoi();
      cidx = ((TObjString *) tokens2->At(2))->GetString().Atoi();
      return true;
    }
    if (tokens2->GetEntries()==4) {
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
  if (par -> CheckPar("num_pads")) { fParCvs.canvasType = 2; return fParCvs; }
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

ejungwoo::TParCanvas2 ejungwoo::getParCanvas2(parContainer *par)
{
  if (par->CheckPar("set_stats")) fParCvs2.setStats = par -> GetParBool("set_stats"); else fParCvs2.setStats = 0;
  if (par->CheckPar("menu_size")) fParCvs2.menuSize = par -> GetParInt("menu_size"); else fParCvs2.menuSize = 25;
  if (par->CheckPar("num_contours")) fParCvs2.numContours = par -> GetParInt("num_contours"); else fParCvs2.numContours = 100;

  fParCvs2.marginCvsL = par -> GetParInt("cvs_mg",0);
  fParCvs2.marginCvsR = par -> GetParInt("cvs_mg",1);
  fParCvs2.marginCvsB = par -> GetParInt("cvs_mg",2);
  fParCvs2.marginCvsT = par -> GetParInt("cvs_mg",3);
  fParCvs2.numPads = par -> GetParInt("num_pads");
  fParCvs2.init();
  for (auto iPad=0; iPad<fParCvs2.numPads; iPad++)
  {
    int jPad = (par->CheckPar(Form("pad_%d",iPad)))?par->GetParInt(Form("pad_%d",iPad)):iPad;
    int lPad = jPad;
    fParCvs2.posPadX[iPad]    = par -> GetParInt(Form("pad_xy_%d",iPad),0);
    fParCvs2.posPadY[iPad]    = par -> GetParInt(Form("pad_xy_%d",iPad),1);
    fParCvs2.sizeFrameX[iPad] = par -> GetParInt(Form("pad_xy_%d",iPad),2);
    fParCvs2.sizeFrameY[iPad] = par -> GetParInt(Form("pad_xy_%d",iPad),3);
    lPad = (par->CheckPar(Form("pad_mg_%d",iPad)))?iPad:jPad;
    fParCvs2.marginPadL[iPad] = par -> GetParInt(Form("pad_mg_%d",lPad),0);
    fParCvs2.marginPadR[iPad] = par -> GetParInt(Form("pad_mg_%d",lPad),1);
    fParCvs2.marginPadB[iPad] = par -> GetParInt(Form("pad_mg_%d",lPad),2);
    fParCvs2.marginPadT[iPad] = par -> GetParInt(Form("pad_mg_%d",lPad),3);
    lPad = (par->CheckPar(Form("pad_lg_%d",iPad)))?iPad:jPad;
    fParCvs2.legendFont[iPad] = par -> GetParInt(Form("pad_lg_%d",lPad),0);
    fParCvs2.legendSize[iPad] = par -> GetParInt(Form("pad_lg_%d",lPad),1);
         if (fParCvs2.posPadX[iPad]==-1) { if (iPad==0) fParCvs2.posPadX[iPad] = 0; else fParCvs2.posPadX[iPad] = fParCvs2.posPadX[iPad-1]; }
    else if (fParCvs2.posPadX[iPad]< -1) { if (iPad==0) fParCvs2.posPadX[iPad] = 0; else fParCvs2.posPadX[iPad] = fParCvs2.posPadX[iPad-1] + fParCvs2.sizePadX[iPad-1]; }
         if (fParCvs2.posPadY[iPad]==-1) { if (iPad==0) fParCvs2.posPadY[iPad] = 0; else fParCvs2.posPadY[iPad] = fParCvs2.posPadY[iPad-1]; }
    else if (fParCvs2.posPadY[iPad]< -1) { if (iPad==0) fParCvs2.posPadY[iPad] = 0; else fParCvs2.posPadY[iPad] = fParCvs2.posPadY[iPad-1] + fParCvs2.sizePadY[iPad-1]; }
    if (fParCvs2.sizeFrameX[iPad]<=0) { if (iPad==0) fParCvs2.sizeFrameX[iPad] = 250; else fParCvs2.sizeFrameX[iPad] = fParCvs2.sizeFrameX[iPad-1]; }
    if (fParCvs2.sizeFrameY[iPad]<=0) { if (iPad==0) fParCvs2.sizeFrameY[iPad] = 200; else fParCvs2.sizeFrameY[iPad] = fParCvs2.sizeFrameY[iPad-1]; }
    fParCvs2.sizePadX[iPad] = fParCvs2.marginPadL[iPad] + fParCvs2.sizeFrameX[iPad] + fParCvs2.marginPadR[iPad];
    fParCvs2.sizePadY[iPad] = fParCvs2.marginPadB[iPad] + fParCvs2.sizeFrameY[iPad] + fParCvs2.marginPadT[iPad];
    for (auto ixyz : {0,1,2}) { int jxyz = ixyz;
      lPad = (par->CheckPar(Form("pad_tt_%d",iPad)))?iPad:jPad;
      if (par->GetParN(Form("pad_tt_%d",lPad))<=ixyz*5) { if (par->GetParN(Form("pad_tt_%d",lPad))<=(ixyz-1)*5) jxyz = ixyz-2; else jxyz = ixyz-1; } else jxyz = ixyz;
      fParCvs2.titleFont    [ixyz][iPad] = par -> GetParDouble(Form("pad_tt_%d",lPad),jxyz*5+0);
      fParCvs2.titleSize    [ixyz][iPad] = par -> GetParDouble(Form("pad_tt_%d",lPad),jxyz*5+1);
      fParCvs2.titleOffset  [ixyz][iPad] = par -> GetParDouble(Form("pad_tt_%d",lPad),jxyz*5+2);
      fParCvs2.titleColor   [ixyz][iPad] = par -> GetParColor(Form("pad_tt_%d",lPad),jxyz*5+3);
      fParCvs2.titleAlign   [ixyz][iPad] = par -> GetParDouble(Form("pad_tt_%d",lPad),jxyz*5+4);
      lPad = (par->CheckPar(Form("pad_lb_%d",iPad)))?iPad:jPad;
      if (par->GetParN(Form("pad_lb_%d",lPad))<=ixyz*5) { if (par->GetParN(Form("pad_lb_%d",lPad))<=(ixyz-1)*5) jxyz = ixyz-2; else jxyz = ixyz-1; } else jxyz = ixyz;
      fParCvs2.labelFont    [ixyz][iPad] = par -> GetParDouble(Form("pad_lb_%d",lPad),jxyz*5+0);
      fParCvs2.labelSize    [ixyz][iPad] = par -> GetParDouble(Form("pad_lb_%d",lPad),jxyz*5+1);
      fParCvs2.labelOffset  [ixyz][iPad] = par -> GetParDouble(Form("pad_lb_%d",lPad),jxyz*5+2);
      fParCvs2.labelColor   [ixyz][iPad] = par -> GetParColor(Form("pad_lb_%d",lPad),jxyz*5+3);
      fParCvs2.labelAlign   [ixyz][iPad] = par -> GetParDouble(Form("pad_lb_%d",lPad),jxyz*5+4);
      lPad = (par->CheckPar(Form("pad_ax_%d",iPad)))?iPad:jPad;
      if (par->GetParN(Form("pad_ax_%d",lPad))<=ixyz*2) { if (par->GetParN(Form("pad_ax_%d",lPad))<=(ixyz-1)*2) jxyz = ixyz-2; else jxyz = ixyz-1; } else jxyz = ixyz;
      fParCvs2.axisNdivision[ixyz][iPad] = par -> GetParDouble(Form("pad_ax_%d",lPad),jxyz*2+0);
      fParCvs2.axisTickSize [ixyz][iPad] = par -> GetParDouble(Form("pad_ax_%d",lPad),jxyz*2+1);
    }
    lPad = (par->CheckPar(Form("pad_mt_%d",iPad)))?iPad:jPad;
    fParCvs2.mtitleFont [iPad] = par -> GetParInt(Form("pad_mt_%d",lPad),0);
    fParCvs2.mtitleSize [iPad] = par -> GetParInt(Form("pad_mt_%d",lPad),1);
    fParCvs2.mtitleAlign[iPad] = par -> GetParInt(Form("pad_mt_%d",lPad),2);
    fParCvs2.mtitleColor[iPad] = par -> GetParColor(Form("pad_mt_%d",lPad),3);
  }

  fParCvs2.sizeCvsX = 0;
  fParCvs2.sizeCvsY = 0;
  for (auto iPad=0; iPad<fParCvs2.numPads; iPad++) {
    double x2 = fParCvs2.posPadX[iPad] + fParCvs2.sizePadX[iPad];
    double y2 = fParCvs2.posPadY[iPad] + fParCvs2.sizePadY[iPad];
    if (fParCvs2.sizeCvsX<x2) fParCvs2.sizeCvsX = x2;
    if (fParCvs2.sizeCvsY<y2) fParCvs2.sizeCvsY = y2;
    fParCvs2.posPadX[iPad] += fParCvs2.marginCvsL;
    fParCvs2.posPadY[iPad] += fParCvs2.marginCvsB;
    for (auto ixyz : {0,1,2}) {
      double sizePad = (ixyz==0?fParCvs2.sizePadY[iPad]:fParCvs2.sizePadX[iPad]);
      double sizeFrame = (ixyz==0?fParCvs2.sizeFrameY[iPad]:fParCvs2.sizeFrameX[iPad]);
      double rfpCross = (ixyz==0?double(fParCvs2.sizeFrameX[iPad])/fParCvs2.sizePadX[iPad]:double(fParCvs2.sizeFrameY[iPad])/fParCvs2.sizePadY[iPad]);
      double rfpCurrent = (ixyz==0?double(fParCvs2.sizeFrameY[iPad])/fParCvs2.sizePadY[iPad]:double(fParCvs2.sizeFrameX[iPad])/fParCvs2.sizePadX[iPad]);
      if (fParCvs2.labelOffset[ixyz][iPad]>1) fParCvs2.labelOffset[ixyz][iPad] = double(fParCvs2.labelOffset[ixyz][iPad])/sizePad;
      if (fParCvs2.axisTickSize[ixyz][iPad]>1) fParCvs2.axisTickSize[ixyz][iPad] = double(fParCvs2.axisTickSize[ixyz][iPad])/sizeFrame;
      fParCvs2.axisTickSize[ixyz][iPad] = fParCvs2.axisTickSize[ixyz][iPad]/rfpCross*rfpCurrent;
    }
  }
  fParCvs2.sizeCvsX += fParCvs2.marginCvsL + fParCvs2.marginCvsR;
  fParCvs2.sizeCvsY += fParCvs2.marginCvsB + fParCvs2.marginCvsT;
  fParCvs2.sizeWindowX = fParCvs2.sizeCvsX;
  fParCvs2.sizeWindowY = fParCvs2.sizeCvsY + fParCvs2.menuSize;
  fParCvs2.marginRatioCvsL = double(fParCvs2.marginCvsL)/fParCvs2.sizeCvsX;
  fParCvs2.marginRatioCvsR = double(fParCvs2.marginCvsR)/fParCvs2.sizeCvsX;
  fParCvs2.marginRatioCvsB = double(fParCvs2.marginCvsB)/fParCvs2.sizeCvsY;
  fParCvs2.marginRatioCvsT = double(fParCvs2.marginCvsT)/fParCvs2.sizeCvsY;
  for (auto iPad=0; iPad<fParCvs2.numPads; iPad++) {
    fParCvs2.marginRatioPadL[iPad] = double(fParCvs2.marginPadL[iPad])/fParCvs2.sizePadX[iPad];
    fParCvs2.marginRatioPadR[iPad] = double(fParCvs2.marginPadR[iPad])/fParCvs2.sizePadX[iPad];
    fParCvs2.marginRatioPadB[iPad] = double(fParCvs2.marginPadB[iPad])/fParCvs2.sizePadY[iPad];
    fParCvs2.marginRatioPadT[iPad] = double(fParCvs2.marginPadT[iPad])/fParCvs2.sizePadY[iPad];
  }

  return fParCvs2;
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
ejungwoo::TParCanvas2    ejungwoo::getParCanvas2  (TString name) { return getParCanvas2  (conf(name)); }
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
  if (fParCvs.canvasType==2)
    return canvas2(nameCvs,nameConf);

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
      //double y1 = y2 - padRatioInnEdge2[1] - padRatioInn[1];
      //double y1 = y2 - padRatioInnEdge2[1] - padRatioSidePad[1];
      double y1 = y2 - padRatioSidePad[1];
      //cout << "asdlfkajwoei ============= " << fParCvs.sidePad[1] << " " << double(fParCvs.sidePad[1])/topSize[1] << endl;

      auto padi = new TPad(namePad,titlePad,x1,y1,x2,y2);
      //cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
      padi -> SetNumber(countPad);
      padi -> SetLeftMargin (0);
      padi -> SetRightMargin (0);
      padi -> SetBottomMargin(0);
      padi -> SetTopMargin   (0);
      //padi -> SetRightMargin (double(fParCvs.marginTop[1]+fParCvs.marginInn[1])/(fParCvs.marginTop[0]+fParCvs.marginTop[1]+innSize[1]));
      //padi -> SetBottomMargin(double(fParCvs.marginTop[2]+fParCvs.marginInn[2])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+fParCvs.sidePad[1]));
      padi -> SetBottomMargin(double(fParCvs.marginTop[3]+fParCvs.marginInn[3])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+2*innSize[1]));
      padi -> SetTopMargin   (double(fParCvs.marginTop[3]+fParCvs.marginInn[3])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+2*innSize[1]));

      //double marginRatioEdg[4]; for (auto img : {0,1,2,3}) marginRatioEdg[img] = double(fParCvs.marginTop[img]+fParCvs.marginInn[img])/(fParCvs.marginTop[img]+innSize[int(img<2?0:1)]);
      double(fParCvs.marginTop[2]+fParCvs.marginInn[2])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1]);
      double(fParCvs.marginTop[3]+fParCvs.marginInn[3])/(fParCvs.marginTop[2]+fParCvs.marginTop[3]+innSize[1]);

      padi -> Draw();
    }
  }

  return cvs;
}

TCanvas *ejungwoo::canvas2(const char *nameCvs, const char *nameConf)
{
  if (strcmp(nameCvs,"")==0)
    nameCvs = Form("canvas_%d",fCanvasIndex);
  fCanvasIndex++;

  if (TString(nameConf).IsNull()) nameConf = fNameCanvasConf2;
  auto par = conf(nameConf);
  getParCanvas2(par);

  const char *titleCvs = Form("%s  ____  %s.%d.%d",nameCvs,nameConf,fParCvs2.numPads,0);
  auto cvs = new TCanvas(nameCvs, titleCvs, fCanvasIndex*10, 25+fCanvasIndex*10, fParCvs2.sizeWindowX, fParCvs2.sizeWindowY);
  cvs -> SetMargin(fParCvs2.marginRatioCvsL,fParCvs2.marginRatioCvsR,fParCvs2.marginRatioCvsB,fParCvs2.marginRatioCvsT);
  if (fCanvasArray==nullptr) fCanvasArray = new TObjArray(50);
  fCanvasArray -> Add(cvs);

  for (auto iPad=0; iPad<fParCvs2.numPads; iPad++)
  {
    const char *namePad = Form("%s_%d",cvs->GetName(),iPad);
    const char *titlePad = Form("____  %s.%d.%d",nameConf,fParCvs2.numPads,iPad);
    double x1 = double(fParCvs2.posPadX[iPad])/fParCvs2.sizeCvsX;
    double y1 = double(fParCvs2.posPadY[iPad])/fParCvs2.sizeCvsY;
    double dx = double(fParCvs2.sizePadX[iPad])/fParCvs2.sizeCvsX;
    double dy = double(fParCvs2.sizePadY[iPad])/fParCvs2.sizeCvsY;
    auto pad = new TPad(namePad,titlePad,x1,y1,x1+dx,y1+dy);
    //cout << setw(5) << x1 << setw(5) << y1 << setw(5) << dx << setw(5) << dy << setw(5) << iPad+1 << endl;;
    //auto pad = new TPad(namePad,titlePad,fParCvs2.posPadX[iPad],fParCvs2.posPadY[iPad],fParCvs2.sizePadX[iPad],fParCvs2.sizePadY[iPad]);
    pad -> SetNumber(iPad+1);
    pad -> SetMargin(fParCvs2.marginRatioPadL[iPad],fParCvs2.marginRatioPadR[iPad],fParCvs2.marginRatioPadB[iPad],fParCvs2.marginRatioPadT[iPad]);
    pad -> Draw();
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
    name = histName();
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

TGraphErrors *ejungwoo::zbgraph(TH2D *hist, double z_boundary, bool include0) {
  binning bnx(hist,1);
  binning bny(hist,2);
  auto graph = new TGraphErrors();

  bnx.reset();
  while (bnx.next(-1)) {
    bny.reset();
    while (bny.next()) {
      auto content1 = hist -> GetBinContent(bnx.bi()  ,bny.bi());
      auto content2 = hist -> GetBinContent(bnx.bi()+1,bny.bi());
      if (!include0&&(content1==0||content2==0)) continue;
      if ((content1<=z_boundary&&content2> z_boundary) ||
          (content1> z_boundary&&content2<=z_boundary)) {
        graph -> SetPoint(graph->GetN(),bnx.high(),bny.val());
        graph -> SetPointError(graph->GetN()-1,bnx.width()/2.,bny.width()/2.);
      }
    }
  }
  bny.reset();
  while (bny.next(-1)) {
    bnx.reset();
    while (bnx.next()) {
      auto content1 = hist -> GetBinContent(bnx.bi(),bny.bi()  );
      auto content2 = hist -> GetBinContent(bnx.bi(),bny.bi()+1);
      if (!include0&&(content1==0||content2==0)) continue;
      if ((content1<z_boundary&&content2>=z_boundary) ||
          (content1>=z_boundary&&content2<z_boundary)) {
        graph -> SetPoint(graph->GetN(),bnx.val(),bny.high());
        graph -> SetPointError(graph->GetN()-1,bnx.width()/2.,bny.width()/2.);
      }
    }
  }
  graph -> SetMarkerStyle(20);
  graph -> SetMarkerSize(0.2);
  return graph;
}

TGraphErrors *ejungwoo::tograph(TH1D *hist, TString opt, int iatt, TString nameConf)
{
  binning bnx(hist);
  auto graph = new TGraphErrors();
  bool tostaterr = false;
  bool ignore0 = false;
  bool reversex = false;
  bool reversey = false;
  bool absy = false;
  bool copyex = true;
  bool copyey = true;
  bool limity = false;
  double ymxlow = -DBL_MAX;
  double ypxlow = -DBL_MAX;
  double ymxhigh = DBL_MAX;
  double ypxhigh = DBL_MAX;
  double xlow  = -DBL_MAX;
  double xhigh =  DBL_MAX;
  double multy = 1;

  vector<double> addatxarray;
  vector<double> addatxiarrayx;
  vector<int>    addatxiarrayi;
  const char *evalopt = "";

  auto optArray = opt.Tokenize(";");
  auto numArray = optArray->GetEntries();
  for (auto i=0; i<numArray; ++i) {
    TString option = tok(optArray,i);
    if (option=="0") ignore0 = true;
    if (option=="-x") reversex = true;
    if (option=="-y") reversey = true;
    if (option=="|y|") absy = true;
    if (option=="ex") copyex = true;
    if (option=="ey") copyey = true;
    if (option=="!ex") copyex = false;
    if (option=="!ey") copyey = false;
    if (option=="staterr") tostaterr = true;
    if (option.Index("evals")==0) evalopt = "S";
    if (option.Index("x<")==0) { option.ReplaceAll("x<",""); xhigh=option.Atof(); }
    if (option.Index("multy=")==0) { option.ReplaceAll("multy=",""); multy=option.Atof(); }
    if (option.Index("x>")==0) { option.ReplaceAll("x>",""); xlow =option.Atof(); }
    if (option.Index("y+x<")==0) { option.ReplaceAll("y+x<",""); ypxhigh=option.Atof(); limity = true; }
    if (option.Index("y+x>")==0) { option.ReplaceAll("y+x>",""); ypxlow =option.Atof(); limity = true; }
    if (option.Index("y-x<")==0) { option.ReplaceAll("y-x<",""); ymxhigh=option.Atof(); limity = true; }
    if (option.Index("y-x>")==0) { option.ReplaceAll("y-x>",""); ymxlow =option.Atof(); limity = true; }
    if (option.Index("+x=")==0) { option.ReplaceAll("+x=",""); addatxarray.push_back(option.Atof()); }
    if (option.Index("+xi=")==0) {
      option.ReplaceAll("+xi=","");
      auto opt12 = option.Tokenize(",");
      opt12->GetEntries();
      auto addat = tok(opt12,0).Atof();
      auto idata = tok(opt12,1).Atoi();
      addatxiarrayx.push_back(addat);
      addatxiarrayi.push_back(idata);
    }
  }

  while (bnx.next()) {
    if (bnx.val()<xlow ) continue;
    if (bnx.val()>xhigh) continue;
    if (bnx.val()<xlow ) continue;
    if (bnx.val()>xhigh) continue;

    auto content = hist -> GetBinContent(bnx.bi());
    auto error = hist -> GetBinError(bnx.bi());
    if (ignore0 && content==0) {}
    else {
      double xvalue = bnx.val();
      double yvalue = content;
      if (reversex) xvalue = -xvalue;
      if (reversey) yvalue = -yvalue;
      if (absy&&yvalue<0) yvalue = -yvalue;
      if (tostaterr) {
        if (content==0) yvalue = 0;
        else yvalue = error / content;
      }
      yvalue = yvalue * multy;
      graph -> SetPoint(graph->GetN(), xvalue, yvalue);
      double xerror = bnx.width()/2.;
      double yerror = hist -> GetBinError(bnx.bi());
      yerror = yerror *multy;
      if (!copyex) xerror = 0;
      if (!copyey) yerror = 0;
      graph -> SetPointError(graph->GetN()-1, xerror, yerror);
    }
  }

  auto npoints = graph -> GetN();
  if (npoints>=2)
  {
    auto numxi = addatxiarrayx.size();
    for (auto i=0; i<numxi; ++i) {
      int    idata = addatxiarrayi[i];
      double addatx = addatxiarrayx[i];
      double yatx = get_y(graph,idata);
      double xerror = graph -> GetErrorX(idata);
      double yerror = graph -> GetErrorY(idata);
      graph -> SetPoint(graph->GetN(),addatx,yatx);
      graph -> SetPointError(graph->GetN()-1, xerror, yerror);
    }
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

      auto yatx = graph -> Eval(addatx,(TSpline*)nullptr,evalopt);
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
    graph -> Sort();
  }

  /*
  if (limity) {
    auto npoints = graph -> GetN();
    for (auto ii=0; ii<npoints; ++ii) {
      auto x0 = get_x(graph,ii);
      if (abs(x0-addatx)<dxmin) {
        ibin2 = ibin1;
        ibin1 = ii;
        dxmin = abs(x0-addatx);
      }
    }
    auto yatx = graph -> Eval(addatx,(TSpline*)nullptr,evalopt);
  }
  */


  if (iatt>=0) {
    nameConf = nameConf + ";" + opt;
    ejungwoo::att(graph,iatt,nameConf);
  }

  return graph;
}

TGraphErrors *ejungwoo::tograph(TF1 *f1, TString opt, int iatt, TString nameConf)
{
  double rmin, rmax;
  f1 -> GetRange(rmin,rmax);
  binning bnx(f1->GetNpx(),rmin,rmax);
  auto graph = new TGraphErrors();

  while (bnx.next()) {
    auto content = f1 -> Eval(bnx.val());
    double xvalue = bnx.val();
    double yvalue = content;
    graph -> SetPoint(graph->GetN(), xvalue, yvalue);
  }

  if (iatt>=0) {
    nameConf = nameConf + ";" + opt;
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

TString ejungwoo::makeNameVersion(TString nameVersion, bool createDirectory) {
  if (nameVersion.IsNull()) {
    TTimeStamp time;
    UInt_t year, month, day;
    time.GetDate(true, 0, &year, &month, &day);
    nameVersion = Form("save_%04d_%02d_%02d", year, month, day);
  }
  if (createDirectory)
    gSystem -> mkdir(nameVersion);
  return nameVersion;
}

void ejungwoo::savePDF(TCanvas *cvs, TString nameVersion) {
  cvs -> Modified();
  cvs -> Update();
  nameVersion = makeNameVersion(nameVersion);
  TString pathToData = nameVersion + "/pdf/";
  gSystem -> mkdir(pathToData);
  cvs -> SaveAs(pathToData+cvs->GetName()+".pdf");
}

void ejungwoo::saveEPS(TCanvas *cvs, TString nameVersion) {
  cvs -> Modified();
  cvs -> Update();
  nameVersion = makeNameVersion(nameVersion);
  TString pathToData = nameVersion + "/eps/";
  gSystem -> mkdir(pathToData);
  cvs -> SaveAs(pathToData+cvs->GetName()+".eps");
}

void ejungwoo::savePNG(TCanvas *cvs, TString nameVersion) {
  cvs -> Modified();
  cvs -> Update();
  nameVersion = makeNameVersion(nameVersion);
  TString pathToData = nameVersion + "/png/";
  gSystem -> mkdir(pathToData);
  cvs -> SaveAs(pathToData+cvs->GetName()+".png");
}

void ejungwoo::savePDF(TString nameVersion) {
  nameVersion = makeNameVersion(nameVersion);
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
  nameVersion = makeNameVersion(nameVersion);
  TString pathToData = nameVersion + "/png/";
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
  nameVersion = makeNameVersion(nameVersion);
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
  nameVersion = makeNameVersion(nameVersion);
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
  if (fParCvs.canvasType==2)
    return make2(hist,vpad,drawOption);

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

TH1 *ejungwoo::make2(TH1 *hist, TVirtualPad *vpad, TString drawOption)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd();

  int nx = 1;
  int ny = 1;
  int iPad = 0;
  TString nameConf = "";
  if (!findConf(padi,nameConf,nx,ny,iPad)) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  getParCanvas2(par);

  hist -> SetStats(fParCvs2.setStats);

  for (auto ixyz : {0,1,2}) {
    TAxis *axis;
         if (ixyz==0) axis = (TAxis *) hist -> GetXaxis();
    else if (ixyz==1) axis = (TAxis *) hist -> GetYaxis();
    else if (ixyz==2) axis = (TAxis *) hist -> GetZaxis();
    if (fParCvs2.titleAlign[ixyz][iPad]==2) axis -> CenterTitle();
    else axis -> CenterTitle(false);
    axis -> SetTitleOffset(fParCvs2.titleOffset[ixyz][iPad]);
    axis -> SetTitleSize(fParCvs2.titleSize[ixyz][iPad]);
    axis -> SetTitleColor(fParCvs2.titleColor[ixyz][iPad]);
    axis -> SetTitleFont(fParCvs2.titleFont[ixyz][iPad]);
    axis -> SetLabelFont(fParCvs2.labelFont[ixyz][iPad]);
    axis -> SetLabelSize(fParCvs2.labelSize[ixyz][iPad]);
    axis -> SetLabelOffset(fParCvs2.labelOffset[ixyz][iPad]);
    axis -> SetTickSize(fParCvs2.axisTickSize[ixyz][iPad]);
    int ndivision = fParCvs2.axisNdivision[ixyz][iPad];
    if (ndivision!=0)
      axis -> SetNdivisions(ndivision);
  }

  hist -> SetContour(fParCvs2.numContours);

  /*
  if (!(nx==1&&ny==1)&&iPad>0)
  {
    if (isLeft&&isRight) {}
    else  {
      if (ix!=0 && fParCvs2.removeInnerPadAxis[1]) {
        hist -> GetYaxis() -> SetLabelOffset(100);
        hist -> GetYaxis() -> SetTitleOffset(100);
      }
      if (isRight&&fParCvs2.removeInnerZaxis) {
        TString drawOptionStr = drawOption;
        if (drawOptionStr.Index("z")>=0)
          drawOptionStr.ReplaceAll("z","");
        drawOption = drawOptionStr;
      }
    }
    if (isTop&&isBottom) {}
    else {
      if (iy!=ny-1&&fParCvs2.removeInnerPadAxis[0]) {
        hist -> GetXaxis() -> SetLabelOffset(100);
        hist -> GetXaxis() -> SetTitleOffset(100);
      }
      if (iy!=0&&fParCvs2.removeInnerMainTitle) {
        hist -> SetTitle("");
        hist -> SetTitle("");
      }
    }
  }
  */

  double marginRatioCur[4] = {padi->GetLeftMargin(), padi->GetRightMargin(), padi->GetBottomMargin(), padi->GetTopMargin()};
  padi -> Modified();
  padi -> Update();
  if (padi!=nullptr) {
    auto list_primitive = padi -> GetListOfPrimitives();
    auto mainTitle = (TPaveText*) list_primitive -> FindObject("title");
    if (mainTitle!=nullptr)
    {
      mainTitle -> SetTextFont(fParCvs2.mtitleFont[iPad]);
      mainTitle -> SetTextAlign(fParCvs2.mtitleAlign[iPad]);
      mainTitle -> SetTextSizePixels(fParCvs2.mtitleSize[iPad]);
      mainTitle -> SetTextColor(fParCvs2.mtitleColor[iPad]);
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
  if (fParCvs.canvasType==2)
    return drawz2(hist,vpad,titlez,logaxis);

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

TGaxis *ejungwoo::drawz2(TH1* hist, TVirtualPad *vpad, const char *titlez, bool logaxis)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd();

  int nx = 1;
  int ny = 1;
  int cidx = 0;
  TString nameConf = "";
  if (!findConf(padi,nameConf,nx,ny,cidx)) nameConf = fNameCanvasConf2;
  auto par = conf(nameConf);
  getParCanvas2(par);

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

TH1 *ejungwoo::draw(TH1 *obj, TVirtualPad *vpad, TString drawOption, int iatt, TString drawConf)
{
  if (iatt>=0) att(obj,iatt,drawConf);
  if (vpad==nullptr)
    vpad = canvas();
  vpad -> cd();
  if (vpad->GetPrimitive("TFrame")!=nullptr) {
    if (drawOption.Index("same")<0) drawOption += "same";
  }
  obj -> Draw(drawOption);
  make(obj,vpad,drawOption);
  return obj;
}

TGraph *ejungwoo::draw(TGraph *obj, TVirtualPad *vpad, TString drawOption, int iatt, TString drawConf)
{
  if (iatt>=0) att(obj,iatt,drawConf);
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
  return obj;
}

TF1 *ejungwoo::draw(TF1 *obj, TVirtualPad *vpad, TString drawOption, int iatt, TString drawConf)
{
  if (iatt>=0) att(obj,iatt,drawConf);
  if (vpad==nullptr)
    vpad = canvas();
  vpad -> cd();
  if (vpad->GetPrimitive("TFrame")!=nullptr) {
    if (drawOption.Index("same")<0) drawOption += "same";
  } else {
    drawOption.ReplaceAll("same","");
    if (drawOption.Index("l")<0) drawOption += "l";
  }
  obj -> Draw(drawOption);
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
TLegend *ejungwoo::make(TLegend *legend, TVirtualPad *vpad, double x1InRatio, double y1InRatio, double dxInRatio, double dyInRatio, double marginObj, int tf, double tz, int ta, int fs, int fc, int bs)
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
  if (fParCvs.canvasType==2)
    return make2(legend,vpad,x1InRatio,y1InRatio,dxInRatio,dyInRatio,marginObj,tf,tz,ta,fs,fc,bs);

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

TLegend *ejungwoo::make2(TLegend *legend, TVirtualPad *vpad, double x1InRatio, double y1InRatio, double dxInRatio, double dyInRatio, double marginObj, int tf, double tz, int ta, int fs, int fc, int bs)
{
  auto pad = (TPad *) vpad;

  double fLegendTextScale = 0.028;
  double fWidthPerLengthLegend = .8;
  double fWidthDefaultLegend = fWidthPerLengthLegend*2.5;
  double fLegendHeightPerEntry = 2.5;
  double xOffset = 0;
  double yOffset = 0;

  int dummy, iPad;
  TString nameConf = "";
  if (!findConf(pad,nameConf,dummy,dummy,iPad)) nameConf = fNameCanvasConf;
  auto par = conf(nameConf);
  getParCanvas2(par);

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
  else legend -> SetTextFont(fParCvs2.legendFont[iPad]);

  if (tz>0) legend -> SetTextSize(tz);
  else if (fParCvs2.legendSize[iPad]>0)
    legend -> SetTextSize(fParCvs2.legendSize[iPad]);

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

void ejungwoo::drawef(TPad *vpad) {
  vpad -> cd();
  vpad -> SetFrameLineColor(0);
  //auto histef = vpad -> DrawFrame(0.,100.,0.,100.);
  auto histef = binning2(0,100,0,100).newHist();
  histef -> SetStats(0);
  histef -> GetXaxis() -> SetLabelOffset(100);
  histef -> GetYaxis() -> SetLabelOffset(100);
  histef -> GetXaxis() -> SetNdivisions(0);
  histef -> GetYaxis() -> SetNdivisions(0);
  vpad -> cd();
  histef -> Draw("a");
  ////histef -> Draw();
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
TObject *ejungwoo::att(TObject *obj, int idx, TString nameConf)
{
  TString sNameConf = nameConf;
  int olz=idx, ols=idx, olc=idx, omz=idx, oms=idx, omc=idx, ofs=idx, ofc=idx, otz=idx, otf=idx, otc=idx;
  double xmz=1, xtz=1;

  int vlz = -1, vls = -1, vlc = -1, vms = -1, vmc = -1, vfs = -1, vfc = -1, vta = -1, vtf = -1, vtc = -1;
  double vmz = -1, vtz = -1; 

  if (sNameConf.IsNull())
    sNameConf = fNameAttributeConf;
  else if (sNameConf.Index(";")==0)
    sNameConf = TString(fNameAttributeConf) + ";" + sNameConf;

  if (sNameConf.Index(";")>0)
  {
    auto elements = sNameConf.Tokenize(";");
    sNameConf = ((TObjString *) elements->At(0))->GetString();

    if(exist_keyv(elements,"vlz")) vlz = get_valuei(elements,"vlz");
    if(exist_keyv(elements,"vls")) vls = get_valuei(elements,"vls");
    if(exist_keyv(elements,"vlc")) vlc = get_valuei(elements,"vlc");
    if(exist_keyv(elements,"vmz")) vmz = get_valued(elements,"vmz");
    if(exist_keyv(elements,"vms")) vms = get_valuei(elements,"vms");
    if(exist_keyv(elements,"vmc")) vmc = get_valuei(elements,"vmc");
    if(exist_keyv(elements,"vfs")) vfs = get_valuei(elements,"vfs");
    if(exist_keyv(elements,"vfc")) vfc = get_valuei(elements,"vfc");
    if(exist_keyv(elements,"vtz")) vtz = get_valued(elements,"vtz");
    if(exist_keyv(elements,"vta")) vta = get_valuei(elements,"vta");
    if(exist_keyv(elements,"vts")) vtf = get_valuei(elements,"vts");
    if(exist_keyv(elements,"vtf")) vtf = get_valuei(elements,"vtf");
    if(exist_keyv(elements,"vtc")) vtc = get_valuei(elements,"vtc");

    if (!exist_key(elements,"ll")) {
      if(exist_key(elements,"lz")) olz=idx; else if(exist_keyv(elements,"lz")) olz=get_valuei(elements,"lz"); else olz=-1; 
      if(exist_key(elements,"ls")) ols=idx; else if(exist_keyv(elements,"ls")) ols=get_valuei(elements,"ls"); else ols=-1; 
      if(exist_key(elements,"lc")) olc=idx; else if(exist_keyv(elements,"lc")) olc=get_valuei(elements,"lc"); else olc=-1; 
    }
    if (!exist_key(elements,"mm")) {
      if(exist_key(elements,"mz")) omz=idx; else if(exist_keyv(elements,"mz")) omz=get_valuei(elements,"mz"); else omz=-1; 
      if(exist_key(elements,"ms")) oms=idx; else if(exist_keyv(elements,"ms")) oms=get_valuei(elements,"ms"); else oms=-1; 
      if(exist_key(elements,"mc")) omc=idx; else if(exist_keyv(elements,"mc")) omc=get_valuei(elements,"mc"); else omc=-1; 
    }
    if (!exist_key(elements,"ff")) {
      if(exist_key(elements,"fs")) ofs=idx; else if(exist_keyv(elements,"fs")) ofs=get_valuei(elements,"fs"); else ofs=-1; 
      if(exist_key(elements,"fc")) ofc=idx; else if(exist_keyv(elements,"fc")) ofc=get_valuei(elements,"fc"); else ofc=-1; 
    }
    if (!exist_key(elements,"tm")) {
      if(exist_key(elements,"tz")) otz=idx; else if(exist_keyv(elements,"tz")) otz=get_valuei(elements,"tz"); else otz=-1; 
      if(exist_key(elements,"ts")) otf=idx; else if(exist_keyv(elements,"ts")) otf=get_valuei(elements,"ts"); else otf=-1; 
      if(exist_key(elements,"tf")) otf=idx; else if(exist_keyv(elements,"tf")) otf=get_valuei(elements,"tf"); else otf=-1; 
      if(exist_key(elements,"tc")) otc=idx; else if(exist_keyv(elements,"tc")) otc=get_valuei(elements,"tc"); else otc=-1; 
    }
    if (exist_keyv(elements,"xmz")) xmz = get_valued(elements,"xmz");
    if (exist_keyv(elements,"xtz")) xtz = get_valued(elements,"xtz");
    if (exist_key(elements,"cc")) { olc=idx; omc=idx; ofc=idx; otc=idx; }

    if (olz==-1&&ols==-1&&olc==-1&&omz==-1&&oms==-1&&omc==-1&&ofs==-1&&ofc==-1&&otz==-1&&otf==-1&&otc==-1) {
      olz=idx; ols=idx; olc=idx; omz=idx; oms=idx; omc=idx; ofs=idx; ofc=idx; otz=idx; otf=idx; otc=idx;
    }

    auto par = conf(sNameConf);
    getParAttribute(par);
  }
  else {
    auto par = conf(sNameConf);
    getParAttribute(par);
  }

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

    if (vlz>=0) hist -> SetLineWidth(  vlz);
    if (vls>=0) hist -> SetLineStyle(  vls);
    if (vlc>=0) hist -> SetLineColor(  vlc);
    if (vmz>=0) hist -> SetMarkerSize( vmz);
    if (vms>=0) hist -> SetMarkerStyle(vms);
    if (vmc>=0) hist -> SetMarkerColor(vmc);
    if (vfs>=0) hist -> SetFillStyle(  vfs);
    if (vfc>=0) hist -> SetFillColor(  vfc);
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

    if (vlz>=0) graph -> SetLineWidth(  vlz);
    if (vls>=0) graph -> SetLineStyle(  vls);
    if (vlc>=0) graph -> SetLineColor(  vlc);
    if (vmz>=0) graph -> SetMarkerSize( vmz);
    if (vms>=0) graph -> SetMarkerStyle(vms);
    if (vmc>=0) graph -> SetMarkerColor(vmc);
    if (vfs>=0) graph -> SetFillStyle(  vfs);
    if (vfc>=0) graph -> SetFillColor(  vfc);
  }

  else if (obj->InheritsFrom(TF1::Class())) {
    auto f1 = (TF1 *) obj;
    if (olz>=0) f1 -> SetLineWidth(fParAtt.lineWidth[olz]);
    if (ols>=0) f1 -> SetLineStyle(fParAtt.lineStyle[ols]);
    if (olc>=0) f1 -> SetLineColor(fParAtt.lineColor[olc]);

    if (vlz>=0) f1 -> SetLineWidth(  vlz);
    if (vls>=0) f1 -> SetLineStyle(  vls);
    if (vlc>=0) f1 -> SetLineColor(  vlc);
  }

  else if (obj->InheritsFrom(TMarker::Class())) {
    auto marker = (TMarker *) obj;
    if (omz>=0) marker -> SetMarkerSize(fParAtt.markerSize[omz]*xmz);
    if (oms>=0) marker -> SetMarkerStyle(fParAtt.markerStyle[oms]);
    if (omc>=0) marker -> SetMarkerColor(fParAtt.markerColor[omc]);

    if (vmz>=0) marker -> SetMarkerSize( vmz);
    if (vms>=0) marker -> SetMarkerStyle(vms);
    if (vmc>=0) marker -> SetMarkerColor(vmc);
  }

  else if (obj->InheritsFrom(TLine::Class())) {
    auto line = (TLine *) obj;
    if (ols>=0) line -> SetLineWidth(fParAtt.lineWidth[ols]);
    if (olz>=0) line -> SetLineStyle(fParAtt.lineStyle[olz]);
    if (olc>=0) line -> SetLineColor(fParAtt.lineColor[olc]);

    if (vlz>=0) line -> SetLineWidth( vlz);
    if (vls>=0) line -> SetLineStyle(vls);
    if (vlc>=0) line -> SetLineColor(vlc);
  }

  else if (obj->InheritsFrom(TText::Class())) {
    auto text = (TText *) obj;
    text -> SetTextAlign(fParAtt.textAlign);
    if (otz>=0) text -> SetTextSize(fParAtt.textSize*xtz);
    if (otf>=0) text -> SetTextFont(fParAtt.textFont);
    if (otc>=0) if (otc<fParAtt.textColor.size()) text -> SetTextColor(fParAtt.textColor[otc]);

    if (vta>=0) text -> SetTextAlign(vta);
    if (vtz>=0) text -> SetTextSize( vtz);
    if (vtf>=0) text -> SetTextFont( vtf);
    if (vtc>=0) text -> SetTextColor(vtc);
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

TString ejungwoo::readline(std::ifstream& infile) {
  std::string line;
  if (std::getline(infile, line))
    return TString(line);
  else
    return TString("EOF");
}

//TString ejungwoo::readline(std::ifstream& infile) {
//  while (getline(file_sys,line_dm))
//    count_line++;
//  istringstream ss(line_st);
//  TString line0;
//  ss >> line0;
//  TString line1 = line0(1,line0.Index(")*",2)-1);
//  return ss;
//}

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

TF1 *ejungwoo::fitg(TGraph *graph, double distSigma=2.5, Option_t *opt="RQ0")
{
  int nn = graph -> GetN();
  double x1 = get_x(graph,0) - graph -> GetErrorX(0);
  double x2 = get_x(graph,nn-1) - graph -> GetErrorX(nn-1);
  auto hist = binning(nn,x1,x2).newHist();
  for (auto i=0; i<nn; ++i) {
    double w = get_y(graph,i);
    if (w!=0) hist -> Fill(get_x(graph,i),w);
  }

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

  auto fit = new TF1(Form("%s_fitg",graph -> GetName()),"gaus(0)",xmax-xerr*distSigma,xmax+xerr*distSigma);
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

  graph -> Fit(fit,opt);
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

    graph -> Fit(fit,opt);
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

TCutG *ejungwoo::cutg(TString f, TString cut_name, TString x, TString y) {
  auto file = new TFile(f.Data(),"read");
  return cutg(file, cut_name, x, y);
}

TCutG *ejungwoo::cutg(TFile *file, TString cut_name, TString x, TString y) {
  auto cg = (TCutG *) file -> Get(cut_name.Data());
  cg -> SetVarX(x.Data());
  cg -> SetVarY(y.Data());
  return cg;
}

TH1 *ejungwoo::cutg(TH1 *h, TCutG *cut)
{
  auto h2 = (TH2 *) h;
  auto hnew = (TH2 *) h2->Clone();

  hnew->Reset();
  auto nx=h2->GetXaxis()->GetNbins();
  auto ny=h2->GetYaxis()->GetNbins();
  for (auto binx=1;binx<=nx;++binx) {
    for (auto biny=1;biny<=ny;++biny) {
      auto x=h2->GetXaxis()->GetBinCenter(binx);
      auto y=h2->GetYaxis()->GetBinCenter(biny);
      if(cut->IsInside(x,y)==1) {
        int n=h2->GetBinContent(binx,biny);
        //hnew->SetBinContent(binx,biny,h2->GetBinContent(binx,biny));
        if (n<1)
          hnew->Fill(x,y,h2->GetBinContent(binx,biny));
        else
          for (auto i=0;i<int(n);++i) hnew->Fill(x,y);
      }
    }
  }

  TString cut_name = cut -> GetName();
  if (!cut_name.IsNull()) {
    hnew -> SetName(Form("%s_%s",hnew->GetName(),cut_name.Data()));
    hnew -> SetTitle(Form("%s_%s",hnew->GetTitle(),cut_name.Data()));
  }

  std::cout<<"cutg->["<<hnew->GetName()<<"],[c:"<<cut->GetName()<<"]->"<<hnew->GetEntries()<<endl;
  return (TH1 *) hnew;
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
