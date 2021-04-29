#ifndef ejungwoo_HH
#define ejungwoo_HH

namespace ejungwoo
{
  class binning; ///< 1-dimensional binning class
  class binning2; ///< 2-dimensional binning class
  typedef KBParameterContainer parContainer; ///< parameter container (for configuration) class.

  class parCanvas { public:
      int font, menuSize, padSize[2], marginInn[4], marginTop[4], marginDiv[2] ,axisNdivision[3], titleAlign[4], titleColor[4];
      bool removeInnerPadAxis, removeInnerMainTitle, removeInnerZaxis, titleRotate[4];
      double tickSize[3], labelSize[3], labelOffset[3], titleSize[4], titleOffset[4];
  };
  class parAttribute { public:
      vector<int> colorArrayMain, colorArrayFaint, colorArrayDark;
  };


  TClonesArray *fParameterArray = nullptr;
  int fCanvasIndex = 0;
  TObjArray *fCanvasArray = nullptr;

  parContainer *conf(const char *nameConf);

  parCanvas fParCvs;
  const char *fNameCanvasConf = "hist1";
  void setCanvasPar(parContainer *par);

  parAttribute fParAtt;
  const char *fNameAttributeConf = "att1";
  void setAttributePar(parContainer *par);


  TCanvas* canvas(const char *nameCvs, int nx, int ny, const char *nameConf="");
  TCanvas* canvas(const char *nameCvs, const char *nameConf="") { return canvas(nameCvs, 1, 1, nameConf); }
  TH1*     make  (TH1*     hist,   TVirtualPad *vpad, int idx=0, const char *drawOption="");
  TH1*     draw  (TH1*     hist,   TVirtualPad *vpad, int idx=0, const char *drawOption="");
  TGaxis*  drawz (TH1*     hist,   TVirtualPad *vpad, int idx=0, const char *titlez="");
  TLegend* make  (TLegend* legend, TVirtualPad *vpad, int idx=0, double x1InRatio=-1, double y1InRatio=-1, double dxInRatio=0, double dyInRatio=0, double marginObj=-1);
  TLegend* draw  (TLegend* legend, TVirtualPad *vpad, int idx=0, double x1InRatio=-1, double y1InRatio=-1, double dxInRatio=0, double dyInRatio=0, double marginObj=-1);
  void saveAll(const char *nameVersion="");


  void att(TObject *obj, const char *nameConf="");


  TString toString(double value);
};

class ejungwoo::binning
{
  public:
    int fN = 0; ///< number of bins
    double fMin = 0; ///< lower bound
    double fMax = 0; ///< upper bound
    double fW = 0; ///< binning space width
    double fValue = 0; ///< value will be set after iteration using next(), back(), nextb()
    int fIdx = 0; ///< index for iteration
    const char *fTitle = "";
    const char *fExpression = "";

  public:
    binning() : fN(0), fMin(0), fMax(0) {}
    binning(int n, double min, double max, const char *ttl="", const char *expression="") : fN(n), fMin(min), fMax(max), fTitle(ttl), fExpression(expression) { init(); }
    ~binning() {}

    void init() { if (fW>0&&fN<1) setW(fW); else if (fN>0&&fW<1) setN(fN); }
    void make(TTree *tree, const char* branchName) {
      if (fN<1) fN = 100;
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

    void setN         (double n)                 { fN = n; fW = (fMax-fMin)/fN; }
    void setMin       (double min)               { fMin = min; fW = (fMax-fMin)/fN; }
    void setMax       (double max)               { fMax = max; fW = (fMax-fMin)/fN; }
    void setW         (double w)                 { fW = w; fN = int((fMax-fMin)/fW); }
    void setTitle     (const char *ttl)          { fTitle = ttl; }
    void setExpression(const char *expression)   { fExpression = expression; }

    void set(double n, double min, double max, const char *ttl="", const char *expression="") {
      fN = n;
      fMin = min;
      fMax = max;
      fW = (fMax-fMin)/fN;
      if (strcmp(ttl,"")!=0) fTitle = ttl;
      if (strcmp(expression,"")!=0) fExpression = expression;
    }

    bool   isNull()               const { if (fN<1||(fMin==0&&fMax==0)) return true; return false; }
    int    getN()                 const { return fN; }
    double getMin()               const { return fMin; }
    double getMax()               const { return fMax; }
    double getW()                 const { return fW; }
    const char *getTitle()        const { return fTitle; }
    const char *getExpression()   const { return fExpression; }

    //iterator
    void reset() { fIdx = 0; }
    bool next()  { if (fIdx>fN-1) return false; fValue = fMin + (fIdx++) * fW + .5 * fW; return true; }
    void end()   { fIdx = fN+1; }
    bool prev()  { if (fIdx<2) return false; fValue = fMin + (fIdx--) * fW + .5 * fW; return true; }

    //iteration values
    int    ii()    const { return fIdx-1; }
    int    bi()    const { return fIdx; }
    double val()   const { return fValue; }
    double low()   const { return lowEdge(fIdx); }
    double high()  const { return highEdge(fIdx); }

    double getFullWidth()         const { return (fMax - fMin); }
    double getFullCenter()        const { return (fMax + fMin)/2.; }
    double lowEdge(int bin= 1)    const { return fMin+(bin-1)*(fMax-fMin)/fN; }
    double highEdge(int bin=-1)   const { if (bin==-1) bin=fN; return fMin+(bin)*(fMax-fMin)/fN; }
    double getCenter(int bin)     const { return (fMin + (bin-.5)*fW); }
    bool   isInside(double value) const { if(value>fMin && value<fMax) return true; return false; }
    int    findBin(double val)    const { return int((val-fMin)/fW); }
    double xByRatio(double ratio) const { return ((1.-ratio)*fMin + ratio*fMax); }

    TH1D *newHist(const char *name, const char *title="") {
      auto titlexy = Form("%s;%s;",title,fTitle);
      return (new TH1D(name,titlexy,fN,fMin,fMax));
    }

    TString print(bool pout=1) const
    {
      TString minString = ejungwoo::toString(fMin);
      TString maxString = ejungwoo::toString(fMax);
      TString wString = ejungwoo::toString(fW);

      TString line = Form("%d,%s,%s",fN,minString.Data(),maxString.Data());
      if (pout)
        cout << line << endl;
      return line;
    }

    void operator=(const binning binn) {
      fN = binn.getN();
      fMin = binn.getMin();
      fMax = binn.getMax();
      fW = binn.getW();
      fTitle = binn.getTitle();
      fExpression = binn.getExpression();
    }

    ejungwoo::binning2 operator*(const binning binn);

    const char *nnmm() { return Form("%s_%d_%.2f_%.2f",fExpression,fN,fMin,fMax); }
};

class ejungwoo::binning2
{
  public :
    const char *fTitleX = "";
    const char *fExpressionX = "";
    int fNX = 0;
    double fMinX = 0;
    double fMaxX = 0;

    const char *fTitleY = "";
    const char *fExpressionY = "";
    int fNY = 0;
    double fMinY = 0;
    double fMaxY = 0;

    binning2(const char *titleX, const char *expX, int nX, double minX, double maxX, const char *titleY,  const char *expY, int nY, double minY, double maxY)
      : fTitleX(titleX), fExpressionX(expX), fNX(nX), fMinX(minX), fMaxX(maxX), fTitleY(titleY) ,fExpressionY(expY), fNY(nY), fMinY(minY), fMaxY(maxY) {}

    TH2D *newHist(const char *nameHist, const char *titleHist="") {
      auto titlexy = Form("%s;%s;%s;",titleHist,fTitleX,fTitleY);
      auto hist = new TH2D(nameHist,titlexy,fNX,fMinX,fMaxX,fNY,fMinY,fMaxY);
      return hist;
    }

    const char *expr()     { return Form("%s:%s",fExpressionY,fExpressionX); }
    const char *namexy()   { return Form("%s%s",fExpressionX,fExpressionY); }
    const char *namexynn() { return Form("%s%d%s%d",fExpressionX,fNX,fExpressionY,fNY); }

    binning bx() { return binning(fNX,fMinX,fMaxX,fTitleX,fExpressionX); }
    binning by() { return binning(fNY,fMinY,fMaxY,fTitleY,fExpressionY); }

    int icvs(int iX, int iY) { return (fNY-iY-1)*fNX+iX+1; }
};

ejungwoo::binning2 ejungwoo::binning::operator*(const binning binn) {
  return binning2(fTitle,fExpression,fN,fMin,fMax, binn.getTitle(),binn.getExpression(),binn.getN(),binn.getMin(),binn.getMax());
}


/// Get configuration parameter container from ejungwoo.[nameConf].conf file.
KBParameterContainer *ejungwoo::conf(const char *nameConf)
{
  KBParameterContainer *par = nullptr;
  if (fParameterArray==nullptr) fParameterArray = new TClonesArray("KBParameterContainer",10);
  else par = (KBParameterContainer *) fParameterArray -> FindObject(nameConf);
  if (par==nullptr) {
    par = (KBParameterContainer *) fParameterArray -> ConstructedAt(fParameterArray -> GetEntriesFast());
    par -> AddFile(Form("ejungwoo.%s.conf",nameConf));
    par -> SetName(nameConf);
  }
  return par;
}

void ejungwoo::setCanvasPar(parContainer *par)
{
  fParCvs.font = par -> GetParInt("font");
  fParCvs.menuSize = par -> GetParInt("menu_size");
  fParCvs.removeInnerZaxis = par -> GetParBool("remove_inn_zaxis");
  fParCvs.removeInnerPadAxis = par -> GetParBool("remove_inn_pad_axis");
  fParCvs.removeInnerMainTitle = par -> GetParBool("remove_inn_main_title");
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
}

void ejungwoo::setAttributePar(parContainer *par)
{
  fParAtt.colorArrayMain = par -> GetParVInt("color_main");
  if (par -> CheckPar("color_faint")) fParAtt.colorArrayFaint = par -> GetParVInt("color_faint");
  if (par -> CheckPar("color_dark"))  fParAtt.colorArrayDark  = par -> GetParVInt("color_dark");
}

/// Create canvas with nameCvs and division number nx and ny.
/// Set nx and ny is 1 for single pad canvas.
/// Configuration of the canvas(+ histogram and legned) is set by ejungwoo.[nameConf].conf file.
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

  const char *titleCvs = Form("%s.%d.%d.%d.   %s",nameConf,nx,ny,0,nameCvs);
  auto cvs = new TCanvas(nameCvs, titleCvs, fullSize[0], fullSize[1]);
  cvs -> SetMargin(marginRatioTop[0],marginRatioTop[1],marginRatioTop[2],marginRatioTop[3]);

  if (fCanvasArray==nullptr) fCanvasArray = new TObjArray();
  fCanvasArray -> Add(cvs);

  if (nx>1||ny>1)
  {
    int countPad = 0;

    for (auto iy=0; iy<ny; ++iy) {
      for (auto ix=0; ix<nx; ++ix)
      {
        cvs -> cd();
        countPad++;

        double y2 = 1 - (iy>0)*(padRatioInnEdge1[3] + (iy-1)*padRatioInn[1]) - iy*marginRatioDiv[1];
        double y1 = y2 - padRatioInn[1];
        if (iy==0&&iy==ny-1) y1 = y2 - padRatioInnEdge2[1];
        else if (iy==0)      y1 = y2 - padRatioInnEdge1[3];
        else if (iy==ny-1)   y1 = y2 - padRatioInnEdge1[2];

        double x1 = (ix>0)*(padRatioInnEdge1[0] + (ix-1)*padRatioInn[0]) + ix*marginRatioDiv[0];
        double x2 = x1 + padRatioInn[0];
        if (ix==0&&ix==nx-1) x2 = x1 + padRatioInnEdge2[0];
        else if (ix==0)      x2 = x1 + padRatioInnEdge1[0];
        else if (ix==nx-1)   x2 = x1 + padRatioInnEdge1[1];

        const char *namePad = Form("%s_%d",cvs->GetName(),countPad);
        const char *titlePad = Form("%d",countPad);
        if (ix==0) titlePad = Form("l%s",titlePad);
        if (ix==nx-1) titlePad = Form("r%s",titlePad);
        if (iy==ny-1) titlePad = Form("t%s",titlePad);
        if (iy==0) titlePad = Form("b%s",titlePad);
        titlePad = Form("%s.%d.%d.%d   %s",nameConf,nx,ny,countPad,titlePad);

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
  }

  return cvs;
}

/// save all canvases created so far in [nameVersion]/figures/[cvs-name].png and [nameVersion]/pdf/[cvs-name].pdf files
void ejungwoo::saveAll(const char *nameVersion) {
  if (strcmp(nameVersion,"")!=0) {
    gSystem -> mkdir(nameVersion);
    gSystem -> mkdir(Form("%s/figures",nameVersion));
    gSystem -> mkdir(Form("%s/pdf",nameVersion));
    auto numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      cvs -> SaveAs(Form("%s/figures/%s.png",nameVersion,cvs->GetName()));
      cvs -> SaveAs(Form("%s/pdf/%s.pdf",nameVersion,cvs->GetName()));
    }
  }
  else {
    gSystem -> mkdir(Form("figures"));
    gSystem -> mkdir(Form("pdf"));
    auto numCanvases = fCanvasArray -> GetEntriesFast();
    for (auto iCanvas=0; iCanvas<numCanvases; ++iCanvas) {
      auto cvs = (TCanvas *) fCanvasArray -> At(iCanvas);
      cvs -> SaveAs(Form("figures/%s.png",cvs->GetName()));
      cvs -> SaveAs(Form("pdf/%s.pdf",cvs->GetName()));
    }
  }
}

/// Make hist fit to the pad. If idx>0, pad where legend is drawn is selected by vpad->cd(idx). 
/// Draw histogram before make(TH1* hist, TVirtualPad *vpad) to apply main title attribute.
TH1 *ejungwoo::make(TH1 *hist, TVirtualPad *vpad, int idx, const char *drawOption)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd(idx);

  int nx = 1;
  int ny = 1;
  int idx0 = idx;
  const char *nameConfCvs = "";
  if (padi!=nullptr) {
    TString titleCvs = padi -> GetTitle();
    auto tokens = titleCvs.Tokenize(".");
    nameConfCvs = ((TObjString *) tokens->At(0))->GetString();
    nx = ((TObjString *) tokens->At(1))->GetString().Atoi();
    ny = ((TObjString *) tokens->At(2))->GetString().Atoi();
    if (idx==0) idx0 = ((TObjString *) tokens->At(3))->GetString().Atoi();
  }

  const char *nameConf = ((padi!=nullptr)?nameConfCvs:fNameCanvasConf);
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
  if (!(nx==1&&ny==1)&&idx0>0)
  {
    for (iy=0; iy<ny; ++iy) {
      bool breaky = false;
      for (ix=0; ix<nx; ++ix) {
        countPad++;
        if (countPad==idx0) {
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

  for (auto ixyz : {0,1,2}) {
    TAxis *axis;
         if (ixyz==0) axis = (TAxis *) hist -> GetXaxis();
    else if (ixyz==1) axis = (TAxis *) hist -> GetYaxis();
    else if (ixyz==2) axis = (TAxis *) hist -> GetZaxis();
    if (fParCvs.titleAlign[ixyz]==2) axis -> CenterTitle();
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

  if (!(nx==1&&ny==1)&&idx0>0)
  {
    if (ix==0&&ix==nx-1) {}
    else  {
      if (ix!=0 && fParCvs.removeInnerPadAxis) {
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
      if (iy!=ny-1&&fParCvs.removeInnerPadAxis) {
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

  return hist;
}

/// Make z axis
TGaxis *ejungwoo::drawz(TH1* hist, TVirtualPad *vpad, int idx, const char *titlez)
{
  auto pad = (TPad *) vpad;
  auto padi = pad -> cd(idx);

  int nx = 1;
  int ny = 1;
  int idx0 = idx;
  const char *nameConfCvs = "";
  if (padi!=nullptr) {
    TString titleCvs = padi -> GetTitle();
    auto tokens = titleCvs.Tokenize(".");
    nameConfCvs = ((TObjString *) tokens->At(0))->GetString();
    nx = ((TObjString *) tokens->At(1))->GetString().Atoi();
    ny = ((TObjString *) tokens->At(2))->GetString().Atoi();
    if (idx==0) idx0 = ((TObjString *) tokens->At(3))->GetString().Atoi();
  }

  const char *nameConf = ((padi!=nullptr)?nameConfCvs:fNameCanvasConf);
  auto par = conf(nameConf);
  setCanvasPar(par);

  int countPad = 0, ix=0, iy=0;
  if (!(nx==1&&ny==1)&&idx0>0)
  {
    for (iy=0; iy<ny; ++iy) {
      bool breaky = false;
      for (ix=0; ix<nx; ++ix) {
        countPad++;
        if (countPad==idx0) {
          breaky = true;
          break;
        }
      }
      if (breaky)
        break;
    }
  }

  auto zaxis = new TGaxis(padi->GetUxmax(),hist->GetYaxis()->GetXmin(),padi->GetUxmax(), hist->GetYaxis()->GetXmax(), 0, hist->GetYaxis()->GetXmax(), fParCvs.axisNdivision[2] ,"+L");
  zaxis -> SetTitle(titlez);
  if (fParCvs.titleAlign[2]==2) zaxis -> CenterTitle();
  zaxis -> SetBit(TAxis::kRotateTitle);
  zaxis -> SetTitleOffset(hist->GetZaxis()->GetTitleOffset());
  zaxis -> SetTitleSize(fParCvs.titleSize[2]);
  zaxis -> SetTitleColor(fParCvs.titleColor[2]);
  zaxis -> SetLabelSize(fParCvs.labelSize[2]);
  zaxis -> SetLabelOffset(fParCvs.labelOffset[2]);
  zaxis -> SetTickSize(fParCvs.tickSize[2]);
  zaxis -> SetNdivisions(fParCvs.axisNdivision[2]);
  zaxis -> SetTitleFont(fParCvs.font);
  zaxis -> SetLabelFont(fParCvs.font);

  if (!(nx==1&&ny==1)&&idx0>0) {
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
TH1 *ejungwoo::draw(TH1 *hist, TVirtualPad *vpad, int idx, const char *drawOption)
{
  vpad -> cd(idx);
  hist -> Draw(drawOption);
  make(hist,vpad,idx,drawOption);
  return hist;
}

/// Make legend fit to the pad. If idx>0, pad where legend is drawn is selected by vpad->cd(idx). 
/// @param legend
/// @param vpad
/// @param idx
/// @param x1InRatio  Set x1 (legend left, x-position) by the ratio inside the histogram frame (0~1). If x1InRatio<0 (by default), the legend position is set to stick to the right most frame.
/// @param y1InRatio  Set y1 (legend bottom, y-position) by the ratio inside the histogram frame (0~1). If y1InRatio<0 (by default), the legend position is set to stick to the top most frame.
/// @param dxInRatio  Set dy (legend width) by the ratio inside the histogram frame (0~1).  If dxInRatio<=0 (by default), the legend width is calcuated from the text length.
/// @param dyInRatio  Set dy (legend height) by the ratio inside the histogram frame (0~1). If dyInRatio<=0 (by default), the legend height is calcuated from the number of legend entries.
/// @param marginObj  Set ratio of the space occupied by object in the legend (compared to the descriptions). By default, it is >0.25
TLegend *ejungwoo::make(TLegend *legend, TVirtualPad *vpad, int idx, double x1InRatio, double y1InRatio, double dxInRatio, double dyInRatio, double marginObj)
{
  auto pad = (TPad *) vpad;

  double fLegendTextScale = 0.028;
  double fWidthPerLengthLegend = 1;
  double fWidthDefaultLegend   = fWidthPerLengthLegend*2.5;
  double fLegendHeightPerEntry = 2.5;
  double xOffset = 0;
  double yOffset = 0;

  const char *nameConfCvs = "";
  if (pad!=nullptr) {
    TString titleCvs = pad -> GetTitle();
    auto tokens = titleCvs.Tokenize(".");
    nameConfCvs = ((TObjString *) tokens->At(0))->GetString();
  }

  const char *nameConf = ((pad!=nullptr)?nameConfCvs:fNameCanvasConf);
  auto par = conf(nameConf);
  setCanvasPar(par);

  auto padi = pad -> cd(idx);

  auto x1Frame = 0. + padi -> GetLeftMargin();
  auto x2Frame = 1. - padi -> GetRightMargin();
  auto y1Frame = 0. + padi -> GetBottomMargin();
  auto y2Frame = 1. - padi -> GetTopMargin();
  auto xUnit = x2Frame - x1Frame;
  auto yUnit = y2Frame - y1Frame;

  double dxTextMax = 0.;
  auto dxNorm = TLatex(0,0,"0").GetXsize();
  TIter next_entry(legend->GetListOfPrimitives());
  while (TLegendEntry *le=(TLegendEntry*)next_entry())
  {
    TString legendString = le->GetLabel();
    auto dxText = TLatex(0,0,legendString).GetXsize()/dxNorm;
    if (dxTextMax<dxText) dxTextMax = dxText;
  }
  auto dxLegend = dxInRatio * xUnit;
  auto dyLegend = dyInRatio * yUnit;
  if (dxInRatio<=0) dxLegend = (fWidthDefaultLegend + dxTextMax * fWidthPerLengthLegend) * fLegendTextScale * xUnit;
  if (dxInRatio<=0) dyLegend = (legend -> GetNRows() * fLegendHeightPerEntry) * fLegendTextScale * yUnit;
  double marginLegend = marginObj;
  if (marginLegend<0) {
    marginLegend = 1.-((dxTextMax*fWidthPerLengthLegend*fLegendTextScale)/(dxLegend/xUnit));
    if (marginLegend<0.25) marginLegend = 0.25;
  }

  double x1Legend = (1.-x1InRatio)*x1Frame + x1InRatio*(x2Frame);
  double y1Legend = (1.-y1InRatio)*y1Frame + y1InRatio*(y2Frame);
  if (x1InRatio<0) x1Legend = x2Frame-dxLegend;
  if (y1InRatio<0) y1Legend = y2Frame-dyLegend;
  double x2Legend = x1Legend + dxLegend;
  double y2Legend = y1Legend + dyLegend;

  legend -> SetX1(x1Legend + xOffset*xUnit);
  legend -> SetX2(x2Legend + xOffset*xUnit);
  legend -> SetY1(y1Legend + yOffset*yUnit);
  legend -> SetY2(y2Legend + yOffset*yUnit);
  legend -> SetFillStyle(0);
  legend -> SetBorderSize(0);
  legend -> SetTextFont(fParCvs.font);
  legend -> SetMargin(marginLegend);

  return legend;
}

/// Make legend fit to the pad and draw.
TLegend *ejungwoo::draw(TLegend* legend, TVirtualPad *vpad, int idx, double x1InRatio, double y1InRatio, double dxInRatio, double dyInRatio, double marginObj)
{
  vpad -> cd(idx);
  make(legend,vpad,idx,x1InRatio,y1InRatio,dxInRatio,dyInRatio,marginObj);
  legend -> Draw();
  return legend;
}

void ejungwoo::att(TObject *obj, const char *nameConf)
{
  if (TString(nameConf).IsNull()) nameConf = fNameAttributeConf;
  auto par = conf(nameConf);
  setAttributePar(par);
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

#endif
