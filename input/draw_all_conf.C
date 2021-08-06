void draw_all_conf(TString selectString="att2")
{
  TIter nextFile(TSystemDirectory("thisDir",".").GetListOfFiles());
  while (auto file = (TSystemFile *) nextFile())
  {
    TString name = file -> GetName();
    if (file -> IsDirectory()) continue;
    if (!name.EndsWith("conf")&&!name.EndsWith("par")) continue;
    name.ReplaceAll(".conf","");
    name.ReplaceAll(".par","");
    if (!selectString.IsNull() && name.Index(selectString.Data())<0) continue;

    if (name.Index("att")==0) {
      auto par = conf(name);
      auto numIndex = par -> GetParVInt("line_color").size();
      auto hist = new TH2D(name,name+";;index",100,0,18,100,-.999,numIndex-.001);
      auto cvs = canvas(name,"wide");
      draw(hist,cvs);
      hist -> GetYaxis() -> SetNdivisions(1.5*numIndex);
      drawz(hist,cvs);
      for (auto idx=0; idx<numIndex; ++idx) {
        auto graph1 = new TGraphErrors(); att(graph1,idx,name);
        auto graph2 = new TGraphErrors(); att(graph2,idx,name);
        for (auto j=1; j<5; ++j) { graph1 -> SetPoint(graph1->GetN(),j+.1*idx,idx); graph1 -> SetPointError(graph1->GetN()-1,0,.4); }
        for (auto j=5; j<9; ++j) { graph2 -> SetPoint(graph2->GetN(),j+.1*idx,idx); graph2 -> SetPointError(graph2->GetN()-1,0,.4); }
        graph1 -> Draw("samepl");
        graph2 -> Draw("samep3");
        auto tt = new TLatex(10,idx,""); att(tt,idx,name);
        //auto pt = newpt(Form("%d) m(%d,%.1f,%d)  t(%d,%d,%d)",idx,
        auto pt = newpt(Form("marker(%d,%.1f,%d)  text(%d,%d,%d)",
              graph1->GetMarkerStyle(), graph1->GetMarkerSize(), graph1->GetMarkerColor(),
              tt->GetTextFont(),int(tt->GetTextSize()),tt->GetTextAlign()),cvs,9,idx-.5,18,idx+.5);
        pt -> Draw();
      }
    }
    else if (name.Index("nn")>=0)
    {
      auto cvs22 = canvas(name,2,2,name);
      auto hist1 = new TH1D(name+1,name+";x;y",10,0,10); draw(hist1,cvs22->cd(1));
      auto hist2 = new TH1D(name+2,name+";x;y",10,0,10); draw(hist2,cvs22->cd(2));
      auto hist3 = new TH1D(name+3,name+";x;y",10,0,10); draw(hist3,cvs22->cd(3));
      auto hist4 = new TH1D(name+4,name+";x;y",10,0,10); draw(hist4,cvs22->cd(4));
    }
    else {
      auto hist = new TH1D(name,name+";x;y",10,0,10);
      draw(hist,canvas(name,name));
    }
  }
}
