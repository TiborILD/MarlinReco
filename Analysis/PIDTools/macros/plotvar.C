/***************************************

  Open root file fn produced by the PIDvarPDF
  processor, plot the distributions of the 
  variable var for all particle types, from 
  lolim to hilim, for comparison.

  S. Lukic, Jan. 2016

****************************************/


void plotvar(const char* fn, const char *var, double lolim, double hilim, double legl=.6, double legr=.9) { 

  double tsiz = 0.05;
  gStyle->SetLineWidth(1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTitleSize(tsiz,"xyz");
  gStyle->SetLabelSize(tsiz,"xyz");
  gStyle->SetTextSize(tsiz);
  gStyle->SetTitleOffset(1.28,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(8.8,6.6);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetOptLogy();

  TFile *file = new TFile(fn);
  if(file->IsOpen()) cout << "Opened file.\n";
  else { cout << "Cannot open file " << fn << " \n"; return ; }

  TTree *tree = NULL;
  file->GetObject("varTree", tree);
  if(!tree)  {
    cout << "Cannot read tree from file " << fn << endl;
    return;
  }

  TCanvas *c = new TCanvas("canvas", var, 800, 600);

  const char *condel = "abs(truePDG)==11";
  const char *condmu = "abs(truePDG)==13";
  const char *condpi = "abs(truePDG)==211";
  const char *condka = "abs(truePDG)==321";
  const char *condpro = "abs(truePDG)==2212";
  const char *condall = Form("(%s||%s||%s||%s||%s)", condel, condmu,
                             condpi, condka, condpro);

  TH1F* hall = new TH1F("hall", Form("%s;%s;", var, var), 
                         200, lolim, hilim);
  tree->Project("hall", var, Form("isReconstructed&&%s", condall));
  hall->SetLineWidth(1);
  hall->SetLineColor(kBlack);

  TH1F* hel = new TH1F("hel", Form("%s;%s;", var, var), 
                         200, lolim, hilim);
  tree->Project("hel", var, Form("isReconstructed&&%s", condel));
  hel->SetLineWidth(1);
  hel->SetLineColor(kBlue);

  TH1F* hmu = new TH1F("hmu", Form("%s;%s;", var, var), 
                         200, lolim, hilim);
  tree->Project("hmu", var, Form("isReconstructed&&%s", condmu));
  hmu->SetLineWidth(1);
  hmu->SetLineColor(kRed+2);

  TH1F* hpi = new TH1F("hpi", Form("%s;%s;", var, var), 
                         200, lolim, hilim);
  tree->Project("hpi", var, Form("isReconstructed&&%s", condpi));
  hpi->SetLineWidth(1);
  hpi->SetLineColor(kYellow+1);

  TH1F* hka = new TH1F("hka", Form("%s;%s;", var, var), 
                         200, lolim, hilim);
  tree->Project("hka", var, Form("isReconstructed&&%s", condka));
  hka->SetLineWidth(1);
  hka->SetLineColor(kGreen+2);

  TH1F* hpro = new TH1F("hpro", Form("%s;%s;", var, var), 
                         200, lolim, hilim);
  tree->Project("hpro", var, Form("isReconstructed&&%s", condpro));
  hpro->SetLineWidth(1);
  hpro->SetLineColor(kGray);

//  hall->SetMinimum(hmu->GetMaximum()/100);
  hall->SetMinimum(1);
  hall->Draw();
  hel->Draw("same");
  hmu->Draw("same");
  hpi->Draw("same");
  hka->Draw("same");
  hpro->Draw("same");

  TLegend *leg = new TLegend(legl, .53, legr, .88);
  leg->SetTextFont(42);
  leg->SetTextSize(tsiz);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4001);

  leg->AddEntry(hall, "all 5 types", "l");
  leg->AddEntry(hel, "electrons", "l");
  leg->AddEntry(hmu, "muons", "l");
  leg->AddEntry(hpi, "pions", "l");
  leg->AddEntry(hka, "kaons", "l");
  leg->AddEntry(hpro, "protons", "l");
  leg->Draw();

  c->Print(Form("%s.%s.pdf", fn, var));

  return;
}
