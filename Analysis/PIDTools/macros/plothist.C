/***************************************

  Open root file fn containing histograms
  of distributions of the variables used by
  the LikelihoodPIDProcessor. Plot histograms
  of variable var for all particle types for 
  comparison.

  S. Lukic, Jan. 2016

****************************************/
void plothist(const char* fn, const char *var, double legl=.6, double legr=.9) { 

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

  const int npart = 5;
  const char *particle[npart] = {"electron", "muon", "pion", "kaon", "proton"};
  TH1F *histos[npart] = { NULL };

  for (int ipart=0; ipart<npart; ipart++) {

    const char *histoname = Form("%s%d", var, ipart+1);
    file->GetObject(histoname, histos[ipart]);

    if(!histos[ipart])  {
      cout << "Cannot read histogram " << histoname 
           << " from file " << fn << endl;
      return;
    }
    histos[ipart]->SetLineColor(ipart+2);
    histos[ipart]->SetMinimum(1.);
  }

  TCanvas *c = new TCanvas("canvas", var, 800, 600);

  TH1F* hall = new TH1F(*histos[0]);
  hall->SetName("hall");
  for (int ipart=1; ipart<npart; ipart++) 
    hall->Add(histos[ipart]);
  hall->SetLineColor(kBlack);
  hall->SetMinimum(1.);
  hall->DrawCopy();

  for (int ipart=0; ipart<npart; ipart++) {
    histos[ipart]->Draw("same");
  }

  TLegend *leg = new TLegend(legl, .53, legr, .88);
  leg->SetTextFont(42);
  leg->SetTextSize(tsiz);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4001);

  leg->AddEntry(hall, "all 5 types", "l");
  for (int ipart=0; ipart<npart; ipart++) {
    leg->AddEntry(histos[ipart], particle[ipart], "l");
  }
  leg->Draw();

  c->Print(Form("%s.%s.pdf", fn, var));

  delete hall;

  return;
}
