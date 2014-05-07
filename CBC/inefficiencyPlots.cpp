
const int nMv = 4;
int mvList[nMv] = {530, 570, 610, 650};
TProfile* mvProfiles[nMv];

void inefficiencyPlots() {
  char filename[1024];
  TFile* myFile;
  TCanvas* myCanvas;
  TObject* myObject;
  for (int i=0; i<nMv; ++i) {
    sprintf(filename, "clusterWidthImpact_TB2012_%dmV_0deg.root", mvList[i]);
    myFile = new TFile(filename, "READ");
    if (myFile) {
      std::cout << filename << " opened"<< std::endl;
      myCanvas = (TCanvas*) myFile->GetObjectChecked("c1", "TCanvas");
      if (myCanvas) {
        myObject = myCanvas->GetPrimitive("clusterSizeHisto0");
        if ((myObject)&&(std::string(myObject->ClassName())=="TProfile")) {
          mvProfiles[i] = new TProfile(*(TProfile*) myObject);
          std::cout << mvProfiles[i]->GetName() << std::endl;
        }
      }
    }
  }

  myCanvas = new TCanvas("myCanvas", "Inefficiency", 600, 600);
  gStyle->SetOptStat(0);
  std::string plotOpt = "l HIST";
  for (int i=0; i<nMv; ++i) {
    for (int iBin=1; iBin<=mvProfiles[i]->GetNbinsX(); ++iBin) {
      mvProfiles[i]->SetBinContent(iBin,(1-mvProfiles[i]->GetBinContent(iBin))*
                                   mvProfiles[i]->GetBinEntries(iBin) );
    }
    std::cout << mvProfiles[i]->GetName() << std::endl;
    if (i<=1) mvProfiles[i]->SetLineColor(kRed);
    else if (i==2) mvProfiles[i]->SetLineColor(kBlue);
    else if (i==3) mvProfiles[i]->SetLineColor(kGreen);
    mvProfiles[i]->SetTitle(";Interstrip Position [pitch];Efficiency");
    mvProfiles[i]->Draw(plotOpt.c_str());
    mvProfiles[i]->GetYaxis()->SetRangeUser(0,1);
    plotOpt="same l HIST";
  }
  myCanvas->SaveAs("eff.png");
}
