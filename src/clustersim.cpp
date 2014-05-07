#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

#include <sensor.h>
#include <Palette.h>

#include <TProfile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TVirtualPad.h>
#include <TLine.h>
#include <TAxis.h>
#include <TFrame.h>

TRandom3 dado(42);


const int middle = 0;
const int leading = 1;
const int gaelle = 2;
const int marchioroInt = 3;
const int marchioroHalf = 4;
const int marchioroQuarter = 5;
const int stefano = 6;
const int nClusterStyles = 7;
// const int nClusterStyles = 4;
// Not used go here
// const int marchioroInt = 5;

const int nNeighbours = 2;
const double positionSigmaCut = 5;

const int maxCW = 3;

class BadConversion : public std::runtime_error {
public:
  BadConversion(const std::string& s)
    : std::runtime_error(s)
  { }
};

inline std::string stringify(double x)
{
  std::ostringstream o;
  if (!(o << x))
    throw BadConversion("stringify(double)");
  return o.str();
}

int willPalette(int clusterWidth) {
  int color;
  if (clusterWidth==0) color = kBlack;
  else if (clusterWidth==1) color = kRed;
  else if (clusterWidth==2) color=kBlue;
  else if (clusterWidth==3) color = kMagenta;
  else color = clusterWidth+10;
  return color;
}

void fillClusterWidthPlot(Sensor& mySensor, TProfile *clusterWidth, double cone,
			  int nSteps, bool crossTalk = true) {
  int count;
  int total;
  double angle;
  double radAngle;
  double refVar;
  //  int sign=0;
  int aCount, aCharge;
  double impactStep =  mySensor.getPitch() / double (nSteps);
  
  for (angle=90+cone; angle>=90-cone ; angle-=1) {
    //     sign=1;
    //     if (angle>90) {
    //       sign=-1;
    //     }
    radAngle=angle/(180./M_PI); // Angle in readiants
    count=0;
    total=0;
    for (double impact=0; impact < mySensor.getPitch() ; impact+=impactStep) {
      // Now we add the cross-talk if needed
      mySensor.trackHit(impact,radAngle);
      if (crossTalk) mySensor.crossTalk();
      mySensor.clusterize(aCount, aCharge);
      total += aCount;
      count++;
    }
    //refVar=cos(radAngle)/sin(radAngle)/(mySensor.getPitch())*(mySensor.getThickness());
    refVar=cos(radAngle)/sin(radAngle);
    clusterWidth->Fill(refVar,double(total)/double(count));
  }
}

void binaryError(Sensor &mySensor, double angleDeg, int nSteps, TH1D* errHisto, int clusterStyle, 
                 double telescopeResolution, bool crossTalk=true,
                 TProfile* myClusterSize = NULL, TProfile* myEfficiency = NULL) {
  double binaryPosition;
  bool onlyLeading = (clusterStyle==leading);
  int nStrip, clusterWidth;
  int aCharge;
  double pitchStep = 4 * mySensor.getPitch() / double(nSteps);
  double centralPosition;
  double baricenter, baricenter_ext, totalCharge;
  bool efficient;

  for (double impact=0; impact<mySensor.getPitch()*4; impact+=pitchStep) {
    mySensor.trackHit(impact,M_PI/180.*(90+angleDeg)); // -1+2*dado.Rndm()));
    mySensor.addNeighbourNoise(nNeighbours);
    centralPosition = impact + nNeighbours*mySensor.getPitch() + mySensor.getThickness()/2.*tan(M_PI/180.*(angleDeg));
    if (telescopeResolution>0) centralPosition+=dado.Gaus(0, telescopeResolution);
    
    // Now we add the cross-talk if needed
    if (crossTalk) mySensor.crossTalk();
    efficient = mySensor.binaryCluster(nStrip, clusterWidth, onlyLeading, baricenter, baricenter_ext, totalCharge);
    if (myEfficiency) {
      int fillValue = (efficient ? 1 : 0);
      //if (!efficient) myEfficiency->Fill(angleDeg, 0);
      //else myEfficiency->Fill(angleDeg, 1);
      myEfficiency->Fill(angleDeg, fillValue);
      //std::cout << "myEfficiency->Fill(" << angleDeg  << ", " << fillValue << ");" << std::endl;
    }
    if (efficient) { // Cluster was found
      
      if (clusterStyle==gaelle) if (clusterWidth%2==1) clusterWidth--;

      if ( (clusterStyle==middle)||(clusterStyle==leading)||(clusterStyle==gaelle) ) {
        binaryPosition = (nStrip-1+clusterWidth/2.);
      } else if ((clusterStyle==marchioroInt)||(clusterStyle==marchioroHalf)||(clusterStyle==marchioroQuarter)) {
        binaryPosition = baricenter;
        if (clusterStyle==marchioroInt) {
          binaryPosition = ceil(baricenter-.5);
        } else if (clusterStyle==marchioroHalf) {
          binaryPosition = ceil(2*baricenter-.5)/2;
        } else if (clusterStyle==marchioroQuarter) {
          binaryPosition = ceil(4*baricenter-.5)/4;
        } else {
          std::cerr << "Error: I should never be here :)" << std::endl;
        }
      } else if (clusterStyle==stefano) {
        binaryPosition = ceil(4*baricenter_ext-.5)/4;
      } else {
        std::cerr << "Warning: unknown cluster style " << clusterStyle << std::endl;
        binaryPosition = 0;
      }

      // The position is now in micrometers
      binaryPosition*=mySensor.getPitch();
    
      // std::cerr << "Marchioro 1/4 = " << ceil(4*baricenter-.5)/4 << ", stefano 1/4 = " <<  ceil(4*baricenter_ext-.5)/4 << std::endl;
      // std::cerr << "Diff " << ceil(4*baricenter-.5)/4 - ceil(4*baricenter_ext-.5)/4 << " ";
      // mySensor.printDigis();

      // if ((angleDeg==91) && (clusterStyle==gaelle) && ((binaryPosition-centralPosition)<-100)) {
      //   std::cout << "nStrip = " << nStrip << ", clusterWidth = " << clusterWidth << std::endl;
      //   std::cout << "  impact = " << impact << " reconstructed = " << binaryPosition << std::endl;
      // }
      // distro->Fill(impact, binaryPosition);
      
      // Debug:
      // if ((binaryPosition-centralPosition<0)&&(clusterStyle==middle)&&(angleDeg<-28)) {
      //   std::cout << angleDeg << "\t";
      //   mySensor.printDigis();
      // }

      errHisto->Fill(binaryPosition-centralPosition);
      if (myClusterSize) myClusterSize->Fill(angleDeg, clusterWidth);
    } // cluster was found
  } // loop over clusters
}


double binarySignal(Sensor &mySensor, double angleDeg, int nSteps, TH1D* signalHisto, TProfile** clusterSizeHisto,
                  TH1D** residuals, int reduction, double telescopeResolution) {
  bool onlyLeading = false;
  int nStrip, clusterWidth;
  int aCharge;
  bool efficient;
  double baricenter, baricenter_ext, totalCharge;
  double pitchStep = 4 * mySensor.getPitch() / double(nSteps);
  double interStrip, centralPosition;
  double binaryPosition;
  int nClusters=0;
  double totalClusterWidth=0;

  for (double impact=0; impact<mySensor.getPitch()*4; impact+=pitchStep) {
    centralPosition = impact + nNeighbours*mySensor.getPitch() + mySensor.getThickness()/2.*tan(M_PI/180.*angleDeg);
    if (telescopeResolution>0) centralPosition+=dado.Gaus(0, telescopeResolution);
    interStrip=centralPosition/mySensor.getPitch(); // Central position in units of pitch
    interStrip+=.5;
    interStrip-=int(interStrip);                    // Remove the integer number of strips

    mySensor.trackHit(impact,M_PI/180.*(90+angleDeg));
    mySensor.addNeighbourNoise(nNeighbours);
    // centralPosition = impact + 5*mySensor.getPitch() + mySensor.getThickness()/2.*tan(M_PI/180.*(angleDeg));
    
    // Now we add the cross-talk
    mySensor.crossTalk();
    efficient = mySensor.binaryCluster(nStrip, clusterWidth, onlyLeading, baricenter, baricenter_ext, totalCharge);
    // std::cout << "totalCharge = " << totalCharge << "\t, efficient = " << efficient << "\t clusterWidth = " << clusterWidth<< std::endl;
    if (efficient) { // Cluster was found
      signalHisto->Fill(totalCharge);
    } // cluster was found
    
    double value;
    for (int i=0; i<=maxCW; ++i) {
      value = (i==clusterWidth) ? 1. : 0.;
      clusterSizeHisto[i]->Fill(interStrip,value);
    }
    binaryPosition = (nStrip-1+clusterWidth/2.)*mySensor.getPitch();
    if ((clusterWidth>0&&clusterWidth<=maxCW)) {
      //std::cout << "centralPosition = " << centralPosition << ", binaryPosition = " << binaryPosition << std::endl;
      residuals[clusterWidth-1]->Fill( binaryPosition-centralPosition );
    }
    if (clusterWidth>0) {
      nClusters++;
      totalClusterWidth+=clusterWidth;
    }
  } // loop over clusters
  return totalClusterWidth/nClusters;
}
                
void fillEtaFunction(Sensor &mySensor,std::vector<TH1D*> etaFunction,
		     double angleDeg, int nSteps, bool crossTalk = true,
		     TH2D* correlation = NULL,  bool zeroSup=false) {
  double eta;
  int nStrip;
  int aCharge;
  double pitchStep = mySensor.getPitch() / double(nSteps);
  for (double impact=0; impact<mySensor.getPitch() ; impact+=pitchStep) {
    mySensor.trackHit(impact,M_PI/180.*(angleDeg+90));
    // Now we add the cross-talk if needed
    if (crossTalk) mySensor.crossTalk();
    mySensor.clusterize(nStrip, aCharge);
    if (zeroSup) {
      eta=mySensor.getEtaZs();
    } else {
      eta=mySensor.getEta();
    }
    if (nStrip>0) {
      if (correlation!=NULL) correlation->Fill(impact/mySensor.getPitch(), eta);
      // std::cerr << impact/mySensor.getPitch()<<"\t"<< eta << std::endl;
      for (int i=0; i<etaFunction.size(); i++) {
	if (nStrip<=(i+1)) {
	  etaFunction[i]->Fill(eta);
	}
      }
    }
  }
}

void createEtaFuncs(double angle, int clusterBins, std::vector<TH1D*> &histos, std::string tag) {
  histos.clear();
  std::string title;
  std::string name;
  int color;
  TH1D* anEta;

  for (int i=0; i < clusterBins; i++) {
    if (i==(clusterBins-1)) {
      color=0;
      name = "EtaFunction" + stringify(angle) + "Clus" + stringify(i+1) + tag;
      title = "Eta Function / " + stringify(angle)+"deg / all-strip clusters " + tag;
    } else {
      color=(i+1)*2;
      name = "EtaFunction" + stringify(angle) + "Clus" + stringify(i+1) + tag;
      title = "Eta Function / " + stringify(angle)+"deg / nStrip<=" + stringify(i+1) + "cluster " + tag;
    }

    anEta = new TH1D(name.c_str(), title.c_str(), 50, -0.2, +1.2);
    anEta->SetFillColor(color);
    histos.push_back(anEta);
  }
}

void plotFuncs(std::vector<TH1D*> &histos) {
  for (int i=histos.size()-1; i>=0; i--) {
    if (i==(histos.size()-1)) {
      histos[i]->Draw();
    } else {
      histos[i]->Draw("same");
    }
  }
}

void createSignalFuncs(double angle, int clusterBins, std::vector<TH1D*> &histos, std::string tag) {
  histos.clear();
  std::string title;
  std::string name;
  int color;
  TH1D* aSign;

  for (int i=0; i < clusterBins; i++) {
    if (i==(clusterBins-1)) {
      color=0;
      name = "SignalDistrib" + stringify(angle) + "Clus" + stringify(i+1) + tag;
      title = "Signal / " + stringify(angle)+"deg / all-strip clusters " + tag;
    } else {
      color=(i+1)*2;
      name = "SignalDistrib" + stringify(angle) + "Clus" + stringify(i+1) + tag;
      title = "Signal / " + stringify(angle)+"deg / nStrip<=" + stringify(i+1) + "cluster " + tag;
    }

    aSign = new TH1D(name.c_str(), title.c_str(), 150, 0., 300.);
    aSign->SetFillColor(color);
    histos.push_back(aSign);
  }
}


void fillSignalFunction(Sensor &mySensor,std::vector<TH1D*> signalFunction,
			double angleDeg, int nSteps, bool crossTalk=true) {
  int signal;
  int nStrip;
  int aSign;
  //  std::cerr << angleDeg << " " << nSteps << std::endl;
  double pitchStep = mySensor.getPitch() / double(nSteps);
  for (double impact=0; impact<mySensor.getPitch() ; impact+=pitchStep) {
    mySensor.trackHit(impact,M_PI/180.*(angleDeg+90));
    if (crossTalk) mySensor.crossTalk();
    mySensor.clusterize(nStrip, aSign);
    //std::cerr << aSign << "\t";
    for (int i=0; i<signalFunction.size(); i++) {
      if (nStrip<=(i+1)) {
	signalFunction[i]->Fill(aSign);
      }
    }   
  }
}

void doAllPlots () {


  //*************************//
  //                         //
  //      SYSTEM SETUP       //
  //                         //
  //*************************//
  
  // A bit of setup
  srand ( 42 );
  std::cout << "Setting up ROOT environment tools" << std::endl;
  TFile *myFile;
  myFile = new TFile("output.root", "RECREATE","Cluster and eta histograms");
  TCanvas* myCanvas = new TCanvas("c1", "c1",0,0,696,472);
  myCanvas->SetFillColor(0);
  myCanvas->SetBorderSize(2);


  //*************************//
  //                         //
  //      SENSOR SETUP       //
  //                         //
  //*************************//
  
  // Set here all the parameters of the sensor. You can later
  // change one of them to build a different plot
  std::cout << "Instantiating the sensor..." << std::endl;  
  Sensor mySensor;
  mySensor.setChargeDensity(0.3);  // ADC per um (90 adc in 300 um) ergo .3
  mySensor.setThickness(290.);     // 290 um silicon thickness
  mySensor.setPitch(80);           // 80 um  strip pitch
  mySensor.setNoiseRms(3.8);       // 3.8 ADC mean noise
  mySensor.setDriftZone(30.);      // 30 um of charge sharing (15 right and 15 left to middle-strip)
  mySensor.setSeedCut(4.);         // 4 Seeding cut
  mySensor.setNeighbourCut(3.);    // 3 Neighbour strip cut
  mySensor.setTotalCut(5.);        // 5 Total charge cut (unused)
  mySensor.setXTalk(0.15);         // 15% of charge shared with neighbours (7.5% right and 7.5% left)
  mySensor.setLandauWidth(0.08);   // landau width in units of MPV. Set to zero to make signal distribution a delta
  //                                  default value is 0.08

  // If you want to switch off the noise from the simulation, set the NoiseRms to 0
  // the "measured" noise used in the clusterize for the threshold will then be the default one (4.)

  // If you want to switch off the Landau fluctuation from the simulation, set the
  // landau width to 0


  //*************************//
  //                         //
  //     CLUSTER WIDTH       //
  //                         //
  //*************************//
  
  std::cout << "Creating cluster width histograms..." << std::endl;
  // Instantiate cluster width plots as you wish
  TProfile* clusterWidthXt = new TProfile("clusterWidthXt", "Cluster Width (cross-talk)", 100, -1, 1);

  // Fill the cluster width plot using tracks with theta uniformely distributed
  // between 50 and 150 deg (cone=50 deg). We will simulate 1000 different impact points
  // uniformely distributed within the pitch. Strip X-talk is switched on (default)
  fillClusterWidthPlot(mySensor, clusterWidthXt, 50., 1000, true);


  //*************************//
  //                         //
  //      ETA FUNCTION       //
  //                         //
  //*************************//
  
  std::cout << "Creating eta function histograms..." << std::endl;
  std::vector<TH1D*> etaFuncsXt;
  std::vector<TH1D*> etaFuncsZsXt;

  // createEtaFuncs is a TH1D Factory which prepares the distribution
  // of eta function for different strip size clusters. (up to 3 in this case)
  // These plots will be used in a lego-plot style, so 2-strip cluster plot
  // couts all the clusters with nStrip<=2. The last one (3-strip in this case)
  // holds all the clusters (even if nStrip>3). The last string is just a tag.
  // [do not put spaces in the tag, as it will be used in the plot 'name'
  createEtaFuncs(90, 3, etaFuncsXt, "Xt");
  createEtaFuncs(90, 3, etaFuncsZsXt, "ZsXt");

  // You may want to know the correlation between impact parameter and eta
  TH2D* myCorrelation = new TH2D("a", "b", 100, -0.2, 1.2, 100, -0.2, 1.2);
 
  // Eta function are filled here for tracks impinging at angleDeg degrees,
  // using nSteps different impact points uniformely distributed inside the pitch.
  // cross-Talk can be switched on and off with crossTalk.
  // You can also get the eta function computed only on clusters (like there was
  // zero suppression) by setting zeroSup=true.
  // By passing a non-null TH2D* as 'correlation' you obtain the correlation
  // between (impactParameter/pitch) and eta function.

  // fillEtaFunction(mySensor, etaFunction, angleDeg, nSteps, crossTalk = true,
  //		  TH2D* correlation = NULL,  bool zeroSup=false);
  fillEtaFunction(mySensor, etaFuncsXt, 0, 1003, true, myCorrelation, false);
  fillEtaFunction(mySensor, etaFuncsZsXt, 0, 1003, true, NULL, true);




  //*************************//
  //                         //
  //  SIGNAL DISTRIBUTION    //
  //                         //
  //*************************//
  
  std::cout << "Creating signal histograms..." << std::endl;
  std::vector<TH1D*> signalXt;

  createSignalFuncs(90, 3, signalXt, "Xt");

  std::cout << "Filling signal histograms..." << std::endl;
  fillSignalFunction(mySensor, signalXt, 0., 10033);



  //*************************//
  //                         //
  // Here we plot everything //
  //                         //
  //*************************//
  

  // Here comes the fun part, where you plot the nice pictures
  std::cout << "Writing histograms to disk..." << std::endl;
  
  clusterWidthXt->Draw();
  myCanvas->SaveAs("clusterWidth.png");

  plotFuncs(etaFuncsXt);
  myCanvas->SaveAs("etaFunctionXt.png");

  plotFuncs(etaFuncsZsXt);
  myCanvas->SaveAs("etaFunctionZsXt.png");

  plotFuncs(signalXt);
  myCanvas->SaveAs("signalXt.png");

  myCorrelation->Draw();
  myCanvas->SaveAs("eta-vs-impact.png");

  myFile->Write();

  return ;
}

bool readParamInt(std::string command, std::string ident, int& param) {
  int position;
  std::string value;
  bool found=false;

  position=command.find(ident);
  if (position==0) {
    found=true;
    value=command.substr(ident.length());
    param=atoi(value.c_str());
  }
  return found;
}

bool readParam(std::string command, std::string ident, double& param) {
  int position;
  std::string value;
  bool found=false;

  position=command.find(ident);
  if (position==0) {
    found=true;
    value=command.substr(ident.length());
    param=atof(value.c_str());
  }
  return found;
}

bool readParamStr(std::string command, std::string ident, std::string& param) {
  int position;
  std::string value;
  bool found=false;

  position=command.find(ident);
  if (position==0) {
    found=true;
    param=command.substr(ident.length());
  }
  return found;
}

double maxValue(double arr[], int asize) {
  double max = arr[0];
  for (int i=1; i<asize; ++i) {
    if (arr[i]>max) max=arr[i];
  }
  return max;
}

void markReferences(TFrame* myFrame, TAxis* xAxis, double thickness, double pitch, bool doResolution = true) {
  double x1 = myFrame->GetX1();
  double x2 = myFrame->GetX2();
  double y1 = myFrame->GetY1();
  double y2 = myFrame->GetY2();

  if (doResolution) {
    double binaryResolution = pitch/sqrt(12);
    TLine* myReferenceBinary = new TLine(x1, binaryResolution, x2, binaryResolution);
    myReferenceBinary->SetLineColor(kGray);
    myReferenceBinary->SetLineStyle(2);
    myReferenceBinary->Draw();
    myReferenceBinary->SetLineWidth(2);
  }

  TLine* myTick;
  double ticky = (y2-y1)*xAxis->GetTickLength();
  double tickx;
  for (double i=-5; i<5.1; i+=1) {
    tickx = atan(i*pitch/thickness)/3.14159265358979312*180;
    if ((tickx>x1)&&(tickx<x2)) {
      myTick = new TLine(tickx, 0, tickx, ticky);
      myTick->SetLineColor(kRed);
      myTick->SetLineWidth(2);
      myTick->Draw();
    }
  }
}

int main (int argc, char*argv[]) {
  Sensor mySensor;
  std::string command;
  double param;
  int paramInt;
  std::string paramStr;
  bool doHelp=false;
  double myAngle = 0;
  std::string myTag = "";
  double maxAngle_=30;
  double telescopeResolution_ = 0;
  
  if (argc==1) doHelp=true;

  for (int i=1; i<argc; i++) {
    if (readParam(argv[i], "--help", param)) {
      doHelp=true;
    }
    if (readParam(argv[i], "help", param)) {
      doHelp=true;
    }
    if (readParam(argv[i], "-h", param)) {
      doHelp=true;
    }
  }
  
  if (doHelp) {
    std::cout << "Syntax " << argv[0] << " [--option=value] [--option=value] ... [command] [--option=value] ..." << std::endl;
    std::cout << std::endl;
    std::cout << "Parameters you can set [default value]:" << std::endl;
    std::cout << "--chargedensity=value         Charge per crossed um in ADC equivalent [" <<  mySensor.getChargeDensity() <<  "]" << std::endl;
    std::cout << "--thickness=value             Silicon thickness in um [" <<  mySensor.getThickness() <<  "]" << std::endl;
    std::cout << "--pitch=value                 Strip pitch in um [" <<  mySensor.getPitch() <<  "]" << std::endl;
    std::cout << "--noiserms=value              Mean noise in ADC [" <<  mySensor.getNoiseRms() <<  "]" << std::endl;
    std::cout << "--driftzone=value             Charge sharing zone between strips in um [" <<  mySensor.getDriftZone() <<  "]" << std::endl;
    std::cout << "--seedcut=value               Seeding cut [" <<  mySensor.getSeedCut() <<  "]" << std::endl;
    std::cout << "--neighbourcut=value          Neighbour strip cut [" <<  mySensor.getNeighbourCut() <<  "]" << std::endl;
    std::cout << "--totalcut=value              Total charge cut (unused) [" <<  mySensor.getTotalCut() <<  "]" << std::endl;
    std::cout << "--crosstalk=value             Fraction of charge shared with neighbours [" <<  mySensor.getXTalk() <<  "]" << std::endl;
    std::cout << "--landauwidth=value           Landau width in units of MPV. Zero to make signal distribution a delta [" <<  mySensor.getLandauWidth() <<  "]" << std::endl;
    std::cout << "--threshold=value             Threshold for binary cluster finding [" <<  mySensor.getBinaryThreshold() <<  "]" << std::endl;
    std::cout << "--angle=value                 Angle of the track impinging the surface [" <<  myAngle <<  "]" << std::endl;
    std::cout << "--saturation=value            maximum value for the signal in local ADC [" <<  mySensor.getSaturation() <<  "]" << std::endl;
    std::cout << "--adcbits=value               number of ADC bits for the signal [" <<  mySensor.getAdcBits() <<  "]" << std::endl;
    std::cout << "--tag=value                   tags the next plot with a string (only works for binErrorScan now)" << std::endl;
    std::cout << "--maxangle=value              maximum track angle for the binErrorScan [" << maxAngle_ << "]" << std::endl;
    std::cout << "--telescoperesolution=value   tracking telescope resolution in um [" << telescopeResolution_ << "] simulating analysis finite resolution" << std::endl;
    std::cout << std::endl;
    std::cout << "Commands you can send" << std::endl;
    std::cout << "help                          Only prints this help" << std::endl;
    std::cout << "clusterwidth                  Creates cluster width plot" << std::endl;
    std::cout << "etafunction                   Creates eta function distribution and impact point correlation" << std::endl;
    std::cout << "binarySignal                  Creates signal distribution for the binary threshold cut" << std::endl;
    std::cout << "binaryError                   Analysis of the error due to binary clusters" << std::endl;
    std::cout << "binErrorScan                  Analysis of the error due to binary clusters vs. incident angle" << std::endl;
    std::cout << std::endl;
    return 0;
  }

  // A bit of setup
  srand ( 42 );
  //TFile *myFile;
  //myFile = new TFile("output.root", "RECREATE","Cluster and eta histograms");
  TCanvas* myCanvas = new TCanvas("c1", "c1",0,0,696,472);
  myCanvas->SetFillColor(0);
  myCanvas->SetBorderSize(2);
  char filename[1000];
  char filename_root[1000];
  int clusterWidthCount=0;
  int etaFunctionCount=0;
  int binaryErrorCount=0;
  int binarySignalCount=0;

  for (int i=1; i<argc; i++) {
    if (readParam(argv[i], "--chargedensity=", param))       mySensor.setChargeDensity(param);
    if (readParam(argv[i], "--thickness=", param))           mySensor.setThickness(param);
    if (readParam(argv[i], "--pitch=", param))               mySensor.setPitch(param);
    if (readParam(argv[i], "--noiserms=", param))            mySensor.setNoiseRms(param);
    if (readParam(argv[i], "--driftzone=", param))           mySensor.setDriftZone(param);
    if (readParam(argv[i], "--seedcut=", param))             mySensor.setSeedCut(param);
    if (readParam(argv[i], "--neighbourcut=", param))        mySensor.setNeighbourCut(param);
    if (readParam(argv[i], "--totalcut=", param))            mySensor.setTotalCut(param);
    if (readParam(argv[i], "--crosstalk=", param))           mySensor.setXTalk(param);
    if (readParam(argv[i], "--rndmcrosstalk=", param))       mySensor.setRndmXTalk(param);
    if (readParam(argv[i], "--landauwidth=", param))         mySensor.setLandauWidth(param);
    if (readParam(argv[i], "--threshold=", param))           mySensor.setBinaryThreshold(param);
    if (readParam(argv[i], "--angle=", param))               myAngle = param;
    if (readParam(argv[i], "--saturation=", param))          mySensor.setSaturation(param);
    if (readParam(argv[i], "--maxangle=", param))            maxAngle_ = param;
    if (readParam(argv[i], "--telescoperesolution=", param)) telescopeResolution_ = param;
    if (readParamInt(argv[i], "--adcbits=", paramInt))       mySensor.setAdcBits(paramInt);
    if (readParamStr(argv[i], "--tag=", paramStr))           myTag = paramStr;

    if (readParam(argv[i], "clusterwidth", param)) {
      TProfile* clusterWidth = new TProfile("clusterWidth", "Cluster Width", 100, -1, 1);
      fillClusterWidthPlot(mySensor, clusterWidth, 50., 1000, true);

      sprintf(filename, "clusterwidth_%d.png", ++clusterWidthCount);
      clusterWidth->Draw();
      myCanvas->SaveAs(filename);
      delete clusterWidth;
    }

    if (readParam(argv[i], "etafunction", param)) {
      std::vector<TH1D*> etaFuncs;
      createEtaFuncs(90, 3, etaFuncs, "");
      TH2D* myCorrelation = new TH2D("correlation", "Correlation", 100, -0.2, 1.2, 100, -0.2, 1.2);
      fillEtaFunction(mySensor, etaFuncs, 0, 1003, true, myCorrelation, false);

      sprintf(filename, "etafunction_%d.png", ++etaFunctionCount);
      plotFuncs(etaFuncs);
      myCanvas->SaveAs(filename);

      sprintf(filename, "etafunctioncorr_%d.png", etaFunctionCount);
      myCorrelation->Draw();
      myCanvas->SaveAs(filename);
      for (std::vector<TH1D*>::iterator it=etaFuncs.begin(); it!=etaFuncs.end(); it++) {
	delete (*it);
      }
      delete myCorrelation;
    }
    
    if (readParam(argv[i] , "binaryError", param)) {
      sprintf(filename, "binaryError_%3.0f_%d.png", myAngle, ++binaryErrorCount);
      TH1D* myResolutionFull = new TH1D("resolutionFull",
                                        "Binary resolution (full cluster);Reconstructed - average [um]",
                                        100, -2*mySensor.getPitch(), 2*mySensor.getPitch());
      binaryError(mySensor, myAngle, 1000, myResolutionFull, false, true);
      
      TH1D* myResolutionLead = new TH1D("resolutionLead",
                                        "Binary resolution (leading strip);Reconstructed - average [um]",
                                        100, -2*mySensor.getPitch(), 2*mySensor.getPitch());
      binaryError(mySensor, myAngle, 1000, myResolutionLead, true, true);
      
      double max[2];
      max[0] = myResolutionLead->GetMaximum();
      max[1] = myResolutionFull->GetMaximum();
      double maximum = maxValue(max, 2);
      
      myResolutionFull->SetMaximum(maximum);
      myResolutionLead->SetMaximum(maximum);
      myResolutionLead->SetFillColor(kRed);
      myResolutionFull->SetFillColor(kBlue);
      myResolutionLead->Draw();
      myResolutionFull->Draw("same");
      myCanvas->SaveAs(filename);
    }

    if (readParam(argv[i] , "binErrorScan", param)) {

      double minAngle_  = -maxAngle_;
      double angleSteps_ = 2*maxAngle_;

      char filename_clustersize[1000];
      char filename_clustersize_root[1000];

      char filename_efficiency[1000];
      char filename_efficiency_root[1000];


      if (myTag=="") {
        binaryErrorCount++;
        sprintf(filename, "binErrorScan_%d.png", binaryErrorCount);
        sprintf(filename_root, "binErrorScan_%d.root", binaryErrorCount);
        sprintf(filename_clustersize, "clusterWidth_%d.png", binaryErrorCount);
        sprintf(filename_clustersize_root, "clusterWidth_%d.root", binaryErrorCount);
        sprintf(filename_efficiency, "efficiency_%d.png", binaryErrorCount);
        sprintf(filename_efficiency_root, "efficiency_%d.root", binaryErrorCount);
      } else {
        sprintf(filename, "binErrorScan_%s.png", myTag.c_str());
        sprintf(filename_root, "binErrorScan_%s.root", myTag.c_str());
        sprintf(filename_clustersize, "clusterWidth_%s.png", myTag.c_str());
        sprintf(filename_clustersize_root, "clusterWidth_%s.root", myTag.c_str());
        sprintf(filename_efficiency, "efficiency_%s.png", myTag.c_str());
        sprintf(filename_efficiency_root, "efficiency_%s.root", myTag.c_str());
        myTag=="";
      }
      
      // middle
      TProfile* myResolutionFull = new TProfile("resolutionFull", "Binary resolution (full cluster);"
						"Incident angle;"
						"Reconstructed - average [um]", angleSteps_, minAngle_, maxAngle_);

      // leading
      TProfile* myResolutionLead = new TProfile("resolutionLead", "Binary resolution (leading cluster);"
						"Incident angle;"
						"Reconstructed - average [um]", angleSteps_, minAngle_, maxAngle_);

      // gaelle
      TProfile* myResolutionGaelle = new TProfile("resolutionGaelle", "Binary resolution (Gaelle-style);"
                                                  "Incident angle;"
                                                  "Reconstructed - average [um]", angleSteps_, minAngle_, maxAngle_);

      // marchioro integer strips
      TProfile* myResolutionMarchioroInt = new TProfile("resolutionMarchioroInt", "TOT 1 strip resolution (Marchioro-style);"
                                                        "Incident angle;"
                                                        "Reconstructed - average [um]", angleSteps_, minAngle_, maxAngle_);
      
      // marchioro half strips
      TProfile* myResolutionMarchioroHalf = new TProfile("resolutionMarchioroHalf", "TOT 1/2 strip resolution (Marchioro-style);"
                                                         "Incident angle;"
                                                         "Reconstructed - average [um]", angleSteps_, minAngle_, maxAngle_);
      
      // marchioro quarter strips
      TProfile* myResolutionMarchioroQuarter = new TProfile("resolutionMarchioroQuarter", "TOT 1/4 strip resolution (Marchioro-style);"
                                                            "Incident angle;"
                                                            "Reconstructed - average [um]", angleSteps_, minAngle_, maxAngle_);
      
      // stefano quarter strips
      TProfile* myResolutionStefano = new TProfile("resolutionStefano", "TOT 1/4 strip resolution (Stefano-style);"
                                                   "Incident angle;"
                                                   "Reconstructed - average [um]", angleSteps_, minAngle_, maxAngle_);

      TProfile* myClusterSize = new TProfile("clusterSize", "Binary cluster size;"
                                             "Incident angle;"
                                             "Cluster size [strips]", angleSteps_, minAngle_, maxAngle_);

      TProfile* myEfficiency = new TProfile("clusterEfficiency", "Binary cluster efficiency;"
                                            "Incident angle;"
                                            "Fraction of clusters found", angleSteps_, minAngle_, maxAngle_);

      TProfile* actualClusterSize;

      TProfile* actualEfficiency;

      TProfile* aProfile;
      
      TH1D* resolutionHisto = new TH1D("resolution", "resolution", 10000, -30*mySensor.getPitch(), 30*mySensor.getPitch());      
      
      double angleStep = (maxAngle_-minAngle_)/double(angleSteps_);
      for (double anAngle=minAngle_; anAngle<maxAngle_; anAngle+=angleStep) {
	for (int clusterStyle=0; clusterStyle<nClusterStyles; ++clusterStyle) {
	  if (clusterStyle==middle) aProfile=myResolutionFull;
	  else if (clusterStyle==leading) aProfile=myResolutionLead;
	  else if (clusterStyle==gaelle) aProfile=myResolutionGaelle;
	  else if (clusterStyle==marchioroInt) aProfile=myResolutionMarchioroInt;
	  else if (clusterStyle==marchioroHalf) aProfile=myResolutionMarchioroHalf;
	  else if (clusterStyle==marchioroQuarter) aProfile=myResolutionMarchioroQuarter;
	  else if (clusterStyle==stefano) aProfile=myResolutionStefano;
	  else aProfile=NULL;
          
          actualEfficiency = NULL;
          actualClusterSize = NULL;
          if (clusterStyle==middle) {
            actualClusterSize = myClusterSize;
            actualEfficiency = myEfficiency;
          } //else if (clusterStyle==leading) {
          //}

	  resolutionHisto->Reset();
          binaryError(mySensor, anAngle, 1000, resolutionHisto, clusterStyle, telescopeResolution_, true, actualClusterSize, actualEfficiency); // onlyLeading if i==1
          // if ((anAngle>=90)&&(anAngle<100)&&(clusterStyle==gaelle)) {
          //   TCanvas* c1 = new TCanvas();
          //   resolutionHisto->Draw();
          //   char nome[256];
          //   sprintf(nome, "debug_%.0f.png", anAngle);
          //   c1->SaveAs(nome);
          //   //return 0;
          // }

          // debug
          // if ((resolutionHisto->GetRMS()>30)&&(clusterStyle==middle)) {
          //   myCanvas->cd();
          //   gStyle->SetOptStat(111);
          //   resolutionHisto->Draw();
          //   char title[1000];
          //   sprintf(title, "reso_%.0fdeg", anAngle);
          //   resolutionHisto->SetTitle(title);
          //   myCanvas->SaveAs("debug.png");
          //   myCanvas->SaveAs("debug.root");
          // }


          resolutionHisto->GetXaxis()->UnZoom();
          double thisRms = resolutionHisto->GetRMS();
          double thisMean = resolutionHisto->GetMean();
          resolutionHisto->GetXaxis()->SetRangeUser(thisMean-positionSigmaCut*thisRms, thisMean+positionSigmaCut*thisRms);
	  aProfile->Fill(anAngle, resolutionHisto->GetRMS());
	}
      }

      gStyle->SetOptStat(0);

      myCanvas->SetCanvasSize(800,400);
      myCanvas->cd();
      
      std::map<int, int> styleColor;
      std::map<int, std::string> styleName;
      styleColor[middle]           =  0; // Black
      styleColor[gaelle]           = 10; // Orange
      styleColor[leading]          =  4; // Green
      styleColor[marchioroHalf]    =  1; // Blue
      styleColor[marchioroQuarter] =  3; // Yellow
      styleColor[stefano]          =  9; // Purple
      styleColor[marchioroInt]     =  6; // Light blue


      char totName[256];
      sprintf(totName, "Baricenter (%d bits over %.0f MIPs, ", mySensor.getAdcBits(), mySensor.getSaturation()/(mySensor.getThickness()*mySensor.getChargeDensity()));
      
      styleName[middle]        = "Middle strip (1/2-strip res.)"; // Black
      styleName[gaelle]        = "Middle strip (1-strip res.)"; // Orange
      styleName[leading]       = "Leading strip (1-strip res.)"; // Green
      styleName[marchioroHalf] = std::string(totName)+"1/2 strip res.)"; // Blue
      styleName[marchioroQuarter] = std::string(totName)+"1/4 strip res.)"; // Blue
      styleName[marchioroInt]  = std::string(totName)+"1 strip res.)"; // Light blue
      styleName[stefano]       = std::string(totName)+"1/4 strip res. incl. neighb)"; // Purple

      TLegend* myLegend = new TLegend(0.7, 0.7, 0.99, 0.99);

      int aColor;
      std::string plotOption="E1 p";
      TAxis* myAxis = NULL;
      for (int clusterStyle=0; clusterStyle<nClusterStyles; ++clusterStyle) {
        aColor=Palette::Color(styleColor[clusterStyle]);

	if (clusterStyle==middle) aProfile=myResolutionFull;
	else if (clusterStyle==leading) aProfile=myResolutionLead;
	else if (clusterStyle==gaelle) aProfile=myResolutionGaelle;
        else if (clusterStyle==marchioroInt) aProfile=myResolutionMarchioroInt;
        else if (clusterStyle==marchioroHalf) aProfile=myResolutionMarchioroHalf;
        else if (clusterStyle==marchioroQuarter) aProfile=myResolutionMarchioroQuarter;
        else if (clusterStyle==stefano) aProfile=myResolutionStefano;
	else aProfile=NULL;

	aProfile->SetMarkerStyle(8);
	aProfile->SetMarkerColor(aColor);
	aProfile->SetLineColor(aColor);
	aProfile->SetMinimum(0);
	aProfile->SetMaximum(50); // maxplotaaa
	
	aProfile->Draw(plotOption.c_str());
	plotOption="E1 p same";
        myAxis=aProfile->GetXaxis();

        myLegend->AddEntry(aProfile, styleName[clusterStyle].c_str(), "p");
      }
      myLegend->Draw();

      myCanvas->Update();
      TFrame* myFrame = myCanvas->GetFrame();
      markReferences(myFrame, myAxis, mySensor.getThickness(), mySensor.getPitch());
      
      myCanvas->SaveAs(filename);
      myCanvas->SaveAs(filename_root);

      myClusterSize->Draw();
      myClusterSize->SetMinimum(0);
      myClusterSize->SetMaximum(4);
      myCanvas->Update();
      myFrame = myCanvas->GetFrame();
      markReferences(myFrame, myAxis, mySensor.getThickness(), mySensor.getPitch(), false);
      myCanvas->SaveAs(filename_clustersize);
      myCanvas->SaveAs(filename_clustersize_root);      


      myEfficiency->SetLineColor(Palette::Color(1));
      myEfficiency->SetLineWidth(2);
      myEfficiency->Draw();
      myEfficiency->SetMinimum(0);
      myEfficiency->SetMaximum(1);
      myCanvas->Update();
      myFrame = myCanvas->GetFrame();
      markReferences(myFrame, myAxis, mySensor.getThickness(), mySensor.getPitch(), false);
      myCanvas->SaveAs(filename_efficiency);
      myCanvas->SaveAs(filename_efficiency_root);      


      // for (int i=1; i<3; ++i)
      // delete myCanvas->GetPad(i);
    }
    

    if (readParam(argv[i] , "binarySignal", param)) {

      int nTracks = 15000;
      int reduction = 200;

      TH1D* signalHisto = new TH1D("signal", "Signal;Charge [ADC];Counts", nTracks/reduction, 0, 4*mySensor.getThickness()*mySensor.getChargeDensity());
      
      TH1D* residualHistos[maxCW];
      TProfile* clusterSizeHisto[maxCW+1];
      for (int i=0; i<=maxCW; ++i) {
        char name[256];
        char title[256];
        sprintf(name, "clusterSizeHisto%d", i);
        sprintf(title, "Clusters with size = %d;Interstrip position;Fraction", i);
        clusterSizeHisto[i] = new TProfile(name, title, nTracks/reduction, 0, 1);
      }
      double sigmaCut=3;
      double residualsCut = sigmaCut*mySensor.getPitch()/sqrt(12);
      for (int i=0; i<maxCW; ++i) {
        char name[256];
        char title[256];
        sprintf(name, "residualHisto%d", i);
        sprintf(title, "Residuals for CW = %d;Residuals [um];Counts", i);
        residualHistos[i] = new TH1D(name, title, nTracks/reduction, -residualsCut, residualsCut);
      }
      
      double avCW = binarySignal(mySensor, myAngle, nTracks, signalHisto, clusterSizeHisto, residualHistos, reduction, telescopeResolution_);
      std::cout << "Thr = " << mySensor.getBinaryThreshold() << " Angle = " << myAngle << " <CW> = " << avCW << std::endl;

      gStyle->SetOptStat(0);
      myCanvas->SetCanvasSize(800,400);
      myCanvas->cd();

      if (myTag=="") {
        binarySignalCount++;
        sprintf(filename, "residuals_%d.png", binarySignalCount);
        sprintf(filename_root, "residuals_%d.root", binarySignalCount);
      } else {
        sprintf(filename, "residuals_%s.png", myTag.c_str());
        sprintf(filename_root, "residuals_%s.root", myTag.c_str());
      }

      char plotOption[30];
      int color;
      sprintf(plotOption, "");
      double maxHisto=0;
      for (int i=0; i<maxCW; ++i) if (residualHistos[i]->GetMaximum()>maxHisto) maxHisto=residualHistos[i]->GetMaximum();
      for (int i=0; i<maxCW; ++i) {
        residualHistos[i]->SetMinimum(0);
        residualHistos[i]->SetMaximum(maxHisto);
        color = willPalette(i+1);
        residualHistos[i]->SetLineColor(color);
        residualHistos[i]->Draw(plotOption);
        sprintf(plotOption, "same");
        std::cout << "Residuals RMS (" <<myTag << ", CW = " << i+1 << ") = " << residualHistos[i]->GetRMS() << std::endl;
      }
      myCanvas->SetCanvasSize(600,600);
      myCanvas->SaveAs(filename);
      myCanvas->SaveAs(filename_root);      

      myCanvas->SetCanvasSize(800,400);
      if (myTag=="") {
        sprintf(filename, "binarySignal_%d.png", binarySignalCount);
        sprintf(filename_root, "binarySignal_%d.root", binarySignalCount);
      } else {
        sprintf(filename, "binarySignal_%s.png", myTag.c_str());
        sprintf(filename_root, "binarySignal_%s.root", myTag.c_str());
      }

      signalHisto->SetLineWidth(2);
      signalHisto->SetLineColor(Palette::Color(1));
      signalHisto->Draw();
      myCanvas->SaveAs(filename);
      myCanvas->SaveAs(filename_root);      

      if (myTag=="") {
        sprintf(filename, "clusterWidthImpact_%d.png", binarySignalCount);
        sprintf(filename_root, "clusterWidthImpact_%d.root", binarySignalCount);
      } else {
        sprintf(filename, "clusterWidthImpact_%s.png", myTag.c_str());
        sprintf(filename_root, "clusterWidthImpact_%s.root", myTag.c_str());
      }

      myCanvas->SetCanvasSize(600, 600);
      sprintf(plotOption, "HIST");
      for (int i=0; i<=maxCW; ++i) {
        clusterSizeHisto[i]->SetMinimum(0);
        clusterSizeHisto[i]->SetMaximum(1.1);
        clusterSizeHisto[i]->SetLineColor(willPalette(i));
        clusterSizeHisto[i]->Draw(plotOption);
        sprintf(plotOption, "HIST same");
      }
      myCanvas->SaveAs(filename);
      myCanvas->SaveAs(filename_root);      

      myCanvas->SetCanvasSize(800,400);
      myTag=="";

    }
    
    
  }
}
