#ifndef _SENSOR_H_
#define _SENSOR_H_

#include <vector>
#include <time.h>
#include <TRandom3.h>
#include <cmath>
#include <cstdlib>

#define RANDOM_SEED 42

typedef std::vector<double> digiV;

extern inline int maxValue(int a, int b) {
  return a > b ? a : b;
}

// The following parameter determines how many tracks left and right the good one
// will be generated to simulate the charge drift
#define DRIFT_WING 1 // n left and n right
#define N_TRACCE ((2*DRIFT_WING)+1)

class Sensor {
 public:
  enum { blockSharing, linearSharing, nSharing };
  
 private:
  TRandom3 myDice;
  int trackHitElement_old(double impact, double theta, digiV &result);
  void trackHitElement_old2(double& x, double& y, double theta, double& prevCharge, double& thisCharge, double& nextCharge);
  void trackHitElement(double& x, double& y, double theta, double& prevCharge, double& thisCharge, double& nextCharge);
  double wallHit(double x, double y, double theta, double width, double height, double& x1, double& y1);
  double landau(double mpv);
  double stripNoise();
  digiV myDigis_;
  double pitch_;
  double chargeDensity_;
  double saturation_;
  int adcBits_;
  double driftZone_;
  double transverseDrift_;
  double thickness_;
  double eta_;
  double etaZs_;
  double seed_cut_;
  double neighbour_cut_;
  double total_cut_;
  double noise_rms_;
  double xTalk_l_;
  double xTalk_c_;
  double xTalk_r_;
  double rndmXTalk_l_;
  double rndmXTalk_c_;
  double rndmXTalk_r_;
  double binaryThreshold_;
  unsigned int sharingMethod_;

  // Landau with parameter expressed
  // in units of MPV
  double landauWidth_;

  // Default parameters in the sensor simulation
  static const double defaultPitch           = 80.;
  static const double defaultTransverseDrift = 0.;   // Obsolete. do not use
  static const double defaultDriftZone       = 10.;
  static const double defaultThickness       = 200.;
  static const double defaultChargeDensity   = 0.3;
  static const double defaultSaturation      = 135;
  static const int defaultAdcBits            = 4;
  static const double defaultNoiseRms        = 4.;
  static const double defaultSeedCut         = 4.;
  static const double defaultNeighbourCut    = 3.;
  static const double defaultTotalCut        = 5.;   // Unused
  static const double defaultRndmXTalkL      = 0.00; // Fraction of charge that can randomly cross talk between strips
  static const double defaultRndmXTalkC      = 1.00; // Fraction of charge that can randomly cross talk between strips
  static const double defaultRndmXTalkR      = 0.00; // Fraction of charge that can randomly cross talk between strips
  static const double defaultXTalkL          = 0.10; // Cross talk between strips
  static const double defaultXTalkC          = 0.80; // each strip shares 20% of its charge
  static const double defaultXTalkR          = 0.10; // with the neighbours TODO: make this adjustable
  static const double defaultLandauWidth     = 0.08; // width is MPV/12.5
  static const double defaultBinaryThreshold = 16; // width is MPV/12.5
  static const unsigned int defaultSharingMethod = linearSharing; // All charge released proportionally to the distance (within charge sharing zone)


  // Simulation parameters
  static const int cycle_gauss   = 5; // Number of flat distributions used to simulate a Gaussian

 public:
  // Constructors and destructor
  ~Sensor();
  Sensor();
  Sensor(double thickness, double pitch, double chargeDensity, double transverseDrift);
  Sensor(double thickness, double pitch, double chargeDensity, double transverseDrift, double, double, double, double);


  // Creates digi correspondent to a specific angle/impact point
  int trackHit_old(double impact, double theta); // old: deprecated
  int trackHit(double impact, double theta);

  // pure binary cluster search
  bool binaryCluster(int& nStrip, int& clusterWidth, const bool& onlyLeading, double& baricenter, double& baricenter_ext, double& totalCharge);


  // Compute cluster on the present digi
  // it reports the number of strips involved and the total
  // (digitised) charge of the cluster
  void clusterize(int& nClus, int& totCharge);

  // Modifies the digi according to the cross-talk
  void crossTalk();

  // Get the eta produced by the last clusterizing
  // Computed on raw data, or zero-suppressed
  double getEta();
  double getEtaZs();

  // Parameter get
  double getNoiseRms()        { return noise_rms_; };
  double getThickness()       { return thickness_; };
  double getPitch()           { return pitch_; };
  double getChargeDensity()   { return chargeDensity_; };
  double getSaturation()      { return saturation_; };
  int getAdcBits()            { return adcBits_; };
  double getTransverseDrift() { return transverseDrift_; };   // obsolete. do not use
  double getDriftZone()       { return driftZone_*2; };
  double getSeedCut()         { return seed_cut_; };
  double getNeighbourCut()    { return neighbour_cut_; };
  double getTotalCut()        { return total_cut_; };
  double getXTalk()           { return (1-xTalk_c_); };
  double getRndmXTalk()       { return (1-rndmXTalk_c_); };
  double getLandauWidth()     { return landauWidth_; }
  double getBinaryThreshold() { return binaryThreshold_; }

  // Parameter set
  bool setSharingMethod(unsigned int newValue);
  void setNoiseRms(double newValue)        { noise_rms_ = newValue; };
  void setThickness(double newValue)       { thickness_ = newValue; };
  void setPitch(double newValue)           { pitch_ = newValue; };
  void setChargeDensity(double newValue)   { chargeDensity_ = newValue; };
  void setSaturation(double newValue)      { saturation_ = newValue; };
  void setAdcBits(int newValue)            { adcBits_ = newValue; };
  void setTransverseDrift(double newValue) { transverseDrift_ = newValue; }; // obsolete. do not use
  void setDriftZone(double newValue)       { driftZone_ = newValue/2; };
  void setSeedCut(double newValue)         { seed_cut_ = newValue; };
  void setNeighbourCut(double newValue)    { neighbour_cut_ = newValue; };
  void setTotalCut(double newValue)        { total_cut_ = newValue; };
  void setXTalk(double newValue) {
    if ((newValue>=0)&&(newValue<=1)) {
      xTalk_c_ = 1 - newValue;
      xTalk_r_ = newValue/2.;
      xTalk_l_ = newValue/2.;
    }
  }
  void setRndmXTalk(double newValue) {
    if ((newValue>=0)&&(newValue<=1)) {
      rndmXTalk_c_ = 1 - newValue;
      rndmXTalk_r_ = newValue/2.;
      rndmXTalk_l_ = newValue/2.;
    }
  }
  void setLandauWidth(double newValue)     { landauWidth_ = newValue; }
  double setBinaryThreshold(const double& newValue ) { binaryThreshold_ = newValue; }
  void addNeighbourNoise(int nStrips);

  void printDigis() {
    for (int i=0; i<myDigis_.size(); ++i) std::cerr << myDigis_[i] <<", ";
    std::cerr << std::endl;
  }
};



#endif
