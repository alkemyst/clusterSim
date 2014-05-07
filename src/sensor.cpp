#include <iostream>
#include <math.h>
#include <sensor.h>

#define _CLUSTERIZE_DEBUG_
#undef _CLUSTERIZE_DEBUG_

Sensor::~Sensor() {
  //  delete etaFunction_;
}

Sensor::Sensor() {
  noise_rms_       = defaultNoiseRms;
  seed_cut_        = defaultSeedCut;
  neighbour_cut_   = defaultNeighbourCut;
  total_cut_       = defaultTotalCut;

  thickness_       = defaultThickness;
  pitch_           = defaultPitch;
  chargeDensity_   = defaultChargeDensity;
  saturation_      = defaultSaturation;
  adcBits_         = defaultAdcBits;
  transverseDrift_ = defaultTransverseDrift;
  driftZone_       = defaultDriftZone;

  xTalk_l_         = defaultXTalkL;
  xTalk_c_         = defaultXTalkC;
  xTalk_r_         = defaultXTalkR;

  rndmXTalk_l_     = defaultRndmXTalkL;
  rndmXTalk_c_     = defaultRndmXTalkC;
  rndmXTalk_r_     = defaultRndmXTalkR;

  myDigis_.clear();
  eta_=-1;
  etaZs_=-1;

  landauWidth_ = defaultLandauWidth;
  binaryThreshold_ = defaultBinaryThreshold;
  sharingMethod_ = defaultSharingMethod;
  myDice.SetSeed(RANDOM_SEED);
}

Sensor::Sensor(double thickness, double pitch, double chargeDensity, double transverseDrift) {
  noise_rms_       = defaultNoiseRms;
  seed_cut_        = defaultSeedCut;
  neighbour_cut_   = defaultNeighbourCut;
  total_cut_       = defaultTotalCut;

  thickness_       = thickness;
  pitch_           = pitch;
  chargeDensity_   = chargeDensity;
  transverseDrift_ = transverseDrift;
  driftZone_       = defaultDriftZone;


  xTalk_l_         = defaultXTalkL;
  xTalk_c_         = defaultXTalkC;
  xTalk_r_         = defaultXTalkR;

  rndmXTalk_l_     = defaultRndmXTalkL;
  rndmXTalk_c_     = defaultRndmXTalkC;
  rndmXTalk_r_     = defaultRndmXTalkR;

  myDigis_.clear();
  eta_=-1;
  etaZs_=-1;

  landauWidth_ = defaultLandauWidth;
  binaryThreshold_ = defaultBinaryThreshold;
  sharingMethod_ = defaultSharingMethod;
  myDice.SetSeed(RANDOM_SEED);
}

Sensor::Sensor(double thickness, double pitch, double chargeDensity, double transverseDrift,
	       double noise_rms, double seed_cut, double neighbour_cut, double total_cut) {
  noise_rms_       = noise_rms;
  seed_cut_        = seed_cut;
  neighbour_cut_   = neighbour_cut;
  total_cut_       = total_cut;

  thickness_       = thickness;
  pitch_           = pitch;
  chargeDensity_   = chargeDensity;
  transverseDrift_ = transverseDrift;
  driftZone_       = defaultDriftZone;

  xTalk_l_         = defaultXTalkL;
  xTalk_c_         = defaultXTalkC;
  xTalk_r_         = defaultXTalkR;

  rndmXTalk_l_     = defaultRndmXTalkL;
  rndmXTalk_c_     = defaultRndmXTalkC;
  rndmXTalk_r_     = defaultRndmXTalkR;

  myDigis_.clear();
  eta_=-1;
  etaZs_=-1;

  landauWidth_ = defaultLandauWidth;
  binaryThreshold_ = defaultBinaryThreshold;
  sharingMethod_ = defaultSharingMethod;
  myDice.SetSeed(RANDOM_SEED);
}

double Sensor::getEta() {
  return eta_;
}

double Sensor::getEtaZs() {
  return etaZs_;
}

double Sensor::landau(double mpv) {
  double tempval;
  if (landauWidth_>0) {
    tempval = myDice.Landau(mpv, mpv*landauWidth_);
    return tempval;
  } else {
    return mpv;
  }
}

double Sensor::stripNoise() {
  if (noise_rms_>0) {
    return (myDice.Gaus(0.,noise_rms_));
  } else {
    return 0;
  }
}

double Sensor::wallHit(double x, double y, double theta, double width, double height, double& x1, double& y1) {
  double deltaX, deltaY;
  
  deltaX=width-x;
  deltaY=height-y;

  // std::cerr << "Track impinging in block" << std::endl;
  // std::cerr << "x: " << x << "\ty: " << y << "\tangle:" << theta/M_PI*180. << "\tdriftzone: " << driftZone_ << std::endl;
  // std::cerr << "deltaX: " << deltaX << "\tdeltaY: " << deltaY << std::endl;

  // TODO: make this nicer in case (theta \simeq \pi)

  if (sin(theta) == 1) {
    x1=x;
    y1=height;
  } else if ((deltaY*cos(theta))==(deltaX*sin(theta))) {
    x1=width;
    y1=height;
  } else if ((deltaY*cos(theta))<(deltaX*sin(theta))) {
    // It will hit the ceiling of the block
    // std::cerr << "It will hit the ceiling of the block" << std::endl;
    x1=x+deltaY*cos(theta)/sin(theta);
    y1=height;
    //} else if ((deltaY*cos(theta))==(deltaX*sin(theta))) {
    //x1=pitch_;
    //y1=thickness_;
  } else if ((deltaY*cos(theta))>(deltaX*sin(theta))) {
    // It will hit the wall of next block
    // std::cerr << "It will hit the wall of next block" << std::endl;
    x1=width;
    y1=y+deltaX*sin(theta)/cos(theta);
  } else {
    std::cerr << "ERROR: something went horribly wrong here!" << std::endl;
  }

  return sqrt( (x1-x)*(x1-x) + (y1-y)*(y1-y) );
  // debug
  // std::cout << x << " " << " " << y << " " << x1 << " " << y1 << std::endl;
}

void Sensor::trackHitElement(double& x0, double& y0, double theta, double& prevCharge, double& thisCharge, double& nextCharge) {
  double x1, y1;
  double wa=(pitch_-driftZone_);
  double L, totalL;

  prevCharge=0;
  nextCharge=0;
  totalL=0;

  if ((driftZone_>0)&&(x0<driftZone_)) { // We start by computing the left-drifting charge
    L = wallHit(x0, y0, theta, driftZone_, thickness_, x1, y1);
    totalL = L;

    // computing the sharing to the left strip
    if (sharingMethod_==blockSharing) {
      prevCharge = L/2;
    } else if (sharingMethod_==linearSharing) {
      prevCharge = L/2;
      prevCharge *= (1 - (x1+x0)/2/driftZone_);
    } else {
      std::cerr << "Unknown sharing method " << sharingMethod_;
    }

    // We move the pointer
    // to the exit point
    x0=x1;
    y0=y1;
  }

  if (y0<thickness_) { // we can go on to the next wall
    if (x0<wa) { // we have a part of the track in the proper bulk
      L = wallHit(x0, y0, theta, wa, thickness_, x1, y1);
      totalL += L;

      // We move the pointer
      // to the exit point
      x0=x1;
      y0=y1;
    }
  }

  if (y0<thickness_) { // we can go on to the next wall
    L = wallHit(x0, y0, theta, pitch_, thickness_, x1, y1);
    totalL += L;

    // computing the sharing to the left strip
    if (sharingMethod_==blockSharing) {
      nextCharge = L/2;
    } else if (sharingMethod_==linearSharing) {
      nextCharge = L/2;
      nextCharge *= (1 - (pitch_-x1+pitch_-x0)/2/driftZone_);
    } else {
      std::cerr << "Unknown sharing method " << sharingMethod_;
    }

    // We move the pointer
    // to the exit point
    x0=x1;
    y0=y1;
  }
    
  thisCharge = totalL - (prevCharge+nextCharge);
  x0-=pitch_;
}

bool Sensor::setSharingMethod(unsigned int newValue) {
  if (newValue<nSharing) {
    sharingMethod_=newValue;
    return true;
  } else {
    std::cerr << "Unknown sharing method " << newValue << std::endl;
    return false;
  }
}


void Sensor::trackHitElement_old2(double& x, double& y, double theta, double& prevCharge, double& thisCharge, double& nextCharge) {
  double x1, y1, yZ, y1Z, dummy;
  double wa=(pitch_-driftZone_);

  // std::cerr << "x1: " << x1 << "\ty1: " << y1 << std::endl;

  wallHit(x,y,theta,pitch_,thickness_,x1,y1);
  
  prevCharge = 0;
  nextCharge = 0;
  if (driftZone_>0) {

    if (x<driftZone_) {
      wallHit(x,y,theta,driftZone_,thickness_,dummy,y1Z);
      if ((dummy != driftZone_)&&(dummy!=x1)) {
	std::cerr << "Error: ((dummy != driftZone_)&&(dummy!=x1))" << std::endl;
	std::cerr << "dummy: " << dummy << "\tdriftZone_: " << driftZone_ << "\tx1: " << x1 << std::endl;
      }
    }

    if (x1>wa) {
      if (x<wa) {
	wallHit(x,y,theta,wa,thickness_,dummy,yZ);
	if (dummy != wa) {
	  std::cerr << "Error: (dummy != wa)" << std::endl;
	  std::cerr << "dummy: " << dummy << "\twa: " << wa << std::endl;
	}
      }
    }

    // WARNING important
    // TODO: change the code in case cos(theta)<<1
    if (x<driftZone_) {
      // We share some charge with previous strip
      // std::cerr << "We share some charge with previous strip" << std::endl;
      if (x1<driftZone_) {
	// x0<a
	// x1<a
	prevCharge = (y1-y)*(2*driftZone_-x-x1)/(4.*driftZone_*sin(theta));
      } else {
	// x0<a
	// x1>a
	prevCharge = (y1Z-y)*(driftZone_-x)/(4*driftZone_*sin(theta));
      }
    }

    if (x1>wa) {
      // std::cerr << "We share some charge with next strip" << std::endl;
      // We share some charge with next strip
      if (x>wa) {
	// x1>w-a
	// x0>w-a
	nextCharge = (y1-y)*(x1+x-2*wa) / (4. * driftZone_ * sin(theta));
      } else {
	// x1>w-a
	// x0<w-a
	nextCharge = (x1-wa) * (y1-yZ)  / (4. * driftZone_ * sin(theta));
      }
    }
  }

  // WARNING important
  // TODO: change the code in case cos(theta)<<1 
  thisCharge = sqrt(pow(y1-y,2)+pow(x1-x,2)); // ((x1-x)/cos(theta));
  // std::cerr << "thisCharge: " << thisCharge << std::endl;
  thisCharge -= (prevCharge+nextCharge);
  // std::cerr << "prevCharge: " << prevCharge << "\tthisCharge: " << thisCharge << "\tnextCharge: " << nextCharge << std::endl;


  // Last: we move the pointer
  // to the exit point
  x=x1;
  y=y1;

  // std::cerr << "Track exiting from block at" << std::endl;
  // std::cerr << "x: " << x << "\ty: " << y << std::endl;

  x-=pitch_;
}



int Sensor::trackHit(double impact, double theta) {
  myDigis_.clear();  
  eta_=-1;
  etaZs_=-1;

  if ((theta<0)||(theta>M_PI)) return 0;
  if (theta>(M_PI/2.)) theta=M_PI-theta;

  digiV tracciaClean;
  tracciaClean.push_back(0);

  double prevCharge, thisCharge, nextCharge;
  double formerNextCharge;
  double x=impact;
  double y=0;
  int lastIndex;
  bool last=false;

  while (x>pitch_) {
    tracciaClean.push_back(0);
    x -= pitch_;
  }

  // std::cerr << std::endl << std::endl;
  // std::cerr << "***************************" << std::endl;
  // std::cerr << "******   NEW TRACK   ******" << std::endl;
  // std::cerr << "***************************" << std::endl;
  // std::cerr << std::endl;
  nextCharge=0;
  while (!last) {
    formerNextCharge=nextCharge;
    trackHitElement(x, y, theta, prevCharge, thisCharge, nextCharge);
    // std::cerr << prevCharge << "\t" << thisCharge << "\t" << nextCharge << std::endl;
    if ((x==pitch_)||(y>=thickness_)) {
      // std::cerr << "quitting loop!"<< std::endl;
      last=true;
    }
    lastIndex=tracciaClean.size()-1;
    tracciaClean[lastIndex]+=prevCharge;
    tracciaClean.push_back(thisCharge+formerNextCharge);
  }
  tracciaClean.push_back(nextCharge);
  // std::cerr << "create loop done" << std::endl;

  // std::cerr << std::endl << "Track lenghts summary:" << std::endl;
  for(int i=0; i<tracciaClean.size(); i++) {
    // std::cerr << tracciaClean[i] << "\t";
    myDigis_.push_back(landau(tracciaClean[i]*chargeDensity_) + stripNoise());
    // Just in case you wanted no noise use the following line
    // instead of previous
    // myDigis_.push_back(tracciaClean[i]*chargeDensity_);
  }
  // std::cerr << std::endl;
  // std::cerr << "******  END OF TRACK  *****" << std::endl;
  // std::cerr << "***************************" << std::endl;
  // std::cerr << std::endl << std::endl;

  return myDigis_.size();
}


int Sensor::trackHitElement_old(double impact, double theta, digiV &result) {

  result.clear();
  if ((theta<0)||(theta>M_PI)) return 0;
  if (theta>(M_PI/2.)) theta=M_PI-theta;

  int nHit = 0;
  double signal;
  double x = impact;
  double y = 0;
  double yStep;
  double yStepPitch = pitch_ * tan(theta);
  bool last = false;

  //double totalSignal=0; // TODO: remove this

  // std::cout << std::endl;
  // std::cout << "*******NEW STRIP***********" << std::endl;

  // std::cout << "Angle: " << theta * 180 / M_PI << std::endl;
  // std::cout << "yStepPitch: " << yStepPitch << std::endl;

  while (x>pitch_) {
    result.push_back(stripNoise());
    x -= pitch_;
  }

  while (!last) {
    // std::cout << "Starting strip " << nHit << std::endl;
    // std::cout << "x: " << x << "; y: " << y << std::endl;
    if (x==0) {
      if ((y + yStepPitch) > thickness_) {
	yStep = thickness_ - y;
	last = true;
      } else {
	yStep = yStepPitch;
      }
    } else {
      yStep = (pitch_-x)/cos(theta)*sin(theta);
      // std::cout << "yStep: " << yStep << std::endl;
      if (yStep<thickness_) {
	x = 0;
      } else {
	yStep = thickness_;
	last = true;
      }
    }
    y += yStep;

    signal = yStep / sin(theta) *  chargeDensity_;
    // std::cout << "Signal[" << nHit << "] = " << signal << "\t";
    signal = landau(signal);
    signal = round (signal + stripNoise());
    // std::cout << std::endl << "###"<<signal << std::endl;
    // std::cout << "with noise = " << signal << std::endl;
    result.push_back(signal);
    // totalSignal += signal;
    nHit++;
  }

  // std::cout << std::endl << "###"<<totalSignal << std::endl;
  // std::cout << "Hit strips: " << nHit << std::endl;
  return nHit;
}

int Sensor::trackHit_old(double impact, double theta) {
  myDigis_.clear();  
  eta_=-1;
  etaZs_=-1;

  if ((theta<0)||(theta>M_PI)) return 0;
  if (theta>(M_PI/2.)) theta=M_PI-theta;

  int maxIndex;
  digiV tracce[N_TRACCE];
  // TODO: maybe here we can add weights to the tracks
  // simualting the drift

  //  double pesi[N_TRACCE];
  double temp;
  double shift;
  double xDrift = transverseDrift_ / sin(theta);
  if (xDrift > (pitch_/2)) std::cerr << "Alarm! Is this track too inclined?" << std::endl;

  for (int i=-DRIFT_WING; i <= DRIFT_WING ; i++ ) {
    //    pesi[i+DRIFT_WING]= [...]
    shift = (xDrift*i/double(DRIFT_WING));
    trackHitElement_old(impact+pitch_+shift, theta, tracce[i+DRIFT_WING]);
    tracce[i+DRIFT_WING].push_back(stripNoise());
  }

  maxIndex = 0;
  for (int i=0; i<N_TRACCE; i++) {
    maxIndex = maxValue(maxIndex, tracce[i].size());
  }

  for (int j=0; j< N_TRACCE; j++) {
    for (int i=tracce[j].size(); i<maxIndex; i++)
      tracce[j].push_back(stripNoise());
  }
  for (int i=myDigis_.size(); i<maxIndex; i++)
    myDigis_.push_back(stripNoise());

  for (int i=0; i<maxIndex; i++) {
    temp = 0;
    for (int j=0; j< N_TRACCE; j++) {
      //temp +=  pow(tracce[j][i],2); // alternative quadratic sum
      temp += tracce[j][i];
    }
    temp /= N_TRACCE;
    //myDigis_[i]=sqrt(temp); // alternative quadratic sum
    myDigis_[i]=temp;
  }
}

void Sensor::crossTalk() {
  digiV newDigi;
  double prev, foll;
  double left_talk;
  double right_talk;
  double unmoved;;
  for (int i=0; i< myDigis_.size(); i++ ) {
    left_talk = xTalk_l_ + rndmXTalk_l_ * myDice.Rndm();
    right_talk = xTalk_r_ + rndmXTalk_r_ * myDice.Rndm();
    unmoved = 1 - left_talk - right_talk;
    if ((i-1)<0) {
      prev=stripNoise();
    } else {
      prev=myDigis_[i-1];
    }
    if ((i+1)>=myDigis_.size()) {
      foll=stripNoise(); 
    } else {
      foll=myDigis_[i+1];
    }
    newDigi.push_back(prev*left_talk+myDigis_[i]*unmoved+foll*right_talk);
  }
  myDigis_.clear();
  for (int i=0;i<newDigi.size(); i++) myDigis_.push_back(newDigi[i]);
}

bool Sensor::binaryCluster(int& nStrip, int& clusterWidth, const bool& onlyLeading, double& baricenter, double& baricenter_ext, double& totalRealCharge) {
  double usedNoise;
  digiV adcDigi;
  baricenter=0;
  baricenter_ext=0;
  totalRealCharge=0;

  if (noise_rms_>0) {
    usedNoise=noise_rms_;
  } else {
    usedNoise=defaultNoiseRms;
  }

  int firstStrip = -1;
  int lastStrip = -2;
  if (onlyLeading) {
    //    int iMax = -1;
    //    double prevVal=0; 
    //    for (int i=0; i< myDigis_.size(); i++ ) {
    //      if (myDigis_[i]>binaryThreshold_) {
    //        if (myDigis_[i]>prevVal) iMax=i;
    //        prevVal = myDigis_[i];
    //      } else {
    //        if (iMax!=-1) break;
    //      }
    double maxVal=0;
    for (int i=0; i< myDigis_.size(); i++ ) {
      // std::cout << "["<<i<<"]="<<myDigis_[i]<< " "; // debug
      if (myDigis_[i]>binaryThreshold_) {
        if (myDigis_[i]>maxVal) {
          maxVal=myDigis_[i];
          firstStrip=i;
        }
      }
    }
    // std::cout << " lead = " <<firstStrip << std::endl; // debug
    lastStrip = firstStrip;
    baricenter = firstStrip;
    baricenter_ext = firstStrip;
    totalRealCharge = maxVal;
  } else { 
    double totalCharge=0;
    double totalCharge_ext=0;
    double signalPrecision=1<<(adcBits_);
    //std::cout << "Signal precision is " << signalPrecision << " ";
    for (int i=0; i< myDigis_.size(); i++ ) {
      //std::cout << "myDigis_["<<i<<"]="<<myDigis_[i]<< " "; // debug
      if (myDigis_[i]>binaryThreshold_) {
        if (firstStrip==-1) {
          if (i!=0) {
            double aValue = myDigis_[i-1];
            if (aValue>saturation_) aValue=.99999999999; else aValue /= saturation_;
            if (aValue<0) aValue=0;
            aValue=(ceil(aValue*signalPrecision)-0.5)/signalPrecision; // n-bit resolution over the saturation value
            baricenter_ext+=(i-1)*aValue;
            totalCharge_ext+=aValue;
          }
          firstStrip=i;
        }
        lastStrip=i;
        double aValue = myDigis_[i];
        //std::cout << "a value " << aValue << " becomes ";
        if (aValue>saturation_) aValue=.99999999999; else aValue /= saturation_;
        //std::cout << aValue;
        aValue=(ceil(aValue*signalPrecision)-0.5)/signalPrecision; // n-bit resolution over the saturation value
        //std::cout << aValue*saturation_ << std::endl;
        baricenter+=i*aValue;
        totalCharge+=aValue;
        totalRealCharge+=myDigis_[i];
      }
      if ((myDigis_[i]<binaryThreshold_)&&(firstStrip!=-1)) {
        if (i+1<myDigis_.size()) {
          double aValue = myDigis_[i+1];
          if (aValue>saturation_) aValue=.99999999999; else aValue /= saturation_;
          if (aValue<0) aValue=0;
          aValue=(ceil(aValue*signalPrecision)-0.5)/signalPrecision; // n-bit resolution over the saturation value
          baricenter_ext+=(i+1)*aValue;
          totalCharge_ext+=aValue;
        }
        break;
      }
    }
    baricenter_ext+=baricenter;
    totalCharge_ext+=totalCharge;
    baricenter/=totalCharge;
    baricenter_ext/=totalCharge_ext;
  }
  //std::cout << std::endl; //debug
  
  nStrip  = firstStrip;
  clusterWidth = lastStrip-firstStrip+1;

  return (clusterWidth>0);
}
 
void Sensor::clusterize(int &nClus, int &totCharge) {
  int seedIndex=-1;
  int minClusterIndex=-1;
  int maxClusterIndex=-1;
  int clusterCount=0;
  digiV adcDigi;
  
  double usedNoise;

  if (noise_rms_>0) {
    usedNoise=noise_rms_;
  } else {
    usedNoise=defaultNoiseRms;
  }

  etaZs_=-1;
  eta_=-1;

#ifdef _CLUSTERIZE_DEBUG_
  std::cerr << std::endl ; 
  std::cerr << "*****************************" <<std::endl;
  std::cerr << "*****START OF NEW DIGI*******" <<std::endl;
  std::cerr << "*****************************" <<std::endl;
  std::cerr <<std::endl;
#endif

  for (int i=0; i< myDigis_.size(); i++ ) {
#ifdef _CLUSTERIZE_DEBUG_
    std::cerr << "digi[" << i << "]=" << myDigis_[i] << "\t";
#endif
    adcDigi.push_back(round(myDigis_[i]));
  }
#ifdef _CLUSTERIZE_DEBUG_
  std::cerr << std::endl;
  std::cerr << "***SEED FIND*****" << std::endl;
#endif
  for (int i = 0; i < adcDigi.size(); i++) {
#ifdef _CLUSTERIZE_DEBUG_
    std::cerr << "digi[" << i << "]=" << adcDigi[i] << "\t" << "Soglia: " << usedNoise * seed_cut_ ;
#endif
    if (adcDigi[i] > (usedNoise * seed_cut_)) {
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << " passed" << std::endl;
#endif
      seedIndex=i;
      clusterCount=1;
      minClusterIndex=i;
      maxClusterIndex=i;
      break;
    } else {
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << std::endl;
#endif
    }
  }
  if (seedIndex != -1) {
#ifdef _CLUSTERIZE_DEBUG_
    std::cerr << "***NEIGHBOUR FIND*****" << std::endl;
#endif
    for (int i=seedIndex-1; i>=0; i--) {
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "digi[" << i << "]=" << adcDigi[i] << "\tSoglia: " << (usedNoise * neighbour_cut_);
#endif
      if ((adcDigi[i])>=(usedNoise * neighbour_cut_)) {
#ifdef _CLUSTERIZE_DEBUG_
	std::cerr << " passed!" << std::endl;
	
#endif
	clusterCount++;
	minClusterIndex=i;
      } else {
#ifdef _CLUSTERIZE_DEBUG_
	std::cerr << std::endl;
#endif	
	break;  
      }
    }
    for (int i=seedIndex+1; i<adcDigi.size(); i++) {
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "digi[" << i << "]=" << adcDigi[i] << "\tSoglia: " << (usedNoise * neighbour_cut_);
#endif
      if ((adcDigi[i])>=(usedNoise * neighbour_cut_)) {
#ifdef _CLUSTERIZE_DEBUG_
	std::cerr << " passed!" << std::endl;
#endif
	clusterCount++; 
	maxClusterIndex=i;
      } else { 
#ifdef _CLUSTERIZE_DEBUG_
	std::cerr << std::endl;
#endif	
	break;  
      }
    }
  }
#ifdef _CLUSTERIZE_DEBUG_
  std::cerr << std::endl << "Now computing the eta, charge and baricenter" << std::endl;
#endif
  double baricenter;
  double totalCharge;
  int sx, dx;
  int iSx, iDx;

  // (((((((((((((((((((
  // )))))))))))))))))))
  // (((((((((((((((((((
  // )))))))))))))))))))
  // (((((((((((((((((((
  // )))))))))))))))))))

  // We compute the eta without the zero suppression
  // first. The baricenter is found using only the
  // cluster, while the eta value is computed using
  // the raw data
  baricenter=0;
  totalCharge=0;
  if (clusterCount!=0) {
    for (int i = minClusterIndex; i <=maxClusterIndex; i++) {
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "digi[" << i << "]=" << adcDigi[i] << "\t";
#endif
      baricenter+=i*adcDigi[i];
      totalCharge+=adcDigi[i];
    }
  }
  baricenter/=totalCharge;
#ifdef _CLUSTERIZE_DEBUG_
  std::cerr << std::endl;
  std::cerr << "  baricenter  = " << baricenter  << std::endl;
  std::cerr << "  totalcharge = " << totalCharge << std::endl;
#endif

  // The strips for the eta function are the ones
  // just right and left with respect to the cluster bariccenter
  iSx = int(floor(baricenter));
  iDx = int(ceil(baricenter));
  
  // If the baricenter of the cluser layes exactly underneath a strip
  // then we shall see which strip (right or left to the baricenter)
  // has the highest signal. We'll do this by moving the iDx and iSx indexes
  if (iDx==iSx) {
    dx = int(adcDigi[iDx+1]);
    sx = int(adcDigi[iDx-1]);
#ifdef _CLUSTERIZE_DEBUG_
    std::cerr << "dx: " << dx << "\tsx: " << sx << std::endl;
#endif
    if (dx==sx) {
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "Symmetric cluster:";
#endif
      // Symmetric cluster
      // we choose random wherther to pick the
      // right or left one
      if (rand()<double(RAND_MAX)/2.) {
	iDx=iDx+1;
      } else {
	iSx=iSx-1;
      }
    } else {
      // Non-symmetric cluster: we pick
      // the highest signal strip
      if (dx>sx) {
#ifdef _CLUSTERIZE_DEBUG_
    std::cerr << "dx>sx" << std::endl;
#endif  
	iDx=iDx+1;
      } else {
	iSx=iSx-1;
      }
    }
  }

  // And now we can compute the eta
  sx = int(adcDigi[iSx]);
  dx = int(adcDigi[iDx]);
  eta_ = double(sx) / (double(dx) + double(sx));
#ifdef _CLUSTERIZE_DEBUG_
  std::cerr << "Sx: " << sx << "\tDx: " << dx << "\tEta: " << eta_ << std::endl;
#endif

  baricenter=0;
  totalCharge=0;

  // We compute the eta with the zero suppression.
  // The baricenter is found using only the
  // cluster. If we get a 1-cluster signal, then we will
  // choose randomly between 0 and 1
  if (clusterCount!=0) {
    if (clusterCount==1) {
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "Single strip" << std::endl;
#endif
      if (rand()<double(RAND_MAX)/2.) {
	etaZs_ = 0;
      } else {
	etaZs_ = 1;
      }
      totalCharge=adcDigi[seedIndex];
      baricenter=seedIndex;
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "  baricenter  = " << baricenter  << std::endl;
      std::cerr << "  totalcharge = " << totalCharge << std::endl;
      std::cerr << "  Eta:Zs      = " << etaZs_ << std::endl;
#endif
    } else {
      for (int i = minClusterIndex; i <=maxClusterIndex; i++) {
#ifdef _CLUSTERIZE_DEBUG_
	std::cerr << "digi[" << i << "]=" << adcDigi[i] << "\t";
#endif
	baricenter+=i*adcDigi[i];
	totalCharge+=adcDigi[i];
      }
      baricenter/=totalCharge;
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "  baricenter  = " << baricenter  << std::endl;
      std::cerr << "  totalcharge = " << totalCharge << std::endl;
#endif
      iSx = int(floor(baricenter));
      iDx = int(ceil(baricenter));
      if (iDx==iSx) {
	dx = int(adcDigi[iDx+1]);
	sx = int(adcDigi[iDx-1]);
	if (dx==sx) {
	  // Symmetric cluster
	  // we choose random whether to pick the
	  // right or left one
	  if (rand()<double(RAND_MAX)/2.) {
	    iDx=iDx+1;
	  } else {
	    iSx=iSx-1;
	  }
	} else {
	  // Non-symmetric cluster: we pick
	  // the highest signal strip
	  if (dx>sx) {
	    iDx=iDx+1;
	  } else {
	    iSx=iSx-1;
	  }
	}
      } 
      sx = int(adcDigi[iSx]);
      dx = int(adcDigi[iDx]);
      etaZs_ = double(sx) / (double(dx) + double(sx));
#ifdef _CLUSTERIZE_DEBUG_
      std::cerr << "Sx: " << sx << "\tDx: " << dx << "\tEta:Zs " << etaZs_ << std::endl;
#endif
    }
    //  etaFunction_->Fill(baricenter/4);
  }
  nClus = clusterCount;
  totCharge = int(totalCharge);
#ifdef _CLUSTERIZE_DEBUG_
  std::cerr << "nCluster  = " << nClus     << std::endl;
  std::cerr << "totCharge = " << totCharge << std::endl;
  std::cerr << std::endl;
#endif
  
}

void Sensor::addNeighbourNoise(int nStrips) {
  digiV newDigi;
  for (int i=0; i<nStrips; ++i) newDigi.push_back(stripNoise());
  for (int i=0; i< myDigis_.size(); i++ ) newDigi.push_back(myDigis_[i]);
  for (int i=0; i<nStrips; ++i) newDigi.push_back(stripNoise());
  myDigis_.clear();
  for (int i=0;i<newDigi.size(); i++) myDigis_.push_back(newDigi[i]);
}
