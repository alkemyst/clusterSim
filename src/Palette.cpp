#include <Palette.h>

Color_t Palette::Color(const unsigned int& plotIndex) {
  std::string colorCode;
  
  if (plotIndex==0) colorCode = "#000000";
  else {
    int nColor=(plotIndex-1) % 12;
    switch (nColor) {
    case 0 :
      colorCode="#004586";
      break;
    case 1 :
      colorCode="#FF420E";
      break;
    case 2 :
      colorCode="#FFD320";
      break;
    case 3 :
      colorCode="#579D1C";
      break;
    case 4 :
      colorCode="#7E0021";
      break;
    case 5 :
      colorCode="#83CAFF";
      break;
    case 6 :
      colorCode="#314004";
      break;
    case 7 :
      colorCode="#AECF00";
      break;
    case 8 :
      colorCode="#4B1F6F";
      break;
    case 9 :
      colorCode="#FF950E";
      break;
    case 10 :
      colorCode="#C5000B";
      break;
    case 11 :
      colorCode="#0084D1";
      break;
    default :
      std::cerr << "ERROR: in Vizard::getNiceColor() n%12 is not an int between 0 and 11! This should not happen." << std::endl;
      colorCode="#000000";
      break;
    }
  }
  
  return TColor::GetColor(colorCode.c_str());
}

