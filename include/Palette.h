#ifndef palette_h
#define palette_h

#include <TColor.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <map>

class Palette {
 public:
  static Color_t Color(const unsigned int& plotIndex);
};

#endif
