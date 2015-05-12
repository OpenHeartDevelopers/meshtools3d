#ifndef _CHRONO_HPP
#define _CHRONO_HPP

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include<iostream>
#include <time.h>

class Chrono
{
  public:
    Chrono();
    ~Chrono();
    void start();
    void stop();
    void reset();
    void printElapsedTime();
    friend std::ostream& operator << (std::ostream& c,const Chrono& cht);
  private:
  bool counting;
  time_t t_start;
  time_t t_end;
  time_t t_elapsed;
  
  
  
};

#endif
