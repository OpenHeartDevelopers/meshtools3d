#include"../include/Chrono.hpp"
#include <iostream>
#include <time.h>
#include<string>

std::ostream& operator << (std::ostream& c,const Chrono& cht)
{

  size_t hours=0, minutes=0,seconds=0;
  if(cht.counting)
  {
     hours=9999999999, minutes=59,seconds=59;
  }
  else
  {
    time_t t_diff=cht.t_elapsed; // I copy the value if I want to look at partial times
    if(t_diff>=3600)
    {
      hours=t_diff/3600;
      t_diff=t_diff-hours*3600;
    }

    if(t_diff>=60)
    {
      minutes=t_diff/60;
      t_diff=t_diff-minutes*60;
    }
    seconds=t_diff;

  }
  c<<hours<<":"<<std::setw(2)<<std::right<<std::setfill('0')<<minutes<<":"<<std::setw(2)<<std::right<<std::setfill('0')<<seconds;
  return c;
};


Chrono::Chrono()
:counting(false),
t_start(0),
t_end(0),
t_elapsed(0)
{}

Chrono::~Chrono()
{
  reset();
}

void Chrono::start()
{
  if(~counting)
  {
    time(&t_start);
    t_end=t_start;
    counting=true;
  }
  else
  {
    std::cout<<"CHRONO: WARNING: NO START SINCE COUNTER WAS PREVIOUSLY STARTED"<<std::endl;
  }

}
void Chrono::stop()
{
  if(counting)
  {
    time(&t_end);
    t_elapsed=t_elapsed+(t_end-t_start);
    counting=false;
  }

}

void Chrono::reset()
{
  this->stop();
  counting=false;
  t_start=0;
  t_end=0;
  t_elapsed=0;
}

void Chrono::printElapsedTime()
{
  if(counting)
  {
    std::cout<<"CHRONO: WARNING: YOU HAVE TO STOP CHRONOMETER BEFORE EVALUATING ELAPSED TIME"<<std::endl;
  }
  else
  {
    time_t t_diff=t_elapsed; // I copy the value if I want to look at partial times
    size_t hours=0, minutes=0,seconds=0;
    if(t_diff>=3600)
    {
      hours=t_diff/3600;
      t_diff=t_diff-hours*3600;
    }

    if(t_diff>=60)
    {
      minutes=t_diff/60;
      t_diff=t_diff-minutes*60;
    }
    seconds=t_diff;
    std::cout<<hours<<":"<<std::setw(2)<<std::right<<std::setfill('0')<<minutes<<":"<<std::setw(2)<<std::right<<std::setfill('0')<<seconds;

  }
}
