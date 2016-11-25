

#include "../include/type.hpp"
#include <PDE.hpp>
#include "Timer.h"


int main(int argc, char** argv) 
{
  siwir::Timer time;
  PDE p(2048,2048);
  p.applyBoundary();
  //p.writeFile("u.out",p.u);
  
  time.reset();
  for(uint i=0;i<500;++i)
  {
  p.RedBlackGaussSeidal();
  }
  real runtime = time.elapsed();
  std::cout<< "Runtime for 500 iterations: "<<runtime<<std::endl;
  p.writeFile("update.txt",p.u);
  p.Residual();
  
  //p.print(p.force);
  
  return 0;
}
