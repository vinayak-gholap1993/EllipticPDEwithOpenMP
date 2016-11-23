

#include "../include/type.hpp"
#include <PDE.hpp>


int main(int argc, char** argv) 
{
  PDE p(1024,1024);
  p.applyBoundary();
  p.writeFile("u.out",p.u);
  
  for(uint i=0;i<500;++i)
  {
  p.RedBlackGaussSeidal();
  }
  p.writeFile("update.out",p.u);
  p.Residual();
  
  //p.print(p.force);
  
  return 0;
}