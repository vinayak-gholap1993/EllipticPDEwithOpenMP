

#pragma once

#include "type.hpp"

class PDE
{
 public:
   
   PDE(){;}
   ~PDE(){;}
   
   PDE( const uint&, const uint& );
   void applyBoundary(void);
   void print(const std::vector<real> &obj);
   void trigTable(void);
   void RedBlackGaussSeidal(void);
   void Residual(void);
   bool writeFile(const std::string& fileName,const std::vector<real>& vec);

  std::vector<double> u;	// displacement vector
  std::vector<double> residual;	// residual vector   
  std::vector<double> force;	// force vector
	
  uint nx , ny;		// number of grid in x(hor.) and y(ver.) direction
  real hx , hy;		//mesh size in x and y direction
     
};	//class
  
