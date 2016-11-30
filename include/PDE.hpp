
#pragma once

#include "../include/type.hpp"

class PDE
{
 public:
   
   PDE(){;}
   ~PDE(){
     delete[] u;
     delete[] residual;
     delete[] force;
  }
   

   PDE( const uint&,const uint& ,const double& hinX, const double& hinY);
   void applyBoundary(void);
   void RedBlackGaussSeidal(const uint& iteration);
   void gaussJacobi(real* src, real* dst, const uint& iteration);
   bool writeFile(const std::string& fileName,const real* vec);
   
   void operator()(real *__restrict src,real *__restrict dst);
   real ResidualNorm(void);

   void clearData(real *__restrict src);
   
  real *__restrict u;	// displacement vector
  real *__restrict residual;	// residual vector   
  real *__restrict force;	// force vector
	
  uint nx , ny;		// number of grid in x(hor.) and y(ver.) direction
  real hx , hy;		//mesh size in x and y direction
     
};	//class
  
