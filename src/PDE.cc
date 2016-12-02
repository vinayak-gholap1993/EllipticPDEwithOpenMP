

#include "../include/PDE.hpp"
#include<omp.h>

#include <fstream>


PDE::PDE( const uint& numx , const uint& numy , const double& hinX ,const double& hinY)
{ 
   this->hx = hinX;
   this->hy = hinY;
  
   uint totalGrid; 
   
   if(numx < 300 || numy < 300)
   {
     totalGrid = (numx + 65) * (numy+1); 
   
   this->nx = (numx  + 1 );	//# of grids in x-direction
   this->ny = (numy + 1);	//# of grids in y-direction
  }
  
   else
   {
    totalGrid = ((numx>>2) + 65) * (numy+1); 
   
   this->nx = ((numx >> 2) + 1);	//# of grids in x-direction (Domain partitioned in X axis in 4 parts)
   this->ny = (numy + 1);
    }

#pragma omp parallel
{
   this->u = new double __attribute__((aligned(64))) [totalGrid];
   this->residual = new double __attribute__((aligned(64))) [totalGrid];
   this->force = new double __attribute__((aligned(64))) [totalGrid];
   
#pragma omp for schedule(static)
   for(uint row = 0;row < totalGrid; ++row)
   {
     this->u[row] = 0.0;
     this->force[row] = 0.0;   
   } //for row
   
  
}
   //std::cout<<BOLD(FBLU(" constructor 2 done "))<<std::endl;
}

void PDE::applyBoundary(void)
{
 //------------------------------- u init ----------------------------------------------------------------------------
  real temp = std::sinh(twoPi);
  for(uint row = this->ny-1 , col = 0 ; col < this->nx ; col++ )
  {
    this->u[row * this->nx + col] = std::sin(twoPi * this->hx * col) * temp;
  }
  
//-------------------------------- f init -----------------------------------------------------------------------------
  for(uint row = 0; row < this->ny ; ++row)
    for(uint col = 0; col < this->nx ; ++col)
    {  
      this->force[ row * this->nx + col ] = twoPi * twoPi * sin(twoPi * this->hx * col) * sinh(twoPi * this->hy * row); 
    }
    //std::cout<<BOLD(FBLU(" applied boundary condition done "))<<std::endl;
}


//---------------------------------- RBGS ---------------------------------------
void PDE::RedBlackGaussSeidal(const uint& iteration)
{
    //uint numberofgridpoints = nx * ny;
real hxSquare = 1.0 /(hx * hx), hySquare= 1.0 /(hy * hy);		/// Inverse of hxSquare and hySquare
const real constant = 1.0/ ((2.0 * hxSquare) +(2.0 * hySquare)+ 4.0 * pi * pi);
uint rowEnd = ny-1 , columnEnd = (nx-1) ,column;

////------------------------------------------------RED UPDATE------------------------------------------------------ 

#pragma omp parallel
{
for(uint iter = 0;iter < iteration;++iter)
{
  for(uint skip = 0;skip < 2; ++skip)
  {  
 #pragma omp for schedule(static) private(column)
  for (uint row= 1; row < rowEnd; ++row)
    {
   for (column =1 + ((row+skip) & 1) ; column < columnEnd; column+=2)
   {
     u[ row * nx + column] =  constant * (((u[row * nx +(column-1)] +  u[ row * nx +(column+1)]) * hxSquare) + ((u[(row - 1) * nx +column] + u[(row + 1) * nx +column])* hySquare) + force[row * nx +column]);  
  } //column
    } //row
   } //skip
  
 
} //iter
}
//std::cout<<BOLD(FBLU(" RBGS done "))<<std::endl;
}


//------------------------------- residual -------------------------------
real PDE::ResidualNorm(void)
{
  real hxSquare = hx * hx, hySquare= hy * hy;
  real hxInv = 1.0 / hxSquare , hyInv = 1.0 / hySquare;
  const real constTerm = ((2.0 *hxInv ) +(2.0 * hyInv)+ 4 * pi * pi);
  real temp = 0.0, norm = 0.0, tgInv = 1.0 / ((nx-1)*(ny-1));

    for(uint row = 1; row < ny-1; ++row)
    {
        for(uint column = 1; column < nx-1 ; ++column)
        {
            temp = force[row * nx + column] + (u[row* nx +(column-1)]  + u[row* nx +(column+1)]) * hxInv +
                                                                     (u[(row-1) * nx +column] + u[(row+1) * nx +column]) * hyInv -
                                                                     constTerm * u[row * nx + column];
            norm += temp * temp;
        }
    }
    const real normValue = sqrt(norm * tgInv);
    
    return normValue;
}

//---------------------------------- File write ----------------------------------
bool PDE::writeFile(const std::string& fileName,const real* vec)
{
  std::ofstream file(fileName);		//object of ofstream
  //file.open(fileName,std::fstream::in | std::fstream::out | std::fstream::app);
  real signBit = 1.0;
  if(file.is_open())
  {   
     
  uint newNX = nx-1 , newNY = ny-1;
 // std::cerr << BOLD(FRED("nx "))<<nx-1 <<BOLD(FRED("\t ny "))<<newNY<<std::endl;
  file << "#  x  y  u(x,y) \n";
  
  if(nx < 100 || ny < 100)
  {
     for(uint row = 0; row <  ny-1 ; ++row)
        for(uint col = 0; col < nx-1 ; ++col)
	  {
	    file << col*hx <<" "<<row*hy <<" "<< vec[row*nx+col] <<"\n";
	 }
	 file<<"\n";
  }
  
  else
  {
    
	for(uint row = 0; row <= newNY ; ++row)
	{
	  for(uint loop = 0;loop<=3 ; ++loop)
	    {
	      if (loop%2 == 1)	{signBit = -1.0;}	//odd
	      else 	{signBit = 1.0;}
	
	      real temps = loop * newNX  ,colS = temps ;
	      //std::cout<<temps<<std::endl;
	  
	    for(uint col = 0 ; col < newNX ; ++col,colS++)
	    {
	      file << colS*hx <<" "<<row*hy <<" "<< ((signBit) * vec[row*nx+col]) <<"\n";
	      if(loop == 3 && col== newNX-1)		{ file << 2.0 <<" "<<row*hy <<" "<<  0 <<"\n";	}	
	    }
	     temps = 0.0;
	    } //loop
	   file<<"\n ";
	}
  }

        
   }//if
     
      //std::cout << BOLD(FBLU(" file write done  ")) <<std::endl;
    return true;
 /*
  else
  {
    std::cerr << BOLD(FRED("\n Unable to open file  ")) <<std::endl;
    return false;
  }
  ---*/
  file.close(); 
}

