

#include "../include/PDE.hpp"
#include <vector>



PDE::PDE( const uint& numx , const uint& numy )
{ 
  hx = 2.0 / numx;
  hy = 1.0 / numy;
  
  uint totalGridPoint = (numx+1) * (numy+1);
  nx = numx+1; 
  ny = numy+1;
  
    u.resize(totalGridPoint,0.0);
    residual.resize(totalGridPoint,0.0);
    force.resize(totalGridPoint,0.0);
    
}

// hor(col) = j = ny = hy
// ver(row) = i = nx = hx
void PDE::applyBoundary(void)
{
 //------------------------------- u init ---------------------
  real temp = sinh(twoPi);
  for(uint row = ny-1 , col = 0 ; col < nx ; col++ )
  {
    u[row * nx + col] = sin(twoPi * hx*col) * temp;
  }
//-------------------------------- f init ---------------------
  for(uint row = 0; row < ny ; ++row)
    for(uint col = 0; col < nx ; ++col)
    {  
      force[ row * nx + col ] = twoPi * twoPi * sin(twoPi * hx * col) * sinh(twoPi * hy * row); 
    }
}

void PDE::print(const std::vector<real> &obj)
{
  /*
  for(auto iter = obj.begin() ; iter < obj.end() ; ++iter)
    std::cout << *iter<< "\t";*/
  
  for(int row = ny-1 ; row >=0 ; row--)
  {
    for(uint col = 0 ; col < nx ; col++)
    {
      std::cout<<obj[ row * nx + col ]<<"\t";
    }
  std::cout<<"\n";
  }
}

void PDE::trigTable(void)
{
  real table[nx*ny]; 
  for(uint row = 0; row < ny ; ++row)
    for(uint col = 0; col < nx ; ++col)
    {
	;
    }
}


void PDE::RedBlackGaussSeidal(void )
{
    //uint numberofgridpoints = nx * ny;
real hxSquare = hx * hx, hySquare= hy * hy;
////------------------------------------------------RED UPDATE-----------------------------------------------------------------

    for (uint row= 1; row< ny-1; ++row)
    {
        if(row & 1)
        {
            for (uint column =1; column < nx-1; column+=2)
            {
               u[row * nx + column] = 1.0/ ((2.0 / hxSquare) +(2.0 /hySquare)+ 4 * pi * pi) * ((u[row* nx +(column-1)]/ hySquare) +
                                                                                        (u[row* nx +(column+1)]/ hySquare) +
                                                                                        (u[(row-1) * nx +column]/ hxSquare) +
                                                                                        (u[(row+1) * nx +column]/ hxSquare) +
                                                                                        force[row* nx +column]);
            }
        }
    else
        {
            for (uint column =2; column < nx-1; column+=2)
            {
               u[row* nx + column] = 1.0/ ((2.0 / hxSquare) +(2.0 /hySquare)+ 4 * pi * pi) * ((u[row* nx +(column-1)]/ hySquare) +
                                                                                        (u[row* nx +(column+1)]/ hySquare) +
                                                                                        (u[(row-1) * nx +column]/ hxSquare) +
                                                                                        (u[(row+1) * nx +column]/ hxSquare) +
                                                                                        force[row* nx +column]);
            }

        }
    }



////------------------------------------------------BLACK UPDATE---------------------------------------------------------------

    for (uint row= 1; row< ny-1; ++row)
    {
        if(row & 1)
        {
            for (uint column =2; column < nx-1; column+=2)
            {
               u[row* nx + column] = 1.0/ ((2.0 / hxSquare) +(2.0 /hySquare)+ 4 * pi * pi) * ((u[row* nx +(column-1)]/ hySquare) +
                                                                                        (u[row* nx +(column+1)]/ hySquare) +
                                                                                        (u[(row-1) * nx +column]/ hxSquare) +
                                                                                        (u[(row+1) * nx +column]/ hxSquare) +
                                                                                        force[row* nx +column]);
            }
        }
    else
        {
            for (uint column =1; column < nx-1; column+=2)
            {
               u[row* nx + column] = 1.0/ ((2.0 / hxSquare) +(2.0 /hySquare)+ 4 * pi * pi) * ((u[row* nx +(column-1)]/ hySquare) +
                                                                                        (u[row* nx +(column+1)]/ hySquare) +
                                                                                        (u[(row-1) * nx +column]/ hxSquare) +
                                                                                        (u[(row+1) * nx +column]/ hxSquare) +
                                                                                        force[row* nx +column]);
            }

        }
    }
}

void PDE::Residual(void)
{
  real hxSquare = hx * hx, hySquare= hy * hy;
  
    for(uint row = 1; row < ny -1; ++row)
    {
        for(uint column = 1; column < nx-1 ; ++column)
        {
            residual[row * nx + column] = force[row * nx + column] + ((u[row* nx +(column-1)]/ hySquare) +
                                                                       (u[row* nx +(column+1)]/ hySquare) +
                                                                       (u[(row-1) * nx +column]/ hxSquare) +
                                                                       (u[(row+1) * nx +column]/ hxSquare) -
                                                                       (((2.0 / hxSquare) +(2.0 /hySquare)+ 4 * pi * pi) * u[row * nx + column]));
        }
    }
}

bool PDE::writeFile(const std::string& fileName,const std::vector<real>& vec)
{
  std::ofstream file(fileName);		//object of ofstream

  if(file.is_open())
  {   
   // file << nx << " " << ny <<"\n";
    /*
    for(auto iter = vec.begin();  iter<vec.end();++iter)
      file << *iter <<"\n";
    */
    
     for(uint row = 0; row < ny ; ++row)
      for(uint col = 0; col < nx ; ++col)
      {
	file << col*hx <<"\t"<<row*hy <<"\t"<<vec[row*nx+col] <<"\n";
      }
      std::cout << BOLD(FGRN("\n file write done  ")) <<std::endl;
    return true;
  }//if
  else
  {
    std::cerr << BOLD(FRED("\n Unable to open file  ")) <<std::endl;
    return false;
  }
  
  file.close(); 
}
