

#include "../include/type.hpp"

#include "../include/Timer.h"
#include "../include/PDE.hpp"


//------------------------------------------------------------------------
// this function is evaluate at commpile time

inline float fraction(size_t A,size_t B)
{
    float f = static_cast<float>(A) / B;
    return f ;
}

//meta-programming
template <int A,int B>
struct Fraction{
  enum { val = A / B };
};

template <uint A,uint B>
struct Product{
    enum { value = A * B };
};

//--------------------------------------------------------

int main(int argc, char** argv) 
{
  if (argc != 4 ){
       std::cerr<< "Not Enough Arguments "<<std::endl;
       return 0;
   }
  else
  {
    siwir::Timer time;
    
   //std::cout<<Product<2,3>::value<<std::endl;
//---------------------------------------------------------------------
    
    uint nx = std::atoi(argv[1]) , ny = std::atoi(argv[2]); 
    real hy = 1.0 / ny;
    real hx = 2.0 / nx;
  
      
    PDE p( nx , ny , hx ,hy);
    p.applyBoundary();

    time.reset();
    //for(uint i= 0 ;i < 10;++i)
    {      
      p.RedBlackGaussSeidal(std::atoi(argv[3]));
    }
      real runtime = time.elapsed();
     
      std::cout<<"runtime "<<runtime<<std::endl;
   
      std::cout<<"L2 Norm "<< p.ResidualNorm()<<std::endl;
    
      p.writeFile("solution.txt",p.u);
  }
  return 0;
}
