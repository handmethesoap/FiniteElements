#include "FEM.h"

double FEM::rbgs(int level)
{
  
  
  return 0;
}

void FEM::initialise_f(double(* u_initialiser)(double))
{
  for(int i = 0; i < ((1<<levels) -1); ++i)
  {
    f[i] = u_initialiser((i+1)*step_size);
  }
}

void FEM::print_f(void){
  for(int i = 0; i < (1<< levels)-1; ++i)
  {
    std::cout << f[i] << " ";
  }
  std::cout << std::endl;
}