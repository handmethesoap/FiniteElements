#include "Colsamm_Source/Colsamm.h"
#include "FEM.h"
#include <stdlib.h>
#include <vector>
#include "initialisers.h"

int main(int argc, char* argv[])
{
  if(argc != 2)
  {
    std::cout << "A grid size is needed" << std::endl;
    return 0;
  }
  
  int L = atoi(argv[1]);
  int array_length = 0;
  int ghost_array_length = 0;
  int cells = (1<<L);
  double h = 1/(cells);
  
  for(int i = 1; i <= L; ++i)
  {
    array_length += (1<<i); 
    ghost_array_length += (1<<i) +1;
  }
  FEM femproblem;
  
  double* v = new double[ghost_array_length];
  double* f = new double[ghost_array_length];
  double* r = new double[ghost_array_length];
  for(int i = 0; i < ghost_array_length; ++i)
  {
    v[i] = 0.0;
    f[i] = 0.0;
    r[i] = 0.0;
  }
  
  double (*stiffness)[4] = new double[array_length][4];
  double (*mass)[4] = new double[array_length][4];
  
  double* b = new double[2];
  b[0] = 1;
  b[1] = 1;
  
  femproblem.init_localstiffness(stiffness, (cells));
  femproblem.init_localmass(mass, (cells));
  int position = 0;
  
  for(int i = L; i > 1; --i)
  {
    femproblem.coarsen_local(stiffness + position, stiffness + position + (1<<i), i-1);
    femproblem.coarsen_local(mass + position, mass + position + (1<<i), i-1);
    position += (1<<i);
  }
  
  ELEMENTS::Interval my_element;
  for( int i = 0; i < 50; i++ )
  {
    femproblem.multigrid(v,f,r,L, stiffness, mass, b);
  }

  
  std::cout << "--------------Result--------------" << std::endl;
  for(int i = 0; i <= (cells); ++i)
  {
    std::cout << "position: " << i*h << ", value: " << v[i] << std::endl;
  }
    
  
  delete[] v;
  delete[] f;
  delete[] r;
  delete[] stiffness;
  delete[] mass;
  delete[] b;
  
  return 0;
}
