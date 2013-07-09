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
  
  FEM feproblem(L);
  feproblem.set_boundaries(0.0,0.0);
  feproblem.initialise_f(ramp);
  feproblem.print_f();
  
  return 0;
}
