#include <vector>
#include <iostream>
#include "Colsamm_Source/Colsamm.h"

using namespace ::_COLSAMM_;

class FEM
{
  private:

  
  public:
    
  FEM()
  {
    
  }
  
  ~FEM()
  {
    
  }
  
  void multigrid( double* vector_v, double* vector_b, double* vector_res, const int level, double localstiffness[][4], double localmass[][4], double vector_bc[]);
  double apply_stencil(double localstiffness[][4],double localmass[][4], double vector_x[], double vector_f[],double vector_bc[],const int cells, int row_index, const double sigma);
  void init_localstiffness(double localstiffness[][4], int cells);
  void init_localmass(double localmass[][4], int cells);
  void coarsen_local( double local_fine[][4], double local_coarse[][4], int level);
  void GSIteration(double* vector_x, double* vector_b, const int cells, double localstiffness[][4], double localmass[][4], double vector_bc[]);
  void Residual(double* vector_v, double* vector_b, double* vector_res, int cells, double localstiffness[][4], double localmass[][4], double vector_bc[]);
  void Restriction(double* vector_res, double* vector_b, int level);
  void Prolongation(double* vector_res, double* vector_v, int level);
  void Correction(double* vector_v, double* vector_res, int cells); //to be implemented
  double sigma(double x); //to be implemented
  
};