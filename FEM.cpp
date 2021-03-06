#include "FEM.h"

void FEM:: multigrid( double* vector_v, double* vector_b, double* vector_res, const int level, double localstiffness[][4], double localmass[][4], double vector_bc[]){
  
  double zero_bc[2] = {0.0, 0.0};
  if(level > 1)
  {
    int cells = (1<<level);
    
    for(int i = 0; i < 2; ++i)
    {
      GSIteration(vector_v, vector_b, cells, localstiffness, localmass, vector_bc);
    }
    Residual(vector_v, vector_b, vector_res, cells, localstiffness, localmass, vector_bc);
    Restriction(vector_res, vector_b+(cells+1), level);
    multigrid(vector_v + (cells+1), vector_b + (cells+1), vector_res + (cells+1), level-1, localstiffness + cells, localmass + cells, zero_bc);
    Prolongation(vector_res, vector_v + (cells + 1), level);
    Correction(vector_v, vector_res, cells);
    for(int i = 0; i<2; ++i)
    {
      GSIteration(vector_v, vector_b, cells, localstiffness, localmass, vector_bc);
    }
  }
  else
  {
    for(int i = 0; i< 1000; ++i)
      GSIteration(vector_v, vector_b, (1<<level), localstiffness, localmass, vector_bc);
  }
  
} 
  
void FEM::GSIteration(double* vector_x, double* vector_b, const int cells, double localstiffness[][4], double localmass[][4], double vector_bc[]){
  const double h = 1/cells;
  vector_x[0] = vector_bc[0];
  vector_x[cells] = vector_bc[1];

  //Gauss Seidel Red
  for(int i = 1; i < cells; i+=2)
  {
    vector_x[i] = apply_stencil(localstiffness, localmass, vector_x, vector_b, vector_bc, cells, i, sigma(i*h));
  }
  
  //Gauss Seidel Black
  for(int i = 2; i < cells; i+=2)
  {
    vector_x[i] = apply_stencil(localstiffness, localmass, vector_x, vector_b, vector_bc,cells,i,sigma(i*h));
  }
}

void FEM::coarsen_local( double local_fine[][4], double local_coarse[][4], int level){
  int n_coarse = 1<<level;
  for(int i = 0; i < n_coarse; ++i)
  {
    local_coarse[i][0] = local_fine[2*i][0] + 0.5*local_fine[2*i][1] + 0.5*local_fine[2*i][2] + 0.25*(local_fine[2*i][3]+local_fine[2*i+1][0]);
    local_coarse[i][1] = 0.5*local_fine[2*i][1] + 0.5*local_fine[2*i+1][2]+ 0.25*(local_fine[2*i][3]+local_fine[2*i+1][0]);
    local_coarse[i][2] = 0.5*local_fine[2*i+1][2] + 0.5*local_fine[2*i][1] +0.25*(local_fine[2*i][3]+local_fine[2*i+1][0]);
    local_coarse[i][3] = local_fine[2*i+1][3] + 0.5*local_fine[2*i+1][2] + 0.5*local_fine[2*i+1][1] + 0.25*(local_fine[2*i][3]+local_fine[2*i+1][0]);
  }
}

void FEM::init_localstiffness(double localstiffness[][4], int cells){
  const double h = 1./cells;
  int dim, num;

  ELEMENTS::Interval Element;
  dim = Element.dimension();
  num = Element.getNumberOfCorners();
  std::vector<std::vector<double> > Stencil;
  std::vector<double> data(num*dim,0.);
  
  for(int i = 0; i < cells; ++i)
  {
    data[0] = i*h;
    data[1] = (i+1)*h;
    Element(data);
    Stencil = Element.integrate(C_(1) * d_dx(v_())*d_dx(w_()));
    localstiffness[i][0] = Stencil[0][0];
    localstiffness[i][1] = Stencil[0][1];
    localstiffness[i][2] = Stencil[1][0];
    localstiffness[i][3] = Stencil[1][1];
  }
}

void FEM::init_localmass(double localmass[][4], int cells){
    const double h = 1./cells;
  int dim, num;

  ELEMENTS::Interval Element;
  dim = Element.dimension();
  num = Element.getNumberOfCorners();
  std::vector<std::vector<double> > Stencil;
  std::vector<double> data(num*dim,0.);
  
  for(int i = 0; i < cells; ++i)
  {
    data[0] = i*h;
    data[1] = (i+1)*h;
    Element(data);
    Stencil = Element.integrate(C_(1) * v_()*w_());
    localmass[i][0] = Stencil[0][0];
    localmass[i][1] = Stencil[0][1];
    localmass[i][2] = Stencil[1][0];
    localmass[i][3] = Stencil[1][1];
  }
}

double FEM::apply_stencil(double localstiffness[][4],double localmass[][4], double vector_x[], double vector_f[],double vector_bc[],const int cells, int row_index, const double sigma){
  double row_value = 0.;
  //right interval: i = P1[]
  // (1|1) ->index 0
  row_value += localmass[row_index][0]*vector_f[row_index];
  // (1|2) ->index 1
  row_value -= localstiffness[row_index][1]* vector_x[row_index+1];
  row_value -= sigma*localmass[row_index][1]*vector_x[row_index+1];
  row_value += localmass[row_index][1]*vector_f[row_index+1];
  //left interval: i = P2
  // (2|1) -> index 2
  row_value -= localstiffness[row_index-1][2]* vector_x[row_index-1];
  row_value -= sigma*localmass[row_index-1][2]*vector_x[row_index-1];
  row_value += localmass[row_index-1][2]*vector_f[row_index-1];
  // (2|2) -> index 3
  row_value += localmass[row_index-1][3]*vector_f[row_index];
  row_value /= (localstiffness[row_index][0] + sigma*localmass[row_index][0] +localstiffness[row_index-1][3] + sigma*localmass[row_index-1][3]);
  return row_value;
}

double FEM::sigma(double x){
  
  return 0.0;
}

void FEM::Residual(double* vector_v, double* vector_b, double* vector_res, int cells, double localstiffness[][4], double localmass[][4], double vector_bc[]){
  const double h = 1/cells;
  vector_v[0] = vector_bc[0];
  vector_v[cells] = vector_bc[1];
  
  for(int i = 1; i < cells; ++i)
  {
    vector_res[i] += localmass[i][0]*vector_b[i];
    // (1|2) ->index 1
    vector_res[i] -= localstiffness[i][1]* vector_v[i+1];
    vector_res[i] -= sigma(i*h)*localmass[i][1]*vector_v[i+1];
    vector_res[i] += localmass[i][1]*vector_b[i+1];
    //left interval: i = P2
    // (2|1) -> index 2
    vector_res[i] -= localstiffness[i-1][2]* vector_v[i-1];
    vector_res[i] -= sigma(i*h)*localmass[i-1][2]*vector_v[i-1];
    vector_res[i] += localmass[i-1][2]*vector_b[i-1];
    // (2|2) -> index 3
    vector_res[i] += localmass[i-1][3]*vector_b[i];
    vector_res[i] -= ((localstiffness[i][0] + sigma(i*h)*localmass[i][0] +localstiffness[i-1][3] + sigma(i*h)*localmass[i-1][3]))*vector_v[i];
  }
}

void FEM::Restriction(double* vector_res, double* vector_b, int level){
  for(int i = 1; i < (1<<(level-1)); ++i)
  {
    vector_b[i] = 0.25*(vector_res[2*i-1] + 2*vector_res[2*i] + vector_res[2*i +1]);  
  }
}

void FEM::Prolongation(double* vector_res, double* vector_v, int level){
  for( int i = 0; i <= (1<<level); ++i)
  {
    vector_res[i] = 0.0;
  }
  for( int i = 1; i < (1 << (level-1)); ++i)
  {
    vector_res[2*i] += vector_v[i];
    vector_res[2*i -1] += 0.5*vector_v[i];
    vector_res[2*i +1] += 0.5*vector_v[i];
  }
}

void FEM::Correction(double* vector_v, double* vector_res, int cells){
  for(int i = 1; i< cells; ++i)
  {
    vector_v[i] += vector_res[i];
  }
}