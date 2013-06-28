#include <vector>
#include <iostream>

class FEM
{
  private:
    
  std::vector<double> u;
  std::vector<double> f;
  std::vector<double> bc;
  std::vector< std::vector<double> > local_stiffness;
  std::vector< std::vector<double> > local_mass;
  int levels;
  double step_size;
  double omega;
  
  public:
    
  FEM(int _levels)
  {
    levels = _levels;
    int vector_length = 0;
    
    for(int i = 1; i<= levels; ++i)
    {
      vector_length += ((1 << i) - 1);
    }
    
    u.resize(vector_length, 1.0);
    f.resize(vector_length);
    local_stiffness.resize(vector_length, std::vector<double>(4));
    local_mass.resize(vector_length), std::vector<double>(4);
    
    step_size = 1.0/(1<<levels);
 
    
    
    #ifdef DEBUG
    
	std::cout << vector_length << std::endl;
	for( std::vector<double>::const_iterator i = u.begin(); i != u.end(); ++i)
	{
	  std::cout << *i << std::endl;
	}
    #endif
  }
  
  ~FEM()
  {
    
  }
  
  void set_boundaries(double left_boundary, double right_boundary)
  {
    bc.push_back(left_boundary);
    bc.push_back(right_boundary);
    
    #ifdef DEBUG
	for( std::vector<double>::const_iterator i = bc.begin(); i != bc.end(); ++i)
	{
	  std::cout << *i << std::endl;
	}
    #endif
  }
  
  void set_omega(double _omega)
  {
    omega = _omega;
  }
  
  void initialise_f(double(* u_initialiser)(double));
  double rbgs(int level);
  void print_f(void);
  
};