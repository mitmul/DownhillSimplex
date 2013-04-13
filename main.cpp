#include <iostream>
#include "downhill.h"

using namespace std;

int main()
{
  DownhillSimplex opt(100, 2);
  opt.converge(1E-10, 10000);

  VectorXd result = opt.getResult();
  for(int i = 0; i < result.rows(); ++i)
  {
    printf("%.5f\t", result[i]);
  }
  printf("\n");

  printf("cost: %.10f\n", opt.getFinalCost());

  return 0;
}

