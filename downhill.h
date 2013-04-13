#ifndef DOWNHILL_H
#define DOWNHILL_H

#include <vector>
#include <iostream>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

class DownhillSimplex
{
  public:
    DownhillSimplex(const int _simplex_num, const int _dimension);
    void setVariableRange(vector<double> mins, vector<double> maxs);
    void converge(const double epsilon, const int iter_limit);

    virtual double cost(const VectorXd &x);

    VectorXd getResult();
    double getFinalCost();

  private:
    int iter_num;
    int simplex_num;
    int dimension;
    double alpha;
    double beta;
    double gamma;
    vector<VectorXd> simplex;
    vector<double> evaluation;

    VectorXd x_h;
    VectorXd x_s;
    VectorXd x_l;
    VectorXd x_g;
    double x_h_eval;
    double x_s_eval;
    double x_l_eval;
    int worst_pos;

    void evalAll();
    void calcParams();

    VectorXd reflect();
    VectorXd expand(const VectorXd x_r);
    VectorXd contract();
    void shrink();

    void optimize();
    double criterion();

    double random(const double min, const double max);
};

#endif // DOWNHILL_H
