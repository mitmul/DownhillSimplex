#include "downhill.h"

DownhillSimplex::DownhillSimplex(const int _simplex_num, const int _dimension)
  : iter_num(0),
    simplex_num(_simplex_num),
    dimension(_dimension),
    alpha(1.0),
    beta(0.5),
    gamma(2.0)
{
  // ランダムに初期simplexを配置
  srand((unsigned int)time(NULL));
  for(int i = 0; i < simplex_num; ++i)
  {
    VectorXd s(dimension);
    for(int j = 0; j < dimension; ++j)
      s[j] = random(-1.0, 1.0);

    simplex.push_back(s);
  }

  x_g.resize(dimension);

  evalAll();
  calcParams();
}

void DownhillSimplex::setVariableRange(vector<double> mins, vector<double> maxs)
{
  // ランダムに初期simplexを配置
  simplex.clear();
  for(int i = 0; i < simplex_num; ++i)
  {
    VectorXd s(dimension);
    for(int j = 0; j < dimension; ++j)
      s[j] = random(mins[j], maxs[j]);

    simplex.push_back(s);
  }

  evalAll();
  calcParams();
}

void DownhillSimplex::converge(const double epsilon, const int iter_limit)
{
  while(criterion() > epsilon && iter_num < iter_limit)
  {
    optimize();
    ++iter_num;
  }
}

double DownhillSimplex::cost(const VectorXd &x)
{
  double a = 10.0;
  double b = 0.0;
  for(int i = 0; i < dimension; ++i)
  {
    b += pow(x[i], 2) - a * cos(2 * M_PI * x[i]);
  }

  return a * dimension + b;
}

VectorXd DownhillSimplex::getResult()
{
  return x_l;
}

double DownhillSimplex::getFinalCost()
{
  return x_l_eval;
}

void DownhillSimplex::evalAll()
{
  evaluation.clear();
  for(int i = 0; i < (int)simplex.size(); ++i)
  {
    double c = cost(simplex[i]);
    evaluation.push_back(c);
  }
}

void DownhillSimplex::calcParams()
{
  worst_pos = 0;
  double worst_val = evaluation[0];
  for(int i = 0; i < simplex_num; ++i)
  {
    // 最も値が大きい位置
    if(evaluation[i] > worst_val)
    {
      worst_pos = i;
      worst_val = evaluation[i];
    }
  }
  x_h = simplex[worst_pos];

  int best_pos = 0;
  double best_val = evaluation[0];
  for(int i = 0; i < simplex_num; ++i)
  {
    // 最も値が小さい位置
    if(evaluation[i] < best_val)
    {
      best_pos = i;
      best_val = evaluation[i];
    }
  }
  x_l = simplex[best_pos];

  int next_worst_pos = 0;
  double next_worst_val = evaluation[best_pos];
  for(int i = 0; i < simplex_num; ++i)
  {
    // worst_posを除いて
    if(i != worst_pos)
    {
      // 最も値が大きい位置
      if(evaluation[i] > next_worst_val)
      {
        next_worst_pos = i;
        next_worst_val = evaluation[i];
      }
    }
  }
  x_s = simplex[next_worst_pos];

  // worstを除いた重心
  x_g.setZero();
  for(int i = 0; i < simplex_num; ++i)
  {
    // worst_posを除いて
    if(i != worst_pos)
    {
      x_g += simplex[i];
    }
  }
  x_g /= (double)(simplex_num - 1);

  x_h_eval = cost(x_h);
  x_s_eval = cost(x_s);
  x_l_eval = cost(x_l);

//  for(int i = 0; i < dimension; ++i)
//    cout << x_l[i] << "\t";
//  cout << " : " << x_l_eval << endl;
}

VectorXd DownhillSimplex::reflect()
{
  return (1 + alpha) * x_g - alpha * x_h;
}

VectorXd DownhillSimplex::expand(const VectorXd x_r)
{
  return gamma * x_r + (1 - gamma) * x_g;
}

VectorXd DownhillSimplex::contract()
{
  return beta * x_h + (1- beta) * x_g;
}

void DownhillSimplex::shrink()
{
  for(int i = 0; i < simplex_num; ++i)
  {
    simplex[i] = 0.5 * (x_l + simplex[i]);
  }
}

void DownhillSimplex::optimize()
{
  VectorXd x_r = reflect();
  double x_r_eval = cost(x_r);

  if(x_r_eval < x_l_eval)
  {
    VectorXd x_e = expand(x_r);
    double x_e_eval = cost(x_e);

    if(x_e_eval < x_r_eval)
    {
      simplex[worst_pos] = x_e;
    }
    else
    {
      simplex[worst_pos] = x_r;
    }
  }
  else if(x_s_eval < x_r_eval && x_r_eval < x_h_eval)
  {
    simplex[worst_pos] = x_r;

    VectorXd x_c = contract();
    double x_c_eval = cost(x_c);

    if(x_c_eval < x_h_eval)
    {
      simplex[worst_pos] = x_c;
    }
    else
    {
      shrink();
    }
  }
  else
  {
    simplex[worst_pos] = x_r;
  }

  evalAll();
  calcParams();
}

double DownhillSimplex::criterion()
{
  return 2 * fabs(x_h_eval - x_l_eval) / (fabs(x_h_eval) + fabs(x_l_eval));
}

double DownhillSimplex::random(const double min, const double max)
{
  // 0〜1の乱数
  double r = (double)rand() / RAND_MAX;

  // 乱数の幅
  double length = fabs(max - min);

  r *= length;
  r += min;

  return r;
}
