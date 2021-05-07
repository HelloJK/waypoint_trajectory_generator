#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <OsqpEigen/OsqpEigen.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}
/*

    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGenerationCloseForm function

    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients

*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGenerationCloseForm(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */

    /*   Produce the dereivatives in X, Y and Z axis directly.  */

    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */


    return PolyCoeff;
}


Eigen::MatrixXd TrajectoryGeneratorWaypoint::GetQk(double t_k, int p_order) {
    MatrixXd Q_k = Eigen::MatrixXd::Zero(p_order+1, p_order+1);
    for (int i = 4; i <= p_order; ++i) {
        for (int l = 4; l <= p_order; ++l) {
            double den = i+l-7;
            Q_k(i, l) = i*(i-1)*(i-2)*(i-3)*l*(l-1)*(l-2)*(l-3) / den * pow(t_k, den);
        }
    }
}



Eigen::MatrixXd TrajectoryGeneratorWaypoint::GetQ(Eigen::VectorXd time, int p_order) {
  auto Q_k = [](double t_k, int p_order) -> Eigen::MatrixXd {
      Eigen::MatrixXd Q_k_t = Eigen::MatrixXd::Zero(p_order + 1, p_order + 1);
      for (int i = 4; i <= p_order; ++i) {
        for (int l = 4; l <= p_order; ++l) {
            double den = i+l-7;
            Q_k_t(i, l) = i*(i-1)*(i-2)*(i-3)*l*(l-1)*(l-2)*(l-3) / den * pow(t_k, den);
        }
      }
      return Q_k_t;
  };

  int dims = time.size() * (p_order + 1);
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(dims, dims);
  int p_num1d = p_order + 1;
  for (int i = 0; i < time.size(); ++i) {
      double time_now = time(i) ;
      Q.block(i * p_num1d, i * p_num1d, p_num1d, p_num1d) = Q_k(time_now, p_order);
  }
  // ROS_INFO("[TRAJ_GEN] Q: ", Q);
  // std::cout << "Q: " << std::endl << Q << std::endl;

  return Q;
}


void TrajectoryGeneratorWaypoint::GetCons(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time,          // time allocation in each segment
            int dim_index, // 0:x, 1:y, 2:z
            Eigen::MatrixXd& Aeq,
            Eigen::MatrixXd& beq) {
  // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
  int p_order   = 2 * d_order - 1;              // the order of polynomial
  int p_num1d   = p_order + 1;                  // the number of variables in each segment
  int seg_num = Time.size();                          // the number of segments
  int wp_num = Path.rows();
  int param_num = p_num1d * seg_num;
  
  // start point
  Eigen::MatrixXd Aeq_start = Eigen::MatrixXd::Zero(4, param_num);
  for (int k = 0; k < 4; ++k) { // p, v, a, j
    for (int i = k; i <= p_order; ++i) {
      Aeq_start(k, i) = Factorial(i) / Factorial(i-k) * pow(0, i-k);
    }
  }
  Eigen::MatrixXd beq_start = Eigen::MatrixXd::Zero(4, 1); // start from (0, 0, 0, 0)
  beq_start(0, 0) = Path(0, dim_index);
  beq_start(1, 0) = Vel(0, dim_index);
  beq_start(2, 0) = Acc(0, dim_index);
  beq_start(3, 0) = 0;

  
  // end point
  Eigen::MatrixXd Aeq_end = Eigen::MatrixXd::Zero(4, param_num);
  for (int k = 0; k < 4; ++k) {  // p, v, a, j
    for (int i = k; i <= p_order; ++i) {
      Aeq_end(k, p_num1d * (seg_num-1)+i) = Factorial(i) / Factorial(i-k) * pow(Time(seg_num-1), i-k);
    }
  }
  Eigen::MatrixXd beq_end = Eigen::MatrixXd::Zero(4, 1);// Path.row(wp_num - 1).transpose();
  beq_end(0, 0) = Path(wp_num - 1, dim_index);  // x
  beq_end(1, 0) = Vel(1, dim_index);
  beq_end(2, 0) = Acc(1, dim_index);
  beq_end(3, 0) = 0;

  // middle way point
  Eigen::MatrixXd Aeq_wp = Eigen::MatrixXd::Zero(seg_num-1, param_num);
  for (int n = 0; n < seg_num-1; ++n) {
    int index = p_num1d * n;
    double T = Time(n);
    for (int i = 0; i < p_num1d; ++i) {
      Aeq_wp(n, index + i) = pow(T, i);
    }
  }
  Eigen::MatrixXd beq_wp = Eigen::MatrixXd::Zero(seg_num-1, 1);
  for (int i = 0; i < seg_num - 1; ++i) {
    beq_wp(i, 0) = Path(i+1, dim_index); // x
  }
  
  // continus constraint
  Eigen::MatrixXd Aeq_cons_p = Eigen::MatrixXd::Zero(seg_num-1, param_num);
  int k = 0;
  for (int n = 0; n < seg_num - 1; ++n) {
    int index1 = p_num1d * n;
    double T1 = Time(n);
    int index2 = p_num1d * (n + 1);
    double T2 = 0; // Time(n+1);
    for (int i = k; i < p_num1d; ++i) {
      Aeq_cons_p(n, index1 + i) = Factorial(i) / Factorial(i-k) * pow(T1, i-k); // pow(T1, i);
      Aeq_cons_p(n, index2 + i) = -Factorial(i) / Factorial(i-k) * pow(T2, i-k);  // -pow(T2, i);
    }
  }
  Eigen::MatrixXd Aeq_cons_v = Eigen::MatrixXd::Zero(seg_num-1, param_num);
  k = 1;
  for (int n = 0; n < seg_num - 1; ++n) {
    int index1 = p_num1d * n;
    double T1 = Time(n);
    int index2 = p_num1d * (n + 1);
    double T2 = 0; // Time(n+1);
    for (int i = k; i < p_num1d; ++i) {
      Aeq_cons_v(n, index1 + i) = Factorial(i) / Factorial(i-k) * pow(T1, i-k); // pow(T1, i);
      Aeq_cons_v(n, index2 + i) = -Factorial(i) / Factorial(i-k) * pow(T2, i-k);  // -pow(T2, i);
    }
  }
  Eigen::MatrixXd Aeq_cons_a = Eigen::MatrixXd::Zero(seg_num-1, param_num);
  k = 2;
  for (int n = 0; n < seg_num - 1; ++n) {
    int index1 = p_num1d * n;
    double T1 = Time(n);
    int index2 = p_num1d * (n + 1);
    double T2 = 0; // Time(n+1);
    for (int i = k; i < p_num1d; ++i) {
      Aeq_cons_a(n, index1 + i) = Factorial(i) / Factorial(i-k) * pow(T1, i-k); // pow(T1, i);
      Aeq_cons_a(n, index2 + i) = -Factorial(i) / Factorial(i-k) * pow(T2, i-k);  // -pow(T2, i);
    }
  }
  Eigen::MatrixXd Aeq_cons_j = Eigen::MatrixXd::Zero(seg_num-1, param_num);
  k = 3;
  for (int n = 0; n < seg_num - 1; ++n) {
    int index1 = p_num1d * n;
    double T1 = Time(n);
    int index2 = p_num1d * (n + 1);
    double T2 = 0; // Time(n+1);
    for (int i = k; i < p_num1d; ++i) {
      Aeq_cons_j(n, index1 + i) = Factorial(i) / Factorial(i-k) * pow(T1, i-k); // pow(T1, i);
      Aeq_cons_j(n, index2 + i) = -Factorial(i) / Factorial(i-k) * pow(T2, i-k);  // -pow(T2, i);
    }
  }
  Eigen::MatrixXd Aeq_con = Eigen::MatrixXd::Zero(4*(seg_num-1), param_num);
  Aeq_con.block(0, 0, seg_num-1, param_num) = Aeq_cons_p;
  Aeq_con.block(seg_num-1, 0, seg_num-1, param_num) = Aeq_cons_v;
  Aeq_con.block(2*(seg_num-1), 0, seg_num-1, param_num) = Aeq_cons_a;
  Aeq_con.block(3*(seg_num-1), 0, seg_num-1, param_num) = Aeq_cons_j;
  Eigen::MatrixXd beq_con = Eigen::MatrixXd::Zero(4 * (seg_num-1), 1);

  int total_constraint = Aeq_start.rows() + Aeq_end.rows() + Aeq_wp.rows() + Aeq_con.rows();
  // Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(total_constraint, param_num);
  Aeq.resize(total_constraint, param_num);
  Aeq << Aeq_start,
          Aeq_end,
          Aeq_wp,
          Aeq_con;
  // Eigen::MatrixXd beq = Eigen::MatrixXd::Zero(total_constraint, 1);
  beq.resize(total_constraint, 1);
  beq << beq_start,
          beq_end,
          beq_wp,
          beq_con; 
}

void TrajectoryGeneratorWaypoint::CastMatrix2SparseMatrix(const Eigen::MatrixXd& mat, Eigen::SparseMatrix<double>& spar_mat) {
  spar_mat.resize(mat.rows(), mat.cols());
  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      if (mat(i, j) != 0) {
        spar_mat.insert(i, j) = mat(i, j);
      }
    }
  }
}

bool TrajectoryGeneratorWaypoint::Solve(const Eigen::MatrixXd& Q,
                              const Eigen::VectorXd& p,
                              const Eigen::MatrixXd& Aeq,
                              const Eigen::VectorXd& beq,
                              Eigen::VectorXd& result) {
  
  solver.settings()->setWarmStart(false);
  solver.data()->setNumberOfVariables(Q.rows());
  solver.data()->setNumberOfConstraints(Aeq.rows());

  Eigen::SparseMatrix<double> hessian;
  CastMatrix2SparseMatrix(Q, hessian);
  solver.data()->clearHessianMatrix();
  if (!solver.data()->setHessianMatrix(hessian)) {
    ROS_ERROR("[TRAJ_GEN] set hessian error");
    return false;
  }
  Eigen::VectorXd gradient = p;
  if (!solver.data()->setGradient(gradient)) {
    ROS_ERROR("[TRAJ_GEN] set gradient error");
    return false;
  }
  Eigen::SparseMatrix<double> linearMat;
  CastMatrix2SparseMatrix(Aeq, linearMat);
  solver.data()->clearLinearConstraintsMatrix();
  if (!solver.data()->setLinearConstraintsMatrix(linearMat)) {
    ROS_ERROR("[TRAJ_GEN] set Aeq error");
    return false;
  }
  Eigen::VectorXd lower_bound = beq;
  if (!solver.data()->setLowerBound(lower_bound)) {
    ROS_ERROR("[TRAJ_GEN] set lower bound error");
    return false;
  }
  Eigen::VectorXd upper_bound = beq;
  if (!solver.data()->setUpperBound(upper_bound)) {
    ROS_ERROR("[TRAJ_GEN] set upper bound error");
    return false;
  }
  if (!solver.initSolver()) {
    ROS_ERROR("[TRAJ_GEN] init solver error");
    return false;
  }
  
  if (!solver.solve()) {
    ROS_ERROR("[TRAJ_GEN] solve error");
    return false;
  }

  Eigen::VectorXd res = solver.getSolution();
  result = res;

  // std::cout << "result:" << std::endl << result << std::endl;
  solver.clearSolver();
  return true;
}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
  int p_order   = 2 * d_order - 1;              // the order of polynomial
  int p_num1d   = p_order + 1;                  // the number of variables in each segment
  int seg_num = Time.size();                          // the number of segments
  int wp_num = Path.rows();
  int param_num = p_num1d * seg_num;
  
  // std::cout << "time: " << std::endl << Time << std::endl;
  
  Eigen::MatrixXd Q = GetQ(Time, p_order);
  Eigen::VectorXd p = Eigen::VectorXd::Zero(param_num);
    
  Eigen::MatrixXd Aeq_x;
  Eigen::MatrixXd beq_x;
  int x_dim_index = 0;
  GetCons(d_order, Path, Vel, Acc, Time, x_dim_index, Aeq_x, beq_x);
  // std::cout << "A:" << std::endl << Aeq_x << std::endl;
  // std::cout << "b:" << std::endl << beq_x << std::endl;
  Eigen::VectorXd x_solution = Eigen::VectorXd::Zero(param_num);
  if (!Solve(Q, p, Aeq_x, beq_x, x_solution)) {
    ROS_ERROR("[TRAJ_GEN] solve x error");
    return Eigen::MatrixXd::Zero(1, 1);
  }
  // std::cout << "solve x succeed" << std::endl;
  // std::cout << "x_param: " << std::endl << x_solution << std::endl;

  Eigen::MatrixXd Aeq_y;
  Eigen::MatrixXd beq_y;
  int y_dim_index = 1;
  GetCons(d_order, Path, Vel, Acc, Time, y_dim_index, Aeq_y, beq_y);
  Eigen::VectorXd y_solution = Eigen::VectorXd::Zero(param_num);
  if (!Solve(Q, p, Aeq_y, beq_y, y_solution)) {
    ROS_ERROR("[TRAJ_GEN] solve y error");
    return Eigen::MatrixXd::Zero(1, 1);
  }
  // std::cout << "y_param: " << std::endl << y_solution << std::endl;

  Eigen::MatrixXd Aeq_z;
  Eigen::MatrixXd beq_z;
  int z_dim_index = 2;
  GetCons(d_order, Path, Vel, Acc, Time, z_dim_index, Aeq_z, beq_z);
  Eigen::VectorXd z_solution = Eigen::VectorXd::Zero(param_num);
  if (!Solve(Q, p, Aeq_z, beq_z, z_solution)) {
    ROS_ERROR("[TRAJ_GEN] solve z error");
    return Eigen::MatrixXd::Zero(1, 1);
  }
  // std::cout << "z_param: " << std::endl << z_solution << std::endl;

  Eigen::MatrixXd result;
  result.resize(param_num, 3);
  result << x_solution, y_solution, z_solution;
  ROS_INFO("[TRAJ_GEN] result size: %d %d %d", result.rows(), result.cols(), x_solution.size());
  return result;
}