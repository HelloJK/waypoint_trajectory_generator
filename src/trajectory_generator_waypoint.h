#ifndef _TRAJECTORY_GENERATOR_WAYPOINT_H_
#define _TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <OsqpEigen/OsqpEigen.h>
#include <vector>

class TrajectoryGeneratorWaypoint {
  private:
		double _qp_cost;
		Eigen::MatrixXd _Q;
		Eigen::VectorXd _Px, _Py, _Pz;
    OsqpEigen::Solver solver;
  public:
    TrajectoryGeneratorWaypoint();

    ~TrajectoryGeneratorWaypoint();

    Eigen::MatrixXd PolyQPGeneration(
        const int order,
        const Eigen::MatrixXd &Path,
        const Eigen::MatrixXd &Vel,
        const Eigen::MatrixXd &Acc,
        const Eigen::VectorXd &Time);
    
    Eigen::MatrixXd PolyQPGenerationCloseForm(
        const int order,
        const Eigen::MatrixXd &Path,
        const Eigen::MatrixXd &Vel,
        const Eigen::MatrixXd &Acc,
        const Eigen::VectorXd &Time);
    
    int Factorial(int x);

    Eigen::MatrixXd GetQk(double t_k, int p_order);
    
    Eigen::MatrixXd GetQ(Eigen::VectorXd time, int p_order);

    void GetCons(
        const int d_order,                    // the order of derivative
        const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
        const Eigen::MatrixXd &Vel,           // boundary velocity
        const Eigen::MatrixXd &Acc,           // boundary acceleration
        const Eigen::VectorXd &Time,          // time allocation in each segment
        int dim_index, // 0:x, 1:y, 2:z
        Eigen::MatrixXd& Aeq,
        Eigen::MatrixXd& beq);
    
    bool Solve(const Eigen::MatrixXd& Q,
                const Eigen::VectorXd& p,
                const Eigen::MatrixXd& Aeq,
                const Eigen::VectorXd& beq,
                Eigen::VectorXd& result);
    void CastMatrix2SparseMatrix(const Eigen::MatrixXd& mat, Eigen::SparseMatrix<double>& spar_mat);
};
        

#endif
