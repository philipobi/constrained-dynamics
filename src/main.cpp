#include <iostream>
#include <utils.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/OrderingMethods>
#include <chrono>
#include <thread>

#define N 3
#define M 3
#define L 1
#define FRAMETIME 500

void grid(Eigen::VectorXd &positions, Eigen::VectorXd &fixpts, Constraints &constraints)
{
    auto p_pos = positions.begin();
    
    for(int i=0; i<N; i++) for(int j=0; j<M; j++) 
    {
        *p_pos++ = i*L; //x
        *p_pos++ = j*L; //y
        *p_pos++ = 0; //z
    }
    
    p_pos = positions.begin();

    //horizontal bonds
    for(int i=0; i<N; i++) for(int j=0; j<(M-1); j++) constraints.add_distance_constraint(
        Distance_Constraint(L, p_pos+3*(i*M+j), p_pos+3*(i*M+j+1))
        );

    //vertical bonds
    for(int j=0; j<M; j++) for(int i=0; i<(N-1); i++) constraints.add_distance_constraint(
        Distance_Constraint(L, p_pos+3*(i*M+j), p_pos+3*((i+1)*M+j))
        );

    //fixpoints
    for(int j=0; j<M; j++)
    {
        fixpts(3*j+1) = j*L;
        constraints.add_fixpoint_constraint(Fixpoint_Constraint(p_pos+3*j, fixpts.begin()+3*j));
        
    }
}


int main()
{    
    int n_coords = 3*N*M;
    Eigen::VectorXd q(n_coords), v(n_coords), Q(n_coords);
    Eigen::VectorXd f(Eigen::VectorXd::Zero(Eigen::Index(3*M)));

    for(int i=0; i<N*M; i++) Q(3*i+2) = -9.81;

    Constraints c;
    
    grid(q,f,c);

    c.reserve();

    Jacobian Jacobian(N,M,c,q);

    JQDeriv Jdq(N,M,c,q,v);

    double dt = double(FRAMETIME) / double(1000);

    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

    Eigen::SparseMatrix<double> J, JT, A;
    
    Eigen::VectorXd y, dq, Qc;

    Screen s(50,50,1,1);

    Eigen::VectorXd::iterator it;

    double x_,y_,z_;

    while(1)
    {
        
        c.eval();

        dq = v * dt;
        
        std::cout << q << std::endl;

        J = Jacobian.eval();

        std::cout << J << std::endl;

        JT = J.transpose();
        A = J*JT;
        A.makeCompressed();

        solver.compute(A);
        
        //Ax = b <=> (J*J^T) * y = -J*dq - J*Q
        y = solver.solve(-Jdq.eval() - J*q);
        
        v += (Q + JT*y) * dt;

        q += dq;

        std::cout << q << std::endl;

        std::this_thread::sleep_for(std::chrono::milliseconds(FRAMETIME));
    }
}