#include <eigen3/Eigen/Dense>

typedef Eigen::VectorXd::iterator_type vector_it;

struct Constraint {
    double value, norm, diff[3], derivative[3];
    vector_it q1, q2;

    Constraint(){};
    Constraint(double value, vector_it &q1, vector_it &q2);
    Constraint& operator=(Constraint&& other);

    void eval();
};

struct Jacobian
{
    int n,m,nConstraints;
    double *zeroes, *matrix, *positions;
    Constraint *constraints;
    
    Jacobian(int n, int m, int nConstraints, Constraint *constraints, double* positions);

    ~Jacobian();

    double *eval();
};


void printMatrix(int n, int m, double *matrix);