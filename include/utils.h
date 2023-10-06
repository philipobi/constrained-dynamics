struct Constraint {
    double value, *p_q1, *p_q2;
    double diff[3], norm, derivative[3];

    Constraint(){};
    Constraint(double value, double *p_q1, double *p_q2);
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