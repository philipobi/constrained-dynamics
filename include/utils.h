#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#include <vector>

typedef Eigen::VectorXd::iterator iterator;

struct Distance_Constraint
{
    double distance, *diff;
    iterator p_q1, p_q2;

    Distance_Constraint() {}
    Distance_Constraint(double d, iterator p_q1, iterator p_q2);
    Distance_Constraint(Distance_Constraint &&other);
    ~Distance_Constraint();
};

struct Fixpoint_Constraint
{
    double *diff;
    iterator p_q, p_fix;

    Fixpoint_Constraint() {}
    Fixpoint_Constraint(iterator p_q, iterator p_fix);
    Fixpoint_Constraint(Fixpoint_Constraint&& other);
    ~Fixpoint_Constraint();
};

struct Constraints
{
    int n;
    Eigen::VectorXd values;
    std::vector<Distance_Constraint> distance_constraints;
    std::vector<Fixpoint_Constraint> fixpoint_constraints;

    Constraints() {n=0;}
    
    void reserve();
    void eval();
    void add_distance_constraint(Distance_Constraint &&c);
    void add_fixpoint_constraint(Fixpoint_Constraint &&c);
};


struct Jacobian
{
    int n,m;
    Eigen::VectorXd &positions;
    Constraints &constraints;
    
    Jacobian(int n, int m, Constraints &constraints, Eigen::VectorXd &positions);

    Eigen::SparseMatrix<double> eval();
};

struct JQDeriv
{
    int n, m;
    Constraints &constraints;
    iterator p_pos, p_vel;
    Eigen::VectorXd result;

    JQDeriv(int n, int m, Constraints &constraints, Eigen::VectorXd &positions, Eigen::VectorXd &velocities);

    Eigen::VectorXd &eval();

};


struct Screen
{
    //private:

    int x, y, x_half, y_half;
    double d1, d2;
    char *front, *back; 
    double *zbuf_front, *zbuf_back;

    //public:

    Screen(int x, int y, double d1, double d2); //d1: eye->screen, d2: screen->object

    ~Screen();
    
    void clear();

    void set(double xi, double yi, double zi, const char c);

    friend std::ostream &operator<<(std::ostream &stream, Screen &screen);
};