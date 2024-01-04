#include <utils.h>
#include <math.h>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <iostream>
#include <utility>


Distance_Constraint::Distance_Constraint(double d, iterator p_q1, iterator p_q2) : distance(d), diff(new double[3]), p_q1(p_q1), p_q2(p_q2) {}
Distance_Constraint::Distance_Constraint(Distance_Constraint &&other) : distance(other.distance), diff(other.diff), p_q1(other.p_q1), p_q2(other.p_q2) {other.diff = nullptr;}
Distance_Constraint::~Distance_Constraint() {delete[] diff;}

Fixpoint_Constraint::Fixpoint_Constraint(iterator p_q, iterator p_fix) : diff(new double[3]), p_q(p_q), p_fix(p_fix) {}
Fixpoint_Constraint::Fixpoint_Constraint(Fixpoint_Constraint &&other) : diff(other.diff), p_q(other.p_q), p_fix(other.p_fix) {other.diff = nullptr;}
Fixpoint_Constraint::~Fixpoint_Constraint() {delete[] diff;}

void Constraints::reserve() {values = Eigen::VectorXd(n);}

void Constraints::eval()
{
    n = distance_constraints.size() + fixpoint_constraints.size();
    iterator p_q1, p_q2;
    double *diff;
    auto p_values = values.begin();
    
    for(auto p_constraint = distance_constraints.begin(); p_constraint!=distance_constraints.end(); p_constraint++)
    {
        p_q1 = p_constraint->p_q1;
        p_q2 = p_constraint->p_q2;
        diff = p_constraint->diff;
        for (int i=0; i<3; i++) diff[i] = *p_q1++ - *p_q2++;
        *p_values++ = pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2) - pow(p_constraint->distance,2);
    }
    
    for(auto p_constraint = fixpoint_constraints.begin(); p_constraint!=fixpoint_constraints.end(); p_constraint++)
    {
        p_q1 = p_constraint->p_q;
        p_q2 = p_constraint->p_fix;
        diff = p_constraint->diff;
        for (int i=0; i<3; i++) diff[i] = *p_q1++ - *p_q2++;
        *p_values++ = pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2);
    }
}

void Constraints::add_distance_constraint(Distance_Constraint &&c) 
{
    distance_constraints.push_back(std::move(c));
    n++;
}
void Constraints::add_fixpoint_constraint(Fixpoint_Constraint &&c)
{
    fixpoint_constraints.push_back(std::move(c));
    n++;
}

Jacobian::Jacobian(int n, int m, Constraints &constraints, Eigen::VectorXd &positions) : 
    n(n), m(m), constraints(constraints), positions(positions) {}

Eigen::SparseMatrix<double> Jacobian::eval()
{   
    Eigen::SparseMatrix<double,Eigen::RowMajor> J(constraints.n, 3*n*m);
    J.reserve(Eigen::VectorXi::Constant(constraints.n,6));
    
    int offset, i=0;
    double *diff;

    //bonds
    for(auto p_constraint = constraints.distance_constraints.begin(); p_constraint!=constraints.distance_constraints.end(); i++, p_constraint++)
    {   
        diff = p_constraint->diff;
        //components for q1
        offset = p_constraint->p_q1 - positions.begin();
        for(int j=0; j<3; j++) J.insert(i,offset+j) = 2*diff[j];
        
        //components for q2
        offset = p_constraint->p_q2 - positions.begin();
        for(int j=0; j<3; j++) J.insert(i,offset+j) = -2*diff[j];
    }

    //fixpoints
    for(auto p_constraint = constraints.fixpoint_constraints.begin(); p_constraint!=constraints.fixpoint_constraints.end(); i++, p_constraint++)
    {
        diff = p_constraint->diff;
        //components for q
        offset = p_constraint->p_q - positions.begin();
        for(int j=0; j<3; j++) J.insert(i,offset+j) = 2*diff[j];
    }

    J.makeCompressed();
    return J;
}

JQDeriv::JQDeriv(int n, int m, Constraints &constraints, Eigen::VectorXd &positions, Eigen::VectorXd &velocities) :
    n(n), m(m), constraints(constraints), p_pos(positions.begin()), p_vel(velocities.begin()), result(constraints.n) {}

Eigen::VectorXd &JQDeriv::eval()
{
    int a,b;
    double *diff, value;
    auto p_result = result.begin();
    iterator p_v1, p_v2;
    
    for(auto p_constraint = constraints.distance_constraints.begin(); p_constraint!=constraints.distance_constraints.end(); p_constraint++)
    {
        value = 0;
        
        p_v1 = p_vel + (p_constraint->p_q1 - p_pos);
        p_v2 = p_vel + (p_constraint->p_q2 - p_pos);
        
        for(int j=0; j<3; j++) value += ( pow(p_v1[j],2) + pow(p_v2[j],2) - 2*p_v1[j]*p_v2[j] );

        *p_result++ = 2*value;
    }

    for(auto p_constraint = constraints.fixpoint_constraints.begin(); p_constraint!=constraints.fixpoint_constraints.end(); p_constraint++)
    {
        value = 0;
        
        p_v1 = p_vel + (p_constraint->p_q - p_pos);
        
        for(int j=0; j<3; j++) value += pow(p_v1[j],2);

        *p_result++ = 2*value;
    }

    return result;
}


Screen::Screen(int x, int y, double d1, double d2) : 
x(x), x_half(x/2), y(y), y_half(y/2), d1(d1), d2(d2), //d1: eye->screen, d2: screen->coordinate origin
front(new char[(x+3)*(y+2)+1]), back(new char[(x+3)*(y+2)+1]), zbuf_front(new double[x*y]), zbuf_back(new double[x*y]{}) 
{
    char *bar = new char[x+3];
    char *ptr = bar;

    for(int j=0; j<x+2; j++) *ptr++ = '-';
    *ptr = '\n';

    char *frame = new char [x+3];
    ptr = frame;
    *ptr++ = '|';
    for(int j=0; j<x; j++) *ptr++ = ' ';
    *ptr++ = '|';
    *ptr = '\n';

    memcpy(back, bar, (x+3)*sizeof(char));
    ptr = back+x+3;
    for(int i=0; i<y; i++, ptr+=(x+3)) memcpy(ptr, frame, (x+3)*sizeof(char));
    memcpy(ptr, bar, (x+3)*sizeof(char));
    back[(x+3)*(y+2)] = '\0';
    
    delete[] bar; delete[] frame;
    clear();
}

Screen::~Screen() {delete[] front; delete[] back; delete[] zbuf_front; delete[] zbuf_back;}
    
void Screen::clear()
{
    memcpy(front, back, ((x+3)*(y+2)+1) * sizeof(char));
    memcpy(zbuf_front,zbuf_back, x*y * sizeof(double));
}

void Screen::set(double xi, double yi, double zi, const char c) {
    
    double factor = d1/(zi+d1+d2);
    
    int xc = x_half + round(xi*factor);
    int yc = y_half - round(yi*factor);
    
    if(0<=yc && yc<y && 0<= xc && xc<x) 
    {
        if(factor > zbuf_front[yc*x + xc]) 
        {
            zbuf_front[yc*x + xc] = factor;
            front[(yc+1)*(x+3) + xc + 1] = c;
        }
    }
}

std::ostream &operator<<(std::ostream &stream, Screen &screen) {
    stream << screen.front;
    return stream;
}