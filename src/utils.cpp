#include <cblas.h>
#include <iostream>
#include <utils.h>
#include <math.h>

Constraint::Constraint(double value, double *p_q1, double *p_q2) : value(value), p_q1(p_q1), p_q2(p_q2) {}

Constraint& Constraint::operator=(Constraint&& other) 
{
    value = other.value; p_q1 = other.p_q1; p_q2 = other.p_q2;
    other.p_q1 = nullptr; other.p_q2=nullptr;
    return *this;
}

void Constraint::eval()
{
    for(int i=0; i<3; i++) diff[i] = p_q1[i]-p_q2[i];
    norm = cblas_dnrm2(3,diff,1);
    for(int i=0; i<3; i++) derivative[i] = diff[i]/norm;
}


Jacobian::Jacobian(int n, int m, int nConstraints, Constraint *constraints, double* positions) : 
    n(n), m(m), nConstraints(nConstraints), constraints(constraints), positions(positions),
    zeroes(new double[nConstraints*3*n*m]{}), matrix(new double[nConstraints*3*n*m]) {}

Jacobian::~Jacobian() {delete[] zeroes; delete[] matrix;}

double *Jacobian::eval()
{
    memcpy(matrix, zeroes, sizeof(double)*nConstraints*3*n*m);
    
    Constraint *p_constraint = constraints;
    double *ptr = nullptr;
    
    for(int i=0; i<nConstraints; i++, p_constraint++)
    {   
        //components for q1
        ptr = matrix + i*nConstraints + (p_constraint->p_q1 - positions);
        for(int j=0; j<3; j++, ptr++) *ptr = p_constraint->derivative[j];
        
        //components for q2
        ptr = matrix + i*nConstraints + (p_constraint->p_q2 - positions);
        for(int j=0; j<3; j++, ptr++) *ptr = -p_constraint->derivative[j];
    }

    return matrix;
}

int pairs[] = {
    1,2,
    1,3,
    2,3
    };
int pairs1[] = {
    1,5,
    1,6,
    2,4,
    2,6,
    3,4,
    3,5
    };
int pairs2[] = {
    4,5,
    4,6,
    5,6
    };

struct JQDeriv
{
    int n, m, nConstraints;
    Constraint *constraints;
    double *positions, *velocities, *result;

    JQDeriv(int n, int m, int nConstraints, Constraint * constraints, double* positions, double* velocities) :
    n(n), m(m), nConstraints(nConstraints), constraints(constraints), positions(positions), velocities(velocities), 
    result(new double[nConstraints]) {}

    ~JQDeriv() {delete[] result;}

    double *eval()
    {
        Constraint *p_constraints = constraints;
        double *p_result = result;
        double *p_q1, *p_q2, *p_v1, *p_v2,*diff;
        int a,b;
        for(int i=0; i<nConstraints; i++, p_constraints++, p_result++)
        {
            p_q1 = p_constraints->p_q1;
            p_q2 = p_constraints->p_q2;
            p_v1 = velocities + (p_q1-positions);
            p_v2 = velocities + (p_q2-positions);
            diff = p_constraints->diff;
            
            double c0 = 1/p_constraints->norm;
            double c03 = pow(c0,3);

            double value = 0.;

            for(int j=0; j<3; j++)
            {
                value += (c0 - c03 * pow(diff[j],2)) * (pow(p_v1[j],2) + pow(p_v2[j],2) -2*p_v1[j]*p_v2[j]);
            }

            for(int *ptr = pairs; ptr < ptr+3*2; ptr++)
            {
                a = *ptr++ -1;
                b = *ptr -1;
                
                value -= 2*(c03 * diff[a] * diff[b] * p_v1[a] * p_v1[b]);
            }
            
            for(int *ptr = pairs1; ptr < ptr+6*2; ptr++)
            {
                a = *ptr++ -1;
                b = *ptr -4;
                
                value += 2*(c03 * diff[a] * diff[b] * p_v1[a] * p_v2[b]);
            }
            
            for(int *ptr = pairs2; ptr < ptr+3*2; ptr++)
            {
                a = *ptr++ -4;
                b = *ptr -4;
                
                value -= 2*(c03 * diff[a] * diff[b] * p_v2[a] * p_v2[b]);
            }

            *p_result = value;
        }

        return result;
    }
    
};

void printMatrix(int n, int m, double *matrix)
{
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++) std::cout << matrix[i*n+j] << '\t';
        std::cout << '\n';
    }
}