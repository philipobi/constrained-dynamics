#include <stdio.h>
#include <cstdlib>

struct Node {
    double pos[3];
    int i,j;
    Node(int i, int j, double x, double y, double z) : i(i), j(j), pos{x,y,z} {}
};

struct DistanceConstraint {
    double value;
    Node *a,*b;
    DistanceConstraint(double value, Node *a, Node *b) : value(value), a(a), b(b) {}
};



void generateLattice(int n, int m, double l, Node *nodes, DistanceConstraint *constraints)
{
    
    for(int i=0; i<n; i++) for(int j=0; j<m; j++) nodes[i*n+j] = Node(i+1, j+1, i*l, j*l, 0);

    //horizontal bonds
    for(int i=0; i<n; i++) for(int j=0; j<(m-1); j++) constraints[i*(m-1)+j] = DistanceConstraint (l, &nodes[i*n+j], &nodes[i*n+j+1]);

    //vertical bonds
    for(int j=0; j<m; j++) for(int i=0; i<(n-1); i++) constraints[n*(m-1)+j*(n-1)+i] = DistanceConstraint (l, &nodes[i*n+j], &nodes[(i+1)*n+j]);

}

void printConstraints(DistanceConstraint *constraints, int n, int m)
{
    for(int i=0; i<((m-1)*n + (n-1)*m); i++) 
    {
        Node *a = constraints[i].a;
        Node *b = constraints[i].b;
        printf("Constraint between (%d,%d) and (%d,%d)\n", a->i, a->j, b->i, b->j);
    }
}

int main()
{
    int n = 3;
    int m = 3;

    Node *nodes = (Node *) malloc(n*m * sizeof(Node));
    DistanceConstraint *constraints = (DistanceConstraint *) malloc(((m-1)*n + (n-1)*m) * sizeof(DistanceConstraint));

    generateLattice(n,m,3.,nodes,constraints);

    printConstraints(constraints, n, m);

    free(nodes); free(constraints);
}