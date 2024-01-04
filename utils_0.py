import numpy as np
from sympy import Symbol, Matrix, simplify
from scipy import linalg

limit = 0.1

distance_constraint = lambda v1,v2,L: (Matrix(v1) - Matrix(v2)).norm()**2 - L**2
fixpoint_constraint = lambda v1,v2: (Matrix(v1) - Matrix(v2)).norm()**2

class Simulation:
    def __init__(self,n,m,L):
        #setup coordinates
        lst = []
        lst1 = []
        for i in range(n):
            for j in range(m):
                lst.append(i*L) #x
                lst.append(j*L) #y
                lst.append(0) #z

                lst1.append(0) #Fx
                lst1.append(0) #Fy
                lst1.append(-9.81) #Fz

        self.q = np.array(lst,dtype=np.float64)
        self.Q = np.array(lst1)

        self.N = 3*n*m

        self.v = np.zeros(self.N)

        self.q_sym = []
        self.dq_sym = []

        for i in range(self.N):
            self.q_sym.append(Symbol("q_{%d}"%i, real=True))
            self.dq_sym.append(Symbol("\dot{q_{%d}}"%i, real=True))

        #setup constraints
        self.constraints = []

        #horizontal bonds
        for i in range(n):
            for j in range(m-1):   
                self.constraints.append(distance_constraint(
                    [self.q_sym[3*(i*m+j)+k] for k in range(3)],
                    [self.q_sym[3*(i*m+j+1)+k] for k in range(3)],
                    L
                ))
    
        #vertical bonds
        for j in range(m):
            for i in range(n-1):
                self.constraints.append(distance_constraint(
                    [self.q_sym[3*(i*m+j)+k] for k in range(3)],
                    [self.q_sym[3*((i+1)*m+j)+k] for k in range(3)],
                    L
                ))

        #fixpoints
        for j in range(m):
            self.constraints.append(fixpoint_constraint(
                [self.q_sym[3*j+k] for k in range(3)],
                [0,j*L,0]
            ))

        #jacobian
        self.J = simplify(Matrix(self.constraints).jacobian(self.q_sym))

        #jacobian derivative
        lst = list()
        for constraint in self.constraints:
            value = 0
            for qi,dqi in zip(self.q_sym,self.dq_sym):
                for qj,dqj in zip(self.q_sym,self.dq_sym):
                    value += constraint.diff(qi).diff(qj) * dqi * dqj
            lst.append(value)
                    
        self.JdQ = simplify(Matrix(lst))


    def progress(self, dt):
        dq = self.v*dt

        if 1:
            subs = {sym:val for sym, val in [*zip(self.q_sym,self.q),*zip(self.dq_sym,self.v)]}
            J_ = np.array(self.J.evalf(subs = subs)).astype(np.float64)
            print(J_.round())
            J_T = J_.T
            A = np.matmul(J_,J_T)
            A[np.abs(A)<limit] = 0 
            
            JdQ_ = np.array(list(self.JdQ.evalf(subs=subs))).astype(np.float64)
        
            A_inv = linalg.pinv(A)

            y = np.matmul(A_inv,-JdQ_ - np.matmul(J_,self.Q))

            self.v += (self.Q + np.matmul(J_T,y)) * dt
        
        #self.v += self.Q
        self.q += dq

    def positions(self):
        return self.q.reshape((-1,3)).T #[[x...],[y...],[z...]]
