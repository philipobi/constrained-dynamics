import numpy as np
from scipy.sparse import bsr_array
from scipy.sparse.linalg import lsmr as lsq_solve
from time import sleep

class Constraint:
    pass

class DistanceConstraint(Constraint):
    def __init__(self, a:int, b:int, L:float):
        self.a = a
        self.b = b
        self.L = L

    def construct_jacobian(self,indptr,indices,data,q):
        indptr.append(indptr[-1] + 2)
        indices.append(int(self.a/3))
        indices.append(int(self.b/3))
        
        self.diff = q[self.a:self.a+3] - q[self.b:self.b+3]
        data.append([2*self.diff])
        data.append([-2*self.diff])

    def construct_djdq(self,dq):
        diff = dq[self.a:self.a+3] - dq[self.b:self.b+3]
        return 2*np.sum(np.square(diff))
    
    def eval(self):
        return np.sum(np.square(self.diff))-self.L**2
    
    def info(self,q):
        return f"Distance Constraint between \n{str(q[self.a:self.a+3])} and {str(q[self.b:self.b+3])}\nlength:{self.L}"

class FixpointConstraint(Constraint):
    def __init__(self,a:int,p):
        self.a = a
        self.p = p

    def construct_jacobian(self,indptr,indices,data,q):
        indptr.append(indptr[-1] + 1)
        indices.append(int(self.a/3))
        
        self.diff = q[self.a:self.a+3] - self.p
        data.append([2*self.diff])

    def construct_djdq(self,dq):
        return 2*np.sum(np.square(dq[self.a:self.a+3]))
    
    def eval(self):
        return np.sum(np.square(self.diff))
    
    def info(self,q):
        return f"Point {str(q[self.a:self.a+3])} fixed to {str(self.p)}"
    
    
class Simulation():
    def __init__(self,n=3,m=3,L=5):
        self.m = m
        self.n = n
        self.L = L
        self.ks = 5
        self.kd = 5
        self.null = np.zeros(3*m*n,dtype=np.float64)

        q = []
        
        for i in range(m):
            for j in range(n):
                q.extend([i*L, j*L, 0]) #x,y,z
        
        self.q = np.array(q,dtype=np.float64)

        self.constraints = []

        for i in range(m):
            for j in range(n-1):
                x = i*3*n + j*3
                self.constraints.append(DistanceConstraint(x,x+3,L))

        for j in range(n):
            for i in range(m-1):
                x = i*3*n + j*3
                self.constraints.append(DistanceConstraint(x,x+3*n,L))

        arr = np.array([0,1,0],dtype=np.float64)
        for j in range(n):
            self.constraints.append(FixpointConstraint(j*3,arr*j*L))

        self.dq = np.zeros(3*m*n, dtype=np.float64)
        arr = np.array([[0,0,1]], dtype=np.float64)
        self.f_ext = np.repeat(arr,m*n,0).flatten() * (-9.81)

    def eval(self):
        a,b = (len(self.constraints),3*self.m*self.n)
        
        indptr,indices,data = [[0],[],[]]
        for constraint in self.constraints:
            constraint.construct_jacobian(indptr,indices,data,self.q)

        self.J = bsr_array((data,indices,indptr),shape=(a,b))
        self.J_T = self.J.transpose()
        
        self.dJdq = np.array([constraint.construct_djdq(self.dq) for constraint in self.constraints], dtype=np.float64)
        self.C = np.array([constraint.eval() for constraint in self.constraints], dtype=np.float64)
        self.dC = self.J.dot(self.dq)
        self.left = self.J @ self.J_T
        self.right = -self.dJdq -self.J.dot(self.f_ext) -self.ks*self.C - self.kd*self.dC

    def solve(self):
        self.x, self.istop, itn, normr = lsq_solve(self.left,self.right)[:4]
        self.f_constraint = self.J_T.dot(self.x) if self.istop < 7 else self.null
        print(self.istop)

    def prop(self,dt):
        self.q += self.dq*dt
        self.dq += (self.f_ext+self.f_constraint)*dt

    def run(self,dt):
        self.eval()
        self.solve()
        self.prop(dt)

    
    





    

    
