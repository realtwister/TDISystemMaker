import numpy as np
from itertools import combinations, product
from sympy import Matrix, latex

# Settings
M= np.array([[1,-3],[-2,-1],[1,4]]); # From question
d= np.array([0,0,7]);	# From question
x= Matrix([[-1,2],[3,1],[0,0]]);	# The polytope vertices

#Get bounds
bound=np.array([[np.amax(M), np.amin(M)] for i in range(M.shape[1])]);
for n in range(M.shape[0]):	
	for i in combinations(range(M.shape[0]),n+1):
		res=np.sum(M[i,:],axis=0);
		for d in range(M.shape[1]):
			bound[d,:]=[min(bound[d,0],res[d]), max(bound[d,1],res[d])];

# Create L
pos = range(bound[0,0],bound[0,1]+1); 
for d in range(1,M.shape[1]):
	pos = product(pos,range(bound[d,0],bound[d,1]+1));

M=Matrix(M);
L=[];
for p in pos:
	for i in combinations(range(M.shape[0]),M.shape[1]):
		res=((M.T)[:,i]).LUsolve(Matrix(p))
		if(min(res) >= 0 and max(res) <= 1):
			break;
	else:
		continue;
	L.append(p);

# Create the final system
A=Matrix(0,2,[]);
b=Matrix(0,1,[]);
for l in L:
	m=[];
	for i in range(x.shape[0]):
		m.append(Matrix(l).dot(x[i,:]));
	A=A.col_join(Matrix(l).T);
	b=b.col_join(Matrix([max(m)]));
print(latex(A.T)+'\\vec{x}\leq'+latex(b.T));
