import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def funye (xe, p1, p2, p3):
    return p1/((p2**2 + (xe-p3)**2))

def qi (ye, xe, p1, p2 ,p3):
    return ye - p1/((p2**2 + (xe-p3)**2))

xe = [0.5, 0.8, 1.1, 1.4, 1.7, 2, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5]
ye = [1.7777778e-001,2.0325203e-001,2.3228804e-001,2.6455026e-001,2.9895366e-001,3.3333333e-001,3.6429872e-001,3.8759690e-001,3.9920160e-001,3.9682540e-001,3.8095238e-001,3.5460993e-001,3.2206119e-001,2.8735632e-001,2.5348542e-001,2.2222222e-001]

p1 = sp.symbols("p1")
p2 = sp.symbols("p2")
p3 = sp.symbols("p3")
p = [p1,p2,p3]

# #Q
Q = 0
for i in range(len (xe)):
    Q = Q + qi(ye[i],xe[i],p1,p2,p3) 

print(Q.subs({p1: 1,p2: 1,p3: 1}))


# #dQ/dp 
dQ = []
for j in range(3):
    dQj = 0
    for i in range(len(xe)):
        dQj = dQj + qi(ye[i],xe[i],p1,p2,p3)*sp.diff(qi(ye[i],xe[i],p1,p2,p3),p[j])
    dQ.append(dQj)

# #d2Q/dp2 
# d2Q = []
# for j1 in range(3):
#     d2Qj = 0
#     for j2 in range(3):
#         for i in range(len(xe)):
#             d2Qj = sp.diff(dQ[j1], p[j1], p[j2])
#     d2Q.append(d2Qj)
       


#QUESTION 2:__________________________________________________________
y = sp.symbols("y")
x = sp.symbols("x")
print("Trois expressions de derivee du premier ordre: ")
for j in range(len(p)):
    dqi = sp.diff(qi(y,x,p1,p2,p3),p[j])
    print("Expression",j+1,":", dqi)
    
        
print("Six expressions de derivee du deuxieme ordre: ")
i = 0
diffqi = []
for pj1 in p:
    for pj2 in p:
        d2qi = sp.diff(qi(y,x,p1,p2,p3),pj1, pj2)
        if d2qi in diffqi:
            pass
        else:  
            i += 1
            diffqi.append(d2qi)
            print("Expression",i,":", d2qi)


F_prime = np.zeros((3,3))
F = np.zeros((3,1))
P_matrix = np.zeros((3,1))
pn = [-1,-1,-1]
# Construciton matrice P 3x1
for i in range(3):
    P_matrix[i][0] = pn[i]

P_matrix_1 = np.zeros((3,1))
tol = 10e-6
print(tol)
delta_p1 = 20
delta_p2 = 20
delta_p3 = 20
d = []
dt = []
t = 0
while abs(delta_p1) > tol or abs(delta_p2) > tol or abs(delta_p3) > tol:
#while t<10:
    
    # Construciton matrice F 3x1
    for i in range(3):
        F[i][0] = dQ[i].subs({p1:P_matrix[0][0],p2:P_matrix[1][0],p3:P_matrix[2][0]})

    # Construction matrice F' 3x3
    for i in range(3):
        for j in range(3):
            derivation = sp.diff(dQ[i], p[j])
            F_prime[i][j] = derivation.subs({p1:P_matrix[0][0],p2:P_matrix[1][0],p3:P_matrix[2][0]})

    F_prime_inv = np.linalg.inv(F_prime)
    delta = F_prime_inv@F

    for i in range(3):
        P_matrix_1[i][0] = P_matrix[i][0] - delta[i][0]
    print(P_matrix)
        
    delta_p1 = abs(P_matrix_1[0][0] - P_matrix[0][0])
    delta_p2 = abs(P_matrix_1[1][0] - P_matrix[1][0])
    delta_p3 = abs(P_matrix_1[2][0] - P_matrix[2][0])
    
        
    for i in range(3):
        P_matrix[i][0] = P_matrix_1[i][0]
        
    t += 1 
    d.append(delta_p1)
    print("itération", t)
    # print(delta_p1,delta_p2,delta_p3)

 
plt.plot(d)
plt.show()
    
print(P_matrix_1)