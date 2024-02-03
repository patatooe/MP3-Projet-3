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
#1.26975848e-05] [1.34968235e-05] [4.99996032e-01]
pn = [1.26975848e-05, 1.34968235e-05,(4.99996032e-01)*-1]
# Construciton matrice P 3x1
for i in range(3):
    P_matrix[i][0] = pn[i]

P_matrix_1 = np.zeros((3,1))
tol = 10e-6
print(tol)
delta_p1 = 20
delta_p2 = 20
delta_p3 = 20
d1, d2, d3= [], [], []
dt = []
p3_list = [P_matrix[2][0]]
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
        P_matrix_1[i][0] = P_matrix[i][0] + delta[i][0]
    # print(P_matrix)
        
    delta_p1 = abs(P_matrix_1[0][0] - P_matrix[0][0])
    delta_p2 = abs(P_matrix_1[1][0] - P_matrix[1][0])
    delta_p3 = abs(P_matrix_1[2][0] - P_matrix[2][0])
    
    p3_list.append(P_matrix_1[2][0])
    for i in range(3):
        P_matrix[i][0] = P_matrix_1[i][0]
        
    t += 1 
    d1.append(delta_p1)
    d2.append(delta_p2)
    d3.append(delta_p3)
    
    print("itération", t)
    # print(delta_p1,delta_p2,delta_p3)

 
print(P_matrix_1)

# print('La valeur de la fonction erreur est : ', Q.subs({p1:P_matrix_1[0][0], p2:P_matrix_1[1][0], p3:P_matrix_1[2][0]}))
 
# plt.plot(d1, label = "dp1")
# plt.plot(d2, label = "dp2")
# plt.plot(d3, label = "dp3")
# plt.title("Valeur de dp en fonction du nombre d'itérations ")
# plt.legend()
# plt.xlabel('Nombre itérations')
# plt.ylabel('dp')
# plt.show()





#Question 2 d)__________________________________________________________________________________________________



print(p3_list)
Qp3 = []
for i in range(len(p3_list)):
    Qp3.append(Q.subs({p1:1.26975848e-05, p2:1.34968235e-05, p3:p3_list[i]}))


# En runnant le code plusieurs fois afin de d'acquerir ces listes (flemme de faire une boucle)
# p3_list_0 = [0.0, 0.40533164180453146, 0.4581773868521386, 0.479736282771327, 0.48996031086095027, 0.49500006325506135, 0.49751053855897354, 0.49876466001128394, 0.4993917780636059, 0.49970578435223223, 0.4998639350833191, 0.49994692520160533, 0.49949630165169834, 0.4997582492183961, 0.49989073449453114, 0.49996321102895197, 0.49996242575104344]
# p3_list_00001 = [0.000499996032, 0.4053638879557266, 0.4581888220320365, 0.47974151554736955, 0.48996286104393433, 0.4950013256932094, 0.49751116413415913, 0.49876496799276543, 0.4993919272013181, 0.49970585371025567, 0.49986396356329826, 0.49994692905859933, 0.49948524820473783, 0.49975269094174624, 0.4998878732568432, 0.4999613282613007, 0.4999580100334851]
# p3_list_2 = [0.999992064, 1.0766730124544903, 1.0881756708522419, 1.0937887297865923, 1.0966168617435474, 1.098080171189832, 1.0988654694429252, 1.0993059702627141, 1.0995640554767785, 1.0997209133126382, 1.0998189597278376, 1.0998815009855, 1.0999219650725278, 1.0999484020340977, 1.0999657893916588, 1.0999772762015323, 1.0999848877228615]
# Qp3_0 = []
# for i in range(len(p3_list_0)):
#     Qp3_0.append(Q.subs({p1:1.26975848e-05, p2:1.34968235e-05, p3:p3_list_0[i]}))
# Qp3_00001 = []
# for i in range(len(p3_list_00001)):
#     Qp3_00001.append(Q.subs({p1:1.26975848e-05, p2:1.34968235e-05, p3:p3_list_00001[i]}))
# Qp3_2 = []
# for i in range(len(p3_list_2)):
#     Qp3_2.append(Q.subs({p1:1.26975848e-05, p2:1.34968235e-05, p3:p3_list_2[i]}))




# plt.plot(p3_list_0, Qp3_0, label = "p3_ini = 0")
# plt.legend()
# plt.xlabel("p3")
# plt.ylabel("Q")
# plt.title("Fonction Q versus p3")
# plt.yscale("log")
# plt.show()
    
# plt.plot(p3_list_00001, Qp3_00001, label = "p3_ini = 0.001")
# plt.legend()
# plt.xlabel("p3")
# plt.ylabel("Q")
# plt.title("Fonction Q versus p3")
# plt.yscale("log")
# plt.show()

# plt.plot(p3_list_2, Qp3_2, label = "p3_ini = 2")
# plt.legend()
# plt.xlabel("p3")
# plt.ylabel("Q")
# plt.title("Fonction Q versus p3")
# plt.yscale("log")
# plt.show()
    
