import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

#Définition de la fonction qi
def qi (ye, xe, p1, p2 ,p3):
    return ye - (p1/((p2**2 + (xe-p3)**2)))


#Données expérimentales du fichier txt
xe = [0.5, 0.8, 1.1, 1.4, 1.7, 2, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1, 4.4, 4.7, 5]
ye = [1.7777778e-001,2.0325203e-001,2.3228804e-001,2.6455026e-001,2.9895366e-001,3.3333333e-001,3.6429872e-001,3.8759690e-001,3.9920160e-001,3.9682540e-001,3.8095238e-001,3.5460993e-001,3.2206119e-001,2.8735632e-001,2.5348542e-001,2.2222222e-001]

#Définition des paramètres p1,p2,p3 en terme de variables avec Sympy
p1 = sp.symbols("p1")
p2 = sp.symbols("p2")
p3 = sp.symbols("p3")
p = [p1,p2,p3]

#Définition de la fonction Q
Q = 0
for i in range(len (xe)):
    Q = Q + qi(ye[i],xe[i],p1,p2,p3) 

#Définition de la dérivée première de Q telle que décrite dans les consignes du premier devoir
dQ = []
for j in range(3):
    dQj = 0
    for i in range(len(xe)):
        dQj = dQj + qi(ye[i],xe[i],p1,p2,p3)*sp.diff(qi(ye[i],xe[i],p1,p2,p3),p[j])
    dQ.append(dQj)


#QUESTION 2:__________________________________________________________
#Définition des variables x et y
y = sp.symbols("y")
x = sp.symbols("x")
#Dérivées premières de dqi en fonction de p1,p2,p3
print("Trois expressions de derivee du premier ordre: ")
for j in range(len(p)):
    dqi = sp.diff(qi(y,x,p1,p2,p3),p[j])
    print("Expression",j+1,":", dqi)

#Dérivées secondes de dqi p1,p2,p3
print("Six expressions de derivee du deuxieme ordre: ")
i = 0
diffqi = []
for w in range(len(p)):
    dqi = sp.diff(qi(y,x,p1,p2,p3),p[w])
    for e in range(len(p)):
        d2qi = sp.diff(dqi, p[e])
        if d2qi in diffqi:
            pass
        else:  
            i += 1
            diffqi.append(d2qi)
            print("Expression",i,":", d2qi)

#Commencement de la méthode de Newton
F_prime = np.zeros((3,3)) #Mise en place de la matrice F'
F = np.zeros((3,1))       #Mise en place du vecteur F
P_matrix = np.zeros((3,1))#Mise en place du vecteur P
pn = [1.99999996,2.23606795,3.47]
# Construction du vecteur P en fonction des paramètres initiaux
for i in range(3):
    P_matrix[i][0] = pn[i]
P_matrix_1 = np.zeros((3,1)) #Mise en place du vecteur Pn+1
tol = 10e-6 #niveau de tolérance
delta_p1, delta_p2, delta_p3 = 20, 20, 20 #delta_p qui est la différence en pn et pn+1 
d1, d2, d3= [], [], [] 
dt = []
p3_list = [P_matrix[2][0]]
t = 0
#Lorsque dp est supérieur à la tolérance, la méthode de Newton continue ses itérations
while abs(delta_p1) > tol or abs(delta_p2) > tol or abs(delta_p3) > tol:
# while t<10: #Test
    # Construction matrice F 3x1
    for i in range(3):
        F[i][0] = dQ[i].subs({p1:P_matrix[0][0],p2:P_matrix[1][0],p3:P_matrix[2][0]})

    # Construction matrice F' 3x3
    for i in range(3):
        for j in range(3):
            derivation = sp.diff(dQ[i], p[j])
            F_prime[i][j] = derivation.subs({p1:P_matrix[0][0],p2:P_matrix[1][0],p3:P_matrix[2][0]})

    #Inversion de la matrice F'
    F_prime_inv = np.linalg.inv(F_prime)
    # F_prime_inv = np.linalg.solve(F_prime, np.eye(3)) #Test
    #Calcul de delta
    delta = F_prime_inv@F
    #Itération de la méthode de Newton pour Pn+1
    for i in range(3):
        P_matrix_1[i][0] = P_matrix[i][0] - delta[i][0]
    print(P_matrix)
        
        
    #Stockage de delta_p à des fins de représentations graphiques
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

print('La valeur de la fonction erreur est : ', Q.subs({p1:P_matrix_1[0][0], p2:P_matrix_1[1][0], p3:P_matrix_1[2][0]}))
 
plt.plot(d1, label = "dp1")
plt.plot(d2, label = "dp2")
plt.plot(d3, label = "dp3")
plt.title("Valeur de dp en fonction du nombre d'itérations ")
plt.legend()
plt.xlabel('Nombre itérations')
plt.ylabel('dp')
plt.show()




#Question 2 d)__________________________________________________________________________________________________

#Tracer la fonction Q en fonction de p3
p3_list = np.linspace(0,10,1000)
Q_list = []
for i in range(len(p3_list)):
    Q_list.append(Q.subs({p1:1.99999996, p2:2.23606795, p3:p3_list[i]}))
    if Q.subs({p1:1.99999996, p2:2.23606795, p3:p3_list[i]}) <= 0:
        print(p3_list[i]) #À des fins d'analyse


plt.plot(p3_list, Q_list, label = "p3_ini")
plt.legend()
plt.xlabel("p3")
plt.ylabel("Q")
plt.title("Fonction Q versus p3")
plt.yscale("log")
plt.show()
    
