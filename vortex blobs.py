"""
Universidad del valle
Escuela de ingenieria mecanica
mecanica de fluidos avanzada

autor: Juan David Contreras, Ruben Dario Aponte
"""

#llamamos librerias a usar
import matplotlib.pyplot as plt
from numpy import *
import math
import numpy as np


# definimos constrantes globales
pi = math.pi  #numero pi

class Cgrid():
    #crear grid para circulo centrado en origen
      def __init__(self, div, dist):  
          self.div = div #numero de viviciones a lo largo del diametro hotizontal del circulo
          self.dist = dist #distrancia ortogonal entre cada nodo o divicion 
          self.radio = self.dist*self.div/2  #radio del circulo
          self.sgmt = linspace(-self.radio, self.radio, div) #arreglo de las diviciones 
          self.r2 = self.radio + self.dist*0.5  #radio de condicion
          
      def N_nodos(self):  #funcion que calcula en numero de nodos
          sgmt, r2=  self.sgmt, self.r2        
          nodos = 0  #conteo inica desde cero
          for i in sgmt: 
              for j in sgmt:
                  L = sqrt(i**2 + j**2)
                  if L < r2:
                      nodos = nodos +1
          return nodos
      def grid(self): 
      #funcion que genera un array con las coordenadas de los puntos del circulo
          n = self.N_nodos()
          hy, r2 =  self.sgmt, self.r2              
          puntos = zeros((n,2))    #crea arreglo para posicionar coordenadas          
          nodo = 0
          for p in hy:                  
              for m in hy:
                  L = sqrt(p**2 + m**2)                    
                  if L < r2:
                      puntos[nodo] = [p , m] #asignar coordenada a cada fila
                      nodo = nodo + 1
          return puntos
          
class XYgrid(Cgrid):  
#generalizacion de la clase Cgrid, para circulo con centro en (x,y)
    def __init__(self, div, dist):
        Cgrid.__init__(self, div, dist)  
        
    def xygrid(self, Cx, Cy): #
        offset = np.array([Cx, Cy]) #crea vertor de desplazamiento
        return Cgrid.grid(self) + offset #desplaza cada punto de el arreglo de coordenadas
                        
          
class VelP():   
#clase que calcula la velocidad inducida por la vorticidad de un punto
    def __init__(self, circulacion):
        self.circulacion = circulacion
       
     
    def vel1(self, Coordenadas, Punto0, Punto1 ):
        circulacion =  self.circulacion
        dx = Coordenadas[Punto0,0] - Coordenadas[Punto1,0]
        dy = Coordenadas[Punto0,1] - Coordenadas[Punto1,1]
        velx =  -circulacion*dy/(2*pi*(dx**2 + dy**2)) #vel en la direccion x
        vely = circulacion*dx/(2*pi*(dx**2 + dy**2))  #vel en la direccion y
        return np.array([velx, vely])
    
class VelT(VelP):

     def __init__(self, circulacion):
         VelP.__init__(self, circulacion)
     
     def vel(self, Coordenadas, Punto, V0=np.array([0,0])):
          # velocidad inducida por una coleccion de puntos a una paeticula
         N = len(Coordenadas) #numero de filas del arreglo de coordenadas
         V = V0               #condcion inicial
         for i in range(0,N):
             if i != Punto:                 
                 V = V + VelP.vel1(self, Coordenadas, Punto, i)
         return V

     def velT(self, Coordenadas, V0=np.array([0,0])):
         #calcular la velocidad de todos las particos inducida por la vorticidad de las demas
        NN = len(Coordenadas)
        VT = zeros((NN,2))
        for j in range(0,NN):
            VT[j] = self.vel(Coordenadas, j, V0)
        return VT
            
class Euler():
     def __init__(self, pasos, x0, m, dt):
         self.pasos, self.x0, self.m, self.dt = pasos, x0, m, dt
         
     def euler1(self):
         #metodod para calcular la paroximacion de un solo paso dt
         X1 = self.x0 + self.m*self.dt
         return X1
         
     def euler(self,  circulacion):
         pasos, x0, m, dt = self.pasos, self.x0, self.m, self.dt
         N = len(x0)
         Cor = x0 #coordenadas iniciales
         V = VelT(circulacion)
         for i in range(0,pasos):       #aplica el la apliximacion i veces      
             for j in range(0, N):      #aplica la aproximacion para cada punto (x,y)
                 Cor[j] = Cor[j] + m[j]*dt
             m = V.velT(Cor)   #recalcula la velocidad en cada nuevo punto para las nuevas corrdenadas
         return Cor   
         
'''
##############################################################################
parametros del problema
##############################################################################
'''

D = 0.51          #diametro del los parches circulares de velocidad
area = pi*(D/2)**2     #area de cada parcher
w = -1                 #vorticidad de cada parche
G = w*area              #Circulacion de cada parche
tiempo = 1              #tiempo total de simulacion
dt = 0.2               #paso de tiempo para la aproximacion
d = 0.4
'''
###############################################################################
Posiciones iniciales de los blood de los parches de velocidad
definicio de geometria del problema
##############################################################################
'''
div = 21              #numero de diviciones para el orden de los 350 bloobs por parche
dis = D/div
#definimos los puntos del ciculo             
Circulo = XYgrid(div, dis) 
    #Puntos del circulo centrado en cero
N_n = Circulo.N_nodos()      
C1 = Circulo.xygrid(0 , d)    #puntos del circulo centrado en (0,4)
C2 = Circulo.xygrid(0 , -d)
C3 = Circulo.xygrid(d, 0)
C4 = Circulo.xygrid(-d, 0)   #opcionales

########################################################
#unir los dos arreglos de coordenadas 
#en uno solo con los puntos de todos los parches
########################################################

C = concatenate((C1,C2)) 


'''
################################################################
grafica de la posicion inicial de los bloods
################################################################
'''


fig1 = plt.figure(1)
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
ax1.axis('equal')
ax1.scatter(C1.T[0], C1.T[1], s=8, color='red')
ax1.scatter(C2.T[0], C2.T[1], s=8, color='red')
ax1.scatter(C3.T[0], C3.T[1], s=8, color='red')
ax1.scatter(C4.T[0], C4.T[1], s=8, color='red')
#ax1.show()



'''
###########################################################################
calcular la velocidad inducida por la vorticidad de los bloobs
############################################################################
'''

g=G/(N_n)

#define un objero VelT --vel total-- con diametro 0.01 y vorticidad -2
V = VelT(g)
vT = V.velT(C)   #calcula la velocidad de cada punto de los 2 circulos

#se calcula el desplazamiento con el metodo de Euler
pasos = int(tiempo/dt)  #numero de iteracion que utilizaca la aproximacion 
#objeto para usar el metodo de euler para el calculo aproximado del desplaxamiento de las particulas de vorticidad
move = Euler(pasos, C, vT, dt)
move1= move.euler(g) #posicion de las particulas en el tiempo t
#
'''
#########################################################################
grafica de la posicion de las perticulas
#########################################################################
'''
fig2 = plt.figure(2)
ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
ax2.axis('equal') #igualar la escala de los ejes
ax2.scatter(move1.T[0], move1.T[1], s=5, color='blue')  # el .T indica transpuesta
ax2.show()



    
    
    
          
                          
                         