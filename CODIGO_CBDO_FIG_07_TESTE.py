import rebound
from math import *
import numpy as np
import random
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1
import matplotlib.pyplot as plt2
from mpl_toolkits.mplot3d import Axes3D
import statistics
import os.path
import os
from decimal import Decimal

random.seed(5)


def modulo(a, b, c):
    return sqrt( a*a + b*b + c*c )




massa_asteroide = 7.329*pow(10, 10) # kg
raio_asteroide = 0.28237  #km
p_asteroide = 1190   # kg/m3

Noutputs = 10000
size = 4000
x_SEM_EFEITO = np.zeros((2, Noutputs))
y_SEM_EFEITO = np.zeros((2, Noutputs))
z_SEM_EFEITO = np.zeros((2, Noutputs))
a_SEM_EFEITO = np.zeros((1, Noutputs))

x_COM_EFEITO = np.zeros((2, Noutputs))
y_COM_EFEITO = np.zeros((2, Noutputs))
z_COM_EFEITO = np.zeros((2, Noutputs))
a_COM_EFEITO = np.zeros((1, Noutputs))

t_SIMULACAO = np.zeros((1, Noutputs))

year = 2. * np.pi  # One year in units where G=1
tempo_simulacao = 100 # anos
times = np.linspace(0., tempo_simulacao * year, Noutputs)

resultado = 0




def forca_yarkovsky(simp):

    sim = simp.contents
    ps = sim.particles

    #D = 1.88752687e-9 * 2 # UA 0.28237*2
    D = 0.28237*2
    #y = 178
    m = 7.329*pow(10, 10) # kg
    Tau = 350

    F = 136.3 # wm-2  Fluxo da Radiação Solar
    #c = 173.1446/(24*60*60)   # UA/s  300000 # km/s
    c = 300000 # km/s
    E = 0.9 # emissividade infravermelha de rocha entre 0.88 - 0.95      OK
    sigma = 5.670400 * 10**-8  # constante de stefan-boltzmam  OK
    n = 1.665454*10**-7     # demora 436,65 dias para rodar o sol
    w = 4*10**-4  #  demora 4.29 para rodar a si mesmo
    #p = 1.3  # g/cm3
    A = 0.4 #4.4 # segundo o site.
    alfa = 1.0 - A   # OK

    Phi = pi * D * D * F / (4 * m * c)

    T_asterisco = 0
    aux =  ((alfa * F) / (E * sigma))
    if aux < 0:
        aux = abs(aux)
        T_asterisco = pow(aux,1/4)*(-1)
    else:
        T_asterisco = pow(aux,1/4)

    Theta_diurno = Tau * sqrt(w) / (E * sigma * pow(T_asterisco, 3))
    Theta_secular = Tau * sqrt(n) / (E * sigma * pow(T_asterisco, 3))

    b = cos(radians(y))
    if(y==90):
        b=0

    fa_diurno = (4 * alfa / 9) * (Phi / n) * b * Theta_diurno / ( 1+Theta_diurno+0.5*(Theta_diurno**2) )
    fa_secular = -(2 * alfa / 9) * (Phi / n) * (sin( radians(y)) ** 2) * Theta_secular / ( 1+Theta_secular+0.5*(Theta_secular**2) )

    resultado = (fa_diurno + fa_secular)
    #print( resultado )

    ps[2].a += resultado
    #return resultado



def simulacao_SEM_EFEITO_():

    print("SEM Y = " + str(y) + "º")
    fn = "CorposSistema.bin"  # file to store the simulation (NASA Horizon queries are slow)
    menorD = 9999999999999999
    menorT = 0


    if os.path.isfile(fn):
        sim = rebound.Simulation(fn)
    else:
        sim = rebound.Simulation()
        sim.units = ("kg", "km", "yr")
        sim.add("Sun")
        sim.add("Earth")
        sim.add("BENNU",  m=massa_asteroide)  # ESTÁ NO GRUPO DOS ASTEROIDES APOLLO  (Bennu é o alvo da missão OSIRIS-REx )
        sim.save(fn)

    sim.dt = 0.1
    sim.integrator = "IAS15"  # IAS15 is the default integrator, so we actually don't need this line
    sim.move_to_com()  # We always move to the center of momentum frame before an integration
    ps = sim.particles  # ps is now an array of pointers and will change as the simulation runs

    for i, time in enumerate(times):
        sim.integrate(time)
        #print(i)

        #if(i >= 6000):
            #print(i)
        ps[2].r = raio_asteroide
        x_SEM_EFEITO[0][i] = ps[1].x  # This stores the data which allows us to plot it later
        y_SEM_EFEITO[0][i] = ps[1].y
        z_SEM_EFEITO[0][i] = ps[1].z
        x_SEM_EFEITO[1][i] = ps[2].x
        y_SEM_EFEITO[1][i] = ps[2].y
        z_SEM_EFEITO[1][i] = ps[2].z
        a_SEM_EFEITO[0][i] = ps[2].a
        t_SIMULACAO[0][i] = time

        aux2 = np.sqrt(
            np.square(ps[1].x - ps[2].x) + np.square(ps[1].y - ps[2].y) + np.square(ps[1].z - ps[2].z))
        if (aux2 < menorD):
            menorD = aux2
            menorT = time

    print("Menor D sem efeito -->" + str(menorD) + "    TEMPO -->> " + str(menorT/year))


def simulacao_COM_EFEITO_():

    print("COM Y = "+str(y)+"º")

    menorD = 9999999999999999
    menorT = 0

    fn = "CorposSistema.bin"  # file to store the simulation (NASA Horizon queries are slow)

    if os.path.isfile(fn):
        sim1 = rebound.Simulation(fn)
    else:
        sim1 = rebound.Simulation()
        sim1.units = ("kg", "km", "yr")
        sim1.add("Sun")
        sim1.add("Earth")
        sim1.add("BENNU",  m=massa_asteroide, r=raio_asteroide )  # ESTÁ NO GRUPO DOS ASTEROIDES APOLLO  (Bennu é o alvo da missão OSIRIS-REx )
        sim1.save(fn)

    sim1.dt = 0.1
    sim1.integrator = "IAS15"  # IAS15 is the default integrator, so we actually don't need this line
    sim1.additional_forces = forca_yarkovsky
    sim1.force_is_velocity_dependent = False
    sim1.move_to_com()  # We always move to the center of momentum frame before an integration
    ps1 = sim1.particles  # ps is now an array of pointers and will change as the simulation runs

    for i, time in enumerate(times):
        sim1.integrate(time)

        #if (i >= 6000):
            #print(i)
        x_COM_EFEITO[0][i] = ps1[1].x  # This stores the data which allows us to plot it later
        y_COM_EFEITO[0][i] = ps1[1].y
        z_COM_EFEITO[0][i] = ps1[1].z
        x_COM_EFEITO[1][i] = ps1[2].x
        y_COM_EFEITO[1][i] = ps1[2].y
        z_COM_EFEITO[1][i] = ps1[2].z
        a_COM_EFEITO[0][i] = ps1[2].a
        t_SIMULACAO[0][i] = time

        aux2 = np.sqrt( np.square(ps1[1].x - ps1[2].x) + np.square(ps1[1].y - ps1[2].y) + np.square( ps1[1].z - ps1[2].z ))
        if(aux2 < menorD):
            menorD = aux2
            menorT = time

    print( "Menor D com efeito -->"+ str(menorD) + "    TEMPO -->> " + str(menorT/year) )



y_vetor = [178]
y = 0;

fig = plt.figure("Figura 01")
ax = plt.subplot(111)


for i in range(0, len(y_vetor)):

    y = y_vetor[i]
    simulacao_SEM_EFEITO_()
    simulacao_COM_EFEITO_()
    #print(resultado)

    aux1 = a_SEM_EFEITO[0]  #np.sqrt(np.square(x_SEM_EFEITO[0] - x_SEM_EFEITO[1]) + np.square(y_SEM_EFEITO[0] - y_SEM_EFEITO[1]) + np.square(z_SEM_EFEITO[0] - z_SEM_EFEITO[1]))
    aux2 = a_COM_EFEITO[0] #np.sqrt(np.square(x_COM_EFEITO[0] - x_COM_EFEITO[1]) + np.square(y_COM_EFEITO[0] - y_COM_EFEITO[1]) + np.square(z_COM_EFEITO[0] - z_COM_EFEITO[1]))
    #diferenca = aux2 - aux1
    diferenca = a_COM_EFEITO[0] - a_SEM_EFEITO[0]
    #print(diferenca)
    ax.set_xlabel("TIME [yrs]", fontsize=18)
    ax.set_ylabel("RELATIVE DISTANCE [Km]", fontsize=18)
    plt.yscale("log")
    #plt.xscale("log")
    ax.set_title("RELATIVE DISTANCE BETWEEN EARTH AND BENNU ASTEROID WITH AND WITHOUT THE YARKOVSKY EFFECT", fontsize=18)
    #ax.plot(t_SIMULACAO[0]/year , aux1   , label="SEM y = " + str(y) + "º");
    ax.plot(t_SIMULACAO[0]/year , aux1 ,  label="WITHOUT YARKOVSKY EFFECT");
    ax.plot(t_SIMULACAO[0]/year, aux2, label="WITH YARKOVSKY EFFECT");
    #ax.plot(t_SIMULACAO[0] / year, diferenca, label="y = " + str(y) + "º");


plt.rcParams.update({'font.size': 18})
ax.legend(loc="best")
plt.tight_layout()
plt.show()



































