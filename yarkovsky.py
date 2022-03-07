import math, rebound, os
from math import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas

y = 0                               # Ângulo e obliquidade do asteroide bennu
diametro_asteroide = 0.28237*2      # Diametro do asteroide em km
massa_asteroide = 7.329*pow(10, 10) # Massa do asteroide em kg
Tau = 350                           # Inercia térmica do asteroide
F = 136.6                           # Fluxo de radiação solar wm
c = 300000                          # Velocidade da luz km/s
E = 0.9                             # Emissividade infravermelha de rocha entre 0.88 - 0.95
sigma = 5.670400 * 10**-8           # Constante de stefan-boltzmam
n = 1.665454*10**-7                 # Demora 436,65 dias para rodar o sol
w = 4*10**-4                        # Demora 4.29 para rodar a si mesmo
albedo = 0.4                        # Segundo Daniel
alfa = 1 - albedo                   # Parametro Alfa

Phi = (np.pi * diametro_asteroide * diametro_asteroide * F) / (4 * massa_asteroide * c)
t_asterisco_aux = ((alfa * F)/(E*sigma))
t_asterisco = pow(t_asterisco_aux,1/4)

# Efeitos diurno e secular do efeito YARKOVSKY
theta_diurno = Tau * sqrt(w) / (E * sigma * pow(t_asterisco, 3))
theta_secular = Tau * sqrt(n) / (E * sigma * pow(t_asterisco, 3))

b = cos(radians(y))

# Verificando que se o angulo for 90, o cos ser 0 exato
if(y==90): b=0

fa_diurno = (4 * alfa / 9) * (Phi / n) * b * theta_diurno / ( 1+theta_diurno+0.5*(theta_diurno**2) )
fa_secular = -(2 * alfa / 9) * (Phi / n) * (sin( radians(y)) ** 2) * theta_secular / ( 1+theta_secular+0.5*(theta_secular**2) )

resultado = (fa_diurno + fa_secular)

simulacao = "corpos.bin"

anos = 100000
tempo_simulacao = 60*60*24*365*anos
times = np.linspace(0,tempo_simulacao, anos)

def forca_yarkovsky(simp):
    sim = simp.contents
    particles = sim.particles

    particles[2].a += resultado

def sem_efeito():

    inicio = []
    db1 = pd.DataFrame({"semi_eixo sem efeito":inicio, "ano": inicio})

    if(os.path.isfile(simulacao)):
        sim = rebound.Simulation(simulacao)

    else:
        sim = rebound.Simulation()
        sim.units = ("kg","km","s")
        sim.add("Sun")
        sim.add("Earth")
        sim.add("BENNU", m=massa_asteroide)
        sim.save(simulacao)

    sim.integrator = "IAS15"
    sim.move_to_com()

    for i, time in enumerate(times):

        sim.integrate(time)
        db1 = db1.append({"semi_eixo sem efeito":sim.particles[2].a, "ano": i},ignore_index=True)
        print("Sem Efeito: {} ano\tSegundos: {}".format(i,time))

    db1.to_csv("sem_efeito.csv",sep=",",header=False, mode='a')

def com_efeito():

    inicio = []
    db2 = pd.DataFrame({"semi_eixo com efeito":inicio, "ano": inicio})

    if(os.path.isfile(simulacao)):
        sim = rebound.Simulation(simulacao)

    else:
        sim = rebound.Simulation()
        sim.units = ("kg","km","s")
        sim.add("Sun")
        sim.add("Earth")
        sim.add("BENNU", m=massa_asteroide)
        sim.save(simulacao)


    sim.integrator = "IAS15"
    sim.additional_forces = forca_yarkovsky
    sim.force_is_velocity_dependent = False
    sim.move_to_com()

    for i, time in enumerate(times):

        sim.integrate(time)
        db2 = db2.append({"semi_eixo com efeito":sim.particles[2].a, "ano": i},ignore_index=True)
        print("Com Efeito: {} ano\tSegundos: {}".format(i,time))

    db2.to_csv("com_efeito.csv",sep=",",header=False, mode='a')

sem_efeito()
com_efeito()
