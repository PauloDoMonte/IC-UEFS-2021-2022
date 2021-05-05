import rebound,sqlite3
import numpy as np
import pandas as pd
import os, sys
import random
import math

def inicializacao(asteroide):
    if(os.path.exists('{}/'.format(asteroide)) == True):
        if(os.path.exists('{}/raio_hill/'.format(asteroide)) == True):
            if(os.path.exists('{}/raio_hill/passo=1ano'.format(asteroide)) == True):
                pass
            else:
                os.mkdir('{}/raio_hill/passo=1ano'.format(asteroide))
            if(os.path.exists('{}/raio_hill/passo=ano10'.format(asteroide)) == True):
                pass
            else:
                os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
            if(os.path.exists('{}/raio_hill/passo=ano100'.format(asteroide)) == True):
                pass
            else:
                os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))
        else:
            os.mkdir('{}/raio_hill/'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=1ano'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))
    else:
        os.mkdir('{}/'.format(asteroide))
        os.mkdir('{}/raio_hill/'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=1ano'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))

AU = 1.498e8

def raio_hills(m1,m2,a,e):
    return(a*(1-e)*(pow(m2/(3*m1),1/3)))

class particula:
    def __init__(self,nome,massa,raio,a,e,inc,omega,Omega,f):
        self.nome = nome
        self.m = massa
        self.r = raio
        self.a = a
        self.e = e
        self.inc = inc*0.0174533
        self.omega = omega*0.0174533    # argumento do perigeu
        self.Omega = Omega*0.0174533    # Longitude do no ascendente
        self.f = f*0.0174533

sol = particula('sol',1.989e+30,696340,0,0,0,0,0,0)
terra = particula('terra',5.973332e+24,6378.1366,149.60e6,0.01671022,0.00005,-11.26064,102.94719,0)

while True:
    tmax = random.uniform(10000,1000000)
    passo = 1
    massa_neo = random.uniform(1,10000)
    raio_neo = random.uniform(0.001,100)
    
    delta_r = random.uniform(1,10)
    delta_v = random.uniform(11,20)

    sim = rebound.Simulation()
    sim.units = ("kg","km","yr")
    sim.integrator = "IAS15"

    sim.add(m=sol.m,r=sol.r)
    sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=terra.f*0.0174533)
    sim.add(m=massa_neo,r=raio_neo,x=sim.particles[1].x+delta_r,y=sim.particles[1].y+delta_r,z=sim.particles[1].z+delta_r,vx=sim.particles[1].vx+delta_v,vy=sim.particles[1].vy+delta_v,vz=sim.particles[1].vz+delta_v)

    tempo_maximo = tmax
    incremento_tempo = passo

    tamanho = int(tempo_maximo/(incremento_tempo))
    times = np.linspace(-tempo_maximo,0,tamanho)
    
    # Elementos keplerianos
    f1      = np.zeros(int(tmax*passo))
    a2      = np.zeros(int(tmax*passo))
    e2      = np.zeros(int(tmax*passo))
    inc2    = np.zeros(int(tmax*passo))
    omega2  = np.zeros(int(tmax*passo))
    Omega2  = np.zeros(int(tmax*passo))
    f2      = np.zeros(int(tmax*passo))

    for i,time in enumerate(times):
        sim.move_to_com()
        sim.integrate(time)

        # Elementos keplerianos
        f1[i]       = sim.particles[1].f
        a2[i]       = sim.particles[2].a
        e2[i]       = sim.particles[2].e
        inc2[i]     = sim.particles[2].inc
        omega2[i]   = sim.particles[2].omega
        Omega2[i]   = sim.particles[2].Omega
        f2[i]       = sim.particles[2].f

    conn = sqlite3.connect('raio_hill.db')
    cursor = conn.cursor()
    cursor.execute("""INSERT INTO RH (f1,m2,r2,a2,e2,inc2,om2,O2,f2,tmax) VALUES (?,?,?,?,?,?,?,?,?,?)""",
    (f1[i],massa_neo,raio_neo,a2[i],e2[i],inc2[i],omega2[i],Omega2[i],f2[i],tmax))
    conn.commit()
    cursor.close()
    conn.close()

    print(tmax)
