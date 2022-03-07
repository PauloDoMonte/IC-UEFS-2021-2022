import rebound
import numpy as np
import pandas as pd

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

_2019pg1    = particula('2019pg1',80000, 400/1000,1.037103597335171*149597871,0.03388638947447908,0.140625575916969,281.6394465120655,13.55855529785005,0)
sol         = particula('sol',1.989e+30,696340,0,0,0,0,0,0)
terra       = particula('terra',5.973332e+24,6378.1366,149.60e6,0.01671022,0.00005,-11.26064,102.94719,0)
neo = _2019pg1

anos = 100
passo = 365*24
tempo_simulacao = 60*60*24*365*anos
times = np.linspace(0,tempo_simulacao, anos*passo)

def alteracao(simp):
    sim = simp.contents
    particles = sim.particles

    particles[2].a += 10

def sem_alteracao():

    inicio = []
    db_s = pd.DataFrame({"a":inicio, "e": inicio, "ano": inicio})
    db_s.to_csv("sem_efeito.csv",sep=",",header=True, mode='a')

    sim = rebound.Simulation()
    sim.units = ("kg","km","s")
    sim.integrator = "IAS15"

    sim.add(m=sol.m,r=sol.r)
    sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=0)
    sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=0)

    sim.move_to_com()

    for i, time in enumerate(times):

        sim.integrate(time)
        db_s = db_s.append({"a":sim.particles[2].a, "e":sim.particles[2].e,"ano": i/passo},ignore_index=True)
        print("Sem efeito\tFalta: {}".format(anos*passo - i))

    db_s.to_csv("sem_efeito.csv",sep=",",header=False, mode='a')

def com_alteracao():

    inicio = []
    db_c = pd.DataFrame({"a":inicio, "e": inicio, "ano": inicio})
    db_c.to_csv("com_efeito.csv",sep=",",header=True)

    sim = rebound.Simulation()
    sim.units = ("kg","km","s")
    sim.integrator = "IAS15"

    sim.add(m=sol.m,r=sol.r)
    sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=0)
    sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=0)

    sim.additional_forces=alteracao
    sim.move_to_com()

    for i, time in enumerate(times):

        sim.integrate(time)
        db_c = db_c.append({"semi_eixo com efeito":sim.particles[2].a, "ano": i},ignore_index=True)
        print("Sem Efeito: {}\tAno:{}".format(sim.particles[2].a,i))

    db_c.to_csv("com_efeito.csv",sep=",",header=False, mode='a')

sem_alteracao()
print("Ok")
#com_alteracao()
