import rebound,math
import numpy as np
import pandas as pd
import os, sys
import sqlite3

AU = 1.498e8

def comprimento_orbita(a,e):
    b = a*math.sqrt(1-(e*e))

    h = ((a-b)**2)/((a+b)**2)
    return(math.pi*(a+b)*(1+((3*h)/(10+math.sqrt(4-(3*h))))))

def inicializacao(asteroide):
    if(os.path.exists('{}/'.format(asteroide)) == True):
        if(os.path.exists('{}/raio_hill/'.format(asteroide)) == True):
            if(os.path.exists('{}/raio_hill/passo=ano1'.format(asteroide)) == True):
                pass
            else:
                os.mkdir('{}/raio_hill/passo=ano1'.format(asteroide))
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
            os.mkdir('{}/raio_hill/passo=ano1'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))
    else:
        os.mkdir('{}/'.format(asteroide))
        os.mkdir('{}/raio_hill/'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano1'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))

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

# Lista de particulas disponiveis
_2019pg1 = particula('2019pg1',80000, 400/1000,1.037103597335171*AU,0.03388638947447908,0.140625575916969,281.6394465120655,13.55855529785005,0)
_2021af8 = particula('2021af8',80000, 400/1000,2.018135760648355*AU,0.5151605500121827,9.697802979250167,168.9511312385866,42.56717817033352,0)
_2005vc = particula('2005vc',80000, 1000/1000, 2.082761852899944*AU,0.5947841108529824,4.483281580745502,290.2042173581335,227.9717174732001,0)
_2005cz36 = particula('2005cz36',80000, 1000/1000,2.238066485936266*AU,0.5747761271925266,16.14636501981154,139.3809136995648,116.8081182427147,0)
_29075 = particula('29075',80000,1.3/2,1.698665726802664*AU,0.5077360606964316,12.16734013772536,224.6831644667536,356.6558403536075,0)

sol = particula('sol',1.989e+30,696340,0,0,0,0,0,0)
terra = particula('terra',5.973332e+24,6378.1366,149.60e6,0.01671022,0.00005,-11.26064,102.94719,0)

def distancia_f(neo,f0,ff,passo):

    terra_f = f0
    neo_f = 1
    rh = raio_hills(sol.m,terra.m,terra.a,terra.e)
    inicio = []
    d = pd.DataFrame({"Distancia":inicio,"Distancia/Rh":inicio,"NEO F":inicio,"TERRA F":inicio})

    while(terra_f <= ff):

        sim = rebound.Simulation()
        sim.units = ("kg","km","s")
        sim.integrator = "IAS15"

        sim.add(m=sol.m,r=sol.r)
        sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=terra_f*0.0174533)
        sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=neo_f*0.0174533)

        x_relativo = sim.particles[2].x - sim.particles[1].x
        y_relativo = sim.particles[2].y - sim.particles[1].y
        z_relativo = sim.particles[2].z - sim.particles[1].z

        distancia_neo_terra = math.sqrt(pow(x_relativo,2)+pow(y_relativo,2)+pow(z_relativo,2)) - (neo.r + terra.r)
        print("Distancia: {}\tNEO F: {}\tTERRA F:{}".format(distancia_neo_terra,neo_f,terra_f))

        d = d.append({"Distancia":distancia_neo_terra,"Distancia/Rh":distancia_neo_terra/rh,"NEO F":neo_f,"TERRA F":terra_f},ignore_index=True)

        neo_f += passo
        if(neo_f > 360):
            neo_f = 1
            terra_f += passo
    caminho = "ast={},passo={},f0={},ff={}.csv".format(neo.nome,passo,f0,ff)
    d.to_csv(caminho,sep=",")

#distancia_f(_29075,f0,ff,passo)
#distancia_f(_2005cz36,f0,ff,passo)
