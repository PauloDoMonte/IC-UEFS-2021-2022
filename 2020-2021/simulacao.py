import rebound,math,sqlite3
import numpy as np
import pandas as pd
import os, sys

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

inicializacao('2005vc')
inicializacao('2005cz36')

if(sys.argv[1] == '-h'):
    print("python3 simulacao.py [ASTEROIDE] [TEMPO MAXIMO] [INCREMENTO TEMPO]")
    sys.exit()

AU = 1.498e8

def raio_hills(m1,m2,a,e):
    return(a*(1-e)*(pow(m2/(3*m1),1/3)))

class particula:
    def __init__(self, massa,raio,a,e,inc,omega,Omega,f):
        self.m = massa
        self.r = raio
        self.a = a
        self.e = e
        self.inc = inc*0.0174533
        self.omega = omega*0.0174533    # argumento do perigeu
        self.Omega = Omega*0.0174533    # Longitude do nó ascendente
        self.f = f*0.0174533

_2019pg1 = particula(80000, 400/1000,1.037103597335171*AU,0.03388638947447908,0.140625575916969,281.6394465120655,13.55855529785005,0)
_2021af8 = particula(80000, 400/1000,2.018135760648355*AU,0.5151605500121827,9.697802979250167,168.9511312385866,42.56717817033352,0)
_2005vc = particula(80000, 1000/1000, 2.082761852899944*AU,0.5947841108529824,4.483281580745502,290.2042173581335,227.9717174732001,0)
_2005cz36 = particula(80000, 1000/1000,2.238066485936266*AU,0.5747761271925266,16.14636501981154,139.3809136995648,116.8081182427147,0)


sol = particula(1.989e+30,696340,0,0,0,0,0,0)
terra = particula(5.973332e+24,6378.1366,149.60e6,0.01671022,0.00005,-11.26064,102.94719,0)

if(sys.argv[1] == '2005vc'):
    neo = _2005vc
elif(sys.argv[1] == '2005cz36'):
    neo = _2005cz36

minuto = 60
hora = minuto*60
dia = hora*24
ano = dia*365

tempo_maximo = ano*int(sys.argv[2])
incremento_tempo = ano/int(sys.argv[3])

tamanho = int(tempo_maximo/(incremento_tempo))
times = np.linspace(0,tempo_maximo,tamanho)

quantidade = 360
angulos = np.linspace(0,quantidade,quantidade)

indice = 0
anos = np.zeros(quantidade*quantidade)
f_terra = np.zeros(quantidade*quantidade)
f_neo = np.zeros(quantidade*quantidade)
n = int(sys.argv[4])   # mutiplicador do raio de hill

if(neo == _2005vc):
    arquivo = '2005vc/raio_hill/passo=ano{}/2005vc-{}x.csv'.format(int(sys.argv[3]),n)
elif(neo == _2005cz36):
    arquivo = '2005cz36/raio_hill/passo=ano{}/2005cz36-{}x.csv'.format(int(sys.argv[3]),n)

for neo.f in angulos:
    for terra.f in angulos:

        sim = rebound.Simulation()
        sim.units = ("kg","km","s")
        sim.integrator = "IAS15"

        sim.add(m=sol.m,r=sol.r)
        sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=terra.f*0.0174533)
        sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=neo.f*0.0174533)

        for i,time in enumerate(times):

            sim.move_to_com()
            sim.integrate(time)

            x_relativo = sim.particles[2].x - sim.particles[1].x
            y_relativo = sim.particles[2].y - sim.particles[1].y
            z_relativo = sim.particles[2].z - sim.particles[1].z

            h_terra = n*raio_hills(sol.m,terra.m,sim.particles[1].a,sim.particles[1].e)
            distancia_neo_terra = math.sqrt(pow(x_relativo,2)+pow(y_relativo,2)+pow(z_relativo,2)) - (neo.r + terra.r)

            if(distancia_neo_terra <= h_terra):
                f_neo[indice] = neo.f
                f_terra[indice] = terra.f
                anos[indice] = time/ano

                indice = indice + 1

                vector = pd.DataFrame({"F NEO":f_neo,"F TERRA":f_terra,"ANO":anos})
                vector.to_csv(arquivo,sep=",")

                break

            if(i >= (len(times)-1)):
                f_neo[indice] = neo.f
                f_terra[indice] = terra.f
                anos[indice] = 99999

                indice = indice + 1

                vector = pd.DataFrame({"F NEO":f_neo,"F TERRA":f_terra,"ANO":anos})
                vector.to_csv(arquivo,sep=',')


        print("NEO:{}\tTERRA:{}".format(neo.f,terra.f))
