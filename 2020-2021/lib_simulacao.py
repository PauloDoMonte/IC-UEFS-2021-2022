import rebound,math
import numpy as np
import pandas as pd
import os, sys
import seaborn as se
import matplotlib.pyplot as plt

AU = 1.498e8
e=math.e

def comprimento_orbita(a,e):
    b = a*math.sqrt(1-(e*e))

    h = ((a-b)**2)/((a+b)**2)
    return(math.pi*(a+b)*(1+((3*h)/(10+math.sqrt(4-(3*h))))))

def inicializacao(asteroide):
    if(os.path.exists('{}/'.format(asteroide)) == True):
        if(os.path.exists('{}/raio_hill/'.format(asteroide)) == True):
            if(os.path.exists('{}/raio_hill/passo=ano1'.format(asteroide)) != True):
                os.mkdir('{}/raio_hill/passo=ano1'.format(asteroide))
            if(os.path.exists('{}/raio_hill/passo=ano10'.format(asteroide)) != True):
                os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
            if(os.path.exists('{}/raio_hill/passo=ano100'.format(asteroide)) != True):
                os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))
        else:
            os.mkdir('{}/raio_hill/'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=ano1'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
            os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))

        if(os.path.exists('{}/estatico/'.format(asteroide)) != True):
            os.mkdir('{}/estatico'.format(asteroide))
            os.mkdir('{}/estatico/passo=1'.format(asteroide))
            os.mkdir('{}/estatico/passo=0.1'.format(asteroide))

        if(os.path.exists('{}/graficos/'.format(asteroide)) != True):
            os.mkdir('{}/graficos/'.format(asteroide))

        if(os.path.exists('{}/colisao/'.format(asteroide)) != True):
            os.mkdir('{}/colisao/'.format(asteroide))
            os.mkdir('{}/colisao/parte1/'.format(asteroide))
            os.mkdir('{}/colisao/parte2/'.format(asteroide))

    else:
        os.mkdir('{}/'.format(asteroide))
        os.mkdir('{}/raio_hill/'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano1'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano10'.format(asteroide))
        os.mkdir('{}/raio_hill/passo=ano100'.format(asteroide))
        os.mkdir('{}/estatico/'.format(asteroide))
        os.mkdir('{}/estatico/passo=1'.format(asteroide))
        os.mkdir('{}/estatico/passo=0.1'.format(asteroide))
        os.mkdir('{}/graficos/'.format(asteroide))
        os.mkdir('{}/colisao/'.format(asteroide))
        os.mkdir('{}/colisao/parte1/'.format(asteroide))
        os.mkdir('{}/colisao/parte2/'.format(asteroide))

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

# Particulas padrões
sol = particula('sol',1.989e+30,696340,0,0,0,0,0,0)
terra = particula('terra',5.973332e+24,6378.1366,149.60e6,0.01671022,0.00005,-11.26064,102.94719,0)

rh = raio_hills(sol.m,terra.m,terra.a,terra.e)

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
    caminho = "{}/estatico/passo={}/ast={},passo={},f0={},ff={}.csv".format(neo.nome,passo,neo.nome,passo,f0,ff)
    d.to_csv(caminho,sep=",")

def grafico_df(neo,passo):

    caminho = "{}/estatico/passo={}/ast={},passo={},f0={},ff={}.csv".format(neo,passo,neo,passo,0,360)
    arquivo = pd.read_csv(caminho)

    dx = arquivo.pivot(index="NEO F",columns="TERRA F",values="Distancia/Rh")

    fig = plt.figure(figsize=(15,10))
    se.heatmap(dx,cbar_kws={'label':'Distância/Raio de Hill da Terra'})
    plt.title("Distancia/Raio de hill da Terra em função das anomalias verdadeiras")
    plt.ylabel("Anomalia verdadeira do {}".format(neo))
    plt.xlabel("Anomalia verdadeira da Terra")
    plt.grid(color='green', linestyle='--', linewidth=0.5)
    fig.savefig('{}/graficos/distanciarhxanomaliasverdadeiras.png'.format(neo))

def colisao_neo_satelite_parte1(neo,passo,tempo,vex,vey,vez,gamma,chi):

    def propulsao(simp):
        sim = simp.contents
        particles = sim.particles

        particles[3].ax += ((-vex*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
        particles[3].ay += ((-vey*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
        particles[3].az += ((-vez*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))

    terra_f = 1
    neo_f = 1

    tempo_simulacao = tempo#60*60*1

    inicio = []
    d = pd.DataFrame({"NEO F":inicio,"TERRA F":inicio,"X0":inicio,"Y0":inicio,"Z0":inicio,"X":inicio,"Y":inicio,"Z":inicio,"VX0":inicio,"VY0":inicio,"VZ0":inicio,"VX":inicio,"VY":inicio,"VZ":inicio,"VC":inicio})
    caminho = "{}/colisao/parte1/t={},vex={},vey={},vez={},gamma={},chi={}.csv".format(neo.nome,tempo_simulacao,vex,vey,vez,gamma,chi)
    d.to_csv(caminho,sep=",")

    while(terra_f <= 360):
        for v in [-1e-1/math.sqrt(3),-1/math.sqrt(3),1/math.sqrt(3),1e1/math.sqrt(3)]:
            sim = rebound.Simulation()
            sim.units=("kg","km","s")
            sim.integrator = "IAS15"

            sim.additional_forces = propulsao
            sim.force_is_velocity_dependent = False

            sim.add(m=sol.m,r=sol.r)
            sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=terra_f*0.0174533)
            sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=neo_f*0.0174533)
            sim.add(m=1000,r=0.01,x=sim.particles[2].x+1e-3,y=sim.particles[2].y,z=sim.particles[2].z,vx=sim.particles[2].vx+v,vy=sim.particles[2].vy+v,vz=sim.particles[2].vz+v)

            r0 = [sim.particles[3].x,sim.particles[3].y,sim.particles[3].z]
            v0 = [sim.particles[3].vx,sim.particles[3].vy,sim.particles[3].vz]

            sim.move_to_com()
            sim.integrate(-tempo_simulacao)

            rf = [sim.particles[3].x,sim.particles[3].y,sim.particles[3].z]
            vf = [sim.particles[3].vx,sim.particles[3].vy,sim.particles[3].vz]

            d = d.append({"NEO F":neo_f,"TERRA F":terra_f,"X0":r0[0],"Y0":r0[1],"Z0":r0[2],"X":rf[0],"Y":rf[1],"Z":rf[2],"VX0":v0[0],"VY0":v0[1],"VZ0":v0[2],"VX":vf[0],"VY":vf[1],"VZ":vf[2],"VC":v*math.sqrt(3)},ignore_index=True)

        print("TERRA:{}\tNEO:{}".format(terra_f,neo_f))

        neo_f += passo
        if(neo_f > 360):
            neo_f = 1
            terra_f += passo

        if(len(d) > 200):
            d.to_csv(caminho,sep=",",header=False,mode='a')
            d = pd.DataFrame({"NEO F":inicio,"TERRA F":inicio,"X0":inicio,"Y0":inicio,"Z0":inicio,"X":inicio,"Y":inicio,"Z":inicio,"VX0":inicio,"VY0":inicio,"VZ0":inicio,"VX":inicio,"VY":inicio,"VZ":inicio,"VC":inicio})

    d.to_csv(caminho,sep=",",header=False,mode='a')

def colisao_neo_satelite_parte2(neo,passo,tempo,vex,vey,vez,gamma,chi):

    sim = rebound.Simulation()
    sim.units=("kg","km","s")
    sim.integrator = "IAS15"

    sim.add(m=sol.m,r=sol.r)
    sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=0*0.0174533)
    sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=0*0.0174533)

    distancia_inicial_neo = np.zeros(360)
    distancia_inicial_terra = np.zeros(360)
    distancia_terra_neo = np.zeros(360*360)

    l0_terra = sim.particles[1].a*(1-sim.particles[1].e**2)
    l0_neo = sim.particles[2].a*(1-sim.particles[2].e**2)

    for f in range(0,360):
        distancia_inicial_neo[f] = l0_neo/(1+(sim.particles[2].e*math.cos(f)))
        distancia_inicial_terra[f] = l0_terra/(1+(sim.particles[1].e*math.cos(f)))

    inicio = []
    distancias = pd.DataFrame({"NEO F":inicio,"TERRA F":inicio,"R_RELATIVO":inicio})
    distancias_caminho = "{}/colisao/parte2/t={},vex={},vey={},vez={},gamma={},chi={}_distancias_iniciais.csv".format(neo.nome,tempo,vex,vey,vez,gamma,chi)
    distancias.to_csv(distancias_caminho,sep=",")

    for neo_f in range(0,360):
        for terra_f in range(0,360):
            distancia_inicial_terra_neo = abs(distancia_inicial_neo[neo_f] - distancia_inicial_terra[terra_f])
            distancias = distancias.append({"NEO F":neo_f,"TERRA F":terra_f,"R_RELATIVO":distancia_inicial_terra_neo},ignore_index=True)

        if(len(distancias)>100):
            distancias.to_csv(distancias_caminho,sep=",",header=False,mode='a')
            distancias = pd.DataFrame({"NEO F":inicio,"TERRA F":inicio,"R_RELATIVO":inicio})

        print("PARTE 1 Falta: {}".format(360-neo_f))
    distancias.to_csv(distancias_caminho,sep=",",header=False,mode='a')

    lista = ["NEO F","TERRA F","VX","VY","VZ","VC"]
    caminho_parte1 = "{}/colisao/parte1/t={},vex={},vey={},vez={},gamma={},chi={}.csv".format(neo.nome,tempo,vex,vey,vez,gamma,chi)
    df = pd.read_csv(caminho_parte1, usecols=lista)

    distancias = pd.DataFrame({"NEO F":inicio,"TERRA F":inicio,"MIN_R":inicio,"MAX_R":inicio,"VC":inicio})
    distancias_caminho = "{}/colisao/parte2/t={},vex={},vey={},vez={},gamma={},chi={}_distancias_finais.csv".format(neo.nome,tempo,vex,vey,vez,gamma,chi)
    distancias.to_csv(distancias_caminho,sep=",")

    for i in range(0,len(df["NEO F"])):

        sim = rebound.Simulation()
        sim.units=("kg","km","s")
        sim.integrator = "IAS15"

        terra_f = df["TERRA F"][i]
        neo_f = df["NEO F"][i]

        sim.add(m=sol.m,r=sol.r)
        sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=terra_f*0.0174533)
        sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=neo_f*0.0174533)

        sim.particles[2].vx = (1000*df["VX"][i] + sim.particles[2].vx*neo.m)/(1000+neo.m)
        sim.particles[2].vy = (1000*df["VY"][i] + sim.particles[2].vy*neo.m)/(1000+neo.m)
        sim.particles[2].vz = (1000*df["VZ"][i] + sim.particles[2].vz*neo.m)/(1000+neo.m)

        distancia_neo = np.zeros(360)
        distancia_terra = np.zeros(360)
        distancia_terra_neo = np.zeros(360*360)
        l_neo = sim.particles[2].a*(1-sim.particles[2].e**2)
        l_terra = sim.particles[1].a*(1-sim.particles[1].e**2)

        aux=0
        for neo_f in range(0,360):
            distancia_neo[neo_f] = l_neo/(1+(sim.particles[2].e*math.cos(neo_f)))
            for terra_f in range(0,360):
                distancia_terra[terra_f] = l_terra/(1+(sim.particles[1].e*math.cos(terra_f)))
                distancia_terra_neo[aux] = abs(distancia_neo[neo_f] - distancia_terra[terra_f])
                aux+=1
                print(aux)

        distancias = distancias.append({"NEO F":df["NEO F"][i],"TERRA F":df["TERRA F"][i],"MIN_R":min(distancia_terra_neo),"MAX_R":max(distancia_terra_neo),"VC":df["VC"][i]},ignore_index=True)

        distancias.to_csv(distancias_caminho,sep=",",header=False,mode='a')
        distancias = pd.DataFrame({"NEO F":inicio,"TERRA F":inicio,"MIN_R":inicio,"MAX_R":inicio,"VC":inicio})

        print("PARTE 2 Falta: {}".format(df["NEO F"][i]))


neos_ = [_2019pg1,_2021af8,_2005vc,_2005cz36,_29075]
neos = ['2019pg1','2021af8','2005vc','2005cz36','29075']

for i in range(0,len(neos)):
    inicializacao(neos[i])

#distancia_f(neos_[int(sys.argv[1])],0,360,1)
#grafico_df(neos[int(sys.argv[1])],1)
#for tempo in [1,10,100]:
#    colisao_neo_satelite_parte1(neos_[int(sys.argv[1])],1,60*60*tempo,2.5,2.5,2.5,1e-2,10)
colisao_neo_satelite_parte2(neos_[int(sys.argv[1])],1,60*60*1,2.5,2.5,2.5,1e-2,10)
