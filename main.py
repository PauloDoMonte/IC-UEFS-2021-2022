import math, rebound, os, sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

"""
Programa feito com os dados do dia 09 de Março de 2022

Só precisa acessar o site https://ssd.jpl.nasa.gov/horizons/app.html#/ Para atulaizar os dados

Ordem das particulas
    0 -> sol
    1 -> terra
    2 -> lua
    3 -> asteroide de estudo
    4 -> veiculo (se existir)
    5 -> outros corpos para deixar mais real a simulacao
"""

e = math.e

class particula:
    def __init__(self,nome,massa,raio,x,y,z,vx,vy,vz):
        self.nome = nome
        self.m = massa
        self.r = raio
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

#Vetores de posição são [x,y,z,vx,vy,vz] [Km, Km/s]

# Corpos maiores
sol         = particula('sol',1.989e+30,696340,0,0,0,0,0,0)
terra       = particula('terra',5.973332e+24,6378.1366,-1.453091711410064E+08,3.065055370063527E+07,-7.925042448136955E+02,-6.623093040898302E+00,-2.926662424847545E+01,4.500915920360171E-04)
lua         = particula("lua",7.349e+22,1737.4,-1.451227069269194E+08,3.100640892359096E+07,4.048109483337030E+03,-7.468662400475338E+00,-2.879215349084939E+01,8.620510923978841E-02)
marte       = particula("marte",6.39e23,3389.5,7.847166899591953E+06,-2.169039184369964E+08,-4.738304891883358E+06,2.512537869754003E+01,2.960113422561584E+00,-5.542873414249656E-01)
venus       = particula("venus",4.867e24,6051.8,-9.851920681168786E+07,-4.392016022776731E+07,5.081886511314837E+06,1.401894735725095E+01,-3.214702470522624E+01,-1.250217084123532E+00)
mercurio    = particula("mercurio",3.285e23,2439.7,1.153980244513476E+07,-6.717383885296583E+07,-6.547869499216489E+06,3.824054643954008E+01,1.074356995196843E+01,-2.629738621953904E+00)
jupiter     = particula("jupiter",1.898e27,69911,7.193793808177494E+08,-1.929473860559401E+08,-1.529344941711865E+07,3.233183651658163E+00,1.324605788582764E+01,-1.274033663158294E-01)
saturno     = particula("saturno",5.683e26,58232,1.077248375084723E+09,-1.017244140420120E+09,-2.518456225300723E+07,6.098570465953345E+00,7.016951543580232E+00,-3.651793127808300E-01)
urano       = particula("urano",8.681e25,25362,2.126454968859567E+09,2.043228090871781E+09,-1.996870870835972E+07,-4.765608297192462E+00,4.606318969096011E+00,7.873895603154546E-02)
netuno      = particula("netuno",1.024e26,24622,4.437085751516443E+09,-5.804937498670350E+08,-9.030853375854236E+07,6.740998684630994E-01,5.436435612720855E+00,-1.279439628294634E-01)
plutao      = particula("plutao",1.307e22,1188.3,2.291946160811365E+09,-4.617640558055246E+09,-1.685448919017437E+08,4.993574173877477E+00,1.253568218382351E+00,-1.558734198633245E+00)

# Luas
ganymede    = particula("ganymede",1.4819e23,2634.1,7.202158835075567E+08,-1.922810341626214E+08,-1.525651807149172E+07,-3.543088673762314E+00,2.176103907656820E+01,1.036583241782791E-01)
titan       = particula("titan",1.3452e23,2574.7,1.076008248198582E+09,-1.017330937315913E+09,-2.501623003992635E+07,6.819835097862795E+00,2.208530978820356E+00,2.042213961357846E+00)

# Asteroides
_53319      = particula("53319",80000,7/2,2.003279550628811E+08,9.215005923352031E+07,-5.139847022506043E+07,2.527386161372272E+00,2.849901843765103E+01,-5.295103263965284E+00)
_2019pg1    = particula("2019pg1",80000,0.270/2,-3.699502640466033E+08,3.633794785576811E+08,-1.071690952445368E+08,-3.017912783426565E+00,-1.128920619033803E+01,-1.939558069024700E+00)
_2021af8    = particula("2021af8",80000,0.310/2,4.022417713216301E+08,2.319215361430710E+07,-4.349310479420440E+07,5.397642054765392E+00,1.360278945736246E+01,1.090835142153088E+00)
#_2005vc     = particula("2005vc")
#_2005cz36   = particula("2005cz36")
#_29075      = particula("29075")

def inicializacao(asteroide):
    if(os.path.exists('{}/'.format(asteroide)) == True):
        if(os.path.exists('{}/grafico'.format(asteroide))== True):
            pass
        else:
            os.mkdir('{}/grafico/'.format(asteroide))
    else:
        os.mkdir('{}/'.format(asteroide))
        os.mkdir('{}/grafico/'.format(asteroide))

def pre_colisao(neo):

    inicio = []
    db = pd.DataFrame({"d_neo_terra":inicio,"d_terra_lua":inicio,"ano":inicio})
    caminho = "{}/pre_colisao.csv".format(neo.nome)
    db.to_csv(caminho,sep=",")

    anos = 10000
    passo = 365
    tempo_simulacao = 60*60*24*365*anos
    times = np.linspace(0, tempo_simulacao, anos*passo)

    sim = rebound.Simulation()
    sim.units = ("kg","km","s")
    sim.integrator = "IAS15"

    sim.add(m=sol.m,r=sol.r,x=sol.x,y=sol.y,z=sol.z,vx=sol.vx,vy=sol.vy,vz=sol.vz)
    sim.add(m=terra.m,r=terra.r,x=terra.x,y=terra.y,z=terra.z,vx=terra.vx,vy=terra.vy,vz=terra.vz)
    sim.add(m=lua.m,r=lua.r,x=lua.x,y=lua.y,z=lua.z,vx=lua.vx,vy=lua.vy,vz=lua.vz)
    sim.add(m=neo.m,r=neo.r,x=neo.x,y=neo.y,z=neo.z,vx=neo.vx,vy=neo.vy,vz=neo.vz)

    # Planetas
    sim.add(m=venus.m,r=venus.r,x=venus.x,y=venus.y,z=venus.z,vx=venus.vx,vy=venus.vy,vz=venus.vz)
    sim.add(m=marte.m,r=marte.r,x=marte.x,y=marte.y,z=marte.z,vx=marte.vx,vy=marte.vy,vz=marte.vz)
    sim.add(m=mercurio.m,r=mercurio.r,x=mercurio.x,y=mercurio.y,z=mercurio.z,vx=mercurio.vx,vy=mercurio.vy,vz=mercurio.vz)
    sim.add(m=jupiter.m,r=jupiter.r,x=jupiter.x,y=jupiter.y,z=jupiter.z,vx=jupiter.vx,vy=jupiter.vy,vz=jupiter.vz)
    sim.add(m=saturno.m,r=saturno.r,x=saturno.x,y=saturno.y,z=saturno.z,vx=saturno.vx,vy=saturno.vy,vz=saturno.vz)
    sim.add(m=urano.m,r=urano.r,x=urano.x,y=urano.y,z=urano.z,vx=urano.vx,vy=urano.vy,vz=urano.vz)
    sim.add(m=netuno.m,r=netuno.r,x=netuno.x,y=netuno.y,z=netuno.z,vx=netuno.vx,vy=netuno.vy,vz=netuno.vz)

    # Luas
    sim.add(m=ganymede.m,r=ganymede.r,x=ganymede.x,y=ganymede.y,z=ganymede.z,vx=ganymede.vx,vy=ganymede.vy,vz=ganymede.vz)


    for i, time in enumerate(times):
        sim.integrate(time)

        d_neo_terra = math.sqrt(((sim.particles[3].x-sim.particles[1].x)**2 + (sim.particles[3].y-sim.particles[1].y)**2 + (sim.particles[3].z-sim.particles[1].z)**2))
        d_terra_lua = math.sqrt(((sim.particles[2].x-sim.particles[1].x)**2 + (sim.particles[2].y-sim.particles[1].y)**2 + (sim.particles[2].z-sim.particles[1].z)**2))

        db = db.append({"d_neo_terra":d_neo_terra,"d_terra_lua":d_terra_lua,"ano":time/(60*60*24*365)},ignore_index=True)
        print(len(times)-i)

        if(len(db)>200):
            db.to_csv(caminho,sep=",", mode='a',header=False)
            db = pd.DataFrame({"d_neo_terra":inicio,"d_terra_lua":d_terra_lua,"ano":inicio})

    db.to_csv(caminho,sep=",", mode='a',header=False)

def colisao(neo):
    inicio = []
    db = pd.DataFrame({"d_sat_neo_0":inicio,"d_sat_terra_0":inicio,"d_neo_terra_0":inicio,
                       "vex":inicio, "vey":inicio, "vez":inicio, "gamma":inicio, "chi":inicio,"t_queima":inicio,
                       "v_colisao":inicio})

    caminho = "{}/colisao.csv".format(neo.nome)
    db.to_csv(caminho,sep=",")

    n_neo = 3

    for ve in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
        vex,vey,vez = ve/math.sqrt(3),ve/math.sqrt(3),ve/math.sqrt(3)

        for gamma in [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10]:
            for chi in [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]:
                for vc in range(1,51):
                    vc = vc/math.sqrt(3) #Colisão no eixo xyz
                    for tempo_simulacao in [24,24*2,24*3,24*4,24*5,24*6,24*7,24*8,24*9,24*10,24*11,24*12,24*13,24*14,24*15,24*16,24*17,24*18,24*19,24*20]:
                        tempo_simulacao = tempo_simulacao*3600

                        def propulsao(simp):
                            sim = simp.contents
                            particles = sim.particles

                            particles[n_neo+1].ax += ((-vex*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
                            particles[n_neo+1].ay += ((-vey*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
                            particles[n_neo+1].az += ((-vez*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))

                        sim = rebound.Simulation()
                        sim.units = ("kg","km","s")
                        sim.integrator = "IAS15"

                        sim.additional_forces = propulsao
                        sim.force_is_velocity_dependent = False

                        sim.add(m=sol.m,r=sol.r,x=sol.x,y=sol.y,z=sol.z,vx=sol.vx,vy=sol.vy,vz=sol.vz)
                        sim.add(m=terra.m,r=terra.r,x=terra.x,y=terra.y,z=terra.z,vx=terra.vx,vy=terra.vy,vz=terra.vz)
                        sim.add(m=lua.m,r=lua.r,x=lua.x,y=lua.y,z=lua.z,vx=lua.vx,vy=lua.vy,vz=lua.vz)
                        sim.add(m=neo.m,r=neo.r,x=neo.x,y=neo.y,z=neo.z,vx=neo.vx,vy=neo.vy,vz=neo.vz)
                        sim.add(m=1000,r=0.01,x=sim.particles[n_neo].x+1e-3,y=sim.particles[n_neo].y+1e-3,z=sim.particles[n_neo].z+1e-3,
                            vx=sim.particles[n_neo].vx+vc,vy=sim.particles[n_neo].vy+vc,vz=sim.particles[n_neo].vz+vc)

                        sim.add(m=venus.m,r=venus.r,x=venus.x,y=venus.y,z=venus.z,vx=venus.vx,vy=venus.vy,vz=venus.vz)
                        sim.add(m=marte.m,r=marte.r,x=marte.x,y=marte.y,z=marte.z,vx=marte.vx,vy=marte.vy,vz=marte.vz)
                        sim.add(m=mercurio.m,r=mercurio.r,x=mercurio.x,y=mercurio.y,z=mercurio.z,vx=mercurio.vx,vy=mercurio.vy,vz=mercurio.vz)
                        sim.add(m=jupiter.m,r=jupiter.r,x=jupiter.x,y=jupiter.y,z=jupiter.z,vx=jupiter.vx,vy=jupiter.vy,vz=jupiter.vz)
                        sim.add(m=saturno.m,r=saturno.r,x=saturno.x,y=saturno.y,z=saturno.z,vx=saturno.vx,vy=saturno.vy,vz=saturno.vz)
                        sim.add(m=urano.m,r=urano.r,x=urano.x,y=urano.y,z=urano.z,vx=urano.vx,vy=urano.vy,vz=urano.vz)
                        sim.add(m=netuno.m,r=netuno.r,x=netuno.x,y=netuno.y,z=netuno.z,vx=netuno.vx,vy=netuno.vy,vz=netuno.vz)

                        sim.integrate(-tempo_simulacao)

                        #print(sim.particles[0])
                        #print(sim.particles[1])
                        #print(sim.particles[2])
                        #print(sim.particles[3])
                        #print(sim.particles[4])
                        #fig = rebound.OrbitPlot(sim)
                        #plt.show()
                        #sys.exit()

                        d_sat_neo_0     = math.sqrt(((sim.particles[4].x-sim.particles[3].x)**2 + (sim.particles[4].y-sim.particles[3].y)**2 + (sim.particles[4].z-sim.particles[3].z)**2))
                        d_sat_terra_0   = math.sqrt(((sim.particles[4].x-sim.particles[1].x)**2 + (sim.particles[4].y-sim.particles[1].y)**2 + (sim.particles[4].z-sim.particles[1].z)**2))
                        d_neo_terra_0   = math.sqrt(((sim.particles[3].x-sim.particles[1].x)**2 + (sim.particles[3].y-sim.particles[1].y)**2 + (sim.particles[3].z-sim.particles[1].z)**2))

                        print("VE:{}\tGAMMA:{}\tCHI:{}\tVC:{}\tT:{}\tNeo:{}".format(ve,gamma,chi,vc,tempo_simulacao,neo.nome))

                        db = db.append({"d_sat_neo_0":d_sat_neo_0 ,"d_sat_terra_0":d_sat_terra_0,"d_neo_terra_0":d_neo_terra_0,
                                           "vex":vex, "vey":vey, "vez":vez, "gamma":gamma, "chi":chi,"t_queima":tempo_simulacao,
                                           "v_colisao":vc*math.sqrt(3)},ignore_index=True)

                if(len(db)> 100):
                    db.to_csv(caminho, sep=",", header=False, mode='a')
                    db = pd.DataFrame({"d_sat_neo_0":inicio,"d_sat_terra_0":inicio,"d_neo_terra_0":inicio,
                                       "vex":inicio, "vey":inicio, "vez":inicio, "gamma":inicio, "chi":inicio,"t_queima":inicio,
                                       "v_colisao":inicio})

    db.to_csv(caminho,sep=",",header=False,mode='a')

def pos_colisao(neo):

    anos = 100
    passo = 365
    tempo_simulacao = 60*60*24*365*anos
    times = np.linspace(0, tempo_simulacao, anos*passo)

    for vc in range(1,51):

        inicio = []
        db = pd.DataFrame({"d_neo_terra":inicio,"d_terra_lua":inicio,"ano":inicio})
        caminho = "{}/pos_colisao_{}km_s.csv".format(neo.nome,vc)
        db.to_csv(caminho,sep=",")

        sim = rebound.Simulation()
        sim.units = ("kg","km","s")
        sim.integrator = "IAS15"

        sim.add(m=sol.m,r=sol.r,x=sol.x,y=sol.y,z=sol.z,vx=sol.vx,vy=sol.vy,vz=sol.vz)
        sim.add(m=terra.m,r=terra.r,x=terra.x,y=terra.y,z=terra.z,vx=terra.vx,vy=terra.vy,vz=terra.vz)
        sim.add(m=lua.m,r=lua.r,x=lua.x,y=lua.y,z=lua.z,vx=lua.vx,vy=lua.vy,vz=lua.vz)
        sim.add(m=neo.m,r=neo.r,x=neo.x,y=neo.y,z=neo.z,vx=neo.vx+vc,vy=neo.vy+vc,vz=neo.vz+vc)

        sim.add(m=venus.m,r=venus.r,x=venus.x,y=venus.y,z=venus.z,vx=venus.vx,vy=venus.vy,vz=venus.vz)
        sim.add(m=marte.m,r=marte.r,x=marte.x,y=marte.y,z=marte.z,vx=marte.vx,vy=marte.vy,vz=marte.vz)
        sim.add(m=mercurio.m,r=mercurio.r,x=mercurio.x,y=mercurio.y,z=mercurio.z,vx=mercurio.vx,vy=mercurio.vy,vz=mercurio.vz)
        sim.add(m=jupiter.m,r=jupiter.r,x=jupiter.x,y=jupiter.y,z=jupiter.z,vx=jupiter.vx,vy=jupiter.vy,vz=jupiter.vz)
        sim.add(m=saturno.m,r=saturno.r,x=saturno.x,y=saturno.y,z=saturno.z,vx=saturno.vx,vy=saturno.vy,vz=saturno.vz)
        sim.add(m=urano.m,r=urano.r,x=urano.x,y=urano.y,z=urano.z,vx=urano.vx,vy=urano.vy,vz=urano.vz)
        sim.add(m=netuno.m,r=netuno.r,x=netuno.x,y=netuno.y,z=netuno.z,vx=netuno.vx,vy=netuno.vy,vz=netuno.vz)

        for i, time in enumerate(times):
            sim.integrate(time)

            d_neo_terra = math.sqrt(((sim.particles[3].x-sim.particles[1].x)**2 + (sim.particles[3].y-sim.particles[1].y)**2 + (sim.particles[3].z-sim.particles[1].z)**2))
            d_terra_lua = math.sqrt(((sim.particles[2].x-sim.particles[1].x)**2 + (sim.particles[2].y-sim.particles[1].y)**2 + (sim.particles[2].z-sim.particles[1].z)**2))

            db = db.append({"d_neo_terra":d_neo_terra,"d_terra_lua":d_terra_lua,"ano":time/(60*60*24*365)},ignore_index=True)
            print("{}\t{}".format(vc,len(times)-i))

            if(len(db)>200):
                db.to_csv(caminho,sep=",", mode='a',header=False)
                db = pd.DataFrame({"d_neo_terra":inicio,"d_terra_lua":d_terra_lua,"ano":inicio})

        db.to_csv(caminho,sep=",", mode='a',header=False)

# Graficos pré colisao
def grafico_pre_colisao(neo):
    caminho = "{}/pre_colisao.csv".format(neo.nome)
    db = pd.read_csv(caminho)

    print(min(db['d_neo_terra']))
    print(max(db['d_terra_lua']))

    fig = plt.figure()

    plt.plot(db['ano'],db['d_neo_terra'], marker="o", markersize=1)
    plt.plot(db['ano'],db['d_terra_lua'], marker="o", markersize=1)
    plt.yscale('log')
    plt.legend(["Distancia Terra-Neo","Distancia Terra-Lua"], loc="best")

    plt.grid()
    plt.show()

def grafico_pos_colisao(neo):

    for vc in [2,4,6,8,10,12,14,16,18,20]:
        caminho = "{}/pos_colisao_{}km_s.csv".format(neo.nome,vc)
        db = pd.read_csv(caminho)

        print(min(db['d_neo_terra']))

        legenda = "VC {}km/s".format(vc)
        plt.plot(db['ano'],db['d_neo_terra'], marker="o", markersize=1,label=legenda)


    plt.grid()
    plt.title("Distância Terra-Neo pós colisão do satélite com Neo")
    plt.ylabel("Distância [KM]")
    plt.xlabel("Anos [ANO]")
    plt.legend(loc="best")
    plt.show()



neo = _53319
inicializacao(neo.nome)

while(True):
    terminal = input("DoMonte >>> ")
    if(terminal == "pre colisao"):
        pre_colisao(neo)
    elif(terminal == "colisao"):
        colisao(neo)
    elif(terminal == "pos colisao"):
        pos_colisao(neo)
    elif(terminal == "grafico pre colisao"):
        grafico_pre_colisao(neo)
    elif(terminal == "grafico colisao"):
        grafico_colisao(neo)
    elif(terminal == "grafico pos colisao"):
        grafico_pos_colisao(neo)
