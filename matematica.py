import math, rebound
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

AU = 149597871
e = math.e

def raio_hills(m1,m2,a,e):
    return(a*(1-e)*(pow(m2/(3*m1),1/3)))

def colisao(sol,terra,neo,nome_simulacao):

    def propulsao(simp):
        sim = simp.contents
        particles = sim.particles

        particles[3].ax += ((-vex*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
        particles[3].ay += ((-vey*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
        particles[3].az += ((-vez*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))

    inicio = []
    db = pd.DataFrame({"neo_f0":inicio,"terra_f0": inicio, "neo_ff": inicio, "terra_ff":inicio,
                       "d_sat_neo_0":inicio,"d_sat_terra_0":inicio,"d_neo_terra_0":inicio,
                       "xis":inicio,"yis":inicio,"zis":inicio,"vxis":inicio,"vyis":inicio,"vzis":inicio,
                       "xin":inicio,"yin":inicio,"zin":inicio,"vxin":inicio,"vyin":inicio,"vzin":inicio,
                       "vx":inicio,"vy":inicio,"vz":inicio,"gamma":inicio,"chi":inicio,"t_queima":inicio,
                       "v_colisao":inicio})

    caminho = "{}/{}_colisao.csv".format(neo.nome,nome_simulacao)
    db.to_csv(caminho, sep=",")


    for terra_f in range(0,360,1):

        for neo_f in range(0,360,1):

            for ve in [5]:
                vex,vey,vez = ve/math.sqrt(3),ve/math.sqrt(3),ve/math.sqrt(3)

                for gamma in [1e-1]:

                    for chi in [10]:

                        for vc in [20]:
                            vc = vc/math.sqrt(2)

                            for tempo_simulacao in [1,10,100]:
                                tempo_simulacao = tempo_simulacao*3600

                                def propulsao(simp):
                                    sim = simp.contents
                                    particles = sim.particles

                                    particles[3].ax += ((-vex*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
                                    particles[3].ay += ((-vey*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))
                                    particles[3].az += ((-vez*gamma*pow(e,(-1*gamma)*(-sim.t)))/(chi+pow(e,(-1*gamma)*(-sim.t))))

                                sim = rebound.Simulation()
                                sim.units = ("kg","km","s")
                                sim.integrator = "IAS15"

                                sim.additional_forces = propulsao
                                sim.force_is_velocity_dependent = False

                                sim.add(m=sol.m,r=sol.r)
                                sim.add(m=terra.m,r=terra.r,a=terra.a,e=terra.e,inc=terra.inc,omega=terra.omega,Omega=terra.Omega,f=terra_f*0.0174533)
                                sim.add(m=neo.m,r=neo.r,a=neo.a,e=neo.e,inc=neo.inc,omega=neo.omega,Omega=neo.Omega,f=neo_f*0.0174533)
                                sim.add(m=1000,r=0.01,x=sim.particles[2].x+1e-3,y=sim.particles[2].y,z=sim.particles[2].z,
                                    vx=sim.particles[2].vx+vc,vy=sim.particles[2].vy+vc,vz=sim.particles[2].vz)

                                rfs = [sim.particles[3].x,sim.particles[3].y,sim.particles[3].z]
                                vfs = [sim.particles[3].x,sim.particles[3].y,sim.particles[3].z]

                                rfn = [sim.particles[2].x,sim.particles[2].y,sim.particles[2].z]
                                vfn = [sim.particles[2].vx,sim.particles[2].vy,sim.particles[2].vz]

                                neo_ff      = neo_f
                                terra_ff    = terra_f

                                sim.move_to_com()
                                sim.integrate(-tempo_simulacao)

                                ris = [sim.particles[3].x,sim.particles[3].y,sim.particles[3].z]
                                vis = [sim.particles[3].x,sim.particles[3].y,sim.particles[3].z]

                                rin = [sim.particles[2].x,sim.particles[2].y,sim.particles[2].z]
                                vin = [sim.particles[2].vx,sim.particles[2].vy,sim.particles[2].vz]

                                neo_f0      = sim.particles[2].f*57.2958
                                terra_f0    = sim.particles[1].f*57.2958

                                d_sat_neo_0     = math.sqrt(((sim.particles[3].x-sim.particles[2].x)**2 + (sim.particles[3].y-sim.particles[2].y)**2 + (sim.particles[3].z-sim.particles[2].z)**2))
                                d_sat_terra_0   = math.sqrt(((sim.particles[3].x-sim.particles[1].x)**2 + (sim.particles[3].y-sim.particles[1].y)**2 + (sim.particles[3].z-sim.particles[1].z)**2))
                                d_neo_terra_0   = math.sqrt(((sim.particles[2].x-sim.particles[1].x)**2 + (sim.particles[2].y-sim.particles[1].y)**2 + (sim.particles[2].z-sim.particles[1].z)**2))


                                print("T_F:{}\tN_F:{}\tVE:{}\tGAMMA:{}\tCHI:{}\tVC:{}\tT:{}\tNeo:{}".format(terra_f,neo_f,ve,gamma,chi,vc,tempo_simulacao,neo.nome))

                                db = db.append({"neo_f0":neo_f0,"terra_f0": terra_f0, "neo_ff":neo_ff, "terra_ff":terra_ff,
                                                "d_sat_neo_0":d_sat_neo_0,"d_sat_terra_0":d_sat_terra_0,"d_neo_terra_0":d_neo_terra_0,
                                                "xis":ris[0],"yis":ris[1],"zis":ris[2],"vxis":vis[0],"vyis":vis[1],"vzis":vis[2],
                                                "xin":rin[0],"yin":rin[1],"zin":rin[2],"vxin":vin[0],"vyin":vin[1],"vzin":vin[2],
                                                "vx":vex,"vy":vey,"vz":vez,"gamma":gamma,"chi":chi,"t_queima":tempo_simulacao,
                                                "v_colisao":vc*math.sqrt(2)},ignore_index=True)

                        if(len(db)> 100):
                            db.to_csv(caminho, sep=",", header=False, mode='a')
                            db = pd.DataFrame({"neo_f0":inicio,"terra_f0": inicio, "neo_ff":inicio, "terra_ff":inicio,
                                               "d_sat_neo_0":inicio,"d_sat_terra_0":inicio,"d_neo_terra_0":inicio,
                                               "xis":inicio,"yis":inicio,"zis":inicio,"vxis":inicio,"vyis":inicio,"vzis":inicio,
                                               "xin":inicio,"yin":inicio,"zin":inicio,"vxin":inicio,"vyin":inicio,"vzin":inicio,
                                               "vx":inicio,"vy":inicio,"vz":inicio,"gamma":inicio,"chi":inicio,"t_queima":inicio,
                                               "v_colisao":inicio})

    db.to_csv(caminho,sep=",",header=False,mode='a')

def grafico(neo,nome_grafico):

    #nome_grafico = ['maiores_valores','menores_valores','variando_chi','variando_gamma','variando_ts',
    #'variando_vc','variando_ve']

    caminho = "{}/{}_colisao.csv".format(neo.nome,nome_grafico)
    db = pd.read_csv(caminho)

    if(nome_grafico == 'variando_ts'):

        plt.plot(db['t_queima'],['d_sat_neo_0'])
        plt.grid()
        plt.show()
