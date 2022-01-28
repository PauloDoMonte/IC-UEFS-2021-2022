import rebound,math,sqlite3
import numpy as np
import pandas as pd
import os, sys
import lib.inicializacao as ini
import lib.matematica as mate

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

_2019pg1    = particula('2019pg1',80000, 400/1000,1.037103597335171*mate.AU,0.03388638947447908,0.140625575916969,281.6394465120655,13.55855529785005,0)
_2021af8    = particula('2021af8',80000, 400/1000,2.018135760648355*mate.AU,0.5151605500121827,9.697802979250167,168.9511312385866,42.56717817033352,0)
_2005vc     = particula('2005vc',80000, 1000/1000, 2.082761852899944*mate.AU,0.5947841108529824,4.483281580745502,290.2042173581335,227.9717174732001,0)
_2005cz36   = particula('2005cz36',80000, 1000/1000,2.238066485936266*mate.AU,0.5747761271925266,16.14636501981154,139.3809136995648,116.8081182427147,0)
_29075      = particula('29075',80000,1.3/2,1.698665726802664*mate.AU,0.5077360606964316,12.16734013772536,224.6831644667536,356.6558403536075,0)
sol         = particula('sol',1.989e+30,696340,0,0,0,0,0,0)
terra       = particula('terra',5.973332e+24,6378.1366,149.60e6,0.01671022,0.00005,-11.26064,102.94719,0)

neos_disponiveis = [_2019pg1,_2021af8,_2005vc,_2005cz36,_29075]
neo = 0

decisao = int(input("1 para rodar todo\t 2 so grafico: "))
if(decisao == 1):
    if(len(sys.argv)>1):
        if(sys.argv[1] == '-h'):
            print("python3 simulacao.py [ASTEROIDE] [TEMPO MAXIMO] [INCREMENTO TEMPO (ANO/n)] [MUTIPLICADOR DO RAIO DE HILL]")
            sys.exit()

    for i in neos_disponiveis:
        ini.inicializacao(i.nome)
        mate.colisao(sol,terra,i)

    mate.grafico(_2019pg1)

if(decisao == 2):
    mate.grafico(_2019pg1)
