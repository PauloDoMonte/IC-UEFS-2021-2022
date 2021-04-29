import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import numpy as np

def graficocalor_raiohill(asteroide,complemento):
    for i in [1,2,3,4,5,6,7,8,9,10]:
        arquivo = '{}/raio_hill/{}/2019pg1-{}x.csv'.format(asteroide,complemento,i)
        db = pd.read_csv(arquivo)

        filtro = db['ANO'] < 99999
        db_filtrado = db[filtro]

        df_new = db_filtrado.pivot(index='F NEO',columns='F TERRA',values='ANO')
        heat_map = sb.heatmap(df_new)
        titulo = "Tempo para atingir {}x Raio de Hill da Terra em função das anomalias verdadeira com o maximo de {} anos".format(i,1000)
        plt.ylabel("Anomalia verdadeira do Neo")
        plt.xlabel("Anomalia verdadeira da Terra")
        plt.title(titulo)
        plt.show()

def histograma_raiohill(asteroide, complemento):
    for rh in [1,2,3,4,5,6,7,8,9,10]:
        controle = [0,0,0,0,0]
        total = 0
        arquivo = '{}/raio_hill/{}/2019pg1-{}x.csv'.format(asteroide,complemento,rh)
        db = pd.read_csv(arquivo)

        filtro = db['ANO'] < 99999
        db_filtrado = db[filtro]

        for i in db_filtrado['ANO']:
            if(i > 0 and i <= 200):
                controle[0] += 1
                total += 1
            elif(i > 200 and i <= 400):
                controle[1] += 1
                total += 1
            elif(i > 400 and i <= 600):
                controle[2] += 1
                total += 1
            elif(i > 600 and i <= 800):
                controle[3] += 1
                total += 1
            elif(i > 800 and i <= 1000):
                controle[4] += 1
                total += 1

        print("Raio de hill: {}".format(rh))
        print("""\t0-200: {}
        200-400: {}
        400-600: {}
        600-800: {}
        800-1000: {}
""".format(controle[0],controle[1],controle[2],controle[3],controle[4]))

histograma_raiohill('2019pg1','passo=ano-10')
