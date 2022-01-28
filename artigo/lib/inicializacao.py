import os, sys

def inicializacao(asteroide):
    if(os.path.exists('{}/'.format(asteroide)) == True):
        if(os.path.exists('{}/grafico'.format(asteroide))== True):
            pass
        else:
            os.mkdir('{}/grafico/'.format(asteroide))
    else:
        os.mkdir('{}/'.format(asteroide))
        os.mkdir('{}/grafico/'.format(asteroide))
