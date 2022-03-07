import pandas as pd
import matplotlib.pyplot as plt

db = pd.read_csv("sem_efeito.csv")
plt.plot(db['ano'],db['semi_eixo sem efeito'])
plt.show()
