import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


data = pd.read_csv("/home/nadzhou/Desktop/bonds.csv")


interaction = np.array(data['interaction'])
value = np.array(data['value'])


sns.barplot(x=interaction, y=value)
plt.show()
