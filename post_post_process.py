import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('settlement.txt').T[0]
plt.hist(data, bins=10)
plt.show()
