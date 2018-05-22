import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-10,10,0.5)
expx = np.exp(-x**2)
conv = np.convolve(expx,expx,"same")

plt.plot(x,conv,"-o",label="convo of expxx")
plt.plot(x,expx,"-*",label="expxx")
plt.show()
