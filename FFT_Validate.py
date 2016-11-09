from matplotlib import pyplot as plt
import numpy as np
import math as mt

x=np.zeros(32)
y=np.copy(x)
yt=np.copy(x)
dTheta=mt.pi*2.0/x.size
#set up x on interval 0:2PI
for i in range(x.size):
	x[i]=dTheta*i
#assign y values
for i in range(x.size):
	#y[i]=3.0*mt.sin(2.0*mt.pi*x[i])+5.0*mt.cos(2.0*mt.pi*x[i]/3.0)
	y[i]=mt.cos(x[i])
yt=np.fft.fft(y)
print(yt)
plt.figure(1)
plt.plot(x,y)
plt.show()
