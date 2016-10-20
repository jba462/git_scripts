import numpy as np
import scipy as sp
import scipy.stats
import matplotlib.pyplot as plt

x_min = -1.5
x_max = 0

x_interval = 1000j

y_min = 0
y_max = 1

y_interval = 1000j

iterations = 10

x,y=np.ogrid[x_min:x_max:x_interval,y_min:y_max:y_interval]

print('')
print('Grid set')
print('')

c=x + 1j*y
z=0

for g in range(iterations):
        print('Iteration number: ',g)
        z=z**2 + c
        z[np.abs(z) > 10] = 3

z[np.abs(z) > 2] = 3

threshold = 2
mask=np.abs(z) < threshold

print('')
print('Plotting using imshow()')
plt.imshow(mask.T,extent=[x_min,x_max,y_min,y_max])

print('')
print('plotting done')
print('')

plt.gray()

print('')
print('Preparing to render')
print('')

plt.show()





