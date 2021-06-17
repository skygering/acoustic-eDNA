from itertools import repeat
import matplotlib.pyplot as plt
import numpy as np

x = np.array([*range(0,29)])
y=np.array([*range(0,5)] + [*repeat(5,5)] + [*range(6,10)] +  [*repeat(10,5)] + [*range(9,-1, -1)])
print(len(x))
print(len(y))

slopes =[]
for i in range(len(x)-1):
    s = (y[i+1] - y[i])/(x[i+1] - x[i])
    slopes.append(s)
slopes.append(slopes[-1]) # add on value for last element
print(slopes)

p1 = []
p2 = []
for j in range(len(slopes)):
    if slopes[j] ==0 or slopes[j-1] == 0: # likely place for an error due to indexing
        p1.append((x[j], y[j]))


plt.plot(x,y)
plt.plot(*zip(*p1))
plt.show()
plt.close()