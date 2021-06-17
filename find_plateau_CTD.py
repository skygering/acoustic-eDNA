from itertools import repeat
import matplotlib.pyplot as plt
import numpy as np

x = np.array([*range(0,29)])
y=np.array([*range(0,5)] + [*repeat(5,5)] + [*range(6,10)] +  [*repeat(10,5)] + [*range(9,-1, -1)])
#print(len(x))
#print(len(y))

def calc_slope(x1, y1, x2, y2):
    return (y2 - y1)/(x2 - x1)

seg_num = 0
segments = [[]]
current_slope = calc_slope(x[0], y[0], x[1], y[1])
for j in range(len(x)-1):
    new_slope = calc_slope(x[j], y[j], x[j+1], y[j+1])
    if new_slope != current_slope:
        if current_slope == 0:
            segments[seg_num].append((x[j], y[j]))
            segments.append([])
        else:
            segments.append([(x[j], y[j])])
        current_slope = new_slope
        seg_num += 1
    else:
        segments[seg_num].append((x[j], y[j]))
#print(segments)

for segs in segments:
    plt.plot(*zip(*segs))
#plt.plot(x,y)
#plt.plot(*zip(*p1))
plt.show()
plt.close()