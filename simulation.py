import time
import matplotlib.pyplot as plt
from matplotlib import axes
from actin import ActinState as AS, Actin as A, ActinPolymer as AP
from dynamics import *

c = 1.38
pol = AP([A(1)] * 100, c)

dt = 0.01
t = 0.0

times = []
length = []
minus = []
plus = []

Running = True
# while Running:
while t < 1000:
    times.append(t)
    length.append(len(pol.monomers))
    minus.append(pol.added['minus'])
    plus.append(pol.added['plus'])

    pol.update(dt)
    t += dt
    # print(pol)
    # time.sleep(dt)

fig, axs = plt.subplots(3, 1)

axs[0].plot(times, length)
axs[1].plot(times, minus)
axs[2].plot(times, plus)

plt.show(block=True)

print(f'Emergent plus reaction speed: {round(plus[-1] / times[-1], 2)}')
print(f'Predicted plus reaction speed: {round(AP.K['plus'][AS.atp]['on'] * pol.C_monomer - AP.K['plus'][AS.adp]['off'], 2)}')