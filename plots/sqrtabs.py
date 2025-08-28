#plots sqrt(abs(x)) function

import numpy as np
import matplotlib.pyplot as plt
import math

x = np.linspace(-20,20,800)
y = [math.sqrt(abs(xi)) for xi in x]

plt.plot(x,y)
plt.xlabel("x")
plt.ylabel("sqrt(abs(x))")
plt.savefig("pngs/sqrtabs.png")
