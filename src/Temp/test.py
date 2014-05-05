import numpy as np
import matplotlib.pyplot as plt
a=np.random.rand(3,3)
plt.pcolor(a)
plt.plot([1,2,3,4])

def plots():
    plt.plot([2 ,5, 3])

plots()
plt.get_backend()
plt.show()

