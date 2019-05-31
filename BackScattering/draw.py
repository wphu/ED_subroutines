import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


rate_n = np.loadtxt("rn.txt")
rate_e = np.loadtxt("re.txt")
print(rate_n.shape)
print(rate_e.shape)

angle, energy = np.mgrid[slice(1.0, 91.0, 1.0), slice(1.0, 501.0, 1.0)]
print(angle.shape, energy.shape)

fig=plt.figure(figsize=(10,8))

'''
# 2d figure
ax = fig.add_subplot(1,1,1)
ax.plot(angle[:, 40], rate[:, 40])
ax.plot(energy[40, :], rate[40, :])
'''

# 3d figure
ax0 = fig.add_subplot(1, 2, 1, projection = '3d')
ax0.plot_surface(angle, energy, rate_n)
ax0.set_xlabel("angle")
ax0.set_ylabel("energy")
ax0.set_zlabel("rn")

ax1 = fig.add_subplot(1, 2, 2, projection = '3d')
ax1.plot_surface(angle, energy, rate_e)
ax1.set_xlabel("angle")
ax1.set_ylabel("energy")
ax1.set_zlabel("re")


plt.show()

file_name = "figures/rn_re_"
fig.savefig(file_name, dpi = 300)