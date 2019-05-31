import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


rate = np.loadtxt("rn.txt")
print(rate.shape)

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
ax = Axes3D(fig)
ax.plot_surface(angle, energy, rate)
ax.set_xlabel("angle")
ax.set_ylabel("energy")
ax.set_zlabel("sputtering rate")

plt.show()

pdf_file_name = "figures/sputtering_rate_"
fig.savefig(pdf_file_name, dpi = 300)