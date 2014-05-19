
from kcorrect import ext

print dir(ext)

z = 0.4
omega = 0.3
omegal = 0.7
print ext.ztor(z, omega, omegal)

a = ext.k_read_ascii_table('kcorrect/data/templates/vmatrix.default.dat')
print a
print a.shape

import matplotlib.pyplot as plt
for i in xrange(a.shape[0]):
    plt.plot(a[i] / max(a[i]))
plt.show()
