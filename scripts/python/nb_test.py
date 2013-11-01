# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import matplotlib.pyplot as plt

x = range(100)
y = range(100,200)
fig = plt.figure()

ax1 = fig.add_subplot(221,title="ax1")
opts = 'c=\'b\', marker=\"s\"'
print opts
#ax1.scatter(x[:4], y[:4], s=10, opts)
ax1.scatter(y[4:], x[4:], s=10, c='r', marker="o")


ax2 = fig.add_subplot(222, title="ax2")
ax2.scatter(x[:4], y[:4], s=1, c='b', marker="s")
ax2.scatter(y[4:], x[4:], s=1, c='r', marker="o")

#plt.show()
plt.savefig("test2.jpg")

# <codecell>


