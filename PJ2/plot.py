# dots from result.txt with format y,x

import matplotlib.pyplot as plt

x = []
y = []

with open('answer2.txt', 'r') as f:

    # 散点
    for line in f:
        line = line.strip()
        # remove ()
        line = line.replace('(', '')
        line = line.replace(')', '')
        xbegin = float(line.split(',')[0])
        xend = float(line.split(',')[1])
        ybegin = float(line.split(',')[2])
        yend = float(line.split(',')[3])
        if (xbegin == xend and ybegin == yend):
            x.append(xbegin)
            y.append(ybegin)
        plt.plot([xbegin, xend], [ybegin, yend])

plt.scatter(x, y)
plt.axis('scaled')
plt.show()
