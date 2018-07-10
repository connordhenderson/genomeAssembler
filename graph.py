import matplotlib.pyplot as plt
import sys

class graph:
    def __init__(self, name, data, labels=["x","y"]):
        self.name = name
        self.ydata = data[0]
        self.xdata = None
        if len(data) == 2:
            self.xdata = data[1]
        self.xlabel = labels[0]
        self.ylabel = None
        if len(labels) == 2:
            self.ylabel = labels[1]

    def show(self):
        plt.plot(self.ydata)
        plt.xlabel = self.xlabel
        plt.ylabel = self.ylabel
        sys.stdout.flush()
        plt.show()
