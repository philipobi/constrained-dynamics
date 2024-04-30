#!.venv/bin/python3
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation

global i
i = 0
positions = np.loadtxt("positions.txt")
n,_ = positions.shape

def main():

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.axes.set_xlim3d(left=-15,right=15)
    ax.axes.set_ylim3d(bottom=-15,top=15)
    ax.axes.set_zlim3d(bottom=-15,top=15)


    data = positions[i].reshape(-1,3)
    graph=ax.scatter(data[:,0], data[:,1], data[:,2])
    dt=.01

    def update_graph(num):
        global i
        i = (i+1) if i < n else 0
        data = positions[i].reshape(-1,3)
        graph._offsets3d = (data[:,0], data[:,1], data[:,2])

    ani = matplotlib.animation.FuncAnimation(fig, update_graph, 2000, 
                                interval=1, blit=False)

    plt.show()


main()
