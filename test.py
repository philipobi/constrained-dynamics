import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from utils import Simulation

#a = np.random.rand(2000, 3)*10
#t = np.array([np.ones(100)*i for i in range(20)]).flatten()
#df = pd.DataFrame({"time": t ,"x" : a[:,0], "y" : a[:,1], "z" : a[:,2]})




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('3D Test')
ax.axes.set_xlim3d(left=-30,right=30)
ax.axes.set_ylim3d(bottom=-30,top=30)
ax.axes.set_zlim3d(bottom=-30,top=30)


sim = Simulation(10,10,3)

data=sim.q.reshape(-1,3)
graph=ax.scatter(data[:,0], data[:,1], data[:,2])
dt=.01

def update_graph(num):
    #input()
    sim.run(dt)
    data=sim.q.reshape(-1,3)
    graph._offsets3d = (data[:,0], data[:,1], data[:,2])

ani = matplotlib.animation.FuncAnimation(fig, update_graph, 2000, 
                               interval=500, blit=False)

plt.show()