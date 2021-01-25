#import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from utils import pipeline
import pdb

def animate(i, xi, yi, G):
    plt.contourf(xi, yi, G[i], 15, vmin=0.0, vmax=1.8, cmap=plt.cm.jet)
    plt.title('Velocity Contour at Time Step: %d' % i)

def main():
    data_dir = '../data/flow_solution'
    airfoil_path = '../data/airfoil'
    z_slice = 0.0125
    scale =  20
    alpha = 8
    Nx = 200
    Ny = 120
    G, xi, yi, x2, y2 = pipeline(alpha, scale, z_slice, Nx, Ny, data_dir, airfoil_path)
    i = 0
    fig = plt.figure(figsize=(12,4),dpi=100)
    plt.contourf(xi, yi, G[i], 15, vmin=0.0, vmax=1.8, cmap=plt.cm.jet)
    plt.fill(x2,y2,fill=True, lw=3)
    plt.colorbar()
    plt.xlim(-10,40)
    plt.ylim(-12,12)

    anim = FuncAnimation(fig, animate, fargs=(xi, yi, G), frames=len(G), interval=1000, repeat=True)
    fig.show()
    anim.save('animation.mp4')

if __name__ == "__main__":
    main()