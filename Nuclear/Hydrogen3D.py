# 3D visulization of the electron density in hydrogen
import numpy as np
import math
from mpl_toolkits import mplot3d

a0 = 53    #Bohr Radius in pm
Z = 1    #Amount of Protons in the atom

#Wave Functions for orbitals up to n = 4
_3d = (1/(81*math.sqrt(6*np.pi)))*((Z/a0)**(3/2))*(((p**2)*e)**(-p/3))*((3*(math.cos(theta)**2))-1)
_3p = (1/81)*(2/math.pi)*((6*r)-(p**2))*((Z/a0)**(3/2))*((e)**(-p/3))*math.cos(45)
_3s = (1/(81*math.sqrt(3*math.pi)))*(27-(18*p)+(2*(p**2)))*((Z/a0)**(3/2))*(e**(-p/3))
_2p = (1/(4*math.sqrt(2*math.pi)))*((Z/a0)**(3/2))*((p*e)**(-p/(2)))*math.cos(theta)
_2s = (1/(4*math.sqrt(2*math.pi)))*((Z/a0)**(3/2))*(2-((Z*r)/a0))*(e**(-p/(2)))
_1s = (1/math.sqrt(math.pi))*((Z/a0)**(3/2))*(e**((-Z*r)/a0))

#function for the probability
def prob(x,y,z):
    r=np.sqrt(np.square(x)+np.square(y)+np.square(z))
    p = (Z*r)/a0
    psi = _2s    #I chose the 2s orbital
    return np.square(psi)
    
#Random coordinates
x=np.linspace(0,1,30)
y=np.linspace(0,1,30)
z=np.linspace(0,1,30)

#finding the probability for each random coordinate
elements = []
probability = []
for ix in x:
    for iy in y:
        for iz in z:
            #Serialize into 1D object
            elements.append(str((ix,iy,iz)))
            probability.append(prob(ix,iy,iz))
            
#Make the sum of the probability 1
probability = probability/sum(probability)

#Getting electron coordinates based on probabiliy
coord = np.random.choice(elements, size=5000, replace=True, p=probability)
elem_mat = [i.split(',') for i in coord]
elem_mat = np.matrix(elem_mat)
x_coords = [float(i.item()[1:]) for i in elem_mat[:,0]] 
y_coords = [float(i.item()) for i in elem_mat[:,1]] 
z_coords = [float(i.item()[0:-1]) for i in elem_mat[:,2]]

#Finding the density for the graph
from scipy.stats import gaussian_kde
from scipy import stats
from mayavi import mlab
xyz = np.vstack([x_coords,y_coords,z_coords])
kde = stats.gaussian_kde(xyz)
density = kde(xyz)

#Plotting the graph
figure = mlab.figure('DensityPlot')
pts = mlab.points3d(x_coords, y_coords, z_coords, density, scale_mode='none', scale_factor=0.03)
mlab.axes()
mlab.show()
