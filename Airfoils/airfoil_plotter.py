# NACA Airfoil Plotter
# Kabir Khwaja, 6/10/24

# INPUTS: NACA Airfoil number, number of points desired
# OUTPUTS: plot, point file
  # point file format: "NACAXXXX_YYY.txt"
  # XXXX = airfoil number, YYY = number of points

# airfoil format: NACA ABCD
  # A = max camber as % of c
  # B = location of max camber as % of c
  # CD = max thickness as % of c
# eg. NACA 6412
  # m = 0.06
  # p = 0.4
  # t = 0.12

import matplotlib.pyplot as plt
import numpy as np

def main():
  airfoil = input("Enter NACA Airfoil number: ")

  # error handling
  if (airfoil.isnumeric() == False):
    print("error: enter numeric characters!!!")
    exit()
  if (len(airfoil) != 4):
    print("error: enter 4 digit number!!!")
    exit()

  # plotting equation parameters
  m = float(airfoil[0]) / 100     # A / 100
  p = float(airfoil[1]) / 10      # B / 10
  t = float(airfoil[-2:]) / 100   # CD / 100

  # error handling
  if (m == 0 and p != 0):
    print("error: missing max camber!!!")
    exit()
  if (m != 0 and p == 0):
    print("error: missing location of max camber!!!")
    exit()
  
  num_points = input("Enter number of points desired: ")
  step_size = 1 / (float(num_points) / 2)

  # set up vectors for point file
  x = np.arange(0, 1 + step_size, step_size)
  yu = np.zeros(x.size)
  yl = np.zeros(x.size)

  # abracadabra
  print("Plotting NACA", airfoil, "...")
  plot_airfoil(m, p, t, airfoil, step_size, yu, yl)
  write_points(airfoil, num_points, x, yu, yl)
##### end main #####

##### plotting function #####
def plot_airfoil(m, p, t, airfoil, s, yu, yl):
  title_str = "NACA " + airfoil + " Airfoil"
    
  ##### uncambered airfoil #####
  if (m == 0 and p == 0):
    x = np.arange(0, 1+s, s)
	  # upper surface
    y1 = 5*t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    # lower surface
    y2 = -5*t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    for i in range(yu.size):
      yu[i] = y1[i]
    for j in range(yl.size):
      yl[j] = y2[j]
        
    # plotting
    plt.plot(x, y1, color="r", linewidth=1)
    plt.plot(x, y2, color="r", linewidth=1)
    ax = plt.gca()
    ax.axis('equal')
    ax.set(xlim=(-0.01, 1.01), ylim=(-0.5, 0.5))
    ax.set(xlabel='x/c', ylabel='y/c', title=title_str) 
    ax.grid()
    plt.show()
	  
  ##### cambered airfoil #####
  else: 
	  # x vectors
    x1 = np.arange(0,p,s)
    x2 = np.arange(p,1+s,s)
    xc = np.concatenate([x1, x2])
        
	  # camber line
    yc1 = (m/(p**2)) * (2*p*x1 - x1**2)
    yc2 = (m/(1-p)**2) * ((1-2*p) + 2*p*x2 -x2**2)
    yc = np.concatenate([yc1, yc2])
            
    # transformations
    dydx1 = (2*m/(p**2)) * (p-x1)
    dydx2 = (2*m/(1-p)**2) * (p-x2)
    theta1 = np.arctan(dydx1)
    theta2 = np.arctan(dydx2)
    yt1 = 5*t * (0.2969*np.sqrt(x1) - 0.1260*x1 - 0.3516*x1**2 + 0.2843*x1**3 - 0.1015*x1**4)
    yt2 = 5*t * (0.2969*np.sqrt(x2) - 0.1260*x2 - 0.3516*x2**2 + 0.2843*x2**3 - 0.1015*x2**4)
    xu1 = x1 - yt1*np.sin(theta1)
    xu2 = x2 - yt2*np.sin(theta2)
    xu_plot = np.concatenate([xu1, xu2])
    xl1 = x1 + yt1*np.sin(theta1)
    xl2 = x2 + yt2*np.sin(theta2)
    xl_plot = np.concatenate([xl1, xl2])
            
    # upper surface
    yu1 = yc1 + yt1*np.cos(theta1)
    yu2 = yc2 + yt2*np.cos(theta2)
    yu_plot = np.concatenate([yu1, yu2])
    for i in range(yu_plot.size):
      yu[i] = yu_plot[i]
          
    # lower surface
    yl1 = yc1 - yt1*np.cos(theta1)
    yl2 = yc2 - yt2*np.cos(theta2)
    yl_plot = np.concatenate([yl1, yl2])
    for i in range(yl_plot.size):
      yl[i] = yl_plot[i]
            
    # plotting
    plt.plot(xc, yc, color="b", linewidth=0.75, label="Camber Line")
    plt.plot(xu_plot, yu_plot, color="r", linewidth=1, label="Airfoil Surface")
    plt.plot(xl_plot, yl_plot, color="r", linewidth=1)
    ax = plt.gca()
    ax.axis('equal')
    ax.set(xlim=(0, 1), ylim=(-0.5, 0.5))
    ax.set(xlabel='x/c', ylabel='y/c', title=title_str) 
    ax.grid()
    ax.legend(loc='upper right', shadow=True)
    plt.show()
##### end plotting #####	

##### point writing function #####
def write_points(airfoil, npoints, xvec, yu_vec, yl_vec):
  filename = "NACA" + airfoil + "_" + npoints + ".txt"
  f = open(filename,"w+")

  # reverse lower surface vectors so points loop around
  reversed_x = np.flipud(xvec)
  reversed_yl = np.flipud(yl_vec)

  # upper surface
  for i in range(xvec.size):
    line = str(round(xvec[i], 3)) + " " + str(round(yu_vec[i], 6)) + " \n"
    f.write(line)
  
  # lower surface
  for j in range(xvec.size-2):
    line = str(round(reversed_x[j+1], 3)) + " " + str(round(reversed_yl[j+1], 6)) + " \n"
    f.write(line)
  
  f.close()
##### end point writing #####
    
if __name__ == "__main__":
  main()
    