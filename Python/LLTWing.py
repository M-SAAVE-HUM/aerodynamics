# class defining wing and related functions to solve lifting line theory on wing

import numpy as np
import matplotlib.pyplot as plt

class Wing:

    # constructor - alpha, AR, Nu are required for initialization
    def __init__(self, alpha, AR, Nu=20, **kwargs):
        self.alpha = np.deg2rad(alpha)                      # angle of attack [input deg, use rad]
        self.AR = AR                                        # aspect ratio [-]
        self.Nu = Nu                                        # number of unknowns/vortices [-]
        self.lam = kwargs.get('taper', 1)                   # taper ratio [c_tip/c_root, -]
        self.alphaL0 = np.deg2rad(kwargs.get('alphaL0', 0)) # zero-lift angle of attack [input deg, use rad]
        self.a0 = kwargs.get('a0', 2*np.pi)                 # lift-curve slope [cl/rad]
        self.show_plots = kwargs.get('plots', False)        # option to show lift distribution plots
        self.b = kwargs.get('b', 1)                         # wing span [usually m]
        self.rho = kwargs.get('rho', 1)                     # air density [usually kg/m^3]
        self.Uinf = kwargs.get('Uinf', 1)                   # freestream velocity [usually m/s]

    # LLT function
    def lifting_line(self):

        print('='*60)
        print(f'Solving LLT for alpha={np.rad2deg(self.alpha):.4f} deg, AR={self.AR:.4f} ...')
        print('Solver options: ')
        print(f'    Number unknowns = {self.Nu}')
        print(f'    Taper = {self.lam:.4f}, alphaL0 = {self.alphaL0:.4f}')
        print(f'    a0 = {self.a0:.4f}, show plots = {self.show_plots}')

        M = np.zeros((self.Nu, self.Nu)) # matrix of coefficients for A modes
        F = np.zeros(self.Nu) # RHS vector (sectional angle - 0 lift angle)

        # spanwise collocation points (pi/2 to pi bc left-right symmetry)
        # basically just doing right side of wing (looking down from top)
        # go up to pi - 2/(pi*N) so last point isn't at wingtip
        # Gamma at wingtip is 0, so don't need to solve for it there
        theta = np.linspace(np.pi/2, np.pi - np.pi/(2*self.Nu), self.Nu) 

        # calculate b and c(y) based on lambda and AR
        c_root = 1 # normalized
        c_tip = self.lam
        b = self.AR * (c_root + c_tip) / 2
        c_theta = c_root * (1 - (1 - self.lam) * (-np.cos(theta)))

        # calculate ag / alphaL0 stuff in terms of theta

        # build matrix system
        for i in range(self.Nu): # loop over spanwise locations (thetas)
            F[i] = self.alpha - self.alphaL0
            for j in range(self.Nu): # loop over unknowns (modes)
                # generalized form of PLL equation
                c_i = float(c_theta[i])
                n = 2*j + 1
                M[i, j] = (4*b/(self.a0*c_i)) * np.sin(n * theta[i]) + n * np.sin(n * theta[i]) / np.sin(theta[i])

        # A = M\F = M^-1 * F
        A = np.linalg.solve(M, F)

        # compute coefficients
        A1 = A[0]
        CL = np.pi * self.AR * A1
        delta = 0
        for k in range(1, self.Nu):  # start at k = 1 (higher modes only)
            n = 2*k + 1
            delta += n * (A[k] / A1)**2
        e = 1/(1+delta)
        CDi = CL**2/(np.pi*self.AR*e)

        print('Solved!')
        print('-'*60)
        print('Results:')
        print('-'*60)
        print(f'    CL = {CL:2f}, CDi = {CDi:2f}, e = {e:2f}')
        print('='*60)

        # plot lift distributions
        if self.show_plots == True:
            self.make_plots(A)

        return CL, CDi, e

    def make_plots(self, A):

        theta = np.linspace(0, np.pi, 2*self.Nu) 
        y = -self.b/2 * np.cos(theta)

        # first plot L'(y)
        # L'(theta) = 2*b*rho*Uinf^2 * (A1*sin(theta) + A3*sin(3*theta)) -- divide out the 2 because LL func uses half span
        L_prime = self.b*self.rho*self.Uinf**2 * (A[0]*np.sin(theta) + A[2]*np.sin(3*theta)) # actual 

        # elliptical (ideal)
        L_prime_elliptical = self.b*self.rho*self.Uinf**2 * A[0]*np.sin(theta)

        fig1, ax1 = plt.subplots()
        ax1.plot(y, L_prime, 'b-', label='Actual')
        ax1.plot(y, L_prime_elliptical, 'r--', label='Elliptical')
        ax1.set_xlabel('y')
        ax1.set_ylabel('L\'(y)')
        ax1.legend()

        # then plot cl(y)
        # cl(theta) = 4*b/(c_root * sin(theta)) * (A1*sin(theta) + A3*sin(3*theta))
        c_root = self.AR/(2*self.b*(1+self.lam))
        c_theta = c_root * (1 - (1 - self.lam) * np.abs(np.cos(theta)))
        cl = L_prime / (0.5 * self.rho * self.Uinf**2 * c_theta)
        fig2, ax2 = plt.subplots()
        ax2.plot(y, cl, label=r'$c_l$(y)')
        ax2.set_xlabel('y')
        ax2.set_ylabel(r'$c_l$(y)')
        ax2.legend()
        plt.show()