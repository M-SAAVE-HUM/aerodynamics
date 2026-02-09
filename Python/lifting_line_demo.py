# Wing class with lifting line solver demo

import numpy as np
import matplotlib.pyplot as plt
import LLTWing

def ar_sweep_rect(ARs):
    Na = 40 # number of alphas to run
    alphas = np.linspace(-10, 10, Na)

    for AR in ARs:
        CLs = np.zeros(Na)
        CDis = np.zeros(Na)
        for i in range(Na):
            a = alphas[i]
            Nu = 40 # num unknowns
            wing_i = LLTWing.Wing(a, AR, Nu)
            CLs[i], CDis[i], _ = wing_i.lifting_line()

        # plot drag polar for that AR
        plt.plot(CDis, CLs, label=f'AR = {AR}')

    plt.xlabel(r'$C_{Di}$')
    plt.ylabel(r'$C_L$')
    plt.legend()
    plt.grid(True, alpha=0.25)
    plt.show()

def taper_sweep(lams, ARs):
    Nt = np.size(lams) # number of ratios to try

    for AR in ARs:
        e = np.zeros(Nt)

        for i in range(Nt):
            lam = lams[i]
            alpha = np.deg2rad(4.0) # at 4 deg aoa
            Nu = 40 # num lifting line unknowns
            wing_i = LLTWing.Wing(alpha, AR, Nu, taper=lam)
            _, _, e[i] = wing_i.lifting_line()

        # plot e vs taper ratio
        plt.plot(lams, e, label=f'AR = {AR}')
    
    plt.xlabel(r'Taper ratio, $\lambda$')
    plt.ylabel('e')
    plt.legend()
    plt.grid(True, alpha=0.25)
    plt.show()


if __name__ == "__main__":

    # example wing with lift distributions plotted
    # inputs: angle of attack, aspect ratio, number of unknowns, other options
    test_wing_opt = LLTWing.Wing(6.25, 7.87, 40, plots=True, rho=1.225, Uinf=20.22, b=2.5)
    test_wing_opt.lifting_line()

    test_wing_1 = LLTWing.Wing(6.25, 7.87, 40, plots=True, rho=1.225, Uinf=20.22, b=2.5)
    test_wing_1.lifting_line()

    # # sweep through aspect ratios for rectangular wings
    # ARs = np.array([4, 6, 8, 10])
    # ar_sweep_rect(ARs)

    # # sweep through taper ratios for various ARs, compare e
    # lams = np.linspace(0.05, 1, 50)
    # taper_sweep(lams, ARs)