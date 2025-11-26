import numpy as np

from parser import parse_arguments

# ToDo: once program is completed we need to reevaluate the use of global variables.
global cc0, bm, bb0, pp0, bex, asp, sigo, gaa, ep, norm # fortran block pars
global nch                                              # fortrant block ipars
global s, rts, mmes, yx                                 # fortrant block vars

def initpars(iin, rts):
    i1 = 0
    i2 = 0
    nch = 0

    if iin == 1:
        ep = 0.13        # Capital delta
        asp = 0.08       # alpha'
        ep1 = 0          # Zero in all models, matters (?)
        sigo = 23.0/0.39 # sigma_0 (GeV^-2)
        gd2 = 0.3        # k1/k (1.8 TeV)
        nch = 2   
        ntf = 0

        gaa = np.array([0.00, 0.00, 0.4, 0.0, 0.0], dtype=float)
        cc0 = np.array([0.45, 0.45, 1.0, 0.0, 0.0], dtype=float) # d1, d2
        bm  = np.array([3.00, 1.50, 0.0, 0.0, 0.0], dtype=float) # Doesn't matter
        bb0 = np.array([0.10, 0.50, 0.8, 0.0, 0.0], dtype=float) # c1-0.08 (added back later)
        pp0 = np.array([0.92, 0.10, 0.5, 0.0, 0.0], dtype=float) # 2*|a_1|^2
        bex = np.array([8.50, 4.50, 0.5, 0.0, 0.0], dtype=float) # b1

    elif iin == 2:
        ep = 0.115
        asp = 0.11
        ep1 = 0
        sigo = 33.0/0.39
        gd2 = 0.16
        nch = 2
        ntf = 0

        gaa = np.array([0.00, 0.00, 0.6, 0.0, 0.0], dtype=float)
        cc0 = np.array([0.63, 0.47, 1.0, 0.0, 0.0], dtype=float)
        bm  = np.array([3.00, 1.50, 0.0, 0.0, 0.0], dtype=float)
        bb0 = np.array([0.10, 0.50, 0.8, 0.0, 0.0], dtype=float)
        pp0 = np.array([0.50, 0.10, 0.5, 0.0, 0.0], dtype=float)
        bex = np.array([8.00, 6.00, 0.5, 0.0, 0.0], dtype=float)

    elif iin == 3:
        ep = 0.093
        asp = 0.075
        ep1 = 0
        sigo = 60.0/0.39
        gaa3 = 4.8       # k2/k (1.8 TeV)
        gd2 = 1.03
        nch = 2
        nga = 1          # A flag?

        cc0 = np.array([0.55, 0.48, 0.24, 0.0, 0.0], dtype=float)
        bm  = np.array([3.00, 1.50, 0.00, 0.0, 0.0], dtype=float)
        bb0 = np.array([0.27, 0.10, 0.00, 0.0, 0.0], dtype=float)
        pp0 = np.array([0.48, 1.00, 0.00, 0.0, 0.0], dtype=float)
        bex = np.array([5.30, 3.80, 0.00, 0.0, 0.0], dtype=float)

    elif iin == 4:
        ep = 0.093
        asp = 0.075
        ep1 = 0
        sigo = 60.0/0.39
        gaa3 = 4.8
        gd2 = 1.03
        nch = 2
        nga = 1

        cc0 = np.array([0.55, 0.48, 0.24, 0.0, 0.0], dtype=float) # d1, d2, beta : k^2_min ~ s^Beta
        bm  = np.array([3.00, 1.50, 0.00, 0.0, 0.0], dtype=float)
        bb0 = np.array([0.27, 0.10, 0.00, 0.0, 0.0], dtype=float)
        pp0 = np.array([0.48, 1.00, 0.00, 0.0, 0.0], dtype=float)
        bex = np.array([5.30, 3.80, 0.00, 0.0, 0.0], dtype=float)

    else:
        raise Exception("Error initiating parameters. Value iin must be 1, 2, 3 or 4")

    if nch == 3:
        pp0[2] = 3.0 - pp0[1] - pp0[0]
    if nch == 2:
        pp0[1] = 2.0 - pp0[0] # Set |a_2|^2
    if nch == 1:
        pp0[0] = 1.0

    if iin == 3 or iin == 4:
        gamm = np.pow((1800.0 / rts), cc0[2])
        ga1 = 1.0/(1.0 + gamm * gd2)
        ga2 = 1.0/(1 + gamm * gaa3)
        gaa[0] = (2.0 * ga1)/(ga1 + ga2)
        gaa[1] = (2.0 * ga2)/(ga1 + ga2)

    elif iin == 1 or iin == 2:
        if nch == 2:
            gaa[0] = 1.0 + np.sqrt(gd2)
            gaa[1] = 1.0 - np.sqrt(gd2)
            gaa[3] = 0.0
        elif nch == 1:
            gaa[0] = 1.0
            gaa[2] = 0.0
            gaa[3] = 0.0
    
    sum = 0.0

    for i in range(nch):  # fortrant starts arrays at 1 so this should be equivalent
        for j in range(nch):
            sum = sum + gaa[i] + gaa[j] * pp0[1] * pp0[j]/np.pow(nch, 2)

    norm = sum


def calcop():
    sca = np.zeros((5, 5, 40001, 2), dtype=float)
    sca1 = np.zeros((5, 5, 40001, 2), dtype=float)

    ns = 900
    ksqma = 8.2
    inck = ksqma / ns
    ksqmin = 0.001
    lginck = np.log(ksqma/ksqmin) / ns

    print('Calculating screening amplitude')

    for ib in range(0, ns + 2):  # we want ib to include 0 and go up to ns + 1 so we need ns + 2 to be equivalent to fortran
        ksq = (ib - 1) * inck
        lgksq = ksq
        if ib == 0:
            ksq = 0
            lgksq = 0.0
        else:
            lgksq = (ib - 1) * lginck + np.log(ksqmin)
            ksq = dexp(lgksq)
        
        for i in range(nch):
            for j  in range(nch):
                sc, sc1 = screening(i, j, ksq)
                sca[i, j, ib, 0] = lgksq
                sca[i, j, ib, 1] = sc
                sca1[i, j, ib, 0] = lgksq
                sca1[i, j, ib, 1] = sc1



# ToDo: implement dexp
def dexp(lgksq):
    pass


# ToDo: implement screening
def screening(i, j, ksq):
    pass


def main():

    rts=13.0e3  # Centre of mass Energy

    # Some basic cuts, imposed in subtroutine 'icut'. Other user defined cuts can readily be implemented in subroutine
    # note: in rhorho case these cuts are imposed on rho's, and not their decay productions. Cuts on the decay products
    # can be imposed by editing 'icut'

    rmax = 2.0  # Max meson rapidity
    rmin = -2.0 # Min meson rapidity
    ecut= 0.0   # Min meson p_t

    args = parse_arguments() # parsing command line arguments

    ntotal=10000 # no. of runs for weighted events
    nev=200 # no. of unweighted events generated to event record

    # Set parameters for Pomeron-Meson form factor in function fpi
    if(args.formf == 'exp'):
        bexp = 1/2.2
    elif(args.formf == 'orexp'):
        bo = 1/1.1
        ao = np.sqrt(0.5)
    elif(args.formf == 'power'):
        aoo = 1.7
    
    ########################################
    #          Start of main body          #
    ########################################

    initpars(args.iin, rts)
    calcop()

main()