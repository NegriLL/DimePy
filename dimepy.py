import numpy as np
from scipy.special import j0

from parser import parse_arguments, parse_meson_parameters

# Initialize global variables
nch = 0
cc0 = np.zeros(5, dtype=float)
bm = np.zeros(5, dtype=float)
bb0 = np.zeros(5, dtype=float)
pp0 = np.zeros(5, dtype=float)
bex = np.zeros(5, dtype=float)
asp = 0.0
sigo = 0.0
gaa = np.zeros(5, dtype=float)
ep = 0.0
norm = 0.0
rts = 0.0

# Initialize opacity arrays
op = np.zeros((5, 5, 10000, 2), dtype=float)
oph = np.zeros((5, 5, 10000, 2), dtype=float)

def initpars(iin, rts):
    global cc0, bm, bb0, pp0, bex, asp, sigo, gaa, ep, norm, nch
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
        raise ValueError(f"Error initiating parameters. Unknown iin value {iin}. Value iin must be 1, 2, 3 or 4")

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

    for i in range(nch):
        for j in range(nch):
            sum = sum + gaa[i] * gaa[j] * pp0[i] * pp0[j]/np.pow(nch, 2)

    norm = sum

def calcscreen():
    sca = np.zeros((5, 5, 40001, 2), dtype=float)
    sca1 = np.zeros((5, 5, 40001, 2), dtype=float)

    ns = 900
    ksqma = 8.2
    inck = ksqma / ns
    ksqmin = 0.001
    lginck = np.log(ksqma/ksqmin) / ns

    print('Calculating screening amplitude')

    for ib in range(ns + 2): 
        if ib == 0:
            ksq = 0.0
            lgksq = 0.0
        else:
            lgksq = (ib - 1) * lginck + np.log(ksqmin)
            ksq = np.exp(lgksq)
        
        for i in range(nch):
            for j in range(nch):
                sc, sc1 = screening(i, j, ksq)
                sca[i, j, ib, 0] = lgksq
                sca[i, j, ib, 1] = sc
                sca1[i, j, ib, 0] = lgksq
                sca1[i, j, ib, 1] = sc1

def screening(i, j, ktsq):
    nb = 5000
    hb = 99.0 / nb

    sc = 0
    sc1 = 0

    for ib in range(1, nb + 1):
        bt = (ib - 1) * hb
        wt = -bt / 2.0 / np.pi * hb

        fr, fr1 = opacityint(i, j, bt)

        sige = sigo * np.exp(np.log(rts) * 2 * ep)
        fr = fr * gaa[i] * gaa[j] * sige
        fr1 = fr1 * gaa[i] * gaa[j] * sige

        sc = sc + wt * (1 - np.exp(-fr/2)) * j0(bt * np.sqrt(ktsq)) * gaa[i] * gaa[j]
        sc1 = sc1 + wt * (1 - np.exp(-fr1/2)) * j0(bt * np.sqrt(ktsq)) * gaa[i] * gaa[j]

    return sc, sc1

def opacityint(i, j, bt):
    global op, oph
    incbt = op[0, 0, 1, 0] - op[0, 0, 0, 0]
    it = int(np.floor(bt/incbt))

    m = (op[i, j, it + 1, 1] - op[i, j, it, 1])/(op[i, j, it + 1, 0] - op[i, j, it, 0])
    delta = bt - op[0, 0, it, 0]
    mh = (oph[i, j, it + 1, 1] - oph[i, j, it, 1])/(oph[i, j, it + 1, 0] - oph[i, j, it, 0])
    deltah = bt - oph[0, 0, it, 0]

    fr = m * delta + op[i, j, it, 1]
    fr1 = mh * deltah + oph[i, j, it, 1]

    return fr, fr1


def calcop():
    global op, oph
    nb = 900
    hb = 100.0/nb

    print("Calculating Opacity")
    with open('output.dat', 'w') as file:

        for ib in range(nb + 1):
            bt = ib * hb

            for i in range(nch):
                for j in range(nch):
                    fr, fr1 = opacity(i, j, bt)

                    op[i, j, ib, 0] = bt
                    op[i, j, ib, 1] = fr
                    oph[i, j, ib, 0] = bt
                    oph[i, j, ib, 1] = fr1
                    file.write(f"{bt:24.16f}  {fr:.16e}  {fr1:.16e}\n")


def opacity(i, j, bt):
    ampi = 0.02
    amro = 1
    a4 = 4 * ampi
    alo = np.log(amro/ampi)
    bpol = 2.4

    nt = 6000
    htt = 6/np.pow(nt, 2)

    fr = 0
    fr1 = 0

    for it in range(nt + 1):
        t = np.pow(it, 2) * htt
        wt = htt * 2 * it/4/np.pi
        if it == 0:
            t = 1e-8
            wt = wt/2
        bes0 = j0(bt*np.sqrt(t))

        ffi = np.exp(-np.pow((t+0.08+bb0[i])*bex[i],cc0[i])
                    +np.pow(bex[i]*(bb0[i]+0.08),cc0[i]))
        ffj = np.exp(-np.pow((t+0.08+bb0[j])*bex[j],cc0[j])
                    +np.pow(bex[j]*(bb0[j]+0.08),cc0[j]))
        
        asp1 = asp
        form1 = np.log(ffi*ffj) - 2 * t * asp1 * np.log(rts)

        ah = np.sqrt(1+a4/t)
        h1pi=2*a4+t*(alo - pow(ah, 3) * np.log((1+ah)/(ah-1))) # Pion loop insertion
        h1pi=h1pi * sigo/(72 * np.pow(np.pi,3) * np.pow((1+t/bpol),2))

        ww = bes0 * np.exp(form1-2*h1pi*np.log(rts))
        aspt=t*asp+h1pi

        fr = fr + ww * wt
        fr1 = fr1 + bes0 * wt * ffi * ffj 

    return fr, fr1

def main():
    global rts
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
    calcscreen()

    meson_parameters = parse_meson_parameters(args.pflag)

    # Other parameters
    ebeam = rts / 2.0
    s = rts**2
    zi = 0 + 1j
    rt2 = np.sqrt(2.0)
    pi = np.arccos(-1.0)
    bp = rts / np.sqrt(2.0)
    mp = 0.93827
    beta = np.sqrt(1.0 - 4.0 * mp**2 / s)
    s0 = 1.0

    # Pomeron + t-slope
    bb = 4.0
    bjac = 6.0
    bjac1 = 2.0

    alphap = 0.25    # D-L 92 fit
    alpha0 = 1.0808
    alphapr = 0.93
    alpha0r = 0.55

    mf127 = 1.275
    mf1525 = 1.525

    cpom = meson_parameters.get("sig0") / 0.389
    aff = -0.860895
    ar = -1.16158

    # Fortran code assigns this in a weird way (lines 406 and 407)
    # where mmes might not be initialized if pflag is set to rho or phi
    # and then it reassigns it later. This should be equivalent in python 
    # while avoiding uninitiated parameters
    mmes1 = mmes2 = meson_parameters.get('mmes') or meson_parameters.get('mmes0')



main()