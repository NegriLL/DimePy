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

    # Initialise RAMBO (rho0 decay)
    if args.pflag == 'rho':
        nparts = 2
        am = np.zeros(4)
        for i in range(4):
            am[i] = 0.13957018
    elif args.pflag == 'phi':
        nparts = 2
        am = np.zeros(4)
        for i in range(4):
            am[i] = 0.493677

    # initialise counters
    nhist = 1
    sum = 0
    sum1 = 0.0
    ncut = 0

    weight = 0.0

    # ToDo: initilise histograms once I know  what kind of values are stored

    num = 0

    # Incoming protons to event record array
    ID = np.zeros(20)
    q = np.zeros((4, 20))
    pup = np.zeros((5, 500))

    ID[0] = 2212
    q[:, 0] = [0, 0, ebeam * beta, ebeam]
    
    meson_parameters["istup"][0] = -1
    meson_parameters["idup"][0] = 2212
    meson_parameters["mothup"][0, 0] = 0 # Some of these lines setting things to zero are irrelevant but I'll keep them here for now
    meson_parameters["mothup"][1, 0] = 0
    meson_parameters["icolup"][0, 0] = 0
    meson_parameters["icolup"][1, 0] = 0
    for i in range(4):
        pup[i, 0] = q[i, 0]
    pup[4,0] = np.sqrt(q[3,0]**2 - q[2, 0]**2 - q[1, 0]**2 - q[0,0]**2)
    meson_parameters["vtimup"][0] = 0
    meson_parameters["spinup"][0] = 9

    q[:, 1] = [0, 0, -ebeam * beta, ebeam]
    meson_parameters["istup"][1] = -1
    if args.ppbar:
        meson_parameters["idup"][1] = -2212
    else:
        meson_parameters["idup"][1] = 221
    meson_parameters["mothup"][0, 1] = 0
    meson_parameters["mothup"][1, 1] = 0
    meson_parameters["icolup"][0, 1] = 0
    meson_parameters["icolup"][1, 1] = 0
    for i in range(4):
        pup[i, 1] = q[i, 1]
    pup[4, 1] = np.sqrt(q[3, 1]**2 - q[2, 1]**2 - q[1, 1]**2 - q[0,1]**2)
    meson_parameters["vtimup"][1] = 0
    meson_parameters["spinup"][1] = 9

    # Outgoing initial info
    meson_parameters['istup'][2] = 1
    if args.ppbar:
        meson_parameters['idup'][1] = -2212
    else:
        meson_parameters['idup'][1] = 2212
    meson_parameters["idup"][2] = 2212
    meson_parameters["mothup"][0, 2] = 1
    meson_parameters["mothup"][1, 2] = 2
    meson_parameters["icolup"][0, 2] = 0
    meson_parameters["icolup"][1, 2] = 0
    meson_parameters["vtimup"][2] = 0
    meson_parameters["spinup"][2] = 9

    meson_parameters["istup"][3] = 1
    meson_parameters["idup"][3] = 2212
    meson_parameters["mothup"][0, 3] = 1
    meson_parameters["mothup"][1, 3] = 2
    meson_parameters["icolup"][0, 3] = 0
    meson_parameters["icolup"][1, 3] = 0
    meson_parameters["vtimup"][3] = 0
    meson_parameters["spinup"][3] = 9

    # HEPEVT
    if args.output == 'hepevt':

        nmxhep = 4000
        isthep = np.zeros(nmxhep, dtype=int)
        idhep = np.zeros(nmxhep, dtype=int)
        jmohep = np.zeros((2, nmxhep), dtype=int)
        jdahep = np.zeros((2, nmxhep), dtype=int)
        phep = np.zeros((5, nmxhep), dtype=float)
        vhep = np.zeros((4, nmxhep), dtype=float)
        
        nhep = meson_parameters['nup']
        
        for k in range(5):
            phep[k, 0] = pup[k, 0]
            phep[k, 1] = pup[k, 1]
        
        for k in range(meson_parameters['nup']):
            isthep[k] = meson_parameters['istup'][k]
            idhep[k] = meson_parameters['idup'][k]
            jmohep[0, k] = meson_parameters['mothup'][0, k]
            jmohep[1, k] = jmohep[0, k]
            jdahep[0, k] = 0
            jdahep[1, k] = 0
            vhep[0, k] = 0.0
            vhep[1, k] = 0.0
            vhep[2, k] = 0.0
            vhep[3, k] = 0.0
        
        if args.pflag == 'rho' or args.pflag == 'phi':
            jdahep[0, 4] = 7
            jdahep[1, 4] = 8
            jdahep[0, 5] = 9
            jdahep[1, 5] = 10
        
        jmohep[1, 4] = 2
        jmohep[1, 5] = 2
    
    # Set lmax based on unw flag
    if args.unw:
        lmax = 2
    else:
        lmax = 1
    
    # Main event generation loop
    for ll in range(1, lmax + 1):
        if ll == 2:
            # ntotal = nev * 10
            ntotal = 1000000000
        
        ip = ntotal + 1
        
        #######################################
        #                                     #
        #       START OF EVENT LOOP           #
        #                                     #
        #######################################
        
        if ll == 1:
            print("Generating weighted events...")
        else:
            print("Generating unweighted events...")
        
        for i in range(1, ntotal + 1):
            weight = 0.0
            
            ran0 = np.random.random()
            ran1 = np.random.random()
            ran2a = np.random.random()
            ran3 = np.random.random()
            ran4 = np.random.random()
            ran5 = np.random.random()
            ran6 = np.random.random()
            
            ID[0] = 2212
            q[:, 0] = [0, 0, ebeam, ebeam]
            
            ID[1] = 2212
            q[:, 1] = [0, 0, -ebeam, ebeam]
    
            phi1 = 2.0 * pi * ran0
            phi2 = 2.0 * pi * ran1
            
            pt1sq = -np.log(ran2a) / bjac
            pt2sq = -np.log(ran3) / bjac
            
            weight = (np.exp(bjac * pt1sq) * np.exp(bjac * pt2sq)) / bjac**2
            
            pt1 = np.array([np.sqrt(pt1sq) * np.sin(phi1), np.sqrt(pt1sq) * np.cos(phi1)])
            pt2 = np.array([np.sqrt(pt2sq) * np.sin(phi2), np.sqrt(pt2sq) * np.cos(phi2)])
            
            pt1x = pt1[0]
            pt1y = pt1[1]
            pt2x = pt2[0]
            pt2y = pt2[1]
            
            if args.pflag == 'rho':
                while True:
                    rm1 = np.random.random()
                    rm2 = np.random.random()
                    
                    mmes0 = meson_parameters['mmes0']
                    mwidth = meson_parameters['mwidth']
                    msmax = mmes0 + 4.0 * mwidth
                    msmin = 2.0 * 0.13957018
                    
                    almin = np.arctan(-(mmes0**2 - msmin**2) / mwidth / mmes0)
                    almax = np.arctan(-(mmes0**2 - msmax**2) / mwidth / mmes0)
                    
                    al1 = almin + (almax - almin) * rm1
                    al2 = almin + (almax - almin) * rm2
                    
                    mmes1 = np.sqrt(np.tan(al1) * mmes0 * mwidth + mmes0**2)
                    mmes2 = np.sqrt(np.tan(al2) * mmes0 * mwidth + mmes0**2)
                    
                    weight = weight * (almax - almin)**2
                    weight = weight * mwidth**2 * mmes0**2
                    weight = weight * (1.0 + np.tan(al1)**2) * (1.0 + np.tan(al2)**2)
                    weight = weight / 4.0 / mmes1 / mmes2
                    
                    mwidth1 = mwidth * ((1.0 - 4.0 * 0.13957018**2 / mmes1**2) /
                                       (1.0 - 4.0 * 0.13957018**2 / mmes0**2))**1.5
                    mwidth2 = mwidth * ((1.0 - 4.0 * 0.13957018**2 / mmes2**2) /
                                       (1.0 - 4.0 * 0.13957018**2 / mmes0**2))**1.5
                    
                    if mmes1 >= 2.0 * 0.13957018 and mmes2 >= 2.0 * 0.13957018:
                        break
                
                weight = weight * 2.0 * mmes0 * mmes1 * mwidth1 / pi
                weight = weight * 2.0 * mmes0 * mmes2 * mwidth2 / pi
                weight = weight * mmes1**2 * mmes2**2 / mmes0**4
                weight = weight / ((mmes0**2 - mmes1**2)**2 + mmes1**4 * mwidth1**2 / mmes0**2)
                weight = weight / ((mmes0**2 - mmes2**2)**2 + mmes2**4 * mwidth2**2 / mmes0**2)
            
            elif args.pflag == 'phi':
                while True:
                    rm1 = np.random.random()
                    rm2 = np.random.random()
                    
                    mmes0 = meson_parameters['mmes0']
                    mwidth = meson_parameters['mwidth']
                    msmax = mmes0 + 4.0 * mwidth
                    msmin = mmes0 - 4.0 * mwidth
                    
                    almin = np.arctan(-(mmes0**2 - msmin**2) / mwidth / mmes0)
                    almax = np.arctan(-(mmes0**2 - msmax**2) / mwidth / mmes0)
                    
                    al1 = almin + (almax - almin) * rm1
                    al2 = almin + (almax - almin) * rm2
                    
                    mmes1 = np.sqrt(np.tan(al1) * mmes0 * mwidth + mmes0**2)
                    mmes2 = np.sqrt(np.tan(al2) * mmes0 * mwidth + mmes0**2)
                    
                    weight = weight * (almax - almin)**2
                    weight = weight * mwidth**2 * mmes0**2
                    weight = weight * (1.0 + np.tan(al1)**2) * (1.0 + np.tan(al2)**2)
                    weight = weight / 4.0 / mmes1 / mmes2
                    
                    mwidth1 = mwidth * ((1.0 - 4.0 * 0.493677**2 / mmes1**2) /
                                       (1.0 - 4.0 * 0.493677**2 / mmes0**2))**1.5
                    mwidth2 = mwidth * ((1.0 - 4.0 * 0.493677**2 / mmes2**2) /
                                       (1.0 - 4.0 * 0.493677**2 / mmes0**2))**1.5
                    
                    if mmes1 >= 2.0 * 0.493677 and mmes2 >= 2.0 * 0.493677:
                        break
                
                weight = weight * 2.0 * mmes0 * mmes1 * mwidth1 / pi
                weight = weight * 2.0 * mmes0 * mmes2 * mwidth2 / pi
                weight = weight * mmes1**2 * mmes2**2 / mmes0**4
                weight = weight / ((mmes0**2 - mmes1**2)**2 + mmes1**4 * mwidth1**2 / mmes0**2)
                weight = weight / ((mmes0**2 - mmes2**2)**2 + mmes2**4 * mwidth2**2 / mmes0**2)
            
            # pi rapidities
            ranx1 = np.random.random()
            ranx2 = np.random.random()
            ranx3 = np.random.random()
            ranx4 = np.random.random()
            
            phix1 = 2.0 * pi * ranx1
            
            ptxsqmin = 0.0
            ptxsqmax = 10.0  # generate linear p_t^2
            
            ptxsq1 = ptxsqmin + (ptxsqmax - ptxsqmin) * ranx2
            
            q[0, 4] = np.sqrt(ptxsq1) * np.cos(phix1)
            q[1, 4] = np.sqrt(ptxsq1) * np.sin(phix1)
            q[0, 5] = -pt1[0] - pt2[0] - q[0, 4]
            q[1, 5] = -pt1[1] - pt2[1] - q[1, 4]
            ptxsq2 = q[0, 5]**2 + q[1, 5]**2
            
            rmx1 = np.sqrt(mmes1**2 + ptxsq1)
            rmx2 = np.sqrt(mmes2**2 + ptxsq2)
            
            ymax1 = np.log(rts / rmx1)
            
            if args.cuts == 'true':
                if rmax < ymax1:
                    ymax1 = rmax
            
            ymin1 = -ymax1
            
            yx1 = ymin1 + (ymax1 - ymin1) * ranx3
            
            ymax2 = np.log((rts - rmx1 * np.exp(yx1)) / rmx2)
            ymin2 = -np.log((rts - rmx1 * np.exp(-yx1)) / rmx2)
            
            if args.cuts == 'true':
                if ymax2 > rmax:
                    ymax2 = rmax
                if ymin2 < -rmax:
                    ymin2 = -rmax
            
            if ymax2 < ymin2:
                weight = 0.0
                sum += weight
                sum1 += weight**2
                continue
            
            yx2 = ymin2 + (ymax2 - ymin2) * ranx4
            x1 = (rmx2 * np.exp(yx2) + rmx1 * np.exp(yx1)) / rts
            x2 = (rmx2 * np.exp(-yx2) + rmx1 * np.exp(-yx1)) / rts
            
            weight = weight * (ptxsqmax - ptxsqmin)
            weight = weight * (ymax1 - ymin1) * (ymax2 - ymin2)
            
            q[2, 4] = rmx1 * (np.exp(yx1) - np.exp(-yx1)) / 2.0
            q[3, 4] = rmx1 * (np.exp(yx1) + np.exp(-yx1)) / 2.0
            
            q[2, 5] = rmx2 * (np.exp(yx2) - np.exp(-yx2)) / 2.0
            q[3, 5] = rmx2 * (np.exp(yx2) + np.exp(-yx2)) / 2.0
            
            # Impose massive on-shell condition
            # p1+ + cc1/p2- = aa1
            # p2- + cc2/p1+ = aa2 
            aa1 = bp * (1.0 - x1)
            aa2 = bp * (1.0 - x2)
            cc1 = 0.5 * (pt2sq + mp**2)
            cc2 = 0.5 * (pt1sq + mp**2)
            
            root1sq = (cc1 - cc2 - aa1 * aa2)**2 - 4.0 * cc2 * aa1 * aa2
            root2sq = (cc2 - cc1 - aa1 * aa2)**2 - 4.0 * cc1 * aa1 * aa2
            
            if root1sq <= 0.0 or root2sq <= 0.0:
                weight = 0.0
                # line 1123 for the GOTO 700 then begin loop again
                sum += weight
                sum1 += weight**2
                continue
            
            p1p = (cc2 - cc1 + aa1 * aa2 + np.sqrt(root1sq)) / (2.0 * aa2)
            p2m = (cc1 - cc2 + aa1 * aa2 + np.sqrt(root2sq)) / (2.0 * aa1)
            p1m = (pt1sq + mp**2) / (2.0 * p1p)
            p2p = (pt2sq + mp**2) / (2.0 * p2m)
            
            if p1p < 0.0 or p1m < 0.0 or p2p < 0.0 or p2m < 0.0:
                weight = 0.0
                sum += weight
                sum1 += weight**2
                continue
            
            t1 = -rts * p1m * rt2
            t2 = -rts * p2p * rt2
            
            q[:, 2] = [pt1[0], pt1[1], (p1p - p1m) / rt2, (p1p + p1m) / rt2]
            q[:, 3] = [pt2[0], pt2[1], (p2p - p2m) / rt2, (p2p + p2m) / rt2]
            
            svec = q[:, 4] + q[:, 5]
            sh = np.dot(svec, svec)
            
            # TODO: Implement from line 854 onwards
            
    

main()