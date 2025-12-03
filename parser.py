import numpy as np

import argparse
import textwrap

def parse_arguments():
    parser = argparse.ArgumentParser(prog='DimeMC Python', usage=preamble(), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pflag",  
                        nargs='?',
                        default='phi',
                        type=str,
                        help='Process generated - see preamble for options')
    parser.add_argument("--fsi",
                        nargs='?',
                        default=True,
                        type=bool,
                        help='phenomenological model for exclusive suppression in Pom Pom --> meson pair. To turn on/off --> "true/false"')
    parser.add_argument("--formf",
                        nargs='?',
                        default='orexp',
                        type=str,
                        help='meson - pomeron form factor.',
                        choices=['exp', 'orexp', 'power'])
    parser.add_argument("--ppbar",
                        nargs='?',
                        default=False,
                        type=bool,
                        help='set true if ppbar collisions')
    parser.add_argument("--output",
                        nargs='?',
                        default='hepevt',
                        type=str,
                        help='Event record style, HEPEVT or LHE')
    parser.add_argument("--cuts",
                        nargs='?',
                        default=True,
                        type=bool,
                        help='Impose cuts or not')
    parser.add_argument("--unw",
                        nargs='?',
                        default=True,
                        type=bool,
                        help='Set true for unweighted events')
    parser.add_argument("--iin",
                        nargs='?',
                        default=1,
                        type=int,
                        choices=[1, 2, 3, 4],
                        help='Model for soft survival factor, as described in arXiv:1306.2149.')
    parser.add_argument("--seed",
                        nargs='?',
                        default=1,
                        type=int,
                        help='Random number seed')

    args = parser.parse_args()
    return args


def parse_meson_parameters(pflag):
    params = {}

    # Initialize event record arrays (MAXNUP=500 in Fortran)
    params['istup'] = np.zeros(500, dtype=int)
    params['idup'] = np.zeros(500, dtype=int)
    params['mothup'] = np.zeros((2, 500), dtype=int)
    params['icolup'] = np.zeros((2, 500), dtype=int)
    params['vtimup'] = np.zeros(500, dtype=float)
    params['spinup'] = np.zeros(500, dtype=float)
    
    if pflag == 'pipm':
        params['mmes'] = 0.13957018 # pi+/- mass, PDG 2011 value
        params['sig0'] = 13.63
        params['alphapm'] = 0.7
        params['alpha0m'] = -0.7 * params['mmes']**2
        params['cf'] = 31.79 / 0.389
        params['crho'] = 4.23 / 0.389

        params['nup'] = 6
        params['istup'][4] = 1
        params['idup'][4] = 211
        params['mothup'][0, 4] = 1
        params['mothup'][1, 4] = 2
        params['icolup'][0, 4] = 0
        params['icolup'][1, 4] = 0
        params['vtimup'][4] = 0
        params['spinup'][4] = 9.

        params['istup'][5] = 1
        params['idup'][5] = -211
        params['mothup'][0, 5] = 1
        params['mothup'][1, 5] = 2
        params['icolup'][0, 5] = 0
        params['icolup'][1, 5] = 0
        params['vtimup'][5] = 0
        params['spinup'][5] = 9.

        
    elif pflag == 'pi0':
        params['mmes'] = 0.1349766  # pi0 mass, PDG 2011 value
        params['sig0'] = 13.63
        params['alphapm'] = 0.7
        params['alpha0m'] = -0.7 * params['mmes']**2
        params['cf'] = 31.79 / 0.389
        params['crho'] = 0.0

        params['nup'] = 6
        params['istup'][4] = 1
        params['idup'][4] = 111
        params['mothup'][0, 4] = 1
        params['mothup'][1, 4] = 2
        params['icolup'][0, 4] = 0
        params['icolup'][1, 4] = 0
        params['vtimup'][4] = 0
        params['spinup'][4] = 9.

        params['istup'][5] = 1
        params['idup'][5] = -111
        params['mothup'][0, 5] = 1
        params['mothup'][1, 5] = 2
        params['icolup'][0, 5] = 0
        params['icolup'][1, 5] = 0
        params['vtimup'][5] = 0
        params['spinup'][5] = 9.
        
    elif pflag == 'kpkm':
        params['mmes'] = 0.493677  # K+/- mass, PDG 2011 value
        params['sig0'] = 11.82
        params['alphapm'] = None   # Not set in Fortran for kaons
        params['alpha0m'] = None
        params['cf'] = 17.255 / 0.389
        params['crho'] = 9.105 / 0.389

        params['nup'] = 6
        params['istup'][4] = 1
        params['idup'][4] = 321
        params['mothup'][0, 4] = 1
        params['mothup'][1, 4] = 2
        params['icolup'][0, 4] = 0
        params['icolup'][1, 4] = 0
        params['vtimup'][4] = 0
        params['spinup'][4] = 9.

        params['istup'][5] = 1
        params['idup'][5] = -321
        params['mothup'][0, 5] = 1
        params['mothup'][1, 5] = 2
        params['icolup'][0, 5] = 0
        params['icolup'][1, 5] = 0
        params['vtimup'][5] = 0
        params['spinup'][5] = 9.
        
    elif pflag == 'ks':
        params['mmes'] = 0.497614  # K_0 mass, PDG 2011 value
        params['sig0'] = 11.82
        params['alphapm'] = None
        params['alpha0m'] = None
        params['cf'] = 17.255 / 0.389
        params['crho'] = 0.0

        params['nup'] = 6
        params['istup'][4] = 1
        params['idup'][4] = 310
        params['mothup'][0, 4] = 1
        params['mothup'][1, 4] = 2
        params['icolup'][0, 4] = 0
        params['icolup'][1, 4] = 0
        params['vtimup'][4] = 0
        params['spinup'][4] = 9.

        params['istup'][5] = 1
        params['idup'][5] = 310
        params['mothup'][0, 5] = 1
        params['mothup'][1, 5] = 2
        params['icolup'][0, 5] = 0
        params['icolup'][1, 5] = 0
        params['vtimup'][5] = 0
        params['spinup'][5] = 9.
        
    elif pflag == 'rho':
        params['mmes0'] = 0.77549  # rho mass, PDG 2013 value
        params['mwidth'] = 0.1491  # rho width, PDG 2013 value
        params['sig0'] = 10.0
        params['alphapm'] = None
        params['alpha0m'] = None
        params['cf'] = 0.0
        params['crho'] = 0.0

        params['nup'] = 10
        params['istup'][4] = 2
        params['idup'][4] = 113
        params['mothup'][0, 4] = 1
        params['mothup'][1, 4] = 2
        params['icolup'][0, 4] = 0
        params['icolup'][1, 4] = 0
        params['vtimup'][4] = 0
        params['spinup'][4] = 9.

        params['istup'][5] = 2
        params['idup'][5] = 113
        params['mothup'][0, 5] = 1
        params['mothup'][1, 5] = 2
        params['icolup'][0, 5] = 0
        params['icolup'][1, 5] = 0
        params['vtimup'][5] = 0
        params['spinup'][5] = 9.

        for k in range(6, 10):
            params['istup'][k] = 1
            params['mothup'][1, k] = 0
            params['icolup'][0, k] = 0
            params['icolup'][1, k] = 0
            params['vtimup'][k] = 0
            params['spinup'][k] = 9.

        params['idup'][6] = 211
        params['idup'][7] = -211
        params['idup'][8] = 211
        params['idup'][9] = -211
        params['mothup'][0, 6] = 5
        params['mothup'][0, 7] = 5
        params['mothup'][0, 8] = 6
        params['mothup'][0, 9] = 6
        
    elif pflag == 'phi':
        params['mmes0'] = 1.01946  # phi mass, PDG 2013 value
        params['mwidth'] = 0.00426 # phi width, PDG 2013 value
        params['sig0'] = 5.7       # Scale by (m_rho/m_phi)^2 for phi pair production
        params['alphapm'] = None
        params['alpha0m'] = None
        params['cf'] = 0.0
        params['crho'] = 0.0

        params['nup'] = 10
        params['istup'][4] = 2
        params['idup'][4] = 333
        params['mothup'][0, 4] = 1
        params['mothup'][1, 4] = 2
        params['icolup'][0, 4] = 0
        params['icolup'][1, 4] = 0
        params['vtimup'][4] = 0
        params['spinup'][4] = 9.

        params['istup'][5] = 2
        params['idup'][5] = 333
        params['mothup'][0, 5] = 1
        params['mothup'][1, 5] = 2
        params['icolup'][0, 5] = 0
        params['icolup'][1, 5] = 0
        params['vtimup'][5] = 0
        params['spinup'][5] = 9.

        for k in range(6, 10):
            params['istup'][k] = 1
            params['mothup'][1, k] = 0
            params['icolup'][0, k] = 0
            params['icolup'][1, k] = 0
            params['vtimup'][k] = 0
            params['spinup'][k] = 9.

        params['idup'][6] = 321
        params['idup'][7] = -321
        params['idup'][8] = 321
        params['idup'][9] = -321
        params['mothup'][0, 6] = 5
        params['mothup'][0, 7] = 5
        params['mothup'][0, 8] = 6
        params['mothup'][0, 9] = 6
        
    else:
        raise ValueError(f"Unknown pflag value {pflag}. pflag must be pipm, pi0, kpkm, ks, rho or phi")
    
    return params


def preamble():
    return textwrap.dedent('''
###########################################################
#                                                         #
#     Dime MC for the central exclusive production of     #
#     meson pairs by double Pomeron exchange:             #
#                                                         #
#     p(1) p(2) --> p(3) + M_1(5) M_2(6) + p(4)           #
#                                                         #
#     Momenta for each event in array q(i,j), where j is  #
#     the particle label and i is the 4-momentum          #
#     component, with:                                    #
#                                                         #
#     1,2 = transverse components                         #
#     3   = beam component                                # 
#     4   = energy                                        #
#                                                         #
#     PDG number for ith particle in arrary ID(i)         #
#                                                         #
#     Also gives event record in HEPEVT or LHE format     #
#     (others are available upon request)                 #
#                                                         #
###########################################################
#                                                         #
#     Particles generated:                                #
#                                                         #
#     pflag='pipm'  :  M_1=pi^+ M_2=pi^-                  #
#     pflag='pi0'   :  M_1=M_2=pi^0                       #
#     pflag='kpkm'  :  M_1=K^+  M_2=K^-                   #
#     pflag='ks'    :  M_1=M_2=K_0                        #
#     pflag='rho'   :  M_1=M_2=rho_0                      #
#     pflag='phi'   :  M_1=M_2=phi                        #
#                                                         #
#     with decay: rho(5) --> pi^+(7)pi^-(8)               #
#                 rho(6) --> pi^+(9)pi^-(10)              #
#     according to phase space, with finite rho           #
#     width included. Similarly for phi to Kaon decay.    #
#                                                         #
###########################################################
#                                                         #
#     User defined cuts can readily be implemented        #
#     in subroutine 'icut'                                #
#                                                         #
###########################################################
#                                                         #
#     Python code based on version 1.07 of DimeMC         # 
#                                                         #
#     For python problems: leonardo.negri@helsinki.fi     #
#                                                         #
#     Comments etc to: l.harland-lang@ucl.ac.uk           #
#                                                         #
#     Please do not bother Dr. Harland-Lang if you find   #
#     issues with the python code, he wrote the fortran   #
#     code. Email leonardo.negri@helsinki.fi instead.     #
#                                                         #
#     If you use the code please reference (and see for   #
#     details of model used):                             #
#                                                         #
#     L.A. Harland-Lang, V.A. Khoze, M.G. Ryskin          #
#     and W.J. Stirling arXiv:1105.1626                   #
#                                                         #
#     L.A. Harland-Lang, V.A. Khoze, M.G. Ryskin          #
#     and W.J. Stirling arXiv:1204.4803                   #
#                                                         #
#     L.A. Harland-Lang, V.A. Khoze and  M.G. Ryskin      #
#     arXiv:1312.4553                                     #
#                                                         # 
###########################################################
''')