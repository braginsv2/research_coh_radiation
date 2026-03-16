#!/usr/bin/env python3
"""
ANB - Coherent Bremsstrahlung Calculator
Port of anb.c to Python

compile : cc anb.c vec.c -g -lm -o anb
input   : 'name'.in      for parameters check comments in example
output  : 'name'.txt     readable output (ascii file)
  contents of each row:
    photon energy, total coherent Intensity, diff coherent I = perp-para,
    incoherent crystal I, incoherent amorphous I,
    electron contrib crystal I, electron amorphous I
  use to calculate:
    Irelexp = (Isum+Iinc)/Iamo  for experimental comparison
    Irel    = (Isum+Iinc)/Iinc
    Pol     = Idif/(Isum+Iinc)
"""

import sys
import math
import time
import os


# -------------------------------------------------------
# DEFINES
# -------------------------------------------------------
TRUE  = 1
FALSE = 0

BEAMX_STEPS   = 20
BEAMY_STEPS   = 20
BEAMXY_REGION = 2
MAX_LATTICEVEC = 1000
MASS_E        = 0.511003373282   # electron mass in MeV
A             = 923.7            # Diamond spacing, lattice constants
M_2PI         = 2.0 * math.pi
GUNIT         = M_2PI / A
DEBYE_ON      = TRUE
DEBYE_OFF     = FALSE
SCREEN_DEBOF  = 111.0
SCREEN_DEBON  = 135.6
TNULL         = 273.15
M_EulerC      = 0.5772156649
X_RANGE       = 1000
USTEPS        = 100
UMAX          = 5
RSTEPS        = 5000
THETARANGE    = 16.73
DQ_STEPS      = 1000
ESTEPS        = 9
A_12C         = 12.0
ALPHA2        = 1.0 / (137.0 * 137.0)
DENS_DIAMOND  = 3.513

LOGTERM  = 1
LOGFILE  = 2
LOGERR   = 4
LOGBOTH  = LOGTERM | LOGFILE
LOGERROR = LOGERR  | LOGFILE

V_ZERO = 1e-8

# -------------------------------------------------------
# VEC operations (3D double vector)
# -------------------------------------------------------

class VEC:
    __slots__ = ('x', 'y', 'z')
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z
    def __repr__(self):
        return f"VEC({self.x}, {self.y}, {self.z})"

class MATRIX:
    __slots__ = ('x', 'y', 'z')
    def __init__(self, x=None, y=None, z=None):
        self.x = x if x else VEC()
        self.y = y if y else VEC()
        self.z = z if z else VEC()

v_null = VEC(0.0, 0.0, 0.0)

def v_set(x, y, z):
    return VEC(x, y, z)

def v_setpol(r, th, ph):
    return v_smul(r, v_set(math.sin(th)*math.cos(ph),
                           math.sin(th)*math.sin(ph),
                           math.cos(th)))

def v_add(a, b):
    return VEC(a.x + b.x, a.y + b.y, a.z + b.z)

def v_sub(a, b):
    return VEC(a.x - b.x, a.y - b.y, a.z - b.z)

def v_minus(a):
    return VEC(-a.x, -a.y, -a.z)

def v_smul(s, a):
    return VEC(a.x * s, a.y * s, a.z * s)

def v_sprd(a, b):
    return a.x*b.x + a.y*b.y + a.z*b.z

def v_sqr(a):
    return a.x*a.x + a.y*a.y + a.z*a.z

def v_kreuz(a, b):
    return VEC(a.y*b.z - a.z*b.y,
               a.z*b.x - a.x*b.z,
               a.x*b.y - a.y*b.x)

def v_norm(a):
    return math.sqrt(v_sqr(a))

def v_vdet(a1, a2, a3):
    return (a1.x*a2.y*a3.z + a1.y*a2.z*a3.x + a1.z*a2.x*a3.y
          - a1.z*a2.y*a3.x - a1.x*a2.z*a3.y - a1.y*a2.x*a3.z)

def v_det(A):
    return (A.x.x*A.y.y*A.z.z + A.x.y*A.y.z*A.z.x + A.x.z*A.y.x*A.z.y
          - A.x.z*A.y.y*A.z.x - A.x.x*A.y.z*A.z.y - A.x.y*A.y.x*A.z.z)

def v_scale(scale, a):
    return VEC(a.x * scale.x, a.y * scale.y, a.z * scale.z)

def v_scadd(scale, a, pos):
    return VEC(a.x * scale.x + pos.x,
               a.y * scale.y + pos.y,
               a.z * scale.z + pos.z)

def v_mmul(A, x):
    return VEC(A.x.x*x.x + A.y.x*x.y + A.z.x*x.z,
               A.x.y*x.x + A.y.y*x.y + A.z.y*x.z,
               A.x.z*x.x + A.y.z*x.y + A.z.z*x.z)

def v_lgs(A, b):
    c = VEC(v_vdet(b, A.y, A.z),
            v_vdet(A.x, b, A.z),
            v_vdet(A.x, A.y, b))
    return v_smul(1.0 / v_det(A), c)

def v_equ(a, b):
    c = v_sub(a, b)
    return v_sprd(c, c) < 1.0e-5

def v_normal(a):
    return v_smul(1.0 / v_norm(a), a)

def v_parallel(a, v):
    v = v_normal(v)
    return v_smul(v_sprd(a, v), v)

def v_perpend(a, v):
    return v_sub(a, v_parallel(a, v))

def v_transverse(a):
    return math.sqrt(a.x*a.x + a.y*a.y)

def v_cart2pol(a):
    rt2 = a.x*a.x + a.y*a.y
    r = math.sqrt(rt2 + a.z*a.z)
    if r < V_ZERO:
        return 0.0, 0.0, 0.0
    t = math.acos(a.z / r)
    if rt2 == 0:
        return r, t, 0.0
    if a.y > 0.0:
        p = math.acos(a.x / math.sqrt(rt2))
    else:
        p = 2.0 * math.pi - math.acos(a.x / math.sqrt(rt2))
    return r, t, p

def v_polrot(th, ph, v):
    if th < 0.0:
        raise ValueError(f"vec.v_polrot > ERROR: polar angle ({th}) negative!!!")
    while ph < 0.0:
        ph += M_2PI
    cth = math.cos(th); sth = math.sin(th)
    cph = math.cos(ph); sph = math.sin(ph)
    return VEC(
        cth*cph*v.x - cth*sph*v.y - sth*v.z,
        sph    *v.x + cph    *v.y,
        sth*cph*v.x - sth*sph*v.y + cth*v.z
    )

def v_polrotinv(th, ph, v):
    if th < 0.0:
        raise ValueError(f"vec.v_polrotinv > ERROR: polar angle ({th}) negative!!!")
    while ph < 0.0:
        ph += M_2PI
    th = M_2PI - th
    ph = M_2PI - ph
    cth = math.cos(th); sth = math.sin(th)
    cph = math.cos(ph); sph = math.sin(ph)
    return VEC(
        cph*cth*v.x - sph*v.y - cph*sth*v.z,
        sph*cth*v.x + cph*v.y - sph*sth*v.z,
        sth    *v.x           + cth    *v.z
    )

def v_zrotate(z, v):
    r, th, ph = v_cart2pol(z)
    return v_polrot(th, ph, v)

# -------------------------------------------------------
# Global state (replaces C global variables)
# -------------------------------------------------------

class State:
    def __init__(self):
        # lattice vectors
        self.lv      = [VEC() for _ in range(MAX_LATTICEVEC)]
        self.S_2     = [0.0] * MAX_LATTICEVEC
        self.g_2     = [0.0] * MAX_LATTICEVEC
        self.aps     = [0.0] * MAX_LATTICEVEC
        self.debye   = [0.0] * MAX_LATTICEVEC
        self.ff      = [0.0] * MAX_LATTICEVEC
        self.Irel    = [0.0] * MAX_LATTICEVEC

        self.sig_colx2 = 0.0
        self.sig_coly2 = 0.0
        self.sig_colr2 = 0.0
        self.sig_colas = 0.0
        self.sig_colli = 0.0
        self.collifct  = [0.0] * (USTEPS + 1)

        self.sigx  = 0.0; self.sigy  = 0.0
        self.sigpx = 0.0; self.sigpy = 0.0
        self.sig_mol = 0.0
        self.sig_bx  = 0.0; self.sig_by = 0.0
        self.del_bx  = 0.0; self.del_by = 0.0

        self.umax        = 0.0
        self.alpha       = 0.0
        self.theta       = 0.0
        self.phi         = 0.0
        self.di_thickness= 0.0

        self.E0MeV    = 0.0
        self.sigE0MeV = 0.0
        self.E0       = 0.0
        self.dist     = 0.0
        self.collength= 0.0
        self.rkoll    = 0.0
        self.C        = 0.0

        self.temperature = 0.0
        self.Zamo  = 0.0
        self.Zcry  = 0.0
        self.A3    = 0.0
        self.Atot  = 0.0
        self.Adebye= 0.0

        self.type           = 0
        self.max_latticevec = 0
        self.NMax           = 0

        self.file_log   = None
        self.input_file = ""
        self.logname    = ""
        self.outname    = ""
        self.comment    = ""
        self.Deb_cmnt   = ""
        self._logbuf    = ""
        self.date       = ""

        # function pointer equivalent
        self.int_inc = None

st = State()

# -------------------------------------------------------
# logging
# -------------------------------------------------------

def logging_(what):
    if LOGTERM == (LOGTERM & what):
        print(st._logbuf, end='')
        sys.stdout.flush()
    if LOGFILE == (LOGFILE & what):
        if st.file_log:
            st.file_log.write(st._logbuf)
    if LOGERR == (LOGERR & what):
        sys.stderr.write(st._logbuf)

# -------------------------------------------------------
# get_par: read one parameter from file
# -------------------------------------------------------

def get_par(file, typ):
    """Read one parameter from file, skip comment/blank lines."""
    while True:
        line = file.readline()
        if line == '':
            st._logbuf = "get_par > FATAL ERROR: couldn't get parameter\n"
            logging_(LOGERROR)
            sys.exit(-1)
        stripped = line.strip()
        if stripped and stripped[0] not in ('#', '\n', ' ', '\t') or \
           (stripped and stripped[0] not in ('#',)):
            # use any non-comment, non-empty line
            if stripped == '' or stripped[0] == '#':
                continue
            break
    if typ == 's':
        return line  # whole line including newline
    elif typ == 'c':
        return stripped[0]
    elif typ == 'f':
        return float(stripped.split()[0])
    elif typ == 'd':
        return float(stripped.split()[0])
    elif typ == 'i':
        return int(stripped.split()[0])
    return stripped

# -------------------------------------------------------
# input_parameter
# -------------------------------------------------------

def input_parameter():
    try:
        con = open(st.input_file, 'r')
    except IOError:
        st._logbuf = f"Fatal Error > unable to open input file: '{st.input_file}'.\n"
        logging_(LOGERROR)
        sys.exit(-1)

    st._logbuf = (f"\n\n********* ANB Version 2.0 ({st.date}) ********* Alexander Natter/11/98\n\n"
                  f"input > Reading input file {st.input_file}\n")
    logging_(LOGBOTH)

    st.comment = get_par(con, 's')
    st._logbuf = f"\ninput > Comment:\n{st.comment}\n"
    logging_(LOGBOTH)

    st.E0MeV = get_par(con, 'd')
    st.file_log.write(f"energy = {st.E0MeV} MeV\n")

    st.sigE0MeV = get_par(con, 'd')
    st.file_log.write(f"energy spread = {st.sigE0MeV} MeV\n")

    st.theta = get_par(con, 'd')
    st.file_log.write(f"Theta = {st.theta} rad\n")

    st.alpha = get_par(con, 'd')
    st.file_log.write(f"Alpha = {st.alpha} rad\n")

    st.phi = get_par(con, 'd')
    st.file_log.write(f"phi = {st.phi} rad\n")

    st.sigx = get_par(con, 'd')
    st.file_log.write(f"sigma_x = {st.sigx} mm\n")

    st.sigy = get_par(con, 'd')
    st.file_log.write(f"sigma_y = {st.sigy} mm\n")

    st.sigpx = get_par(con, 'd')
    st.file_log.write(f"sigma_px = {st.sigpx} mrad\n")

    st.sigpy = get_par(con, 'd')
    st.file_log.write(f"sigma_py = {st.sigpy} mrad\n")

    st.di_thickness = get_par(con, 'd')
    st.file_log.write(f"Di-thickness = {st.di_thickness} mm\n")

    st.dist = get_par(con, 'd')
    st.file_log.write(f"Distance target - collimator = {st.dist} m\n")

    st.collength = get_par(con, 'd')
    st.file_log.write(f"Length of collimator = {st.collength} m\n")

    st.rkoll = get_par(con, 'd')
    st.file_log.write(f"collimator radius = {st.rkoll} mm,\n")

    st.type = get_par(con, 'i')
    st.file_log.write(f"{'in' if st.type < 9 else ''}coherente Bremsstrahlung (type:{st.type})\n")

    st.max_latticevec = get_par(con, 'i')
    st.file_log.write(
        f"Max Lattice vector index in Lattice vector table i={st.max_latticevec} in [0..{MAX_LATTICEVEC}]\n")

    st.Zcry = get_par(con, 'd')
    st.file_log.write(f"using Z={st.Zcry:3.0f} for 'crystal' radiator\n")

    st.Zamo = get_par(con, 'd')
    st.file_log.write(f"using Z={st.Zamo:3.0f} for 'amorphous' radiator\n")

    con.close()
    st.file_log.flush()

# -------------------------------------------------------
# calc_collifct
# -------------------------------------------------------

def collifctr(u, r, s):
    return r * math.exp(r*r / (-2.0*s*s)) / (s*s) * math.acos(
        (r*r + u*u - st.C*st.C) / (2.0*r*u + 1e-10)) / math.pi

def calc_collifct(sig):
    st._logbuf = "calc_collifct > calculating collimation function\n"
    logging_(LOGBOTH)

    umax   = st.umax
    ustep  = umax / USTEPS
    C      = st.C

    if sig == 0.0:
        u = 0.0
        for i in range(USTEPS + 1):
            st.collifct[i] = 1.0 if u < C else 0.0
            u += ustep
    else:
        u = 0.0
        for i in range(USTEPS + 1):
            sum_ = 0.0
            if u > 0.0:
                rstep = (C + u - abs(C - u)) / RSTEPS
                r1    = abs(C - u) + rstep
                r2    = C + u - rstep

                sum_ = collifctr(u, r1, sig) * rstep / 2.0
                r = r1
                while r <= r2:
                    sum_ += collifctr(u, r, sig) * rstep
                    r += rstep
                sum_ += collifctr(u, r2, sig) * rstep / 2.0

            if C > u:
                sum_ += 1.0 - math.exp(-(C - u)*(C - u) / (2.0*sig*sig))

            st.collifct[i] = sum_
            u += ustep

    st._logbuf = "calc_collifct > done\n"
    logging_(LOGFILE)

def collifunction(x, xd):
    umax  = st.umax
    ustep = umax / USTEPS
    C     = st.C

    if x == 0.0 or x > xd:
        return 0.0

    denom = 1.0 - xd
    if denom == 0.0:
        return 0.0
    u = math.sqrt((xd / x - 1.0) / denom)
    if u >= umax:
        return 0.0

    i  = int(u / ustep)
    if i >= USTEPS:
        i = USTEPS - 1
    f1 = st.collifct[i]
    f2 = st.collifct[i + 1]

    return f1 + (f2 - f1) / ustep * (u - i * ustep)

# -------------------------------------------------------
# int_Hub
# -------------------------------------------------------

def int_Hub(Z, debye_switch):
    iinc = [0.0] * (X_RANGE + 1)
    iiel = [0.0] * (X_RANGE + 1)

    C  = st.C
    E0 = st.E0

    fc = C*C / (1.0 + C*C)
    Z3 = Z ** (1.0/3.0)
    v  = 1.0 / (1.0 + C*C)
    v2 = v * v

    screening = SCREEN_DEBON if debye_switch == DEBYE_ON else SCREEN_DEBOF

    for i in range(1, X_RANGE):
        x   = float(i) / X_RANGE
        emx = 1.0 - x
        b   = 2.0 * Z3 * E0 / screening * emx / x
        at  = math.atan(C*C*b / (1.0 + C*C + b*b)) * 2.0 / b
        dmy = (x / (2.0 * E0 * emx)) ** 2.0
        Mv  = -math.log(dmy + (v   * Z3 / screening) ** 2.0)
        M1  = -math.log(dmy + (      Z3 / screening) ** 2.0)
        psi1 = 2.0 * (1.0 + M1 - (1.0 + Mv)*v - at)
        psi2 = (-40.0/3.0*v*v2 + 18.0*v2 - (8.0/(b*b) + 6.0)*v + 8.0/(b*b) + 4.0/3.0
                + (4.0*v*v2 - 6.0*v2)*Mv + 2.0*M1
                - 6.0/(b*b)*(Mv - M1 + 2.0/3.0*at))
        iinc[i] = (1.0 + emx*emx) * psi1 - 2.0/3.0 * emx * psi2

        eps = 100.0 / (E0 * Z3*Z3) * x / emx
        if eps <= 0.88:
            eps    = 0.88 - eps
            psieps = 19.70 + eps*(4.177 + eps*(-3.806 + eps*(31.84 + eps*(-58.63 + eps*40.77))))
        else:
            psieps = 19.9 - 4.0 * math.log(eps)

        psi1_e  = psieps - 4.0 - (8//3) * math.log(Z)   # C integer division: 8/3 = 2
        psi2_e  = psi1_e + 0.6666666
        iiel[i] = ((1.0 + emx*emx)*psi1_e - 2.0/3.0*emx*psi2_e) * fc / Z

    iinc[X_RANGE] = 0.0
    iiel[X_RANGE] = 0.0
    return iinc, iiel

# -------------------------------------------------------
# int_HubBd
# -------------------------------------------------------

def int_HubBd(Z, debye_switch):
    iinc = [0.0] * (X_RANGE + 1)
    iiel = [0.0] * (X_RANGE + 1)

    sig_colli = st.sig_colli
    U0   = st.C
    norm = 0.0
    drho = sig_colli / 11.0
    sg2  = 1.0 / (2.0 * sig_colli * sig_colli)
    dph  = M_2PI / 31.0

    st._logbuf = (f"int_HubBd > folding electron divergence and hubbel intensity "
                  f"{'with' if debye_switch else 'without'} debye-factor\n")
    logging_(LOGBOTH)

    rho = drho / 2.0
    while rho < U0:
        ph = dph / 2.0
        while ph < M_2PI:
            st.C = math.sqrt(U0*U0 - (rho*math.sin(ph))**2.0) - rho*math.cos(ph)
            weight = (dph*drho / math.pi
                      * math.exp(-rho*rho*sg2)
                      * rho / (M_2PI * sig_colli * sig_colli))
            norm += weight
            ii, ie = int_Hub(Z, debye_switch)
            for i in range(1, X_RANGE):
                iinc[i] += ii[i] * weight
                iiel[i] += ie[i] * weight
            ph += dph
        rho += drho

    st.C  = U0
    drho2 = sig_colli / 6.0
    rho   = U0 + drho2 / 2.0
    while rho < 3.0 * U0:
        phc = math.pi - math.asin(U0 / rho)
        dph2 = (M_2PI - 2.0*phc) / 10.0
        ph = phc + dph2 / 2.0
        while ph <= M_2PI - phc:
            sq = math.sqrt(U0*U0 - (rho*math.sin(ph))**2.0)
            st.C = sq - rho*math.cos(ph)
            Cm   = -sq - rho*math.cos(ph)
            weight = (dph2*drho2 / math.pi
                      * math.exp(-rho*rho*sg2)
                      * rho / (M_2PI * sig_colli * sig_colli))
            norm += weight
            ii, ie = int_Hub(Z, debye_switch)
            for i in range(1, X_RANGE):
                iinc[i] += ii[i] * weight
                iiel[i] += ie[i] * weight
            st.C = Cm
            ii, ie = int_Hub(Z, debye_switch)
            for i in range(1, X_RANGE):
                iinc[i] -= ii[i] * weight
                iiel[i] -= ie[i] * weight
            ph += dph2
        rho += drho2

    for i in range(1, X_RANGE):
        iinc[i] /= norm
        iiel[i] /= norm

    st.C = U0
    st._logbuf = "int_HubBd > done\n"
    logging_(LOGBOTH)
    return iinc, iiel

# -------------------------------------------------------
# int_BH
# -------------------------------------------------------

def int_BH(Z, debye_switch):
    iinc = [0.0] * (X_RANGE + 1)
    iiel = [0.0] * (X_RANGE + 1)

    psi1e = 4.05
    psi2e = 3.94

    C      = st.C
    E0     = st.E0
    Adebye = st.Adebye
    fc     = C*C / (1.0 + C*C)

    st._logbuf = (f"int_BH > starting to calculate the incoherent (Bethe Heitler) intensities\n"
                  f"int_inc > Z of radiator: {Z}\t{'using' if debye_switch else 'not using'} Debye Factor"
                  f"\tColli: {C}\n")
    logging_(LOGBOTH)

    for i in range(1, X_RANGE):
        x   = float(i) / 1000.0
        emx = 1.0 - x
        del_ = x / (2.0 * E0 * emx)
        dq   = (1.0 - del_) / float(DQ_STEPS)
        psi1 = 0.0
        psi2 = 0.0
        q = del_
        while q < 1.0:
            q2 = q*q
            q3 = q*q2
            if debye_switch == DEBYE_ON:
                debfact = 1.0 - math.exp(-Adebye * q2)
            else:
                debfact = 1.0
            ff = 1.0 - (0.2283 + 1.8359*math.exp(-10528*q2)
                        + 1.8119*math.exp(-4678*q2)
                        + 1.5809*math.exp(-239*q2)
                        + 0.5426*math.exp(-27116*q2)) / 6.0
            dmy   = debfact * ff*ff / q3
            psi1 += dmy * (q - del_)**2.0
            psi2 += dmy * (q2 + del_*del_*(3.0 - 6.0*math.log(q/del_) - 4.0*del_/q))
            q += dq

        psi1 = 4.0     + 4.0 * psi1 * dq
        psi2 = 10.0/3.0 + 4.0 * psi2 * dq
        iinc[i] = ((1.0 + emx*emx)*psi1 - 2.0/3.0*emx*psi2) * fc
        iiel[i] = ((1.0 + emx*emx)*psi1e - 2.0/3.0*emx*psi2e) * fc

    iinc[X_RANGE] = 0.0
    iiel[X_RANGE] = 0.0
    st._logbuf = "int_BH > done\n"
    logging_(LOGFILE)
    return iinc, iiel

# -------------------------------------------------------
# cohbrems
# -------------------------------------------------------

def cohbrems(alpha, theta):
    isum = [0.0] * (X_RANGE + 1)
    idif = [0.0] * (X_RANGE + 1)

    E0 = st.E0

    for l in range(st.max_latticevec):
        g    = st.lv[l]
        dmy  = math.sin(theta) * (g.y * math.cos(alpha) + g.z * math.sin(alpha))
        gl   = (math.cos(theta)*g.x + dmy) * GUNIT
        gt_2 = (g.y*g.y + g.z*g.z + (math.sin(theta)*g.x)**2.0 - dmy*dmy) * GUNIT*GUNIT
        xd   = 2.0 * E0 * gl / (1.0 + 2.0 * E0 * gl)
        xc   = xd / (1.0 + st.C*st.C*(1.0 - xd))

        if st.g_2[l] == 0.0 or gl == 0.0:
            continue

        dmy2 = M_2PI * st.ff[l] / (st.g_2[l] * gl*gl)
        Gg   = st.A3 / 8.0 * dmy2*dmy2 * st.S_2[l] * st.debye[l]

        for i in range(1, X_RANGE):
            x = float(i) / 1000.0
            if x > xd:
                break
            del_ = x / (2.0 * E0 * (1.0 - x))
            psi1 =  4.0 * Gg * del_       * gt_2 * gl*gl
            psi2 = 24.0 * Gg * del_*del_  * gt_2 * (gl - del_)
            psi3 = -4.0 * Gg * del_**3.0  * st.aps[l]
            cf   = collifunction(x, xd)
            isum[i] += ((1.0 + (1.0-x)**2.0)*psi1 - 2.0/3.0*(1.0-x)*psi2) * cf
            idif[i] -= 2.0*(1.0-x)*psi3 * cf

    return isum, idif

# -------------------------------------------------------
# int_coh
# -------------------------------------------------------

def int_coh():
    intsum = [0.0] * (X_RANGE + 1)
    intdif = [0.0] * (X_RANGE + 1)

    st._logbuf = "int_coh > calculating the coherent (sum, dif) intensities\n"
    logging_(LOGBOTH)

    # precompute per-lattice-vector factors
    for l in range(st.max_latticevec):
        g = st.lv[l]
        st.g_2[l]   = v_sprd(g, g) * GUNIT*GUNIT
        st.aps[l]   = ((g.y*g.y - g.z*g.z)*math.cos(2.0*st.phi)
                       + 2.0*g.y*g.z*math.sin(2.0*st.phi)) * GUNIT*GUNIT
        st.debye[l] = math.exp(-st.Adebye * st.g_2[l])
        st.ff[l]    = (0.2283 + 1.8359*math.exp(-10528*st.g_2[l])
                       + 1.8119*math.exp(-4678*st.g_2[l])
                       + 1.5809*math.exp(-239*st.g_2[l])
                       + 0.5426*math.exp(-27116*st.g_2[l]))
        st.ff[l]    = 1.0 - st.ff[l] / 6.0

    norm = 0.0
    beam_x = -BEAMXY_REGION * st.sig_bx
    while beam_x <= BEAMXY_REGION * st.sig_bx:
        beam_y = -BEAMXY_REGION * st.sig_by
        while beam_y <= BEAMXY_REGION * st.sig_by:
            th_beam  = math.sqrt(beam_x*beam_x + beam_y*beam_y)
            phi_beam = (math.atan(beam_y / beam_x) if beam_x != 0.0 else 0.0)
            discr    = (st.theta*st.theta + th_beam*th_beam
                        + 2.0*st.theta*th_beam*math.cos(phi_beam))
            if discr <= 0.0:
                st._logbuf = (f"inc_coh > FATAL ERROR: beam divergence too large: "
                              f"phi:{phi_beam}\tth:{th_beam}\tdiscr:{discr}\n")
                logging_(LOGERROR)
                beam_y += st.del_by
                continue
            th = math.sqrt(discr)
            al = st.alpha + math.asin(math.sin(phi_beam)*th_beam/th) if th > 0 else st.alpha
            weight = ((math.exp(-(beam_x*beam_x)/(2.0*st.sig_bx*st.sig_bx)) if st.sig_bx > 0.0 else 1.0) *
                      (math.exp(-(beam_y*beam_y)/(2.0*st.sig_by*st.sig_by)) if st.sig_by > 0.0 else 1.0))
            norm += weight
            st._logbuf = (f"int_coh > calculating intensities for alpha:{al}\ttheta:{th}"
                          f"\tphib:{phi_beam}\tthb:{th_beam}\tweight:{weight}\n")
            logging_(LOGFILE)

            isum, idif = cohbrems(al, th)
            for i in range(1, X_RANGE + 1):
                intsum[i] += isum[i] * weight
                intdif[i] += idif[i] * weight

            beam_y += st.del_by
        beam_x += st.del_bx

    st._logbuf = f"int_coh > sum of theta/phi-crystal weight: {norm}\n"
    logging_(LOGFILE)

    for i in range(1, X_RANGE + 1):
        intsum[i] /= norm
        intdif[i] /= norm

    st._logbuf = "int_coh > done\n"
    logging_(LOGFILE)
    return intsum, intdif

# -------------------------------------------------------
# sort lattice vectors
# -------------------------------------------------------

def sort_lv(n):
    """Straight insertion sort (piksrt) by Irel descending."""
    for j in range(1, n):
        Ii = st.Irel[j]
        Si = st.S_2[j]
        gi = VEC(st.lv[j].x, st.lv[j].y, st.lv[j].z)
        i  = j - 1
        while i >= 0 and st.Irel[i] < Ii:
            st.Irel[i+1] = st.Irel[i]
            st.S_2[i+1]  = st.S_2[i]
            st.lv[i+1]   = st.lv[i]
            i -= 1
        st.Irel[i+1] = Ii
        st.S_2[i+1]  = Si
        st.lv[i+1]   = gi

# -------------------------------------------------------
# calc_latticevec
# -------------------------------------------------------

def calc_latticevec(lvmax_in):
    E0     = st.E0
    Adebye = st.Adebye
    theta  = st.theta
    alpha  = st.alpha

    hklmax = float(int((lvmax_in**(1.0/3.0) - 1.0) / 2.0))
    lvmax  = int((2.0*hklmax + 1.0)**3.0)

    st._logbuf = (f"calc_latticevec > calculating relative contributions of {lvmax} lattice vectors\n")
    logging_(LOGBOTH)
    st._logbuf = (f"calc_latticevec > checking {lvmax} lattice vectors, maximal miller index:{hklmax}\n")
    logging_(LOGFILE)

    if 3 * int(max(hklmax - 0.5, 0)**2.5) > MAX_LATTICEVEC:
        st._logbuf = f"calc_latticevec > Error: Too many lattice vects requested\nMaximal {MAX_LATTICEVEC} possible\n"
        logging_(LOGERROR)
        sys.exit(-1)

    lvidx = 0
    gx = -hklmax
    while gx <= hklmax:
        gy = -hklmax
        while gy <= hklmax:
            gz = -hklmax
            while gz <= hklmax:
                mod = abs(int(gx))%2 + abs(int(gy))%2 + abs(int(gz))%2
                if mod == 0 and int(gx+gy+gz) % 4 == 0:
                    S2 = 64
                elif mod == 3:
                    S2 = 32
                else:
                    S2 = 0

                if S2 > 0:
                    dmy_ = math.sin(theta) * (gy*math.cos(alpha) + gz*math.sin(alpha))
                    gl   = (math.cos(theta)*gx + dmy_) * GUNIT
                    if gl > 0:
                        xd = 2.0 * E0 * gl / (1.0 + 2.0 * E0 * gl)
                        if 0.0 < xd < 1.0:
                            g     = VEC(gx, gy, gz)
                            g_2   = v_sprd(g, g) * GUNIT*GUNIT
                            gt_2  = (gy*gy + gz*gz + (math.sin(theta)*gx)**2.0 - dmy_*dmy_) * GUNIT*GUNIT
                            ff    = (0.2283 + 1.8359*math.exp(-10528*g_2)
                                     + 1.8119*math.exp(-4678*g_2)
                                     + 1.5809*math.exp(-239*g_2)
                                     + 0.5426*math.exp(-27116*g_2))
                            ff    = 1.0 - ff/6.0
                            del_  = xd / (2.0*E0*(1.0-xd))
                            psi1  = st.A3 * M_2PI*M_2PI / 2.0 * ff*ff * S2 * math.exp(-Adebye*g_2) / del_ / g_2
                            I_val = (1.0 + (1.0-xd)**2.0) * psi1
                            st.Irel[lvidx] = I_val
                            st.S_2[lvidx]  = S2
                            st.lv[lvidx]   = g
                            lvidx += 1
                            if lvidx > MAX_LATTICEVEC:
                                st._logbuf = f"\n\ncalc_latticevec > Error: MAX_LATTICEVEC exceeded (lvidx:{lvidx})\n\n"
                                logging_(LOGERROR)
                                sys.exit(-1)
                gz += 1.0
            gy += 1.0
        gx += 1.0

    sort_lv(lvidx)
    st._logbuf = f"calc_latticevec > sorting {lvidx} contributing lattice vectors\n"
    logging_(LOGFILE)

    for i in range(lvidx):
        gl = (math.cos(theta)*st.lv[i].x
              + math.sin(theta)*(st.lv[i].y*math.cos(alpha) + st.lv[i].z*math.sin(alpha))) * GUNIT
        xd = 2.0*E0*gl / (1.0 + 2.0*E0*gl)
        st.file_log.write(f"{st.Irel[i]}\t{st.S_2[i]}\t{xd}\tg: {st.lv[i].x:10.5f} {st.lv[i].y:10.5f} {st.lv[i].z:10.5f}\n")

    st._logbuf = (f"calc_latticevec > using the {st.max_latticevec} most contributing lattice vectors "
                  f"out of {lvidx} for further calculation\n")
    logging_(LOGFILE)

# -------------------------------------------------------
# moliere_scattering
# -------------------------------------------------------

def solveeq(b):
    re = 40.0; li = 1.2
    f_r = re - math.log(re) - b
    f_l = li - math.log(li) - b
    next_ = (re + li) / 2.0
    while re - li > 1e-9:
        next_   = (re + li) / 2.0
        f_next  = next_ - math.log(next_) - b
        if f_next <= 0.0:
            li = next_; f_l = f_next
        else:
            re = next_; f_r = f_next
    return next_

def sig_moliere(s, pe):
    s_ = s / 10.0
    om = (7800.0 * (1.0 + st.Zcry) * st.Zcry**(1.0/3.0) * DENS_DIAMOND * s_
          / (A_12C * (1.0 + 3.35*(st.Zcry*st.Zcry*ALPHA2))))
    b  = math.log(om) + 1.0 - 2.0*M_EulerC
    B  = solveeq(b)
    thsq = 0.157 * st.Zcry*(1.0+st.Zcry) * DENS_DIAMOND*s_ / (A_12C*pe*pe)
    return math.sqrt(thsq * (B - 1.2))

def moliere_scattering():
    MSSTEPS = 100
    step = st.di_thickness / MSSTEPS
    st._logbuf = (f"moliere > starting to calculate moliere variance, "
                  f"sigmoliere(0,thickness): {sig_moliere(st.di_thickness, st.E0MeV)*1e3} (mrad)\n")
    logging_(LOGFILE)
    sm = 0.0
    s  = step
    while s < st.di_thickness:
        sm += sig_moliere(s, st.E0MeV)
        s  += step
    st._logbuf = f"moliere > mean sigma_moliere: {sm/MSSTEPS*1e3} mrad\n"
    logging_(LOGFILE)
    return sm / MSSTEPS

# -------------------------------------------------------
# energy_spread
# -------------------------------------------------------

def energy_spread(E0, sigE0MeV_val, I):
    Ie    = [0.0] * (X_RANGE + 1)
    weight = [0.0] * (ESTEPS + 2)

    st._logbuf = "energy_spread > fold intensity with primary beam energy spread \n"
    logging_(LOGBOTH)

    Emin   = E0 - 2.0 * sigE0MeV_val / MASS_E
    Emax   = E0 + 2.0 * sigE0MeV_val / MASS_E
    Estep  = (Emax - Emin) / ESTEPS
    sig2E0 = 2.0 * sigE0MeV_val*sigE0MeV_val / (MASS_E*MASS_E)
    st._logbuf = (f"energy_spread > E0:{E0}\tsigE0MeV:{sigE0MeV_val}\tEstep:{Estep}"
                  f"\tEmin:{Emin}\tEmax:{Emax}\n")
    logging_(LOGFILE)

    sumweight = 0.0
    for Eidx in range(ESTEPS + 1):
        Es = Emin + Eidx * Estep
        weight[Eidx] = math.exp(-(Es - E0)**2.0 / sig2E0)
        sumweight += weight[Eidx]

    for i in range(1, X_RANGE):
        x = float(i) / X_RANGE
        for Eidx in range(ESTEPS + 1):
            Es = Emin + Eidx * Estep
            x1 = x * E0 / Es
            i0 = int(X_RANGE * x1)
            if i0 >= X_RANGE:
                continue
            x0 = float(i0) / X_RANGE
            Ie[i] += (I[i0] + (I[i0+1] - I[i0]) * (x1 - x0)) * weight[Eidx]

    for i in range(1, X_RANGE):
        I[i] = Ie[i] / sumweight

    return I

# -------------------------------------------------------
# calc_debyered (stub -- not called in main flow of anb.c)
# -------------------------------------------------------

def calc_debyered():
    pass

# -------------------------------------------------------
# main
# -------------------------------------------------------

def main():
    # init arrays
    intsum = [0.0] * (X_RANGE + 1)
    intdif = [0.0] * (X_RANGE + 1)
    intinc = [0.0] * (X_RANGE + 1)
    intamo = [0.0] * (X_RANGE + 1)
    intiel = [0.0] * (X_RANGE + 1)
    intael = [0.0] * (X_RANGE + 1)

    # get time
    st.date = time.ctime().rstrip('\n')

    for i in range(1, X_RANGE + 1):
        intsum[i] = 0.0
        intdif[i] = 0.0

    # set up input/output filenames
    if len(sys.argv) == 2:
        st.input_file = sys.argv[1]
    else:
        st.input_file = input("Please enter the name of the input file (without the ending 'in')?\n").strip()

    # strip possible .in suffix
    if st.input_file.endswith('.in'):
        st.input_file = st.input_file[:-3]

    st.outname    = st.input_file + ".txt"
    st.logname    = st.input_file + ".log"
    st.input_file = st.input_file + ".in"

    try:
        st.file_log = open(st.logname, 'w')
    except IOError:
        sys.stderr.write("\n\nFatal Error: couldn't open log file\n\n")
        sys.exit(-1)

    input_parameter()

    st._logbuf = (f"anb > Input file read, please wait approx 5 minutes...\n"
                  f"anb > Output is written in file {st.outname} and log file is {st.logname}\n")
    logging_(LOGBOTH)

    st.E0 = st.E0MeV / MASS_E

    if st.max_latticevec > MAX_LATTICEVEC:
        st.max_latticevec = MAX_LATTICEVEC
        st._logbuf = f"anb > Warning: highest value for 'MAX_LATTICEVEC' = {MAX_LATTICEVEC} using\n"
        logging_(LOGFILE)

    if st.di_thickness > 0.0:
        st.sig_mol = moliere_scattering()

    st.sig_bx = math.sqrt(st.sig_mol**2.0 + st.sigpx**2.0 * 1e-6)
    st.sig_by = math.sqrt(st.sig_mol**2.0 + st.sigpy**2.0 * 1e-6)
    st._logbuf = f"anb > Mean electron divergence: sig_bx/y: {st.sig_bx*1e3}\t{st.sig_by*1e3} mrad\n"
    logging_(LOGFILE)

    if st.dist == 0.0:
        st.C    = st.E0 * M_2PI
        st.umax = st.C
        if st.type == 2:
            st.type = 1
    else:
        st.C        = st.rkoll * 1e-3 / (st.dist + st.collength) * st.E0
        st.umax     = UMAX * st.C
        st.sig_colx2 = ((st.sigx**2.0 * 1e-6 / (st.dist**2.0)
                         + st.sigpx**2.0 * 1e-6
                         + st.sig_mol**2.0) * st.E0*st.E0)
        st.sig_coly2 = ((st.sigy**2.0 * 1e-6 / (st.dist**2.0)
                         + st.sigpy**2.0 * 1e-6
                         + st.sig_mol**2.0) * st.E0*st.E0)
        st.sig_colr2 = st.sig_colx2
        st.sig_colas = 1.0/(2.0*st.sig_coly2) - 1.0/(2.0*st.sig_colx2)
        st.sig_colli = st.sig_colx2 + st.sig_coly2

    if st.type == 0:
        st.int_inc = int_BH
        st._logbuf = "anb > Calculating incoherent part via Bethe Heitler and trivial Electron contrib\n"
    elif st.type == 1:
        if st.sig_colli == 0:
            st.int_inc = int_Hub
        else:
            st.int_inc = int_HubBd
        st._logbuf = "anb > Calculating incoherent part via Hubbell and Electron contrib from [Owens]\n"
    else:
        st._logbuf = "anb > FATAL ERROR: process type doesn't exist\n"
        logging_(LOGERROR)
        sys.exit(-1)
    logging_(LOGFILE)

    st._logbuf = (f"anb > Collimation angle (E0*theta_c): {st.C}\t"
                  f"colli-divergence: {st.sig_colli}\n")
    logging_(LOGFILE)

    st.A3     = A**(-3.0)
    # Debye factor -- set a reasonable default (T=0 -> no Debye reduction)
    # The original C code reads Adebye from somewhere; it seems it's set via calc_debyered
    # In the original code, Adebye is a global set via input or calc_debyered.
    # Searching anb.c: Adebye is used but never explicitly set in input_parameter.
    # It appears to default to 0 (C global initialisation). We keep same behaviour.
    # If needed, the user can extend this.
    if not hasattr(st, 'Adebye') or st.Adebye == 0.0:
        st.Adebye = 0.0

    st.del_bx = (2.0*BEAMXY_REGION*st.sig_bx/BEAMX_STEPS if st.sig_bx > 0.0 else 100.0)
    st.del_by = (2.0*BEAMXY_REGION*st.sig_by/BEAMY_STEPS if st.sig_by > 0.0 else 100.0)

    calc_collifct(st.sig_colli * math.sqrt(0.5))
    calc_latticevec(10000)

    intsum, intdif = int_coh()

    intinc, intiel = st.int_inc(st.Zcry, DEBYE_ON)
    intamo, intael = st.int_inc(st.Zamo, DEBYE_OFF)

    if st.sigE0MeV > 0.0:
        intsum = energy_spread(st.E0, st.sigE0MeV, intsum)
        intdif = energy_spread(st.E0, st.sigE0MeV, intdif)
        intinc = energy_spread(st.E0, st.sigE0MeV, intinc)
        intamo = energy_spread(st.E0, st.sigE0MeV, intamo)

    # write output
    try:
        fhdl = open(st.outname, 'w')
    except IOError:
        st._logbuf = "Fatal Error: couldn't open output file\n"
        logging_(LOGERROR)
        sys.exit(-1)

    import struct
    def _f(x):
        """Cast to float32 like C's (float) cast before fprintf %g"""
        return struct.unpack('f', struct.pack('f', x))[0]

    for i in range(1, X_RANGE + 1):
        energy = _f(float(i) / 1000.0 * st.E0MeV)
        fhdl.write(f"{energy:g}\t{_f(intsum[i]):g}\t{_f(intdif[i]):g}\t"
                   f"{_f(intinc[i]):g}\t{_f(intamo[i]):g}\t{_f(intiel[i]):g}\t{_f(intael[i]):g}\n")
    fhdl.close()

    st.date = time.ctime().rstrip('\n')
    st._logbuf = (f"anb > intensities I written in file {st.outname}\n"
                  f"contents: photon energy(MeV) Isum Idif Iinccrystal Iincamorphous Ieleccry Ielecamo\n"
                  f"use paw, xmgr or plotdata to calc Iinc,amo=(Iinc+Ielec)(cry,amo), "
                  f"Irel_exp=(Isum+Iinc)/Iamo, Pol=Idif/(Isum+Iinc)\n"
                  f"anb > Program finished ({st.date})\n")
    logging_(LOGBOTH)
    st.file_log.close()


if __name__ == '__main__':
    main()
