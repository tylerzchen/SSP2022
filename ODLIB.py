import numpy as np
import math
import random
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D
from astropy.io import fits
import ODLIB as od

def dot(vec1, vec2):
    """
    Computes dot product of 2 1x3 vectors
    
    Args:

        vec1 (list): The first vector [x1, y1, z1]

        vec2 (list): The second vector [x2, y2, z2]

    Returns:
    
        Dot product of 2 1x3 vectors (float)
    """
    return vec1[0]*vec2[0] + vec1[1]*vec2[1]+vec1[2]*vec2[2]

def cross(vec1, vec2):
    """
    Computes cross product of 2 1x3 vectors
    
    Args:

        vec1 (list): The first vector [x1, y1, z1]

        vec2 (list): The second vector [x2, y2, z2]

    Returns:
    
        Cross product of 2 1x3 vectors (list)
    """
    return [(vec1[1] * vec2[2] - vec1[2] * vec2[1]), (vec1[2] * vec2[0] - vec1[0] * vec2[2]), (vec1[0] * vec2[1] - vec1[1] * vec2[0])]

def mag(vec):
    """
    Computes magnitude of a vector
    
    Args:
    
        vec (array): A n x 3 vector [an, bn, cn...]
    
    Returns:
    
        Magnitude of vector (float)
    """
    sum = 0
    for i in range(len(vec)):
        sum += vec[i]**2
    return np.sqrt(sum)

def rotateZ(vec, a):
    """
    Rotates a given 3D vector around the z-axis
    
    Args:
    
        vec (array): A 3D vector [x1, y1, z1]
    
        a (float): angle of rotation
    
    Returns:
    
        Rotated vector [x, y, z] (array)
    """
    rotMatrix = np.array([
        [np.cos(a), -np.sin(a), 0],
        [np.sin(a), np.cos(a), 0],
        [0, 0, 1]
    ])
    rotVec = np.matmul(rotMatrix, vec)
    return rotVec

def rotateX(vec, a):
    """
    Rotates a given 3D vector around the x-axis
    
    Args:
    
        vec (array): A 3D vector [x1, y1, z1]
    
        a (float): angle of rotation
    
    Returns:
    
        Rotated vector [x, y, z] (array)
    """
    rotMatrix = np.array([
        [1, 0, 0],
        [0, np.cos(a), -np.sin(a)],
        [0, np.sin(a), np.cos(a)]
    ])
    rotVec = np.matmul(rotMatrix, vec)
    return rotVec

def error(calc, real):
    """
    Calculates percent error between two numbers
    
    Args:

        calc (float): Calculated value

        real (float): Expected/Real value
    
    Returns:
    
        Percent error (float)
    """
    return abs((calc - real)/real) * 100

def listError(lcalc, lreal):
    """
    Calculates percent error between two lists of equal length
    
    Args:
    
        lcalc (list): Calculated values
    
        lreal (list): Expected/Real values
    
    Returns:
    
        Percent error (float)
    """
    lerrors = []
    for i in range(len(lcalc)):
        lerrors.append(abs(lcalc[i] - lreal[i])/lreal[i] * 100)
    return lerrors

def NewtonMethod(f, fprime, xi, error):
    """
    Approximates root value of real-valued function f
    
    Args:
    
        f (function): Real-valued function to be approximated
        
        fprime (function): Derivative of real-valued function
        
        xi (float): Initial guess passed as a root
    
        error (float): Approximation precision
    
    Returns:
    
        Approximated root value (float)
    """
    while True:
        x = xi - f(xi)/fprime(xi)
        if abs(x - xi) < (error):
            return x
        xi = x

def gaussdays(num):
    """
    Converts days into gaussian days
    
    Args:
    
        num (float): Number to be converted
    
    Returns: 
    
        Number converted into gaussian days (float)
    """
    return num * (365.2568983/(2 * np.pi))

def HMStoDeg(hours, mins, sec, isRad):
    """
    Decimalizes right acension
    
    Args:
    
        hours (int): hour value of right acension

        mins (int): minute value of right acension

        sec (float): second value of right acension
    
        isRad (boolean): Determines if returned value is in radians or degrees (true --> radians, false --> degrees)
    
    Returns:
    
        right acension converted into degrees or radians (float)
    """
    hours += (mins/60 + sec/3600)
    hours %= 24
    if(isRad):
        return hours * (15 * np.pi/180)
    return hours * 15

def DMStoDeg(degs, arcmin, arcsec):
    """
    Decimalizes declination
    
    Args:
    
        degs (int): degree value of declination

        arcmin (int): arcmin value of declination 

        arcsec (float): arcsec value of declination 
    
    Returns:
    
        declination converted into degrees (float) 
    """
    degs += (math.copysign(arcmin, (degs)))/60 + (math.copysign(arcsec, (degs)))/3600
    return degs

def RAdecimaltoHMS(degs):
    """
    Converts right acension in degrees to HMS format
    
    Args:
    
        degs (float): degree value of right acension 
    
    Returns:
    
        hours (int), minutes (int), seconds (float) of right acension in HMS format (tuple)
    """
    hours = degs / 15
    hoursR = hours // 1
    mins = (hours - hoursR) * 60
    minsR = mins // 1
    secs = (mins - minsR) * 60
    return ((hoursR, minsR, secs))

def DECdecimaltoDMS(deg):
    """
    Converts declination in degrees to DMS format
    
    Args:
    
        degs (float): degree value of declination 
    
    Returns:
    
        degrees (int), arcminutes (int), arcseconds (float) of declination in DMS format (tuple)
    """
    return((int(deg), int((abs(deg) % 1) * 60), (((abs(deg)%1)*60)%1)*60))

def valAng(sin, cos):
    """
    Converts a given sine and cosine into value of an angle in radians in the correct quadrant
    
    Args:
    
        sin (float): given sine value 

        cos (float): given cosine value 
    
    Returns:
    
        Value of angle in radians (float)
    """
    if np.arcsin(sin) >= 0:
        return np.arccos(cos)
    return 2*np.pi - np.arccos(cos)

def calcTaus(obsTimes):
    """
    Calculates time differences of observation in gaussian days
    
    Args:
    
        obsTimes (array): Time of observations in julian days
    
    Returns:
    
        Time differences in gaussian days [tau1, tau3, tau] (list)
    """
    k = 0.0172020989484
    tau1 = k * (obsTimes[0] - obsTimes[1])
    tau3 = k * (obsTimes[2] - obsTimes[1])
    tau = k * (obsTimes[2] - obsTimes[0])
    return [tau1, tau3, tau]

def calcD(ra1, ra2, ra3, dec1, dec2, dec3, suns):
    """
    Calculates D coefficients and rhohats
    
    Args:
    
        ra1 (list): right acension of asteroid for first observation in HMS [hours, mins, secs]

        ra2 (list): right acension of asteroid for second observation in HMS [hours, mins, secs]

        ra3 (list): right acension of asteroid for third observation in HMS [hours, mins, secs]

        dec1 (list): declination of asteroid for first observation in DMS [degrees, arcmins, arcsecs]

        dec2 (list): declination of asteroid for second observation in DMS [degrees, arcmins, arcsecs]

        dec3 (list): declination of asteroid for third observation in DMS [degrees, arcmins, arcsecs]

        suns (array): Sun to Earth (SBO) vectors in the equitorial plane in AU [[R1], [R2], [R3]]
    
    Returns:
    
        D0 (float): D0 coefficient
        
        D1s (list): List of D1 coefficients
        
        D2s (list): List of D2 coefficients
        
        D3s (list): List of D3 coefficients
        
        rhohats (array): List of rhohats
    """
    ra1deg = HMStoDeg(ra1[0], ra1[1], ra1[2], True) 
    ra2deg = HMStoDeg(ra2[0], ra2[1], ra2[2], True) 
    ra3deg = HMStoDeg(ra3[0], ra3[1], ra3[2], True)
    ras = np.array([ra1deg, ra2deg, ra3deg])
    dec1 = DMStoDeg(dec1[0], dec1[1], dec1[2]) * (np.pi/180)
    dec2 = DMStoDeg(dec2[0], dec2[1], dec2[2]) * (np.pi/180)
    dec3 = DMStoDeg(dec3[0], dec3[1], dec3[2]) * (np.pi/180)
    decs = np.array([dec1, dec2, dec3])
    rhohats = np.array([[np.cos(ras[i]) * np.cos(decs[i]), np.sin(ras[i]) * np.cos(decs[i]), np.sin(decs[i])] for i in range(3)])
    D0 = dot(rhohats[0], cross(rhohats[1], rhohats[2]))
    D1s = []
    D2s = []
    D3s = []
    for i in range(3):
        D1s.append(dot(cross(suns[i], rhohats[1]), rhohats[2]))
        D2s.append(dot(cross(rhohats[0], suns[i]), rhohats[2]))
        D3s.append(dot(rhohats[0], cross(rhohats[1], suns[i])))
    return D0, D1s, D2s, D3s, rhohats

def SEL(taus, Sun2, rhohat2, Ds):
    """
    Calculates real roots of possible rho2s and r2s using scalar equation of lagrange
    
    Args:
    
        taus (list): time differences of observation in julian days [tau1, tau2, tau3]

        Sun2 (list): sun vector of observation 2 [x, y, z]

        rhohat2 (list): rhohat vector of observation 2 [x, y, z]

        Ds (list): D coefficients for observation 2 [D0, D21, D22, D32]
    
    Returns:
    
        roots (list): possible values of R2 [float1, float2, float3]
        
        rhos (list): corresponding possible range values rho2 [float1, float2, float3]
    """
    roots = []
    rhos = []
    possRoots = []

    A1 = taus[1]/taus[2]
    A3 = -taus[0]/taus[2]
    B1 = (A1/6) * (taus[2]**2 - taus[1]**2)
    B3 = (A3/6) * (taus[2]**2 - taus[0]**2)
    
    A = (A1 * Ds[1] - Ds[2] + A3 * Ds[3])/(-Ds[0])
    B = (B1 * Ds[1] + B3 * Ds[3])/(-Ds[0])
    E = -2 * (dot(Sun2, rhohat2))
    F = dot(Sun2, Sun2)
    
    a = -(A**2 + A * E + F)
    b = -(2 * A * B + B * E)
    c = -B**2
    preroots = np.polynomial.polynomial.polyroots([c, 0, 0, b, 0, 0, a, 0, 1])
    for i in preroots:
        if np.isreal(i) and i > 0:
            possRoots.append(np.real(i))
    for i in range(len(possRoots)):
        possibleRho = A + B/possRoots[i]**3
        if possibleRho > 0:
            roots.append(possRoots[i])
            rhos.append(possibleRho)
    return roots, rhos

def fg(tau, r2, r2dot, flag):
    """
    Calculates f and g coefficients
    
    Args:
    
        tau (float): gaussian time interal for current iteration

        r2 (list): position vector for asteroid at observation 2 [x, y, z]

        r2dot (list): velocity vector for asteroid at observation 2 [vx, vy, vz]

        flag (int): determines whether second, third, or fourth degree taylor series is used

    Returns:
    
        f (float): f coefficient
        
        g (float): g coefficient
    """
    if flag == 2:
        u = 1/r2**3
        f = 1 - .5 * u * tau**2
        g = tau
    else:
        u = 1/mag(r2)**3 
        z = dot(r2, r2dot)/mag(r2)**2
        q = dot(r2dot, r2dot)/mag(r2)**2 - u
        if flag == 4:
            f = 1 - .5 * u * tau**2 + .5 * u * z * tau**3 + (1/24) * (3 * u * q - 15 * u * z ** 2 + u**2) * tau**4
            g = tau - (1/6) * u * tau**3 + (1/4) * u * z * tau**4
        elif flag == 3:
            f = 1 - .5 * u * tau**2 + .5 * u * z * tau**3
            g = tau - (1/6) * u * tau**3
    return f, g

def calcC(f1, f3, g1, g3):
    """
    Calculates C coefficients
    
    Args:
    
        f1 (float): f1 coefficient
        
        f3 (float): f3 coefficient
        
        g1 (float): g1 coefficient
        
        g3 (float): g3 coefficient
    
    Returns: 
        
        C coefficients (list)
    """
    c1 = g3/(f1 * g3 - g1 * f3)
    c2 = -1
    c3 = -g1/(f1 * g3 - g1 * f3)
    return [c1, c2, c3]

def rhoMags(Cs, D0, D1s, D2s, D3s):
    """
    Calculates magnitude of rho vectors
    
    Args:
        
        Cs (list): C coefficients
        
        D0 (float): D0 coefficient
        
        D1s (list): D1 coefficients
        
        D2s (list): D2 coefficients
        
        D3s (list): D3 coefficients
    
    Returns: 
        
        Magnitudes of rho vectors (list)
    """
    rho1Mag = (Cs[0] * D1s[0] + Cs[1] * D1s[1] + Cs[2] * D1s[2])/(Cs[0] * D0)
    rho2Mag = (Cs[0] * D2s[0] + Cs[1] * D2s[1] + Cs[2] * D2s[2])/(Cs[1] * D0)
    rho3Mag = (Cs[0] * D3s[0] + Cs[1] * D3s[1] + Cs[2] * D3s[2])/(Cs[2] * D0)
    return [rho1Mag, rho2Mag, rho3Mag]

def rhoVectors(rhomags, rhohats):
    """
    Calculates rho vectors
    
    Args:
        
        rhomags (list): Magnitude of rho vectors
        
        rhohats (array): Rhohats for corresponding rho vectors [[rhohat1], [rhohat2], [rhohat3]]
    
    Returns:
    
        rho vectors (array)
    """
    rho1 = rhohats[0] * rhomags[0]
    rho2 = rhohats[1] * rhomags[1]
    rho3 = rhohats[2] * rhomags[2]
    return np.array([rho1, rho2, rho3])

def rVectors(rhovectors, Rs):
    """
    Calculates position vectors
    
    Args:
        
        rhovectors (array): array of rhovectors [[rho1], [rho2], [rho3]]
        
        Rs (array): Sun vectors
    
    Returns:
    
        position vectors (array)
    """
    r1 = rhovectors[0] - Rs[0]
    r2 = rhovectors[1] - Rs[1]
    r3 = rhovectors[2] - Rs[2]
    return np.array([r1, r2, r3])

def calcr2Dot(rs, f1, f3, g1, g3):
    """
    Calculates velocity vector for observation 2
    
    Args:
    
        rs (array): position vectors ([r1, r2, r3])
        
        f1 (float): f1 coefficient
        
        f3 (float): f3 coefficient
        
        g1 (float): g1 coefficient
        
        g3 (float): g3 coefficient
    
    Returns:
        
        velocity vector for observation 2 (list)
    """
    return (-f3/(-f3 * g1 + g3 * f1)) * rs[0] + (f1/(f1 * g3 - f3 * g1)) * rs[2]

def timeCorrection(obsTimes, rhomags):
    """
    Corrects observation times for light-travel time
    
    Args:
    
        obsTimes (list): list of observation times of asteroid [t1, t2, t3]
        
        
        rhomags (array): array of magnitudes of rhos
    
    Returns:
    
        Corrected observation times [t1C, t2C, t3C] (list)
    """
    cAU = 173.144643267 
    t1 = obsTimes[0] - rhomags[0]/cAU
    t2 = obsTimes[1] - rhomags[1]/cAU
    t3 = obsTimes[2] - rhomags[2]/cAU
    return [t1, t2, t3]

def listDiff(l1, l2):
    """
    Calculates differences between elements of two lists of length 3
    
    Args:
    
        l1 (list): first list
        
        l2 (list): second list
    
    Returns:
    
        Differences of each element in list [diff1, diff2, diff3] (list)
    """
    return [l1[0] - l2[0], l1[1] - l2[1], l1[2] - l2[2]]

def rotateEcliptic(vec):
    """
    Rotates a given 3D vector around the ecliptic
    
    Args:
    
        vec (array): A 3D vector [x1, y1, z1]
    
    Returns:
    
        Rotated vector [x, y, z] (array)
    """
    negEps = (23.4374) * (np.pi/180)
    rotMatrix = np.array([
        [1, 0, 0],
        [0, np.cos(negEps), np.sin(negEps)],
        [0, -np.sin(negEps), np.cos(negEps)]
    ])
    rotVec = np.matmul(rotMatrix, vec)
    return rotVec

def MOG(obsTimes, ras, decs, Rs):
    """
    Calculates asteroid's position and velocity vector for central observation in ecliptic rectangular coordinates in AU and AU/day
    
    Args:
        
        obsTimes (array): time of middle observations in julian days [t1, t2, t3]
        
        ras (array): right acensions in HMS [[ra1], [ra2], [ra3]]
        
        decs (array): declinations of asteroid in DMS [[dec1], [dec2], [dec3]]
        
        Rs (array): Sun vector from JPL Horizons (AU, equatorial cartesian, J200, apparent states)
    
    Returns:
    
        position and velocity vectors for asteroids seconds observation ([r2, r2dot]) (array)
    """
    taus = calcTaus(obsTimes)
    D0, D1s, D2s, D3s, rhohats = calcD(ras[0], ras[1], ras[2], decs[0], decs[1], decs[2], Rs)
    DSELs = [D0, D2s[0], D2s[1], D2s[2]]
    r2MagGuesses, rho2MagGuesses = SEL(taus, Rs[1], rhohats[1], DSELs)
    if len(r2MagGuesses) > 1:
        print(r2MagGuesses)
        print(rho2MagGuesses)
        val = int(input("Enter index of desired position for guesses: "))
        while True:
            if val > len(r2MagGuesses) - 1:
                val = int(input("Please enter valid index: "))
            else:
                break
        r2Mag = r2MagGuesses[val]
    else:
        r2Mag = r2MagGuesses[0]
    f1, g1 = fg(taus[0], r2Mag, 0, 2)
    f3, g3 = fg(taus[1], r2Mag, 0, 2)
    Cs = calcC(f1, f3, g1, g3)
    rhomags = rhoMags(Cs, D0, D1s, D2s, D3s)
    rhos = rhoVectors(rhomags, rhohats)
    rs = rVectors(rhos, Rs)
    r2 = rs[1]
    r2dot = calcr2Dot(rs, f1, f3, g1, g3)
    r2prev = [0, 0, 0]
    r2dotprev = [0, 0, 0]
    while mag(listDiff(r2, r2prev)) > 1e-10 and mag(listDiff(r2dot, r2dotprev)) > 1e-10:
        r2prev = r2
        r2prevdot = r2dot
        timesCorrected = timeCorrection(obsTimes, rhomags)
        taus = calcTaus(timesCorrected)
        f1, g1 = fg(taus[0], r2, r2dot, 4)
        f3, g3 = fg(taus[1], r2, r2dot, 4)
        Cs = calcC(f1, f3, g1, g3)
        rhomags = rhoMags(Cs, D0, D1s, D2s, D3s)
        rhos = rhoVectors(rhomags, rhohats)
        rs = rVectors(rhos, Rs)
        r2 = rs[1]
        r2dot = calcr2Dot(rs, f1, f3, g1, g3)
    rho2 = rotateEcliptic(rhos[1])
    r2 = rotateEcliptic(r2)
    r2dot = rotateEcliptic(r2dot)
    data = np.array([r2, r2dot])
    for i in range(3):
        data[1, i] = data[1, i] * ((2 * np.pi)/365.2568983)
    return data, rho2
    
def angMomentum(r, rdot):
    """
    Computes specific angular momentum
    
    Args:
    
        r (list): position vector [x, y, z]

        rdot (list): velocity vector [vx, vy, vz]
    
    Returns:
    
        Specific angular momentum vector [hx, hy, hz] (list)
    """
    angMomentum = cross(
        [r[0], r[1], r[2]],
        [rdot[0], rdot[1],rdot[2]]
    )
    return angMomentum

def semimajor(r, rdot):
    """
    Calculates the semimajor axis of elliptical orbit (vis visa equation)
    
    Args:
    
        r (list): position vector [x, y, z]

        rdot (list): velocity vector [vx, vy, vz]
    
    Returns:
    
        Semimajor axis of ellipitical orbit (float)
    """
    rmag = mag(r)
    vmagSqr = dot(rdot, rdot)
    a = 1/((2/rmag) - vmagSqr)
    return a

def eccentricity(r, rdot):
    """
    Calculates the eccentricity of elliptical orbit
    
    Args:
    
        r (list): position vector [x, y, z]

        rdot (list): velocity vector [vx, vy, vz]
    
    Returns:
    
        eccentricity of elliptical orbit (float)
    """
    h = angMomentum(r, rdot)
    e = np.sqrt(1 - (mag(h)**2/semimajor(r, rdot)))
    return e

def inclination(r, rdot):
    """
    Calculates inclination 
    
    Args:
    
        r (list): position vector [x, y, z]

        rdot (list): velocity vector [vx, vy, vz]
    
    Returns:
    
        Inclination in degrees (float)
    """
    h = angMomentum(r, rdot)
    i = np.arctan((np.sqrt(h[0]**2 + h[1]**2))/h[2])
    return i * (180/np.pi)

def longNode(r, rdot):
    """
    Calculates longitude of ascending nodes
    
    Args:
    
        r (list): position vector [x, y, z]
    
        rdot (list): velocity vector [vx, vy, vz]
    
    Returns:
    
        Longitude of ascending nodes in degrees (float)
    """
    h = angMomentum(r, rdot)
    i = inclination(r, rdot) * (np.pi/180)
    cos = -h[1]/(mag(h) * np.sin(i))
    sin = h[0]/(mag(h) * np.sin(i))
    return valAng(sin, cos) * (180/np.pi)

def trueAnomaly(data):
    """
    Calculates true anomaly (for poliastro visualization)
    
    Args:
    
        data (array): position and velocity vector for asteroids orbit in AU and days        
        
    Returns:
    
        True anomaly (float) in degrees
    """
    r = data[0].copy()
    rdot = data[1].copy()
    for i in range(len(rdot)):
        rdot[i] = gaussdays(rdot[i])
    h = angMomentum(r, rdot)
    om = longNode(r, rdot) * (np.pi/180)
    i = inclination(r, rdot) * (np.pi/180)
    e = eccentricity(r, rdot)
    a = semimajor(r, rdot)
    cos = (1/e)*((a * (1 - e**2))/mag(r) - 1)
    sin = (a * (1 - e**2))/(mag(h) * e) * (dot(r, rdot)/mag(r))
    V = valAng(sin, cos) * (180/np.pi)
    return V

def argPeri(r, rdot):
    """
    Calculates argument of perihelion
    
    Args:
    
        r (list): position vector [x, y, z]
    
        rdot (list): velocity vector [vx, vy, vz]
    
    Returns:
    
        Argument of perihelion in degrees (float)
    """
    h = angMomentum(r, rdot)
    om = longNode(r, rdot) * (np.pi/180)
    i = inclination(r, rdot) * (np.pi/180)
    e = eccentricity(r, rdot)
    a = semimajor(r, rdot)
    cos = (r[0] * np.cos(om) + r[1] * np.sin(om))/(mag(r))
    sin = (r[2])/(mag(r) * np.sin(i))
    U = valAng(sin, cos) * (180/np.pi)
    cos = (1/e)*((a * (1 - e**2))/mag(r) - 1)
    sin = (a * (1 - e**2))/(mag(h) * e) * (dot(r, rdot)/mag(r))
    V = valAng(sin, cos) * (180/np.pi)
    peri = U - V
    if peri < 0:
        return 360 + peri
    return peri

def meanAnomaly(r, rdot):
    """
    Calculates mean anomaly
    
    Args:
    
        r (list): position vector [x, y, z]
    
        rdot (list): velocity vector [vx, vy, vz]
    
    Returns:
    
        Mean anomaly in degrees (float)
    """
    h = angMomentum(r, rdot)
    e = eccentricity(r, rdot)
    a = semimajor(r, rdot)
    cos = (1/e)*((a * (1 - e**2))/mag(r) - 1)
    sin = (a * (1 - e**2))/(mag(h) * e) * (dot(r, rdot)/mag(r))
    V = valAng(sin, cos)
    cos = (e + np.cos(V))/(1 + e * np.cos(V))
    sin = (mag(r) * np.sin(V))/(a * np.sqrt(1 - e**2))
    E = valAng(sin, cos)
    M = E - e * np.sin(E)
    return M * (180/np.pi)

def periPassage(r, rdot, obsTime):
    """
    Calculates time of perihelion passage
    
    Args:
    
        r (list): position vector [x, y, z]
    
        rdot (list): velocity vector [vx, vy, vz]
        
        obsTime (float): time of observation in julian days
    
    Returns:
    
        Time of perihelion passage (float)
    """
    M = meanAnomaly(r, rdot)
    n = 0.0172020989484/(semimajor(r, rdot))**(3/2)
    T = obsTime - (M * np.pi/180)/n
    return T

def orbitalElements(data, obsTime):
    """
    Calculates all orbital elements of asteroids orbit
    
    Args:
        
        data (array): position and velocity vector for asteroids orbit in AU and days
        
        obsTime (float): time of observation in julian days
        
    Returns:
        
        h (list): specific angular momentum vector [hx, hy, hz]
        
        a (float): semimajor axis of elliptical orbit
        
        e (float): eccentricity of elliptical orbit
        
        i (float): inclination in degrees
        
        om (float): longitude of ascending nodes in degrees
        
        w (float): argument of perihelion in degrees
        
        m (float): mean anomaly in degrees
        
        T (float): time of perihelion passage
    """
    r = data[0].copy()
    rdot = data[1].copy()
    for i in range(len(rdot)):
        rdot[i] = gaussdays(rdot[i])
    h = angMomentum(r, rdot)
    a = semimajor(r, rdot)
    e = eccentricity(r, rdot)
    i = inclination(r, rdot)
    om = longNode(r, rdot)
    w = argPeri(r, rdot)
    m = meanAnomaly(r, rdot)
    T = periPassage(r, rdot, obsTime)
    
    return h, a, e, i, om, w, m, T

def compareElements(data, obsTime, JPLVals):
    """
    Prints calculated, JPL, and % error of asteroid elements with a period correction for the time of perihelion passage JPL
    
    Args:
    
        data (array): position and velocity vector for asteroids orbit in AU and days
        
        obsTime (float): time of observation in julian days
        
        JPLVals (list): JPL Horizon asteroid element values
    """
    h, a, e, i , om, w, m, T = orbitalElements(data, obsTime)
    print("Semimajor Axis Calculated: ", a)
    print("Semimajor Axis JPL: ", JPLVals[0])
    print("Semimajor Axis % Error: ", error(a, JPLVals[0]))
    print()
    print("Eccentricity Calculated: ", e)
    print("Eccentricity JPL: ", JPLVals[1])
    print("Eccentricity % Error: ", error(e, JPLVals[1]))
    print()
    print("Inclination Calculated: ", i)
    print("Inclination JPL: ", JPLVals[2])
    print("Inclination % Error: ", error(i, JPLVals[2]))
    print()
    print("Longitude of Ascending Nodes Calculated: ", om)
    print("Longitude of Ascending Nodes JPL: ", JPLVals[3])
    print("Longitude of Ascending Nodes % Error: ", error(om, JPLVals[3]))
    print()
    print("Argument of Perihelion Calculated: ", w)
    print("Argument of Perihelion JPL: ", JPLVals[4])
    print("Argument of Perihelion % Error: ", error(w, JPLVals[4]))
    print()
    print("Mean Anomaly Calculated: ", m)
    print("Mean Anomaly JPL: ", JPLVals[5])
    print("Mean Anomaly % Error: ", error(m, JPLVals[5]))
    print()
    Tast =  (JPLVals[0]**(3/2) * 365.2568983) 
    JPLVals[6] = obsTime - ((obsTime - JPLVals[6])%Tast)
    print("Time of Perihelion Passage Calculated: ", T)
    print("Time of Perihelion Passage Corrected JPL: ", JPLVals[6])
    print("Time of Perihelion Passage % Error: ", error(T, JPLVals[6]))
    
def ephemeris(data, R, obsTime, dataObsTime):
    """
    Returns ephemeris of asteroid at desired observation time
    
    Args:
    
        data (array): known position and velocity vectors for asteroids orbit in AU and days
        
        R (array): sun distance vector in equitorial plane of desired ephemeris date
        
        obsTime (float): desired observation time in julian days
        
        dataObsTime (float): data observation time in julian days
        
    Returns:
        
        Right acension (HMS format) and declination (DMS) of asteroid at desired observation time (array)
    """
    r = data[0].copy()
    rdot = data[1].copy()
    for i in range(len(rdot)):
        rdot[i] = gaussdays(rdot[i])
    a = semimajor(r, rdot)
    e = eccentricity(r, rdot)
    i = inclination(r, rdot) * (np.pi/180)
    om = longNode(r, rdot) * (np.pi/180)
    w = argPeri(r, rdot) * (np.pi/180)
    eps = (23.4374) * (np.pi/180)
    T = periPassage(r, rdot, dataObsTime)
    k = 0.0172020989484
    n = k * np.sqrt(1/a**3)
    M = n * (obsTime - T)
    f = lambda E: M - E + e * np.sin(E)
    fprime = lambda E: -1 + e * np.cos(E)
    E = NewtonMethod(f, fprime, M, 10**(-10))
    pos = np.array([a * np.cos(E) - a * e, a * np.sqrt(1 - e**2) * np.sin(E), 0])
    pos = rotateZ(pos, w)
    pos = rotateX(pos, i)
    pos = rotateZ(pos, om)
    pos = rotateX(pos, eps)
    rho = pos + R
    rhohat = rho/mag(rho)
    dec = np.arcsin(rhohat[2])
    cos = (rhohat[0]/np.cos(dec))
    sin = (rhohat[1]/np.cos(dec))
    ra = valAng(sin, cos) * (180/np.pi)
    dec = DECdecimaltoDMS(dec * (180/np.pi))
    ra = RAdecimaltoHMS(ra)
    eph = np.array([ra, dec])
    return eph

def ephemerides(data, Rs, obsTimes, dataObsTime):
    """
    Prints ephemerides for desired observation times
    
    Args:
    
        data (array): known position and velocity vectors for asteroids orbit in AU and days
        
        Rs (array): sun distance vectors in equitorial plane of desired ephemeris dates
      
        obsTimes (float): desired observation times in julian days
        
        dataObsTime (float): data observation time in julian days
    
    Returns: 
    
        Prints out ephemerides' time (julian days), RAs (HMS), and Decs (DMS)
    """
    for i in range(len(obsTimes)):
        eph = ephemeris(data, Rs[i], obsTimes[i], dataObsTime)
        print("Time: ", obsTimes[i])
        print("RA and Dec: ")
        print(eph)

def rmsRADec(path):
    """
    Calculates the root mean square of an images ra and dec utilizing corr.fits file
    
    Args:
    
        path (String): corr.fits file path
        
    Returns:
    
        raRMS (float): right acension root mean square (degrees)
        
        decRMS (float): declination root mean square (degrees)
    """
    table = fits.open(path)[1].data
    field_dec, index_dec = table.field_dec[:], table.index_dec[:]
    field_ra, index_ra = table.field_ra[:], table.index_ra[:]
    decRMS = np.sqrt((1/int(table.shape[0])) * sum([(field_dec[i] - index_dec[i])**2 for i in range(int(table.shape[0]))]))
    raRMS = np.sqrt((1/int(table.shape[0])) * sum([(field_ra[i] - index_ra[i])**2 for i in range(int(table.shape[0]))]))
    return decRMS, raRMS

def monteHistogram(dataList, realVal, graphLine, units, xlabel, title):
    """
    Displays formatted histogram for a monte carlo distribution and prints percent error of mean compared to real value. Period corrects JPL time of perihelion passage. 
    
    Args:
    
        dataList (list): list of values generated through monte carlo simulation
        
        realVal (float): real value to be compared to monte carlo data mean
        
        graphLine (boolean): determines if line is graphed on histogram for real and mean value
        
        units (String): units of value being simulated
        
        xlabel (String): histogram x-label
        
        title (String): histogram title
    """
    fig1, axs1 = plt.subplots(1, 1, figsize = (10, 7), tight_layout = True)
    for s in ['top', 'bottom', 'left', 'right']:
        axs1.spines[s].set_visible(False)
    axs1.xaxis.set_ticks_position('none')
    axs1.yaxis.set_ticks_position('none')
    axs1.xaxis.set_tick_params(pad = 5)
    axs1.yaxis.set_tick_params(pad = 10)
    axs1.grid(b = True, color ='grey', linestyle ='-.', linewidth = 0.5, alpha = 0.6)   
    N, bins, patches = axs1.hist(dataList, bins = 20)
    fracs = ((N**(1 / 5)) / N.max())
    norm = colors.Normalize(fracs.min(), fracs.max())
    for thisfrac, thispatch in zip(fracs, patches):
        color = plt.cm.viridis(norm(thisfrac))
        thispatch.set_facecolor(color)
    if graphLine:
        axs1.axvline(realVal, color = 'k', linestyle = 'dashed', linewidth = 1)
        axs1.axvline(sum(dataList)/len(dataList), color = 'b', linewidth = 1)
        customLines = [Line2D([0], [0], color = 'black', linestyle = 'dashed', lw = 1), Line2D([0], [0], color = 'blue', lw = 1)]
        axs1.legend(customLines, ["True Value", "Mean"])
    axs1.set_xlabel(xlabel + " " + units)
    axs1.set_ylabel("Count")
    plt.title(title)
    plt.show()
    mean = sum(dataList)/len(dataList)
    sumSq = sum([(dataList[i] - mean)**2 for i in range(len(dataList))])
    if realVal != 0:
        print(xlabel, " Mean: ", mean)
        print(xlabel, " SDOM: ", np.sqrt(sumSq/len(dataList)))
        print(xlabel, " % error w/JPL : ", od.error(mean, realVal))

def monteCarlo(path1, path2, path3, obsTimes, ras, decs, Rs, JPLVals, graphLine, n):
    """
    Generates 6 monte carlo distribution histograms for orbital elements of asteroid and for time of perihelion passage (period corrected JPL)
    
    Args:
    
        path1 (String): corr.fits file path for observation 1
        
        path2 (String): corr.fits file path for observation 2
        
        path3 (String: corr.fits file path for observation 3
        
        obsTimes (array): time of middle observations in julian days [t1, t2, t3]
        
        ras (array): right acensions in HMS [[ra1], [ra2], [ra3]]
        
        decs (array): declinations of asteroid in DMS [[dec1], [dec2], [dec3]]
        
        Rs (array): Sun vector from JPL Horizons (AU, equatorial cartesian, J200, apparent states)
        
        JPLVals (list): JPL Horizon asteroid element values
        
        graphLine (boolean): determines if JPL and mean vertical lines are displayed on histogram graph
        
        n (int): number of iterations for monte carlo simulation
    """
    decRMS1, raRMS1 = rmsRADec(path1)
    decRMS2, raRMS2 = rmsRADec(path2)
    decRMS3, raRMS3 = rmsRADec(path3)
    rasCopy = np.copy(ras)
    decsCopy = np.copy(decs)
    rasDeg = [HMStoDeg(rasCopy[i, 0], rasCopy[i, 1], rasCopy[i, 2], False) for i in range(3)]
    decsDeg = [DMStoDeg(decsCopy[i, 0], decsCopy[i, 1], decsCopy[i, 2]) for i in range(3)]
    rasHMS = np.zeros((3, 3))
    hList, aList, eList, iList, omList, wList, mList, TList = [[] for i in range(8)]
    for i in range(n):
        decsNew = [decsDeg[0] + np.random.normal() * decRMS1, decsDeg[1] + np.random.normal() * decRMS2, decsDeg[2] + np.random.normal() * decRMS3]
        rasNew = [rasDeg[0] + np.random.normal() * raRMS1, rasDeg[1] + np.random.normal() * raRMS2, rasDeg[2] + np.random.normal() * raRMS3]
        decsDMS = [list(DECdecimaltoDMS(decsNew[0])), list(DECdecimaltoDMS(decsNew[1])), list(DECdecimaltoDMS(decsNew[2]))]
        rasHMS = [list(RAdecimaltoHMS(rasNew[0])), list(RAdecimaltoHMS(rasNew[1])), list(RAdecimaltoHMS(rasNew[2]))]
        dataNew, rho2 = MOG(obsTimes, rasHMS, decsDMS, Rs)
        h, a, e, inc, om, w, m, T = orbitalElements(dataNew, obsTimes[1])
        aList.append(a)
        eList.append(e)
        iList.append(inc)
        omList.append(om)
        wList.append(w)
        mList.append(m)
        TList.append(T)
    Tast =  (JPLVals[0]**(3/2) * 365.2568983) 
    JPLVals[6] = obsTimes[1] - ((obsTimes[1] - JPLVals[6])%Tast)
    monteHistogram(aList, JPLVals[0], graphLine, "(AU)", "Semi-major Axis", "Monte Carlo Distribution of Semi-major Axis")
    monteHistogram(eList, JPLVals[1], graphLine, "", "Eccentricity", "Monte Carlo Distribution of Eccentricity")
    monteHistogram(iList, JPLVals[2], graphLine, "(degrees)", "Inclination", "Monte Carlo Distribution of Inclination")
    monteHistogram(omList, JPLVals[3], graphLine, "(degrees)", "Longitude of Ascending Nodes", "Monte Carlo Distribution of Longitude of Ascending Nodes")
    monteHistogram(wList, JPLVals[4], graphLine, "(degrees)", "Argument of Perihelion", "Monte Carlo Distribution of Argument of Perihelion")
    monteHistogram(mList, JPLVals[5], graphLine, "(degrees)", "Mean Anomaly", "Monte Carlo Distribution of Mean Anomaly")
    monteHistogram(TList, JPLVals[6], graphLine, "(JD)", "Time of Perihelion Passage", "Monte Carlo Distribution of Time of Perihelion Passage")
def futureMean(a, T, obsTime):
    """
    Calculates future mean anomaly for desired date
    
    Args:
        
        a (float): Semimajor axis of orbit
        
        T (float): Time of perihelion passage of orbit
        
        obsTime (float): Desired date in julian days
    
    Returns:
    
        Future mean anomaly of asteroid in degrees (float)
    """
    k = 0.0172020989484
    n = k * np.sqrt(1/a**3)
    M = n * (obsTime - T)
    return M * (180/np.pi)

def parseODTxt(path):
    """
    Parses .txt file for observation times, right acensions, declinations, and sun vectors
    
    Args:
    
        path (String): .txt file path
    
    Returns:
    
        obsTimes (list): array of observation times
        
        ras (array): right acensions of asteroid
        
        decs (array): declinations of asteroid
        
        Rs (array): sun vectors
    """
    fileData = np.loadtxt(path, dtype = float, delimiter =", ")
    obsTimes = [fileData[0, 0], fileData[1, 0], fileData[2, 0]]
    ras = np.array([
        [fileData[0, 1], fileData[0, 2], fileData[0, 3]],
        [fileData[1, 1], fileData[1, 2], fileData[1, 3]],
        [fileData[2, 1], fileData[2, 2], fileData[2, 3]]
    ])
    decs = np.array([
        [fileData[0, 4], fileData[0, 5], fileData[0, 6]],
        [fileData[1, 4], fileData[1, 5], fileData[1, 6]],
        [fileData[2, 4], fileData[2, 5], fileData[2, 6]]
    ])
    Rs = np.array([
        [fileData[0, 7], fileData[0, 8], fileData[0, 9]],
        [fileData[1, 7], fileData[1, 8], fileData[1, 9]],
        [fileData[2, 7], fileData[2, 8], fileData[2, 9]]
    ])
    return obsTimes, ras, decs, Rs

def parseSunTxt(path):
    """
    Parses .txt file for observation time and sun vector data
    
    Args:
    
        path (string): .txt file path
    
    Returns:
    
        futureTimes (list): list of observation times for asteroid (future for ephemeris generation)
        
        futureRs (array): array of sun vectors (future for ephemeris generation) formatted as [[r1], [r2],...[rn]]
    """
    fileData = np.loadtxt(path, dtype = float, delimiter =", ")
    futureTimes = []
    [futureTimes.append(time) for time in fileData[:,0]]
    fileShape = fileData.shape
    futureRs = np.zeros((fileShape[0], fileShape[1] - 1))
    for i in range(fileData.shape[0]):
        for j in range (1, fileData.shape[1]):
            futureRs[i, j - 1] = fileData[i, j]
    return futureTimes, futureRs