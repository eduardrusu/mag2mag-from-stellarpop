"""
A module to compute cosmological distances, including:
    comoving_distance (Dc)
    angular_diameter_distance (Da)
    luminosity_distance (Dl)
    comoving_volume (volume)

"""

import warnings
from math import exp, log10, pi, sin, sinh, sqrt

from scipy import integrate

warnings.warn("Default cosmology is Om=0.3,Ol=0.7,h=0.7,w=-1 and distance units are Mpc!", ImportWarning)

c = 299792458.0
G = 4.3e-6


class Distance:
    def __init__(self, cosmo=[0.3, 0.7, 0.7]):
        self.OMEGA_M = cosmo[0]
        self.OMEGA_L = cosmo[1]
        self.h = cosmo[2]
        self.w = -1.0
        self.wpars = None
        self.Dc = self.comoving_distance
        self.Dt = self.comoving_transverse_distance
        self.Dm = self.comoving_transverse_distance
        self.Da = self.angular_diameter_distance
        self.Dl = self.luminosity_distance
        self.dm = self.distance_modulus
        self.volume = self.comoving_volume

    def set(self, cosmo):
        self.OMEGA_M = cosmo[0]
        self.OMEGA_L = cosmo[1]
        self.h = cosmo[2]

    def reset(self):
        self.OMEGA_M = 0.3
        self.OMEGA_L = 0.7
        self.h = 0.7
        self.w = -1.0

    def age(self, z):
        om = self.OMEGA_M
        ol = self.OMEGA_L
        ok = 1.0 - om - ol

        def f(zp, om, ol, ok):
            return (om / zp + ok + ol * zp**2) ** -0.5

        return (9.778 / self.h) * integrate.romberg(f, 1e-300, 1 / (1.0 + z), (om, ol, ok))

    def comoving_distance(self, z1, z2=0.0):
        if z2 < z1:
            z1, z2 = z2, z1

        def fa(z):
            def wa(z, wpars):
                return (1.0 + self.w(z, wpars)) / (1.0 + z)

            return exp(3.0 * integrate.romberg(wa, 0, z))

        if isinstance(self.w, type(self.comoving_distance)) or isinstance(self.w, type(fa)):

            def f(z, om, ol, ok):
                return (om * (1.0 + z) ** 3 + ok * (1.0 + z) ** 2 + ol * fa(z)) ** -0.5

        elif self.w != -1.0:

            def f(z, om, ol, ok):
                return (om * (1.0 + z) ** 3 + ok * (1.0 + z) ** 2 + ol * (1.0 + z) ** (3.0 * (1.0 + self.w))) ** -0.5

        else:

            def f(z, om, ol, ok):
                return (om * (1.0 + z) ** 3 + ok * (1.0 + z) ** 2 + ol) ** -0.5

        om = self.OMEGA_M
        ol = self.OMEGA_L
        ok = 1.0 - om - ol
        #        return (c/self.h)*integrate.romberg(f,z1,z2,(om,ol,ok))/1e5
        return (c / self.h) * integrate.quad(f, z1, z2, (om, ol, ok))[0] / 1e5

    def comoving_transverse_distance(self, z1, z2=0.0):
        dc = 1e5 * self.comoving_distance(z1, z2) / (c / self.h)
        ok = 1.0 - self.OMEGA_M - self.OMEGA_L
        if ok > 0:
            dtc = sinh(sqrt(ok) * dc) / sqrt(ok)
        elif ok < 0:
            ok *= -1.0
            dtc = sin(sqrt(ok) * dc) / sqrt(ok)
        else:
            dtc = dc
        return (c / self.h) * dtc / 1e5

    def angular_diameter_distance(self, z1, z2=0.0):
        if z2 < z1:
            z1, z2 = z2, z1
        return self.comoving_transverse_distance(z1, z2) / (1.0 + z2)

    def luminosity_distance(self, z):
        return (1.0 + z) * self.comoving_transverse_distance(z)

    def comoving_volume(self, z1, z2=0.0):
        if z2 < z1:
            z1, z2 = z2, z1

        def f(z, om, ol, ok):
            return (self.comoving_distance(0.0, z) ** 2) / ((om * (1.0 + z) ** 3 + ok * (1.0 + z) ** 2 + ol) ** 0.5)

        om = self.OMEGA_M
        ol = self.OMEGA_L
        ok = 1.0 - om - ol
        return 4 * pi * (c / self.h) * integrate.romberg(f, z1, z2, (om, ol, ok)) / 1e5

    def rho_crit(self, z):
        h2 = (self.OMEGA_M * (1 + z) ** 3 + self.OMEGA_L) * (self.h / 10.0) ** 2
        return 3 * h2 / (8.0 * pi * G)

    def distance_modulus(self, z):
        return 5 * log10(self.luminosity_distance(z) * 1e5)
