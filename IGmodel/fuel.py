import numpy as np
from IGmodel.mixture import Ru, Tref

class Fuel:
    def __init__(self, T, P) -> None:
        self.T = T
        self.P = P

        self.xeth = 0.1 # ethanol %
        self.xoct = 0.9 # octane %

        self.LHV = self.xeth*26810 + self.xoct*44430
        self.Cp = self.xeth*2.44 + self.xoct*2.23
        self.hfg = self.xeth*919 + self.xoct*363

        Neth = self.xeth/46.069
        Noct = self.xoct/114.231
        N = Neth + Noct

        self.yeth = Neth/N
        self.yoct = Noct/N

        self.hfo = self.yeth*-235310 + self.yoct*-208450

        self.M = self.yeth*46.069 + self.yoct*114.231
        self.R = Ru/self.M

    @property
    def enthalpy(self):
        return self.hfo/self.M + self.Cp*(self.T-Tref)

    @property
    def internal_energy(self):
        return self.enthalpy - self.R*self.T


