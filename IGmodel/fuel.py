import numpy as np
from IGmodel.mixture import Ru, Tref

class Fuel:
    def __init__(self, T, P) -> None:
        self.T = T
        self.P = P
        self._yi = [0, 0, 0, 0, 0, 0, 0]

        self.xeth = 0.1 # ethanol %
        self.xoct = 0.9 # octane %

        self.LHV = self.xeth*26810 + self.xoct*44430
        self.Cp = self.xeth*2.44 + self.xoct*2.23
        self.hfg = self.xeth*919 + self.xoct*363
        self.density = self.xeth*790 + self.xoct*703

        Neth = self.xeth/44.0535
        Noct = self.xoct/114.232
        N = Neth + Noct

        self.yeth = Neth/N
        self.yoct = Noct/N

        self.hfo = self.yeth*-235310 + self.yoct*-208450

        self.M = self.yeth*44.0535 + self.yoct*114.232
        self.R = Ru/self.M

        self.v = self.R*self.T/self.P
        self.h = self.hfo/self.M + self.Cp*(self.T-Tref)
        self.u = self.h - self.R*self.T

        self.Cv = self.Cp - self.R
        self.k = self.Cp/(self.Cp-self.R)

    def calculate_properties(self, T, P):
        self.__init__(T, P)
