import numpy as np
from IGmodel.mixture import Mixture, Ru

class DryAir(Mixture):
    def __init__(self) -> None:
        super().__init__(N2=0.78, O2=0.21, Ar=0.01)

class HumidAir(Mixture):
    def __init__(self, T, P, RH=1.00) -> None:
        """
        Args:
            - T (float)  : temperature          [K]
            - P (float)  : pressure             [kPa]
            - RH (float) : relative humidity    [%]
        """
        PH2O = RH * 0.61094 * np.exp((17.625*(T-273.15))/(T-273.15+243.04)) # using Magnus Eq. with temperature in [C]
        NH2O = PH2O/(Ru*T)
        Nair = (P-PH2O)/(Ru*T)

        N2 = 0.78 * Nair/(NH2O + Nair)
        O2 = 0.21 * Nair/(NH2O + Nair)
        Ar = 0.01 * Nair/(NH2O + Nair)
        H2O = NH2O/(NH2O + Nair)

        super().__init__(N2=N2, O2=O2, Ar=Ar, H2O=H2O)
        self.RH = RH
        self.calculate_properties(T=T, P=P)
