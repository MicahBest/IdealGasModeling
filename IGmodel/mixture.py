import numpy as np
import pandas as pd

Ru = 8.31447 # kJ/kmol-K

class Mixture:
    def __init__(self, N2=0.0, O2=0.0, H2O=0.0, H2=0.0, CO2=0.0, CO=0.0, Ar=0.0) -> None:
        self._yi = np.array([N2, O2, H2O, H2, CO2, CO, Ar])
        assert self._yi.sum() == 1.0

        df = pd.read_csv("IGmodel/CpModelTable.csv")

        # Molar mass
        _Mi = df["M"]
        self.M = np.dot(self._yi, _Mi)      # kg/kmol
        self.R = Ru/self.M                  # kJ/kg-K

        # enthalpy of formation
        _hfoi = df["hfo"]
        self.hfo = np.dot(self._yi, _hfoi)  # kJ/kmol

        # specific heat model coefficients
        df = df.drop(["Gas", "M", "hfo"], axis=1)
        self.Cp_coefs = np.zeros(7)
        for i, coef in enumerate(df.columns):
            self.Cp_coefs[i] = np.dot(self._yi, df[coef])

        # initialized as 0.0 until calculation is performed
        self.Cp = 0.0   # kJ/kg-K
        self.v = 0.0    # m^3/kg
        self.h = 0.0    # kJ/kg
        self.u = 0.0    # kJ/kg
        self.s = 0.0    # kJ/kg-K

    def calculate_properties_explicit(self, T, P):
        """Explicitly calculates properties from known temperature and pressure."""
        # specific heat calculation
        order = np.array([0, 1, 2, 3, 4, 5, 6])
        T_vec = np.array([T**i for i in order])
        self.Cp = np.dot(self.Cp_coefs, T_vec)*self.R

        # specific volume calculation
        self.v = self.R*T/P

        # specific enthalpy calculation
        FT = np.array([T**(i+1)/(i+1) for i in order])
        hT = np.dot(self.Cp_coefs, FT)*self.R
        self.h = hT + self.hfo/self.M

        # specific internal energy calculation
        self.u = self.h - P*self.v

        # specific entropy calculation
        GT = np.array([np.log(T)] + list(FT[:-1]))
        self.s = self.R*np.dot(self.Cp_coefs, GT) - self.R*np.log(P)


    def __repr__(self):
        """Prints mixture info."""
        print_str = "Mixture Composition\n"
        print_str += f"N2:  {self._yi[0]*100:>10.5f} %\n"
        print_str += f"O2:  {self._yi[1]*100:>10.5f} %\n"
        print_str += f"H2O: {self._yi[2]*100:>10.5f} %\n"
        print_str += f"H2:  {self._yi[3]*100:>10.5f} %\n"
        print_str += f"CO2: {self._yi[4]*100:>10.5f} %\n"
        print_str += f"CO:  {self._yi[5]*100:>10.5f} %\n"
        print_str += f"Ar:  {self._yi[6]*100:>10.5f} %\n"
        print_str += "\n"

        print_str += f"Molar Mass (M)               : {self.M:12f} kg/kmol  \n"
        print_str += f"Gas constant (R)             : {self.R:12f} kJ/kg-K  \n"
        print_str += f"Specific heat (Cp)           : {self.Cp:12f} kJ/kg-K \n"
        print_str += f"Specific volume (v)          : {self.v:12f} m^3/kg   \n"
        print_str += f"Specific enthalpy (h)        : {self.h:12f} kJ/kg    \n"
        print_str += f"Specific internal energy (u) : {self.u:12f} kJ/kg    \n"
        print_str += f"Specific entropy (s)         : {self.s:12f} kJ/kg-K  \n"

        return print_str

if __name__ == "__main__":
    mixture = Mixture(N2=0.79, O2=0.21)
    mixture.calculate_properties_explicit(T=1300, P=3000)
    print(mixture)
