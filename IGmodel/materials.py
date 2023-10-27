import abc
import pandas as pd
from six import add_metaclass
from scipy.interpolate import interp1d

PROPERTY_TABLE_PATH = "IGmodel/property-tables/"

@add_metaclass(abc.ABCMeta)
class Material:
    def __init__(self) -> None:
        pass

    def _interpolate(self, x, x_known, y_known) -> float:
        f = interp1d(x_known, y_known, kind="linear")
        return f(x)

class H2O(Material):
    def __init__(self, phase: str="liquid") -> None:
        super().__init__()
        self._phase = phase
        self._M = 18.015 # kg/kmol
        self._R = 0.4615 # kJ/kg-K
        self._property_table = "H2O-sat-temperature.csv"
        df = pd.read_csv(PROPERTY_TABLE_PATH+self._property_table)

        self._T = df["T"].to_numpy()
        self._Psat = df["Psat"].to_numpy()
        self._ufg = df["ufg"].to_numpy()
        self._hfg = df["hfg"].to_numpy()
        self._sfg = df["sfg"].to_numpy()
        if self._phase == "liquid":
            self._v = df["vf"].to_numpy()
            self._h = df["hf"].to_numpy()
            self._u = df["uf"].to_numpy()
            self._s = df["sf"].to_numpy()
        elif self._phase in ["vapor", "gas"]:
            self._v = df["vg"].to_numpy()
            self._h = df["hg"].to_numpy()
            self._u = df["ug"].to_numpy()
            self._s = df["sg"].to_numpy()

    def enthalpy(self, T):
        return self._interpolate(T, self._T, self._h)

class IdealGas(Material):
    def __init__(self, property_file: str) -> None:
        super().__init__()
        df = pd.read_csv(PROPERTY_TABLE_PATH+property_file)
        self._T = df["T"].to_numpy()
        self._h_bar = df["h-bar"].to_numpy()
        self._u_bar = df["u-bar"].to_numpy()
        self._s_bar = df["s-bar"].to_numpy()

    def specific_heat(self, T) -> float:
        """
        Ideal gas specific heat [kJ/kmol-K]

        Args:
            - T (float): temperature [K]
        """
        a, b, c, d = self._Cp_coefs
        return a + b*T + c*T*T + d*T*T*T

class O2(IdealGas):
    def __init__(self) -> None:
        self._property_table = "O2-ideal-gas.csv"
        self._M = 31.999 # kg/kmol
        self._R = 0.2598 # kJ/kg-K
        self._Cp_coefs = [25.48, 1.520e-2, -0.7155e-5, 1.312e-9]
        super().__init__(self._property_table)

class N2(IdealGas):
    def __init__(self) -> None:
        self._property_table = "N2-ideal-gas.csv"
        self._M = 28.013 # kg/kmol
        self._R = 0.2968 # kJ/kg-K
        self._Cp_coefs = [28.90, -0.1571e-2, 0.8081e-5, -2.873e-9]
        super().__init__(self._property_table)

class CO2(IdealGas):
    def __init__(self) -> None:
        self._property_table = "CO2-ideal-gas.csv"
        self._M = 44.01  # kg/kmol
        self._R = 0.1889 # kJ/kg-K
        self._Cp_coefs = [22.26, 5.981e-2, -3.501e-5, 7.469e-9]
        super().__init__(self._property_table)

class CO(IdealGas):
    def __init__(self) -> None:
        self._property_table = "CO-vapor-ideal-gas.csv"
        self._M = 28.011 # kg/kmol
        self._R = 0.2968 # kJ/kg-K
        self._Cp_coefs = [28.16, 0.1675e-2, 0.5375e-5, -2.222e-9]
        super().__init__(self._property_table)

class H2O_vapor_IG(IdealGas):
    def __init__(self) -> None:
        self._property_table = "H2O-vapor-ideal-gas.csv"
        self._M = 18.015 # kg/kmol
        self._R = 0.4615 # kJ/kg-K
        self._Cp_coefs = [32.24, 0.1923e-2, 1.055e-5, -3.595e-9]
        super().__init__(self._property_table)

if __name__ == "__main__":
    print(H2O(phase="vapor").enthalpy(25))
