import abc
import pandas as pd
from six import add_metaclass

PROPERTY_TABLE_PATH = "IGmodel/property-tables/"

@add_metaclass(abc.ABCMeta)
class Material:
    def __init__(self) -> None:
        pass

    def _interpolate(self, x, known_y, known_x) -> float:
        pass

class H2O(Material):
    def __init__(self) -> None:
        super().__init__()
        self._property_table = "H2O-sat-temperature.csv"
        df = pd.read_csv(PROPERTY_TABLE_PATH+self._property_table)
        self._T = df["T"].to_numpy()
        self._Psat = df["Psat"].to_numpy()
        self._vf = df["vf"].to_numpy()
        self._vg = df["vg"].to_numpy()
        self._uf = df["uf"].to_numpy()
        self._ufg = df["ufg"].to_numpy()
        self._ug = df["ug"].to_numpy()
        self._hf = df["hf"].to_numpy()
        self._hfg = df["hfg"].to_numpy()
        self._hg = df["hg"].to_numpy()
        self._sf = df["sf"].to_numpy()
        self._sfg = df["sfg"].to_numpy()
        self._sg = df["sg"].to_numpy()

class IdealGas(Material):
    def __init__(self, property_file) -> None:
        super().__init__()
        df = pd.read_csv(PROPERTY_TABLE_PATH+property_file)
        self._T = df["T"].to_numpy()
        self._h_bar = df["h-bar"].to_numpy()
        self._u_bar = df["u-bar"].to_numpy()
        self._s_bar = df["s-bar"].to_numpy()

class Air(IdealGas):
    def __init__(self) -> None:
        self._property_table = "Air-ideal-gas.csv"
        self._M = 28.97  # kg/kmol
        self._R = 0.2870 # kJ/kg-K
        super().__init__(self._property_table)

class O2(IdealGas):
    def __init__(self) -> None:
        self._property_table = "O2-ideal-gas.csv"
        self._M = 31.999 # kg/kmol
        self._R = 0.2598 # kJ/kg-K
        super().__init__(self._property_table)

class N2(IdealGas):
    def __init__(self) -> None:
        self._property_table = "N2-ideal-gas.csv"
        self._M = 28.013 # kg/kmol
        self._R = 0.2968 # kJ/kg-K
        super().__init__(self._property_table)

if __name__ == "__main__":
    H2O()
