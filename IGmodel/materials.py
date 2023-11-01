import abc
from six import add_metaclass
from scipy.interpolate import interp1d

PROPERTY_TABLE_PATH = "IGmodel/property-tables/"

@add_metaclass(abc.ABCMeta)
class Material:
    def __init__(self, M) -> None:
        self.M = M # kg/kmol

    @staticmethod
    def _interpolate(x, x_known, y_known) -> float:
        f = interp1d(x_known, y_known, kind="linear")
        return f(x)

class IdealGas(Material):
    def __init__(self, M: float, Cp_coefs: list) -> None:
        super().__init__(M=M)
        self._Cp_coefs = Cp_coefs

    def specific_heat(self, T) -> float:
        """
        Ideal gas specific heat [kJ/kmol-K]

        Args:
            - T (float): temperature [K]
        """
        a, b, c, d, e, f, g = self._Cp_coefs
        return a + b*T + c*T**2 + d*T**3 + e*T**4 + f*T**5 + g*T**6

class N2(IdealGas):
    def __init__(self) -> None:
        M = 28.013 # kg/kmol
        Cp_coefs = [3.675,	-1.21E-03,	2.32E-06,	-6.32E-10,	-2.26E-13,	0.00E+00,	0.00E+00]
        super().__init__(M, Cp_coefs)

class O2(IdealGas):
    def __init__(self) -> None:
        M = 31.999 # kg/kmol
        Cp_coefs = [3.626,	-1.88E-03,	7.06E-06,	-6.76E-09,	2.16E-12,	0.00E+00,	0.00E+00]
        super().__init__(M, Cp_coefs)

class H2O(IdealGas):
    def __init__(self) -> None:
        M = 18.0155 # kg/kmol
        Cp_coefs = [4.07,	-1.11E-03,	4.15E-06,	-2.96E-09,	8.07E-13,	0.00E+00,	0.00E+00]
        super().__init__(M, Cp_coefs)

class H2(IdealGas):
    def __init__(self) -> None:
        M = 2.016 # kg/kmol
        Cp_coefs = [3.057,	2.68E-03,	-5.81E-06,	5.52E-09,	-1.81E-12,	0.00E+00,	0.00E+00]
        super().__init__(M, Cp_coefs)

class CO2(IdealGas):
    def __init__(self) -> None:
        M = 44.01  # kg/kmol
        Cp_coefs = [2.401,	8.74E-03,	-6.61E-06,	2.00E-09,	0.00E+00,	0.00E+00,	0.00E+00]
        super().__init__(M, Cp_coefs)

class CO(IdealGas):
    def __init__(self) -> None:
        M = 28.0105 # kg/kmol
        Cp_coefs = [3.71,	-1.62E-03,	3.69E-06,	-2.03E-09,	2.40E-13,	0.00E+00,	0.00E+00]
        super().__init__(M, Cp_coefs)

class Ar(IdealGas):
    def __init__(self) -> None:
        M = 39.948 # kg/kmol
        Cp_coefs = [2.5,	0.00E+00,	0.00E+00,	0.00E+00,	0.00E+00,	0.00E+00,	0.00E+00]
        super().__init__(M, Cp_coefs)
