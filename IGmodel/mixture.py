import numpy as np
import IGmodel.materials as mat

class Mixture:
    def __init__(self, T, P, N2=0.0, O2=0.0, H2O=0.0, H2=0.0, CO2=0.0, CO=0.0, Ar=0.0) -> None:
        self.T = T
        self.P = P
        self._yi = np.array([N2, O2, H2O, H2, CO2, CO, Ar])
        self._materials = np.array([mat.N2(), mat.O2(), mat.H2O(), mat.H2(), mat.CO2(), mat.CO(), mat.Ar()])

        Mi = []
        for comp in self._materials:
            Mi.append(comp.M)

        self.M = np.dot(self._yi, Mi)

        assert self._yi.sum() == 1.0

    def __repr__(self):
        print_str = "Mixture Composition\n"
        print_str += f"N2:  {self._yi[0]}\n"
        print_str += f"O2:  {self._yi[1]}\n"
        print_str += f"H2O: {self._yi[2]}\n"
        print_str += f"H2:  {self._yi[3]}\n"
        print_str += f"CO2: {self._yi[4]}\n"
        print_str += f"CO:  {self._yi[5]}\n"
        print_str += f"Ar:  {self._yi[6]}\n"
        print_str += "\n"
        print_str += f"Molar Mass: {self.M:.3f} kg/kmol\n"

        return print_str

if __name__ == "__main__":
    mixture = Mixture(T=1300, P=3000, N2=0.79, O2=0.21)
    print(mixture)
    