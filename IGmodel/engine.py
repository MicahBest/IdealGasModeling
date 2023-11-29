import numpy as np

class Engine:
    def __init__(self, cylinders, bore, stroke, compression_ratio, intake_diameter) -> None:
        self.ncyl = cylinders
        self.bore = bore*0.001      # mm to m
        self.stroke = stroke*0.001  # mm to m
        self.rc = compression_ratio
        self.intake_diameter = intake_diameter
        self.TDC = self.cyl_volume/(self.rc - 1)
        self.BDC = self.TDC*self.rc

    @property
    def cyl_volume(self):
        return np.pi/4 * self.bore**2 * self.stroke

    @property
    def displaced_volume(self):
        return self.cyl_volume * self.ncyl
