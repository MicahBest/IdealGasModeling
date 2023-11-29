import numpy as np

class Engine:
    def __init__(self, cylinders, bore, stroke, compression_ratio) -> None:
        self.ncyl = cylinders
        self.bore = bore*0.001      # mm to m
        self.stroke = stroke*0.001  # mm to m
        self.rc = compression_ratio

    @property
    def volume(self):
        return np.pi/4 * self.bore**2 * self.stroke * self.ncyl
