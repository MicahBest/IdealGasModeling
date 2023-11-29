import numpy as np
from IGmodel.air import DryAir, HumidAir

def test_dry_air():
    air = DryAir()

def test_humid_air():
    air = HumidAir(T=300, P=100, RH=0.7)
    np.testing.assert_almost_equal(air.RH, 0.700000, decimal=6)
