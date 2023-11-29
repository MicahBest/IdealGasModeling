from IGmodel.cycle import CombustionCycle

# Case 1
T = 5+273.15 # K
P = 102 # kPa
RH = 0.50 # %
omega = 1500 # crank speed, RPM
AFR = 10

case_1 = CombustionCycle(T, P, RH, omega, AFR)
case_1.run_cycle()
