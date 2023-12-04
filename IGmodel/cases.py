from IGmodel.cycle import CombustionCycle

# Case 0
print("\n---------- Case 0 Results ----------")
T = 20+273.15 # K
P = 100 # kPa
RH = 0.50 # %
omega = 2800 # crank speed, RPM
AFR = 15

case_1 = CombustionCycle(T, P, RH, AFR, omega=omega)
case_1.run_cycle_case_1()

# # Case 1
# print("\n---------- Case 1 Results ----------")
# T = 5+273.15 # K
# P = 102 # kPa
# RH = 0.50 # %
# omega = 1500 # crank speed, RPM
# AFR = 10

# case_1 = CombustionCycle(T, P, RH, AFR, omega=omega)
# case_1.run_cycle_case_1()

# # Case 2
# print("\n---------- Case 2 Results ----------")
# T = 35+273.15 # K
# P = 98 # kPa
# RH = 0.60 # %
# car_speed = 20 # crank speed, RPM
# AFR = 14

# case_2 = CombustionCycle(T, P, RH, AFR, car_speed=0)
# case_2.run_cycle_case_2()
