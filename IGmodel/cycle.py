import numpy as np
from IGmodel.air import HumidAir
from IGmodel.engine import Engine
from IGmodel.fuel import Fuel
from IGmodel.mixture import Mixture, Tref
from IGmodel.state import State

class CombustionCycle:
    def __init__(self, T, P, RH, omega, AFR) -> None:
        self.T = T
        self.P = P
        self.RH = RH
        self.omega = omega
        self.AFR = AFR
        self.flamespeed = 50 # m/s
        self.states = []
        self.engine = Engine(cylinders=3, bore=74.0, stroke=77.0, compression_ratio=9.5, intake_diameter=35.0)

    def run_cycle(self):
        # state 0 : ambient air
        air = HumidAir(self.T, self.P, self.RH)
        T0 = air.T
        P0 = air.P
        h0 = air.h

        # state 2 : fuel vapor charge
        fuel = Fuel(T=Tref, P=self.P)
        T2 = fuel.T
        P2 = fuel.P
        h2 = fuel.h
        v2 = fuel.v

        # state 1 : air charge-cooled
        P1 = P0
        h1 = h0 - fuel.hfg/self.AFR
        air.calculate_properties(P=P1, h=h1)
        T1 = air.T
        v1 = air.v

        # throttle pressure loss calcs
        xair = self.AFR/(self.AFR+1)
        xfuel = 1/(1+self.AFR)
        vi = xair*v1 + xfuel*v2
        rho = 1/vi
        Vrms = np.pi/4*self.engine.bore**2 * self.omega
        Vi = (Vrms/2)/(np.pi/4*self.engine.intake_diameter**2)
        kloss = 5
        deltaP = 0.5*rho*Vi**2*kloss

        # state 3 : air charge after throttle
        T3 = T1
        P3 = P1 - deltaP
        air.calculate_properties(T=T3, P=P3)
        v3 = air.v
        u3 = air.u
        s3 = air.s

        # state 4 : fuel charge after throttle
        T4 = T2
        P4 = P2 - deltaP
        fuel.calculate_properties(T=T4, P=P4)
        v4 = fuel.v
        u4 = fuel.u

        # state 5 : dead air mass post combustion
        T5 = 700 # K
        P5 = 105 # kPa
        deadair = Mixture(N2=0.7, CO2=0.15, H2O=0.15)
        deadair.calculate_properties(T=T5, P=P5)
        v5 = deadair.v
        s5 = deadair.s

        # state 6 : dead air mass, expanded
        P6 = P4
        s6 = s5
        deadair.calculate_properties(P=P6, s=s6)
        T6 = deadair.T
        v6 = deadair.v
        u6 = deadair.u
        V6 = self.engine.TDC
        mdead = V6/v6

        # state 7 : dead air mass, shift position
        T7 = T6
        P7 = P6
        v7 = v6
        u7 = u6
        V7 = V6
        s7 = s6

        # state 8 : fuel vapor from state 4
        T8 = T4
        P8 = P4
        v8 = v4
        u8 = u4

        # state 9 : fresh air from state 3
        T9 = T3
        P9 = P3
        v9 = v3
        u9 = u3
        s9 = s3

        # mass calcs
        Vcharge = self.engine.BDC - V7
        mair = Vcharge/v9
        mfuel = mair/self.AFR

        # pump work
        W_intake = Vcharge * deltaP * 1000 # kJ to J

        # state 10 : dead air after compression
        s10 = s7
        v10 = v7/self.engine.rc
        deadair.calculate_properties(s=s10, v=v10)
        T10 = deadair.T
        P10 = deadair.P
        u10 = deadair.u

        N_H2O_10 = deadair._yi[2]
        N_H2_10 = deadair._yi[3]
        N_CO2_10 = deadair._yi[4]

        # state 11 : fuel vapor after compression
        T11 = T8 * self.engine.rc**(fuel.k - 1)
        v11 = v8/self.engine.rc
        P11 = fuel.R*T11/v11
        fuel.calculate_properties(T=T11, P=P11)
        P11 = fuel.P
        h11 = fuel.h
        u11 = fuel.u

        # state 12 : fresh air after compression
        s12 = s9
        v12 = v9/self.engine.rc
        air.calculate_properties(s=s12, v=v12)
        T12 = air.T
        P12 = air.P
        u12 = air.u

        N_N2_12 = air._yi[0]
        N_O2_12 = air._yi[1]
        N_H2O_12 = air._yi[2]
        N_Ar_12 = air._yi[6]

        # state 13 : combustion products, end of isometric

        # state 14 : combustion products, end of isobaric
        delta_t = self.engine.bore/(2*self.flamespeed)
        print(f"Combustion time: {delta_t}")
        print(f"Crank pos. at cutoff: {self.omega*360/60*delta_t}")
        delta = (self.engine.stroke/2) * (1 - np.sin(self.omega*delta_t + np.pi/2))
        V14 = self.engine.TDC + delta*np.pi/4*self.engine.bore**2
        v14 = V14 / (mair + mfuel + mdead)
        print(f"Cutoff ratio: {V14/self.engine.TDC}")
        print(f"Piston RMS: {Vrms}")

        # heat loss to cylinder walls
        h = 50 # W/m2-K
        Ts = 100+273.15 # K
        As = (V14/(np.pi/4*self.engine.bore**2)*np.pi*self.engine.bore + np.pi/4*self.engine.bore**2)
        Ti = (mdead*T10 + mfuel*T11 + mair*T12)/(mdead + mfuel + mair)
        Te = mfuel*fuel.LHV/((mdead + mfuel + mair)*fuel.Cv) + Ti
        Tinf = (Te + Ti)/2
        Qdot_out = h*As*(Ts-Tinf)
        Q_out = Qdot_out*delta_t

        # combustion mass balance
        Ndead = mdead/deadair.M
        Nfuel = mfuel/fuel.M
        Nair = mair/air.M
        print(f"Working charge mass: {mair+mfuel}")

        Neth = Nfuel*fuel.yeth
        Noct = Nfuel*fuel.yoct

        N_O2_theoretical = 2.5*Neth + 12.5*Noct
        phi = N_O2_12/N_O2_theoretical
        O2_def = N_O2_theoretical*(1 - phi)

        N_CO = np.max([4/3*O2_def, 0.1*(2*Neth + 8*Noct)])
        N_H2 = np.max([4/3*O2_def, 0.1*(2*Neth + 8*Noct)])

        N_CO2_add = 2*Neth + 8*Noct - N_CO
        N_H2O_add = 2*Neth + 9*Noct - N_H2
        N_O2_add = N_O2_12 - N_CO2_add - N_CO/2 - N_H2O_add/2 + Neth/2

        products = Mixture(
            N2=N_H2_10+N_N2_12,
            O2=N_O2_add,
            H2O=N_H2O_10+N_H2O_12+N_H2O_add,
            H2=N_H2,
            CO2=N_CO2_10+N_CO2_add,
            CO=N_CO,
            Ar=N_Ar_12
        )

        # state 15 : combustion products, end of power

        # state 16 : combustion products, expanded
