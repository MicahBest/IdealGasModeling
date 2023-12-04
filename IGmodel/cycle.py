import numpy as np
from IGmodel.air import HumidAir
from IGmodel.engine import Engine
from IGmodel.fuel import Fuel
from IGmodel.mixture import Mixture, Tref
from IGmodel.state import State

class CombustionCycle:
    def __init__(self, T, P, RH, AFR, omega=None, car_speed=None) -> None:
        self.T = T
        self.P = P
        self.RH = RH
        self.omega_rpm = omega
        self.carspeed = car_speed
        self.AFR = AFR
        self.flamespeed = 50 # m/s
        self.states = []
        self.engine = Engine(cylinders=3, bore=74.0, stroke=77.0, compression_ratio=9.5, intake_diameter=35.0, exhaust_diameter=28.0)

    def run_cycle_case_1(self):
        self.omega_deg = self.omega_rpm*360/60 # RPM to deg/s
        self.omega_rad = self.omega_rpm*2*np.pi/60 # RPM to rad/s

        print(f"Crankspeed: {self.omega_rad:.3f} rad/s")

        # state 0 : ambient air
        air = HumidAir(self.T, self.P, self.RH)
        T0 = air.T
        P0 = air.P
        h0 = air.h
        v0 = air.v

        # state 2 : fuel vapor charge
        fuel = Fuel(T=Tref, P=self.P)
        T2 = fuel.T
        P2 = fuel.P
        h2 = fuel.h
        v2 = fuel.v
        print(f"Spec. Vol. Fuel: {v2:.3f}")

        # state 1 : air charge-cooled
        P1 = P0
        h1 = h0 - fuel.hfg/self.AFR
        air.calculate_properties(P=P1, h=h1)
        T1 = air.T
        v1 = air.v
        print(f"Spec. Vol. Air: {v1:.3f}")

        # throttle pressure loss calcs
        xair = self.AFR/(self.AFR+1)
        xfuel = 1/(1+self.AFR)
        vi = xair*v1 + xfuel*v2
        print(f"Spec. Vol. Intake: {vi:.6f}")
        rho = 1/vi
        print(f"Density Intake Air: {rho:.3f}")
        Vrms = self.engine.stroke*self.omega_rad/(2*np.sqrt(2))
        print(f"V RMS: {Vrms:.3f}")
        Vdot_rms = np.pi/4*self.engine.bore**2 * Vrms
        print(f"Vdot RMS: {Vdot_rms:.3f}")
        Vi = (Vdot_rms/2)/(np.pi/4*self.engine.intake_diameter**2)
        print(f"Intake Velocity: {Vi:.3f}")
        kloss = 5
        deltaP_in = 0.5*rho*Vi**2*kloss/1000 # kPa
        print(f"Delta P Intake: {deltaP_in:.3f} kPa")

        # state 3 : air charge after throttle
        T3 = T1
        P3 = P1 - deltaP_in
        air.calculate_properties(T=T3, P=P3)
        v3 = air.v
        u3 = air.u
        s3 = air.s

        # state 4 : fuel charge after throttle
        T4 = T2
        P4 = P2 - deltaP_in
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
        T11 = fuel.T
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

        # state 14 : combustion products, end of isobaric
        delta_t = self.engine.bore/(2*self.flamespeed)
        delta = (self.engine.stroke/2) * (1 - np.sin(self.omega_rad*delta_t + np.pi/2))
        V14 = self.engine.TDC + delta*np.pi/4*self.engine.bore**2
        v14 = V14 / (mair + mfuel + mdead)
        rcutoff = V14/self.engine.TDC

        # heat loss to cylinder walls
        h = 50 # W/m2-K
        Ts = 100+273.15 # K
        As = (V14/(np.pi/4*self.engine.bore**2)*np.pi*self.engine.bore + np.pi/4*self.engine.bore**2)
        Ti = (mdead*T10 + mfuel*T11 + mair*T12)/(mdead + mfuel + mair)
        Te = mfuel*fuel.LHV/((mdead + mfuel + mair)*fuel.Cv) + Ti
        Tinf = (Te + Ti)/2
        print(f"Tinf: {Tinf:.3f}")
        Qdot_out = -h*As*(Ts-Tinf)
        print(f"Convection Q-dot-out: {Qdot_out:.3f}")
        Q_out = Qdot_out*delta_t
        print(f"Convection Q-out: {Q_out:.3f} J")

        # combustion mass balance
        Ndead = mdead/deadair.M
        Nfuel = mfuel/fuel.M
        Nair = mair/air.M

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
        m14 = mair + mfuel + mdead

        u14 = (mair*u12 + mfuel*u11 + mdead*u10 - Q_out)/m14
        products.calculate_properties(v=v14, u=u14)
        T14 = products.T
        P14 = products.P
        s14 = products.s

        # state 13 : combustion products, end of isometric
        T13 = T14 / rcutoff
        v13 = v14 / rcutoff
        P13 = P14
        print(f"T13: {T13:.3f}")
        print(f"T14: {T14:.3f}")

        # state 15 : combustion products, end of power
        v15 = v14 * self.engine.BDC / V14
        s15 = s14
        products.calculate_properties(v=v15, s=s15)
        T15 = products.T
        P15 = products.P
        u15 = products.u

        # state 16 : combustion products, expanded
        rho16 = 1/v15
        Ve = (Vdot_rms/2)/(np.pi/4*self.engine.exhaust_diameter**2)
        deltaP_ex = 0.5*rho16*Ve**2*kloss/1000
        s16 = s15
        P16 = P0 + deltaP_ex
        products.calculate_properties(s=s16, P=P16)
        T16 = products.T
        v16 = products.v
        u16 = products.u
        h16 = products.h
        V16 = m14*v16
        Vexpel = self.engine.TDC - V16

        # computing Wnet,cycle
        W_intake  = deltaP_in * (self.engine.BDC - V6) # kJ
        W_exhaust = deltaP_ex * (Vexpel) # kJ
        if W_exhaust < 0:
            W_exhaust = 0
        viscosity = 0.01 # N-s/m^2
        tau = viscosity*Vrms/0.22e-3
        A = np.pi*self.engine.bore*0.01e-3
        W_stroke  = tau*A*self.engine.stroke
        print(f"Mass Products: {m14:.6f} kg")
        Wout_comb = m14*products.R*(T13-T14)
        Wout_power = m14*(u14-u15)
        Win_comp = mdead*(u10-u7) + mfuel*(u11-u8) + mair*(u12-u9)
        Wnet = Wout_comb + Wout_power - Win_comp - W_intake - W_exhaust - 4*W_stroke

        Qdot_in = mfuel*fuel.LHV
        eta = Wnet/Qdot_in
        MEP = Wnet/(self.engine.BDC-self.engine.TDC)
        Power = self.engine.ncyl * Wnet * self.omega_rpm/120
        mdot_fuel = self.engine.ncyl * mfuel * self.omega_rpm/120 # kg/s
        Vdot_fuel = mdot_fuel / fuel.density * 951019.3885 # m^3/s to gal/hr

        print()
        print(f"Combustion time: {delta_t*1000:.3f} ms")
        print(f"Crank pos. at cutoff: {self.omega_deg*delta_t:.3f} deg")
        print(f"Cutoff ratio: {rcutoff:.3f}")
        print(f"Piston RMS: {Vrms:.3f} m/s")
        print(f"Working charge mass: {(mair+mfuel)*1000:.3f} g")
        print(f"Heat Loss: {Q_out:.3f} kJ")
        print(f"Stroke Work: {W_stroke:.3f} kJ")
        print(f"Intake Pump Work: {W_intake:.3f} kJ")
        print(f"Exhaust Pump Work: {W_exhaust:.3f} kJ")
        print(f"Compression Work: {Win_comp:.3f} kJ")
        print(f"Combustion Work: {Wout_comb:.3f} kJ")
        print(f"Power Stroke Work: {Wout_power:.3f} kJ")
        print(f"Net Work Cycle: {Wnet:.3f} kJ")
        print(f"Fuel Rate: {Vdot_fuel/3600:.6f} g/s")
        print(f"Power: {Power:.3f} kW")
        print(f"Power: {Power*1.341022:.3f} hp")
        print(f"Cycle Efficiency: {eta:.3f}")
        print(f"MEP: {MEP:.3f} kPa")

    def run_cycle_case_2(self):

        self.omega_rpm = 1500
        update_rate = 100

        nIter = 0
        maxIter = 100
        while nIter < maxIter:
            self.omega_deg = self.omega_rpm*360/60 # RPM to deg/s
            self.omega_rad = self.omega_rpm*2*np.pi/60 # RPM to rad/s

            nIter += 1
            # state 0 : ambient air
            air = HumidAir(self.T, self.P, self.RH)
            T0 = air.T
            P0 = air.P
            h0 = air.h
            v0 = air.v

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
            Vrms = self.engine.stroke*self.omega_rad/(2*np.sqrt(2))
            Vdot_rms = np.pi/4*self.engine.bore**2 * Vrms
            Vi = (Vdot_rms/2)/(np.pi/4*self.engine.intake_diameter**2)
            kloss = 5
            deltaP_in = 0.5*rho*Vi**2*kloss/1000 # kPa

            # state 3 : air charge after throttle
            T3 = T1
            P3 = P1 - deltaP_in
            air.calculate_properties(T=T3, P=P3)
            v3 = air.v
            u3 = air.u
            s3 = air.s

            # state 4 : fuel charge after throttle
            T4 = T2
            P4 = P2 - deltaP_in
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
            T11 = fuel.T
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

            # state 14 : combustion products, end of isobaric
            delta_t = self.engine.bore/(2*self.flamespeed)
            delta = (self.engine.stroke/2) * (1 - np.sin(self.omega_rad*delta_t + np.pi/2))
            V14 = self.engine.TDC + delta*np.pi/4*self.engine.bore**2
            v14 = V14 / (mair + mfuel + mdead)
            rcutoff = V14/self.engine.TDC

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
            m14 = mair + mfuel + mdead

            u14 = (mair*u12 + mfuel*u11 + mdead*u10 - Q_out)/m14
            products.calculate_properties(v=v14, u=u14)
            T14 = products.T
            P14 = products.P
            s14 = products.s

            # state 13 : combustion products, end of isometric
            T13 = T14 / rcutoff
            v13 = v14 / rcutoff
            P13 = P14

            # state 15 : combustion products, end of power
            v15 = v14 * self.engine.BDC / V14
            s15 = s14
            products.calculate_properties(v=v15, s=s15)
            T15 = products.T
            P15 = products.P
            u15 = products.u

            # state 16 : combustion products, expanded
            rho16 = 1/v15
            Ve = (Vdot_rms/2)/(np.pi/4*self.engine.exhaust_diameter**2)
            deltaP_ex = 0.5*rho16*Ve**2*kloss/1000
            s16 = s15
            P16 = P0 + deltaP_ex
            products.calculate_properties(s=s16, P=P16)
            T16 = products.T
            v16 = products.v
            u16 = products.u
            h16 = products.h
            V16 = m14*v16
            Vexpel = self.engine.TDC - V16

            # computing Wnet,cycle
            W_intake  = deltaP_in * (self.engine.BDC - V6) # kJ
            W_exhaust = deltaP_ex * (Vexpel) # kJ
            if W_exhaust < 0:
                W_exhaust = 0
            viscosity = 0.01 # N-s/m^2
            tau = viscosity*Vrms/0.22e-3
            A = np.pi*self.engine.bore*0.01e-3
            W_stroke  = tau*A*self.engine.stroke
            Wout_comb = m14*products.R*(T13-T14)
            Wout_power = m14*(u14-u15)
            Win_comp = mdead*(u10-u7) + mfuel*(u11-u8) + mair*(u12-u9)
            Wnet = Wout_comb + Wout_power - Win_comp - W_intake - W_exhaust - 4*W_stroke

            Qdot_in = mfuel*fuel.LHV
            eta = Wnet/Qdot_in
            MEP = Wnet/(self.engine.BDC-self.engine.TDC)
            Power = self.engine.ncyl * Wnet * self.omega_rpm/120
            mdot_fuel = self.engine.ncyl * mfuel * self.omega_rpm/120 # kg/s
            Vdot_fuel = mdot_fuel / fuel.density * 951019.3885 # m^3/s to gal/hr

            CD = 0.34
            Af = 1.86 # m^2
            mcar = 1000 # kg
            g = 9.81 # m/s^2
            Cr = 0.015
            Vcar = self.carspeed # mph
            eta_transmission = 0.85
            FE = np.abs(Vcar / Vdot_fuel) # mpg
            Vcar *= 0.44704 # mph to m/s
            Drag = 0.5/v0 * Vcar**2 * CD * Af
            RollingResistance = mcar*g * Cr
            Wdot_mech = Vcar * (Drag + RollingResistance)
            Wdot_crank = Wdot_mech/eta_transmission

            self.omega_rpm += update_rate*(Wdot_crank-Power)

        print(f"Speed: {self.carspeed:.3f} mph")
        print(f"Speed: {self.carspeed*0.44704:.3f} m/s")
        print(f"Aerodynamic Drag: {Drag:.3f} N")
        print(f"Rolling Resistance: {RollingResistance:.3f} N")
        print(f"Mechanical Power Req.: {Wdot_mech:.3f} kW")
        print(f"Crank Power Req.: {Wdot_crank:.3f} kW")
        print(f"Fuel Economy: {FE:.3f} mpg\n")

        print(f"Crankspeed: {self.omega_rpm:.3f} RPM")
        print(f"Combustion time: {delta_t*1000:.3f} ms")
        print(f"Crank pos. at cutoff: {self.omega_deg*delta_t:.3f} deg")
        print(f"Cutoff ratio: {rcutoff:.3f}")
        print(f"Piston RMS: {Vrms:.3f} m/s")
        print(f"Working charge mass: {(mair+mfuel)*1000:.3f} g")
        print(f"Heat Loss: {Q_out:.3f} kJ")
        print(f"Stroke Work: {W_stroke:.3f} kJ")
        print(f"Intake Pump Work: {W_intake:.3f} kJ")
        print(f"Exhaust Pump Work: {W_exhaust:.3f} kJ")
        print(f"Compression Work: {Win_comp:.3f} kJ")
        print(f"Combustion Work: {Wout_comb:.3f} kJ")
        print(f"Power Stroke Work: {Wout_power:.3f} kJ")
        print(f"Net Work Cycle: {Wnet:.3f} kJ")
        print(f"Fuel Rate: {Vdot_fuel/3600:.6f} g/s")
        print(f"Power: {Power:.3f} kW")
        print(f"Power: {Power*1.341022:.3f} hp")
        print(f"Cycle Efficiency: {eta:.3f}")
        print(f"MEP: {MEP:.3f} kPa")