import numpy as np
import matplotlib.pyplot as plt
import iapws 
from iapws import IAPWS97 as gas

MPa = 10 ** 6
kPa = 10 ** 3
MW  = 10 ** 6
unit = 1 / MPa
to_kelvin = lambda x: x + 273.15 if x else None



class Params:

    def __init__(self, N, P0, t0, P_pp, t_pp, t_fw, Pk, internal_efficiency, mechanical_efficiency, generator_efficiency):	

        self.N = N
        self.internal_efficiency = internal_efficiency
        self.mechanical_efficiency = mechanical_efficiency
        self.generator_efficiency = generator_efficiency
        T_0 = to_kelvin(t0)
        T_pp = to_kelvin(t_pp)
        T_fw = to_kelvin(t_fw)
    
        P_fw = 1.33 * P0 
        delta_p0 = 0.05 * P0
        real_p0 = P0 - delta_p0
        delta_p1 = 0.1 * P_pp
        real_p1 = P_pp + delta_p1
        delta_p_pp = 0.02 * real_p1
        real_p_pp = P_pp - delta_p_pp
    
        # Нахождение параметров
    
        self._point_0 = gas(P = P0 * unit, T=T_0) # Теоретическая точка 0
        self.point_0_true = gas(P = real_p0 * unit, h=self._point_0.h) # Действительня точка 0 с учетом потерь клапанов ЦВД
        self._point_1 = gas (P = real_p1 * unit, s=self._point_0.s) # Теоретическая точка 1 
        self.hp_heat_drop = (self._point_0.h - self._point_1.h) * self.internal_efficiency # Действительный теплоперепад
        self.h_1_true = self._point_0.h - self.hp_heat_drop # Нахождение энтальпии дейстительной  точки 1  
        self.point_1_true = gas(P=real_p1 * unit, h=self.h_1_true) # Нахождение действительной точки 1
        self._point_pp = gas(P=P_pp * unit, T=T_pp) # Нахождение точки ПП 
        self.point_pp_true = gas (P = real_p_pp * unit , h=self._point_pp.h) # Нахождение действительной  точки ПП 
        self._point_k = gas(P=Pk * unit, s=self._point_pp.s) # Нахождение точки конденсатора 
        self.lp_heat_drop = (self._point_pp.h - self._point_k.h) * self.internal_efficiency # Нахождение Теплоперепада после ПП 
        self.h_k_true = self._point_pp.h - self.lp_heat_drop # Нахождение действительной точки конденсатора 
        self.point_k_true = gas(P=Pk * unit, h=self.h_k_true)
        self.point_k_water = gas(P=Pk * unit, x=0) # Нахождение точки воды в состоянии насыщении
        self.point_feed_water = gas(P=P_fw * unit, T=T_fw) # Нахождение точки питательной воды
        
        self.numenator_without = self.point_k_true.T * (self._point_pp.s - self.point_k_water.s)
        self.denumenator_without = (self.point_0_true.h - self._point_1.h) + (self.point_pp_true.h - self.point_k_water.h)
        self.without_part = 1 - (self.numenator_without / self.denumenator_without)

        self.numenator_infinity = self.point_k_true.T * (self._point_pp.s - self.point_feed_water.s)
        self.denumenator_infinity = (self.point_0_true.h - self._point_1.h) + (self.point_pp_true.h - self.point_feed_water.h)
        self.infinity_part = 1 - (self.numenator_infinity / self.denumenator_infinity)
        self.ksi_infinity = 1 - (self.without_part / self.infinity_part)
                                
        self.coeff = (self.point_feed_water.T - self.point_k_true.T) / (to_kelvin(374.2) - self.point_k_true.T) # Коэфф по оси х

        

    def print_coeff (self):
        print("Коэффициент отношений разности температур = ", self.coeff)

    def get_input (self,ksi_input):
        
        ksi = ksi_input * self.ksi_infinity
        return ksi

    def get_mas_flow (self,ksi):
        self.eff_num = self.hp_heat_drop + self.lp_heat_drop
        self.eff_denum = self.hp_heat_drop + (self.point_pp_true.h - self.point_k_water.h)
        self.efficiency = (self.eff_num / self.eff_denum) * (1 / (1 - ksi))
        self.estimated_heat_drop = self.efficiency * ((self.point_0_true.h - self.point_feed_water.h) + (self.point_pp_true.h - self.point_1_true.h))
        self.G0 = (self.N * MW) / (self.estimated_heat_drop * 1000 * self.mechanical_efficiency * self.generator_efficiency)
        self.Gk = (
        (self.N * MW) /
        ((self.point_k_true.h - self.point_k_water.h) * 1000 * self.mechanical_efficiency * self.generator_efficiency) * ((1 / self.efficiency) - 1))
        return self.G0, self.Gk