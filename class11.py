import numpy as np
import matplotlib.pyplot as plt
import iapws 
from iapws import IAPWS97 as gas
import math as mp

class CalcExpenditurePar:
    def __init__(self, p_o, t_o, p_nn, t_nn, t_nv, p_k, effect_oi, effect_mex, effect_eg,n):

        self.n = n
        delt_p_o = 0.05 * p_o
        delt_p_nn = 0.1 * p_nn
        delt_p_1 = 0.03 * p_nn
        p_1 = p_nn + delt_p_nn
        p_nv = 1.85 * p_o
        p_to_iapws = 10 ** (-6)
        self.effect_oi = effect_oi
        self.effect_mex = effect_mex
        self.effect_eg = effect_eg 
        self.o_teor = gas(P = p_o * p_to_iapws,T =  t_o + 273.15)
        self.nn_teor = gas(P = (p_nn - delt_p_nn)  * p_to_iapws,T =  t_nn + 273.15)
        self.nn= gas(P = p_nn * p_to_iapws,T =  t_nn + 273.15)
        self.k_teor = gas(P = p_k * p_to_iapws,s = self.nn.s )
        self._1_teor = gas(P = p_1 * p_to_iapws,s = self.o_teor.s)
        self.k_o = gas(P = p_k * p_to_iapws,x = 0 )
        self.nv = gas(P = p_nv * p_to_iapws,T =  t_nv + 273.15)
        self.o_x_o = gas(P = p_o * p_to_iapws,x = 0)
        h_1 = self.o_teor.h - (self.o_teor.h - self._1_teor.h) * self.effect_oi
        h_k = self.nn.h - (self.nn.h - self.k_teor.h) * self.effect_oi
        self.o = gas(P = (p_o - delt_p_o)  * p_to_iapws,h =  self.o_teor.h)
        self.k = gas(P = p_k * p_to_iapws,h = h_k )
        self._1_ = gas(P = (p_1 + delt_p_1 )  * p_to_iapws,h = h_1)
        
    def calc_xi(self):

        chisl = 1 - ((self.k_teor.T) * (self.nn.s 
                                               - self.k_o.s)/((self.o.h
                                                                      - self._1_teor.h) +(self.nn.h - self.k_o.h)))
        znam = 1 - ((self.k_teor.T) * (self.nn.s 
                                              - self.nv.s)/((self.o.h 
                                                                    - self._1_teor.h) +(self.nn.h - self.nv.h)))
        
        xi_dly_opred = 1 - chisl/znam
        x = ((self.nv.T) - self.k_teor.T)/(self.o_x_o.T - self.k_teor.T)
        xi = (0.5469 * (x ** 3)  -  1.9057 * (x ** 2) + 2.2727 * x + 0.0235) * xi_dly_opred
        # Интерполяция производилась по 20-ти точкам, полиномом 3-й степени 
        return xi

    def calc_expenditure(self):
        
        N_to_G = 10 ** (-3)
        self.effect_ip = ((self.o.h - self._1_teor.h) * 
                          self.effect_oi + (self.nn.h - self.k_teor.h) * 
                          self.effect_oi)/((self.o.h - self._1_teor.h) * self.effect_oi + 
                                           (self.nn.h - self.k_o.h))/(1 - self.calc_xi())
        h_1 = self.o.h - (self.o.h - self._1_teor.h) * self.effect_oi 
        h_i = self.effect_ip * ((self.o.h - self.nv.h) + (self.nn.h - h_1)) 
        self.g_o = self.n * N_to_G /( h_i * self.effect_mex * self.effect_eg)
        self.g_k = self.n * N_to_G /( (self.k.h - self.k_o.h ) * self.effect_mex * self.effect_eg) * ((1/self.effect_ip) - 1)
        return self.g_o, self.g_k