import numpy as np
import iapws 
from iapws import IAPWS97 as par
import pandas as pd
import matplotlib.pyplot as plt
class RegularTurbin:
    
    def __init__(self, tochka_o, n, H_o, ro, G, d, k):
        
        self.k = 1.33
        self.n = n
        self.ro = ro
        self.G = G
        self.H_o = H_o
        self.tochka_o = tochka_o
        self.k = k
        self.d = d
        self.H_oc = (1 - self.ro) * self.H_o 
        self.H_or = self.H_o * self.ro
        self.h_1_t = self.tochka_o.h - self.H_oc 
        self.tochka_1_t = par(h = self.h_1_t, s = self.tochka_o.s)
        self.c_1_t = (2000 * self.H_oc) ** 0.5 
        self.a_1_t = (self.k * self.tochka_1_t.P * self.tochka_1_t.v * (10 ** 6)) ** 0.5
        self.M_1_t = self.c_1_t / self.a_1_t
        self.mu_1_pred = 0.97 
        self.F_1_pred = (self.G  * self.tochka_1_t.v)/(self.mu_1_pred * self.c_1_t)
        self.u = self.d * np.pi * self.n
        self.c_f = (2000 * H_o) ** 0.5
        self.u_c_f = self.u/self.c_f
        
    def calc_sopl_lop(self, b_1, alpha_1_e, t_1opt):
        
        self.b_1 = b_1 * (10 ** (-3)) 
        self.t_1opt = t_1opt
        self.alpha_1_e = alpha_1_e
        self.el_1 = self.F_1_pred/(np.pi * self.d * np.sin(np.deg2rad(self.alpha_1_e)))
        self.e_opt = 6 * (self.el_1 ** 0.5)
        self.e_opt = np.where(self.e_opt < 0.85,self.e_opt , 0.85) 
        self.l_1 = self.el_1/self.e_opt
        self.z_1 = np.round(( np.pi * self.d * self.e_opt)/( self.b_1 * self.t_1opt),0)
        self.z_1 = np.where(self.z_1 % 2 == 0,self.z_1, self.z_1 +1)
        self.t_1 = (np.pi * self.d * self.e_opt) / (self.b_1 * self.z_1)
        self.mu_1 = 0.982 - 0.005 * (self.b_1/self.l_1)
        self.b_1_l_1 = self.b_1 / self.l_1
        self.zit_sum = -29.63 * ((self.t_1opt)**3) + 84.444 * ((self.t_1opt)**2) - 76.434 * (self.t_1opt) + 24.457
        self.zit_prof = 19.394 * ((self.mu_1_pred)**3) - 40.997 * ((self.mu_1_pred)**2) + 26.639 * (self.mu_1_pred) - 2.8748
        self.alpha_yst = self.alpha_1_e - 16 * (self.t_1opt - 0.75) + 23.1
        self.fi = (1 - (self.zit_sum)/100) ** 0.5
        self.c_1 = self.c_1_t * self.fi
        self.alpha_1 = np.arcsin((self.mu_1/self.fi) * np.sin(np.deg2rad(self.alpha_1_e)))
        self.w_1 = ((self.c_1 ** 2) + (self.u ** 2) - 2 * self.c_1 * self.u * np.cos(self.alpha_1))**0.5
        self.tan_bett_1 = np.sin((self.alpha_1))/(np.cos((self.alpha_1)) - (self.u/self.c_1))
        self.bett_1 = np.arctan(self.tan_bett_1)
        return self.fi, self.mu_1, self.l_1 , self.alpha_yst
        
    def calc_for_find_koef_rab_lop(self, pere):
        self.pere = pere
        self.delt_H_c = self.c_1_t ** 2 * ( 1 - self.fi ** 2 ) / 2000    
        h_1 = self.h_1_t + self.delt_H_c
        self.tochka_1 = par(P = self.tochka_1_t.P, h = h_1) 
        h_2_t = self.tochka_1.h - self.H_or
        self.tochka_2_t = par(h = h_2_t, s = self.tochka_1.s) 
        self.w_2_t =(2000 * self.H_or + self.w_1 ** 2) ** 0.5
        self.l_2 = self.l_1 + self.pere
        self.a_2_t = (self.k * self.tochka_2_t.P * self.tochka_2_t.v * (10 ** 6)) ** 0.5
        self.m_2_pred = self.w_2_t/ self.a_2_t
        
        return self.m_2_pred, self.delt_H_c, self.H_oc, self.mu_1
        
    def calc_rab_lop(self, b_2, t_opt_2):
        
        self.b_2 = b_2 * (10 ** (-3)) 
        
        
        self.t_opt_2 = t_opt_2 
        self.b_2_l_2 = self.b_2/self.l_2
        self.mu_2 = 0.965 - 0.01*(self.b_2_l_2)
        self.F_2 = (self.G * self.tochka_2_t.v)/(self.mu_2 * self.w_2_t)
        self.a_2_t = (self.k * self.tochka_2_t.P * self.tochka_2_t.v * (10 ** 6)) ** 0.5 
        self.betta_2e = np.arcsin(self.F_2/(self.e_opt * np.pi * self.d * self.l_2 ))
        self.z_2 = np.round((np.pi * self.d)/(self.b_2 * self.t_opt_2),0) 
        self.z_2 = np.where(self.z_2 % 2 == 0,self.z_2, self.z_2 +1)
        self.t_2 = (np.pi * self.d)/(self.b_2 * self.z_2)
        self.bet_yst = np.rad2deg(self.betta_2e) - 3.25 * ( self.t_opt_2 - 0.7) + 63.8
        self.zit_sum = -434.49 * ((self.m_2_pred)**3) + 1061.5 * ((self.m_2_pred)**2) - 814.76 * (self.m_2_pred) + 205.42
        self.zit_prof = 4.8786 * (self.b_2_l_2) + 4.9714
        self.ksi = (1 - (self.zit_sum /100)) ** 0.5 
        self.w_2 = self.w_2_t * self.ksi
        self.bett_2 = np.arcsin((self.mu_2/self.ksi)* np.sin(self.betta_2e))
        self.c_2 = ((self.w_2**2) + (self.u**2) - 2 * self.w_2 * self.u * np.cos(self.bett_2)) ** 0.5
        self.alpha_2 = np.arctan((np.sin(self.bett_2))/((np.cos(self.bett_2)) - (self.u/self.w_2)))
        
        return self.w_2, self.c_2, self.bett_2, self.alpha_2, self.F_2, self.z_2, self.l_2, self.bet_yst

    def calc_eff(self,xi_vs):

        self.xi_vs = xi_vs
        self.delt_H_p = self.w_2_t ** 2 * (1 - (self.ksi ** 2)) /2000
        self.delt_H_vc = (self.c_2 ** 2)/(2 * 1000)
        self.E_o = self.H_o - self.xi_vs * self.delt_H_vc
        self.eff_oi_loss = (self.E_o - self.delt_H_c - self.delt_H_p - (1 - self.xi_vs) *  self.delt_H_vc) / self.E_o
        self.eff_oi_velocity = self.u * (self.c_1 * (np.cos(self.alpha_1)) + self.c_2 * (np.cos(self.alpha_2))) / (self.E_o * 1000)
        self.u_c_f_opt = (self.fi * np.cos(self.alpha_1)) /(2 * (1 - self.ro) ** 0.5) 
        h_2 = self.tochka_2_t.h + self.delt_H_p
        self.tochka_2 = par( P = self.tochka_2_t.P, h = h_2)
        return self.eff_oi_loss, self.eff_oi_velocity
        
    def calc_int_power(self,z,i):
        
        self.d_p = self.d + self.l_2
        self.mu_a = 0.5
        self.delt_a = 0.0025
        self.mu_r = 0.75
        self.delt_g = 0.001 * self.d_p
        self.z = z
        self.i = i
        self.delt_e = ((1/(self.mu_a * self.delt_a)**2) + (self.z/(self.mu_r * self.delt_g)**2)) ** (-0.5)
        self.ksi_b_y = (np.pi * self.d_p * self.delt_e * self.eff_oi_loss) * ((self.ro + 1.8 * (self.l_2/self.d)) ** 0.5) / self.F_1_pred
        self.delt_H_y = self.ksi_b_y * self.E_o
        self.k_tr = 0.7 * (10 ** -3)
        self.ksi_tr = (self.k_tr * ((self.d) ** 2)) * ((self.u_c_f) ** 3) / self.F_1_pred
        self.delt_H_tr = self.ksi_tr * self.E_o
        self.k_v = 0.065
        self.m = 1
        self.zitta_v = (self.k_v * (1 - self.e_opt) * ((self.u_c_f) ** 3) * self.m) / ((np.sin(np.rad2deg(self.alpha_1_e))) * self.e_opt)
        self.B_2 = self.b_2 * np.sin(np.deg2rad(self.bet_yst))
        self.ksi_segm = 0.25 * (self.B_2 * self.l_2 * self.u_c_f * self.i * self.eff_oi_loss) / self.F_1_pred
        self.ksi_parc = self.zitta_v + self.ksi_segm
        self.delt_H_parc = self.ksi_parc * self.E_o
        self.H_i = self.E_o - self.delt_H_c - self.delt_H_p - (1 - self.xi_vs) *  self.delt_H_vc - self.delt_H_y - self.delt_H_tr - self.delt_H_parc
        self.kpd_oi = self.H_i/self.E_o
        self.N_i = self.G * self.H_i
        
        return self.N_i, self.kpd_oi
        
    def calc_stength(self, W_atl, b_2_atl, sigma_bend_dop):

        self.sigma_bend_dop = sigma_bend_dop
        self.W_atl = W_atl * (10 ** (-6))
        self.b_2_atl = b_2_atl * (10 **(-3))
        self.W_min = self.W_atl * ((self.b_2/self.b_2_atl) ** 3)
        self.sigma_bend = (self.G * self.H_o * self.eff_oi_velocity * self.l_2)/(2 * self.u * self.z_2 * self.W_min * self.e_opt)/1000
        self.b_2 = np.where(self.sigma_bend > self.sigma_bend_dop,self.b_2 ,self.b_2 * (self.sigma_bend / self.sigma_bend_dop) ** 0.5) 
        self.w = 2 * np.pi * self.n 
        self.dens = 7800
        self.sigma_streth = 0.5 * self.dens * ((self.w) ** 2) * self.d * self.l_2 * (10 ** (-6))
        
        return self.sigma_bend, self.sigma_streth, self.b_2,self.betta_2e
        
    def tabl(self):

        self.tabl_one = pd.DataFrame({
            "Абсолютная скорость (с), м/с": [np.round(self.c_1,3), np.round(self.c_2,3)],
            "Абсолютный угол (α) , °": [np.round(np.rad2deg(self.alpha_1),3), np.round(np.rad2deg(self.alpha_2))],
            "Относительная скорость выхода (w) , м/с": [np.round(self.w_1,3), np.round(self.w_2,3)],
            "Относительный угол (β) , °": [np.round(np.rad2deg(self.bett_1),3), np.round(np.rad2deg(self.bett_2),3)],
            "": ["За сопловой решеткой", "За рабочей решеткой"],
        }).set_index("")
        
        self.tabl_two = pd.DataFrame({
            "Потери, %": [np.round(self.delt_H_c/self.E_o * 100,3), np.round(self.delt_H_p/self.E_o * 100,3),np.round(self.delt_H_vc/self.E_o * 100 ,3),np.round(self.delt_H_y/self.E_o * 100,3),np.round(self.delt_H_tr/self.E_o * 100,3),np.round(self.delt_H_parc/self.E_o * 100,3)],
            "": ["В сопловой решетке", "В рабочей решетке","Энергия выходной скорости","Абсолютные потери от утечек через периферийное уплотнение ступени", "Абсолютные потери от трения диска","Абсолютные потери от парциальности" ],
        }).set_index("")

        self.tabl_three =  pd.DataFrame({
            "P, МПа": [np.round(self.tochka_o.P,3), np.round(self.tochka_1_t.P,3),np.round(self.tochka_1.P,3), np.round(self.tochka_2_t.P,3), np.round(self.tochka_2.P,3) ],
            "T, К, °": [np.round(self.tochka_o.T,3), np.round(self.tochka_1_t.T,3),np.round(self.tochka_1.T,3), np.round(self.tochka_2_t.T,3), np.round(self.tochka_2.T,3) ],
            "h, кДж/кг": [np.round(self.tochka_o.h,3), np.round(self.tochka_1_t.h,3),np.round(self.tochka_1.h,3), np.round(self.tochka_2_t.h,3), np.round(self.tochka_2.h,3)],
            "s,кДж/(кг*К)": [np.round(self.tochka_o.s,3), np.round(self.tochka_1_t.s,3),np.round(self.tochka_1.s,3), np.round(self.tochka_2_t.s,3), np.round(self.tochka_2.s,3)],
            "": ["Параметры перед ступенью","Теоретичекие параметры за сопловой решеткой", "Параметры за сопловой решеткой","Теоретичекие параметры за ступенью", "Параметры за ступенью"],
        }).set_index("")
        
        self.tabl_four = pd.DataFrame({
            "Длины лопаток, м": [np.round(self.l_1,3), np.round(self.l_2,3)],
            "Число лопаток, шт": [self.z_1, self.z_2],
            "": ["Cопловая решетка", "Рабочая решетка"],
        }).set_index("")
        
        return self.tabl_one, self.tabl_two, self.tabl_three, self.tabl_four

    def points(self):
        
        h_0 = self.tochka_2_t.h + self.delt_H_vc + self.delt_H_y + self.delt_H_tr + self.delt_H_parc + self.delt_H_p 
        self.tocka_0 = par(h = h_0, P = self.tochka_2_t.P)

        return self.tochka_1_t, self.tochka_1, self.tochka_2_t, self.tochka_2, self.tocka_0

    def triangle(self):

        fig, ax = plt.subplots(1, 1, figsize=(15,5))
        
        c_1_u = self.c_1 * np.cos(self.alpha_1)
        c_1_a = self.c_1 * np.sin(self.alpha_1)
        w_1_u = self.w_1 * np.cos(self.bett_1)
        w_1_a = self.w_1 * np.sin(self.bett_1)
        
        c_2_u = self.c_2 * np.cos(self.alpha_2)
        c_2_a = self.c_2 * np.sin(self.alpha_2)
        w_2_u = self.w_2 * np.cos(self.bett_2)
        w_2_a = self.w_2 * np.sin(self.bett_2)

        ax.plot([0,-c_1_u], [0, -c_1_a], label='c_1', color='red')
        ax.plot([0,-w_1_u], [0, - w_1_a], label='w_1', color='green')
        ax.plot([-c_1_u, -c_1_u + self.u], [- w_1_a, - w_1_a], label='u', color='blue')

        ax.plot([0,c_2_u], [0, - c_2_a], label='c_2', color='red')
        ax.plot([0, w_2_u], [0, - w_2_a], label='w_2', color='green')
        ax.plot([c_2_u, c_2_u + self.u], [- w_2_a, - w_2_a], label='u', color='blue')
        ax.legend()
        ax.grid()
        
    def for_razb(self):
        return self.d, self.l_1, self.l_2, self.fi, self.mu_1, self.mu_2, self.ro, self.alpha_1_e, self.n