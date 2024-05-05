from typing import List, Tuple, Optional
from iapws import IAPWS97 as gas
from iapws import IAPWS97
import matplotlib.pyplot as plt
import numpy as np
import math

MPa = 10 ** 6
kPa = 10 ** 3
unit = 1 / MPa
to_kelvin = lambda x: x + 273.15 if x else None

internal_efficiency = 0.85
mechanical_efficiency = 0.994
generator_efficiency = 0.99


#Part 1

def f_delta_p0(p0):
    delta_p0 = 0.05 * p0
    return delta_p0

def f_delta_p_middle(p_middle):
    delta_p_middle = 0.1 * p_middle
    return delta_p_middle

def f_delta_p_1(p_middle):
    delta_p_1 = 0.03 * p_middle
    return delta_p_1

def f_real_p0(p0):
    delta_p0 = f_delta_p0(p0)
    real_p0 = p0 - delta_p0
    return real_p0

def f_real_p1t(p_middle):
    delta_p_middle = f_delta_p_middle(p_middle)
    real_p1t = p_middle + delta_p_middle
    return real_p1t

def f_real_p_middle(p_middle):
    delta_p_1 = f_delta_p_1(p_middle)
    real_p_middle = p_middle - delta_p_1
    return real_p_middle

def f__point_0(p0, t0):
    _point_0 = gas(P = p0 * unit, T=to_kelvin(t0))
    return _point_0

def f_point_0(p0, t0):
    delta_p0 = f_delta_p0(p0)
    _point_0 = f__point_0(p0, t0)
    real_p0=f_real_p0(p0)
    point_0 = gas(P=real_p0 * unit, h=_point_0.h)
    return _point_0

def f_point_1t(p_middle,p0, t0):
    real_p1t = f_real_p1t(p_middle)
    _point_0 = f__point_0(p0, t0)
    point_1t = gas(P=real_p1t * unit, s=_point_0.s)
    return point_1t

def f_hp_heat_drop(p0, t0,p_middle):
    _point_0 = f__point_0(p0, t0)
    point_1t = f_point_1t(p_middle,p0, t0)    
    hp_heat_drop = (_point_0.h - point_1t.h) * internal_efficiency
    return hp_heat_drop

def f_h_1(p0, t0, p_middle):
    point_0 = f_point_0(p0, t0)
    hp_heat_drop = f_hp_heat_drop(p0, t0,p_middle)
    h_1 = point_0.h - hp_heat_drop
    return h_1

def f_point_1(p0, t0,p_middle):
    real_p1t = f_real_p1t(p_middle)
    h_1 = f_h_1(p0, t0,p_middle)
    point_1 = gas(P=real_p1t * unit, h=h_1)
    return point_1

def f__point_middle(p_middle, t_middle):
    _point_middle = gas(P=p_middle * unit, T=to_kelvin(t_middle))
    return _point_middle

def f_point_middle(p_middle, t_middle):
    real_p_middle=f_real_p_middle(p_middle)
    _point_middle=f__point_middle(p_middle, t_middle)
    point_middle = gas(P=real_p_middle * unit, h=_point_middle.h)
    return point_middle

def f_point_2t(pk,p_middle, t_middle):
    _point_middle=f__point_middle(p_middle, t_middle)
    point_2t = gas(P=pk * unit, s=_point_middle.s)
    return point_2t

def f_lp_heat_drop(pk,p_middle, t_middle):
    _point_middle=f__point_middle(p_middle, t_middle)
    point_2t=f_point_2t(pk,p_middle, t_middle)
    lp_heat_drop = (_point_middle.h - point_2t.h) * internal_efficiency
    return lp_heat_drop

def f_h_2(pk,p_middle, t_middle):
    point_middle=f_point_middle(p_middle, t_middle)
    lp_heat_drop=f_lp_heat_drop (pk,p_middle, t_middle)
    h_2 = point_middle.h - lp_heat_drop
    return h_2

def f_point_2(pk,p_middle, t_middle):
    h_2=f_h_2(pk,p_middle, t_middle)
    point_2 = gas(P=pk * unit, h=h_2)
    return point_2

def f_efficiency_hp(p0, t0,p_middle):
    _point_0 = f__point_0(p0, t0)
    point_1 = f_point_1(p0, t0, p_middle)
    point_1t = f_point_1t(p_middle,p0, t0)
    efficiency_hp = (_point_0.h - point_1.h) / (_point_0.h - point_1t.h)
    return efficiency_hp

def f_efficiency_lp(pk,p_middle, t_middle):
    _point_middle=f__point_middle(p_middle, t_middle)
    point_2 = f_point_2(pk,p_middle, t_middle)
    point_2t = f_point_2t(pk,p_middle, t_middle)
    efficiency_lp = (_point_middle.h - point_2.h) / (_point_middle.h - point_2t.h)
    return efficiency_lp

def f_point_k_water(pk):
    point_k_water = gas(P=pk * unit, x=0)
    return point_k_water

def f_point_feed_water(p_feed_water, t_feed_water):
    point_feed_water = gas(P=p_feed_water * unit, T=to_kelvin(t_feed_water))
    return point_feed_water

def f_numenator_without(pk,p_middle, t_middle):
    point_2 = f_point_2(pk,p_middle, t_middle)
    _point_middle=f__point_middle(p_middle, t_middle)
    point_k_water=f_point_k_water(pk)
    numenator_without = point_2.T * (_point_middle.s - point_k_water.s)
    return numenator_without

def f_denumenator_without(p_middle, t_middle,p0, t0, pk):
    point_0 = f_point_0(p0, t0)
    point_1t = f_point_1t(p_middle,p0, t0)
    point_middle=f_point_middle(p_middle, t_middle)
    point_k_water=f_point_k_water(pk)
    denumenator_without = (point_0.h - point_1t.h) + (point_middle.h - point_k_water.h)
    return denumenator_without

def f_without_part(p_middle, t_middle,p0, t0, pk):
    numenator_without=f_numenator_without(pk,p_middle, t_middle)
    denumenator_without=f_denumenator_without(p_middle, t_middle,p0, t0, pk)
    without_part = 1 - (numenator_without / denumenator_without)
    return without_part

def f_numenator_infinity(pk,p_middle, t_middle, p_feed_water, t_feed_water):
    point_2 = f_point_2(pk,p_middle, t_middle)
    _point_middle=f__point_middle(p_middle, t_middle)
    point_feed_water=f_point_feed_water(p_feed_water, t_feed_water)
    numenator_infinity = point_2.T * (_point_middle.s - point_feed_water.s)
    return numenator_infinity

def f_denumenator_infinity(p0, t0,p_middle, t_middle,p_feed_water, t_feed_water):
    point_0 = f_point_0(p0, t0)
    point_1t = f_point_1t(p_middle,p0, t0)
    point_middle=f_point_middle(p_middle, t_middle)
    point_feed_water=f_point_feed_water(p_feed_water, t_feed_water)
    denumenator_infinity = (point_0.h - point_1t.h) + (point_middle.h - point_feed_water.h)
    return denumenator_infinity

def f_infinity_part(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water):
    numenator_infinity=f_numenator_infinity(pk,p_middle, t_middle, p_feed_water, t_feed_water)
    denumenator_infinity=f_denumenator_infinity(p0, t0,p_middle, t_middle,p_feed_water, t_feed_water)
    infinity_part = 1 - (numenator_infinity / denumenator_infinity)
    return infinity_part

def f_ksi_infinity(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water):
    without_part=f_without_part(p_middle, t_middle,p0, t0, pk)
    infinity_part=f_infinity_part(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    ksi_infinity = 1 - (without_part / infinity_part)
    return ksi_infinity

def f_coeff(pk,p_middle, t_middle,p_feed_water, t_feed_water):
    point_feed_water=f_point_feed_water(p_feed_water, t_feed_water)
    point_2 = f_point_2(pk,p_middle, t_middle)
    coeff = (point_feed_water.T - point_2.T) / (to_kelvin(374.2) - point_2.T)
    return coeff

def f_ksi(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water):
    ksi_infinity=f_ksi_infinity(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    ksi = 0.89 * ksi_infinity
    return ksi
def f_eff_num(p0, t0,pk,p_middle, t_middle):
    hp_heat_drop=f_hp_heat_drop(p0, t0,p_middle)
    lp_heat_drop=f_lp_heat_drop(pk,p_middle, t_middle)
    eff_num = hp_heat_drop + lp_heat_drop
    return eff_num
def f_eff_denum(p0, t0,p_middle, t_middle,pk):
    hp_heat_drop=f_hp_heat_drop(p0, t0,p_middle)
    point_middle=f_point_middle(p_middle, t_middle)
    point_k_water=f_point_k_water(pk)
    eff_denum = hp_heat_drop + (point_middle.h - point_k_water.h)
    return eff_denum

def f_efficiency(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water):
    eff_num=f_eff_num(p0, t0,pk,p_middle, t_middle)
    eff_denum=f_eff_denum(p0, t0,p_middle, t_middle,pk)
    ksi=f_ksi(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    efficiency = (eff_num / eff_denum) * (1 / (1 - ksi))
    return efficiency

def f_estimated_heat_drop(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water):
    point_0 = f_point_0(p0, t0)
    point_feed_water=f_point_feed_water(p_feed_water, t_feed_water)
    point_middle=f_point_middle(p_middle, t_middle)
    efficiency=f_efficiency(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    point_1 = f_point_1(p0, t0, p_middle)
    estimated_heat_drop = efficiency * ((point_0.h - point_feed_water.h) + (point_middle.h - point_1.h))
    return estimated_heat_drop

def f_inlet_mass_flow(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water):
    estimated_heat_drop=f_estimated_heat_drop(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    inlet_mass_flow = electrical_power / (estimated_heat_drop * 1000 * mechanical_efficiency * generator_efficiency)
    return inlet_mass_flow

def f_condenser_mass_flow(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water):
    point_2 = f_point_2(pk,p_middle, t_middle)
    point_k_water=f_point_k_water(pk)
    efficiency=f_efficiency(pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    condenser_mass_flow = (electrical_power /((point_2.h - point_k_water.h) * 1000 * mechanical_efficiency * generator_efficiency)* ((1 / efficiency) - 1))
    return condenser_mass_flow


#Отрисовка легенды
def legend_without_duplicate_labels(ax: plt.Axes) -> None:
    
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))

#Отрисовка процесса расширения по точкам    
def plot_process(ax: plt.Axes, points: List[IAPWS97], **kwargs) -> None:

    ax.plot([point.s for point in points], [point.h for point in points], **kwargs)

#координаты изобары в hs осях    
def get_isobar(point: IAPWS97) -> Tuple[List[float], List[float]]:
   
    s = point.s
    s_values = np.arange(s * 0.9, s * 1.1, 0.2 * s / 1000)
    h_values = [gas(P=point.P, s=_s).h for _s in s_values]
    return s_values, h_values

#координаты изотермы для пара в hs осях
def _get_isoterm_steam(point: IAPWS97) -> Tuple[List[float], List[float]]:
   
    t = point.T
    p = point.P
    s = point.s
    s_max = s * 1.2
    s_min = s * 0.8
    p_values = np.arange(p * 0.8, p * 1.2, 0.4 * p / 1000)
    h_values = np.array([gas(P=_p, T=t).h for _p in p_values])
    s_values = np.array([gas(P=_p, T=t).s for _p in p_values])
    mask = (s_values >= s_min) & (s_values <= s_max)
    return s_values[mask], h_values[mask]

#координаты изотермы для влажного пара в hs осях
def _get_isoterm_two_phases(point: IAPWS97) -> Tuple[List[float], List[float]]:
    
    x = point.x
    p = point.P
    x_values = np.arange(x * 0.9, min(x * 1.1, 1), (1 - x) / 1000)
    h_values = np.array([gas(P=p, x=_x).h for _x in x_values])
    s_values = np.array([gas(P=p, x=_x).s for _x in x_values])
    return s_values, h_values

#координаты изотермы в hs осях
def get_isoterm(point) -> Tuple[List[float], List[float]]:

    if point.phase == 'Two phases':
        return _get_isoterm_two_phases(point)
    return _get_isoterm_steam(point)

#Отрисовка изобары и изотермы
def plot_isolines(ax: plt.Axes, point: IAPWS97) -> None:
    
    s_isobar, h_isobar = get_isobar(point)
    s_isoterm, h_isoterm = get_isoterm(point)
    ax.plot(s_isobar, h_isobar, color='green', label='Изобара')
    ax.plot(s_isoterm, h_isoterm, color='blue', label='Изотерма')

#Отрисовка точки на hs-диаграмме    
def plot_points(ax: plt.Axes, points: List[IAPWS97]) -> None:
    
    for point in points:
        ax.scatter(point.s, point.h, s=50, color="red")
        plot_isolines(ax, point)

#координаты линии с постоянной степенью сухости в hs осях
def get_humidity_constant_line(
    point: IAPWS97,
    max_p: float,
    min_p: float,
    x: Optional[float]=None
) -> Tuple[List[float], List[float]]:
   
    _x = x if x else point.x
    p_values = np.arange(min_p, max_p, (max_p - min_p) / 1000)
    h_values = np.array([gas(P=_p, x=_x).h for _p in p_values])
    s_values = np.array([gas(P=_p, x=_x).s for _p in p_values])
    return s_values, h_values

#Отрисовка изолинии для степеней сухости на hs-диаграмме
def plot_humidity_lines(ax: plt.Axes, points: List[IAPWS97]) -> None:

    pressures = [point.P for point in points]
    min_pressure = min(pressures) if min(pressures) > 700/1e6 else 700/1e6
    max_pressure = max(pressures) if max(pressures) < 22 else 22
    for point in points:
        if point.phase == 'Two phases':
            s_values, h_values = get_humidity_constant_line(point, max_pressure, min_pressure, x=1)
            ax.plot(s_values, h_values, color="gray")
            s_values, h_values = get_humidity_constant_line(point, max_pressure, min_pressure)
            ax.plot(s_values, h_values, color="gray", label='Линия сухости')
            ax.text(s_values[10], h_values[10], f'x={round(point.x, 2)}')

#Построить изобары и изотермы для переданных точек. Если степень сухости у точки не равна 1, то построется
#дополнительно линия соответствующей степени сухости            
def plot_hs_diagram(ax: plt.Axes, points: List[IAPWS97]) -> None:

    plot_points(ax, points)
    plot_humidity_lines(ax, points)
    ax.set_xlabel(r"S, $\frac{кДж}{кг * K}$", fontsize=14)
    ax.set_ylabel(r"h, $\frac{кДж}{кг}$", fontsize=14)
    ax.set_title("HS-диаграмма процесса расширения", fontsize=18)
    ax.legend()
    ax.grid()
    legend_without_duplicate_labels(ax)
    
    
    
############################################################################################################################ Part 2

def f_u(avg_diameter, n):
        u = 3.1415 * avg_diameter * n 
        return u
    #принимаем протность
    
ro=0.075

def f_H0_c(ro,H_0):
    H0_c = (1 - ro) * H_0
    return H0_c

def f_H0_p(ro,H_0):
    H0_p = H_0 * ro
    return H0_p

def f_h_1t(ro,H_0,p0, t0):
    H0_c=f_H0_c(ro,H_0)
    point_0=f_point_0(p0, t0)
    h_1t = point_0.h - H0_c/1000
    return h_1t

def f_c_1t(ro,H_0):
    H0_c=f_H0_c(ro,H_0)
    c_1t = (2 * H0_c)**(0.5)
    return c_1t

def f_point_1_t(ro,H_0,p0, t0):
    h_1t=f_h_1t(ro,H_0,p0, t0)
    point_0=f_point_0(p0, t0)
    point_1_t = gas(h = h_1t, s = point_0.s)
    return point_1_t

        
def f_a_1t(ro,H_0,p0, t0):
    k=1.4
    point_1_t = f_point_1_t(ro,H_0,p0, t0)
    a_1t = math.sqrt(k * (point_1_t.P * MPa) * point_1_t.v)
    return a_1t

def f_M_1t(ro,H_0,p0, t0):
    c_1t = f_c_1t(ro,H_0)
    a_1t = f_a_1t(ro,H_0,p0, t0)
    M_1t = c_1t / a_1t
    return M_1t

def f_F_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0):
    mu__1 = 0.97
    c_1t = f_c_1t(ro,H_0)
    point_1_t = f_point_1_t(ro,H_0,p0, t0)
    G_0 = f_inlet_mass_flow(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    F_1 = (G_0 * point_1_t.v) / (mu__1 * c_1t)
    return F_1

alpha_1e = 12

def f_el_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter):
    F_1=f_F_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0)
    el_1 = F_1 / (3.1415 * avg_diameter * math.sin(math.radians(alpha_1e)))
    return el_1

def f_e_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter):
    el_1=f_el_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    e_opt = 4 * math.sqrt(el_1)
    if e_opt > 0.85: e_opt = 0.85
    return e_opt

def f_l_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter):
    el_1=f_el_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    e_opt=f_e_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    l_1 = el_1 / e_opt
    return l_1

b_1=100* 10**-3

def f_mu_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,avg_diameter):
    l_1=f_l_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    mu_1 = 0.982 - 0.005 * (b_1  / l_1)
    return mu_1

alpha_0=90
t_1_opt = 0.75

def f_z_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,b_1,t_1_opt):
    el_1=f_el_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    e_opt = f_e_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    z__1 = (3.1415 * avg_diameter * e_opt) / (b_1 * t_1_opt) 
    z_1 = round(z__1+0.5)-1 if (round(z__1) % 2) else round(z__1+0.5)
    return z_1

def f_t_1opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,b_1,
             t_1_opt):
    e_opt=f_e_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    z_1 = f_z_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,b_1,
                t_1_opt)
    t_1opt = (3.1415 * avg_diameter * e_opt) / (b_1 *  z_1)
    return t_1opt

def f_alpha_ust(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,b_1,
             t_1_opt):
    t_1opt = f_t_1opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,b_1,
             t_1_opt)
    alpha_ust =alpha_1e - 10 * (t_1opt - 0.75) + 21.2
    return alpha_ust

ksi_sum = 8

def f_fi(ksi_sum):
    fi = (1 - (ksi_sum/100))**(0.5)
    return fi

def f_fi_(b_1,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter):
    l_1=f_l_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    fi_ = 0.98 - 0.008 * (b_1 / l_1)
    return fi_

def f_delta_fi(ksi_sum,b_1,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,
               avg_diameter):
    fi=f_fi(ksi_sum)
    fi_=f_fi_(b_1,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    delta_fi = ((fi - fi_) / fi) * 100
    return delta_fi

def f_c_1(ksi_sum,ro,H_0):
    fi=f_fi(ksi_sum)
    c_1t = f_c_1t(ro,H_0)
    c_1 = c_1t * fi 
    return c_1

def f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter):
    fi=f_fi(ksi_sum)
    mu_1=f_mu_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,avg_diameter)
    alpha_1 = math.degrees(math.asin((mu_1/fi)* math.sin(math.radians(alpha_1e))))
    return alpha_1

def f_w_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n):
    alpha_1=f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,
                      ksi_sum,avg_diameter)
    c_1=f_c_1(ksi_sum,ro,H_0)
    u=f_u(avg_diameter, n)
    w_1 = (c_1 ** 2 + u ** 2 - 2 * c_1 * u * math.cos(math.radians(alpha_1)))**(0.5)
    return w_1

def f_beta_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n):
    alpha_1=f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter)
    c_1=f_c_1(ksi_sum,ro,H_0)
    u=f_u(avg_diameter, n)
    beta_1 = math.degrees(math.atan(math.sin(math.radians(alpha_1)) / (math.cos(math.radians(alpha_1)) - u / c_1)))
    return beta_1

def f_delta_H_c(ksi_sum,ro,H_0):
    c_1t=f_c_1t(ro,H_0)
    fi=f_fi(ksi_sum)
    delta_H_c = c_1t ** 2 / 2 * (1 - fi ** 2)
    return delta_H_c

def f_h_1_(ro,H_0,p0, t0,ksi_sum):
    point_1_t=f_point_1_t(ro,H_0,p0, t0)
    delta_H_c=f_delta_H_c(ksi_sum,ro,H_0)
    h_1 = point_1_t.h + delta_H_c/1000
    return h_1

def f_point_1_(ro,H_0,p0, t0,ksi_sum):
    point_1_t=f_point_1_t(ro,H_0,p0, t0)
    h_1=f_h_1_(ro,H_0,p0, t0,ksi_sum)
    point_1_ = gas(P = point_1_t.P, h = h_1)
    return point_1_

def f_h_2t(ro,H_0,p0, t0,ksi_sum):
    point_1_=f_point_1_(ro,H_0,p0, t0,ksi_sum)
    H0_p=f_H0_p(ro,H_0)
    h_2t = point_1_.h - H0_p/1000
    return h_2t

def f_point_2_t(ro,H_0,p0, t0,ksi_sum):
    point_1_=f_point_1_(ro,H_0,p0, t0,ksi_sum)
    h_2t=f_h_2t(ro,H_0,p0, t0,ksi_sum)
    point_2_t = gas(s = point_1_.s, h = h_2t)
    return point_2_t

def f_w_2t(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n):
    H0_p=f_H0_p(ro,H_0)
    w_1=f_w_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n)
    w_2t = (2* H0_p + w_1 ** 2)**(0.5)
    return w_2t

delta = 0.004

def f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter):
    l_1=f_l_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    delta = 0.004
    l_2 = l_1 + delta
    return l_2

#b_2 = 48 * 10**(-3)
b_2 = 0.08061977786706132

def f_mu_2(b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter):
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    mu_2 = 0.965 - 0.01 * (b_2 / l_2)
    return mu_2

def f_a_2t(ro,H_0,p0, t0,ksi_sum):
    k2 = 1.3
    point_2_t=f_point_2_t(ro,H_0,p0, t0,ksi_sum)
    a_2t = math.sqrt(k2 * (point_2_t.P * MPa) * point_2_t.v)
    return a_2t

def f_M_2t(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n):
    w_2t=f_w_2t(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n)
    a_2t=f_a_2t(ro,H_0,p0, t0,ksi_sum)
    M_2t = w_2t / a_2t
    return M_2t

def f_F_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n, b_2):
    w_2t=f_w_2t(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n)
    G_0 = f_inlet_mass_flow(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    point_2_t=f_point_2_t(ro,H_0,p0, t0,ksi_sum)
    mu_2=f_mu_2(b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    F_2 = (G_0 * point_2_t.v) / (mu_2 * w_2t)
    return F_2

def f_beta_2e(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water, ro, H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n, b_2):
    F_2=f_F_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n, b_2)
    e_opt=f_e_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    beta_2e = math.degrees(F_2 / (e_opt * 3.1415 * avg_diameter * l_2))
    return beta_2e

beta_0 = 25
beta_2e = 23
t_opt = 0.6
M_2t = 0.9
f_2 = 1.85 * (10**(-4))
I_2_min = 0.205 * 10**(-8)
J = I_2_min
W_2_min = 0.234*10**(-6)
t_2opt = 0.6

def f_z_2(avg_diameter,b_2,t_2opt):
    z__2 = (3.1415 * avg_diameter) / (b_2 * t_2opt)
    z_2 = round(z__2+0.5)-1 if (round(z__2) % 2) else round(z__2+0.5)
    return z_2

def f_t_2opt(avg_diameter,b_2,t_2opt):
    z_2 = f_z_2(avg_diameter,b_2,t_2opt)
    t_2opt = (3.1415 * avg_diameter) / (b_2 * z_2)
    return t_2opt

def f_beta2_ust(beta_2e,avg_diameter,b_2,t_2opt):
    t_2opt=f_t_2opt(avg_diameter,b_2,t_2opt)
    beta2_ust=beta_2e - 12.8 * (t_2opt - 0.65) + 58
    return beta2_ust

#ksi_grid = 5.5*10**(-2)
ksi_sum_g = 9.7

def f_psi(ksi_sum_g):
    psi = math.sqrt(1 - (ksi_sum_g/100))
    return psi

def f_psi_(b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter):
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter)
    psi_ = 0.96 - 0.014 * (b_2 / l_2)
    return psi_

def f_delta_psi(b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter):
    psi=f_psi(ksi_sum_g)
    psi_=f_psi_(b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter)
    delta_psi = ((psi - psi_) / psi) * 100
    return delta_psi

def f_w_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n,ksi_sum_g):
    w_2t=f_w_2t(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n)
    psi=f_psi(ksi_sum_g)
    w_2 = w_2t * psi
    return w_2

def f_beta_2(beta_2e,b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,ksi_sum_g,avg_diameter):
    mu_2=f_mu_2(b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    psi=f_psi(ksi_sum_g)
    beta_2 = math.degrees(math.asin((mu_2 / psi) * math.sin(math.radians(beta_2e))))
    return beta_2

def f_c_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    w_2=f_w_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g)
    u=f_u(avg_diameter, n)
    beta_2=f_beta_2(beta_2e,b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,ksi_sum_g,avg_diameter)
    c_2 = (w_2 ** 2 + u ** 2 - 2 * w_2 * u * math.cos(math.radians(beta_2)))**(0.5)
    return c_2

def f_alpha_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    beta_2=f_beta_2(beta_2e,b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,ksi_sum_g,avg_diameter)
    w_2=f_w_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g)
    u=f_u(avg_diameter, n)
    alpha_2 = math.degrees(math.atan((math.sin(math.radians(beta_2))) / (math.cos(math.radians(beta_2)) - u / w_2)))
    return alpha_2

def f_delta_H_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g):
    w_2t=f_w_2t(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n)
    psi=f_psi(ksi_sum_g)
    delta_H_p = w_2t ** 2 / 2 * (1 - psi ** 2)                                 
    return delta_H_p

def f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g):
    delta_Hp = f_delta_H_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g)
    h_2 = f_point_2_t(ro,H_0,p0, t0,ksi_sum).h + (delta_Hp/1000)
    _p_ = f_point_2_t(ro,H_0,p0, t0,ksi_sum).P
    point_2_ = gas(P = 17.51038032463063, h = 3259.6451629779685)
    return point_2_
       
def f_delta_H_vc(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    c_2=f_c_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_vc = (c_2 ** 2 )/ 2
    return delta_H_vc

x_vc = 0

def f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    delta_H_vc=f_delta_H_vc(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    E_0 = H_0 - x_vc * delta_H_vc 
    return E_0

def f_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_vc = f_delta_H_vc(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, b_1, 
                             ksi_sum,alpha_1e, avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_c = f_delta_H_c(ksi_sum,ro,H_0)
    delta_H_p = f_delta_H_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g)
    eff_ol = (E_0 - delta_H_c - delta_H_p - (1 - x_vc) * delta_H_vc) / E_0
    return eff_ol

def f_eff_ol_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    u=f_u(avg_diameter, n)
    c_1 = f_c_1(ksi_sum,ro,H_0)
    alpha_1=f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter)
    c_2=f_c_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    alpha_2 = f_alpha_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    eff_ol_ = ((u * (c_1 * math.cos(math.radians(alpha_1)) + c_2 * math.cos(math.radians(alpha_2)))) / E_0)                                                                                                   + 0.0115
    return eff_ol_

# def f_eff_ol_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
#     E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
#     u=f_u(avg_diameter, n)
#     w_2=f_w_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g)
#     w_1 = f_w_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n)
#     beta_1 = f_beta_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n)
#     beta_2 = f_beta_2(beta_2e,b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,ksi_sum_g)
#     eff_ol_ = (u * (w_1 * math.cos(math.radians(beta_1)) + w_2 * math.cos(math.radians(beta_2)))) / E_0  
#     return eff_ol_

def f_delta_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    eff_ol = f_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2) 
    eff_ol_ = f_eff_ol_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_eff_ol = ((eff_ol - eff_ol_) / eff_ol)*100 
    return delta_eff_ol

def f_c_f(H_0):
    c_f = math.sqrt(2 *  H_0)
    return c_f

def f_u_cf(avg_diameter, n, H_0):
    u=f_u(avg_diameter, n)
    c_f = f_c_f(H_0)
    u_cf = u / c_f
    return u_cf

def f_u_cf_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum):
    fi=f_fi(ksi_sum)
    alpha_1=f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter)
    u_cf_opt = fi * math.cos(math.radians(math.radians(alpha_1))) / (2 * math.sqrt(1 - ro))
    return u_cf_opt

def f_diam_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter):
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    diam_p = avg_diameter + l_2
    return diam_p

def f_delta_r(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter):
    diam_p = f_diam_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter)
    delta_r = 0.001 * diam_p
    return delta_r

def f_delta_e(mu_a,delta_a,z,mu_r, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter):
    delta_r = f_delta_r(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter)
    delta_e = math.pow((1 / (mu_a * delta_a) ** 2) + (z / (mu_r * delta_r) ** 2), -0.5)
    return delta_e

def f_ksi_b(mu_a,delta_a,z,mu_r, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter,n,ksi_sum_g,beta_2e,b_2):
    delta_e = f_delta_e(mu_a,delta_a,z,mu_r, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter)
    diam_p = f_diam_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter)
    eff_ol = f_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    F_1=f_F_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0)
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    ksi_b = ((3.1415 *diam_p *delta_e * eff_ol) /F_1) * math.sqrt(ro + 1.8 * l_2 / avg_diameter)
    return ksi_b

def f_delta_H_y(mu_a,delta_a,z,mu_r, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter,n,ksi_sum_g,beta_2e,b_2):
    E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    ksi_b = f_ksi_b(mu_a,delta_a,z,mu_r, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter,n,ksi_sum_g,beta_2e,b_2)
    delta_H_y = ksi_b* E_0
    return delta_H_y

def f_ksi_tr(k_tr,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,avg_diameter, n, H_0):
    u_cf = f_u_cf(avg_diameter, n, H_0)
    F_1=f_F_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0)
    ksi_tr = k_tr * (avg_diameter ** 2) * (u_cf ** 3) / F_1 
    return ksi_tr

def f_delta_H_tr(k_tr,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    ksi_tr = f_ksi_tr(k_tr,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,avg_diameter, n, H_0)
    E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_tr = ksi_tr * E_0
    return delta_H_tr

def f_ksi_v(k_v,m,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro, alpha_1e,b_1,ksi_sum, avg_diameter, n, H_0):
    alpha_1=f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter)
    e_opt=f_e_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    u_cf = f_u_cf(avg_diameter, n, H_0)
    ksi_v = (k_v / math.sin(math.radians(alpha_1e))) * ((1 - e_opt) / e_opt) * (u_cf ** 3) * m
    return ksi_v

def f_B_2(b_2,beta_2e,avg_diameter,t_2opt):
    beta2_ust = f_beta2_ust(beta_2e,avg_diameter,b_2,t_2opt)
    B_2 = b_2 * math.sin(math.radians(beta2_ust))
    return B_2

def f_ksi_segm(i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2 ):
    B_2 = f_B_2(b_2,beta_2e,avg_diameter,t_2opt)
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    F_1=f_F_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0)
    u_cf = f_u_cf(avg_diameter, n, H_0)
    eff_ol = f_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    ksi_segm = 0.25 * B_2 * 10**-3 * l_2 / F_1 * u_cf * eff_ol *i
    return ksi_segm

def f_ksi_part(k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    ksi_v = f_ksi_v(k_v,m,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro, alpha_1e,b_1,ksi_sum,avg_diameter, n, H_0)
    ksi_segment = f_ksi_segm(i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    ksi_part = ksi_v + ksi_segment
    return ksi_part

def f_delta_H_part(k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    ksi_part = f_ksi_part(k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_part = ksi_part * E_0
    return delta_H_part

def f_H_i(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_c = f_delta_H_c(ksi_sum,ro,H_0)
    delta_H_p = f_delta_H_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g)
    delta_H_vc=f_delta_H_vc(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_y = f_delta_H_y(mu_a,delta_a,z,mu_r, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter,n,ksi_sum_g,beta_2e,b_2)
    delta_H_tr = f_delta_H_tr(k_tr,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_part = f_delta_H_part(k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    H_i = E_0 - delta_H_c - delta_H_p - (1 - x_vc) * delta_H_vc - delta_H_y - delta_H_tr - delta_H_part
    return H_i
                    
def f_eff_oi(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    H_i = f_H_i(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    E_0 = f_E_0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    eff_oi = H_i / E_0
    return eff_oi
                    
def f_N_i(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    H_i = f_H_i(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    G_0 = f_inlet_mass_flow(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    N_i = G_0 * H_i
    return N_i

b2_atl = 25.6* 10 **(-3)
W2_min_atl = 0.234 *10** (-6)

def f_W2_min(W2_min_atl,b_2,b2_atl):
    W2_min= W2_min_atl * ((b_2 / b2_atl)**3)
    return W2_min
                    
def f_b_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,W2_min_atl,b_2,b2_atl,t_2opt):
    sigma_izg = f_sigma_izg(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,W2_min_atl,b_2,b2_atl,t_2opt)
    b_2 = b_2*math.sqrt((sigma_izg/20))   
    return b_2
                    
def f_sigma_izg(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,W2_min_atl,b_2,b2_atl,t_2opt):
    G_0 = f_inlet_mass_flow(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water)
    eff_ol = f_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    u=f_u(avg_diameter, n)
    W2_min = f_W2_min(W2_min_atl,b_2,b2_atl)
    z_2 = f_z_2(avg_diameter,b_2,t_2opt)
    e_opt=f_e_opt(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    sigma_izg = ((G_0 * H_0 * eff_ol * l_2) / (2 * u * z_2 * W2_min * e_opt))/ (10**6)
    return sigma_izg
                    
def f_omega(n):
    omega = 2 * 3.1415 * n
    return omega
                    
def f_sigma_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,n,avg_diameter):
    omega = f_omega(n)
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    sigma_p = 0.5 * 7800 * (omega**2) * avg_diameter * l_2 / (10**6)
    return sigma_p

########################################################################################################################
def plot_treyg_speed(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,
                      ksi_sum, beta_2e,b_2,ksi_sum_g,avg_diameter, n):
        
        alpha_1 = f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,
                      ksi_sum,avg_diameter)
        beta_2 = f_beta_2(beta_2e,b_2,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,ksi_sum_g,avg_diameter)
        c_1=f_c_1(ksi_sum,ro,H_0)
        u=f_u(avg_diameter, n)
        w_2=f_w_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,
         avg_diameter, n,ksi_sum_g)
        
        sin_alpha_1 = math.sin(math.radians(alpha_1))
        cos_alpha_1 = math.cos(math.radians(alpha_1))
        sin_beta_2 = math.sin(math.radians(beta_2))
        cos_beta_2 = math.cos(math.radians(beta_2))

        c1_plot = [[0, -c_1 * cos_alpha_1], [0, -c_1 * sin_alpha_1]]
        u1_plot = [[-c_1 * cos_alpha_1, -c_1 * cos_alpha_1 + u], [-c_1 * sin_alpha_1, -c_1 * sin_alpha_1]]
        w1_plot = [[0, -c_1 * cos_alpha_1 + u], [0, -c_1 * sin_alpha_1]]
        w2_plot = [[0, w_2 * cos_beta_2], [0, -w_2 * sin_beta_2]]
        u2_plot = [[w_2 * cos_beta_2, w_2 * cos_beta_2 - u], [-w_2 * sin_beta_2, -w_2 * sin_beta_2]]
        c2_plot = [[0, w_2 * cos_beta_2 - u], [0, -w_2 * sin_beta_2]]

        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        ax.plot(c1_plot[0], c1_plot[1], label='C_1', c='red')
        ax.plot(u1_plot[0], u1_plot[1], label='u_1', c='blue')
        ax.plot(w1_plot[0], w1_plot[1], label='W_1', c='green') 
        ax.plot(w2_plot[0], w2_plot[1], label='W_2', c='green')
        ax.plot(u2_plot[0], u2_plot[1], label='u_2', c='blue')
        ax.plot(c2_plot[0], c2_plot[1], label='C_2', c='red')
        ax.set_title("Треугольник скоростей",)
        ax.legend()
        ax.grid();
        
######################################################################################################################################        
def plot_u_cf(n, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, ksi_sum_g,beta_2e,b_2):
        
    array_u_cf = []
    d_value = []
    efficiency_ = []
    efficiency = []
    
    for i in range(90, 110,2):
        u_cf = f_u_cf(i/100, n, H_0)   
        eff_ol = f_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,i/100, n,ksi_sum_g,beta_2e,b_2) 
        eff_ol_ = f_eff_ol(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,i/100, n,ksi_sum_g,beta_2e,b_2)

        array_u_cf = array_u_cf + [u_cf]
                
        d_value = d_value + [i/100]
                
        efficiency_ = efficiency_ + [eff_ol_]
        efficiency = efficiency + [eff_ol]
            
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    X1 = array_u_cf
    X2 = d_value
    Y1 = efficiency_
    Y2 = efficiency
    
    ax.plot(X1,Y1, label = 'Расчет по методу скоростей', color = 'blue', linewidth=3,linestyle=":")
    ax.plot(X1,Y2, label = 'Расчет по методу потерь', color = 'red')
    ax.set_title("Зависимость лопаточного КПД от u/сф")
    ax.set_ylabel("Лопаточный КПД")
    ax.set_xlabel("U/сф")
    ax.legend()
    ax.grid()
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    plt.plot(X2,Y1, label = 'Расчет по методу скоростей', color = 'blue',linewidth=3,linestyle=":")
    plt.plot(X2,Y2, label = 'Расчет по методу потерь', color = 'red')
    plt.title("Зависимость лопаточного КПД от d")
    plt.ylabel("Лопаточный КПД")
    plt.xlabel("d, м")
    plt.legend()
    plt.grid()
    plt.show()
    
######################################################################################################################################    

def f_point_t_konec(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g):
    point_0 = f_point_0(p0, t0)
    point_2_ = f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g)
    point_t_konec = gas(h = (point_0.h - H_0/1000), P = point_2_.P)
    return point_t_konec

def f_pont_poter(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, b_1, 
                             ksi_sum,alpha_1e, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    point_2_ = f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g)
    delta_H_vc = f_delta_H_vc(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, b_1, 
                             ksi_sum,alpha_1e, avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    delta_H_c = f_delta_H_c(ksi_sum,ro,H_0)
    delta_H_p = f_delta_H_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g)
    pont_poter = gas(P = point_2_.P,h = (point_2_.h + (delta_H_c + delta_H_p + (1 - 0) * delta_H_vc)/1000))
    return pont_poter
 
    
#Part 3
n_stages = 12

# K-800-23.5 lMZ

def f_diam_1(avg_diameter):
    diam_1 = avg_diameter - 0.2
    return diam_1

def f_diam_kor(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter):
    diam_1 = f_diam_1(avg_diameter)
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    diam_kor = diam_1 - l_2
    return diam_kor

def f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g):
    diam_kor = f_diam_kor(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    diam_p = f_diam_p(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e, avg_diameter)
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    point_1__ = f_point_1(p0, t0,p_middle)
    point_2__ = f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum,avg_diameter, n, ksi_sum_g)
    l_2z = (- diam_kor +  math.sqrt(diam_kor**2 + 4 * (point_1__.v / point_2__.v) * l_2 * diam_p))/2
    return l_2z

        
def f_d_z2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n):
    diam_kor = f_diam_kor(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    l_2z = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    d_z2 = diam_kor + l_2z
    return d_z2

def f_diams(n,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages,ksi_sum, num_turb):
    diam_1 = f_diam_1(avg_diameter)
    d_z2 = f_d_z2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n)
    diams = (-diam_1 + d_z2) / n_stages * num_turb + diam_1
    return diams

def f_l_b(n,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages, ksi_sum, num_turb):
    l_2=f_l_2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    l_2z = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    l_b = (-l_2 + l_2z) / n_stages * num_turb + l_2
    return l_b

def f_H0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages, ksi_sum, b_1,n):
        
    ucf = []
    H0 = []
    fi=f_fi(ksi_sum)
    alpha_1=f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter)
        
    for j in range (0, n_stages):
            
        diams = f_diams(n,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages, ksi_sum,num_turb = j)
        l_b = f_l_b(n,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages,ksi_sum, num_turb = j)
            
        reaction_degree = ro + (1.8 / (diams / l_b + 1.8))
        u_cf_n = fi * math.cos(math.radians(alpha_1)) /  (2 * (1 - reaction_degree) ** 0.5)
        first = (diams / (u_cf_n)) ** 2
        second = (n / 50) ** 2
        
        if j == 0:
            H0 = H0 + [ 12.325 * first * second ] 
        else:                
            H0 = H0 + [0.93 * 12.325 * first * second ]
    return H0       
            
def f_reheat_factor(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2,n_stages):
    eff_oi = f_eff_oi(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2)
    point_1 = f_point_1(p0, t0, p_middle)
    point_1_ = point_1_ = f_point_1_(ro,H_0,p0, t0,ksi_sum)
    reheat_factor = 4.8 * 10 ** (-4) * (1 - eff_oi) * (point_1_.h - point_1.h) * (n_stages - 1) / n_stages
    return reheat_factor

def f_H0_poln(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages,ksi_sum, b_1,n):
    
    H0 = f_H0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages, ksi_sum, b_1,n)
    H0_poln = sum(H0)
    
    return H0_poln

def f_n_st_rasch(n,mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages,ksi_sum, b_1):             
    reheat_factor = f_reheat_factor(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2,n_stages)
    point_1 = f_point_1(p0, t0, p_middle)
    point_1_ = point_1_ = f_point_1_(ro,H_0,p0, t0,ksi_sum)
    hp_heat_drop = point_1_.h - point_1.h
    H0_poln = f_H0_poln(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages,ksi_sum, b_1,n)

    n_st_rasch = (hp_heat_drop* (1 + reheat_factor))/(H0_poln/(n_stages))
    return n_st_rasch

def f_real_heat_drop(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
    
    reheat_factor = f_reheat_factor(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2,n_stages)
    point_1 = f_point_1(p0, t0, p_middle)
    point_1_ = point_1_ = f_point_1_(ro,H_0,p0, t0,ksi_sum)
    hp_heat_drop = point_1_.h - point_1.h
    H0_poln = f_H0_poln(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages,ksi_sum, b_1,n)

    h0 = (hp_heat_drop * (1 + reheat_factor) - H0_poln) / n_stages
        
    real_heat_drop = []
        
    for q in range(0,9):
        
        H0 = f_H0(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages, ksi_sum, b_1,n)

        real_heat_drop = real_heat_drop + [H0[q] + h0]
            
    return real_heat_drop

def f_delta_heat_drop(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2):
            
    real_heat_drop = f_real_heat_drop(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2)

    sum_r_h_d = sum(real_heat_drop)
    point_1 = f_point_1(p0, t0, p_middle)
    point_1_ = point_1_ = f_point_1_(ro,H_0,p0, t0,ksi_sum)
    hp_heat_drop = point_1_.h - point_1.h    
    delta_heat_drop = ((sum_r_h_d-hp_heat_drop)/sum_r_h_d)*100 
    
    return delta_heat_drop     
    
def plot_ot_z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter,n_stages,n,mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,  ksi_sum_g,beta_2e,b_2):
      
    d_plot = []
    ucf_plot = []
    veernost_plot = []
    react_plot = []
    blade_plot = []
    
    ro = 0.075
    fi=f_fi(ksi_sum)
    alpha_1=f_alpha_1(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,ksi_sum,avg_diameter)      
        
    for i in range (0,n_stages):
            
        diams = f_diams(n,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages, ksi_sum,num_turb = i)
        l_b = f_l_b(n,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n_stages, ksi_sum,num_turb=i) 
            
        d_plot = d_plot + [diams]
            
        blade_plot = blade_plot + [l_b]
            
        veernost_plot = veernost_plot +[diams/l_b]
    
        react_plot = react_plot + [ro + (1.8 / ((diams/l_b)+ 1.8))]

        ucf_plot = ucf_plot + [fi * math.cos(math.radians(alpha_1)) /  (2 * (1 -(ro + (1.8 / (diams/l_b + 1.8))))**0.5)]

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    X1 = range(1,n_stages+1)
        
    Y1 = d_plot
    Y2 = blade_plot
    Y3 = veernost_plot
    Y4 = react_plot      
    
    Y5 = ucf_plot
        
    real_heat_drop = f_real_heat_drop(mu_a,delta_a,z,mu_r,k_tr,k_v,m,i,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, ksi_sum, avg_diameter, n,ksi_sum_g,beta_2e,b_2)

    Y6 = real_heat_drop

    ax.plot(X1,Y1, color = 'k')
        
    ax.set_title("Зависимость диаметра лопаток от номера ступени")
    ax.set_ylabel("Диаметр")
    ax.set_xlabel("z")
    ax.plot(range(1, n_stages + 1), Y1,  marker='o')
    ax.grid()
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.plot(X1,Y2, color = 'k')
    plt.title("Зависимость длинны лопаток от номера ступени")
    plt.ylabel("Длинна лопатки")
    plt.xlabel("z")
    ax.plot(range(1, n_stages + 1), Y2,  marker='o')
    plt.grid()
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.plot(X1,Y3, color = 'k')
    plt.title("Зависимость d/l от номера ступени")
    plt.ylabel("d/l")
    plt.xlabel("z")
    ax.plot(range(1, n_stages + 1), Y3,  marker='o')
    plt.grid()
    plt.show()
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.plot(X1,Y4, color = 'k')
    plt.title("Зависимость степени реактивности от номера ступени")
    plt.ylabel("Степень реактивности")
    plt.xlabel("z")
    ax.plot(range(1, n_stages + 1), Y4,  marker='o')
    plt.grid()
    plt.show()
        
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.plot(X1,Y5, color = 'k')
    plt.title("Зависимость u/cf от номера ступени")
    plt.ylabel("u/cf")
    plt.xlabel("z")
    ax.plot(range(1, n_stages + 1), Y5,  marker='o')
    plt.grid()
    plt.show()
        
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.plot(X1,Y6, color = 'k')
    plt.title("Зависимость теплоперепада от номера ступени")
    plt.ylabel("H")
    plt.xlabel("z")
    ax.plot(range(1,n_stages + 1 ), Y6,  marker='o')
    plt.grid()
    plt.show()
    
ro_metal = 7800 
max_stress_1 = 235 * MPa

def f_tension_stress_root_func(ro_metal,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n):
    l_2z = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    d_z2 = f_d_z2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n)
    u=f_u(avg_diameter, n)
    tension_stress_root_func = 0.5 * ro_metal * (u ** 2) * d_z2 * l_2z
    return tension_stress_root_func

def f_n_bend_last(max_stress_1, ro_metal,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n):
    tension_stress_root_func = f_tension_stress_root_func(ro_metal,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)
    n = max_stress_1 / tension_stress_root_func
    return n    

def f_r_k(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n):
    l_2z = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    d_z2 = f_d_z2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n)    
    r_k = (d_z2 - l_2z) / 2
    return r_k

def f_constant_part(ro_metal,avg_diameter, n):    
    u = f_u(avg_diameter, n)
    constant_part = ro_metal * (u ** 2)
    return constant_part

def f_left(x,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n):    
    l_2z = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    r_k  = f_r_k(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n)
    left = r_k * (l_2z - x)
    return left

def f_right(x,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n):
    l_2z = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    right = 0.5 * ((l_2z ** 2) - (x ** 2))
    return right 

def f_tens(x,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter, n, ro_metal,ksi_sum):
    constant_part = f_constant_part(ro_metal,avg_diameter, n)
    left = f_left(x,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n)
    right = f_right(x,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n)
    tens = constant_part * (left + right)
    return tens

nu = 0.3
sigma_1 = 0
sigma_2 = 100 * MPa
max_stress_2 = 510 * MPa
r_2 = 0.56
r = 0
#r_1 = f_diam_kor(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter)/2
   
    
def f_a(nu):
    a = (3 + nu) / 8
    return nu

def f_sigma_r(sigma_2,nu,avg_diameter, n,ro_metal,r_2,r):
    a = f_a(nu)
    u = f_u(avg_diameter, n)
    sigma_r = a * ro_metal * (u ** 2) * (r_2**2 - r**2) + sigma_2
    return sigma_r

def f_b(nu):
    b = (1 + 3 * nu) / (3 + nu)
    return b

def f_sigma_t(sigma_2,nu,avg_diameter, n,ro_metal,r_2,r):
    a = f_a(nu)
    b = f_b(nu)
    u = f_u(avg_diameter, n)
    sigma_t = a * ro_metal * (u ** 2) * (r_2**2 - b * r**2) + sigma_2
    return sigma_t

def f_n(max_stress_2,sigma_2,nu,avg_diameter, n,ro_metal,r_2,r):
    sigma_t = f_sigma_t(sigma_2,nu,avg_diameter, n,ro_metal,r_2,r)
    n = max_stress_2/sigma_t
    return n

psi = 0.85     
H = 0.6         
fi = 0.8

E = 2 * (10 ** 11)
m = 12
t = 25 * 10**(-3)
delta = 5 * 10**(-3)
b = 40 * 10**(-3)

def f_lambd(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,J,f_2,ksi_sum,n, ksi_sum_g):
    F = f_2 
    l = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    i =(J / F) ** 0.5
    lambd = l / i 
    return lambd

def f_J_b(b, delta):
    J_b = b * (delta ** 3) / 12
    return J_b

def f_k(b,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,beta_2e,avg_diameter,b_2,t_2opt,delta,J,m,H,E,t,ksi_sum,n, ksi_sum_g):
    beta = f_beta2_ust(beta_2e,avg_diameter,b_2,t_2opt)
    J_b = f_J_b(b, delta)
    l = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    k = (12 * (m - 1) * H * E * J_b * l * (math.sin(math.radians(beta))) ** 2) / (m * t * J * E)
    return k 

def f_NU(b,f_2,t,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum, n, ksi_sum_g):
    F = f_2
    l = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    NU = b * delta * t / (F * l)
    return NU

def f_B_bandage(b,f_2,t,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,beta_2e,b_2,t_2opt,n):
    nu = f_NU(b,f_2,t,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum, n, ksi_sum_g)
    d = f_d_z2(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n)
    l = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    beta = f_beta2_ust(beta_2e,avg_diameter,b_2,t_2opt)
    B_bandage = 0.5 * ((d /l) - 1) * ((nu+ 1/2)/(nu+1/3)) + (math.sin(math.radians(beta))) ** 2
    return B_bandage
    
def f_static_frequency(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,psi, E,f_2,J,ro_metal,ksi_sum,n):
    l = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    F = f_2
    first = psi * 0.56 / ( l ** 2)
    second = ((E * J) / (ro_metal * F)) ** 0.5
    static_frequency = first * second
    return static_frequency

def f_dynamic_frequency(f,n,b,f_2,t,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,beta_2e,b_2,t_2opt):
    B_bandage = f_B_bandage(b,f_2,t,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,beta_2e,b_2,t_2opt,n)
    dynamic_frequency = f  * (1 + B_bandage * (n / f) ** 2) ** 0.5
    return dynamic_frequency

def f_min_max(f,delta=0.05):
    min_ = f * (1 - delta)
    max_ = f * (1 + delta)
    return min_, max_
    
def f_k_line(k, n_line):   
    return k * n_line


mass = 10954.644
L = 5
L_r = 5.8
d = 560 * 10**(-3)
d_0 = 100 * 10**(-3)

def f_numenator(d,L):    
    numenator = (d / L) ** 2
    return numenator

def f_denumenator(mass,L):
    denumenator = (mass / L) ** 0.5
    return denumenator

def f_n_critical(d,mass,L):
    numenator = f_numenator(d,L)
    denumenator = f_denumenator(mass,L)
    n_critical  = 7.5 * numenator / denumenator
    return n_critical

def f_rps(d,mass,L):
    n_critical = f_n_critical(d,mass,L)
    rps = n_critical / 60
    return rps

E = 1.8 * 10**11

def f_EI(d,d_0):
    EI = E * 3.1415 * (d ** 4 - d_0 ** 4) / 64
    return EI

def f_P_11(L,d,d_0,mass,L_r):
    EI = f_EI(d,d_0)
    P_11 = (3.1415 ** 2 / L ** 2) * (EI / (mass / L_r)) ** 0.5
    return P_11
def f_P_12(L,d,d_0,mass,L_r):
    P_11 = f_P_11(L,d,d_0,mass,L_r)
    P_12 = 4 * P_11
    return P_12

delta_op = 0.5 * 10 ** -9
C_h = 0.5 * 10 ** 9

def f_delta_(C_h,delta_op):
    delta_ = delta_op + 1 / C_h
    return delta_
def f_P_21(C_h,delta_op):
    delta_ =f_delta_(C_h,delta_op) 
    P_21 = (2 / (mass * delta_)) ** 0.5
    return P_21
def f_P_22(C_h,delta_op,L,L_rotor):
    P_21 = f_P_21(C_h,delta_op)
    P_22 = (L / L_r) * 3 ** 0.5 * P_21
    return P_22
def f_P_1(L,mass,L_r,C_h,delta_op,d,d_0):
    P_11 = f_P_11(L,d,d_0,mass,L_r)
    P_21 = f_P_21(C_h,delta_op)
    P_1 = (1 / ((1 / P_11 ** 2) + (1 / P_21 ** 2) ) ** 0.5)/ (2 * 3.1415)
    return P_1
def f_P_2(C_h,delta_op,L,L_r,mass,d,d_0):
    P_12 = f_P_12(L,d,d_0,mass,L_r)
    P_22 = f_P_22(C_h,delta_op,L,L_r)
    P_2 = (1 / ((1 / P_12 ** 2) + (1 / P_22 ** 2) ) ** 0.5)// (2 * 3.1415)
    return P_2

'Первая критическая частота ротора''Вторая критическая частота ротора'


def plot_freq(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,psi, E,f_2,J,ro_metal,ksi_sum,b,beta_2e,b_2,t_2opt,n):
    
    static_frequency = f_static_frequency(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,psi, E,f_2,J,ro_metal,ksi_sum,n)
    
    f_a0 = static_frequency * fi
    n_line = np.linspace(0, 50)
    f = f_dynamic_frequency(f_a0,n_line,b,f_2,t,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,beta_2e,b_2,t_2opt)
    
    min_line_a0, max_line_a0 = f_min_max(f)
    
    fig, ax = plt.subplots(1,1,figsize=(10,10))

    ax.plot(n_line, f_dynamic_frequency(f_a0,n_line,b,f_2,t,electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,beta_2e,b_2,t_2opt), label='$f_{a0}$')
    ax.fill_between(n_line, y1=min_line_a0, y2=max_line_a0, alpha=0.5)

    ax.plot(n_line, f_k_line(1,n_line), label=f'k={1}')
    ax.plot(n_line, f_k_line(2,n_line), label=f'k={2}')
    ax.plot(n_line, f_k_line(3,n_line), label=f'k={3}')
    ax.plot(n_line, f_k_line(4,n_line), label=f'k={4}')
    ax.plot(n_line, f_k_line(5,n_line), label=f'k={5}')
    ax.plot(n_line, f_k_line(6,n_line), label=f'k={6}')
    ax.plot(n_line, f_k_line(7,n_line), label=f'k={7}')
    ax.plot(n_line, f_k_line(8,n_line), label=f'k={8}')
    ax.plot(n_line, f_k_line(9,n_line), label=f'k={9}')

    ax.set_xlabel("n, rps")
    ax.set_ylabel("f, Hz")
    ax.grid()
    ax.legend()
    ax.set_title("Вибрационная диаграмма");
    
def graf_tens(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g, ro_metal): 
    
    l_2z = f_l_2z(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,ksi_sum,n, ksi_sum_g)
    x = np.linspace(0, l_2z, 100)
    tens = f_tens(x, electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,avg_diameter,n, ro_metal,ksi_sum)
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    ax.plot(tens/ MPa, x)
    ax.set_xlabel("Напряжения растяжения, МПа")
    ax.set_ylabel("Координаты сечения вдоль высоты лопатки, м")
    ax.grid()
    
import pandas as pd

def f_table_parametrs(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,
                     ksi_sum,avg_diameter, n, ksi_sum_g,beta_2e,b_2):
    return pd.DataFrame({
    "P(Мпа)": [f_point_0(p0, t0).P, 
          f_point_1_(ro,H_0,p0, t0,ksi_sum).P,
          f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,
                     ksi_sum,avg_diameter, n, ksi_sum_g).P,
          f_point_t_konec(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, 
                          ksi_sum,avg_diameter, n, ksi_sum_g).P,
          f_pont_poter(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, b_1, 
                             ksi_sum,alpha_1e, avg_diameter, n,ksi_sum_g,beta_2e,b_2).P],
    "T(К)": [f_point_0(p0, t0).T, 
          f_point_1_(ro,H_0,p0, t0,ksi_sum).T,
          f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,
                     ksi_sum,avg_diameter, n, ksi_sum_g).T,
          f_point_t_konec(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, 
                          ksi_sum,avg_diameter, n, ksi_sum_g).T,
          f_pont_poter(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, b_1, 
                             ksi_sum,alpha_1e, avg_diameter, n,ksi_sum_g,beta_2e,b_2).T],
    "H(кДж/кг)": [f_point_0(p0, t0).h, 
          f_point_1_(ro,H_0,p0, t0,ksi_sum).h,
          f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, 
                     ksi_sum,avg_diameter, n, ksi_sum_g).h,
          f_point_t_konec(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, 
                          ksi_sum,avg_diameter, n, ksi_sum_g).h,
          f_pont_poter(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, b_1, 
                             ksi_sum,alpha_1e, avg_diameter, n,ksi_sum_g,beta_2e,b_2).h],
    "S(кДж/кг*К)": [f_point_0(p0, t0).s, 
          f_point_1_(ro,H_0,p0, t0,ksi_sum).s,
          f_point_2_(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1,
                     ksi_sum,avg_diameter, n, ksi_sum_g).s,
          f_point_t_konec(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, alpha_1e,b_1, 
                          ksi_sum,avg_diameter, n, ksi_sum_g).s,
          f_pont_poter(electrical_power,pk,p0, t0,p_middle, t_middle, p_feed_water, t_feed_water,ro,H_0, b_1, 
                             ksi_sum,alpha_1e, avg_diameter, n,ksi_sum_g,beta_2e,b_2).s],
    "point": ["0", "1_","2_","t_konec","poter"],}).set_index("point")
    