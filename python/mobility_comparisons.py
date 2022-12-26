import matplotlib.pyplot as plt
import math
import numpy as np
from collections import OrderedDict

linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

convert = 1e-12/3600
kb = 8.617e-5

def gg_calc(d0,df,t):
    return (df**2 - d0**2)/t * convert

def frazier(T):
    if isinstance(T,list):
        return [7.66e-6*np.exp(-1.78/(kb*t)) for t in T]
    else:
        return 7.66e-6*np.exp(-1.78/(kb*T))

def mobility_calc(m0,q,T):
    if isinstance(T, list):
        return [m0*np.exp(-q/(kb*t)) for t in T]
    else:
        return m0*np.exp(-q/(kb*T))

iltis = gg_calc(2.5,33,5,1)
iltis_T_converted = 1/(kb*(1000+273))
prabhakaran = [gg_calc(125,159,4), gg_calc(79,111,4), gg_calc(66,74,4), gg_calc(44,55,4),
               gg_calc(111,125,48), gg_calc(74,100,48), gg_calc(55,100,48)]
prabhakaran_T_converted = [1/(kb*(1210+273)), 1/(kb*(1210+273)), 1/(kb*(1210+273)), 1/(kb*(1210+273)),
                           1/(kb*(900+273)), 1/(kb*(900+273)), 1/(kb*(900+273))]

frazier_T = [273+i for i in range(700,901)]
frazier_T_converted = [1/(kb*T) for T in frazier_T]

T_range_EAM = list(range(1050,1401))
T_range_EAM_converted = [1/(kb*T) for T in T_range_EAM]
meam_100_20 = mobility_calc(667.3e-9,0.27,T_range_EAM)
meam_100_30 = mobility_calc(1.7e-6,0.39,T_range_EAM)
meam_100_45 = mobility_calc(3.1e-6,0.47,T_range_EAM)
meam_110_20 = mobility_calc(279.6e-9,0.17,T_range_EAM)
meam_110_30 = mobility_calc(616e-9,0.24,T_range_EAM)
meam_110_45 = mobility_calc(421.3e-9,0.21,T_range_EAM)
meam_111_20 = mobility_calc(1.1e-6,0.33,T_range_EAM)
meam_111_30 = mobility_calc(438.5e-9,0.22,T_range_EAM)
meam_111_sigma7 = mobility_calc(557.4e-9,0.26,T_range_EAM)
meam_111_45 = mobility_calc(250.5e-9,0.18,T_range_EAM)

eam_100_20 = mobility_calc(421.3e-9,0.31,T_range_EAM)
eam_100_30 = mobility_calc(2.4e-6,0.46,T_range_EAM)
eam_100_45 = mobility_calc(256.1e-6,1.07,T_range_EAM)
eam_110_20 = mobility_calc(848.3e-9,0.29,T_range_EAM)
eam_110_30 = mobility_calc(3.2e-6,0.48,T_range_EAM)
eam_110_45 = mobility_calc(2e-6,0.42,T_range_EAM)
eam_111_20 = mobility_calc(93.3e-6,0.88,T_range_EAM)
eam_111_30 = mobility_calc(236.4e-6,1.05,T_range_EAM)
eam_111_sigma7 = mobility_calc(21.9e-6,0.73,T_range_EAM)
eam_111_45 = mobility_calc(109.5e-6,0.93,T_range_EAM)

T_range_ADP = list(range(1050,1301))
T_range_ADP_converted = [1/(kb*T) for T in T_range_ADP]
adp_100_20 = mobility_calc(88.5e-9,0.13,T_range_ADP)
adp_100_30 = mobility_calc(156.5e-9,0.12,T_range_ADP)
adp_100_45 = mobility_calc(5.6e-6,0.51,T_range_ADP)
adp_110_20 = mobility_calc(530.2e-9,0.22,T_range_ADP)
adp_110_30 = mobility_calc(187.4e-9,0.14,T_range_ADP)
adp_110_45 = mobility_calc(140.2e-9,0.11,T_range_ADP)
adp_111_20 = mobility_calc(282.4e-9,0.16,T_range_ADP)
adp_111_30 = mobility_calc(1.8e-6,0.37,T_range_ADP)
adp_111_sigma7 = mobility_calc(288.1e-9,0.19,T_range_ADP)
adp_111_45 = mobility_calc(966.1e-9,0.33,T_range_ADP)

T = T_range_EAM_converted
plt.semilogy(T,meam_100_20, color = 'red', linestyle = linestyles['solid'], label = "MEAM 100 20$^\circ$")
plt.semilogy(T,meam_100_30, color = 'red', linestyle = linestyles['loosely dotted'], label = "MEAM 100 30$^\circ$")
plt.semilogy(T,meam_100_45, color = 'red', linestyle = linestyles['densely dotted'], label = "MEAM 100 45$^\circ$")
plt.semilogy(T,meam_110_20, color = 'red', linestyle = linestyles['loosely dashed'], label = "MEAM 110 20$^\circ$")
plt.semilogy(T,meam_110_30, color = 'red', linestyle = linestyles['dashed'], label = "MEAM 110 30$^\circ$")
plt.semilogy(T,meam_110_45, color = 'red', linestyle = linestyles['densely dashed'], label = "MEAM 110 45$^\circ$")
plt.semilogy(T,meam_111_20, color = 'red', linestyle = linestyles['loosely dashdotted'], label = "MEAM 111 20$^\circ$")
plt.semilogy(T,meam_111_30, color = 'red', linestyle = linestyles['dashdotted'], label = "MEAM 111 30$^\circ$")
plt.semilogy(T,meam_111_sigma7, color = 'red', linestyle = linestyles['densely dashdotted'], label = "MEAM 111 38.20$^\circ$")
plt.semilogy(T,meam_111_45, color = 'red', linestyle = linestyles['dashdotdotted'], label = "MEAM 111 45$^\circ$")

plt.semilogy(T,eam_100_20, color = 'blue', linestyle = linestyles['solid'], label = "EAM 100 20$^\circ$")
plt.semilogy(T,eam_100_30, color = 'blue', linestyle = linestyles['loosely dotted'], label = "EAM 100 30$^\circ$")
plt.semilogy(T,eam_100_45, color = 'blue', linestyle = linestyles['densely dotted'], label = "EAM 100 45$^\circ$")
plt.semilogy(T,eam_110_20, color = 'blue', linestyle = linestyles['loosely dashed'], label = "EAM 110 20$^\circ$")
plt.semilogy(T,eam_110_30, color = 'blue', linestyle = linestyles['dashed'], label = "EAM 110 30$^\circ$")
plt.semilogy(T,eam_110_45, color = 'blue', linestyle = linestyles['densely dashed'], label = "EAM 110 45$^\circ$")
plt.semilogy(T,eam_111_20, color = 'blue', linestyle = linestyles['loosely dashdotted'], label = "EAM 111 20$^\circ$")
plt.semilogy(T,eam_111_30, color = 'blue', linestyle = linestyles['dashdotted'], label = "EAM 111 30$^\circ$")
plt.semilogy(T,eam_111_sigma7, color = 'blue', linestyle = linestyles['densely dashdotted'], label = "EAM 111 38.20$^\circ$")
plt.semilogy(T,eam_111_45, color = 'blue', linestyle = linestyles['dashdotdotted'], label = "EAM 111 45$^\circ$")

T = T_range_ADP_converted
plt.semilogy(T,adp_100_20, color = 'green', linestyle = linestyles['solid'], label = "ADP 100 20$^\circ$")
plt.semilogy(T,adp_100_30, color = 'green', linestyle = linestyles['loosely dotted'], label = "ADP 100 30$^\circ$")
plt.semilogy(T,adp_100_45, color = 'green', linestyle = linestyles['densely dotted'], label = "ADP 100 45$^\circ$")
plt.semilogy(T,adp_110_20, color = 'green', linestyle = linestyles['loosely dashed'], label = "ADP 110 20$^\circ$")
plt.semilogy(T,adp_110_30, color = 'green', linestyle = linestyles['dashed'], label = "ADP 110 30$^\circ$")
plt.semilogy(T,adp_110_45, color = 'green', linestyle = linestyles['densely dashed'], label = "ADP 110 45$^\circ$")
plt.semilogy(T,adp_111_20, color = 'green', linestyle = linestyles['loosely dashdotted'], label = "ADP 111 20$^\circ$")
plt.semilogy(T,adp_111_30, color = 'green', linestyle = linestyles['dashdotted'], label = "ADP 111 30$^\circ$")
plt.semilogy(T,adp_111_sigma7, color = 'green', linestyle = linestyles['densely dashdotted'], label = "ADP 111 38.20$^\circ$")
plt.semilogy(T,adp_111_45, color = 'green', linestyle = linestyles['dashdotdotted'], label = "ADP 111 45$^\circ$")

plt.semilogy(frazier_T_converted, frazier(frazier_T), color = 'black', linestyle = linestyles['dotted'], label = "Frazier (2018)")
plt.semilogy(iltis_T_converted,iltis, 'k^', label = "Iltis (2017)")
plt.semilogy(prabhakaran_T_converted, prabhakaran, 'ko', label = "Prabhakaran (2019)")
# plt.legend()
plt.savefig("mobility_comparisons.png")

plt.show()
