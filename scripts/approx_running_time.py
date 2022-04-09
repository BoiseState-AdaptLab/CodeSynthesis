from cProfile import label
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

from sklearn.metrics import r2_score


# Data created from runs of bench_harness on files created using make_fake_mtx_files.cpp running COO -> CSR naive
nnzs = []
microseconds_run_times = []
nnzs.append(10000) 
microseconds_run_times.append(59901883)
nnzs.append(1000) 
microseconds_run_times.append(526899)
nnzs.append(100) 
microseconds_run_times.append(4967)
nnzs.append(1100) 
microseconds_run_times.append(630081)
nnzs.append(1200) 
microseconds_run_times.append(752015)
nnzs.append(1300) 
microseconds_run_times.append(881128)
nnzs.append(1400) 
microseconds_run_times.append(1042300)
nnzs.append(1500) 
microseconds_run_times.append(1193379)
nnzs.append(1600) 
microseconds_run_times.append(1365543)
nnzs.append(1700) 
microseconds_run_times.append(1550952)
nnzs.append(1800) 
microseconds_run_times.append(1748754)
nnzs.append(1900) 
microseconds_run_times.append(1947850)
nnzs.append(2000) 
microseconds_run_times.append(2161284)
nnzs.append(200) 
microseconds_run_times.append(19864)
nnzs.append(2100) 
microseconds_run_times.append(2395371)
nnzs.append(2200) 
microseconds_run_times.append(2620501)
nnzs.append(2300) 
microseconds_run_times.append(2863369)
nnzs.append(2400) 
microseconds_run_times.append(3117081)
nnzs.append(2500) 
microseconds_run_times.append(3419000)
nnzs.append(2600) 
microseconds_run_times.append(3717010)
nnzs.append(2700) 
microseconds_run_times.append(4009391)
nnzs.append(2800) 
microseconds_run_times.append(4335064)
nnzs.append(2900) 
microseconds_run_times.append(4695769)
nnzs.append(3000) 
microseconds_run_times.append(5059554)
nnzs.append(300) 
microseconds_run_times.append(45357)
nnzs.append(3100) 
microseconds_run_times.append(5360481)
nnzs.append(3200) 
microseconds_run_times.append(5728426)
nnzs.append(3300) 
microseconds_run_times.append(6097572)
nnzs.append(3400) 
microseconds_run_times.append(6546016)
nnzs.append(3500) 
microseconds_run_times.append(6906683)
nnzs.append(3600) 
microseconds_run_times.append(7424588)
nnzs.append(3700) 
microseconds_run_times.append(7737168)
nnzs.append(3800) 
microseconds_run_times.append(8265198)
nnzs.append(3900) 
microseconds_run_times.append(8658820)
nnzs.append(4000) 
microseconds_run_times.append(9172049)
nnzs.append(400) 
microseconds_run_times.append(83758)
nnzs.append(4100) 
microseconds_run_times.append(9533207)
nnzs.append(4200) 
microseconds_run_times.append(10099795)
nnzs.append(4300) 
microseconds_run_times.append(10515877)
nnzs.append(4400) 
microseconds_run_times.append(11585258)
nnzs.append(4500) 
microseconds_run_times.append(11405829)
nnzs.append(4600) 
microseconds_run_times.append(12479830)
nnzs.append(4700) 
microseconds_run_times.append(12411264)
nnzs.append(4800) 
microseconds_run_times.append(13576534)
nnzs.append(4900) 
microseconds_run_times.append(13559601)
nnzs.append(5000) 
microseconds_run_times.append(14710133)
nnzs.append(500) 
microseconds_run_times.append(129652)
nnzs.append(5100) 
microseconds_run_times.append(15237089)
nnzs.append(5200) 
microseconds_run_times.append(15365133)
nnzs.append(5300) 
microseconds_run_times.append(16622591)
nnzs.append(5400) 
microseconds_run_times.append(17080541)
nnzs.append(5500) 
microseconds_run_times.append(17900086)
nnzs.append(5600) 
microseconds_run_times.append(17953819)
nnzs.append(5700) 
microseconds_run_times.append(19222316)
nnzs.append(5800) 
microseconds_run_times.append(19431468)
nnzs.append(5900) 
microseconds_run_times.append(21743282)
nnzs.append(6000) 
microseconds_run_times.append(21723426)
nnzs.append(600) 
microseconds_run_times.append(184770)
nnzs.append(6100) 
microseconds_run_times.append(22122286)
nnzs.append(6200) 
microseconds_run_times.append(22324288)
nnzs.append(6300) 
microseconds_run_times.append(24598974)
nnzs.append(6400) 
microseconds_run_times.append(23600027)
nnzs.append(6500) 
microseconds_run_times.append(25114384)
nnzs.append(6600) 
microseconds_run_times.append(25468923)
nnzs.append(6700) 
microseconds_run_times.append(27272187)
nnzs.append(6800) 
microseconds_run_times.append(26882361)
nnzs.append(6900) 
microseconds_run_times.append(28966554)
nnzs.append(7000) 
microseconds_run_times.append(29136779)
nnzs.append(700) 
microseconds_run_times.append(249334)
nnzs.append(7100) 
microseconds_run_times.append(29449688)
nnzs.append(7200) 
microseconds_run_times.append(31602943)
nnzs.append(7300) 
microseconds_run_times.append(31131206)
nnzs.append(7400) 
microseconds_run_times.append(33174313)
nnzs.append(7500) 
microseconds_run_times.append(33020227)
nnzs.append(7600) 
microseconds_run_times.append(34982006)
nnzs.append(7700) 
microseconds_run_times.append(34828395)
nnzs.append(7800) 
microseconds_run_times.append(38537058)
nnzs.append(7900) 
microseconds_run_times.append(36388263)
nnzs.append(8000) 
microseconds_run_times.append(39384003)
nnzs.append(800) 
microseconds_run_times.append(328540)
nnzs.append(8100) 
microseconds_run_times.append(40271621)
nnzs.append(8200) 
microseconds_run_times.append(39515530)
nnzs.append(8300) 
microseconds_run_times.append(42189123)
nnzs.append(8400) 
microseconds_run_times.append(41719063)
nnzs.append(8500) 
microseconds_run_times.append(44131074)
nnzs.append(8600) 
microseconds_run_times.append(43729260)
nnzs.append(8700) 
microseconds_run_times.append(46457121)
nnzs.append(8800) 
microseconds_run_times.append(45760987)
nnzs.append(8900) 
microseconds_run_times.append(48376912)
nnzs.append(9000) 
microseconds_run_times.append(48388860)
nnzs.append(900) 
microseconds_run_times.append(423030)
nnzs.append(9100) 
microseconds_run_times.append(49981434)
nnzs.append(9200) 
microseconds_run_times.append(49856072)
nnzs.append(9300) 
microseconds_run_times.append(52774376)
nnzs.append(9400) 
microseconds_run_times.append(52555247)
nnzs.append(9500) 
microseconds_run_times.append(55452615)
nnzs.append(9600) 
microseconds_run_times.append(53817654)
nnzs.append(9700) 
microseconds_run_times.append(55104703)
nnzs.append(9800) 
microseconds_run_times.append(61069560)
nnzs.append(9900) 
microseconds_run_times.append(57343813)


filename_to_nnz = {
"atmosmodd.mtx" : 8814880,
"Baumann.mtx" : 760631,
"cant.mtx" : 2034917,
"chem_master1.mtx" : 201201,
"consph.mtx" : 3046907,
"cop20k_A.mtx" : 1362087,
"denormal.mtx" : 622812,
"dixmaanl.mtx" : 179999,
"ecology1.mtx" : 2998000,
"jnlbrng1.mtx" : 119600,
"Lin.mtx" : 1011200,
"mac_econ_fwd500.mtx" : 1273389,
"majorbasis.mtx" : 1750416,
"obstclae.mtx" : 118804,
"pdb1HYS.mtx" : 2190591,
"pwtk.mtx" : 5926171,
"rma10.mtx" : 2374001,
"scircuit.mtx" : 958936,
"shipsec1.mtx" : 3977139,
"shyy161.mtx" : 329762,
"webbase_1M.mtx" : 3105536,
}

nnzs = np.array(nnzs)
microsecond_run_times = np.array(microseconds_run_times)

# microsecond to second
microsecond_run_times =  microsecond_run_times / 1_000_000

# second to minute
microsecond_run_times =  microsecond_run_times / 60

# minute to hour
microsecond_run_times =  microsecond_run_times / 60

x_for_plot = np.linspace(100, 100_000, 100)


model = Polynomial.fit(nnzs, microsecond_run_times, 2)
plt.scatter(nnzs, microsecond_run_times, c="r", label="data")
plt.plot(x_for_plot, model(x_for_plot), label="model")
plt.xlabel("hours running time")
plt.ylabel("number non zeros")

print(f"r2: {r2_score(microsecond_run_times, model(nnzs))}")
for fn, nnz in filename_to_nnz.items():
    print(f"{fn}, nnz: {nnz}, Approx. running time (hrs) {int(model(nnz))}")


plt.legend()
plt.show()


