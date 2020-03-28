# c(z)=c_0-∆c tanh⁡〖((z-z_0)/σ),〗 where c_0=1490 m/s, z_0=25 m, ∆c=30 m/s, σ=10 m
# c(z)=1490-30*tanh((z-25)/10)

import numpy as np 
import math 

depths = [i for i in range(0,51) if i % 5 == 0]
print(depths)
print(len(depths))
reduced_depths = [(x-25)/10 for x in depths]
print(reduced_depths)
tanh_val = np.tanh(reduced_depths)
print(tanh_val)
c=[round(1490-30*x,2) for x in tanh_val]
print(c)