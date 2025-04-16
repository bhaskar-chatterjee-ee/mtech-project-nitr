""" Created on Fri Apr 14 19:59:55 2023 @author: chatt """

import numpy as np
import math

# IEEE 3 Bus System.
Z = np.array([[ 0.0027696 +0.04556139j, -0.00020388+0.02940507j, 0.00288494+0.03286259j],
              [-0.00020388+0.02940507j,  0.00737818+0.06012434j, 0.00685829+0.04439028j],
              [ 0.00288494+0.03286259j,  0.00685829+0.04439028j, 0.03910555+0.09174853j]])

V_0 = np.array([1+0j, 1.019988055892558-0.004936176305600952j, 0.9353015881671138-0.060116488813617856j ])

# Fault at Bus-2.
fault_bus = 1
V_f = np.array([5.08664088e-01-6.13398315e-02j, 1.97003599e-07-1.67228975e-06j, 1.79456679e-01-3.28784714e-02j])
I_f = 1.9700544108738791-16.722888282454033j

delta_I = [0, 1.9700544108738791-16.722888282454033j, 0 ]

delta_V = V_f - V_0
print(delta_V)
# print(f'\n {Z @ delta_I}\n\n')

# Two PMUs are installed at bus 1 and 3.
pmu_location = np.array([0,2])
a = len(pmu_location)

K_D = np.zeros((a, a), dtype = complex)
for i in range(a):
    K_D[i] = delta_V[pmu_location[i]]/Z[pmu_location[i],fault_bus]

# print(K_D)

K = (K_D[0] + K_D[1])/a
print(K)

for i in range(a):   
    k = (K_D[i]- K).transpose()*(K_D[i]- K)
    
MD = math.sqrt((np.real(k[0])+np.real(k[1]))/a) 
print(f'\n\n {MD}')

