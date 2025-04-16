import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt

from NetworkData import *
from load_flow_combined import *
#from short_circuit_analysis import *
#from step_distance_protection import *
#from state_estimation import *
#from fault_detection import *
#from load_encroachment import *


plt.close("all")
print("\033[H\033[J") 

# =============================================================================
bus = 9
net = Network()
data_file = 'IEEE' + str(bus) + ".xlsx"
read_network_data(net,data_file)
# =============================================================================

system_bus_matrix(net)

# gauss_siedal(net)
FDLF(net)
# alternative_FDLF(net)

# line_outage(net, 1)
# load_change(net)
# FDLF(net)


# net.pmu_location = np.array([0, 2])           #Bus locations of the PMUs in IEEE-3 bus system
# net.pmu_location = np.array([3, 6, 8])   #Bus locations of the PMUs in IEEE-9 bus system
# StateEstimation(net)


# net.Z_f = 0
# fault_on_bus(net, 1)
# fault_on_line(net, 2, 0.5)
# StateEstimation(net)

# net.outage_line = 1
# StateEstimation(net)

# KCL_mismatch(net)
# print_feeder_result(net)

#RelaySensing(net)
#MhoRelayCharacteristics(net)

# fault_detection(net, 0, 1)  #put 0 on the second argument if fault has occured on bus and 1 if the fault has occured on line.


