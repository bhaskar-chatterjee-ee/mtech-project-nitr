""" Created on Mon Mar 13 17:12:52 2023 @author: chatt """

import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt

from NetworkData import *
from load_flow_combined import *
from short_circuit_analysis import *


# =============================================================================
def form_msrmnt_vector(net):
    net.no_V_msrmnt = len(net.pmu_location)         #No of voltage measurements possible with the given PMUs
    
    a = []
    for i in range(len(net.pmu_location)):
        for j in range(net.no_fd):
            if net.fd[j].end_bus[0] == net.pmu_location[i] or net.fd[j].end_bus[1] == net.pmu_location[i]:
                a.append(j)
                
    for i in range(len(net.pmu_location)):
        for j in range(net.no_tr):
            if net.tr[j].end_bus[0] == net.pmu_location[i] or net.tr[j].end_bus[1] == net.pmu_location[i]:
                a.append(j)
    # print(a)            
    net.no_I_msrmnt = len(a)                        #No of current measurements possible with the given PMUs
    
    net.no_msrmnt = net.no_V_msrmnt + net.no_I_msrmnt
    
    if net.outage_line != None:
        a = copy.deepcopy(net.fd[net.outage_line].I[0])
        b = copy.deepcopy(net.fd[net.outage_line].I[1])
        net.fd[net.outage_line].I[0] = 0
        net.fd[net.outage_line].I[1] = 0
    # print(net.fd[net.outage_line].I[0], net.fd[net.outage_line].I[1])
    
    for i in range(len(net.pmu_location)):
        m = Measurement()
        j = net.pmu_location[i] 
        m.type = 1      #voltage measurement on buses
        m.msringV_onbus = j
        m.reading = net.bus[j].V
        net.add_msrmnt(m)
        del m
    
    for i in range(len(net.pmu_location)):
        j = net.pmu_location[i]
        for k in range(net.no_fd):
            if net.fd[k].end_bus[0] == j:
                m = Measurement()                
                m.type = 2      #current measurement through feeders
                m.msringI_onfd = k
                m.reading = net.fd[k].I[0]
                # print(m.reading)
                net.add_msrmnt(m)
                del m
            elif net.fd[k].end_bus[1] == j:
                m = Measurement()                
                m.type = 2      #current measurement through feeders
                m.msringI_onfd = k
                m.reading = net.fd[k].I[1]
                # print(m.reading)
                net.add_msrmnt(m)
                del m
    
    for i in range(len(net.pmu_location)):
        j = net.pmu_location[i]
        for k in range(net.no_tr):
            if net.tr[k].end_bus[0] == j:
                m = Measurement()
                m.type = 3      #current measurement through transformers
                m.msringI_ontr = k
                m.reading = net.tr[k].I[0]
                # print(m.reading)
                net.add_msrmnt(m)
                del m
            elif net.tr[k].end_bus[1] == j:
                m = Measurement()
                m.type = 3      #current measurement through transformers
                m.msringI_ontr = k
                m.reading = net.tr[k].I[1]
                # print(m.reading)
                net.add_msrmnt(m)
                del m
    
    if net.outage_line != None:
        net.fd[net.outage_line].I[0] = a
        net.fd[net.outage_line].I[0] = b           
    
# =============================================================================


# =============================================================================    
def StateEstimation(net):
    
    form_msrmnt_vector(net)
    
    msrmnt_Vector = np.zeros((net.no_msrmnt), dtype = complex)
    for i in range(net.no_msrmnt):
        msrmnt_Vector[i] = net.msrmnt[i].reading   
    print(f'\n\n\nThe measurement vector is:\n {msrmnt_Vector}\n')
        
    #Measurement Bus Incidence Matrix
    A = np.zeros((net.no_msrmnt, net.no_bus), dtype = complex)
    
    for i in range(net.no_msrmnt):
        if net.msrmnt[i].type == 1:             #voltage measurement on buses
            j = net.msrmnt[i].msringV_onbus
            A[i,j] = 1
        elif net.msrmnt[i].type == 2:           #current measurement through feeders
            k = net.msrmnt[i].msringI_onfd
            if net.fd[k].I[0] == net.msrmnt[i].reading:
                m = net.fd[k].end_bus[0]
                n = net.fd[k].end_bus[1]
                A[i,m] = A[i,m] + net.fd[k].y + net.fd[k].ys/2
                A[i,n] = A[i,n] - net.fd[k].y
            elif net.fd[k].I[1] == net.msrmnt[i].reading:
                m = net.fd[k].end_bus[0]
                n = net.fd[k].end_bus[1]
                A[i,m] = A[i,m] - net.fd[k].y
                A[i,n] = A[i,n] + net.fd[k].y + net.fd[k].ys/2       
        elif net.msrmnt[i].type == 3:           #current measurement through transformers
            k = net.msrmnt[i].msringI_ontr
            if net.tr[k].I[0] == net.msrmnt[i].reading:
                m = net.tr[k].end_bus[0]
                n = net.tr[k].end_bus[1]
                A[i,m] = A[i,m] + np.conj(net.tr[k].alpha) * net.tr[k].alpha * net.tr[k].y
                A[i,n] = A[i,n] - np.conj(net.tr[k].alpha) * net.tr[k].y
            elif net.tr[k].I[1] == net.msrmnt[i].reading:
                m = net.tr[k].end_bus[0]
                n = net.tr[k].end_bus[1]
                A[i,m] = A[i,m] - np.conj(net.tr[k].alpha) * net.tr[k].y
                A[i,n] = A[i,n] + np.conj(net.tr[k].alpha) * net.tr[k].alpha * net.tr[k].y
                
            
    print(f'\n\nThe measurement bus incidence matrix:\n {A}\n')
            
    net.state_Vector = np.zeros((net.no_bus), dtype = complex)
    net.state_Vector = np.linalg.inv(A.transpose() @ A) @ A.transpose() @ msrmnt_Vector
    # print(f'\n\nThe system state of IEEE {net.no_bus} Bus system obtained after state estimation is:\n {net.state_Vector}\n')
        
    print_state_estimation_result(net)
    
    error_Vector = np.zeros((net.no_bus), dtype = complex)
    for i in range(net.no_bus):
        error_Vector[i] = net.bus[i].V - net.state_Vector[i]
    print(f'\n\nError in State Estimation:\n {error_Vector}\n')
    print(f'\n\nError squared is:\n {error_Vector.transpose() @ error_Vector}\n')
    
# =============================================================================

# =============================================================================
def print_state_estimation_result(net):
    bus = ["Bus"]
    Vm = ["Vm (pu)"]
    delta = ["delta (deg)"]
    delta_rad = ["delta (rad)"]
    
    for i in range(net.no_bus):
        bus.append(str(net.bus[i].actual_code))
        Vm.append(np.abs(net.state_Vector[i]))
        delta.append(np.angle(net.state_Vector[i]) * 180/np.pi)
        delta_rad.append(np.angle(net.state_Vector[i]))
    
    data = np.transpose([bus, Vm, delta_rad, delta])
    
    print("\n---------------System State after performing State Estimation---------------\n")
    print(tabulate(data,  headers="firstrow",numalign="right", floatfmt=".5f", tablefmt="fancy_grid"))                    
# =============================================================================

