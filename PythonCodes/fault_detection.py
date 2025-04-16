""" Created on Wed May  3 10:53:35 2023 @author: chatt """

import pandas as pd
import numpy as np
import copy
import math

from NetworkData import *
from load_flow_combined import *
from short_circuit_analysis import *

#==============================================================================
def fault_detection(net, F, A):
    if F == 0:  #i.e., fault ocuured on bus
        V_0 = np.zeros((net.no_bus), dtype = complex)
        for i in range(net.no_bus):
            V_0[i] = copy.deepcopy(net.bus[i].V)
        # print(V_0)
        
        global fault_bus
        fault_bus = A
        system_bus_matrix(net)
        mod_Y = copy.deepcopy(net.Y)
        
        for i in range(net.no_bus):
            if net.bus[i].load_status == 1:
                net.bus[i].yload = net.bus[i].load_I / net.bus[i].V 
                
        for i in range(net.no_bus): 
            if net.bus[i].load_status == 1:
                mod_Y[i,i] = mod_Y[i,i] + net.bus[i].yload
                
        for k in range(net.no_gen):
            i = net.gen[k].bus
            mod_Y[i,i] = mod_Y[i,i] + net.gen[k].y
            
        net.Z = np.linalg.inv(mod_Y)
        
        for i in range(net.no_gen):
            j = net.gen[i].bus
            #Calculation of prefault generator emf
            net.gen[i].E = net.bus[j].V + (net.gen[i].Ra + 1j * net.gen[i].Xd_d) * net.bus[j].gen_I
        
        mod_I = np.zeros((net.no_bus), dtype = complex)
        
        for i in range(net.no_gen):
            k = net.gen[i].bus
            mod_I[k] = net.gen[i].E * net.gen[i].y
            
        mod_Y[fault_bus,fault_bus] = 10e6
        
        V_f = np.linalg.inv(mod_Y) @ mod_I
        # print(V_f)
        
        delta_V = V_f - V_0
        # print(f'\n {delta_V}')
        
        I_f = V_0[fault_bus] /(net.Z_f + net.Z[fault_bus,fault_bus])
        # print(f'\n {I_f}')
            
        
    elif F == 1:
        fault_line = A
        fault_dist = 0.5
        net.no_bus = net.no_bus + 1
        net.no_fd = net.no_fd + 1
        
        b = Bus()
        b.code = net.no_bus - 1
        b.zone = None
        b.actual_code = net.no_bus
        b.Pg = 0
        b.Pd = 0    #i.e., fault ocuured on line
        b.Qg = 0
        b.Qd = 0
        b.type = 3
        b.Vm = 0.0
        b.delta = 0.0
        b.V = 0.0 + 1j * 0.0
        b.P = 0.0
        b.Q = 0.0
        net.add_bus(b)
        del b
        
        fd = Feeder()
        fd.code = net.no_fd
        fd.end_bus[0] = net.no_bus - 1
        fd.end_bus[1] = net.fd[fault_line].end_bus[1]
        fd.r = (1 - fault_dist) * net.fd[fault_line].r
        fd.x = (1 - fault_dist) * net.fd[fault_line].x
        fd.z = fd.r + 1j * fd.x
        fd.y = 1/fd.z
        fd.ys = 0.5 * net.fd[fault_line].ys 
        fd.reac_frombus = 0
        fd.reac_tobus = 0
        fd.status = 1
        fd.I = [0.0 + 1j * 0.0, 0.0 + 1j * 0.0]
        fd.S = [0.0 + 1j * 0.0, 0.0 + 1j * 0.0]
        net.add_feeder(fd)
        del fd
        
        for i in range(net.no_fd):
            if net.fd[i].code == fault_line + 1:
                net.fd[i].end_bus[0] = net.fd[i].end_bus[0]
                net.fd[i].end_bus[1] = net.no_bus - 1
                net.fd[i].r = fault_dist * net.fd[i].r
                net.fd[i].x = fault_dist * net.fd[i].x
                net.fd[i].z = net.fd[i].r + 1j * net.fd[i].x
                net.fd[i].y = 1/net.fd[i].z
                net.fd[i].ys = 0.5 * net.fd[fault_line].ys
        
        system_bus_matrix(net)
        FDLF(net)
        
        V_0 = np.zeros((net.no_bus), dtype = complex)
        for i in range(net.no_bus):
            V_0[i] = copy.deepcopy(net.bus[i].V)
        
        mod_Y = copy.deepcopy(net.Y)
        
        for i in range(net.no_bus):
            if net.bus[i].load_status == 1:
                net.bus[i].yload = net.bus[i].load_I / net.bus[i].V 
                
        for i in range(net.no_bus):
            if net.bus[i].load_status == 1:
                mod_Y[i,i] = mod_Y[i,i] + net.bus[i].yload
                
        for k in range(net.no_gen):
            i = net.gen[k].bus
            mod_Y[i,i] = mod_Y[i,i] + net.gen[k].y
            
        net.Z = np.linalg.inv(mod_Y)
            
        for i in range(net.no_gen):
            j = net.gen[i].bus
            #Calculation of prefault generator emf
            net.gen[i].E = net.bus[j].V + (net.gen[i].Ra + 1j * net.gen[i].Xd_d) * net.bus[j].gen_I
        
        mod_I = np.zeros((net.no_bus), dtype = complex)
        
        for i in range(net.no_gen):
            k = net.gen[i].bus
            mod_I[k] = net.gen[i].E * net.gen[i].y
        
        fault_bus = net.no_bus - 1
            
        mod_Y[fault_bus,fault_bus] = 10e6
        
        V_f = np.linalg.inv(mod_Y) @ mod_I
        
        delta_V = V_f - V_0
        # print(f'\n {delta_V}')
        
        I_f = V_0[fault_bus] /(net.Z_f + net.Z[fault_bus,fault_bus])
        # print(f'\n {I_f}')
        
    m = len(net.pmu_location)   #no of PMUs in the network
    
    KD = np.zeros(((net.no_bus),(m)),dtype = float)
    for i in range(net.no_bus):
        for j in range(m):
            k = net.pmu_location[j]
            KD[i,j] = np.abs(delta_V[k]/ net.Z[k,i])
    
    # print(KD)
    
    K_bar = np.zeros((net.no_bus),dtype = float)
    

    for i in range(net.no_bus):
        for j in range(m):
            K_bar[i] = K_bar[i] + KD[i,j]
    
    K_bar = K_bar/m
    # print(K_bar)
    
    delta_k = np.zeros((net.no_bus),dtype = float)       
    for i in range(net.no_bus):
        for j in range(m):
            delta_k[i] = delta_k[i] + (KD[i,j] - K_bar[i])**2
        delta_k[i] = math.sqrt(delta_k[i]/m)
        # print(delta_k[i])
    print(delta_k)
    
    temp = np.min(delta_k)
    for i in range(net.no_bus):
        if temp == delta_k[i]:
            print(f"\nfault has occured on bus {i}.")
        
    
            
        
        
     





























###############################################################################

    