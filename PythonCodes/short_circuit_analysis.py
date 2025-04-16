import pandas as pd
import numpy as np
import copy as copy
import math

from NetworkData import *
from load_flow_combined import *

# =============================================================================
def modified_Y_bus(net):
    net.mod_Y = copy.deepcopy(net.Y)
    
    for i in range(net.no_bus):
        if net.bus[i].load_status == 1:
            net.bus[i].yload = net.bus[i].load_I / net.bus[i].V 
    
    for i in range(net.no_bus):
        if net.bus[i].load_status == 1:
            net.mod_Y[i,i] = net.mod_Y[i,i] + net.bus[i].yload  
    
    for k in range(net.no_gen):
        i = net.gen[k].bus
        net.mod_Y[i,i] = net.mod_Y[i,i] + net.gen[k].y
    
    # print(f"\n\n The modified Admittance Matrix:\n{net.mod_Y}\n\n")
    
    net.Z = np.linalg.inv(net.mod_Y)
    # print(f"\n\n The Impedance matrix:\n {net.Z}\n")
    
# =============================================================================  
 

# =============================================================================
def fault_on_bus(net, fault_bus):
    
    V_prefault = copy.deepcopy(net.bus[fault_bus].V) 
    
    for i in range(net.no_gen):
        j = net.gen[i].bus
        #Calculation of prefault generator emf
        net.gen[i].E = net.bus[j].V + (net.gen[i].Ra + 1j * net.gen[i].Xd_d) * net.bus[j].gen_I
              
    modified_Y_bus(net)
    
    mod_I = np.zeros((net.no_bus), dtype = complex)
    
    for i in range(net.no_gen):
        k = net.gen[i].bus
        mod_I[k] = net.gen[i].E * net.gen[i].y
    
    # print(f"\nModified current Injection Vector = {mod_I}\n")
    
    net.mod_Y[fault_bus,fault_bus] = 10e6
    
    V_f = np.linalg.inv(net.mod_Y) @ mod_I
    # print(f'\nThe post fault voltages  when fault occurs at bus {fault_bus} are: \n{V_f}\n')
    
    # print(f'The post fault voltages of the IEEE {net.no_bus} bus system after fault at bus {fault_bus + 1} are:\n')
    for i in range(net.no_bus):
        net.bus[i].V = V_f[i]
        net.bus[i].Vm = np.abs(net.bus[i].V)
        net.bus[i].delta = np.angle(net.bus[i].V)
        # print(f'{net.bus[i].Vm},{net.bus[i].delta}')
        
    calculate_P(net)
    calculate_Q(net)
    
    for i in range(net.no_bus):
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)
        
    calculate_PgandQg(net)
        
    current_injection(net)
        
    print_postfault_result(net)
    
    # print(f'Prefault Voltage at {fault_bus} is {V_prefault}')
    # print(f'\n{net.Z[fault_bus,fault_bus]}')
    I_f = V_prefault /(net.Z_f + net.Z[fault_bus,fault_bus])
    print(f'\nThe fault current(in cartesian) is: {I_f}\n')
    print(f'The fault current(in polar) is: {np.abs(I_f)}mag_{np.angle(I_f)}rad\n')
    print(f'The fault current(in polar) is: {np.abs(I_f)}mag_{math.degrees(np.angle(I_f))}deg\n')
        
    load_currents(net)
    generator_currents(net)
    feeder_currents(net)
    transformer_currents(net)
    
    print_feeder_current(net)
    
    print_feeder_powerflow(net)
        
    KCL_mismatch(net)
    
# =============================================================================


# =============================================================================
def fault_on_line(net, fault_line, fault_dist):
    # The variable fault_dist is the distance of fault location from the end_bus[0] of the feeder.
    net.no_bus = net.no_bus + 1
    net.no_fd = net.no_fd + 1
    
    b = Bus()
    b.code = net.no_bus - 1
    b.zone = None
    b.actual_code = net.no_bus
    b.Pg = 0
    b.Pd = 0
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
    # fault_on_bus(net, net.no_bus - 1)
    
    fault_bus = net.no_bus - 1
    V_prefault = copy.deepcopy(net.bus[fault_bus].V) 
    
    for i in range(net.no_gen):
        j = net.gen[i].bus
        #Calculation of prefault generator emf
        net.gen[i].E = net.bus[j].V + (net.gen[i].Ra + 1j * net.gen[i].Xd_d) * net.bus[j].gen_I
              
    modified_Y_bus(net)
    
    mod_I = np.zeros((net.no_bus), dtype = complex)
    
    for i in range(net.no_gen):
        k = net.gen[i].bus
        mod_I[k] = net.gen[i].E * net.gen[i].y
    
    # print(f"\nModified current Injection Vector = {mod_I}\n")
    
    net.mod_Y[fault_bus,fault_bus] = 10e6
    
    V_f = np.linalg.inv(net.mod_Y) @ mod_I
    # print(f'\nThe post fault voltages  when fault occurs at bus {f} are: \n{V_f}\n')
    
    # print(f'The post fault voltages of the IEEE {net.no_bus} bus system after fault at bus {fault_bus + 1} are:\n')
    for i in range(net.no_bus):
        net.bus[i].V = V_f[i]
        net.bus[i].Vm = np.abs(net.bus[i].V)
        net.bus[i].delta = np.angle(net.bus[i].V)
        # print(f'{net.bus[i].Vm},{net.bus[i].delta}')
        
    calculate_P(net)
    calculate_Q(net)
    
    for i in range(net.no_bus):
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)
        
    calculate_PgandQg(net)
        
    current_injection(net)
        
    print_postfault_result(net)
    
    # print(f'Prefault Voltage at {fault_bus} is {V_prefault}')
    # print(f'\n{net.Z[fault_bus,fault_bus]}')
    I_f = V_prefault /(net.Z_f + net.Z[fault_bus,fault_bus])
    print(f'\nThe fault current is: \n{I_f}\n')
    print(f'The fault current magnitude is {np.abs(I_f)}')
    print(f'The fault current magnitude is {np.angle(I_f)}')
        
    load_currents(net)
    generator_currents(net)
    feeder_currents(net)
    transformer_currents(net)
    
    print_feeder_current(net)
    
    print_feeder_powerflow(net)
        
    KCL_mismatch(net)
    
# =============================================================================


# =============================================================================
#------------------------------------------------------------------------------
def load_currents(net):   
    #Load Currents
    for i in range(net.no_bus):
        if net.bus[i].load_status == 1:
            net.bus[i].load_I = net.bus[i].V * net.bus[i].yload

#------------------------------------------------------------------------------
def generator_currents(net):
    #Generator Currents
    for i in range(net.no_gen):
        k = net.gen[i].bus
        if net.gen[i].bus == net.bus[k].code:
            net.bus[i].gen_I = (net.gen[i].E - net.bus[k].V) * net.gen[i].y
    else:   
        net.bus[i].gen_I = 0

#------------------------------------------------------------------------------
def feeder_currents(net):
    #Feeder Currents
    for k in range(net.no_fd):
        i = net.fd[k].end_bus[0]
        j = net.fd[k].end_bus[1]
        
        net.fd[k].I[0] = net.fd[k].y * (net.bus[i].V - net.bus[j].V) + (net.fd[k].ys/2) * net.bus[i].V
        net.fd[k].I[1] = net.fd[k].y * (net.bus[j].V - net.bus[i].V) + (net.fd[k].ys/2) * net.bus[j].V
        
        net.fd[k].S[0] = net.bus[i].V * np.conj(net.fd[k].I[0])
        net.fd[k].S[1] = net.bus[j].V * np.conj(net.fd[k].I[1])
        
        net.fd[k].S_loss = net.fd[k].S[0] +  net.fd[k].S[1]

#------------------------------------------------------------------------------
def transformer_currents(net):
    #Transformer Currents
    for k in range(net.no_tr):
        i = net.tr[k].end_bus[0]
        j = net.tr[k].end_bus[1]
        
        # print(f"Bus {i + 1} voltage is {net.bus[i].V}")
        # print(f"Bus {j + 1} voltage is {net.bus[j].V}")
        
        net.tr[k].I[0] = np.conj(net.tr[k].alpha) * net.tr[k].alpha * net.tr[k].y * net.bus[i].V - np.conj(net.tr[k].alpha) * net.tr[k].y * net.bus[j].V
        net.tr[k].I[1] = - net.tr[k].alpha * net.tr[k].y * net.bus[i].V + net.tr[k].y * net.bus[j].V
        
        net.fd[k].S[0] = net.bus[i].V * np.conj(net.tr[k].I[0])
        net.fd[k].S[1] = net.bus[j].V * np.conj(net.tr[k].I[1])
        
        # print(f'The {k + 1} Transformer current is {net.tr[k].I}')
#------------------------------------------------------------------------------
# =============================================================================


# =============================================================================
def KCL_mismatch(net):
    sum_I = np.zeros((net.no_bus),dtype = complex)
    
    current_injection(net)
    
    for i in range(net.no_bus):
        sum_I[i]= net.bus[i].I
    # print(f'\nInjection currents of all the buses = {sum_I}\n')
        
    for i in range(net.no_bus):
        for j in range (net.no_fd):
            if net.bus[i].code == net.fd[j].end_bus[0]:
                sum_I[i] = sum_I[i] - net.fd[j].I[0]
            elif net.bus[i].code == net.fd[j].end_bus[1]:
                sum_I[i] = sum_I[i] - net.fd[j].I[1]
                
        for k in range(net.no_tr):
            if net.bus[i].code == net.tr[k].end_bus[0]:
                sum_I[i] = sum_I[i] - net.tr[k].I[0]
            elif net.bus[i].code == net.tr[k].end_bus[1]:
                sum_I[i] = sum_I[i] - net.tr[k].I[1]
                
    # print(f'\nsum_I = {sum_I}\n')
           
    temp = np.max(np.abs(sum_I))
    # for i in range(net.no_bus):
    #     print(f"Bus {i + 1} of type {net.bus[i].type} has KCL mismatch of {np.abs(sum_I[i])}")
        
    print(f"\nThe maximum mismatch in KCL checking is {temp}") 
     
# =============================================================================


# =============================================================================
def print_postfault_result(net):
    bus = ["Bus"]
    Vm = ["Vm (pu)"]
    delta = ["delta (deg)"]
    delta_rad = ["delta (rad)"]
    Pg = ["Pg (pu)"]
    Pd = ["Pd (pu)"]
    Qg = ["Qg (pu)"]
    Qd = ["Qd (pu)"]
    for i in range(net.no_bus):
        bus.append(str(net.bus[i].actual_code))
        Vm.append(net.bus[i].Vm)
        delta.append(net.bus[i].delta * 180/np.pi)
        delta_rad.append(net.bus[i].delta)
        Pg.append(net.bus[i].Pg)
        Qg.append(net.bus[i].Qg)
        Pd.append(net.bus[i].Pd)
        Qd.append(net.bus[i].Qd)

    data = np.transpose([bus, Vm, delta_rad, delta, Pg, Qg, Pd, Qd])

    print("\n------------Post Fault System State-------------\n")
    print(tabulate(data,  headers="firstrow",numalign="right", floatfmt=".5f", tablefmt="fancy_grid"))
# =============================================================================    