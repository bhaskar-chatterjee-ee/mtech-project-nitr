import pandas as pd
import numpy as np
import copy
from tabulate import tabulate

from NetworkData import *

# =============================================================================
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(complex(x, y))
# =============================================================================

# =============================================================================
def Y_bus_matrix(net):
    net.Y = np.zeros((net.no_bus, net.no_bus), dtype = complex)
    
    for k in range(net.no_fd):
        i= net.fd[k].end_bus[0]
        j= net.fd[k].end_bus[1]
    
        net.Y[i,i] = net.Y[i,i] + net.fd[k].y + net.fd[k].ys/2
        net.Y[i,j] = net.Y[i,j] - net.fd[k].y 
        net.Y[j,i] = net.Y[j,i] - net.fd[k].y
        net.Y[j,j] = net.Y[j,j] + net.fd[k].y + net.fd[k].ys/2
        
    for k in range(net.no_tr):
        i = net.tr[k].end_bus[0]
        j = net.tr[k].end_bus[1]
        a = net.tr[k].alpha
        
        net.Y[i,i] = net.Y[i,i] + np.conj(a) * a * net.tr[k].y
        net.Y[i,j] = net.Y[i,j] - np.conj(a) * net.tr[k].y
        net.Y[j,i] = net.Y[j,i] - a * net.tr[k].y
        net.Y[j,j] = net.Y[j,j] + net.tr[k].y
        
    net.G = np.real(net.Y)
    net.B = np.imag(net.Y)
# =============================================================================

# =============================================================================
def alter_Y_bus_matrix(net):
    net.alter_Y = np.zeros((net.no_bus,net.no_bus),dtype = complex)
    
    for k in range(net.no_fd):
        i = net.fd[k].end_bus[0]
        j = net.fd[k].end_bus[1]
        
        net.alter_Y[i,i] = net.alter_Y[i,i] + 1 / (1j *net.fd[k].x)
        net.alter_Y[i,j] = net.alter_Y[i,j] - 1 / (1j *net.fd[k].x)
        net.alter_Y[j,i] = net.alter_Y[j,i] - 1 / (1j *net.fd[k].x)
        net.alter_Y[j,j] = net.alter_Y[j,j] + 1 / (1j *net.fd[k].x)
    
    for k in range(net.no_tr):
        i = net.tr[k].end_bus[0]
        j = net.tr[k].end_bus[1]
        a = net.tr[k].alpha
    
        net.alter_Y[i,i] = net.alter_Y[i,i] + (np.conj(pol2cart(1, np.angle(a))))*(pol2cart(1, np.angle(a)))*(1 / (1j * net.tr[k].x))
        net.alter_Y[i,j] = net.alter_Y[i,j] - (np.conj(pol2cart(1, np.angle(a))))*(1 / (1j * net.tr[k].x))
        net.alter_Y[j,i] = net.alter_Y[j,i] - (pol2cart(1, np.angle(a)))*(1 / (1j * net.tr[k].x))
        net.alter_Y[j,j] = net.alter_Y[j,j] + (1 / (1j * net.tr[k].x))
        
# =============================================================================

# =============================================================================
def Bd_matrix(net):
    net.Bd = np.zeros(((net.no_bus - 1),(net.no_bus - 1)),dtype = float)
    
    net.Np = []
    for i in range(net.no_bus):
        if net.bus[i].type == 2 or net.bus[i].type == 3:
            net.Np.append(net.bus[i].code)
    
    z = len(net.Np)
    
    for i in range(z):
        for j in range(z):
            a = net.Np[i]
            b = net.Np[j]
            net.Bd[i,j] = - np.imag(net.alter_Y[a,b])
            
    net.Bd_inv = np.linalg.inv(net.Bd)
                    
# =============================================================================

# =============================================================================
def Bdd_matrix(net):
    net.no_pq = net.no_bus - (net.no_pv + 1)
    net.Bdd = np.zeros(((net.no_pq),(net.no_pq)),dtype = float)
    
    net.Npq = []
    for i in range(net.no_bus):
        if net.bus[i].type == 3:
            net.Npq.append(net.bus[i].code)
     
    for i in range(net.no_pq):
        for j in range(net.no_pq):
            a = net.Npq[i]
            b = net.Npq[j]
            net.Bdd[i,j] = - net.B[a,b]
            
    net.Bdd_inv = np.linalg.inv(net.Bdd)

# =============================================================================

# =============================================================================
def opt_Bd_matrix(net):
    net.Bd = np.zeros(((net.no_bus - 1),(net.no_bus - 1)),dtype = float)
    
    net.Bd = - np.imag(net.alter_Y)
    
    i = net.slack_bus.code
    
    net.Bd[i,i] = 1e10
    
    net.Bd_inv = np.linalg.inv(net.Bd)

                    
# =============================================================================

# =============================================================================
def opt_Bdd_matrix(net):
    net.no_pq = net.no_bus - (net.no_pv + 1)
    net.Bdd = np.zeros(((net.no_pq),(net.no_pq)),dtype = float)
    
    net.Bdd = - net.B
    
    i = net.slack_bus.code
    
    net.Bdd[i,i] = 1e10
    
    for i in range(net.no_bus):
        if net.bus[i].type == 3:
            net.Npq.append(net.bus[i].code)
     
    for i in range(net.no_bus):
        if net.bus[i].type == 2:
            net.Bdd[i,i] = 1e10
            
    net.Bdd_inv = np.linalg.inv(net.Bdd)

# =============================================================================

# =============================================================================    
def system_bus_matrix(net):
    Y_bus_matrix(net)
    alter_Y_bus_matrix(net)
    
    Bd_matrix(net)
    # opt_Bd_matrix(net)
    
    Bdd_matrix(net)
    # opt_Bdd_matrix(net)
# =============================================================================

# =============================================================================
def power_injection(net):
    for i in range(net.no_bus):
        net.bus[i].P = net.bus[i].Pg - net.bus[i].Pd
        net.bus[i].Q = net.bus[i].Qg - net.bus[i].Qd
        
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)

        net.bus[i].Psp = net.bus[i].P
        net.bus[i].Qsp = net.bus[i].Q
        
def current_injection(net):
    for i in range(net.no_bus):
        net.bus[i].I = np.conj(net.bus[i].S / net.bus[i].V)
        # print(f'\n{net.bus[i].I}')
# =============================================================================

# =============================================================================
def bus_initialization(net):
     
    #PQ Bus Initialization
    for i in range(net.no_bus):
        if net.bus[i].type == 3:
            net.bus[i].Vm = 1.0
            net.bus[i].delta = 0.0
            net.bus[i].V = pol2cart(net.bus[i].Vm, net.bus[i].delta)
    

    #PV Bus Initialization
    for i in range(net.no_pv):
        k=net.pv_bus[i].code
        net.bus[k].Vm = net.pv_bus[i].Vsp
        net.bus[k].delta = 0.0
        net.bus[k].V = pol2cart(net.bus[k].Vm, net.bus[k].delta)

    #Slack Bus Initialization
    i = net.slack_bus.code
    net.bus[i].Vm = net.slack_bus.Vsp
    net.bus[i].delta = 0.0
    net.bus[i].V = pol2cart(net.bus[i].Vm, net.bus[i].delta)
    
# =============================================================================

# =============================================================================
def calculate_P(net):
    for i in range(net.no_bus):
        net.bus[i].P = 0.0
        for j in range(net.no_bus):
            net.bus[i].P += (net.bus[i].Vm * net.bus[j].Vm*(net.G[i,j] * np.cos(net.bus[i].delta - net.bus[j].delta) + net.B[i,j] * np.sin(net.bus[i].delta - net.bus[j].delta)))


def calculate_Q(net):
    for i in range(net.no_bus):
        net.bus[i].Q = 0.0
        for j in range(net.no_bus):
            net.bus[i].Q += (net.bus[i].Vm * net.bus[j].Vm*(net.G[i,j] * np.sin(net.bus[i].delta - net.bus[j].delta) - net.B[i,j] * np.cos(net.bus[i].delta - net.bus[j].delta)))

# =============================================================================


# =============================================================================
def gauss_siedal(net):
    E = 10e-5
    ALPHA = 1.6
    
    itr = 1
    itrmax = 100
    
    power_injection(net)
    bus_initialization(net)
         
    for itr in range(1,itrmax):
        for i in range(net.no_bus):
            if net.bus[i].type == 2:
                s = 0.0
                for k in range(net.no_bus):
                    s = s + (net.Y[i][k])*(net.bus[k].V)
                    
                a = np.conj(net.bus[i].V) * s
                net.bus[i].Q = -np.imag(a)
                
                s = 0.0
                for k in range(net.no_bus):
                    if i != k:
                        s = s + (net.Y[i][k])*(net.bus[k].V)
                        
                b = complex(net.bus[i].P , -net.bus[i].Q)
                c = b /(np.conj(net.bus[i].V))
                    
                V_new = (1/net.Y[i][i])*(c-s)
                
                V_old = net.bus[i].V
                delta_old = net.bus[i].delta
                
                V_new = net.bus[i].V + ALPHA * (V_new - V_old)
                
                net.bus[i].delta = np.angle(V_new)
                net.bus[i].V = pol2cart(net.bus[i].Vm,net.bus[i].delta)
                
                net.bus[i].del_V = np.abs(net.bus[i].V - V_old)
                net.bus[i].del_delta = np.abs(net.bus[i].delta - delta_old)
                
                
            if net.bus[i].type == 3:
                s = 0.0
                for k in range(net.no_bus):
                    if i != k:
                        s = s + net.Y[i][k] * net.bus[k].V
                    
                b = complex(net.bus[i].P, -net.bus[i].Q)
                # print(f'b = {b}')
                c = b /(np.conj(net.bus[i].V))
                # print(f'c = {c}')
                
                V_new = (1/net.Y[i][i])*(c-s) 
                
                Vm_old = net.bus[i].Vm
                delta_old = net.bus[i].delta
                
                net.bus[i].V = net.bus[i].V + ALPHA * (V_new - net.bus[i].V)
                
                net.bus[i].Vm = np.abs(net.bus[i].V)
                net.bus[i].delta = np.angle(net.bus[i].V)
                
                net.bus[i].del_V = np.abs(net.bus[i].Vm - Vm_old)
                net.bus[i].del_delta = np.abs(net.bus[i].delta - delta_old)
                
        
        max_del_V = np.abs(net.bus[0].del_V)
        max_del_delta = np.abs(net.bus[0].del_delta)
        
        for i in range(1, net.no_bus):
            if max_del_V < np.abs(net.bus[i].del_V):
                max_del_V = np.abs(net.bus[i].del_V)
            if max_del_delta < np.abs(net.bus[i].del_delta):
                max_del_delta = np.abs(net.bus[i].del_delta)
        
        if max_del_V <= E and max_del_delta <= E:
            print(f"The iteration converges in {itr}th iteration.")
            break
    
    for i in range(net.no_bus):
        calculate_P(net)
        calculate_Q(net)
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)
    
    current_injection(net)
    
    print_loadflow_result(net)    
    
    feeder_current(net)
    
    transformer_current(net)
    
    gen_and_load_currents(net)
    
    print_feeder_powerflow(net)
    
    KCL_mismatch(net)
# =============================================================================

# =============================================================================
def FDLF(net):
    
    del_delta = np.zeros((net.no_bus - 1), dtype = float)
    del_PbyV = np.zeros((net.no_bus - 1), dtype = float)
    del_V = np.zeros((net.no_pq), dtype = float)
    del_QbyV = np.zeros((net.no_pq), dtype = float)
    
    E = 0.00001
    itr = 1
    itrmax = 100
    
    power_injection(net)
    bus_initialization(net)
    
    for itr in range(1, itrmax):
        m = 0
        n = 0
        
        calculate_P(net)
        
        for i in range(net.no_bus):
            if net.bus[i].type != 1:
                net.bus[i].delP = net.bus[i].Psp - net.bus[i].P
                del_PbyV[m] = net.bus[i].delP / net.bus[i].Vm
                m += 1
                
        max_delP = np.abs(net.bus[0].delP)
        
        for i in range(1,net.no_bus):
            if max_delP < np.abs(net.bus[i].delP):
                max_delP = np.abs(net.bus[i].delP)
            
        
        if max_delP > E:
            
            del_delta = (net.Bd_inv)@(del_PbyV)
            
            m = 0
            for i in range(net.no_bus):
                if net.bus[i].type != 1:
                    net.bus[i].delta = net.bus[i].delta + del_delta[m]
                    m += 1
        
        calculate_Q(net)
        
        for i in range(net.no_bus):
            if net.bus[i].type == 3:
                net.bus[i].delQ = net.bus[i].Qsp - net.bus[i].Q
                # print(f"Bus {i + 1} Qsp is = {net.bus[i].Qsp}")
                # print(f"Bus {i + 1} Q is = {net.bus[i].Q}")
                del_QbyV[n] = net.bus[i].delQ / net.bus[i].Vm
                n += 1
        
        max_delQ = np.abs(net.bus[0].delQ)
        
        for i in range(1,net.no_bus):
            if max_delQ < np.abs(net.bus[i].delQ):
                max_delQ = np.abs(net.bus[i].delQ)
                
        if max_delQ > E:
            
            del_V = (net.Bdd_inv)@(del_QbyV)
            
            n = 0
            for i in range(net.no_bus):
                if net.bus[i].type == 3:
                    net.bus[i].Vm = net.bus[i].Vm + del_V[n]
                    n += 1
        
        
        for i in range(net.no_bus):
            net.bus[i].V = pol2cart(net.bus[i].Vm, net.bus[i].delta)
           
        if max_delP <= E and max_delQ <= E:
            print(f"\nThe iteration is converged in {itr}th iteration")
            break
    
    current_injection(net)
    
    for i in range(net.no_bus):
        calculate_P(net)
        calculate_Q(net)
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)
            
    calculate_PgandQg(net)
     
    print_loadflow_result(net)    
    
    feeder_current(net)
    
    transformer_current(net)
    
    gen_and_load_currents(net)
    
    print_feeder_current(net)
    
    print_feeder_powerflow(net)
    
    KCL_mismatch(net)
# =============================================================================
    
# =============================================================================
def alternative_FDLF(net):
    
    del_delta = np.zeros(net.no_bus, dtype = float)
    del_PbyV = np.zeros(net.no_bus, dtype = float)
    del_V = np.zeros(net.no_bus, dtype = float)
    del_QbyV = np.zeros(net.no_bus, dtype = float)
    
    E = 0.00001
    itr = 1
    itrmax = 100

    power_injection(net)
    bus_initialization(net)
    
    for itr in range(1, itrmax):
        
        calculate_P(net)
        
        for i in range(net.no_bus):
            if net.bus[i].type != 1:
                net.bus[i].delP = net.bus[i].Psp - net.bus[i].P
                del_PbyV[i] = net.bus[i].delP / net.bus[i].Vm
        
                
        max_delP = np.abs(net.bus[0].delP)
        
        for i in range(1,net.no_bus):
            if max_delP < np.abs(net.bus[i].delP):
                max_delP = np.abs(net.bus[i].delP)
            
        if max_delP > E:
            del_delta = (net.Bd_inv)@(del_PbyV)
            
            for i in range(net.no_bus):
                net.bus[i].delta += del_delta[i]
                    
        calculate_Q(net)
        
        for i in range(net.no_bus):
            if net.bus[i].type == 3:
                net.bus[i].delQ = net.bus[i].Qsp - net.bus[i].Q
                del_QbyV[i] = net.bus[i].delQ / net.bus[i].Vm
            
        
        max_delQ = np.abs(net.bus[0].delQ)
        
        for i in range(1,net.no_bus):
            if max_delQ < np.abs(net.bus[i].delQ):
                max_delQ = np.abs(net.bus[i].delQ)
                
        if max_delQ > E:
            del_V = (net.Bdd_inv)@(del_QbyV)
            
            for i in range(net.no_bus):
                net.bus[i].Vm += del_V[i]
                      
        
        for i in range(net.no_bus):
            net.bus[i].V = pol2cart(net.bus[i].Vm, net.bus[i].delta)
        
            
        if max_delP <= E and max_delQ <= E:
            print(f"\nThe load flow is converged in {itr}th iteration.")
            break
    current_injection(net)
    
    for i in range(net.no_bus):
        calculate_P(net)
        calculate_Q(net)
        net.bus[i].S = complex(net.bus[i].P, net.bus[i].Q)
        
    print_loadflow_result(net)    
    
    feeder_current(net)
    
    transformer_current(net)
    
    gen_and_load_currents(net)
    
    print_feeder_powerflow(net)
    
    KCL_mismatch(net)
# =============================================================================


# =============================================================================
def print_loadflow_result(net):
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

    print("\n---------------Load Flow Result---------------\n")
    print(tabulate(data,  headers="firstrow",numalign="right", floatfmt=".5f", tablefmt="fancy_grid"))
# =============================================================================

# =============================================================================
def feeder_current(net):
    for k in range(net.no_fd):
        i = net.fd[k].end_bus[0]
        j = net.fd[k].end_bus[1]
        
        net.fd[k].I[0] = net.fd[k].y * (net.bus[i].V - net.bus[j].V) + (net.fd[k].ys/2) * net.bus[i].V
        net.fd[k].I[1] = net.fd[k].y * (net.bus[j].V - net.bus[i].V) + (net.fd[k].ys/2) * net.bus[j].V
        
        net.fd[k].S[0] = net.bus[i].V * np.conj(net.fd[k].I[0])
        net.fd[k].S[1] = net.bus[j].V * np.conj(net.fd[k].I[1])
        
        net.fd[k].S_loss = net.fd[k].S[0] +  net.fd[k].S[1]
# =============================================================================

# =============================================================================
def transformer_current(net):
    for k in range(net.no_tr):
        i = net.tr[k].end_bus[0]
        j = net.tr[k].end_bus[1]
        a = net.tr[k].alpha
        
        net.tr[k].I[0] = np.conj(a) * a * net.tr[k].y * net.bus[i].V - np.conj(a) * net.tr[k].y * net.bus[j].V
        net.tr[k].I[1] = - a * net.tr[k].y * net.bus[i].V + net.tr[k].y * net.bus[j].V
        
        net.fd[k].S[0] = net.bus[i].V * np.conj(net.tr[k].I[0])
        net.fd[k].S[1] = net.bus[j].V * np.conj(net.tr[k].I[1])
        
        # print(f'The {k + 1} Transformer current is {net.tr[k].I}')
# =============================================================================


# =============================================================================
def gen_and_load_currents(net):
    
    calculate_P(net)
    calculate_Q(net)
    
    for i in range(net.no_bus):
        net.bus[i].Pg = net.bus[i].P + net.bus[i].Pd
        net.bus[i].Qg = net.bus[i].Q + net.bus[i].Qd

    
        if net.bus[i].Pg != 0 or net.bus[i].Qg != 0:
            Sg = complex(net.bus[i].Pg, net.bus[i].Qg)
            net.bus[i].gen_I = np.conj(Sg/net.bus[i].V)
        else:
            net.bus[i].gen_I = 0.0
        
        if net.bus[i].Pd != 0 or net.bus[i].Qd != 0:
            Sd = complex(net.bus[i].Pd, net.bus[i].Qd)
            net.bus[i].load_I = np.conj(Sd/net.bus[i].V)
        else:
            net.bus[i].load_I = 0.0

# =============================================================================

# =============================================================================
def print_feeder_current(net):
    feeder = ["Feeder"]
    from_bus = ["FromBus"]
    to_bus = ["ToBus"]
    from_to_Bus_I_Mag = ["F-T I Mag(pu)"]
    from_to_Bus_I_Ang = ["F-T I Angle(deg)"]
    to_from_Bus_I_Mag = ["T-F I Mag(pu)"]
    to_from_Bus_I_Ang = ["T-F I Angle(deg)"]
    for i in range(net.no_fd):
        feeder.append(str(net.fd[i].code))
        from_bus.append(net.fd[i].end_bus[0] + 1)
        to_bus.append(net.fd[i].end_bus[1] + 1)
        from_to_Bus_I_Mag.append(round(np.abs(net.fd[i].I[0]), 5))
        from_to_Bus_I_Ang.append(round(np.angle(net.fd[i].I[0])* 180/np.pi, 5))
        to_from_Bus_I_Mag.append(round(np.abs(net.fd[i].I[1]), 5))
        to_from_Bus_I_Ang.append(round(np.angle(net.fd[i].I[1])* 180/np.pi, 5))

    data = np.transpose([feeder, from_bus, to_bus, from_to_Bus_I_Mag, from_to_Bus_I_Ang, to_from_Bus_I_Mag, to_from_Bus_I_Ang])
    
    print("\nFeeder Currents\n")
    print(tabulate(data,  headers="firstrow",numalign="right", floatfmt=".5f", tablefmt="fancy_grid"))
# =============================================================================

# =============================================================================
def print_feeder_powerflow(net):
    feeder = ["Feeder"]
    from_bus = ["FromBus"]
    to_bus = ["ToBus"]
    from_to_Bus_P = ["F-T P(pu)"]
    from_to_Bus_Q = ["F-T Q(pu)"]
    to_from_Bus_P = ["T-F P(pu)"]
    to_from_Bus_Q = ["T-F Q(pu)"]
    for i in range(net.no_fd):
        feeder.append(str(net.fd[i].code))
        from_bus.append(net.fd[i].end_bus[0] + 1)
        to_bus.append(net.fd[i].end_bus[1] + 1)
        from_to_Bus_P.append(round(np.real(net.fd[i].S[0]), 5))
        from_to_Bus_Q.append(round(np.imag(net.fd[i].S[0]), 5))
        to_from_Bus_P.append(round(np.real(net.fd[i].S[1]), 5))
        to_from_Bus_Q.append(round(np.imag(net.fd[i].S[1]), 5))

    data = np.transpose([feeder, from_bus, to_bus, from_to_Bus_P, from_to_Bus_Q, to_from_Bus_P, to_from_Bus_Q])
    
    print("\nPower flowing through Feeders\n")
    print(tabulate(data,  headers="firstrow",numalign="right", floatfmt=".5f", tablefmt="fancy_grid"))
# =============================================================================

# =============================================================================
def calculate_PgandQg(net):
    net.Ng = []
    for i in range(net.no_gen):
        m = net.gen[i].bus
        net.Ng.append(m)
    
    for i in range(len(net.Ng)):
        j = net.Ng[i]
        net.bus[j].Pg = net.bus[j].P + net.bus[j].Pd
        net.bus[j].Qg = net.bus[j].Q + net.bus[j].Qd
        
# =============================================================================

# # =============================================================================
# def print_feeder_result(net):
#     feeder = ["Feeder"]
#     from_bus = ["FromBus"]
#     to_bus = ["ToBus"]
#     from_to_Bus_I = ["F-T Current(pu)"]
#     to_from_Bus_I = ["T-F Current(pu)"]
#     from_to_Bus_S = ["F-T Power Flow(pu)"]
#     to_from_Bus_S = ["T-F Power Flow(pu)"]
#     for i in range(net.no_fd):
#         feeder.append(str(net.fd[i].code))
#         from_bus.append(net.fd[i].end_bus[0] + 1)
#         to_bus.append(net.fd[i].end_bus[1] + 1)
#         from_to_Bus_I.append(round(net.fd[i].I[0], 5))
#         to_from_Bus_I.append(round(net.fd[i].I[1], 5))
#         from_to_Bus_S.append(round(net.fd[i].S[0], 5))
#         to_from_Bus_S.append(round(net.fd[i].S[1], 5))

#     data = np.transpose([feeder, from_bus, to_bus, from_to_Bus_I, to_from_Bus_I, from_to_Bus_S, to_from_Bus_S])
    
#     print("\nFeeder Current and Power Flow\n")
#     print(tabulate(data,  headers="firstrow",numalign="right", floatfmt=".5f", tablefmt="rst"))
# # =============================================================================


# # =============================================================================
# def KCL_checking(net):
#     sum_I = np.zeros((net.no_bus),dtype = complex)
    
#     for i in range(net.no_bus):
#         sum_I[i] = net.bus[i].load_I
        
#         for j in range (net.no_fd):
#             if net.bus[i].code == net.fd[j].end_bus[0]:
#                 sum_I[i] = sum_I[i] + net.fd[j].I[0]        #outgoing current is taken as +ve.
#             elif net.bus[i].code == net.fd[j].end_bus[1]:
#                 sum_I[i] = sum_I[i] + net.fd[j].I[1]
                
#         for k in range(net.no_tr):
#             if net.bus[i].code == net.tr[k].end_bus[0]:
#                 sum_I[i] = sum_I[i] + net.tr[k].I[0]
#             elif net.bus[i].code == net.tr[k].end_bus[1]:
#                 sum_I[i] = sum_I[i] + net.tr[k].I[1]
            
#         for k in range(net.no_gen):
#             if net.bus[i].code == net.gen[k].bus:
#                 sum_I[i] = sum_I[i] - net.bus[i].gen_I      #incoming current is taken as -ve.
                
#     temp = np.max(np.abs(sum_I))       
    
#     # for i in range(net.no_bus):
#     #     print(f"Bus {i + 1} of type {net.bus[i].type} KCL error is = {np.abs(sum_I[i])}")
        
#     print(f"The maximum error in KCL checking is {temp}")
# # =============================================================================

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
        
    print(f"\nThe maximum mismatch in KCL checking is {temp}\n") 
     
# =============================================================================
    
