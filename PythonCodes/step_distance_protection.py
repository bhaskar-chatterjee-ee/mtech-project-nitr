import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
from tabulate import tabulate

from NetworkData import *
from load_flow_combined import *

# =============================================================================
def MhoRelayCharacteristics(net):
    for i in range(net.no_mhorelay):
        if net.mhorelay[i].status == 1:
            plt.figure()
            x = 0; y = 0
            zeta = np.linspace(0, 2 * np.pi, 100 )
            x = (net.mhorelay[i].Zset1/2) * np.cos(zeta) + (net.mhorelay[i].Zset1/2)*np.cos(net.mhorelay[i].MTA)
            y = (net.mhorelay[i].Zset1/2) * np.sin(zeta) + (net.mhorelay[i].Zset1/2)*np.sin(net.mhorelay[i].MTA)
            plt.plot(x, y, color = 'blue')
            x = 0; y = 0
            zeta = np.linspace(0, 2 * np.pi, 100 )
            x = (net.mhorelay[i].Zset2/2) * np.cos(zeta) + (net.mhorelay[i].Zset2/2)*np.cos(net.mhorelay[i].MTA)
            y = (net.mhorelay[i].Zset2/2) * np.sin(zeta) + (net.mhorelay[i].Zset2/2)*np.sin(net.mhorelay[i].MTA)
            plt.plot(x, y, color = 'black')
            x = 0; y = 0
            zeta = np.linspace(0, 2 * np.pi, 100 )
            x = (net.mhorelay[i].Zset3/2) * np.cos(zeta) + (net.mhorelay[i].Zset3/2)*np.cos(net.mhorelay[i].MTA)
            y = (net.mhorelay[i].Zset3/2) * np.sin(zeta) + (net.mhorelay[i].Zset3/2)*np.sin(net.mhorelay[i].MTA)
            plt.plot(x, y, color = 'red')
            plt.gca().set_aspect('equal')
            plt.grid()
            plt.title(f'Characteristics of Mho Relay {net.mhorelay[i].relaycode} near Bus {net.mhorelay[i].neartobus + 1} and on Feeder {net.mhorelay[i].onfeeder + 1}')
            x = net.mhorelay[i].Zseen.real
            y = net.mhorelay[i].Zseen.imag
            plt.scatter(x,y)
            plt.savefig(f"Relay{i+1}.png")
            # plt.clf()
# =============================================================================        

# =============================================================================        
def RelaySensing(net):
    for i in range(net.no_mhorelay):
        if net.mhorelay[i].status == 1:
            j = net.mhorelay[i].neartobus
            k = net.mhorelay[i].onfeeder
            # print(i,j,k)
            net.mhorelay[i].Zseen = net.bus[j].V/net.fd[k].I[0]
            # print(net.mhorelay[i].Zseen)
        
# =============================================================================


    
        