""" Created on Wed May 10 14:45:10 2023 @author: chatt """
import pandas as pd
import numpy as np
import copy as copy
import math

from NetworkData import *
from load_flow_combined import *
from short_circuit_analysis import *
from step_distance_protection import *

#==============================================================================
def line_outage(net, outage_fd):
    x = outage_fd
    i= net.fd[x].end_bus[0]
    j= net.fd[x].end_bus[1]
    
    net.Y[i,i] = net.Y[i,i] - net.fd[x].y - net.fd[x].ys/2
    net.Y[i,j] = net.Y[i,j] + net.fd[x].y 
    net.Y[j,i] = net.Y[j,i] + net.fd[x].y
    net.Y[j,j] = net.Y[j,j] - net.fd[x].y - net.fd[x].ys/2
    
#------------------------------------------------------------------------------        

#==============================================================================
def load_change(net):
    for i in range(net.no_bus):
        if net.bus[i].type == 3:
            net.bus[i].Pd = 3.25 * net.bus[i].Pd
            net.bus[i].Qd = 2 * net.bus[i].Qd
    
    for i in range(net.no_bus):
        if net.bus[i].type == 2:
            net.bus[i].Pg = 2 * net.bus[i].Pg
            
#------------------------------------------------------------------------------


#==============================================================================
def participation_factor(net):
    sum_Pg = 0
    for i in range(net.no_gen):
        j = net.gen[i].bus
        sum_Pg = sum_Pg + net.bus[i].Pg
    
    for i in range(net.no_gen):
        net.gen[i].participation_factor = net.bus[i].Pg / sum_Pg
        
#------------------------------------------------------------------------------


#==============================================================================
# def load_distribution(net):
#     delta_Pd = 
          