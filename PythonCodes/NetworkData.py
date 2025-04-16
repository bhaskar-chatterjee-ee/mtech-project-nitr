import pandas as pd
import copy
import numpy as np


# ================================== Bus Data =================================
class Bus:
    def __init__(self):
        self.code = 0    #busNo
        self.zone = 0
        self.Pg = 0.0
        self.Pd = 0.0
        self.Qg = 0.0
        self.Qd = 0.0
        self.type = 0       #slack or PV or PQ
        self.Vm = 0.0    #voltage mag
        self.delta = 0.0
        self.V = 0.0 + 1j * 0.0      #voltage
        self.del_V = 0.0 + 1j * 0.0
        self.del_delta = 0.0
        self.actual_code = 0.0
        self.P = 0.0               #Pcal
        self.Q = 0.0               #Qcal
        self.delP = 0.0
        self.delQ = 0.0
        self.Psp = 0.0
        self.Qsp = 0.0
        self.gen = -1
        self.gen_I = 0.0 + 1j * 0.0
        self.load_I = 0.0 + 1j * 0.0
        self.yload =  0.0 + 1j * 0.0
        self.load_status = 0
        self.S = 0.0 + 1j * 0.0
# =============================================================================


# ================================ Line Data ==================================
class Feeder:
    def __init__(self):
        self.code = 0
        self.end_bus = [0, 0]
        self.r = 0.0
        self.x = 0.0
        self.ys = 0.0
        self.reac_frombus = 0.0
        self.reac_tobus = 0.0
        self.status = 0
        self.y = 0.0
        self.z = 0.0
        self.I = [0.0 + 1j * 0.0, 0.0 + 1j * 0.0]
        self.S = [0.0 + 1j * 0.0, 0.0 + 1j * 0.0]
        self.Plimit = 0.0
        # self.Pijmax = 0.0
        # self.Pijmin = 0.0
# =============================================================================        


# ================================ Transformer Data ===========================
class Transformer:
    def __init__(self):
        self.code = 0
        self.end_bus = [0, 0]
        self.r = 0.0
        self.x = 0.0
        self.alpha = 0.0
        self.status = 0
        self.y = 0.0
        self.z = 0.0
        self.I = [0.0 + 1j*0.0, 0.0 +1j*0.0]
        self.S = [0.0, 0.0]
# =============================================================================       

# ===============================  PV Bus Data ================================
class PVBUs:
    def __init__(self):
        self.code = 0
        self.pmin = 0.0
        self.pmax = 0.0
        self.qmin = 0.0
        self.qmax = 0.0
        self.Vsp = 0.0
# =============================================================================

# ==============================  Slack Bus Data ==============================
class SlackBus:
    def __init__(self):
        self.code = 0
        self.pmin = 0.0
        self.pmax = 0.0
        self.qmin = 0.0
        self.qmax = 0.0
        self.Vsp = 0.0
# =============================================================================

# =============================  Shunt load Data ==============================
class ShuntLoad:
    def __init__(self):
        self.code = 0
        self.g = 0.0
        self.b = 0.0
        self.status = 0
        self.y = 0.0 + 1j* 0.0
        self.I = 0.0
        self.S = 0.0
# =============================================================================

# ==============================  GeneratorData ===============================
class Generator:
    def __init__(self):
        self.bus = -1
        self.rated_mva = 0.0
        self.H = 0.0
        self.D = 0.0
        self.Ra = 0.0
        self.Xd = 0.0
        self.Xq = 0.0
        self.Xd_d = 0.0
        self.Xq_d = 0.0
        self.I = 0.0 + 1j*0.0
        self.status = 1
        self.sync_condenser = 0
        self.Pe = 0.0
        self.E = 0.0 + 1j*0.0
        self.Em = 0.0
        self.delta = 0.0
        self.w = 0.0
        self.f = 0.0
        self.y = 0.0 + 1j*0.0
        self.Pm = 0.0
        self.delta_dot_1 = 0.0
        self.delta_dot_2 = 0.0  # Predicted slope
        self.w_dot_1 = 0.0
        self.w_dot_2 = 0.0      #Predicted slope
        self.delta_0 = 0.0      # Predicted delta
        self.w_0 = 0.0          # Predicted w
        self.R = 0              # Regulation
        self.Pmax = 0.0
        self.Pmin = 0.0
        self.Pset = 0.0
        self.delta_Pset = 0.0
        self.participation_factor = 0.0
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.Qmin = 0.0
        self.Qmax = 0.0
        self.Pset_max = 0.0
        self.Pset_min = 0.0
        self.IG = 0.0
        self.gamma = 0.0
        self.Id = 0.0
        self.Iq = 0.0
        self.Vd = 0.0
        self.Vq = 0.0
        self.Ed_d = 0.0
        self.Eq_d = 0.0
        self.Efd = 0.0   
        self.VR = 0.0   
        self.Rf = 0.0   
        self.Vref = 0.0   
        self.TM = 0.0   
        self.M = 0.0   
        self.Tr = 0.0   
        self.Ka = 0.0   
        self.Ta = 0.0   
        self.VRmax = 0.0   
        self.VRmin = 0.0   
        self.Ke = 0.0   
        self.Te = 0.0   
        self.Kf = 0.0   
        self.Tf = 0.0   
        self.E1 = 0.0   
        self.SE1 = 0.0   
        self.E2 = 0.0   
        self.SE2 = 0.0  
        self.Tch = 0.0
        self.Tsv = 0.0
        self.Rd = 0.0
        self.Psv = 0.0
        self.Pc = 0.0
        self.Qg = 0.0
        self.Td0_d = 0.0
        self.Tq0_d = 0.0
        self.Id_0 = 0.0
        self.Iq_0 = 0.0
# =============================================================================

# ================================  Disturbance ===============================   
class Disturbance:
    def __init__(self):
        self. initiation_time = 0.0
        self.clearance_time = 0.0
        self.type = 0
        self.bus = 0
        self.feeder = 0
# =============================================================================

# =============================================================================
class MhoRelay:
    def __init__(self):
        self.code = 0
        self.neartobus = 0
        self.onfeeder = 0
        self.Zset1 = 0.0
        self.Zset2 = 0.0
        self.Zset3 = 0.0
        self.MTA = 0.0
        self.status = 0
        self.Zseen = 0.0
# =============================================================================

# =============================================================================
class OCRelay:
    def __init__(self):
        self.code = 0
        self.generator = 0
        self.Ipk = 0
        self.status = 0
# =============================================================================

# =============================================================================
class Measurement:
    def __init__(self):
        self.type = 0
        self.onfeeder = 0
        self.onbus = 0
        self.ontr = 0
        self.reading = 0

# =============================================================================


# =============================================================================
class Network:
    def __init__(self):
        self.no_bus = 0
        self.no_fd = 0
        self.no_tr = 0
        self.no_pv = 0
        self.no_pq = 0
        self.no_shunt = 0
        self.no_gen = 0
        self.no_exciter = 0
        self.base_MVA = 0.0
        self.bus = []
        self.fd = []
        self.tr = []
        self.pv_bus = []
        self.slack_bus = []
        self.shunt = []
        self.gen = []
        self.disturbance = Disturbance()
        self.frequency = 50.0
        self.ws = 2 * np.pi * 50
        self.Y = []
        #Y-Bus Formation
        self.alter_Y = []
        self.Np =[]
        self.Npq = []
        self.Bd = []
        self.Bdd = []
        self.Bd_inv = []
        self.Bdd_inv = []
        self.G = []
        self.B = []
        self.Ng = []     #Buses at which generators are present
        #Fault on Bus
        self.fault_bus = None       #Bus on which fault has occured
        self.mod_Y =[]
        self.Z = []
        self.Z_f = 0
        #Fault on Feeder
        self.fault_line = None      #Feeder on which fault has occured
        self.fault_dist = 0         #Distance of the fault from end_bus[0] of the feeder
        self.new_Y = []
        self.new_G = []
        self.new_B = []
        #Relay Characteristics
        self.no_mhorelay = 0
        self.mhorelay = []
        self.no_ocrelay = 0
        self.ocrelay = []
        #State Estimation
        self.no_msrmnt = 0
        self.msrmnt = []
        self.no_V_msrmnt = 0
        self.V_msrmnt = []
        self.no_I_msrmnt = 0
        self.I_msrmnt = []
        self.outage_line = None
        #Phasor Measurement Units
        self.pmu_location = 0
        #State Estimation
        self.state_Vector = []
        

    def add_bus(self, b):
        self.bus.append(b)

    def add_feeder(self, f):
        self.fd.append(f)

    def add_transformer(self, tr):
        self.tr.append(tr)

    def add_PV_bus(self, pv):
        self.pv_bus.append(pv)

    def add_slack_bus(self, s):
        self.slack_bus.append(s)

    def add_shuntload(self, s):
        self.shunt.append(s)

    def add_generator(self, g):
        self.gen.append(g)
        
    def add_mhorelay(self, m):
        self.mhorelay.append(m)
    
    def add_ocrelay(self, o):
        self.ocrelay.append(o)
        
    def add_msrmnt(self, m):
        self.msrmnt.append(m)
    
# =============================================================================     
    

# =============================================================================
def read_bus_data(net, data_file):
    # =========== Read Bus Data ========================
    df = pd.read_excel(data_file, sheet_name='Bus Data')

    code = df['Code']
    Zone = df['Zone']
    Pg = df['Pg']
    Qg = df['Qg']
    Pd = df['Pd']
    Qd = df['Qd']
    bustype = df['Type']

    for i in range(net.no_bus):
        b = Bus()
        b.code = code[i]
        b.zone = Zone[i]
        b.Pg = Pg[i]
        b.Pd = Pd[i]
        b.Qg = Qg[i]
        b.Qd = Qd[i]
        b.type = bustype[i]
        net.add_bus(b)
        del b

    del Pg
    del Pd
    del Qg
    del Qd
    del bustype
    del code
    
# =============================================================================

# =============================================================================
def read_feeder_data(net, data_file):
    # =========== Read Feeder Data ========================
    if net.no_fd > 0:
        df = pd.read_excel(data_file, sheet_name='Feeder Data')

        code = df['Code']
        frombus = df['From Bus']
        tobus = df['To Bus']
        r = df['Resistance']
        x = df['Reactance']
        ys = df['Charging Admittance']
        reac_frombus = df['Reactance From Bus']
        reac_tobus = df['Reactance To Bus']
        status = df['Status']
        

        for i in range(net.no_fd):
            fd = Feeder()
            fd.code = code[i]
            fd.end_bus[0] = frombus[i]
            fd.end_bus[1] = tobus[i]
            fd.frombus = frombus[i]
            fd.tobus = tobus[i]
            fd.r = r[i]
            fd.x = x[i]
            fd.ys = ys[i]
            fd.reac_frombus = reac_frombus[i]
            fd.reac_tobus = reac_tobus[i]
            fd.status = status[i]
            net.add_feeder(fd)
            del fd

        del frombus
        del tobus
        del r
        del x
        del ys
        del code
        del reac_frombus
        del reac_tobus
        del status
# =============================================================================

# =============================================================================
def read_transformer_data(net, data_file):
    # =========== Read Transformer Data ========================
    if net.no_tr > 0:
        df = pd.read_excel(data_file, sheet_name='Transformer Data')

        code = df['Code']
        frombus = df['From Bus']
        tobus = df['To Bus']
        r = df['Resistance']
        x = df['Reactance']
        alpha = df['Off-nominal Tap Ratio']
        status = df['Status']

        for i in range(net.no_tr):
            tr = Transformer()
            tr.code = code[i]
            tr.end_bus[0] = frombus[i]
            tr.end_bus[1] = tobus[i]
            tr.frombus = frombus[i]
            tr.tobus = tobus[i]
            tr.r = r[i]
            tr.x = x[i]
            #tr.alpha = 1/alpha[i]   ============================ change 1
            # tr.alpha = 1.0  ==================================== change 2
            tr.alpha = alpha[i]
            tr.status = status[i]
            net.add_transformer(tr)
            del tr

        del frombus
        del tobus
        del r
        del x
        del code
        del alpha
        del status
        
# =============================================================================

# =============================================================================
def read_PVbus_data(net, data_file):
    # =========== Read PV Bus Data ========================
    if net.no_pv > 0:
        df = pd.read_excel(data_file, sheet_name='PV Bus Data')

        no = df['Number']
        code = df['Bus Code']
        pmin = df['Pmin']
        pmax = df['Pmax']
        qmin = df['Qmin']
        qmax = df['Qmax']
        Vsp = df['Specified Voltage']

        for i in range(net.no_pv):
            pv = PVBUs()
            pv.code = code[i]
            pv.pmin = pmin[i]
            pv.pmax = pmax[i]
            pv.qmin = qmin[i]
            pv.qmax = qmax[i]
            pv.Vsp = Vsp[i]
            net.add_PV_bus(pv)
            del pv

        del code
        del pmin
        del pmax
        del qmin
        del qmax
        del Vsp

# =============================================================================


# =============================================================================
def read_slackbus_data(net, data_file):
    # =========== Read Slack Bus Data ========================
    df = pd.read_excel(data_file, sheet_name='Slack Bus Data')

    code = df['Bus Code']
    pmin = df['Pmin']
    pmax = df['Pmax']
    qmin = df['Qmin']
    qmax = df['Qmax']
    Vsp = df['Specified Voltage']

    s = SlackBus()
    i = 0
    s.code = code[i]
    s.pmin = pmin[i]
    s.pmax = pmax[i]
    s.qmin = qmin[i]
    s.qmax = qmax[i]
    s.Vsp = Vsp[i]

    net.slack_bus = copy.deepcopy(s)
    # net.add_slack_bus(s)
    del s

    del code
    del pmin
    del pmax
    del qmin
    del qmax
    del Vsp
    
# =============================================================================


# =============================================================================
def read_shuntload_data(net, data_file):
    # =========== Read Shunt Load Data ========================
    if net.no_shunt > 0:
        df = pd.read_excel(data_file, sheet_name='Shunt Load Data')

        no = df['Number']
        code = df['Bus Code']
        G = df['G']
        B = df['B']
        status = df['Status']

        for i in range(net.no_shunt):
            sh = ShuntLoad()
            sh.code = code[i]
            sh.g = G[i]
            sh.b = B[i]
            sh.status = status[i]
            net.add_shuntload(sh)
            del sh

        del code
        del G
        del B
        del status
        
# =============================================================================


# =============================================================================
def read_generator_data(net, data_file):
    # =========== Read Generator Data ========================
    if net.no_gen > 0:
        df = pd.read_excel(data_file, sheet_name='Generator Data')


        bus = df['Bus Code']
        MVA = df['Rated MVA']
        H = df['H']
        D = df['D']
        Ra = df['Ra']
        Xd = df['Xd']
        Xq = df['Xq']
        Xd_d = df['Xd_d']
        Xq_d = df['Xq_d']
        status = df['Status']
        Td0_d = df['Td0_d']
        Tq0_d = df['Tq0_d']

        for i in range(net.no_gen):
            g = Generator()
            g.bus = bus[i]
            g.rated_mva = MVA[i]
            g.H = H[i]
            g.D = D[i]
            g.Ra = Ra[i]
            g.Xd = Xd[i]
            g.Xq = Xq[i]
            g.Xd_d = Xd_d[i]
            g.Xq_d = Xq_d[i]
            g.status = status[i]
            g.Td0_d = Td0_d[i]
            g.Tq0_d = Tq0_d[i]  

            net.add_generator(g)
            del g

        del bus
        del H
        del D
        del Ra
        del Xd
        del Xq
        del Xd_d
        del Xq_d
        
# =============================================================================

    # =========== Read exciter Data ========================                   
    if net.no_exciter > 0:
        df = pd.read_excel(data_file, sheet_name='Exciter Data')


        Tr = df['Tr']
        Ka = df['Ka']
        Ta = df['Ta']
        VRmax = df['VRmax']
        VRmin = df['VRmin']
        Ke = df['Ke']
        Te = df['Te']
        Kf = df['Kf']
        Tf = df['Tf']
        E1 = df['E1']
        SE1 = df['SE1']
        E2 = df['E2']
        SE2 = df['SE2']
        

        for i in range(net.no_gen):
            net.gen[i].Tr = Tr[i]
            net.gen[i].Ka = Ka[i]
            net.gen[i].Ta = Ta[i]
            net.gen[i].VRmax = VRmax[i]
            net.gen[i].VRmin = VRmin[i]
            net.gen[i].Ke = Ke[i]
            net.gen[i].Te = Te[i]
            net.gen[i].Kf = Kf[i]
            net.gen[i].Tf = Tf[i]
            net.gen[i].E1 = E1[i]
            net.gen[i].SE1 = SE1[i]
            net.gen[i].E2 = E2[i]
            net.gen[i].SE2 = SE2[i]
            

        del Tr
        del Ka
        del Ta
        del VRmax
        del VRmin
        del Ke
        del Te
        del Kf
        del Tf
        del E1
        del E2
        del SE1
        del SE2
    # ======================================================
    
# =============================================================================                    
 
       
# =============================================================================
def read_mhorelay_data(net, data_file):
    df = pd.read_excel(data_file, sheet_name = 'Mho Relay Settings')
    
    relaycode = df['Code']
    neartobus = df['Bus']
    onfeeder = df['Feeder']
    Zset1 = df['Zset1']
    Zset2 = df['Zset2']
    Zset3 = df['Zset3']
    MTA = df['MTA']
    status = df['Status']
    
    for i in range(net.no_mhorelay):
        m = MhoRelay()
        m.relaycode = relaycode[i]
        m.neartobus = neartobus[i]
        m.onfeeder = onfeeder[i]
        m.Zset1 = Zset1[i]
        m.Zset2 = Zset2[i]
        m.Zset3 = Zset3[i]
        m.MTA = MTA[i]
        m.status = status[i]
        
        net.add_mhorelay(m)
        del m
    
    del relaycode
    del neartobus
    del onfeeder
    del Zset1
    del Zset2
    del Zset3
    del MTA
    del status
    
# =============================================================================    
  
  
# =============================================================================
def read_overcurrentrelay_data(net, data_file):
    df = pd.read_excel(data_file, sheet_name = 'OC Relay Settings')
    
    relaycode = df['Code']
    gen = df['Generator']
    Ipk = df['Pickup Current']
    status = df['Status']
    
    for i in range(net.no_ocrelay):
        o = OCRelay()
        o.relaycode = relaycode[i]
        o.gen = gen[i]
        o.Ipk = Ipk[i]
        o.status = status[i]
        net.add_ocrelay(o)
        del o
        
    del relaycode
    del gen
    del Ipk
    del status
    
# =============================================================================    
    
# =============================================================================
def read_general_data(net, data_file):
    # =========== Read general Information ========================
    df = pd.read_excel(data_file, sheet_name='General Info')

    net.no_bus = df['no of bus'][0]
    net.no_fd = df['no of feeder'][0]
    net.no_tr = df['no of transformer'][0]
    net.no_pv = df['no of pv bus'][0]
    net.no_shunt = df['no of shunt load'][0]
    net.base_MVA = df['base mva'][0]
    net.no_gen = df['no of generator'][0]
    net.no_exciter = df['no of exciter'][0]
# =============================================================================


# =============================================================================
def read_protection_info(net, data_file):
    df = pd.read_excel(data_file, sheet_name = 'Protection General Info')
    
    net.no_mhorelay = df['no of mho relay'][0]
    net.no_ocrelay = df['no of oc relay'][0]
# =============================================================================

# =============================================================================
def read_network_data(net, data_file):
    read_general_data(net, data_file)
    read_bus_data(net, data_file)
    read_feeder_data(net, data_file)
    read_transformer_data(net, data_file)
    read_PVbus_data(net, data_file)
    read_slackbus_data(net, data_file)
    read_shuntload_data(net, data_file)
    read_generator_data(net, data_file)
    read_protection_info(net, data_file)
    read_mhorelay_data(net, data_file)
    read_overcurrentrelay_data(net, data_file)

    # ================= Renumber ==============================
    renumber_bus_code(net)

    for i in range(net.no_fd):
        net.fd[i].end_bus[0] = code_to_number(net, net.fd[i].end_bus[0])
        net.fd[i].end_bus[1] = code_to_number(net, net.fd[i].end_bus[1])

    for i in range(net.no_tr):
        net.tr[i].end_bus[0] = code_to_number(net, net.tr[i].end_bus[0])
        net.tr[i].end_bus[1] = code_to_number(net, net.tr[i].end_bus[1])

    for i in range(net.no_pv):
        net.pv_bus[i].code = code_to_number(net, net.pv_bus[i].code)

    for i in range(net.no_shunt):
        net.shunt[i].code = code_to_number(net, net.shunt[i].code)

    net.slack_bus.code = code_to_number(net, net.slack_bus.code)

    for i in range(net.no_gen):
        net.gen[i].bus = code_to_number(net, net.gen[i].bus)
        
    for i in range(net.no_mhorelay):
        net.mhorelay[i].code = code_to_number(net, net.mhorelay[i].code)
        net.mhorelay[i].neartobus = code_to_number(net, net.mhorelay[i].neartobus)
        net.mhorelay[i].onfeeder = code_to_number(net, net.mhorelay[i].onfeeder)

    process_network_data(net)
    
# =============================================================================


# =============================================================================
def process_network_data(net):
    # =========== Process Network Data ========================

    for i in range(net.no_bus):
        net.bus[i].Pg /= net.base_MVA
        net.bus[i].Pd /= net.base_MVA
        net.bus[i].Qg /= net.base_MVA
        net.bus[i].Qd /= net.base_MVA

    for i in range(net.no_pv):
        net.pv_bus[i].pmin /= net.base_MVA
        net.pv_bus[i].pmax /= net.base_MVA
        net.pv_bus[i].qmin /= net.base_MVA
        net.pv_bus[i].qmax /= net.base_MVA

    net.slack_bus.pmin /= net.base_MVA
    net.slack_bus.pmax /= net.base_MVA
    net.slack_bus.qmin /= net.base_MVA
    net.slack_bus.qmax /= net.base_MVA
    
    for i in range(net.no_gen):
        net.gen[i].a *= net.base_MVA ** 2
        net.gen[i].b *= net.base_MVA
        net.gen[i].Pmin /= net.base_MVA
        net.gen[i].Pmax /= net.base_MVA
        
        
    for i in range(net.no_fd):
        net.fd[i].z = net.fd[i].r + 1j * net.fd[i].x
        net.fd[i].y = 1 / net.fd[i].z
        net.fd[i].ys = 1j * net.fd[i].ys
        net.fd[i].tempz = 1j * net.fd[i].x
        net.fd[i].tempy = 1 / net.fd[i].tempz
        # net.fd[i].Pijmax /= net.base_MVA

    for i in range(net.no_tr):
        net.tr[i].z = net.tr[i].r + 1j * net.tr[i].x
        net.tr[i].y = 1 / net.tr[i].z
        # net.tr[i].aplha = 1/net.tr[i].alpha

    for i in range(net.no_shunt):
        net.shunt[i].y = (net.shunt[i].g + 1j * net.shunt[i].b)

    for i in range(net.no_bus):
        if abs(net.bus[i].Pd) != 0 or abs(net.bus[i].Qd) != 0:
            net.bus[i].load_status = 1

    for i in range(net.no_gen):
        net.gen[i].H *= net.gen[i].rated_mva / net.base_MVA
        net.gen[i].D *= net.gen[i].rated_mva / net.base_MVA
        net.gen[i].Ra *= net.base_MVA / net.gen[i].rated_mva
        
        net.gen[i].Xd_d *= net.base_MVA / net.gen[i].rated_mva
        net.gen[i].Xd *= net.base_MVA / net.gen[i].rated_mva
        net.gen[i].Xq_d *= net.base_MVA / net.gen[i].rated_mva
        net.gen[i].Xq *= net.base_MVA / net.gen[i].rated_mva
        
        b = net.gen[i].bus
        net.bus[b].gen = i
        z = net.gen[i].Ra + 1j * net.gen[i].Xd_d
        net.gen[i].y = 1 / z
        net.gen[i].D *= 5
        # mark synchronous condenser
        if net.bus[b].Pg == 0:
            net.gen[i].sync_condenser = 1
        # change D
        # print(net.gen[i].D)
            
        if net.gen[i].D == 0:
            if net.no_bus != 9:
                net.gen[i].D = 2 # 2
        if net.no_bus == 9:
            net.gen[i].D *= 5
            
        # net.gen[i].D *= 5


    net.frequency = 50
    net.ws = 2 * np.pi * net.frequency
    
# =============================================================================

# =============================================================================
def code_to_number(net, i):
    n = 0
    for k in range(net.no_bus):
        if net.bus[k].actual_code == i:
            n = k
            break
    return n
# =============================================================================

# =============================================================================
def renumber_bus_code(net):
    for i in range(net.no_bus):
        net.bus[i].actual_code = net.bus[i].code
        net.bus[i].code = i
# =============================================================================


