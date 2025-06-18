import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as ip
from datetime import datetime 
from scipy.special import gammaincc


 #parameters & constants

N_a = 6.02214e23 # Avogadro const.
m_e = 511e3# in eV
gA = 1.27 # Axial vector weak coupling constant

### Xenon 136 parameters & constants
ia = 0.089 # Natural abundance of 136Xe -- this is the abundance per amount (mole), not per mass!
M_Xe136 = 135.907e-6 # t per mol 
M_Xe = 131.293e-6 # t per mol 

E_Qbb = 2457.8 
G_0nu_Xe136 = 14.58e-15 #1/yr ([https://arxiv.org/pdf/1209.5722.pdf])
NME_Xe136 = [1.11, 4.77]

eff = 0.683 # = signal coverage of FWHM-ROI 

### Background Models ###
### Intrinsic Model 
# assuming flat spectra in the ROI, a priori no MS-rejection

def Intr_Bkg(Xe137_Key, BiPo_eff, enrichment, flux):
    
    ### NONE OF THESE NUMBERS INCLUDE SS/MS DISCRIMINATION, those factors are applied later
    b8 =     1.77e-4  ### LZ number with updated Flux = 5.25e-6 cm-2 s-1, and complete P_ee 
                      ### corresponding to 1.36e-4 (=3.72e-10 dru) after SS acceptance with 3mm
    #b8 =     1.84e-4  ### NEW DARWIN number with updated Flux = 5.25e-6 cm-2 s-1, and complete P_ee 
    #xe136 =  7.05e-6  ### DARWIN number corresponding to 0.8% E_res
    rn222 =  1.90e-4 * (100-BiPo_eff)/(100-99.9)  ### assuming 99.9% BiPo, 0.1mBq/t ==> Rescaled

    xe136 = 5.04e-6   ### Using smeared Kotila&Iachello spectrum (10 keV bins), 0.65% E_res, 2.165e21 yr HL
    Ffactor = 2.05e-4 # Fabian's sim results corrected for the new Xe-137 production: 0.822/0.79
    if enrichment == 0: # scaling with the new production rates from Jose for different enrichment scenarios
        Ffactor = Ffactor * 0/0.822
        xe136 = 0
    elif enrichment == 5:
        Ffactor = Ffactor * 0.551/0.822
        xe136 = 2.84e-06
    elif enrichment == 10:
        Ffactor = Ffactor * 0.872/0.822
        xe136 = 5.68e-06
    elif enrichment == 15:
        Ffactor = Ffactor * 1.29/0.822
        xe136 = 8.52e-06
    elif enrichment == 20:
        Ffactor = Ffactor * 1.38/0.822
        xe136 = 1.14e-05
    elif enrichment == 25:
        Ffactor = Ffactor * 1.84/0.822
        xe136 = 1.42e-05
    
    if Xe137_Key == "LNGS": xe137 =  Ffactor
    if Xe137_Key == "SURF": xe137 =  Ffactor * 0.142/0.822      # scaled using DARWIN Xe-137 projections for SURF
    if Xe137_Key == "SNOLab": xe137 =  Ffactor * 0.00675/0.822  # scaled using DARWIN Xe-137 projections for SNOLab
    if Xe137_Key == "Boulby-1300": xe137 =  Ffactor * 14.6/29.7 # this and the others are scaled using the muon flux
    if Xe137_Key == "Boulby": xe137 =  Ffactor * 32.3/29.7 
    if Xe137_Key == "Kamioka": xe137 =  Ffactor * 128/29.7
    if Xe137_Key == "generic": 
        xe137 = Ffactor * flux/29.7
        #print('Not one of the labs, using a flux of ',flux)
    
    return np.array([(xe137+xe136+rn222+b8), b8, rn222, xe136, xe137])
    
    
## External background model (DARWIN, 0.8%)
# assuming 0.8% energy res. and 15% gamma acceptance (= 15mm in xyz)
def Bkg_fit_func(m, exponent, power1, power2,  scale):
    return np.exp(m**(power1) * exponent) * m**(power2) * scale 

    

### External background model by A. Lindote (Jan 2023)def EresRateNormalisation(rate_array, eres):
def EresRateNormalisation(rate_array, eres):
    """
    Normalises the gamma rate for a particular energy resolution value (within
    the ones studied/supported). For that, need to remove the /keV dependency 
    from the ROI width first. If the energy resolution is not supported, don't 
    change the normalisation.
    """
    norm = 1
    Q = 2457.8 # keV
    ROI = 2*(Q*eres/100)   # ROI width (in keV) of the energy resolution being used
    ROI_0 =  2*(Q*1.0/100) # ROI with (in keV) with 1% energy resolution
    # scaling factors for the respective energy resolutions, from LZ sims
    resolution = [0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.64, 0.8, 1.0, 2.0]
    scale = [0.001, 0.043, 0.178, 0.486, 0.686, 0.804, 0.898, 0.957, 1.0, 1.075]
    # get an interpolating function
    f = ip.interp1d(resolution, scale, kind='cubic')
    if eres<0.1:
        norm = 0
        print('WARNING: Energy resolution too low')
    elif eres>2.0:
        norm = 1.075
        print('WARNING: Energy resolution too high')
    else:
        norm = f(eres)
#    print('For an energy resolution of ', eres, ' the normalization will be ', norm)
    return rate_array*ROI_0*norm/ROI

def Extr_Bkg(m, Model = 'XLZD_1:1', LXeMass = 80, eres = 1.0, gamma_frac = 1.0):
    if Model == 'DARWIN':
        popt_DARWIN = np.asarray([6.57488161e-02, 1.12434327e+00, 1.96618120e+00, 5.27331167e-05]) ## 4 params
        return 1 / 0.15 * Bkg_fit_func(m, *popt_DARWIN) * gamma_frac
    if Model == 'XLZD': 
        Model = 'XLZD_1:1'
        print('WARNING: Old Model String used, converted to _1:1')
    if Model == 'XLZD_1:1':
        #filename = './ExtBackground_files/ext_gammas_' + str(LXeMass) + 't.npz'       
        filename = './ExtBackground_files/ext_gammas_' + str(LXeMass) + 't_15cmRFR.npz'       
        npzfile = np.load(filename)
        mass_array, rate_array = npzfile['mass'], npzfile['rate']
        npzfile.close()   
        rate_array = EresRateNormalisation(rate_array, eres) # normalise rate from 1% to the given energy resolution
        rate_array *= 1000*365  # normalise to events / ton / year / keV
        rate_array *= gamma_frac / 0.10  # normalise to no SS/MS discrimination and rescaled gamma flux
        func = ip.interp1d(mass_array, rate_array, kind='linear'
                           #, bounds_error = False, fill_value="extrapolate"
                          )
        return(func(m))
    if Model == 'XLZD_shallow' and LXeMass == 40:
#        filename = './ExtBackground_files/ext_gammas_' + str(LXeMass) + 't_shallow.npz'       
        filename = './ExtBackground_files/ext_gammas_' + str(LXeMass) + 't_shallow_15cmRFR.npz'       
        npzfile = np.load(filename)
        mass_array, rate_array = npzfile['mass'], npzfile['rate']
        npzfile.close()   
        rate_array = EresRateNormalisation(rate_array, eres) # normalise rate from 1% to the given energy resolution
        rate_array *= 1000*365  # normalise to events / ton / year / keV
        rate_array *= gamma_frac / 0.10  # normalise to no SS/MS discrimination and rescaled gamma flux
        func = ip.interp1d(mass_array, rate_array, kind='cubic'
                           #, bounds_error = False, fill_value="extrapolate"
                          )
        return(func(m))   
    
    if Model == 'XLZD_tall' and LXeMass == 80:
        #filename = './ExtBackground_files/ext_gammas_' + str(LXeMass) + 't_tall.npz'       
        filename = './ExtBackground_files/ext_gammas_' + str(LXeMass) + 't_tall_15cmRFR.npz'       
        npzfile = np.load(filename)
        mass_array, rate_array = npzfile['mass'], npzfile['rate']
        npzfile.close()   
        rate_array = EresRateNormalisation(rate_array, eres) # normalise rate from 1% to the given energy resolution
        rate_array *= 1000*365  # normalise to events / ton / year / keV
        rate_array *= gamma_frac / 0.10  # normalise to no SS/MS discrimination and rescaled gamma flux
        func = ip.interp1d(mass_array, rate_array, kind='cubic'
                           #, bounds_error = False, fill_value="extrapolate"
                          )
        return(func(m))   
    else: 
        print('No valid Background model')
        return 0
    
def Extr_Bkg_shell(m_lower, m_higher, Model = 'XLZD_1:1', LXeMass = 80, eres = 1.0, gamma_frac = 1.0):
    Mean_Bkg = ((m_higher * Extr_Bkg(m_higher, Model = Model, LXeMass = LXeMass, eres=eres, gamma_frac=gamma_frac)) 
                - (m_lower * Extr_Bkg(m_lower, Model = Model, LXeMass = LXeMass, eres=eres, gamma_frac=gamma_frac))
               ) / (m_higher-m_lower)
    return Mean_Bkg



### Sensitivity Calculation

### Relying on the figure of merit approach
# for 90% C.L. Exclusion Limit and for 3sigma (99.73% C.L.) Discovery Potential 
   #S3σ(B) denotes the Poisson signal expectation at which 50% of the measurements in an 
   #ensemble of identical experiments would report a 3σ positive fluctuation above B.        

def sensitivity(mass,bkg,time, ROI_width, enrichment = 0, metric = 'Excl_90p'):
    ia = enrichment/100
    if isinstance(mass, list):
        mass = mass[1]-mass[0]
    if metric == 'Excl_90p':
        return np.log(2)*(eff*ia*N_a/1.64/M_Xe)*(np.sqrt(mass*time/bkg/ROI_width))
    elif metric == 'Disc_3sigma': 
        S = S_3sigma(bkg *mass *time * ROI_width)
        return np.log(2)*(N_a*ia/M_Xe)* eff* mass *time / S
    else: 
        print('Sensitivity metric not defined')
        return 0

def S_3sigma(bkg, CL = 0.9973):
    if isinstance(bkg, int) or isinstance(bkg, float): bkg = [bkg]
    S_list = []
    for B in bkg: 
        if B > 200: S_list.append(np.inf) 
        else:
            for C in np.arange(0,np.max([100, 2*B]), 0.01):
                if gammaincc(C+1, B) > CL: break
            for S in np.arange(0,100, 0.01):
                if (1-gammaincc(C +1 , S+B)) > 0.5: 
                    S_list.append(S) 
                    break
    if len(bkg) == 1: return S_list[0]
    else: return S_list
    
    

   
def Sensitivity_curve(FVmass_range, signal_acc, beta_acc, gamma_acc, eres = 1.0, 
                      gamma_args = ['XLZD', 80, 1.0], beta_args = ["LNGS", 99.9],  
                      exposuretime = 10, output= "short", enrichment = 8.9, flux = 4.6,
                      metric = "Excl_90p" ):
    
                      # gamma_args = [DetectorModel, LXeMass, gamma_frac], beta_args = [Lab, BiPoTagging_eff]
    
    ROI_width = round(2 * eres / 100 * E_Qbb, 1)
    Sensitivity = signal_acc * sensitivity(FVmass_range, 
                                           gamma_acc * Extr_Bkg(FVmass_range, Model = gamma_args[0], LXeMass = gamma_args[1], 
                                           eres = eres, gamma_frac = gamma_args[2]) +
                                           beta_acc * Intr_Bkg(*beta_args, enrichment, flux)[0], time = exposuretime, 
                                           ROI_width = ROI_width, enrichment = enrichment, metric = metric)
    
    max_index = np.argmax(Sensitivity)
    #print('Initial sensitivity in optimal fiducial:', FVmass_range[max_index], 't', Sensitivity[max_index],'yr')
    FV_Points_list = [FVmass_range[max_index]]
    Sensitivities_all = [Sensitivity]
    Sensitivity_summed =[Sensitivity[max_index]]
    """
    # Some debugging info
    Internal = beta_acc * Intr_Bkg(*beta_args, enrichment, flux)
    bkg = gamma_acc * Extr_Bkg(FVmass_range, Model = gamma_args[0], LXeMass = gamma_args[1], eres = eres, 
                               gamma_frac = gamma_args[2]) + Internal[0]
    print(beta_args[0], gamma_args[1], 'Total background: ', bkg[max_index])
    print('B-8: ', Internal[1], 'Rn-222: ', Internal[2], 'Xe-136: ', Internal[3], 'Xe-137: ', Internal[4])
    print('Fiducial mass: ', FV_Points_list, ' Sensitivity: ', Sensitivity_summed)
    """
    
    for i in range(10): 
        Sensitivity_here = signal_acc * sensitivity([FVmass_range[max_index], FVmass_range[max_index+1:]], 
                                                    gamma_acc * Extr_Bkg_shell(FVmass_range[max_index], 
                                                                               FVmass_range[max_index+1:], 
                                                                               Model = gamma_args[0], LXeMass = gamma_args[1],
                                                                               eres = eres, gamma_frac = gamma_args[2]) +
                                                    beta_acc * Intr_Bkg(*beta_args, enrichment, flux)[0], time = exposuretime, 
                                                    ROI_width = ROI_width, enrichment = enrichment, metric = metric)
        if len(Sensitivity_here) < 2: break
        Sensitivity_here = np.concatenate((np.array([0 for i in range(max_index)]), Sensitivity_here, [Sensitivity_here[-1]]))
        Sensitivities_all.append(Sensitivity_here)
        max_index = np.argmax(Sensitivity_here)
        FV_Points_list.append(FVmass_range[max_index])
        Sensitivity_summed.append(np.sqrt(Sensitivity_summed[-1]**2 + np.max(Sensitivity_here)**2))
        #print('Added sensitivity:',FVmass_range[max_index],'t',Sensitivity_here[max_index],'yr',Sensitivity_summed[-1])

    if output == "long": return([Sensitivity, Sensitivity_summed, FV_Points_list, Sensitivities_all])
    else: return ([Sensitivity, Sensitivity_summed])
       
        
def return_Sensitivity(idx0, idx1, idx2, idx3, idx4, idx5, metric):
    gamma_frac,BiPo_eff, LXeMass = idx0 / 100, idx1, idx3
    if idx4 == 1 or idx4 == idx3: Model = 'XLZD_1:1'
    elif (idx4 == 60 and idx3 == 40): Model = 'XLZD_shallow'
    elif ((idx4 == 60) and idx3 == 80): Model = 'XLZD_tall'
    else: 
        print("WARNING: Invalid Background Model")
        return 0,0
    
    if idx2 == 2: signal_acc, beta_acc, gamma_acc, eres = 0.853, 0.767, 0.098, 0.65  ### XLZD new baseline (3mm in z-only)
    elif idx2 == 1: signal_acc, beta_acc, gamma_acc, eres = .825, 0.739, 0.070, 0.60 ### XLZD new progressive (2mm in z-only)
        
    if idx5 == 1: Lab = 'Kamioka' 
    elif idx5 == 2: Lab = 'LNGS'
    elif idx5 == 3: Lab = 'Boulby-1300'
    elif idx5 == 4: Lab = 'SURF'
    elif idx5 == 5: Lab = 'SNOLab' 
       
    gamma_args = [Model, LXeMass, gamma_frac]
    beta_args = [Lab, BiPo_eff]  
    FVmass_range = np.arange(0.1, 0.8*LXeMass, 0.1)
    [Sensitivity, Sensitivity_summed] = Sensitivity_curve(FVmass_range, signal_acc, beta_acc, gamma_acc, eres=eres, gamma_args = gamma_args , beta_args = beta_args,  exposuretime =  10, metric = metric)
    if np.max(Sensitivity) == Sensitivity[-1]:
        print("WARNING: Sensitivity maximizes at largest allowed FV")
    return np.max(Sensitivity) / 1e27, np.max(Sensitivity_summed) / 1e27, 
    
    
    
### Conversion functions m_bb -> T_half


def m_bb(T_half, M_0nu = NME_Xe136[0], m_e = m_e, G_0nu = G_0nu_Xe136):
    return m_e / np.sqrt(G_0nu * T_half * M_0nu**2) / gA**2 * 1000


def T_half(m_bb, M_0nu = NME_Xe136[0], m_e = m_e, G_0nu = G_0nu_Xe136):
    return (m_e / m_bb * 1000)**2 / G_0nu / M_0nu**2 / gA**4
    



