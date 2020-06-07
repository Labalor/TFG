
# -*- coding: utf-8 -*-
########################################################################################################################
########################################################################################################################

#                           				TFG dirigido por NICOLÁS CARDIEL 
# 
#                   				ELAVORACIÓN DE UN MODELO NUMÉRICO DE INTERIOR ESTELAR
# 
#                              				Autor: LUIS ABALO RODRÍGUEZ
# 
#                              				FACULTAD DE CIENCIAS FÍSICAS
#                           			 UNIVERSIDAD COMPLUTENSE DE MADRID
# 
#                               				CURSO ACADÉMICO 2019-2020
#
########################################################################################################################
########################################################################################################################

#------------------------------------------------------------------------------------------------------------------------
# PAQUETES NECESARIOS
#------------------------------------------------------------------------------------------------------------------------
from A_parametros_iniciales import * 	# Initial parameters are imported
from numpy import sqrt 					# The square root of the numpy package is imported
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# 									FUNCTIONS TO CALCULATE POWER GENERATION
#------------------------------------------------------------------------------------------------------------------------
def pp_fun(T):
	'''   
    Study the value of energy generation for a certain temperature for the proton-proton chain.

    It turns out that the easiest nuclear reaction is the reaction between a deuteron 
    and a proton. This reaction can happen with temperatures around 2 million degrees. 
    
    Parameters
    ----------
    T : float
        The tempeture that we want to study.

    Returns
    -------
    float
        Value of the energy rate production.
    '''
	global X

	T *= 10

	if T >= 4.0 and T <= 6.0:
		eps_1_pp = 10**(-6.84) 
		nu_pp = 6.0
	elif T >= 6.0 and T <= 9.5:
		eps_1_pp = 10**(-6.04) 
		nu_pp = 5.0
	elif T >= 9.5 and T <= 12.0:
		eps_1_pp = 10**(-5.56) 
		nu_pp = 4.5
	elif T >= 12.0 and T <= 16.5:
		eps_1_pp = 10**(-5.02) 
		nu_pp = 4.0
	elif T >= 16.5 and T <= 24.0:
		eps_1_pp = 10**(-4.40) 
		nu_pp = 3.5
	else:
		eps_1_pp = 0
		nu_pp = 0

	return eps_1_pp * (X**2) * (T**nu_pp)

def CN_fun(T):
	'''   
    Study the value of energy generation for a certain temperature for the carbon-nitrogen cycle.

    In the CNO cycle, four protons fuse, using carbon, nitrogen, and oxygen isotopes as catalysts, 
    to produce one alpha particle, two positrons and two electron neutrinos. Although there are various 
    paths and catalysts involved in the CNO cycles, all these cycles have the same net result.
    
    Parameters
    ----------
    T : float
        The tempeture that we want to study.

    Returns
    -------
    float
        Value of the energy rate production.
    '''
	global X, Z

	T *= 10

	if T >= 12.0 and T <= 16.0:
		eps_1_CN = 10**(-22.2) 
		nu_CN = 20.0
	elif T >= 16.0 and T <= 22.5:
		eps_1_CN = 10**(-19.8) 
		nu_CN = 18.0
	elif T >= 22.5 and T <= 27.5:
		eps_1_CN = 10**(-17.1) 
		nu_CN = 16.0
	elif T >= 27.5 and T <= 36.0:
		eps_1_CN = 10**(-15.6) 
		nu_CN = 15.0
	elif T >= 36.0 and T <= 50.0:
		eps_1_CN = 10**(-12.5) 
		nu_CN = 13.0
	else:
		eps_1_CN = 0
		nu_CN = 0 

	return eps_1_CN * (X*Z/3) * (T**nu_CN)

def generacion_energia(T): 
	'''   
    Selection of the highest energy production value: p-p chain or CNO cycle.
    
    Parameters
    ----------
    T : float
        The tempeture that we want to study.

    Returns
    -------
    float
        Value of the energy rate production.
    '''
	global  nu, pp_fun, CN_fun 

	eps = max(pp_fun(T), CN_fun(T))
	return eps

def ciclo(T):
	'''   
    Selection of the highest energy production method.
    
    Parameters
    ----------
    T : float
        The tempeture that we want to study.

    Returns
    -------
    string
        Name of the method of the energy rate production.
    '''
	global generacion_energia, pp_fun

	eps = generacion_energia(T)
	eps1 = pp_fun(T)

# Below 0.017, power generation is considered to be nil.
	if eps < 0.017: 						  
		return "--"
	elif eps == eps1:
		return "PP"
	else:
		return "CN"
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# FUNDAMENTAL DIFFERENTIAL EQUATIONS: RADIOACTIVE CASE OF TABLE 3
#------------------------------------------------------------------------------------------------------------------------
# Equation (18): Mass
def dM(r, T, P):
	'''   
    Differential Equation of Mass Study.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.
    T : float
        The tempeture of the layer to study.
	P : float
        The pressure of the layer to study.

    Returns
    -------
    float
        Value for the mass that accumulates in the list f_i.
    '''
	global mu
	C_m = 0.01523 * mu
	return C_m * (P / T) * (r ** 2)

# Equation (19): Pressure
def dP(r, T, P, M): 
	'''   
    Differential Equation of Pressure Study.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.
    T : float
        The tempeture of the layer to study.
	P : float
        The pressure of the layer to study.
    M : float
        The pressure of the layer to study.

    Returns
    -------
    float
        Value for the pressure that accumulates in the list f_i.
    '''
	global mu
	C_p = 8.084 * mu
	return  - C_p * (P / T) * M / (r ** 2)

# Equation (20): Temperature
def dT(r, T, P, L):
	'''   
    Differential Equation of Temperature Study.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.
    T : float
        The tempeture of the layer to study.
	P : float
        The pressure of the layer to study.
    L : float
        The luminosity of the layer to study.

    Returns
    -------
    float
        Value for the temperature that accumulates in the list f_i.
    '''
	global mu, Z, X
	C_t = 0.01679 * Z * (1 + X) * (mu) ** 2
	return  - C_t * (P ** 2) * L * (T ** (-17 / 2)) * (r ** (-2))

# Equation (21): Luminosity
def dL(r, T, P):
	'''   
    Differential Equation of Luminosity Study.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.
    T : float
        The tempeture of the layer to study.
	P : float
        The pressure of the layer to study.

    Returns
    -------
    float
        Value for the luminosity that accumulates in the list f_i.
    '''
	global mu, generacion_energia
	eps = generacion_energia(T) 

# 	Below 0.017, power generation is considered to be nil.
	if eps < 0.017: 							   
		C = 0
	else:
		C = 0.01845 
	return C * eps * (P ** 2) * (T ** -2) * (r ** 2) * (mu ** 2)
#------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------
# OUTER LAYER EQUATIONS: TABLE 4
#------------------------------------------------------------------------------------------------------------------------
# Equation (35): Temperature
def T(r, R_tot): 
	'''   
    Study of the temperature in a given layer.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.
    R_tot : float
        The initial total radius. 

    Returns
    -------
    float
        Layer temperature value.
    '''
	global mu, M_tot 
	A_1 = 1.9022 * mu * M_tot 
	return A_1 * ((1 / r) - (1 / R_tot))

# Equation (36): Pressure
def P(T, L_tot): 
	'''   
    Study of the temperature in a given layer.
    
    Parameters
    ----------
    T : float
        he tempeture of the layer to study.
    L_tot : float
        The initial total luminosity. 

    Returns
    -------
    float
        Layer temperature value.
    '''
	global M_tot, mu, Z, X
	A_2 = 10.645 * sqrt((M_tot / (mu * Z * (1 + X) * L_tot)))
	return A_2 * T ** (17/4)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# ESTIMATED PRESSURE AND TEMPERATURE 
#------------------------------------------------------------------------------------------------------------------------
# Estimated pressure
def presion_estimada(P_pres, f_i, h, i): # P_pres = integracion[i][4]
	'''   
    Study of the pressure in a given layer.
    
    Parameters
    ----------
    P_pres : float
        The pressure to study.
    f_i : array
        List of first derivatives 
    h : float
        The step
    i : float
        The number of the layer	

    Returns
    -------
    float
        Estimated pressure value
    '''
	fi = f_i[i][4] 
	fi_1 = f_i[i-1][4]
	fi_2 = f_i[i-2][4]

	delta1 = h * fi - h * fi_1
	delta2 = h * fi - 2 * h * fi_1 + h * fi_2
	return P_pres + h * fi + (1/2) * delta1 + (5/12) * delta2

# Estimated temperature
def temperatura_estimada(T_temp, f_i, h, i): # T_temp = integracion[i][5]
	'''   
    Study of the temperature in a given layer.
    
    Parameters
    ----------
    T_temp : float
        The temperature to study.
    f_i : array
        List of first derivatives. 
    h : float
        The step. 
    i : float
        The number of the layer	

    Returns
    -------
    float
        Estimated temperature value.
    '''
	fi = f_i[i][5]
	fi_1 = f_i[i-1][5]

	delta1 = h * fi - h * fi_1
	return T_temp + h * fi + (1/2) * delta1
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# CALCULATED PRESSURE, TEMPERATURE, MASS AND LUMINOSITY 
#------------------------------------------------------------------------------------------------------------------------
# Calculated pressure, temperature and mass (PTM) 
def PTM_calculada(ptm, dptm, fi, h, i): # P = integracion[i][4], T = integracion[i][5], M = integracion[i][7]
# P: fi = f_i[i][4], T: fi = f_i[i][5], M: fi = f_i[i][7]
	'''   
    Study of the pressure, temperature or mass in a given layer.
    
    Parameters
    ----------
    ptm : float
        Value of PTM to study.
    dptm : float
        Last first derivative of PTM. 
    fi : array
        List of first derivatives of PTM. 
    h : float
        The step. 
    i : float
        The number of the layer.	

    Returns
    -------
    float
        Calculated PTM value.
    '''    
	delta1 = h * dptm - h * fi
	return ptm + h * dptm - (1/2) * delta1

# Calculated luminosity
def L_calculada(L, dLu, f_i, h, i): # L = L_tot
	'''   
    Study of the calculated temperature in a given layer.
    
    Parameters
    ----------
    L : float
        The luminosity of the layer to study.
    dLu : float
        Last first derivative. 
    f_i : array
        List of first derivatives of Luminosity. 
    h : float
        The step. 
    i : float
        The number of the layer.	

    Returns
    -------
    float
        Calculated luminosity value.
    ''' 
	fi = f_i[i][6]
	fi_1 = f_i[i-1][6]

	delta1 = h * dLu - h * fi
	delta2 = h * dLu - 2 *  h * fi + h * fi_1
	return L + h * dLu - (1/2) * delta1 - (1/12) * delta2
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE THE RELATIVE ERROR
#------------------------------------------------------------------------------------------------------------------------
Erel_max = 0.0001
def error_relativo(calculada, estimada):
	'''   
    Study of the relative error between the calculated value and the estimated value.
    
    Parameters
    ----------
    calculada : float
        Calculated value.
    estimada : float
        Estimated value. 

    Returns
    -------
    float
        Value of the relative error.
    ''' 
	return abs(calculada - estimada) / calculada
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# n+1 PARAMETER
#------------------------------------------------------------------------------------------------------------------------
def n_fun(Tcal, Pcal, dPr, dTe):
	'''   
    Study of the value of the parameter n+1.
    
    Parameters
    ----------
    Tcal : float
        Calculated temperature.
    Pcal : float
        Calculated pressure. 
    dPr : float
        Pressure of the first derivative.
    dTe : float
        Temperature of the first derivative. 

    Returns
    -------
    float
        Value of the n+1 parameter.
    ''' 
	return (Tcal / Pcal) * (dPr / dTe)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# INTERPOLATION
#------------------------------------------------------------------------------------------------------------------------
def interpolation(y0,y1,n0,n1):
	'''   
    Linear interpolation study.
    
    Parameters
    ----------
    y0 : float
        First value of Y.
    y1 : float
        Second value of Y. 
    n0 : float
        First value of X.
    n1 : float
        Second value of X. 

    Returns
    -------
    float
        The interpolated value for n+1 = 2.5.
    ''' 
	n = 2.5
	return y0 + ((y1 - y0) / (n1 - n0)) * (n - n0)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# FUNDAMENTAL DIFFERENTIAL EQUATIONS: CONVECTIVE CASE OF TABLE 3
#------------------------------------------------------------------------------------------------------------------------
# Equation (22): mass
def dMc(T,r, K): 
	'''   
    Differential Equation of Mass Study.
    
    Parameters
    ----------
    T : float
        The temperature of the layer to study.
    r : float
        The radius of the layer to study.
    K : float
        Polytrope constant.

    Returns
    -------
    float
        Value for the mass that accumulates in the list f_i.
    ''' 
	global mu
	C_m = 0.01523 * mu
	return C_m * K * (T ** 1.5) * (r ** 2)

# Equation (23): pressure
def dPc(T,r, M, K):
	'''   
    Differential Equation of Pressure Study.
    
    Parameters
    ----------
    T : float
        The temperature of the layer to study.
    r : float
        The radius of the layer to study. 
    M : float
        The mass of the layer to study.
    K : float
        Polytrope constant.

    Returns
    -------
    float
        Value for the pressure that accumulates in the list f_i.
    ''' 
	global mu
	C_p = 8.084 * mu
	return - C_p * K * (T ** 1.5) * M * (r ** (-2))

# Equation (24): luminosity
def dLc(r, T, K):
	'''   
    Differential Equation of Luminosity Study.
    
    Parameters
    ----------
    T : float
        The temperature of the layer to study.
    r : float
        The radius of the layer to study. 
    K : float
        Polytrope constant.

    Returns
    -------
    float
        Value for the luminosity that accumulates in the list f_i.
    '''
	global mu, generacion_energia
	eps = generacion_energia(T)   
# 	So that dL = 0 and, thus, dLc = 0 (I consider that energy is generated above 0.017).	
	if eps < 0.017: 
		C = 0
	else:
		C = 0.01845 
	return C * eps * (K ** 2) * (T ** 3) * (r ** 2) * (mu ** 2)

# Equation (25): temperature
def dTc(r, M):
	'''   
    Differential Equation of Temperature Study.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study. 
    M : float
        The mass of the layer to study.

    Returns
    -------
    float
        Value for the temperature that accumulates in the list f_i..
    '''
	global mu
	C_t = 3.234 * mu
	return - C_t * M / (r ** 2)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# CONVECTIVE CASE EQUATIONS: TABLE 5
#------------------------------------------------------------------------------------------------------------------------
# Equation (43): mass
def Mcentro(r, K, T_c):
	'''   
    Study of the mass in a determined central layer.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study. 
    K : float
        Polytrope constant.
    T_c : float
        The core temperature of the layer to study.

    Returns
    -------
    float
        Value of the mass of a central layer.
    '''
	global mu
	return 0.005077 * mu * K * (T_c ** 1.5) * (r ** 3)

# Equation (44): luminosity
def Lcentro(r, K, T_c):
	'''   
    Study of the luminosity in a determined central layer.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.  
    K : float
        Polytrope constant.
    T_c : float
        The core temperature of the layer to study.	

    Returns
    -------
    float
        Value of the luminosity of a central layer.
    '''
	global mu, generacion_energia
	eps = generacion_energia(T_c)
	return 0.006150 * eps * (K ** 2) * (T_c ** 3) * (r ** 3) * (mu ** 2)

# Equation (45): temperature
def Tcentro(r, K, T_c):
	'''   
    Study of the temperature in a determined central layer.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.  
    K : float
        Polytrope constant.
    T_c : float
        The core temperature of the layer to study.	

    Returns
    -------
    float
        Value of the temperature of a central layer.
    '''
	global mu
	return  T_c - 0.008207 * (mu ** 2) * K * (T_c ** 1.5) * (r ** 2)

# Equation (46): pressure
def Pcentro(r, K, T_centro):
	'''   
    Study of the pressure in a determined central layer.
    
    Parameters
    ----------
    r : float
        The radius of the layer to study.  
    K : float
        Polytrope constant.
    T_centro : float
        The core temperature of the layer to study.	

    Returns
    -------
    float
        Value of the pressure of a central layer.
    '''
	return K * (T_centro ** 2.5)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# TOTAL RELATIVE ERROR
#------------------------------------------------------------------------------------------------------------------------
def TotalRelEror(P_Interp, P_Interp_cen, T_Interp, T_Interp_cen, L_Interp, L_Interp_cen,
 M_Interp, M_Interp_cen):
	'''   
    Study of the relative total error in the radioactive-convective boundary.
    
    Parameters
    ----------
    P_Interp : float
        Pressure value interpolated in integration from outside. 
    P_Interp_cen :  float
        Pressure value interpolated in integration from inside.
    T_Interp : float
        Temperature value interpolated in integration from outside.
    T_Interp_cen : float
        Temperature value interpolated in integration from inside. 
    L_Interp : float
        Luminosity value interpolated in integration from outside.
    L_Interp_cen : float
        Luminosity value interpolated in integration from inside.
    M_Interp : float
        Mass value interpolated in integration from outside. 
    M_Interp_cen : float
        Mass value interpolated in integration from inside.

    Returns
    -------
    float
        Total relative error value.
    '''
	a = ((P_Interp - P_Interp_cen) / P_Interp) ** 2
	b = ((T_Interp - T_Interp_cen) / T_Interp) ** 2
	c = ((L_Interp - L_Interp_cen) / L_Interp) ** 2
	d = ((M_Interp - M_Interp_cen) / M_Interp) ** 2
	return sqrt(a + b + c + d)
#------------------------------------------------------------------------------------------------------------------------ 
