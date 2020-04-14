
########################################################################################################################
########################################################################################################################

#                           				TFG dirigido por NICOLÁS CARDIEL 
# 
#                   				ELAVORACIÓN DE UN MODELO NUMÉRICO DE INTERIOR ESTELAR
# 
#                              				Autor: LUIS ABALO RODRÍGUEZ
# 
#                              				FACULTAD DE CIENCIAS FÍSICAS
#                           				UNIVERSIDAD COMPLUTENSE DE MADRID
# 
#                               				CURSO ACADÉMICO 2019-2020
#
########################################################################################################################
########################################################################################################################

from A_parametros_iniciales import * 												# Importamos los parámetros iniciales
from numpy import sqrt 													  # Importamos la raíz cuadrada del paquete numpy

# FUNCIONES PARA CALCULAR LA GENERACIÓN DE ENERGÍA
#------------------------------------------------------------------------------------------------------------------------
def pp_fun(T):
	'''
	Estudio de la generación de energía para la cadena p-p
	Input: Temperatura que se va a estudiar
	Output: Valor de la generación de energía
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
	Estudio de la generación de energía para el ciclo CNO
	Input: Temperatura que se va a estudiar
	Output: Valor de la generación de energía
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
	Selección del valor mayor de producción de energía: cadena p-p o ciclo CNO
	Input: Temperatura que se va a estudiar
	Output: Valor de la generación de energía
	'''
	global  nu, pp_fun, CN_fun 

	eps = max(pp_fun(T), CN_fun(T))
	return eps

def ciclo(T):
	'''
	Selección del ciclo de producción de energía mayor
	Input: Temperatura que se va a estudiar
	Output: Valor de la generación de energía
	'''
	global generacion_energia, pp_fun

	eps = generacion_energia(T)
	eps1 = pp_fun(T)

	if eps < 0.017: 						  # Por debajo de 0.017, se considera que la generación de la energía es nula
		return "--"
	elif eps == eps1:
		return "PP"
	else:
		return "CN"
#------------------------------------------------------------------------------------------------------------------------


# FUNCIÓN PARA CALCULAR EL ERROR RELATIVO
#------------------------------------------------------------------------------------------------------------------------
Erel_max = 0.0001
def error_relativo(calculada, estimada):
	'''
	Estudio del error relativo entre el valor calculado y el valor estimado
	Input: valor calculado y valor estimado
	Output: error relativo

	'''
	return abs(calculada - estimada) / calculada
#------------------------------------------------------------------------------------------------------------------------


# ECUACIONES DIFERENCIALES FUNDAMENTALES: CASO RADIACTIVO DE LA TABLA 3
#------------------------------------------------------------------------------------------------------------------------
# Ecuación (18): Masa
def dM(r, T, P):
	'''
	Estudio de la ecuación diferencial de la masa
	Input: radio, temperatura y presión 
	Output: valor para la masa que se acumula en la lista f_i
	'''
	global mu
	C_m = 0.01523 * mu
	return C_m * (P / T) * (r ** 2)

# Ecuación (19): Presión
def dP(r, T, P, M): 
	'''
	Estudio de la ecuación diferencial de la presión
	Input: radio, temperatura, presión y masa
	Output: valor para la presión que se acumula en la lista f_i
	'''
	global mu
	C_p = 8.084 * mu
	return  - C_p * (P / T) * M / (r ** 2)

# Ecuación (20): Temperatura
def dT(r, T, P, L):
	'''
	Estudio de la ecuación diferencial de la temperatura
	Input: radio, temperatura, presión y luminosidad
	Output: valor para la temperatura que se acumula en la lista f_i
	'''
	global mu, Z, X
	C_t = 0.01679 * Z * (1 + X) * (mu) ** 2
	return  - C_t * (P ** 2) * L * (T ** (-17 / 2)) * (r ** (-2))

# Ecuación (21): Luminosidad
def dL(r, T, P):
	'''
	Estudio de la ecuación diferencial de la luminosidad
	Input: radio, temperatura y presión 
	Output: valor para la luminosidad que se acumula en la lista f_i
	'''
	global mu, generacion_energia
	eps = generacion_energia(T) 

	if eps < 0.017: 							   # Se considera que se genera energía por encima de 0.017 de producción
		C = 0
	else:
		C = 0.01845 
	return C * eps * (P ** 2) * (T ** -2) * (r ** 2) * (mu ** 2)
#------------------------------------------------------------------------------------------------------------------------

# ECUACIONES PARA LAS CAPAS EXTERIORES: TABLA 4
#------------------------------------------------------------------------------------------------------------------------

# Ecuación (35): Temperatura
def T(r, R_tot): 
	'''
	Estudio de la temperatura en una capa determinada
	Input: radio de la capa y radio total
	Output: valor de la temperatura de la capa
	'''
	global mu, M_tot 
	A_1 = 1.9022 * mu * M_tot 
	return A_1 * ((1 / r) - (1 / R_tot))

# Ecuación (36): Presión
def P(T, L_tot): 
	'''
	Estudio de la presión en una capa determinada
	Input: temperatura de la capa y luminosidad total 
	Output: valor de la presión de la capa
	'''
	global M_tot, mu, Z, X
	A_2 = 10.645 * sqrt((M_tot / (mu * Z * (1 + X) * L_tot)))
	return A_2 * T ** (17/4)
#------------------------------------------------------------------------------------------------------------------------

# PRESIÓN Y TEMPERATURA ESTIMADAS
#------------------------------------------------------------------------------------------------------------------------
# Presión estimada
def presion_estimada(P_pres, f_i, h, i): # P_pres = integracion[i][4]
	'''
	Estudio de la presión estimada en una capa determinada 
	Input: presión, lista de las primeras derivadas, paso y número de la capa 
	Output: Valor de la presión estimada
	'''
	fi = f_i[i][4] 
	fi_1 = f_i[i-1][4]
	fi_2 = f_i[i-2][4]

	delta1 = h * fi - h * fi_1
	delta2 = h * fi - 2 * h * fi_1 + h * fi_2
	return P_pres + h * fi + (1/2) * delta1 + (5/12) * delta2

# Temperatura estimada
def temperatura_estimada(T_temp, f_i, h, i): #T_temp = integracion[i][5]
	'''
	Estudio de la temperatura estimada en una capa determinada 
	Input: temperatura, lista de las primeras derivadas, paso y número de la capa 
	Output: Valor estimado de la temperatura 
	'''
	fi = f_i[i][5]
	fi_1 = f_i[i-1][5]

	delta1 = h * fi - h * fi_1
	return T_temp + h * fi + (1/2) * delta1
#------------------------------------------------------------------------------------------------------------------------

# PRESIÓN, TEMPERATURA, MASA Y LUMINOSIDAD CALCULADAS
#------------------------------------------------------------------------------------------------------------------------
# Presión, temperatura y masa (PTM) calculadas 
def PTM_calculada(ptm, dptm, fi, h, i): # P = integracion[i][4], T = integracion[i][5], M = integracion[i][7]
# P: fi = f_i[i][4], T: fi = f_i[i][5], M: fi = f_i[i][7]    
	'''
	Estudio de la PTM calculada de una capa determinada 
	Input: PTM, última primera derivada de PTM, lista de
	primeras derivadas de PTM, paso y número de la capa
	Output: Valor de la PTM calculada 
	'''
	delta1 = h * dptm - h * fi
	return ptm + h * dptm - (1/2) * delta1

# Luminosidad calculada
def L_calculada(L, dLu, f_i, h, i): # L = L_tot
	'''
	Estudio de la luminosidad calculada de una determinada capa 
	Input: luminosidad, última primera derivada y lista de pri-
	meras derivadas de la luminosidad, paso y número de la capa 
	Output: Valor de la luminosidad calculada 
	'''
	fi = f_i[i][6]
	fi_1 = f_i[i-1][6]

	delta1 = h * dLu - h * fi
	delta2 = h * dLu - 2 *  h * fi + h * fi_1
	return L + h * dLu - (1/2) * delta1 - (1/12) * delta2
#------------------------------------------------------------------------------------------------------------------------





# 
#------------------------------------------------------------------------------------------------------------------------
# We are going to define n+1:
def n_fun(Tcal, Pcal, dPr, dTe):
	'''
	Study the value of the parameter that gives us information about the radioactive or convective behavior
	Input: Calculated Temperature and Pressure of the shell, pressure and temperature first derivative of the shell
	Output: Value of the n-parameter 
	'''
	return (Tcal / Pcal) * (dPr / dTe)
#------------------------------------------------------------------------------------------------------------------------

# 
#------------------------------------------------------------------------------------------------------------------------
def interpolation(y0,y1,n0,n1):
	'''
	Generates a linear interpolation
	Input: two pairs of values (x,y) to make the linear intepolation
	Output: The value interpolated for n=2.5
	'''
	n = 2.5
	return y0 + ((y1 - y0) / (n1 - n0)) * (n - n0)
#------------------------------------------------------------------------------------------------------------------------

#---------------------------FUNCTIONS FOR PHYICAL PARAMETERS (CONVECTIVE PHASE)--------------------------
# We are going to define the equations to obtain the parameters in the convective core
# First the functions for the first derivative of the parameters

# ECUACIONES FUNDAMENTALES: CASO CONVECTIVO DE LA TABLA 3
#------------------------------------------------------------------------------------------------------------------------
C_m = 0.01523*mu
C_p = 8.084*mu
C_t = 0.01679*Z*(1+X)*(mu)**2
C_t_c = 3.234 * mu

# Ecuación (22)
def dMc(T,r, K): 
	'''
	Study of the mass first derivative
	Input: Temperature, Radius and Polytrope model constant
	Output: Value of fi_M
	'''
	global mu
	C_m = 0.01523 * mu
	return C_m * K * (T ** 1.5) * (r ** 2)

# Ecuación (23)
def dPc(T,r, M, K):
	'''
	Study of the pressure first derivative
	Input: Temperature, Radius, Mass and Polytrope model constant
	Output: Value of fi_P
	'''
	global mu
	C_p = 8.084 * mu
	return - C_p * K * (T ** 1.5) * M * (r ** (-2))

# Ecuación (24)
def dLc(r, T, K):
	'''
	Study of the luminosity first derivative
	Input: Temperature, Radius and the Polytrope model constant
	Output: Value of fi_L
	'''
	global mu, generacion_energia
	eps = generacion_energia(T)   
	if eps < 0.017: # Para que dL = 0 (considero que se genera energía por encima de 0.017)
		C = 0
	else:
		C = 0.01845 
	return C * eps * (K ** 2) * (T ** 3) * (r ** 2) * (mu ** 2)

# Ecuación (25)
def dTc(r, M):
	'''
	Study of the temperature first derivative
	Input: Mass and Radius
	Output: Value of fi_T
	'''
	global mu
	C_t = 3.234 * mu
	return - C_t * M / (r ** 2)
#------------------------------------------------------------------------------------------------------------------------

# ECUACIONES FUNDAMENTALES: CASO RADIACTIVO DE LA TABLA 5
#------------------------------------------------------------------------------------------------------------------------
# Then we define the functions to calculate the T(r) and P(T)
# We will use this functions for the first three shells

# Ecuación (43)
# Mass(r)
def Mcentro(r, K, T_c):
	'''
	Study of the mass in a certain shell
	Input: Radius ,Polytrope model constant and center Temperature
	Output: Value of M(r)
	'''
	global mu
	return 0.005077 * mu * K * (T_c ** 1.5) * (r ** 3)

# Ecuación (44)
# Luminosity(r)
def Lcentro(r, K, T_c):
	'''
	Study of the luminosity in a certain shell
	Input: Radius, the Polytrope model constant and center Temperature
	Output: Value of L(r)
	'''
	global mu, generacion_energia
	eps = generacion_energia(T_c)
	return 0.006150 * eps * (K ** 2) * (T_c ** 3) * (r ** 3) * (mu ** 2)

# Ecuación (45)
# Temperature(r)
def Tcentro(r, K, T_c): 
	'''
	Study of the temperature in a certain shell
	Input: Radius, the Polytrope model constant and center Temperature
	Output: Value of T(r)
	'''
	global mu
	return  T_c - 0.008207 * (mu ** 2) * K * (T_c ** 1.5) * (r ** 2)

# Ecuación (46)
# Pressure(r)
def Pcentro(r, K, T_centro):
	'''
	Study of the Mass in a certain shell
	Input: Radius, Temperature and the Polytrope model constant 
	Output: Value of P(r)
	'''    
	return K * (T_centro ** 2.5)
#------------------------------------------------------------------------------------------------------------------------

# 
#------------------------------------------------------------------------------------------------------------------------
# We can define a function to calculate the total relative error
def TotalRelEror(P_Interp, P_Interp_cen, T_Interp, T_Interp_cen, L_Interp, L_Interp_cen,
 M_Interp, M_Interp_cen):
	'''
	Study the total relative error in the frontier between the radioative and convective part
	Inputs: values of pressure, temperature, luminosity and mass in the frontier, in each part
	Output: value of the total relative error
	'''    
	a = ((P_Interp - P_Interp_cen) / P_Interp) ** 2
	b = ((T_Interp - T_Interp_cen) / T_Interp) ** 2
	c = ((L_Interp - L_Interp_cen) / L_Interp) ** 2
	d = ((M_Interp - M_Interp_cen) / M_Interp) ** 2
	return sqrt(a + b + c + d)
#------------------------------------------------------------------------------------------------------------------------ 
