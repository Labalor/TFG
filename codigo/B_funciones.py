
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

#------------------------------------------------------------------------------------------------------------------------
# PAQUETES NECESARIOS
#------------------------------------------------------------------------------------------------------------------------
from A_parametros_iniciales import * 	# Se importan los parámetros iniciales
from numpy import sqrt 					# Se importan la raíz cuadrada del paquete numpy
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# 									FUNCIONES PARA CALCULAR LA GENERACIÓN DE ENERGÍA
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

# Por debajo de 0.017, se considera que la generación de la energía es nula
	if eps < 0.017: 						  
		return "--"
	elif eps == eps1:
		return "PP"
	else:
		return "CN"
#------------------------------------------------------------------------------------------------------------------------

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

# 	Se considera que se genera energía por encima de 0.017 de producción
	if eps < 0.017: 							   
		C = 0
	else:
		C = 0.01845 
	return C * eps * (P ** 2) * (T ** -2) * (r ** 2) * (mu ** 2)
#------------------------------------------------------------------------------------------------------------------------


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

#------------------------------------------------------------------------------------------------------------------------
# PARÁMETRO n+1
#------------------------------------------------------------------------------------------------------------------------
def n_fun(Tcal, Pcal, dPr, dTe):
	'''
	Estudio del valor n+1
	Input: temperatura y presión calculadas, temperatura y presión de las primeras derivadas
	Output: valor de n+1
	'''
	return (Tcal / Pcal) * (dPr / dTe)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# INTERPOLACIÓN
#------------------------------------------------------------------------------------------------------------------------
def interpolation(y0,y1,n0,n1):
	'''
	Estudio de la interpolación lineal
	Input: un par de valores X-Y 
	Output: el valor interpolado para n+1 = 2.5 
	'''
	n = 2.5
	return y0 + ((y1 - y0) / (n1 - n0)) * (n - n0)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# ECUACIONES DIFERENCIALES FUNDAMENTALES: CASO CONVECTIVO DE LA TABLA 3
#------------------------------------------------------------------------------------------------------------------------
# Ecuación (22): masa
def dMc(T,r, K): 
	'''
	Estudio de la ecuación diferencial de la masa
	Input: temperatura, radio y constante del polítropo 
	Output: valor para la masa que se acumula en la lista f_i
	'''
	global mu
	C_m = 0.01523 * mu
	return C_m * K * (T ** 1.5) * (r ** 2)

# Ecuación (23): presión
def dPc(T,r, M, K):
	'''
	Estudio de la ecuación diferencial de la presión
	Input: temperatura, radio, masa y constante del polítropo 
	Output: valor para la presión que se acumula en la lista f_i
	'''
	global mu
	C_p = 8.084 * mu
	return - C_p * K * (T ** 1.5) * M * (r ** (-2))

# Ecuación (24): luminosidad
def dLc(r, T, K):
	'''
	Estudio de la ecuación diferencial de la luminosidad
	Input: radio, temperatura y constante del polítropo 
	Output: valor para la luminosidad que se acumula en la lista f_i
	'''
	global mu, generacion_energia
	eps = generacion_energia(T)   
# 	Para que dL = 0 y, así, dLc = 0 (considero que se genera energía por encima de 0.017)	
	if eps < 0.017: 
		C = 0
	else:
		C = 0.01845 
	return C * eps * (K ** 2) * (T ** 3) * (r ** 2) * (mu ** 2)

# Ecuación (25): temperatura
def dTc(r, M):
	'''
	Estudio de la ecuación diferencial de la temperatura
	Input: radio y masa 
	Output: valor para la temperatura que se acumula en la lista f_i
	'''
	global mu
	C_t = 3.234 * mu
	return - C_t * M / (r ** 2)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# ECUACIONES DEL CASO CONVECTIVO: TABLA 5
#------------------------------------------------------------------------------------------------------------------------
# Ecuación (43): masa
def Mcentro(r, K, T_c):
	'''
	Estudio de la masa en una capa central determinada
	Input: radio, constante de polítropo y temperatura central 
	Output: valor de la masa de una capa central
	'''
	global mu
	return 0.005077 * mu * K * (T_c ** 1.5) * (r ** 3)

# Ecuación (44): luminosidad
def Lcentro(r, K, T_c):
	'''
	Estudio de la luminosidad en una capa central determinada
	Input: radio, constante de polítropo y temperatura central 
	Output: valor de la luminosidad de una capa central
	'''
	global mu, generacion_energia
	eps = generacion_energia(T_c)
	return 0.006150 * eps * (K ** 2) * (T_c ** 3) * (r ** 3) * (mu ** 2)

# Ecuación (45): temperatura
def Tcentro(r, K, T_c): 
	'''
	Estudio de la temperatura en una capa central determinada
	Input: radio, constante de polítropo y temperatura central 
	Output: valor de la temperatura de una capa central
	'''
	global mu
	return  T_c - 0.008207 * (mu ** 2) * K * (T_c ** 1.5) * (r ** 2)

# Ecuación (46): presión
def Pcentro(r, K, T_centro):
	'''
	Estudio de la presión en una capa central determinada
	Input: radio, constante de polítropo y temperatura central 
	Output: valor de la presión de una capa central
	'''  
	return K * (T_centro ** 2.5)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# ERROR RELATIVO TOTAL
#------------------------------------------------------------------------------------------------------------------------
def TotalRelEror(P_Interp, P_Interp_cen, T_Interp, T_Interp_cen, L_Interp, L_Interp_cen,
 M_Interp, M_Interp_cen):
	'''
	Estudio del error total relativo en la frontera radioactiva-convectiva
	Inputs: valores de presión, temperatura, luminosidad y masa interpolados, para ambas integraciones
	Output: valor del error relativo total
	'''    
	a = ((P_Interp - P_Interp_cen) / P_Interp) ** 2
	b = ((T_Interp - T_Interp_cen) / T_Interp) ** 2
	c = ((L_Interp - L_Interp_cen) / L_Interp) ** 2
	d = ((M_Interp - M_Interp_cen) / M_Interp) ** 2
	return sqrt(a + b + c + d)
#------------------------------------------------------------------------------------------------------------------------ 
