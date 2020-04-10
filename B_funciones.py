from A_parametros_iniciales import *
from numpy import sqrt

# FUNCIONES PARA CALCULAR LA GENERACIÓN DE ENERGÍA
#
# Instrucción: 
# Calculamos eps mediante:
# generacion_energia(0.5)
# y luego clasificamos por qué medio se produce la generación de energía
# print(ciclo(0.5))

def generacion_energia(T): # T / 10**6 
    global nu_pp, nu_CN, eps_1_pp, eps_1_CN, eps_pp, eps_CN, eps, C_l, nu
    T *= 10
# Parámetros para la cadena protón-protón

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

# Parámetros para el ciclo CN

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

    eps_pp = eps_1_pp*(X**2)*(T**nu_pp)
    eps_CN = eps_1_CN*(X*Z/3)*(T**nu_CN)

    #Se elige el ciclo que produzca más energía    
    eps = max(eps_pp, eps_CN)
    nu = max(nu_pp, nu_CN)

    Cl1 = 0.01845 * eps_1_pp * (X**2) * (10 ** nu_pp) * (mu ** 2)
    Cl2 = 0.01845 * eps_1_CN * (X*Z/3) * (10 ** nu_CN) * (mu ** 2)
    C_l = max(Cl1, Cl2)
    return eps

def ciclo(T):
    generacion_energia(T)
    if eps < 0.017: # Por debajo de 0.017, considero que la generación de la energía es nula
        return "--"
    elif eps == eps_pp:
        return "PP"
    else:
        return "CN"

# FUNCIÓN PARA CALCULAR EL ERROR RELATIVO

Erel_max = 0.0001
def error_relativo(calculada, estimada):
    return abs(calculada - estimada) / calculada

#ECUACIONES FUNDAMENTALES: CASO RADIACTIVO DE LA TABLA 3
#----------------------------------------------------------------------------
#-------------------------FUNCTIONS FOR PHYICAL PARAMETERS (RADIOACTIVITE PHASES)------------------------
# We are going to define the functions to calculate the different physical parameters
# First the functions for the first derivative of the parameters

# Ecuación (18)

# Mass first derivative
def dM(r, T, P):
    '''
    Study of the mass first derivative
    Input: Pressure, Temperature and Radius
    Output: Value of fi_M
    '''
    global mu
    C_m = 0.01523 * mu
    return C_m * (P / T) * (r ** 2)

# Ecuación (19)

# Presure first derivative
def dP(r, T, P, M): 
    '''
    Study of the pressure first derivative
    Input: Pressure, Temperature, Radius and Mass
    Output: Value of fi_P
    '''
    global mu
    C_p = 8.084 * mu
    return  - C_p * (P / T) * M / (r ** 2)

# Ecuación (20)

# Temperture first derivative
def dT(r, T, P, L):
    '''
    Study of the temperature first derivative
    Input: Pressure, Temperature, Radius and Luminosity
    Output: Value of fi_T
    '''
    global mu, Z, X
    C_t = 0.01679 * Z * (1 + X) * (mu) ** 2
    return  - C_t * (P ** 2) * L * (T ** (-17 / 2)) * (r ** (-2))

# Ecuación (21)

# Luminosity first derivative
def dL(r, T, P):
    '''
    Study of the luminosity first derivative
    Input: Pressure, Temperature and Radius
    Output: Value of fi_L
    '''
    global mu, C_l, nu
    ciclo(T)
    if eps < 0.017: # Para que dL = 0 (considero que se genera energía por encima de 0.017)
        C_l = 0
    return C_l * (P ** 2) * (T ** (nu - 2)) * (r ** 2)


#ECUACIONES (35) Y (36) DE LA TABLA 4
#----------------------------------------------------------------------------
# Constantes de las expresiónes matemáticas

# Then we define the functions to calculate the T(r) and P(T)
# We will use this functions for the first three shells

# Ecuación (35)

# R_tot lo pongo en la función porque va a ir cambiando para miminizar mi error
# Temperature in a certain shell
def T(r, R_tot): # Unidades: K
    '''
    Study of the temperature in a certain shell
    Input: Radius of the shell and the total Radius of the star
    Output: Value of T(r)
    '''
    global mu, M_tot
    A_1 = 1.9022 * mu * M_tot 
    return A_1 * ((1 / r) - (1 / R_tot))

# Ecuación (36)

# Pressure in a certain shell
def P(T, L_tot): #Unidades: 10**(15) din cm-2
    '''
    Study of the pressure in a certain shell
    Input: Temperature of the shell and the total Luminosity of the star
    Output: Value of P(r)
    '''
    global M_tot, mu, Z, X
    A_2 = 10.645 * sqrt((M_tot / (mu * Z * (1 + X) * L_tot)))
    return A_2 * T ** (17/4)


# First of all we have to define the lists where we will host the parameters
# Lists that host the fi or each parameter
fi_P = [] # Pressure fi
fi_T = [] # Temperature fi
fi_M = [] # Mass fi
fi_L = [] # Luminosity fi
 
# Lists that host the value or each parameter
R_rad_i = [] # Radius in the radioactive surface
P_rad_i = [] # Pressure in the radioactive surface
T_rad_i = [] # Temperature in the radioactive surface
M_rad_i = [] # Mass in the radioactive surface
L_rad_i = [] # Luminosity in the radioactive surface
E_rad_i = [] # Energy rate production in the radioactive surface
n_rad_i = [] # n-parameter in the radioactive surface
f_rad_i = [] # phase in the radioative surface

# Estimated Pressure
def presion_estimada(P_pres, f_i, h, i): # P_pres = integracion[i][4]
    '''
    Study of the estimated pressure in a certain shell
    Input: Pressure of the shell, list of pressure first derivative and the step between shells 
    Output: Value of Estimated Pressure 
    '''
    fi = f_i[i][4]
    fi_1 = f_i[i-1][4]
    fi_2 = f_i[i-2][4]

    delta1 = h * fi - h * fi_1
    delta2 = h * fi - 2 * h * fi_1 + h * fi_2
    return P_pres + h * fi + (1/2) * delta1 + (5/12) * delta2


# Estimated Temperatere
def temperatura_estimada(T_temp, f_i, h, i): #T_temp = integracion[i][5]
    '''
    Study of the estimated temperature in a certain shell
    Input: Temperature of the shell, list of temperature first derivative and the step between shells 
    Output: Value of Estimated Temperature
    '''
    fi = f_i[i][5]
    fi_1 = f_i[i-1][5]

    delta1 = h * fi - h * fi_1
    return T_temp + h * fi + (1/2) * delta1

# Calculated Pressure, Temperature and Mass
def PTM_calculada(ptm, dptm, fi, h, i): # P = integracion[i][4], T = integracion[i][5], M = integracion[i][7]
# P: fi = f_i[i][4], T: fi = f_i[i][5], M: fi = f_i[i][7]    
    '''
    Study of the calculated pressure, temperature or mass in a certain shell
    Input: Parameter value of the shell, list of parameter values first derivative and the step between shells 
    Output: Value of Calculated Pressure, Temperature or Mass
    '''
    delta1 = h * dptm - h * fi
    return ptm + h * dptm - (1/2) * delta1

# Calculated Luminosity
def L_calculada(L, dLu, f_i, h, i): # L = L_tot
    '''
    Study of the calculated luminosity in a certain shell
    Input: Luminosity of the shell, list of luminosity first derivative and the step between shells 
    Output: Value of Calculated Luminosity
    '''
#    fi1 = dLu
#    fi = f_i[i][6]
#   fi_1 = f_i[i-1][6]

#    return L + (h/12) * (5*fi1 + 8*fi - fi_1)


    fi = f_i[i][6]
    fi_1 = f_i[i-1][6]
    fi_2 = f_i[i-2][6]

    return L + (h/12) * (5*fi + 8*fi_1 - fi_2)

# Calculated Luminosity
def L_calculadaa(L, dLu, f_i, h, i): # L = L_tot
    '''
    Study of the calculated luminosity in a certain shell
    Input: Luminosity of the shell, list of luminosity first derivative and the step between shells 
    Output: Value of Calculated Luminosity
    '''
#    fi1 = dLu
#    fi = f_i[i][6]
#   fi_1 = f_i[i-1][6]

#    return L + (h/12) * (5*fi1 + 8*fi - fi_1)


    fi = f_i[i][6]
    fi_1 = f_i[i-1][6]

    delta1 = h * dLu - h * fi
    delta2 = h * dLu - 2 *  h * fi + h * fi_1
    return L + h * dLu - (1/2) * delta1 - (1/12) * delta2



def L_cal_fun(Li, fi1_L, fi_L, h):
    '''
    Study of the calculated luminosity in a certain shell
    Input: Luminosity of the shell, list of luminosity first derivative and the step between shells 
    Output: Value of Calculated Luminosity
    '''
    fi1 = fi1_L
    fi = fi_L[-1]
    fi_1 = fi_L[-2]
    return Li + (h/12) * (5*fi1 + 8*fi - fi_1)

# We are going to define n+1:
def n_fun(Tcal, Pcal, dPr, dTe):
    '''
    Study the value of the parameter that gives us information about the radioactive or convective behavior
    Input: Calculated Temperature and Pressure of the shell, pressure and temperature first derivative of the shell
    Output: Value of the n-parameter 
    '''
    return (Tcal / Pcal) * (dPr / dTe)


def interpolation(y0,y1,n0,n1):
    '''
    Generates a linear interpolation
    Input: two pairs of values (x,y) to make the linear intepolation
    Output: The value interpolated for n=2.5
    '''
    n = 2.5
    return y0 + ((y1 - y0) / (n1 - n0)) * (n - n0)

#---------------------------FUNCTIONS FOR PHYICAL PARAMETERS (CONVECTIVE PHASE)--------------------------
# We are going to define the equations to obtain the parameters in the convective core
# First the functions for the first derivative of the parameters

#ECUACIONES FUNDAMENTALES: CASO CONVECTIVO DE LA TABLA 3
#----------------------------------------------------------------------------
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
    global mu, C_l, nu
    ciclo(T)   
    return C_l * (K ** 2) * (T ** (3 + nu)) * (r ** 2)

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

#ECUACIONES FUNDAMENTALES: CASO RADIACTIVO DE LA TABLA 5
#----------------------------------------------------------------------------

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
    global mu, eps
    ciclo(T_c)
    return 0.006150 * eps * (K ** 2) * (T_c ** (3)) * (r ** 3) * (mu ** 2)


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
def Pcentro(r, T_centro, K):
    '''
    Study of the Mass in a certain shell
    Input: Radius, Temperature and the Polytrope model constant 
    Output: Value of P(r)
    '''    
    return K * (T_centro ** 2.5)


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
 
