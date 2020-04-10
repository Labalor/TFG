#####################################################################################
#####################################################################################

#                           TFG dirigido por NICOLÁS CARDIEL 
# 
#                   ELAVORACIÓN DE UN MODELO NUMÉRICO DE INTERIOR ESTELAR
# 
#                              Autor: LUIS ABALO RODRÍGUEZ
# 
#                              FACULTAD DE CIENCIAS FÍSICAS
#                           UNIVERSIDAD COMPLUTENSE DE MADRID
# 
#                               CURSO ACADÉMICO 2019-2020
#
#####################################################################################
#####################################################################################

#                                     COMENTARIOS: 
# 
#                                         (1) 
# Las ecuaciónes que se enumeran en el presente código siguen las equiquetas del guión 

# ___________________________________________________________________________________
# ___________________________________________________________________________________
# 
#                                DESARROLLO DEL CÓDIGO 
# ___________________________________________________________________________________
# ___________________________________________________________________________________

#VALORES INICIALES
#----------------------------------------------------------------------------
M_tot = 5 #Masa total de la Tierra ; Unidades: e33 g
L_tot = 70 
R_tot = 11.5
T_c = 2.0
Erel_max = 0.0001
#----------------------------------------------------------------------------
#
#
#
#COMPOSICIÓN QUÍMICA
#----------------------------------------------------------------------------
X = 0.75 #Fracción en masa del Hidrógeno
Y = 0.22 #Fracción en masa del Helio
Z = 1 - X - Y #Fracción en masa de elementos pesados
mu = 1 / (2*X + .75*Y + .5*Z) #Peso molecular medio
#----------------------------------------------------------------------------
#
#
#
#CARGA DE LIBRERÍAS
#----------------------------------------------------------------------------
import numpy as np 
from tabulate import tabulate
#----------------------------------------------------------------------------
#
#
#
#RITMO DE GENERACIÓN DE ENERGÍA
#----------------------------------------------------------------------------
# En función de la temperatura, se calculan los parámetros del ritmo de generación de energía 
# para los dos casos
print("_" * 78)
print("GENERACIÓN DE ENERGÍA")

def gener_energia(temp): # T / 10**6 
    eps_pp,eps_CN = 1,1
    global nu_pp
    global nu_CN

# Parámetros para la cadena protón-protón

    if temp >= 4.0 and temp <= 6.0:
        eps_1_pp = 10**(-6.84) 
        nu_pp = 6.0
    elif temp >= 6.0 and temp <= 9.5:
        eps_1_pp = 10**(-6.04) 
        nu_pp = 5.0
    elif temp >= 9.5 and temp <= 12.0:
        eps_1_pp = 10**(-5.56) 
        nu_pp = 4.5
    elif temp >= 12.0 and temp <= 16.5:
        eps_1_pp = 10**(-5.02) 
        nu_pp = 4.0
    elif temp >= 16.5 and temp <= 24.0:
        eps_1_pp = 10**(-4.40) 
        nu_pp = 3.5
    else:
        eps_pp = 0
        nu_pp = 0

# Parámetros para el ciclo CN

    if temp >= 12.0 and temp <= 16.0:
        eps_1_CN = 10**(-22.2) 
        nu_CN = 20.0
    elif temp >= 16.0 and temp <= 22.5:
        eps_1_CN = 10**(-19.8) 
        nu_CN = 18.0
    elif temp >= 22.5 and temp <= 27.5:
        eps_1_CN = 10**(-17.1) 
        nu_CN = 16.0
    elif temp >= 27.5 and temp <= 36.0:
        eps_1_CN = 10**(-15.6) 
        nu_CN = 15.0
    elif temp >= 36.0 and temp <= 50.0:
        eps_1_CN = 10**(-12.5) 
        nu_CN = 13.0
    else:
        eps_CN = 0
        nu_CN = 0    

    if eps_pp != 0:
        eps_pp = eps_1_pp*(X**2)*(temp**nu_pp)
        global Cl1
        Cl1 = 0.01845 * eps_1_pp * (X**2) * (10 ** nu_pp) * (mu ** 2)
        print("Cl_PP",Cl1)

    if eps_CN != 0:
        eps_CN = eps_1_CN*(X*Z/3)*(temp**nu_CN)
        global Cl2
        Cl2 = 0.01845 * eps_1_CN * (X*Z/3) * (10 ** nu_CN) * (mu ** 2)
        print("Cl_CN",Cl2)

    #Se elige el ciclo que produzca más energía    
    global eps
    if eps_pp > eps_CN:
        eps = eps_pp
        print("Para una temperatura de",temp,"K se elige la cadena protón-protón.")
        print('RITMO DE GENERACIÓN DE ENERGÍA: ',eps)
        global k
        k = 0
    else:
        eps = eps_CN
        print("Para una temperatura de ",temp,"K se elige el ciclo CN.")
        print('RITMO DE GENERACIÓN DE ENERGÍA: ',eps)
#
# ++++ No copiado este caso para T = 17
gener_energia(17)
print("Cl1",Cl1, "Cl2", Cl2)
k = 1
if k == 0:
    Cl = Cl1
    nu = nu_pp
else:
    Cl = Cl2
    nu = nu_CN
print("Cl definitivo", Cl, "nu definitivo", nu)
#
#----------------------------------------------------------------------------
#
#
#
#RADIO INICIAL Y RADIO DE CADA CAPA
#----------------------------------------------------------------------------
R_ini = 0.9 * R_tot     #Para evitar problemas de convergencia 
h = - R_ini / 100 #Paso de integración
r = np.linspace(R_ini,0,101) #Radio 
#----------------------------------------------------------------------------
#
#
#
#ECUACIONES (35) Y (36) DE LA TABLA 4
#----------------------------------------------------------------------------
# Constantes de las expresiónes matemáticas

A_1 = 1.9022 * mu * M_tot 
A_2 = 10.645 * np.sqrt((M_tot/(mu*Z*(1 + X)*L_tot)))

# Ecuación (35)
def T(j): # Unidades: K
    return A_1 * ((1/j) - (1/R_tot))

# Ecuación (36)

def P(T_temp): #Unidades: 10**(15) din cm-2
    return A_2 * T_temp ** (17/4)
#----------------------------------------------------------------------------
#
#
#
#CONSTANTES DE ECUACIONES DE LA TABLA 3
#----------------------------------------------------------------------------
C_m = 0.01523*mu
C_p = 8.084*mu
C_t = 0.01679*Z*(1+X)*(mu)**2
C_t_c = 3.234 * mu
#----------------------------------------------------------------------------
#
#
#
#ECUACIONES FUNDAMENTALES: CASO RADIACTIVO DE LA TABLA 3
#----------------------------------------------------------------------------
# Ecuación (18)

def dM(r, T_temp, P_pres):
    return C_m*(P_pres/T_temp)*(r**2)

# Ecuación (19)
def dP(r, T_temp, P_pres, M): 
    return  -C_p*(P_pres/T_temp)*M/(r**2)

# Ecuación (20)
def dT(r, T_temp, P_pres, L):
    return  -C_t * (P_pres ** 2) * L * (T_temp ** (-17/2)) * (r ** (-2))

# Ecuación (21)
def dL(r, T_temp, P_pres, nu, C_l):
    return C_l * (P_pres ** 2) * (T_temp ** (nu-2)) * (r **2)
#----------------------------------------------------------------------------
#
#
#
#ECUACIONES FUNDAMENTALES: CASO CONVECTIVO DE LA TABLA 3
#----------------------------------------------------------------------------
# Ecuación (22)
def dMc(T_temp,r):
    return C_m * K * (T_temp ** 1.5) * (r ** 2)

# Ecuación (23)
def dPc(T_temp,r, M):
    return C_p * K * (T_temp ** 1.5) * M * (r ** (-2))

# Ecuación (24)
def dLc(r, T_temp, nu, C_l):
    return C_l * (K ** 2) * (T_temp ** (3 + nu)) * (r ** 2)

# Ecuación (25)
def dTc(r, M):
    return - C_t_c * M / (r ** 2)
#----------------------------------------------------------------------------
#
# ++++ continua aquí
#
#PRIMERAS TRES CAPAS
#----------------------------------------------------------------------------
integracion = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]
f_i = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]

for i in np.arange(3): # Crea [0, 1, 2]
    T_temp = T(r[i])
    P_pres = P(T_temp)
    f_i_dM = dM(r[i], T_temp, P_pres)
    f_i_dP = dP(r[i], T_temp, P_pres, M_tot)
    f_i_dT = dT(r[i], T_temp, P_pres, L_tot)
    f_i_dL = 0

    integracion.append(['--', 'INICIO',  i, r[i], P_pres, T_temp, L_tot, M_tot])
    f_i.append(['--', 'INICIO', i, r[i], f_i_dP, f_i_dT, f_i_dL, f_i_dM])

# ++++ no copio las siguientes dos líneas que printean
print(tabulate(integracion, headers='firstrow', tablefmt='fancy_grid'))
print(tabulate(f_i, headers='firstrow', tablefmt='fancy_grid'))
# print(f_i[3][4])
#----------------------------------------------------------------------------
#
#
#
#ALGORITMO FASE A.1.1. Envoltura radiativa: luminosidad y masa constantes
#----------------------------------------------------------------------------
i = 3 #número de capa
integracionA11 = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]

print("_" * 78)
print("Algoritmo A.1.1.")

loop1 = True
while loop1:
    dLu = 0 #No hay generación de energía
    print("\n * Calculando capa número",i,"\n")
    print("Ejecutando paso 1")
    print("Ejecutando paso 2")
#   Cuando selecciono la fila i, como tiene firstrow = headers, me selecciona la anterior
#   La columna 4 es la presión (en f_i es la calculada por dP)
    delta_1_iP = h * f_i[i][4] - h * f_i[i-1][4]
    delta_2_iP = h * f_i[i][4] - 2 * h * f_i[i-1][4] + h * f_i[i-2][4]
    Pest = integracion[i][4] + h * f_i[i][4] + (1/2) * delta_1_iP + (5/12) * delta_2_iP

    delta_1_iT = h * f_i[i][5] - h * f_i[i-1][5]
#    delta_2_i = h * f_i[i][5] - 2 * h * f_i[i-1][5] + h * f_i[i-2][5]
    Test = integracion[i][5] + h * f_i[i][5] + (1/2) * delta_1_iT
 
    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
            print("Ejecutando paso 4")
            dPr = dP(r[i], Test, Pest, M_tot)
            delta_1_i1P = h * dPr - h * f_i[i][4]
            Pcal = integracion[i][4] + h * dPr - (1/2) * delta_1_i1P

            ErelP = abs(Pcal - Pest) / Pcal

            if ErelP < Erel_max:
                loop3 = False
            else:
                Pest = Pcal

        print("Ejecutando paso 7")
        dTe = dT(r[i], Test, Pcal, L_tot)
        delta_1_i1T = h * dTe - h * f_i[i][5]
        Tcal = integracion[i][5] + h * dTe - (1/2) * delta_1_i1T

        ErelT = abs(Tcal - Test) / Tcal

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal


    print("Ejecutando paso 3")
    dMa = dM(r[i], Tcal, Pcal)
    delta_1_i1M = h * dMa - h * f_i[i][7]
    Mcal = M_tot + h * dMa - (1/2) * delta_1_i1M 

    ErelM = abs(Mcal - M_tot) / M_tot

    if ErelM > Erel_max :
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.1.2. -- \n")
        print(tabulate(integracionA11, headers='firstrow', tablefmt='fancy_grid'))         
    else:
        #pasamos a la siguiente capa
        Mcal = M_tot
        integracionA11.append(['--', 'A.1.1.',  i, r[i], Pcal, Tcal, L_tot, Mcal])
        integracion.append(['--', 'A.1.1.',  i, r[i], Pcal, Tcal, L_tot, Mcal])
        f_i.append(['--', 'A.1.1.',  i, r[i], dPr, dTe, dLu, dMa])

        i += 1

print("_" * 78)
print("Algoritmo A.1.2.")
print("_" * 78)
#----------------------------------------------------------------------------
#
#
#
#ALGORITMO FASE A.1.2. Envoltura radiativa: luminosidad constantes y masa variable
#----------------------------------------------------------------------------
integracionA12 = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]

loop1 = True
while loop1:
    print("\n * Calculando capa número",i,"\n")
    print("Ejecutando paso 1")
    print("Ejecutando paso 2")
    delta_1_iP = h * f_i[i][4] - h * f_i[i-1][4]
    delta_2_iP = h * f_i[i][4] - 2 * h * f_i[i-1][4] + h * f_i[i-2][4]
    Pest = integracion[i][4] + h * f_i[i][4] + (1/2) * delta_1_iP + (5/12) * delta_2_iP

    delta_1_iT = h * f_i[i][5] - h * f_i[i-1][5]
#    delta_2_i = h * f_i[i][5] - 2 * h * f_i[i-1][5] + h * f_i[i-2][5]
    Test = integracion[i][5] + h * f_i[i][5] + (1/2) * delta_1_iT

    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
            print("Ejecutando paso 3")
            dMa = dM(r[i], Test, Pest)
            delta_1_i1M = h * dMa - h * f_i[i][7]
            Mcal = integracion[i][7] + h * dMa - (1/2) * delta_1_i1M  
            print("Ejecutando paso 4")
            dPr = dP(r[i], Test, Pest, Mcal)
            delta_1_i1P = h * dPr - h * f_i[i][4]
            Pcal = integracion[i][4] + h * dPr - (1/2) * delta_1_i1P

            ErelP = abs(Pcal - Pest) / Pcal

            if ErelP < Erel_max:
                loop3 = False
            else:
                Pest = Pcal

        print("Ejecutando paso 7")
        dTe = dT(r[i], Test, Pcal, L_tot)
        delta_1_i1T = h * dTe - h * f_i[i][5]
        Tcal = integracion[i][5] + h * dTe - (1/2) * delta_1_i1T

        ErelT = abs(Tcal - Test) / Tcal

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal

    print("Ejecutando paso 6")
    Tgen = Tcal * 10 #Para que esté en K/e6, unidades de la funcion de gener_energia (aquí K/e7)
    k = 1
    # gener_energia(Tgen) #Me genera un Tgen = 0.8, por eso no me coge ningun Cl2. No puedo borrar
    # de 136 a 145 porque necesita un valor para Cl2
    # if k == 0:
    #     Cl = Cl1
    #     nu = nu_pp
    # else:
    #     Cl = Cl2
    #     nu = nu_CN
    # print("Cl definitivo", Cl, "nu definitivo", nu)
    if eps < 0.017:
        dLu = 0
        Lcal = L_tot
    else:
        dLu = dL(r[i], Tgen, Pcal, nu, Cl)
        
        delta_1_i1L = h * dLu - h * f_i[i][6]
        delta_2_i1L = h * dLu - 2 *  h * f_i[i][6] + h * f_i[i-1][6]

        Lcal = L_tot + h * dLu - (1/2) * delta_1_i1L - (1/12) * delta_2_i1L

    ErelL = abs(Lcal - L_tot) / L_tot


    if ErelL > Erel_max :
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.1.3. -- \n")
        print(tabulate(integracionA12, headers='firstrow', tablefmt='fancy_grid'))         
    else:
        #pasamos a la siguiente capa
        integracionA12.append(['--', 'A.1.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        integracion.append(['--', 'A.1.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i.append(['--', 'A.1.2.',  i, r[i], dPr, dTe, dLu, dMa])

        i += 1

print("_" * 78)
print("Algoritmo A.1.3.")
print("_" * 78)
#----------------------------------------------------------------------------
#
#
#
#ALGORITMO FASE A.1.3. Envoltura radiativa: masa y luminosidad variable
#----------------------------------------------------------------------------
integracionA13 = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M', 'n+1']]

loop1 = True
while loop1:
    print("\n * Calculando capa número",i,"\n")
    print("Ejecutando paso 1")
    print("Ejecutando paso 2")
    delta_1_iP = h * f_i[i][4] - h * f_i[i-1][4]
    delta_2_iP = h * f_i[i][4] - 2 * h * f_i[i-1][4] + h * f_i[i-2][4]
    Pest = integracion[i][4] + h * f_i[i][4] + (1/2) * delta_1_iP + (5/12) * delta_2_iP

    delta_1_iT = h * f_i[i][5] - h * f_i[i-1][5]
#    delta_2_i = h * f_i[i][5] - 2 * h * f_i[i-1][5] + h * f_i[i-2][5]
    Test = integracion[i][5] + h * f_i[i][5] + (1/2) * delta_1_iT



    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
            print("Ejecutando paso 3")
            dMa = dM(r[i], Test, Pest)
            delta_1_i1M = h * dMa - h * f_i[i][7]
            Mcal = integracion[i][7] + h * dMa - (1/2) * delta_1_i1M  
            print("Ejecutando paso 4")
            dPr = dP(r[i], Test, Pest, Mcal)
            delta_1_i1P = h * dPr - h * f_i[i][4]
            Pcal = integracion[i][4] + h * dPr - (1/2) * delta_1_i1P

            ErelP = abs(Pcal - Pest) / Pcal

            if ErelP < Erel_max:
                loop3 = False
            else:
                Pest = Pcal
        print("Ejecutando paso 6")
        Tgen = Test * 10 #Para que esté en K/e6, unidades de la funcion de gener_energia (aquí K/e7)
        k = 1
        gener_energia(Tgen)
        if k == 0:
            Cl = Cl1
            nu = nu_pp
            E = 'PP'
        else:
            Cl = Cl2
            nu = nu_CN
            E = 'CN'
        print("Cl definitivo", Cl, "nu definitivo", nu)

        dLu = dL(r[i], Test, Pcal, nu, Cl)
        delta_1_i1L = h * dLu - h * f_i[i][6]
        delta_2_i1L = h * dLu - 2 *  h * f_i[i][6] + h * f_i[i-1][6]

        Lcal = integracion[i][6] + h * dLu - (1/2) * delta_1_i1L - (1/12) * delta_2_i1L

        print("Ejecutando paso 7")
        dTe = dT(r[i], Test, Pcal, Lcal)
        delta_1_i1T = h * dTe - h * f_i[i][5]
        Tcal = integracion[i][5] + h * dTe - (1/2) * delta_1_i1T

        ErelT = abs(Tcal - Test) / Tcal

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal

    print("Ejecutando paso 9")
    n = (Tcal / Pcal) * (dPr / dTe) 
    print(n)

    if n < 2.5:
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.2. -- \n")
        print(tabulate(integracionA13, headers='firstrow', tablefmt='fancy_grid'))
        global K
        K = Pcal / (Tcal ** 2.5) # Constante del polítropo         
    else:
        #pasamos a la siguiente capa
        integracionA13.append([E, 'A.1.3.',  i, r[i], Pcal, Tcal, Lcal, Mcal, n])
        integracion.append([E, 'A.1.3.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i.append([E, 'A.1.3.',  i, r[i], dPr, dTe, dLu, dMa])

        i += 1

print("_" * 78)
print("Algoritmo A.2.")
print("_" * 78)
print("Núcleo convectivo")
#----------------------------------------------------------------------------
#
#
#
#ALGORITMO FASE A.2. Núcleo convectivo
#----------------------------------------------------------------------------
integracionA2 = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M', 'n+1']]
jj = 1

loop1 = True
while loop1:
    print("\n * Calculando capa número",i,"\n")
    print("Ejecutando paso 1")
    print("Ejecutando paso 2bis")
    delta_1_iT = h * f_i[i][5] - h * f_i[i-1][5]
    Test = integracion[i][5] + h * f_i[i][5] + (1/2) * delta_1_iT

    loop2 = True
    while loop2:
        print("Ejecutando Polítropo")
        Pest = K * (Test ** 2.5)

        print("Ejecutando paso 3")
        print("Test", Test, "Pest", Pest)
        dMa = dM(r[i], Test, Pest)
        print("dMa", dMa)
        delta_1_i1M = h * dMa - h * f_i[i][7]
        Mcal = integracion[i][7] + h * dMa - (1/2) * delta_1_i1M  
        print("Mcal", Mcal)

        print("Ejecutando paso 7bis")
        if r[i] == 0:
            Tcal = Test
        else:
            dTe = dTc(r[i], Mcal)
            delta_1_i1T = h * dTe - h * f_i[i][5]
            Tcal = integracion[i][5] + h * dTe - (1/2) * delta_1_i1T
        print("Tcal", Tcal)

        ErelT = abs(Tcal - Test) / Tcal
        print("ErelT",ErelT)
        jj += 1
        if jj == 100:
            loop2 = False

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal
    Pcal = K * (Tcal ** 2.5)
    

    print("Ejecutando paso 6")
    Tgen = Tcal * 10 #Para que esté en K/e6, unidades de la funcion de gener_energia (aquí K/e7)
    k = 1
    gener_energia(Tgen)
    if k == 0:
        Cl = Cl1
        nu = nu_pp
        E = 'PP'
    else:
        Cl = Cl2
        nu = nu_CN
        E = 'CN'
    print("Cl definitivo", Cl, "nu definitivo", nu)

    dLu = dL(r[i], Tcal, Pcal, nu, Cl)
    delta_1_i1L = h * dLu - h * f_i[i][6]
    delta_2_i1L = h * dLu - 2 *  h * f_i[i][6] + h * f_i[i-1][6]

    Lcal = integracion[i][6] + h * dLu - (1/2) * delta_1_i1L - (1/12) * delta_2_i1L

    K = Pcal / (Tcal ** 2.5) # Constante del polítropo         
    print("K válida", K)

    if r[i] == 0:
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.2. -- \n")
        integracionA2.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal, n])
        integracion.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        n = 0
        i += 1
        print(tabulate(integracionA2, headers='firstrow', tablefmt='fancy_grid'))         
    else:
        #pasamos a la siguiente capa
        integracionA2.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal, n])
        integracion.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        n = 0
        i += 1

# print("_" * 78)
# print("Fase B --integracion desde el centro--")
# print("_" * 78)
#----------------------------------------------------------------------------
#
print(tabulate(integracion, headers='firstrow', tablefmt='fancy_grid'))      
#
# VALORES EN LA FRONTERA (interpolando linealmente para cuando n + 1 = 2.5)
#----------------------------------------------------------------------------
# Valores en la zona radiactiva
radio_Antes25 = r[81]
P_Antes25 = integracionA13[22][4]
T_Antes25 = integracionA13[22][5]
L_Antes25 = integracionA13[22][6] 
M_Antes25 = integracionA13[22][7]
n_Antes25 = integracionA13[22][8]

# Valores en la zona convectiva
radio_Despues25 = r[82]
P_Despues25 = integracionA2[1][4]
T_Despues25 = integracionA2[1][5]
L_Despues25 = integracionA2[1][6]
M_Despues25 = integracionA2[1][7]
n_Despues25 = integracionA2[1][8]

from scipy.optimize import curve_fit

def poldeg1(x,a,b):
    return a * x + b

equis = [n_Antes25,n_Despues25]

val, cov = curve_fit(poldeg1,equis,[radio_Antes25,radio_Despues25])
a,b = val
r_Interp = poldeg1(2.5,a,b)
#print("r_Interp",r_Interp)

val, cov = curve_fit(poldeg1,equis,[P_Antes25,P_Despues25])
a,b = val
P_Interp = poldeg1(2.5,a,b)
# print("P_Interp",P_Interp)

val, cov = curve_fit(poldeg1,equis,[T_Antes25,T_Despues25])
a,b = val
T_Interp = poldeg1(2.5,a,b)
# print("T_Interp",T_Interp)

val, cov = curve_fit(poldeg1,equis,[L_Antes25,L_Despues25])
a,b = val
L_Interp = poldeg1(2.5,a,b)
# print("L_Interp",L_Interp)

val, cov = curve_fit(poldeg1,equis,[M_Antes25,M_Despues25])
a,b = val
M_Interp = poldeg1(2.5,a,b)
# print("M_Interp",M_Interp)



interpolacion = [['r', 'P', 'T', 'L', 'M']]
interpolacion.append([r_Interp, P_Interp, T_Interp, L_Interp, M_Interp])

print(tabulate(interpolacion, headers='firstrow', tablefmt='fancy_grid')) 

#----------------------------------------------------------------------------
# Método dos de interpolación
#
# def interpolation(y0,y1):
#     
#     n = 2.5
#     n0 = 2.500934
#     n1 = 2.379468 
#     return y0 + ((y1 - y0) / (n1 - n0)) * (n - n0)
# 
# 
# We have obtained the values manually
# R_rad = interpolation(1.96650, 1.86300)
# P_rad = interpolation(39.8420132, 42.3165500)
# T_rad = interpolation(1.4562414, 1.4917669)
# L_rad = interpolation(66.582873, 66.097434)
# M_rad = interpolation(0.702943, 0.607118)
#----------------------------------------------------------------------------
#
#
#
# INTEGRANDO DESDE EL CENTRO
#----------------------------------------------------------------------------
#
#
#
#RADIO INICIAL Y RADIO DE CADA CAPA
#----------------------------------------------------------------------------
R_ini = 0.9 * R_tot #Para evitar problemas de convergencia 
h = R_ini / 100 #Paso de integración
r = np.linspace(0,R_ini,101) #Radio 
#----------------------------------------------------------------------------
#
#
# +++++ Copiadas en "funciones.py"
#ECUACIONES FUNDAMENTALES: CASO RADIACTIVO DE LA TABLA 5
#----------------------------------------------------------------------------
def Mcentro(r):
    return 0.005077 * mu * K * (T_c ** 1.5) * (r ** 3)

def Lcentro(r):
    return 0.006150 * eps * (K ** 2) * (T_c ** (3)) * (r ** 3) * (mu ** 2)

def Tcentro(r): 
    return  T_c - 0.008207 * (mu ** 2) * K * (T_c ** 1.5) * (r ** 2)
#----------------------------------------------------------------------------
#
#
#
#PRIMERAS TRES CAPAS
#----------------------------------------------------------------------------
centro = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]
f_i_centro = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]
# f_i = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]

gener_energia(T_c * 10) 

print("mu", mu, "K", K, "T_c", T_c, "eps", eps, "nu", nu)
k = 1
if k == 0:
    Cl = Cl1
    nu = nu_pp
else:
    Cl = Cl2
    nu = nu_CN
# K = 15.547208668705919

Temp_centro = T_c
for i in np.arange(3):
    M_centro = Mcentro(r[i])
    L_centro = Lcentro(r[i])
    T_centro = Tcentro(r[i])
    P_centro = K * (T_centro ** 2.5)

    f_i_dM = dMc(T_centro,r[i])
    f_i_dL = dLc(r[i], T_centro, nu, Cl)
    if r[i] == 0:
        f_i_dT = 0
    else:
        f_i_dT = dTc(r[i], M_centro)    

    centro.append(['--', 'INICIO',  i, r[i], P_centro, T_centro, L_centro, M_centro])
    f_i_centro.append(['--', 'INICIO', i, r[i], " - ", f_i_dT, f_i_dL, f_i_dM])

print(tabulate(centro, headers='firstrow', tablefmt='fancy_grid'))
print(tabulate(f_i_centro, headers='firstrow', tablefmt='fancy_grid'))
#
#
# Para poder continuar: 
# f_i_centro[1][5] = 0 
#
#
#----------------------------------------------------------------------------
#
#
# 
#----------------------------------------------------------------------------
#ALGORITMO FASE A.2. Núcleo convectivo
#----------------------------------------------------------------------------
centro_integracionA2 = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]
jj = 1
i = 3
print("K buena", K)


loop1 = True
while loop1:
    print("\n * Calculando capa número",i,"\n")
    print("Ejecutando paso 1")
    print("Ejecutando paso 2bis")
    delta_1_iT = h * f_i_centro[i][5] - h * f_i_centro[i-1][5]
    Test = centro[i][5] + h * f_i_centro[i][5] + (1/2) * delta_1_iT

    loop2 = True
    while loop2:
        print("Ejecutando Polítropo")
        Pest = K * (Test ** 2.5)

        print("Ejecutando paso 3")
        print("Test", Test, "Pest", Pest)
        dMa = dMc(Test, r[i])
        print("dMa", dMa)
        delta_1_i1M = h * dMa - h * f_i_centro[i][7]
        Mcal = centro[i][7] + h * dMa - (1/2) * delta_1_i1M  
        print("Mcal", Mcal)

        print("Ejecutando paso 7bis")
        dTe = dTc(r[i], Mcal)
        delta_1_i1T = h * dTe - h * f_i_centro[i][5]
        Tcal = centro[i][5] + h * dTe - (1/2) * delta_1_i1T
        print("Tcal", Tcal)

        ErelT = abs(Tcal - Test) / Tcal
        print("ErelT",ErelT)
        jj += 1
        if jj == 100:
            loop2 = False

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal
#         loop2 = False #AÑADIDO AHORA
    Pcal = K * (Tcal ** 2.5)
    

    print("Ejecutando paso 6")
    Tgen = Tcal * 10 #Para que esté en K/e6, unidades de la funcion de gener_energia (aquí K/e7)
    k = 1
    gener_energia(Tgen)
    if k == 0:
        Cl = Cl1
        nu = nu_pp
        E = 'PP'
    else:
        Cl = Cl2
        nu = nu_CN
        E = 'CN'
    print("Cl definitivo", Cl, "nu definitivo", nu)

    dLu = dLc(r[i], Tcal, nu, Cl)
    delta_1_i1L = h * dLu - h * f_i_centro[i][6]
    delta_2_i1L = h * dLu - 2 *  h * f_i_centro[i][6] + h * f_i_centro[i-1][6]

    Lcal = centro[i][6] + h * dLu - (1/2) * delta_1_i1L - (1/12) * delta_2_i1L

    K = Pcal / (Tcal ** 2.5) # Constante del polítropo         
    print("K válida", K)

    if r[i] > r_Interp:
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.2. -- \n")
        centro_integracionA2.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        centro.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i_centro.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        n = 0
        i += 1
        print(tabulate(centro_integracionA2, headers='firstrow', tablefmt='fancy_grid'))         
    else:
        #pasamos a la siguiente capa
        centro_integracionA2.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        centro.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i_centro.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        n = 0
        i += 1
# loop1 = False # AÑADIDO AHORA 
#
#
#
# VALORES EN LA FRONTERA (interpolando linealmente para cuando n + 1 = 2.5)
#----------------------------------------------------------------------------
# Valores en la zona radiactiva
radio_Antes25_cen = r[18]
P_Antes25_cen = centro_integracionA2[16][4]
T_Antes25_cen = centro_integracionA2[16][5]
L_Antes25_cen = centro_integracionA2[16][6] 
M_Antes25_cen = centro_integracionA2[16][7]
n_Antes25_cen = integracionA13[22][8]

# Valores en la zona convectiva
radio_Despues25_cen = r[19]
P_Despues25_cen = centro_integracionA2[17][4]
T_Despues25_cen = centro_integracionA2[17][5]
L_Despues25_cen = centro_integracionA2[17][6]
M_Despues25_cen = centro_integracionA2[17][7]
n_Despues25_cen = integracionA2[1][8]

from scipy.optimize import curve_fit

def poldeg1(x,a,b):
    return a * x + b

equis_cen = [n_Antes25_cen,n_Despues25_cen]

val, cov = curve_fit(poldeg1,equis_cen,[radio_Antes25_cen,radio_Despues25_cen])
a,b = val
r_Interp_cen = poldeg1(2.5,a,b)
#print("r_Interp",r_Interp)

val, cov = curve_fit(poldeg1,equis_cen,[P_Antes25_cen,P_Despues25_cen])
a,b = val
P_Interp_cen = poldeg1(2.5,a,b)
# print("P_Interp",P_Interp)

val, cov = curve_fit(poldeg1,equis_cen,[T_Antes25_cen,T_Despues25_cen])
a,b = val
T_Interp_cen = poldeg1(2.5,a,b)
# print("T_Interp",T_Interp)

val, cov = curve_fit(poldeg1,equis_cen,[L_Antes25_cen,L_Despues25_cen])
a,b = val
L_Interp_cen = poldeg1(2.5,a,b)
# print("L_Interp",L_Interp)

val, cov = curve_fit(poldeg1,equis_cen,[M_Antes25_cen,M_Despues25_cen])
a,b = val
M_Interp_cen = poldeg1(2.5,a,b)
# print("M_Interp",M_Interp)



interpolacion_cen = [['r', 'P', 'T', 'L', 'M']]
interpolacion_cen.append([r_Interp_cen, P_Interp_cen, T_Interp_cen, L_Interp_cen, M_Interp_cen])

print(tabulate(interpolacion_cen, headers='firstrow', tablefmt='fancy_grid')) 

#----------------------------------------------------------------------------
#
#
#
#----------------------------------------------------------------------------
# Error relativo total

def TotalRelEror(P_Interp, P_Interp_cen, T_Interp, T_Interp_cen, L_Interp, L_Interp_cen,
 M_Interp, M_Interp_cen):
    
    a = ((P_Interp - P_Interp_cen) / P_Interp) ** 2
    b = ((T_Interp - T_Interp_cen) / T_Interp) ** 2
    c = ((L_Interp - L_Interp_cen) / L_Interp) ** 2
    d = ((M_Interp - M_Interp_cen) / M_Interp) ** 2
    return np.sqrt(a + b + c + d)

TotRE = TotalRelEror(P_Interp, P_Interp_cen, T_Interp, T_Interp_cen, L_Interp, L_Interp_cen,
 M_Interp, M_Interp_cen)
print("Error relativo total = {:.4f}".format(100*TotRE))

#----------------------------------------------------------------------------
#
#
#
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
# Ajuste de las soluciones a un radio intermedio

error_relativo = []
Temp_c_plot = []

#PRIMERAS TRES CAPAS
#----------------------------------------------------------------------------

for i in np.linspace(1.9,2.0,1000):
    T_c = i
    centro = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]
    f_i_centro = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]
# f_i = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]

    gener_energia(T_c * 10) 

    print("mu", mu, "K", K, "T_c", T_c, "eps", eps, "nu", nu)
    k = 1
    if k == 0:
        Cl = Cl1
        nu = nu_pp
    else:
        Cl = Cl2
        nu = nu_CN
    # K = 15.547208668705919

    Temp_centro = T_c
    for i in np.arange(3):
        M_centro = Mcentro(r[i])
        L_centro = Lcentro(r[i])
        T_centro = Tcentro(r[i])
        P_centro = K * (T_centro ** 2.5)

        f_i_dM = dMc(T_centro,r[i])
        f_i_dL = dLc(r[i], T_centro, nu, Cl)
        if r[i] == 0:
            f_i_dT = 0
        else:
            f_i_dT = dTc(r[i], M_centro)    

        centro.append(['--', 'INICIO',  i, r[i], P_centro, T_centro, L_centro, M_centro])
        f_i_centro.append(['--', 'INICIO', i, r[i], " - ", f_i_dT, f_i_dL, f_i_dM])

    print(tabulate(centro, headers='firstrow', tablefmt='fancy_grid'))
    print(tabulate(f_i_centro, headers='firstrow', tablefmt='fancy_grid'))
#
#
# Para poder continuar: 
# f_i_centro[1][5] = 0 
#
#
#----------------------------------------------------------------------------
#
#
# 
#----------------------------------------------------------------------------
#ALGORITMO FASE A.2. Núcleo convectivo
#----------------------------------------------------------------------------
    centro_integracionA2 = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]
    jj = 1
    i = 3
    print("K buena", K)


    loop1 = True
    while loop1:
        print("\n * Calculando capa número",i,"\n")
        print("Ejecutando paso 1")
        print("Ejecutando paso 2bis")
        delta_1_iT = h * f_i_centro[i][5] - h * f_i_centro[i-1][5]
        Test = centro[i][5] + h * f_i_centro[i][5] + (1/2) * delta_1_iT

        loop2 = True
        while loop2:
            print("Ejecutando Polítropo")
            Pest = K * (Test ** 2.5)

            print("Ejecutando paso 3")
            print("Test", Test, "Pest", Pest)
            dMa = dMc(Test, r[i])
            print("dMa", dMa)
            delta_1_i1M = h * dMa - h * f_i_centro[i][7]
            Mcal = centro[i][7] + h * dMa - (1/2) * delta_1_i1M  
            print("Mcal", Mcal)

            print("Ejecutando paso 7bis")
            dTe = dTc(r[i], Mcal)
            delta_1_i1T = h * dTe - h * f_i_centro[i][5]
            Tcal = centro[i][5] + h * dTe - (1/2) * delta_1_i1T
            print("Tcal", Tcal)

            ErelT = abs(Tcal - Test) / Tcal
            print("ErelT",ErelT)
            jj += 1
            if jj == 100:
                loop2 = False

            if ErelT < Erel_max:
                loop2 = False
            else:
                Test = Tcal
#         loop2 = False #AÑADIDO AHORA
        Pcal = K * (Tcal ** 2.5)
    

        print("Ejecutando paso 6")
        Tgen = Tcal * 10 #Para que esté en K/e6, unidades de la funcion de gener_energia (aquí K/e7)
        k = 1
        gener_energia(Tgen)
        if k == 0:
            Cl = Cl1
            nu = nu_pp
            E = 'PP'
        else:
            Cl = Cl2
            nu = nu_CN
            E = 'CN'
        print("Cl definitivo", Cl, "nu definitivo", nu)

        dLu = dLc(r[i], Tcal, nu, Cl)
        delta_1_i1L = h * dLu - h * f_i_centro[i][6]
        delta_2_i1L = h * dLu - 2 *  h * f_i_centro[i][6] + h * f_i_centro[i-1][6]

        Lcal = centro[i][6] + h * dLu - (1/2) * delta_1_i1L - (1/12) * delta_2_i1L

        K = Pcal / (Tcal ** 2.5) # Constante del polítropo         
        print("K válida", K)

        if r[i] > r_Interp:
            loop1 = False  
            print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.2. -- \n")
            centro_integracionA2.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
            centro.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
            f_i_centro.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
            n = 0
            i += 1
            print(tabulate(centro_integracionA2, headers='firstrow', tablefmt='fancy_grid'))         
        else:
        #pasamos a la siguiente capa
            centro_integracionA2.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
            centro.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
            f_i_centro.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
            n = 0
            i += 1
# loop1 = False # AÑADIDO AHORA 
#
#
#
# VALORES EN LA FRONTERA (interpolando linealmente para cuando n + 1 = 2.5)
#----------------------------------------------------------------------------
# Valores en la zona radiactiva
    radio_Antes25_cen = r[18]
    P_Antes25_cen = centro_integracionA2[16][4]
    T_Antes25_cen = centro_integracionA2[16][5]
    L_Antes25_cen = centro_integracionA2[16][6] 
    M_Antes25_cen = centro_integracionA2[16][7]
    n_Antes25_cen = integracionA13[22][8]

# Valores en la zona convectiva
    radio_Despues25_cen = r[19]
    P_Despues25_cen = centro_integracionA2[17][4]
    T_Despues25_cen = centro_integracionA2[17][5]
    L_Despues25_cen = centro_integracionA2[17][6]
    M_Despues25_cen = centro_integracionA2[17][7]
    n_Despues25_cen = integracionA2[1][8]

    from scipy.optimize import curve_fit

    def poldeg1(x,a,b):
        return a * x + b

    equis_cen = [n_Antes25_cen,n_Despues25_cen]

    val, cov = curve_fit(poldeg1,equis_cen,[radio_Antes25_cen,radio_Despues25_cen])
    a,b = val
    r_Interp_cen = poldeg1(2.5,a,b)
#print("r_Interp",r_Interp)

    val, cov = curve_fit(poldeg1,equis_cen,[P_Antes25_cen,P_Despues25_cen])
    a,b = val
    P_Interp_cen = poldeg1(2.5,a,b)
# print("P_Interp",P_Interp)

    val, cov = curve_fit(poldeg1,equis_cen,[T_Antes25_cen,T_Despues25_cen])
    a,b = val
    T_Interp_cen = poldeg1(2.5,a,b)
# print("T_Interp",T_Interp)

    val, cov = curve_fit(poldeg1,equis_cen,[L_Antes25_cen,L_Despues25_cen])
    a,b = val
    L_Interp_cen = poldeg1(2.5,a,b)
# print("L_Interp",L_Interp)

    val, cov = curve_fit(poldeg1,equis_cen,[M_Antes25_cen,M_Despues25_cen])
    a,b = val
    M_Interp_cen = poldeg1(2.5,a,b)
# print("M_Interp",M_Interp)



    interpolacion_cen = [['r', 'P', 'T', 'L', 'M']]
    interpolacion_cen.append([r_Interp_cen, P_Interp_cen, T_Interp_cen, L_Interp_cen, M_Interp_cen])

    print(tabulate(interpolacion_cen, headers='firstrow', tablefmt='fancy_grid')) 

#----------------------------------------------------------------------------
#
#
#
#----------------------------------------------------------------------------
# Error relativo total

    def TotalRelEror(P_Interp, P_Interp_cen, T_Interp, T_Interp_cen, L_Interp, L_Interp_cen,
    M_Interp, M_Interp_cen):
    
        a = ((P_Interp - P_Interp_cen) / P_Interp) ** 2
        b = ((T_Interp - T_Interp_cen) / T_Interp) ** 2
        c = ((L_Interp - L_Interp_cen) / L_Interp) ** 2
        d = ((M_Interp - M_Interp_cen) / M_Interp) ** 2
        return np.sqrt(a + b + c + d)

    TotRE = TotalRelEror(P_Interp, P_Interp_cen, T_Interp, T_Interp_cen, L_Interp, L_Interp_cen,
    M_Interp, M_Interp_cen)
    print("Error relativo total = {:.4f}".format(100*TotRE))

    error_relativo.append(TotRE)
    Temp_c_plot.append(T_c)
#----------------------------------------------------------------------------
#
#
#
import matplotlib.pyplot as plt

plt.plot(Temp_c_plot, error_relativo)
x_min = np.argmin(error_relativo)
print(x_min, Temp_c_plot[x_min], error_relativo[x_min])

plt.plot(Temp_c_plot[x_min], error_relativo[x_min],'ro')


print(tabulate(integracion, headers='firstrow', tablefmt='fancy_grid'))         
# plt.show()