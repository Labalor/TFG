from B_funciones import *
from numpy import linspace
import warnings
from tabulate import tabulate
from copy import copy

warnings.filterwarnings("ignore")

R_tot_inicial = 11.629113924050632 
L_tot_inicial = 45.2506329113924
T_centro = 1.8603603603603602  

# 1 iteracion:  
# 2 iteracion: 
# 3 iteracion: 3.1178905511325166 
# 4 iteracion: 2.31011772281779 
# 5 iteracion: 1.7780332999027155 
# 6 iteracion: 1.2374477877105639
# 7 iteracion: 0.9653478031581336 
# 8 iteracion: 0.6966399145461489
# 9 iteracion: 0.4217028845084499  

def desde_fuera(R_tot, L_tot):

    #RADIO INICIAL Y RADIO DE CADA CAPA
    #----------------------------------------------------------------------------
    R_ini = 0.9 * R_tot       #Para evitar problemas de convergencia 
    h = - R_ini / 100         #Paso de integración
    r = linspace(R_ini,0,101) #Radio 
    #----------------------------------------------------------------------------

    #PRIMERAS TRES CAPAS
    #----------------------------------------------------------------------------
    # Creo por primera vez los dos vectores: 
    integracion = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M', 'n+1']]
    f_i = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]

    i = 0
    while i <= 2: 
        T_temp = T(r[i], R_tot)
        P_pres = P(T_temp,L_tot)
        E = ciclo(T_temp)

        f_i_dT = dT(r[i], T_temp, P_pres, L_tot)
        f_i_dP = dP(r[i], T_temp, P_pres, M_tot)
        f_i_dM = 0
        f_i_dL = 0

        integracion.append([E, 'INICIO',  i, r[i], P_pres, T_temp, L_tot, M_tot, 0])
        f_i.append([E, 'INICIO', i, r[i], f_i_dP, f_i_dT, f_i_dL, f_i_dM])

        i += 1
    #----------------------------------------------------------------------------
    #
    # A.1.1.
    #
    loop1 = True
    while loop1:
        
        Pest = presion_estimada(integracion[i][4], f_i, h, i)
        Test = temperatura_estimada(integracion[i][5], f_i, h, i)

        loop2 = True
        while loop2:
            loop3 = True
            while loop3:

                dPr = dP(r[i], Test, Pest, M_tot)
                Pcal = PTM_calculada(integracion[i][4], dPr, f_i[i][4], h, i)

                ErelP = error_relativo(Pcal, Pest)

                if ErelP < Erel_max:
                    loop3 = False
                
                Pest = Pcal

            dTe = dT(r[i], Test, Pcal, L_tot)
            Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)

            if ErelT < Erel_max:
                loop2 = False
            
            Test = Tcal

        dMa = dM(r[i], Tcal, Pcal)
        Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)

        ErelM = error_relativo(M_tot, Mcal)

        if ErelM > Erel_max :
            loop1 = False          
      
        else:
            #pasamos a la siguiente capa
            integracion.append([ciclo(Tcal), 'A.1.1.',  i, r[i], Pcal, Tcal, L_tot, Mcal, 0])
            f_i.append([ciclo(Tcal), 'A.1.1.',  i, r[i], dPr, dTe, 0, dMa])

            i += 1

    # A.1.2.

    loop1 = True
    while loop1:

        Pest = presion_estimada(integracion[i][4], f_i, h, i)
        Test = temperatura_estimada(integracion[i][5], f_i, h, i)

        loop2 = True
        while loop2:
            loop3 = True
            while loop3:

                dMa = dM(r[i], Test, Pest)
                Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)

                dPr = dP(r[i], Test, Pest, Mcal)
                Pcal = PTM_calculada(integracion[i][4], dPr, f_i[i][4], h, i)

                ErelP = error_relativo(Pcal, Pest)
                if ErelP < Erel_max:
                    loop3 = False
                
                Pest = Pcal

            dTe = dT(r[i], Test, Pcal, L_tot)
            Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)

            if ErelT < Erel_max:
                loop2 = False
            
            Test = Tcal

    #    print("Ejecutando paso 6")
        dLu = dL(r[i], Tcal, Pcal)
        Lcal = L_calculada(L_tot, dLu, f_i, h, i)

        ErelL = error_relativo(Lcal, L_tot)

        if ErelL > Erel_max :
            loop1 = False  

        else:
            #pasamos a la siguiente capa
            integracion.append([ciclo(Tcal), 'A.1.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal, 0])
            f_i.append([ciclo(Tcal), 'A.1.2.',  i, r[i], dPr, dTe, dLu, dMa])

            i += 1


    # A.1.3.
    loop1 = True
    while loop1:

        Pest = presion_estimada(integracion[i][4], f_i, h, i)
        Test = temperatura_estimada(integracion[i][5], f_i, h, i)


        loop2 = True
        while loop2:
            loop3 = True
            while loop3:

                dMa = dM(r[i], Test, Pest)
                Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)

                dPr = dP(r[i], Test, Pest, Mcal)
                Pcal = PTM_calculada(integracion[i][4], dPr, f_i[i][4], h, i)

                ErelP = error_relativo(Pcal, Pest)
                if ErelP < Erel_max:
                    loop3 = False
                
                Pest = Pcal

            dLu = dL(r[i], Test, Pcal)
            Lcal = L_calculada(integracion[i][6], dLu, f_i, h, i)

            dTe = dT(r[i], Test, Pcal, Lcal)
            Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)

            if ErelT < Erel_max:
                loop2 = False
            
            Test = Tcal


        dPr = dP(r[i], Tcal, Pcal, Mcal)
        dTe = dT(r[i], Tcal, Pcal, Lcal)
        n = n_fun(Tcal, Pcal, dPr, dTe)

        if n < 2.5:
            loop1 = False  
            K = Pcal / (Tcal ** 2.5) # Constante del polítropo         
        else:
            #pasamos a la siguiente capa
            dLu = dL(r[i], Tcal, Pcal)
            dMa = dM(r[i], Tcal, Pcal)
            integracion.append([ciclo(Tcal), 'A.1.3.',  i, r[i], Pcal, Tcal, Lcal, Mcal, n])
            f_i.append([ciclo(Tcal), 'A.1.3.',  i, r[i], dPr, dTe, dLu, dMa])

            i += 1


    #ALGORITMO FASE A.2. Núcleo convectivo
    #----------------------------------------------------------------------------
    tam = len(integracion) - 1

    loop1 = True
    while loop1:

        Test = temperatura_estimada(integracion[i][5], f_i, h, i)

        loop2 = True
        while loop2:

            Pest = Pcentro(r[i], K, Test)

            dMa = dM(r[i], Test, Pest)
            Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)

            if r[i] == 0:
                Tcal = Test
            else:
                dTe = dTc(r[i], Mcal)
                Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)
     

            if ErelT < Erel_max:
                loop2 = False
            
            Test = Tcal

        Pcal = Pcentro(r[i], K, Tcal)

        dLu = dL(r[i], Tcal, Pcal)
        Lcal = L_calculada(Lcal, dLu, f_i, h, i)


        if r[i] == 0:
            loop1 = False  

        integracion.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal, n])
        f_i.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        n = 0

        i += 1

    return integracion, tam, K



def desde_dentro(R_tot, T_centro, K, r_inter):

    R_ini = 0.9 * R_tot     #Para evitar problemas de convergencia 
    h = R_ini / 100
    r = linspace(R_ini,0,101)

    centro = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]
    f_i_centro = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]

    r = r[::-1]

    i = 0
    while i <= 2: 
        Mcal = Mcentro(r[i], K, T_centro)
        Lcal = Lcentro(r[i], K, T_centro)
        Tcal = Tcentro(r[i], K, T_centro)
        Pcal = Pcentro(r[i], K, Tcal)

        f_i_dM = dMc(Tcal, r[i], K)
        f_i_dL = dLc(r[i], Tcal, K)
        f_i_dP = dPc(Tcal,r[i], Mcal, K)
        f_i_dT = dTc(r[i], Mcal)    

        centro.append(['--', 'INICIO',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i_centro.append(['--', 'INICIO', i, r[i], f_i_dP, f_i_dT, f_i_dL, f_i_dM])

        i += 1

    loop1 = True
    while loop1:

        Test = temperatura_estimada(centro[i][5], f_i_centro, h, i)

        loop2 = True
        while loop2:

            Pest = Pcentro(r[i], K, Test)

            dMa = dMc(Test, r[i], K)
            Mcal = PTM_calculada(centro[i][7], dMa, f_i_centro[i][7], h, i)

            dTe = dTc(r[i], Mcal)
            Tcal = PTM_calculada(centro[i][5], dTe, f_i_centro[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)

            if ErelT < Erel_max:
                loop2 = False

            Test = Tcal

        Pcal = Pcentro(r[i], K, Tcal)   

        dLu = dLc(r[i], Tcal, K)
        Lcal = L_calculada(Lcal, dLu, f_i_centro, h, i)

        dPr = dPc(Tcal,r[i], Mcal, K)
        dTe = dTc(r[i], Mcal)
        dMa = dMc(Tcal, r[i], K)
          
        if r[i] > r_inter:
            loop1 = False  

        centro.append([ciclo(Tcal), 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i_centro.append([ciclo(Tcal), 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        i += 1

    return centro


calculos = desde_fuera(R_tot_inicial, L_tot_inicial)

def error_total_funcion(R_tot, L_tot, Tc):

    global desde_fuera, desde_dentro, calculos, R_tot_inicial, L_tot_inicial

    calculos_2 = copy(calculos)
    if L_tot != L_tot_inicial or R_tot != R_tot_inicial:
        calculos_2 = desde_fuera(R_tot, L_tot)

    integracion = calculos_2[0]
    tam = calculos_2[1]
    K = calculos_2[2]

    n_rad = integracion[tam][8]
    n_conv = integracion[tam+1][8]

    r_rad = integracion[tam][3]
    r_conv = integracion[tam+1][3]
    r_inter = interpolation(r_rad,r_conv,n_rad,n_conv)

    P_rad = integracion[tam][4]
    P_conv = integracion[tam+1][4]
    P_inter = interpolation(P_rad,P_conv,n_rad,n_conv)

    T_rad = integracion[tam][5]
    T_conv = integracion[tam+1][5]
    T_inter = interpolation(T_rad,T_conv,n_rad,n_conv)

    L_rad = integracion[tam][6]
    L_conv = integracion[tam+1][6]
    L_inter = interpolation(L_rad,L_conv,n_rad,n_conv)

    M_rad = integracion[tam][7]
    M_conv = integracion[tam+1][7]
    M_inter = interpolation(M_rad,M_conv,n_rad,n_conv)


    centro = desde_dentro(R_tot, Tc, K, r_inter)

    r_rad = centro[-1][3]
    r_conv = centro[-2][3]
    r_inter_centro = interpolation(r_rad,r_conv,n_rad,n_conv)
    
    P_rad = centro[-1][4]
    P_conv = centro[-2][4]
    P_inter_centro = interpolation(P_rad,P_conv,n_rad,n_conv)
    
    T_rad = centro[-1][5]
    T_conv = centro[-2][5]
    T_inter_centro = interpolation(T_rad,T_conv,n_rad,n_conv)
    
    L_rad = centro[-1][6]
    L_conv = centro[-2][6]
    L_inter_centro = interpolation(L_rad,L_conv,n_rad,n_conv)
    
    M_rad = centro[-1][7]
    M_conv = centro[-2][7]
    M_inter_centro = interpolation(M_rad,M_conv,n_rad,n_conv)
    

    return 100 * TotalRelEror(P_inter, P_inter_centro, T_inter, T_inter_centro, L_inter, L_inter_centro, M_inter, M_inter_centro)

