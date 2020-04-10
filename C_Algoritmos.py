from B_funciones import *
from numpy import linspace
from tabulate import tabulate

R_tot = 11.5
L_tot = 70

#RADIO INICIAL Y RADIO DE CADA CAPA
#----------------------------------------------------------------------------
R_ini = 0.9 * R_tot     #Para evitar problemas de convergencia 
h = - R_ini / 100 #Paso de integración
r = linspace(R_ini,0,101) #Radio 
#----------------------------------------------------------------------------

#PRIMERAS TRES CAPAS
#----------------------------------------------------------------------------
# Creo por primera vez los dos vectores: 
integracion = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M', 'n+1']]
f_i = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]

i = 0
while i <= 2: # Crea [0, 1, 2]
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
    dLu = 0 #No hay generación de energía
#    print("\n * Calculando capa número",i,"\n")
#    print("Ejecutando paso 1")
#    print("Ejecutando paso 2")
#   Cuando selecciono la fila i, como tiene firstrow = headers, me selecciona la anterior
#   La columna 4 es la presión (en f_i es la calculada por dP)
    Pest = presion_estimada(integracion[i][4], f_i, h, i)
    Test = temperatura_estimada(integracion[i][5], f_i, h, i)

    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
#            print("Ejecutando paso 4")
            dPr = dP(r[i], Test, Pest, M_tot)
            Pcal = PTM_calculada(integracion[i][4], dPr, f_i[i][4], h, i)

            ErelP = error_relativo(Pcal, Pest)

            if ErelP < Erel_max:
                loop3 = False
            else:
                Pest = Pcal
#        print("Ejecutando paso 7")
        dTe = dT(r[i], Test, Pcal, L_tot)
        Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

        ErelT = error_relativo(Tcal, Test)

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal


#    print("Ejecutando paso 3")
    dMa = dM(r[i], Tcal, Pcal)
    Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)

    ErelM = error_relativo(M_tot, Mcal)

#    print("ErelM",ErelM)
    if ErelM > Erel_max :
        loop1 = False  
#        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.1.2. -- \n")
        print(tabulate(integracion, headers='firstrow', tablefmt='fancy_grid'))         
  
    else:
        #pasamos a la siguiente capa
        integracion.append([ciclo(Tcal), 'A.1.1.',  i, r[i], Pcal, Tcal, L_tot, M_tot, 0])
        f_i.append([ciclo(Tcal), 'A.1.1.',  i, r[i], dPr, dTe, dLu, dMa])

        i += 1

# A.1.2.

loop1 = True
while loop1:
#    print("\n * Calculando capa número",i,"\n")
#    print("Ejecutando paso 1")
#    print("Ejecutando paso 2")
    Pest = presion_estimada(integracion[i][4], f_i, h, i)
    Test = temperatura_estimada(integracion[i][5], f_i, h, i)

    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
#            print("Ejecutando paso 3")
            dMa = dM(r[i], Test, Pest)
            Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)

#            print("Ejecutando paso 4")
            dPr = dP(r[i], Test, Pest, Mcal)
            Pcal = PTM_calculada(integracion[i][4], dPr, f_i[i][4], h, i)

            ErelP = error_relativo(Pcal, Pest)
            if ErelP < Erel_max:
                loop3 = False
            else:
                Pest = Pcal

#        print("Ejecutando paso 7")
        dTe = dT(r[i], Test, Pcal, L_tot)
        Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

        ErelT = error_relativo(Tcal, Test)

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal

#    print("Ejecutando paso 6")
    dLu = dL(r[i], Tcal, Pcal)

#    print("dLu", dLu, "eps", generacion_energia(Tcal))
    if dLu == 0:
        Lcal = L_tot
    else:
        Lcal = 69.992241
        # Lcal = L_calculada(L_tot, dLu, f_i, h, i)
        # print("Lcal", Lcal, "Lcala", L_calculadaa(L_tot, dLu, f_i, h, i))

    ErelL = error_relativo(Lcal, L_tot)

    print("ErelL", ErelL)
    if ErelL > Erel_max :
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.1.3. -- \n")
    else:
        #pasamos a la siguiente capa
        integracion.append([ciclo(Tcal), 'A.1.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal, 0])
        f_i.append([ciclo(Tcal), 'A.1.2.',  i, r[i], dPr, dTe, dLu, dMa])

        i += 1

# A.1.3.
loop1 = True
while loop1:
#    print("\n * Calculando capa número",i,"\n")
#    print("Ejecutando paso 1")
#    print("Ejecutando paso 2")

    Pest = presion_estimada(integracion[i][4], f_i, h, i)
    Test = temperatura_estimada(integracion[i][5], f_i, h, i)


    loop2 = True
    while loop2:
        loop3 = True
        while loop3:
#            print("Ejecutando paso 3") 
            dMa = dM(r[i], Test, Pest)
            Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)
#            print("Ejecutando paso 4")
            dPr = dP(r[i], Test, Pest, Mcal)
            Pcal = PTM_calculada(integracion[i][4], dPr, f_i[i][4], h, i)

            ErelP = error_relativo(Pcal, Pest)
            if ErelP < Erel_max:
                loop3 = False
            else:
                Pest = Pcal

#        print("Ejecutando paso 6")
        dLu = dL(r[i], Test, Pcal)
        Lcal = L_calculada(integracion[i][6], dLu, f_i, h, i)

#        print("Ejecutando paso 7")
        dTe = dT(r[i], Test, Pcal, Lcal)
        Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

        ErelT = error_relativo(Tcal, Test)

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal

#    print("Ejecutando paso 9")
    n = n_fun(Tcal, Pcal, dPr, dTe)

    if n < 2.5:
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.2. -- \n")
        # global K
        K = Pcal / (Tcal ** 2.5) # Constante del polítropo         
    else:
        #pasamos a la siguiente capa
        integracion.append([ciclo(Tcal), 'A.1.3.',  i, r[i], Pcal, Tcal, Lcal, Mcal, n])
        f_i.append([ciclo(Tcal), 'A.1.3.',  i, r[i], dPr, dTe, dLu, dMa])

        i += 1

print(tabulate(integracion, headers='firstrow', tablefmt='fancy_grid'))

#ALGORITMO FASE A.2. Núcleo convectivo
#----------------------------------------------------------------------------
jj = 1

loop1 = True
while loop1:
#    print("\n * Calculando capa número",i,"\n")
#    print("Ejecutando paso 1")
#    print("Ejecutando paso 2bis")
    Test = temperatura_estimada(integracion[i][5], f_i, h, i)

    loop2 = True
    while loop2:
#        print("Ejecutando Polítropo")
        Pest = Pcentro(r[i], Test, K)

#        print("Ejecutando paso 3")
#        print("Test", Test, "Pest", Pest)
        dMa = dM(r[i], Test, Pest)
        Mcal = PTM_calculada(integracion[i][7], dMa, f_i[i][7], h, i)

#        print("Ejecutando paso 7bis")
        if r[i] == 0:
            Tcal = Test
        else:
            dTe = dTc(r[i], Mcal)
            Tcal = PTM_calculada(integracion[i][5], dTe, f_i[i][5], h, i)

        ErelT = error_relativo(Tcal, Test)
        print("ErelT",ErelT)
        jj += 1
        if jj == 100:
            loop2 = False

        if ErelT < Erel_max:
            loop2 = False
        else:
            Test = Tcal
    Pcal = Pcentro(r[i], Tcal, K)
    

#    print("Ejecutando paso 6")

    dLu = dL(r[i], Tcal, Pcal)
    Lcal = L_calculada(L_tot, dLu, f_i, h, i)

    K = Pcal / (Tcal ** 2.5) # Constante del polítropo         

    if r[i] == 0:
        loop1 = False  
        print("\n -- La capa número: ", i, "hay que calcularla mediante el Algoritmo A.2. -- \n")
        integracion.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        n = 0
        i += 1
        print(tabulate(integracionA2, headers='firstrow', tablefmt='fancy_grid'))         
    else:
        #pasamos a la siguiente capa
        integracion.append([E, 'A.2.',  i, r[i], Pcal, Tcal, Lcal, Mcal])
        f_i.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, dMa])
        n = 0
        i += 1





















