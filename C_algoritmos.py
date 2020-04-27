
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
from B_funciones import *
from numpy import linspace
import warnings
from tabulate import tabulate
from copy import copy
#------------------------------------------------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------------------------------------------------
# INTEGRACIÓN DESDE LA SUPERFICIE
#------------------------------------------------------------------------------------------------------------------------
def desde_fuera(R_tot, L_tot):
	'''
	Cálculo de la integración desde la superficie: zona radioactiva
	Input: radio total y luminosidad total
	Output: lista de las magnitudes calculadas, su tamaño y la constante del polítropo
	'''
	#--------------------------------------------------------------------------------------------------------------------
    # RADIO INICIAL Y RADIO DE CADA CAPA
	#--------------------------------------------------------------------------------------------------------------------
    R_ini = 0.9 * R_tot 		# Radio inicial, para evitar problemas de convergencia
    h = - R_ini / 100         	# Paso de integración
    r = linspace(R_ini,0,101) 	# Radio relativo a cada capa 
	#--------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------
	# ALMACENAMIENTO DE DATOS
	#--------------------------------------------------------------------------------------------------------------------
    # Se crean las listas donde se almacenarán los valores de las magnitudes y las ec. diferenciales: 
    integracion = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M', 'n+1']]	# Valores de las magnitudes
    f_i = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]				# Valores de las ec. diferenciales
    #--------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------
    # PRIMERAS TRES CAPAS
    #--------------------------------------------------------------------------------------------------------------------
    i = 0
    while i <= 2: 
   		# El valor de la luminosidad y la masa es el valor total
        T_temp = T(r[i], R_tot)		# Temperatura
        P_pres = P(T_temp,L_tot)	# Presión
        E = ciclo(T_temp)			# Ciclo de generación de energía

        # Se calculan las ecuaciones diferenciales 
        f_i_dT = dT(r[i], T_temp, P_pres, L_tot)	# Temperatura
        f_i_dP = dP(r[i], T_temp, P_pres, M_tot)	# Presión
        f_i_dM = 0									# Masa
        f_i_dL = 0									# Luminosidad

        # Se almacenan los valores generados en las listas 
        # Nota: el 0 (último elemento de 'integración') se añade por completitud
        integracion.append([E, 'INICIO',  i, r[i], P_pres, T_temp, L_tot, M_tot, 0])
        f_i.append([E, 'INICIO', i, r[i], f_i_dP, f_i_dT, f_i_dL, f_i_dM])

        # Pasa a la siguiente capa
        i += 1
    #--------------------------------------------------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------------------------------------------------
    # ALGORITMO FASE A.1.1. (envoltura radioactiva)
    #--------------------------------------------------------------------------------------------------------------------
    loop1 = True 												
    while loop1:
        Pest = presion_estimada(integracion[i][4], f_i, 		# Presión estimada 
        	h, i)
        Test = temperatura_estimada(integracion[i][5], 			# Temperatura estimada
        	f_i, h, i)

        loop2 = True
        while loop2:
            loop3 = True
            while loop3:

                dPr = dP(r[i], Test, Pest, M_tot)				# Ec. diferencial de la presión
                Pcal = PTM_calculada(integracion[i][4], 		# Presión calculada
                	dPr, f_i[i][4], h, i)

                ErelP = error_relativo(Pcal, Pest)				# Error relativo entre presiones

                if ErelP < Erel_max:							# Se repite el loop3 hasta conseguir un 
                    loop3 = False								# error relativo menor que el máximo
                
                Pest = Pcal 									# Para el error del siguiente loop

            dTe = dT(r[i], Test, Pcal, L_tot)					# Ec. diferencial de la temperatura
            Tcal = PTM_calculada(integracion[i][5], dTe, 		# Temperatura calculada
            	f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)					# Error relativo entre temperaturas

            if ErelT < Erel_max:								# Se repite el loop2 hasta conseguir un 
                loop2 = False									# error relativo menor que el máximo
            
            Test = Tcal 										# Para el error del siguiente loop

        dMa = dM(r[i], Tcal, Pcal)								# Ec. diferencial de la masa
        Mcal = PTM_calculada(integracion[i][7], dMa, 			# Masa calculada
        	f_i[i][7], h, i)

        ErelM = error_relativo(M_tot, Mcal)						# Error relativo entre masas
        														# La masa estimada es la masa total				
        
        if ErelM > Erel_max :									# se repite el loop1 hasta conseguir un 
            loop1 = False          								# error relativo menor que el máximo
      
        else:													# Si se consigue, se pasa a la siguiente capa

            integracion.append([ciclo(Tcal), 'A.1.1.',  		# Almacenamiento de valores calculados
            	i, r[i], Pcal, Tcal, L_tot, Mcal, 0])
            f_i.append([ciclo(Tcal), 'A.1.1.',  i, r[i], 
            	dPr, dTe, 0, dMa])

            i += 1												# Se pasa a la siguiente capa
    #--------------------------------------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------------------------------------------
    # ALGORITMO FASE A.1.2. (envoltura radioactiva)
    #--------------------------------------------------------------------------------------------------------------------
    loop1 = True
    while loop1:

        Pest = presion_estimada(integracion[i][4], f_i, 		# Presión estimada
        	h, i)
        Test = temperatura_estimada(integracion[i][5], 			# Temperatura estimada
        	f_i, h, i)

        loop2 = True
        while loop2:
            loop3 = True
            while loop3:

                dMa = dM(r[i], Test, Pest)						# Ec. diferencial de la masa
                Mcal = PTM_calculada(integracion[i][7], 		# Masa calculada
                	dMa, f_i[i][7], h, i)

                dPr = dP(r[i], Test, Pest, Mcal)				# Ec. diferencial de la presión
                Pcal = PTM_calculada(integracion[i][4], 		# Presión calculada
                	dPr, f_i[i][4], h, i)

                ErelP = error_relativo(Pcal, Pest)				# Error relativo entre presiones
                
                if ErelP < Erel_max:							# Se repite el loop3 hasta conseguir un
                    loop3 = False								# error relativo menor que el máximo
                
                Pest = Pcal 									# Para el error del siguiente loop

            dTe = dT(r[i], Test, Pcal, L_tot)					# Ec. diferencial de la temperatura
            Tcal = PTM_calculada(integracion[i][5], dTe, 		# Temperatura calculada
            	f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)					# Error relativo entre temperaturas

            if ErelT < Erel_max:								# Se repite el loop2 hasta conseguir un 
                loop2 = False									# error relativo menor que el máximo
            
            Test = Tcal 										# Para el error del siguiente loop

        dLu = dL(r[i], Tcal, Pcal)								# Ec. diferencial para la luminosidad
        Lcal = L_calculada(L_tot, dLu, f_i, h, i)				# Luminosidad calculada

        ErelL = error_relativo(Lcal, L_tot)						# Error relativo entre luminosidades

        if ErelL > Erel_max :									# Se repite el loop1 hasta conseguir un
            loop1 = False  										# error relativo menor que el máximo

        else:													# Si se consigue, se pasa a la siguiente capa
            
            integracion.append([ciclo(Tcal), 'A.1.2.',  		# Almacenamiento de valores calculados  
            	i, r[i], Pcal, Tcal, Lcal, Mcal, 0])
            f_i.append([ciclo(Tcal), 'A.1.2.',  i, r[i], 
            	dPr, dTe, dLu, dMa])

            i += 1												# Se pasa a la siguiente capa
    #--------------------------------------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------------------------------------------
    # ALGORITMO FASE A.1.3. (envoltura radioactiva)
    #--------------------------------------------------------------------------------------------------------------------
    loop1 = True
    while loop1:

        Pest = presion_estimada(integracion[i][4], f_i, 		# Presión estimada
        	h, i)
        Test = temperatura_estimada(integracion[i][5], 			# Temperatura estimada
        	f_i, h, i)


        loop2 = True
        while loop2:
            loop3 = True
            while loop3:

                dMa = dM(r[i], Test, Pest)						# Ec. diferencial de la masa 
                Mcal = PTM_calculada(integracion[i][7], 		# Masa calculada
                	dMa, f_i[i][7], h, i)

                dPr = dP(r[i], Test, Pest, Mcal)				# Ec. diferencial de la presión
                Pcal = PTM_calculada(integracion[i][4], 		# Presión calculada
                	dPr, f_i[i][4], h, i)

                ErelP = error_relativo(Pcal, Pest)				# Error relativo entre presiones

                if ErelP < Erel_max:							# Se repite el loop3 hasta conseguir un 
                    loop3 = False								# error relativo menor que el máximo
                
                Pest = Pcal 									# Para el error del siguiente loop

            dLu = dL(r[i], Test, Pcal)							# Ec. diferencial de la luminosidad
            Lcal = L_calculada(integracion[i][6], dLu, 			# Luminosidad calculada
            	f_i, h, i)

            dTe = dT(r[i], Test, Pcal, Lcal)					# Ec. diferencial de la temperatura
            Tcal = PTM_calculada(integracion[i][5], dTe, 		# Temperatura calculada
            	f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)					# Error relativo entre temperaturas

            if ErelT < Erel_max:								# Se repite el loop2 hasta conseguir un 
                loop2 = False									# error relativo menor que el máximo
            
            Test = Tcal 										# Para el error del siguiente loop


        #dPr = dP(r[i], Tcal, Pcal, Mcal)
        #dTe = dT(r[i], Tcal, Pcal, Lcal)
        n = n_fun(Tcal, Pcal, dPr, dTe) 						# Cálculo del parámetro n+1

        if n < 2.5: 											# Se repite el loop1 hasta conseguir un 
            loop1 = False  										# error relativo menor que el máximo
            K = Pcal / (Tcal ** 2.5) 							# Se calcula la constante del polítropo         
        
        else:													# Si se consigue, se pasa a la siguiente capa
            # dLu = dL(r[i], Tcal, Pcal)
            # dMa = dM(r[i], Tcal, Pcal)
            integracion.append([ciclo(Tcal), 'A.1.3.',  		# Almacenamiento de valores calculados
            	i, r[i], Pcal, Tcal, Lcal, Mcal, n])
            f_i.append([ciclo(Tcal), 'A.1.3.',  i, r[i], 
            	dPr, dTe, dLu, dMa])

            i += 1												# Se pasa a la siguiente capa
    #--------------------------------------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------------------------------------------
    # ALGORITMO FASE A.2 (núcleo convectivo)
    #--------------------------------------------------------------------------------------------------------------------
    tam = len(integracion) - 1 									# Calculo del tamaño de las capas radioactivas

    loop1 = True
    while loop1:

        Test = temperatura_estimada(integracion[i][5], 			# Temperatura estimada
        	f_i, h, i)

        loop2 = True
        while loop2:

            Pest = Pcentro(r[i], K, Test)						# Presión estimada

            dMa = dM(r[i], Test, Pest)							# Ec. diferencial de la masa
            Mcal = PTM_calculada(integracion[i][7], 			# Masa calculada
            	dMa, f_i[i][7], h, i)

            if r[i] == 0:										# Para el caso de radio nulo
                Tcal = Test 									# La temperatura calculada es la estimada
            
            else:												# Si no, se procede con normalidad:
                
                dTe = dTc(r[i], Mcal)							# Ec. diferencial de la temperatura
                Tcal = PTM_calculada(integracion[i][5], 		# Temperatura calculada
                	dTe, f_i[i][5], h, i)

            ErelT = error_relativo(Tcal, Test) 					# Error relativo entre temperaturas

            if ErelT < Erel_max:								# Se repite el loop2 hasta conseguir un
                loop2 = False									# error relativo menor que el máximo
            
            Test = Tcal 										# Para el error del siguiente loop

        Pcal = Pcentro(r[i], K, Tcal) 							# Presión calculada

        dLu = dL(r[i], Tcal, Pcal) 								# Ec. diferencial de la luminosidad
        Lcal = L_calculada(Lcal, dLu, f_i, h, i)				# Luminosidad calculada

        if r[i] == 0: 											# Cuando se llega a la última capa de la 
            loop1 = False  										# estrella, se detiene el loop1

        integracion.append([E, 'A.2.',  i, r[i], Pcal, 			# Almacenamiento de valores calculados
        	Tcal, Lcal, Mcal, n])
        f_i.append([E, 'A.2.',  i, r[i], dPr, dTe, dLu, 
        	dMa])
        n = 0 											

        i += 1													# Se pasa a la siguiente capa

    return integracion, tam, K
    #--------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# INTEGRACIÓN DESDE EL CENTRO: ZONA CONVECTIVA
#------------------------------------------------------------------------------------------------------------------------
def desde_dentro(R_tot, T_centro, K, r_inter):
	'''
	Cálculo de la integración desde el centro: zona convectiva
	Input: radio total, temperatura centra, constante del polítropo y radio interpolado de la zona radioactiva
	Output: lista de las magnitudes calculadas 
	'''

	#--------------------------------------------------------------------------------------------------------------------
	# RADIO INICIAL Y RADIO DE CADA CAPA
	#--------------------------------------------------------------------------------------------------------------------
    R_ini = 0.9 * R_tot     	# Radio inicial, para evitar problemas de convergencia
    h = R_ini / 100				# Paso de integración, invertido
    r = linspace(R_ini,0,101)	# Radio relativo a cada capa
    r = r[::-1]					# Inversión del radio relativo a cada capa
    #--------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------
	# ALMACENAMIENTO DE DATOS
	#--------------------------------------------------------------------------------------------------------------------
    centro = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M']]
    f_i_centro = [['E', 'FASE', 'i', 'r', 'dP', 'dT', 'dL', 'dM']]
	#--------------------------------------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------------------------------------------
    #ALGORITMO FASE A.2. Núcleo convectivo
    #--------------------------------------------------------------------------------------------------------------------	
    i = 0 														# Se establece la capa número cero

    while i <= 2: 												# Cálculo de las tres primeras capas
        Mcal = Mcentro(r[i], K, T_centro)						# Masa calculada
        Lcal = Lcentro(r[i], K, T_centro)						# Luminosidad calculada
        Tcal = Tcentro(r[i], K, T_centro)						# Temperatura calculada
        Pcal = Pcentro(r[i], K, Tcal)							# Presión calculada

        f_i_dM = dMc(Tcal, r[i], K) 							# Ec. diferencial de la masa
        f_i_dL = dLc(r[i], Tcal, K) 							# Ec. diferencial de la luminosidad
        f_i_dP = dPc(Tcal,r[i], Mcal, K)						# Ec. diferencial de la presión
        f_i_dT = dTc(r[i], Mcal)    							# Ec. diferencial de la temperatura

        centro.append(['--', 'CENTRO',  i, r[i], Pcal, 			# Almacenamiento de valores calculados
        	Tcal, Lcal, Mcal])
        f_i_centro.append(['--', 'CENTRO', i, r[i], 
        	f_i_dP, f_i_dT, f_i_dL, f_i_dM])

        i += 1													# Se pasa a la siguiente capa

    loop1 = True
    while loop1:

        Test = temperatura_estimada(centro[i][5], 				# Temperatura estimada
        	f_i_centro, h, i)

        loop2 = True
        while loop2:

            Pest = Pcentro(r[i], K, Test)						# Presión estimada

            dMa = dMc(Test, r[i], K)							# Ec. diferencial de la masa
            Mcal = PTM_calculada(centro[i][7], dMa, 			# Masa calculada
            	f_i_centro[i][7], h, i)

            dTe = dTc(r[i], Mcal)								# Ec. diferencial de la temperatura
            Tcal = PTM_calculada(centro[i][5], dTe, 			# Temperatura calculada
            	f_i_centro[i][5], h, i)

            ErelT = error_relativo(Tcal, Test)					# Error relativo entre temperaturas

            if ErelT < Erel_max:								# Se repite el loop2 hasta conseguir un
                loop2 = False									# error mínimo menor que el máximo

            Test = Tcal 										# Para el error del siguiente loop

        Pcal = Pcentro(r[i], K, Tcal)   						# Presión calculada

        dLu = dLc(r[i], Tcal, K) 								# Ec. diferencial de la luminosidad
        Lcal = L_calculada(Lcal, dLu, f_i_centro, 				# Luminosidad calculada
        	h, i)

        # dPr = dPc(Tcal,r[i], Mcal, K) 
        # dTe = dTc(r[i], Mcal)
        # dMa = dMc(Tcal, r[i], K)
          
        if r[i] > r_inter: 										# Si se llega a la zona radioactiva
            loop1 = False  										# se detiene el loop1

        centro.append([ciclo(Tcal), 'A.2.',  i, 				# Almacenamiento de valores calculados
        	r[i], Pcal, Tcal, Lcal, Mcal])
        f_i_centro.append([ciclo(Tcal), 'A.2.',  
        	i, r[i], dPr, dTe, dLu, dMa])

        i += 1													# Se pasa a la siguiente capa

    return centro
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# Almacenamiento de outputs de 'INTEGRACIÓN DESDE LA SUPERFICIE: ZONA RADIOACTIVA' en variable local
#------------------------------------------------------------------------------------------------------------------------
# BORRAR: se calcula fuera de la función 'desde_fuera' debido a que no depende de la temperatura central. Sin embargo, 
# cuando introducimos los valores totales de radio y luminosidad, tengo que recalcular las capas. Por ello, crearé
# un condicional en la línea 402
calculos = desde_fuera(R_tot_inicial, L_tot_inicial)
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# INTERPOLACIÓN Y ERROR RELATIVO TOTAL
#------------------------------------------------------------------------------------------------------------------------
def error_total_funcion(R_tot, L_tot, Tc):
	'''
	Cálculo de la interpolación en los valores frontera y estudio del error relativo total
	Input: radio total, luminosidad total y temperatura central
	Output: valor del error relativo total
	'''
    global desde_fuera, desde_dentro, calculos, R_tot_inicial, L_tot_inicial

    # Se crea una variable local con los datos de la variable global para no modificar los de la global. Si no, estaría 
    # modificando los datos de la variable local por datos que no son los iniciales
    calculos_2 = copy(calculos) 

    # Cuando se introducen los valores totales de radio y luminosidad, se recalculan las capas de la zona radiactiva
    if L_tot != L_tot_inicial or R_tot != R_tot_inicial:
        calculos_2 = desde_fuera(R_tot, L_tot)

	#--------------------------------------------------------------------------------------------------------------------
    # Almacenamiento de los outputs de la función 'desde_fuera' (ZONA RADIOACTIVA) en variables locales
    #--------------------------------------------------------------------------------------------------------------------
    integracion = calculos_2[0]					# Lista de las magnitudes calculadas
    tam = calculos_2[1]							# Tamaño de la lista de magnitudes calculadas
    K = calculos_2[2]							# Constante del polítropo

    n_rad = integracion[tam][8]					# n+1 frontera-radioactiva
    n_conv = integracion[tam+1][8]				# n+1 frontera-convectiva

    r_rad = integracion[tam][3]					# Radio frontera-radioactiva
    r_conv = integracion[tam+1][3]				# n+1 frontera-convectiva
    r_inter = interpolation(r_rad,				# Radio interpolado 
    	r_conv,n_rad,n_conv)

    P_rad = integracion[tam][4]					# Presión frontera-radioactiva
    P_conv = integracion[tam+1][4]				# Presión frontera-convectiva
    P_inter = interpolation(P_rad,				# Presión interpolado 
    	P_conv,n_rad,n_conv)

    T_rad = integracion[tam][5]					# Temperatura frontera-radioactiva
    T_conv = integracion[tam+1][5]				# Temperatura frontera-convectiva
    T_inter = interpolation(T_rad,				# Temperatura interpolado 
    	T_conv,n_rad,n_conv)

    L_rad = integracion[tam][6]					# Luminosidad frontera-radioactiva
    L_conv = integracion[tam+1][6]				# Luminosidad frontera-convectiva
    L_inter = interpolation(L_rad,				# Luminosidad interpolado 
    	L_conv,n_rad,n_conv)

    M_rad = integracion[tam][7]					# Masa frontera-radioactiva
    M_conv = integracion[tam+1][7]				# Masa frontera-convectiva
    M_inter = interpolation(M_rad,				# Masa interpolado 
    	M_conv,n_rad,n_conv)
	#--------------------------------------------------------------------------------------------------------------------

	#--------------------------------------------------------------------------------------------------------------------
    # Almacenamiento del output de la función 'desde_dentro' (ZONA CONVECTIVA) en una variable local
    #--------------------------------------------------------------------------------------------------------------------
    centro = desde_dentro(R_tot, 
    	Tc, K, r_inter)

    r_rad = centro[-1][3]						# Radio frontera-radioactiva
    r_conv = centro[-2][3]						# Radio frontera-convectiva
    r_inter_centro = interpolation(r_rad,		# Radio interpolado 
    	r_conv,n_rad,n_conv)
    
    P_rad = centro[-1][4]						# Presión frontera-radioactiva
    P_conv = centro[-2][4]						# Presión frontera-convectiva
    P_inter_centro = interpolation(P_rad,		# Presión interpolad
    	P_conv,n_rad,n_conv)
    
    T_rad = centro[-1][5]						# Temperatura frontera-radioactiva
    T_conv = centro[-2][5]						# Temperatura frontera-convectiva
    T_inter_centro = interpolation(T_rad,		# Temperatura interpolado
    	T_conv,n_rad,n_conv)
    
    L_rad = centro[-1][6]						# Luminosidad frontera-radioactiva
    L_conv = centro[-2][6]						# Luminosidad frontera-convectiva
    L_inter_centro = interpolation(L_rad,		# Luminosidad interpolado
    	L_conv,n_rad,n_conv)
    
    M_rad = centro[-1][7]						# Masa frontera-radioactiva
    M_conv = centro[-2][7]						# Masa frontera-convectiva
    M_inter_centro = interpolation(M_rad,		# Masa interpolado
    	M_conv,n_rad,n_conv)
    #--------------------------------------------------------------------------------------------------------------------

    return 100 * TotalRelEror(P_inter, P_inter_centro, T_inter, T_inter_centro, 
    	L_inter, L_inter_centro, M_inter, M_inter_centro)
#------------------------------------------------------------------------------------------------------------------------
