
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
#                           			  UNIVERSIDAD COMPLUTENSE DE MADRID
# 
#                               				CURSO ACADÉMICO 2019-2020
#
########################################################################################################################
########################################################################################################################

#------------------------------------------------------------------------------------------------------------------------
# PAQUETES NECESARIOS
#------------------------------------------------------------------------------------------------------------------------
# from C_algoritmos import * 				# Se importan los algoritmos de C_Algoritmos
from D_estudio import * 				# Se importa el estudio completo 
import matplotlib.pyplot as plt 		# Se importa la biblioteca de matplotlib.pyplot como 'plt'
from numpy import arange				# Se importa arange del paquete numpy
#------------------------------------------------------------------------------------------------------------------------

# (*) Se carga C_Algoritmos y se introducen los valores de R, L y Tc manualmente 
'''
R = 11.6263157895 
L = 45.2526315789 
Tc = 1.86036036036 
'''

#------------------------------------------------------------------------------------------------------------------------
# CREACIÓN DE UNA ÚNICA LISTA DE INTEGRACIÓN POR CAPAS
#------------------------------------------------------------------------------------------------------------------------
def ayuda1(R, L):
	'''
	Estudio de las 10 primeras capas del modelo
	Input: Radio total y luminosidad total de la estrella
	Output: 10 primeras capas del modelo
	'''
	herramienta = [['E', 'FASE', 'i', 'r', 'P', 'T', 'L', 'M', 'n+1']] 
	r_inicio = linspace(R, 0.9 * R,11)

	i = 0
	while i <= 9: 
	    T_temp = T(r_inicio[i], R)
	    P_pres = P(T_temp, L)
	    E = ciclo(T_temp)

	    herramienta.append([E, '^^^^^^', -10+i, r_inicio[i], P_pres, T_temp, L, M_tot, 0])
	    i += 1
	return herramienta

def ayuda2(lista):
	'''
	Cambio del índice 1 sobre la integración desde dentro
	Input: lista de las magnitudes calculadas desde la integración desde dentro
	Output: lista ordenada de manera consecuente respecto a la integración desde fuera
	'''
	for pos in range(len(lista)):
		lista[pos][2] = 101 - len(lista) + pos
	return lista

# Se calcula la integración desde fuera a través de la función definida en 'C_Algoritmos'
valores_desde_fuera = desde_fuera(R, L)

# Se almacenan los outputs de la función 'desde_fuera' en variables locales
integracion = valores_desde_fuera[0] 	# Valores de la integración
tam = valores_desde_fuera[1]			# Tamaño de la lista 'integración'
K = valores_desde_fuera[2]				# Constante del polítropo

n_rad = integracion[tam][8]				# Valor n+1 en la capa radiactiva
n_conv = integracion[tam+1][8]			# Valor n+1 en la capa convectiva

r_rad = integracion[tam][3]				# Radio en la capa radiactiva
r_conv = integracion[tam+1][3]			# Radio en la capa convectiva
r_inter = interpolation(r_rad,			# Radio interpolado 
	r_conv,n_rad,n_conv)

# Se calcula la integración desde dentro a través de la función definida en 'C_Algoritmos'
# Y se almacena el único output de la función 'desde_dentro' en una variable local
centro = desde_dentro(R, Tc, K, r_inter)


# Se preparan las listas calculadas en ambas integraciones
integracion_plot = integracion[1:tam+1] 	# Se excluye la primera fila (el encabezado)
centro_plot = centro[::-1][1:-1] 			# Se invierte la lista & se excluye la primera fila 
											# (el encabezado) y la última (es radioactiva)

# Se unen las tres listas en una definitiva: primeras diez capas + integración desde fuera + integración desde dentro
plot_final = ayuda1(R,L) + integracion_plot + ayuda2(centro_plot)

# Se redondean los valores de la lista definitiva: 'plot_final'
for i in range(len(plot_final)-1):
	i += 1
	for j in arange(3,8):
		plot_final[i][j] = round(plot_final[i][j], 5)

# Se muestra por pantalla la lista definitiva: 'plot_final'
# print(tabulate(plot_final, headers='firstrow', tablefmt='fancy_grid'))
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# CREACIÓN DEL GRÁFICO DE LAS MAGNITUDES EN FUNCIÓN DEL RADIO NORMALIZADO
#------------------------------------------------------------------------------------------------------------------------
# Índices del input de la función 'ayuda3':
# 0 = 'Generación de energía'; 1 = 'FASE'; 2 = número de capa; 3 = radio; 
# 4 = presión; 5 = temperatura; 6 = luminosidad; 7 = masa; 8 = n+1
def ayuda3(plot_final, pos):
	'''
	Generación de una lista con los valores en cada capa de la magnitud que se indique
	Input: lista definitiva de integración y la posición equivalente a la magnitud
	Output: lista de la magnitud deseada
	'''
	herramienta = []
	for i in range(len(plot_final)-1):
		herramienta.append(plot_final[i+1][pos])
	return herramienta

# Se almacenan en variables locales el output de la función 'ayuda3'
r_fin = ayuda3(plot_final, 3) 	# Radio 
P_fin = ayuda3(plot_final, 4)	# Presión 
T_fin = ayuda3(plot_final, 5)	# Temperatura
L_fin = ayuda3(plot_final, 6)	# Luminosidad
M_fin = ayuda3(plot_final, 7)	# Masa

# Se normalizan las magnitudes mediante funciones compactas
r_visualizar = [r / max(r_fin) for r in r_fin] 	# Radio normalizado
P_visualizar = [P / max(P_fin) for P in P_fin]	# Presión normalizada
T_visualizar = [T / max(T_fin) for T in T_fin]	# Temperatura normalizada
L_visualizar = [L / max(L_fin) for L in L_fin]	# Luminosidad normalizada
M_visualizar = [M / max(M_fin) for M in M_fin]	# Masa normalizada
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# GRÁFICO DE LAS MAGNITUDES EN FUNCIÓN DEL RADIO NORMALIZADO
#------------------------------------------------------------------------------------------------------------------------
# Se grafican las magnitudes integradas en fundion del radio normalizado
fig, ax = plt.subplots(figsize = (7,4))
ax.plot(r_visualizar, P_visualizar, linewidth = 3, label = "Presion")
ax.plot(r_visualizar, T_visualizar, linewidth = 3, label = "Temperatura")
ax.plot(r_visualizar, L_visualizar, linewidth = 3, label = "Luminosidad")
ax.plot(r_visualizar, M_visualizar, linewidth = 3, label = "Masa")
ax.set_title("Magnitudes a lo largo del radio", fontsize = 16)
ax.set_xlabel("Radio normalizado", fontsize = 11)
ax.set_ylabel("Magnitudes normalizadas", fontsize = 11)
ax.legend()
ax.grid(True)

plt.show()
#------------------------------------------------------------------------------------------------------------------------


