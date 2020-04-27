
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
#                           				UNIVERSIDAD COMPLUTENSE DE MADRID
# 
#                               				CURSO ACADÉMICO 2019-2020
#
########################################################################################################################
########################################################################################################################

#------------------------------------------------------------------------------------------------------------------------
# PAQUETES NECESARIOS
#------------------------------------------------------------------------------------------------------------------------
from C_algoritmos import * 				# Se importan los algoritmos de C_Algoritmos
import matplotlib.pyplot as plt 		# Se importa la biblioteca de matplotlib.pyplot como 'plt'
from numpy import argmin, zeros, where 	# Se importan las funciones argmin, zeros y where
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# TEMPERATURA CENTRAL CON MÍNIMO ERROR RELATIVO
#------------------------------------------------------------------------------------------------------------------------
# Se estudia la temperatura central con el mínimo error relativo 
T_centro_lista = linspace(1.5, 2.0, 1000) 	# Valores de la temperatura central
error_relativo_lista = []					# Lista-almacén de valores de la temperatura central

# Se calcula el error relativo para cada temperatura con los valores de Radio-Luminosidad iniciales
for n in T_centro_lista:
	error_relativo_lista.append(error_total_funcion(R_tot_inicial, L_tot_inicial, n))


x_min = argmin(error_relativo_lista) 	# Se selecciona la posición del error mínimo
Tc = T_centro_lista[x_min] 				# Se selecciona la temperatura con el error mínimo

# Se imprimen la temperatura central y su error relativo mínimo
print("\n --- Valor de la temperatura central --- \n",  "\n     Temperatura centro", round(Tc, 5), 
	"\n     Error relativo    ", round(error_relativo_lista[x_min], 5), "\n")
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# PAR RADIO-LUMINOSIDAD CON MÍNIMO ERROR RELATIVO
# (La temperatura central que se utiliza es la de menor error relativo, calculada previamente)
#------------------------------------------------------------------------------------------------------------------------
# Se procede al cálculo de la matriz de errores mínimos para cada par Radio-Luminosidad
num = 20							# Número de divisiones
L_lista = linspace(45, 45.4, num) 	# Límite de la luminosidad
R_lista = linspace(11.5, 11.7, num)	# Límite del radio 

# Introduce el valor del tanto por ciento calculado
print("\n--- Tanto por ciento del modelo calculado ---\n") 		

# Calculo y almacenamiento del error relativo de cada par Radio-Luminosidad en una matriz
error_relativo_matriz = zeros([num, num])
for i in range(num):

	# Muestra en pantalla el tanto por ciento calculado
	if (100*(i+1)/num) % 20 == 0: 
		print("    {}%".format(round(100*(i+1)/num),0))

	for j in range(num):
		error_relativo_matriz[i,j] = error_total_funcion(R_lista[i], L_lista[j], Tc)

# Se localizan el par Radio-Luminosidad con menor error relativo
pos_min = where(error_relativo_matriz == error_relativo_matriz.min())
R = R_lista[pos_min[0][0]] 	# Radio total con menor error relativo
L = L_lista[pos_min[1][0]]	# Luminosidad total con menor error relativo

# Se imprimen los resultados en pantalla en una tabla
print("\n", "--- RESULTADOS ---", "\n") 
valores = [['Magnitud', 'Valor'], ['Radio total', R], ['Luminosidad total', L], 		# Lista de resultados
['Temperatura centro', Tc], ['Error relativo total', error_relativo_matriz.min()] ]
print(tabulate(valores, headers='firstrow', tablefmt='fancy_grid', numalign="center"))	
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# GRÁFICO DEL MÍNIMO ERROR RELATIVO PARA CADA PAR RADIO-LUMINOSIDAD 
#------------------------------------------------------------------------------------------------------------------------
# Gráfico - mapa de colores del mínimo error relativo para cada par Radio-Luminosidad
fig, ax = plt.subplots(figsize = (7,4)) 						# Creando la figura
color_map = plt.imshow(error_relativo_matriz, 					# Mapa de mínimos errores
	origin="lower", cmap="Spectral", interpolation="bilinear")
color_bar = fig.colorbar(color_map, ax = ax) 					# Barra de colores
color_bar.set_label('Error relativo', fontsize = 11)           	# Etiqueta de la barra de colores	
ax.set_xlabel("Luminosidad", fontsize = 11)						# Etiqueta del eje X
ax.set_ylabel("Radio", fontsize = 11)							# Etiqueta del eje Y
ax.set_title("Mínimo error relativo", fontsize = 16) 			# Título
texto = plt.text(10, 1, "\nError relativo = {:.4f} %\n".format( # Cuadro de texto 
	error_relativo_matriz.min()), fontsize=10,ha = "center", va = "center")
                   
plt.show() 
#------------------------------------------------------------------------------------------------------------------------
