
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
from E_magnitudes_estudio import * 		# Se importa el estudio completo 
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
# GRÁFICO DE LA VARIACION DE LAS MAGNITUDES Y EL ERROR EN CADA UNA DE LAS ITERACIONES
#------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------
# Valores de las magnitudes y errores relativas a cada una de las siete iteraciones
#-------------------------------------------------------------------------------------------
# Radio total
total_radio = [12, 11.35897435897436, 11.430769230769231, 11.502564102564103, 
11.538461538461538, 11.574358974358974, 11.574358974358974, 11.602564102564102]

# Luminosidad total
total_luminosidad = [40, 42.717948717948715, 43.38461538461539, 43.8974358974359, 
44.333333333333336, 44.58974358974359, 44.84615384615385, 44.98461538461538]

# Temperatura central
centro_temperatura = [1.5, 1.8508508508508508, 1.8533533533533533, 1.8553553553553552, 
1.8568568568568569, 1.8578578578578577, 1.8588588588588588, 1.8593593593593594]

# Error temperatura central
error_centro_temperatura = [19.241971065636857, 4.79205108412574, 3.613497795049584, 
2.6452965387842884, 1.9401136362716924, 1.4590071718887634, 1.161278501762381]

# Error relativo total
total_error_relativo = [5.525840074514378, 4.190215495204249, 3.137465377639427, 
2.3192895968709157, 1.813116116966384, 1.264149149902196, 0.9652464210041786]

# Número de iteraciones
iteraciones = [0, 1, 2, 3, 4, 5, 6, 7]			# Se añade el 0 relativo al valor inicial
iteraciones_error = [1, 2, 3, 4, 5, 6, 7]		# Número de iteraciones

#-------------------------------------------------------------------------------------------
# Construcción de los gráficos
#-------------------------------------------------------------------------------------------
# Gráfico de la variación del radio total
fig, ax = plt.subplots(figsize = (7,4))
ax.plot(iteraciones, total_radio, linewidth = 2, label = "R", color = "blue")
ax.plot(0,12,'co')
ax.plot(7,11.602564102564102,'ro')
ax.set_title("Variación radio total en 7 iteraciones", fontsize = 12)
ax.set_xlabel("Iteraciones", fontsize = 9)
ax.set_ylabel("Radio total (cm)", fontsize = 9)
ax.legend()
ax.grid(True)

# Gráfico de la variación de la luminosidad total
fig, ax = plt.subplots(figsize = (7,4))
ax.plot(iteraciones, total_luminosidad, linewidth = 2, label = "L", color = "green")
ax.plot(0,40,'co')
ax.plot(7,44.98461538461538,'ro')
ax.set_title("Variación luminosidad total en 7 iteraciones", fontsize = 12)
ax.set_xlabel("Iteraciones", fontsize = 9)
ax.set_ylabel("Luminosidad total (erg s-1)", fontsize = 9)
ax.legend()
ax.grid(True)

# Gráfico de la variación de la temperatura central
fig, ax = plt.subplots(figsize = (7,4))
ax.plot(iteraciones, centro_temperatura, linewidth = 2, label = "Tc", color = "orange")
ax.plot(0,1.5,'co')
ax.plot(7,1.8593593593593594,'ro')
ax.set_title("Variación temperatura del centro en 7 iteraciones", fontsize = 12)
ax.set_xlabel("Iteraciones", fontsize = 9)
ax.set_ylabel("Temperatura centro (K)", fontsize = 9)
ax.legend()
ax.grid(True)

# Gráfico de la variación del error relativo total y del error de la temperatura central
fig, ax = plt.subplots(figsize = (7,4))
ax.plot(iteraciones_error, error_centro_temperatura, linewidth = 2, label = "Error Tc", color = "orange")
ax.plot(iteraciones_error, total_error_relativo, linewidth = 2, label = "Error total", color = "black")
ax.plot(1,5.525840074514378,'co')
ax.plot(7,0.9652464210041786,'ro')
plt.axhline(y=1, color='r', linestyle='-')
ax.set_title("Variación error relativo total", fontsize = 12)
ax.set_xlabel("Iteraciones", fontsize = 9)
ax.set_ylabel("Error relativo total (%)", fontsize = 9)
ax.legend()
ax.grid(True)

plt.show()
#------------------------------------------------------------------------------------------------------------------------