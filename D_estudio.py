from C_Algoritmos import * 
import matplotlib.pyplot as plt
from numpy import argmin, zeros, where



T_centro_lista = linspace(1.5, 2.0, 1000)
error_relativo_lista = []

for n in T_centro_lista:
	error_relativo_lista.append(error_total_funcion(R_tot_inicial, L_tot_inicial, n))




x_min = argmin(error_relativo_lista)
# print(x_min, T_centro_lista[x_min], error_relativo_lista[x_min])

Tc = T_centro_lista[x_min]

#plt.plot(T_centro_lista, error_relativo_lista)
#plt.plot(T_centro_lista[x_min], error_relativo_lista[x_min],'ro')
#plt.show()


num = 80
L_lista = linspace(44, 45.3, num)
R_lista = linspace(11.5, 11.7, num)
error_relativo_matriz = zeros([num, num])

avanze = 20	
print("\nTanto por ciento calculado:\n")
for i in range(num):

	if i == 8 or i == 16 or i == 24 or i == 32 or i == 40:
		print(avanze ,"%")
		avanze += 20

	for j in range(num):
		error_relativo_matriz[i,j] = error_total_funcion(R_lista[i], L_lista[j], Tc)

pos_min = where(error_relativo_matriz == error_relativo_matriz.min())
R = R_lista[pos_min[0][0]] 
L = L_lista[pos_min[1][0]]

print("\n", "  Radio total", R, "\n")
print("\n", "  Luminosidad total", L, "\n")
print("\n", "  Temperatura centro", Tc, "\n")
print("\n", "  Error relativo total", error_relativo_matriz.min(), "\n")


plt.imshow(error_relativo_matriz, interpolation="bilinear")
plt.show()

'''
variable = error_total_funcion(R, L, Tc)
integracion_final = variable[1]
centro_final = variable[2]

presion = integracion_final[:][4]
luminosidad = integracion_final[:][5]
temperatura = integracion_final[:][6]
masa = integracion_final[:][7] 
'''