from C_Algoritmos import *
import matplotlib.pyplot as plt
from numpy import arange


R = 11.6263157895 
L = 45.2526315789 
Tc = 1.86036036036 


def ayuda1(R, L):
	'''
	Calcula las primeras 10 capas
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
	Cambia el índice 1 sobre la integración desde dentro
	'''

	for pos in range(len(lista)):
		lista[pos][2] = 101 - len(lista) + pos

	return lista


valores_desde_fuera = desde_fuera(R, L)

integracion = valores_desde_fuera[0]
tam = valores_desde_fuera[1]
K = valores_desde_fuera[2]

n_rad = integracion[tam][8]
n_conv = integracion[tam+1][8]

r_rad = integracion[tam][3]
r_conv = integracion[tam+1][3]
r_inter = interpolation(r_rad,r_conv,n_rad,n_conv)

centro = desde_dentro(R, Tc, K, r_inter)


integracion_plot = integracion[1:tam+1] # No cogemos el encabezado
centro_plot = centro[::-1][1:-1] # No cogemos ni el encabezado ni el1o de la capa convectiva porque es de la radiactiva. 
# Primero, le damos la vuelta a la lista

# juntamos las tres partes
plot_final = ayuda1(R,L) + integracion_plot + ayuda2(centro_plot)

# redondeamos lso valores de plot_final
for i in range(len(plot_final)-1):
	i += 1
	for j in arange(3,8):
		plot_final[i][j] = round(plot_final[i][j], 5)


print(tabulate(plot_final, headers='firstrow', tablefmt='fancy_grid'))


def ayuda3(plot_final, pos):
	'''
	Genera una lista de la magnitud que le indico con los valores en cada capa 
	'''
	herramienta = []
	for i in range(len(plot_final)-1):
		herramienta.append(plot_final[i+1][pos])

	return herramienta

r_fin = ayuda3(plot_final, 3)
P_fin = ayuda3(plot_final, 4)
T_fin = ayuda3(plot_final, 5)
L_fin = ayuda3(plot_final, 6)
M_fin = ayuda3(plot_final, 7)

'''
r_rad = ayuda(integracion, tam, 3)
P_rad = ayuda(integracion, tam, 4)
T_rad = ayuda(integracion, tam, 5)
L_rad = ayuda(integracion, tam, 6)
M_rad = ayuda(integracion, tam, 7)

tam_2 = len(centro)-2
r_conv = ayuda(centro, tam_2, 3)[::-1]
P_conv = ayuda(centro, tam_2, 4)[::-1]
T_conv = ayuda(centro, tam_2, 5)[::-1]
L_conv = ayuda(centro, tam_2, 6)[::-1]
M_conv = ayuda(centro, tam_2, 7)[::-1]

r_fin = r_rad + r_conv
P_fin = P_rad + P_conv
T_fin = T_rad + T_conv
L_fin = L_rad + L_conv
M_fin = M_rad + M_conv
'''

# Normalizar las magnitudes 
r_visualizar = [r / max(r_fin) for r in r_fin]
P_visualizar = [P / max(P_fin) for P in P_fin]
T_visualizar = [T / max(T_fin) for T in T_fin]
L_visualizar = [L / max(L_fin) for L in L_fin]
M_visualizar = [M / max(M_fin) for M in M_fin]

fig, ax = plt.subplots(figsize = (7,4))
ax.plot(r_visualizar, P_visualizar, linewidth = 3, label = "Presion")
ax.plot(r_visualizar, T_visualizar, linewidth = 3, label = "Temperatura")
ax.plot(r_visualizar, L_visualizar, linewidth = 3, label = "Luminosidad")
ax.plot(r_visualizar, M_visualizar, linewidth = 3, label = "Masa")

ax.set_title("Magnitudes a lo largo del radio")
ax.set_xlabel("Radio normalizado")
ax.set_ylabel("Magnitudes normalizado")
ax.legend()
ax.grid(True)

plt.show()



