import numpy as np
import iapws as ip
import matplotlib.pyplot as plt

#Seleção dos casos
Caso = int(input('Insira o caso (1 ou 2): '))
if Caso == 1:
    fator_potencia = 1
    fator_vazao = 1
if Caso == 2:
    tempo = int(input('Insira o tempo (32 ou 100) em segundos: '))
    if tempo == 32:
        fator_potencia = 0.05
        fator_vazao = 0.5273
    if tempo == 100:
        fator_potencia = 0.0389 
        fator_vazao = 0.1353 
            

##Dados

T_entrada = 30 + 273 #em Kelvin
P_reator = 5000000*fator_potencia #em W
Vazao_vol_EC = 22.8*fator_vazao #m³/h
FR = 1.914 #maior do esquema fornecido
N_EC = 24 #dado
N_controles = 1 #dado
N_placas_combustiveis = 18 #dado
N_placas_controles = 0 #dado
L_placa = 62.6 #mm
H_placa = 600 #mm
L_canal = 67.1 #mm
E_combustivel = 0.76 #mm
E_canal_interno = 2.89 #mm
E_revestimento = 0.38 #mm
k_revestimento = 180 #W/mK
P_entrada = 0.160 #MPa
P_saida = 0.150 #P_entrada-Perda de carga de 10kPa



##Cálculo Placas

N_total_placa = (N_EC*N_placas_combustiveis) + (N_placas_controles*N_controles)
A_total_troca = N_total_placa*H_placa*(L_placa*10**-6)*2
A_troca_1_placa = A_total_troca/N_total_placa
q_2linha_medio = P_reator/A_total_troca
D_h = 4*(E_canal_interno*L_canal)/10**3/(2*E_canal_interno + 2*L_canal)
A_escoamento_1_canal = ((E_canal_interno*L_canal)*10**-6)
v_escoamento = Vazao_vol_EC/((N_placas_combustiveis - 1)*3600*A_escoamento_1_canal)
V_total_combustivel = N_total_placa*E_combustivel*H_placa*(L_placa*10**-9)
q_3linha = P_reator/V_total_combustivel




##Propriedades da Água

def mi1_H2O(T_entrada, P_reator):
    H2O = ip.IAPWS97(T = T_entrada, P = P_reator)
    G = H2O.mu
    return G

def cp1_H2O(T_entrada, P_reator):
    H2O = ip.IAPWS97(T = T_entrada, P = P_reator)
    G = H2O.cp*1000
    return G

def rho1_H2O(T_entrada, P_reator):
    H2O = ip.IAPWS97(T = T_entrada, P = P_reator)
    G = H2O.rho
    return G

def k1_H2O(T_entrada, P_reator):
    H2O = ip.IAPWS97(T = T_entrada, P = P_reator)
    G = H2O.k
    return G

def h1(k_H2O, Re, Pr):
    G = 0.023*(k_H2O/D_h)*(Re**0.8)*(Pr**0.4)
    return G

def Pr(mi_H2O, cp_H2O, k_H2O):
    G = mi_H2O*cp_H2O/k_H2O
    return G

def Re(v_escoamento,rho_H2O, mi_H2O, D_h):
    G = v_escoamento*rho_H2O*D_h/mi_H2O
    return G

def k1_combustivel(T_entrada):
    T = (9/5*T_entrada) + 32
    k = 3978.1/(692.61 + T) + (6.02366*10**-12)*((T + 460)**3)
    G = 1.73*k
    return G

def q_1linha(Y):
    G = FR*q_2linha_medio*np.cos((np.pi/2)*(Y/(H_placa/2)))
    return G

def kelvin_to_celsius(temps):
    return [temp - 273.15 for temp in temps]


##Criando a malha

X_pontos = int(input('Insira o Tamanho da Malha: '))
Ponto_inicial = - (H_placa/2) + (H_placa/X_pontos)/2
Ponto_final = Ponto_inicial + (X_pontos - 1)*(H_placa/X_pontos)
Pontos_Y = np.linspace(Ponto_inicial, Ponto_final, X_pontos)
print()
print('-------------------------------------------------')
print()

Lista_Pontos_Y = []
for i in range(len(Pontos_Y)):
    Lista_Pontos_Y.append(Pontos_Y[i])

for i in range(len(Lista_Pontos_Y)):
    if i == 0:
        print('Centros utilizados:', Lista_Pontos_Y[i], 'mm')
    else:
        print('                   ', Lista_Pontos_Y[i], 'mm')
print()
print('-------------------------------------------------')
print()




##Cálculo do fluxo de calor por unidade de área

q_2linha_especifico = np.array([q_1linha(Pontos_Y[i]) for i in range(X_pontos)])

Lista_q_2linha = []
for i in range(len(q_2linha_especifico)):
    Lista_q_2linha.append(q_2linha_especifico[i])
for i in range(len(q_2linha_especifico)):
    if i == 0:
        print('q" =', Lista_q_2linha[i], 'W/m**2')
    else:
        print('    ', Lista_q_2linha[i], 'W/m**2')
print()
print('-------------------------------------------------')
print()
print('q_medio =', q_2linha_medio, 'W/m**2')
print()
print('-------------------------------------------------')
print()
print('A_total_troca =', A_total_troca, 'm**2')
print()
print('-------------------------------------------------')
print()




##Obtenção dos perfis de temperaturas (T4, T3, T2, T1)

Vazao_massica_canal = Vazao_vol_EC*rho1_H2O(T_entrada, P_entrada)/(3600*(N_placas_combustiveis - 1))
Gradiente_de_pressao = np.linspace(P_entrada, P_saida, X_pontos)
A_unidade_de_controle = A_troca_1_placa/X_pontos

#Solução para o Fluído (T4)
T4 = []
for i in range(X_pontos):
    T4.append(0)
T4[0] = T_entrada + q_2linha_especifico[0]*A_unidade_de_controle/(Vazao_massica_canal*cp1_H2O(T_entrada, Gradiente_de_pressao[0]))
for i in range(1, X_pontos):
    cp_H2O = cp1_H2O(T4[i - 1], Gradiente_de_pressao[i])
    T4[i] = T4[i - 1] + q_2linha_especifico[i]*A_unidade_de_controle/(Vazao_massica_canal*cp_H2O)
T4C = kelvin_to_celsius(T4)
for i in range(len(T4)):	
    if i == 0:
        print('T4 =', T4C[i], '°C')
    else:
        print('    ', T4C[i], '°C')
print()
print('-------------------------------------------------')
print()


#Solução para a interface Fluído/Revestimento (T3)
T3 = []
H = []
for i in range(X_pontos):
    T3.append(0)
    H.append(0)
mi_H2O = mi1_H2O(T4[0], Gradiente_de_pressao[0])
rho_H2O = rho1_H2O(T4[0], Gradiente_de_pressao[0])
cp_H2O = cp1_H2O(T4[0], Gradiente_de_pressao[0])
k_H2O = k1_H2O(T4[0], Gradiente_de_pressao[0])
v_escoamento = (Vazao_massica_canal/rho_H2O)/A_escoamento_1_canal
re = Re(v_escoamento, rho_H2O, mi_H2O, D_h)
pr = Pr(mi_H2O, cp_H2O, k_H2O)
h = h1(k_H2O, re, pr)
H[0] = h
T3[0] = T4[0] + q_2linha_especifico[0]/h
for i in range(1, X_pontos):
    mi_H2O = mi1_H2O(T3[i - 1], Gradiente_de_pressao[i - 1])
    rho_H2O = rho1_H2O(T3[i - 1], Gradiente_de_pressao[i - 1])
    cp_H2O = cp1_H2O(T3[i - 1], Gradiente_de_pressao[i - 1])
    k_H2O = k1_H2O(T3[i - 1], Gradiente_de_pressao[i - 1])
    v_escoamento = (Vazao_massica_canal/rho_H2O)/A_escoamento_1_canal
    re = Re(v_escoamento, rho_H2O, mi_H2O, D_h)
    pr = Pr(mi_H2O, cp_H2O, k_H2O)
    h = h1(k_H2O, re, pr)
    T3[i] = T4[i] + q_2linha_especifico[i]/h
    H[i] = h
T3C = kelvin_to_celsius(T3)
for i in range(len(T3)):	
    if i == 0:
        print('T3 =', T3C[i], '°C')
    else:
        print('    ', T3C[i], '°C')
print()
print('-------------------------------------------------')
print()
for i in range(len(T3)):	
    if i == 0:
        print('h =', H[i], 'W/m²·K')
    else:
        print('    ', H[i], 'W/m²·K')
print()
print('-------------------------------------------------')
print()

#Solução para o Revestimento (T2)
T2 = []
for i in range(X_pontos):
    T2.append(0)
for i in range(X_pontos):
    T2[i] = T3[i] + q_2linha_especifico[i]*(E_revestimento/1000/k_revestimento)
T2C = kelvin_to_celsius(T2)
for i in range(len(T2)):	
    if i == 0:
        print('T2 =', T2C[i], '°C')
    else:
        print('    ', T2C[i], '°C')
print()
print('-------------------------------------------------')
print()


#Solução para o Combustível (T1)
T1 = []
for i in range(X_pontos):
    T1.append(0)
T1[0] = T2[0] + q_2linha_especifico[0]*(E_combustivel/(2*1000))/(2*k1_combustivel(T2[0]))
for i in range(1, X_pontos):
    k_combustivel = k1_combustivel(T1[i - 1])
    T1[i] = T2[i] + q_2linha_especifico[i]*(E_combustivel/(2*1000))/(2*k_combustivel)
T1C = kelvin_to_celsius(T1)
for i in range(len(T1)):	
    if i == 0:
        print('T1 =', T1C[i], '°C')
    else:
        print('    ', T1C[i], '°C')



##Gráfico dos perfis de temperaturas

plt.plot(Pontos_Y, T1C, "--o", label="T1")
plt.plot(Pontos_Y, T2C, "--o", label="T2")
plt.plot(Pontos_Y, T3C, "--o", label="T3")
plt.plot(Pontos_Y, T4C, "--o", label="T4")

plt.legend()

plt.title('Distribuição de Temperaturas para ' + str(X_pontos) + ' pontos')
plt.xlabel('Posição X (mm)')
plt.ylabel('Temperatura (°C)')
plt.show()

##Gráfico do coeficiente de troca de calor por convecção

plt.plot(Pontos_Y, H, "--o", label="h")
plt.legend()
plt.title('Coeficiente de Transferência de Calor')
plt.xlabel('Posição X (mm)')
plt.ylabel('h (W/m²·K)')
plt.show()