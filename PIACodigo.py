import random
import time
import sys
import pandas as pd
import seaborn as se
from math import e

def ins_BuenosAires(lon, mod):
    instancia = {}
    for i in range((lon + 1)**2):
        instancia[str(i)] = []
        
    if lon % 2 != 0:
        print('La longuitud debe ser par')
        return instancia
    elif lon < 5:
        print('La longitud debe de ser mayor a 4')
        return instancia
    
    
        
    ## Las diagonales hacia el cuadrado
    if mod >= 2:
        for i in range(lon):
            instancia[str((lon + 2) * i)].append((lon + 2) * (i + 1))
            instancia[str((lon + 1) * (lon - i) + i)].append((lon + 1) * (lon - (i + 1)) + (i + 1))
            if mod > 3:
                instancia[str((lon + 2) * (i + 1))].append((lon + 2) * i)
                instancia[str((lon + 1) * (lon - (i + 1)) + (i + 1))].append((lon + 1) * (lon - i) + i)
    iz_der = True
    for i in range(lon + 1):
        for j in range(lon):
            if iz_der:
                ## Vertical hacia abajo
                instancia[str((lon + 1) * ( j + 1) + i)].append((lon + 1) * j + i)
            else:
                ## Vertical hacia arriba
                instancia[str((lon + 1) * j + i)].append((lon + 1) * (j + 1) + i)
        iz_der = not iz_der
        
    ## El rombo
    if mod >= 3:
        for i in range(int(lon / 2)):
            ## Parte superior izquierda del rombo
            instancia[str(int(lon / 2) + (lon * i))].append(int(lon / 2) + (lon * (i + 1)))
            ## Parte superior derecha del rombo
            instancia[str(((int(lon / 2) + 1) * (lon + 1) - 1)-(lon + 2) * i)].append(((int(lon / 2) + 1) * (lon + 1) - 1)-(lon + 2) * (i+1))
            ## Parte inferior izquierda del rombo 
            instancia[str(((int(lon / 2) + 1) * (lon + 1) - 1) + (lon * (i + 1)))].append(((int(lon / 2) + 1) * (lon + 1) - 1) + (lon * i))
            ## Parte inferior derecha
            instancia[str(int(lon / 2) * (lon + 1) + (lon + 2) * i)].append(int(lon / 2) * (lon + 1) + (lon + 2) * (i + 1))
            if mod > 3:
                instancia[str(int(lon / 2) + (lon * (i + 1)))].append(int(lon / 2) + (lon * i))
                instancia[str(((int(lon / 2) + 1) * (lon + 1) - 1)-(lon + 2) * (i+1))].append(((int(lon / 2) + 1) * (lon + 1) - 1)-(lon + 2) * i)
                instancia[str(((int(lon / 2) + 1) * (lon + 1) - 1) + (lon * i))].append(((int(lon / 2) + 1) * (lon + 1) - 1) + (lon * (i + 1)))
                instancia[str(int(lon / 2) * (lon + 1) + (lon + 2) * (i + 1))].append(int(lon / 2) * (lon + 1) + (lon + 2) * i)
                
  
    ## Sentido de las carreteras
    iz_der = True
    for i in range(lon + 1):
        for j in range(lon): 
            if iz_der:
                ## Horizontal a la izquierda
                instancia[str(lon + (lon + 1) * i - j)].append(lon + ((lon + 1) * i) - (j + 1))
            else:
                ## Horizontal a la derecha
                instancia[str((lon + 1) * i + j)].append((lon + 1) * i + (j + 1))
        iz_der = not iz_der
        
    if mod == 1:
        ## invertir el ultimo nodo
        del instancia[str(((lon+1)**2)-1)][1]
        instancia[str(((lon+1)**2)-2)].append(((lon+1)**2)-1)
        
        instancia['0'].append(1)
        del instancia['1'][1]
        
    return instancia

def dis_prob(lon, mod):
    lista = [1] * ((lon + 1) ** 2)
    prob = 1 / ((lon / 2) + 1) 
    if mod == 1 or mod == 3:
        for i in range(int(lon/2)):
            for j in range(i):
                ## Hacia arriba
                lista[(lon + 1) * (lon - j) + i] = prob / (4 * (2 * i + 1))
                lista[(lon + 1) * ((lon + 1) - j) - (i + 1)] = prob / (4 * (2 * i + 1))
                ## abajo
                lista[i + (lon + 1) * j] = prob / (4 * (2 * i + 1))
                lista[(lon - i) + (lon + 1) * j] = prob / (4 * (2 * i + 1))
                ## Hacia la izquierda o a la derecha
                lista[(lon + 1) * (lon - i) + j] = prob / (4 * (2 * i + 1))
                lista[(lon + 1)* ((lon + 1)- i) - (j + 1)] = prob / (4 * (2 * i + 1))
                lista[i * (lon + 1) + j] = prob / (4 * (2 * i + 1))
                lista[(lon - j) + (lon + 1) * i] = prob / (4 * (2 * i + 1))
            ## Diagonal arriba a la izquierda
            lista[(lon + 2) * i] = prob / (4 * (2 * i + 1))
            ## Diagonal arriaba a la derecha
            lista[lon * (i + 1)] = prob / (4 * (2 * i + 1))
            ## Diagonal abajo izquierda
            lista[(lon + 1) * (lon - i) + i] = prob / (4 * (2 * i + 1))
            ## Diagonal abajo a la derecha
            lista[(lon + 1) * ((lon + 1) - i) - (i + 1)] = prob / (4 * (2 * i + 1))
        
        for i in range(lon + 1):
            lista[int(lon / 2) * (lon + 1) + i] = prob / (2 * lon + 1)
            lista[int(lon / 2) + (lon + 1) * i] = prob / (2 * lon + 1)
            
        if mod == 3:
            random.shuffle(lista)
            
    if mod == 2:
        for i in range(int(lon / 2)):
            ## Diagonal arriba a la izquierda
            lista[(lon + 2) * i] = prob / (4 *(lon - (2 * i + 1)) + 4)
            ## Diagonal arriaba a la derecha
            lista[lon * (i + 1)] = prob / (4 *(lon - (2 * i + 1)) + 4)
            ## Diagonal abajo izquierda
            lista[(lon + 1) * (lon - i) + i] = prob / (4 *(lon - (2 * i + 1)) + 4)
            ## Diagonal abajo derecha
            lista[(lon + 1) * ((lon + 1) - i) - (i + 1)] = prob / (4 *(lon - (2 * i + 1)) + 4)
            
            for j in range(2 * (int(lon / 2) - (i + 1)) + 1): 
                ## Lado superior
                lista[(lon + 2) * i + (j + 1)] =  prob / (4 *(lon - (2 * i + 1)) + 4)
                ## Lado izquierdo
                lista[lon * (i + 1) + (lon + 1) * (j + 1)] = prob / (4 *(lon - (2 * i + 1)) + 4)
                ## Lado inferior
                lista[(lon + 1) * (lon - i) + i + (j + 1)] = prob / (4 *(lon - (2 * i + 1)) + 4)
                ## Lado derecho 
                lista[(lon + 2) * i + (lon + 1) * (j + 1)] = prob / (4 *(lon - (2 * i + 1)) + 4)
        ## Punto del medio
        lista[int((((lon+1)**2) - 1) / 2) ] = prob
    return lista

def cal_prob(lista):
    prob = random.random()
    cot_inf = 0.0
    cot_sup = lista[0]
    for i in range(len(lista)-1):
        if (cot_inf <= prob) and (prob <= cot_sup):
            return i
        cot_inf = cot_inf + lista[i]
        cot_sup = cot_sup + lista[i + 1]
    i = random.randint(0,len(lista))
    if i > len(lista)-1:
        i = i - 1
    return i

def bfs(grafo,nodo,objetivo):
    tiempo = time.perf_counter()
    if nodo == objetivo:
        return [nodo],0.0
    visitados = [False] * len(grafo)
    reg_caminos =[]
    for i in range(len(grafo)):
        reg_caminos.append([]) 
    cola = [nodo]
    while len(cola)>0:
        nodo = cola.pop(0)
        if not visitados[nodo]:
            visitados[nodo] = True
            reg_caminos[nodo].append(nodo)
            for llave in grafo[str(nodo)]:
                if not llave == objetivo:
                    if not visitados[llave]:
                        cola.append(llave)
                        reg_caminos[llave] = reg_caminos[nodo][:]
                else:
                    reg_caminos[nodo].append(llave)
                    return reg_caminos[nodo],  time.perf_counter()-tiempo
def datos(lon):
    lista_cost = []
    for i in range(((lon + 1)**2)):
        lista_cost.append([])
    for k in range(((lon + 1)**2)):
        reg = k // (lon + 1)
        col = k % (lon + 1) 
        for i in range(lon + 1):
            for j in range(lon + 1):
                lista_cost[reg * (lon + 1) + col].append((((i - reg)**2)+((j - col)**2))**0.5)
    return lista_cost
    
    
def A_est(grafo, nodo, objetivo, datos ):
    tiempo = time.perf_counter()
    if nodo == objetivo:
        return [nodo],0.0
    visitados = [False] * len(grafo)
    pasos = 0
    reg_caminos = []
    for i in range(len(grafo)):
        reg_caminos.append([])
    sucesores = [nodo]
    while len(sucesores) > 0:
        
        nodo = sucesores.pop(0)
        if not visitados[nodo]:
            visitados[nodo] = True
            reg_caminos[nodo].append(nodo)
            
        for llave in grafo[str(nodo)]:
            if not llave == objetivo:
                if not visitados[llave]:
                    sucesores.append(llave)
                    reg_caminos[llave] = reg_caminos[nodo][:]
            else:
                reg_caminos[nodo].append(llave)
                return reg_caminos[nodo],  time.perf_counter()-tiempo
        
        pasos = pasos + 1
        costos = [1] * len(sucesores)
        for i in range(len(sucesores)):
            costos[i] = pasos + datos[objetivo][sucesores[i]]
        sucesores = [x for (i, x) in sorted(zip(costos, sucesores))]

def dfs(grafo,nodo,objetivo):
    tiempo = time.perf_counter()
    if nodo == objetivo:
        return [nodo],0.0
    visitados = [False] * len(grafo)
    reg_caminos =[]
    for i in range(len(grafo)):
        reg_caminos.append([]) 
    pila = [nodo]
    while len(pila)>0:
        nodo = pila.pop()
        if not visitados[nodo]:
            visitados[nodo] = True
            reg_caminos[nodo].append(nodo)
            grafo[str(nodo)].reverse()
            for llave in grafo[str(nodo)]:
                if not llave == objetivo:
                    if not visitados[llave]:
                        pila.append(llave)
                        reg_caminos[llave] = reg_caminos[nodo][:]
                else:
                    reg_caminos[nodo].append(llave)
                    return reg_caminos[nodo],  time.perf_counter()-tiempo


def solucion(grafo,nodo, obj, lon):
    ## Identificaremos hacia donde tenemos que orientar el sentido del perimetro segun un sistema de cuadrantes donde si
    ## la diferencia entre las y del los nodos (que representara lo ubicacion en renglones de la matriz) es positiva eso quiere
    ## decir que el nodo objetivo se encuentra por debajo del nodo de salida, si es negativo se encontrara arriba
    sum_res = [0,0]
    y_del = (obj // (lon + 1) - nodo // (lon + 1))
    if y_del >= 0:
        sum_res[1] = 1 
    ## En caso del componente de la diferencia del componente de las x si es positiva eso quiere decir que se encuentra a la
    ## derecha, en el caso contrario se encontraria a la izquierda
    x_del = ( obj % (lon + 1) - nodo % (lon + 1))
    if x_del >= 0:
        sum_res[0] = 1

    reg_camino = [nodo]
    while not (obj in grafo[str(nodo)]): 
        pat_dis = [[],[]]
        for i in grafo[str(nodo)]:
            
            if i - nodo == (lon + 1) or i - nodo == -(lon + 1):  
                if (i - nodo > 0) == sum_res[1]:
                    pat_dis[0].append(1/len(grafo[str(nodo)]))
                    pat_dis[1].append(len(pat_dis[0]) - 1)
                else:
                    pat_dis[0].append(1/(10*len(grafo[str(nodo)])))
            elif i - nodo == 1 or i - nodo -1:
                if (i - nodo > 0) == sum_res[0]:
                    pat_dis[0].append(1/len(grafo[str(nodo)]))
                    pat_dis[1].append(len(pat_dis[0]) - 1)
                else:
                    pat_dis[0].append(1/(10*len(grafo[str(nodo)])))
            
            elif i- nodo != -7 or i - nodo != -1 or i- nodo != 7 or i - nodo != 1 :
                if  (obj // (lon + 1) - nodo // (lon + 1) > 0) == sum_res[1]:
                    if (obj // (lon + 1) - nodo // (lon + 1) > 0) == sum_res[0]:
                        pat_dis[0].append(1/(len(grafo[str(nodo)])))  
                        pat_dis[1].append(len(pat_dis[0]) - 1)
                    else:
                         pat_dis[0].append(1/(10*len(grafo[str(nodo)])))
                else:
                    pat_dis[0].append(1/(10*len(grafo[str(nodo)])))

        if len(pat_dis[1]) != len(pat_dis[0]):
            restos = (1-sum(pat_dis[0]))/(len(pat_dis[0])-len(pat_dis[1]))
            for i in pat_dis[1]:
                pat_dis[0][i] = pat_dis[0][i] + restos
        
        if len(pat_dis[1]) == 0:
            con = 0
            for i in pat_dis[0]:
                pat_dis[0][con] = i*10
                con = con + 1
        nodo = grafo[str(nodo)][cal_prob(pat_dis[0])]
        reg_camino.append(nodo)
    reg_camino.append(obj)
    return reg_camino

def recocido_sim(grafo, nodo, obj, lon, tem, enf):
    tiempo = time.perf_counter()
    recorrido_sal = solucion(grafo,nodo, obj, lon)    ## Realizamos una solucion inicial
    while tem>0:                                      ## Detenemos el ciclo hasta parar la temperatura
        recorrido_p = solucion(grafo,nodo, obj, lon)  ## Realizamos una solucion probicional 
        if len(recorrido_sal) > len(recorrido_p):     ## Es la solucion probicional mejor que la actual
            recorrido_sal = recorrido_p               ## Si, asignamos la solucion probicionala la actual
        else: 
            if e**((len(recorrido_sal) - len(recorrido_p))/tem) > .05:  ## No, decidimos si tomar la solucion con la prbabilidad 
            ## ^-------------------------------------------------- en base a esta funcion
                recorrido_sal= recorrido_p         ## Si toca entonces cambiamos la solucion
        tem = tem - enf                            ## Disminuimos la temperatira en cada iteracion
    return recorrido_sal, time.perf_counter() - tiempo
l, t = recocido_sim(ins_BuenosAires(6,4), 8 , 3, 6, 10, 1)

if __name__ == __main__:
    def problema_UberDidi(viajes, lon, mod_ciudad, mod_prob):
        ubicacion = random.randint(0,((lon + 1)**2 - 1))
        pasajero = 0
        destino = 0
        ins = ins_BuenosAires(lon,mod_ciudad)
        if mod_prob == 0:
            dis_pasajeros = dis_prob(lon, 1)
            dis_destinos = dis_prob(lon, 2)
        
        elif mod_prob == 1:
            dis_pasajeros = dis_prob(lon, 2)
            dis_destinos = dis_prob(lon, 1)
        
        elif mod_prob == 2:
            dis_pasajeros = dis_prob(lon, 3)
            dis_destinos = dis_prob(lon, 3)
        
        temp = sys.stdout
        sys.stdout = open('ReporteDeResultados_ProblemaDidi_UberBFS.csv','a')
        for i in range(viajes):
            pasajero = cal_prob(dis_pasajeros)
            destino = cal_prob(dis_destinos)
            lista, tiempo = bfs(ins,ubicacion,pasajero)
            print(ubicacion,',',pasajero,',',len(lista),',',tiempo,',',lon,',',mod_ciudad,',',mod_prob)
            lista, tiempo = bfs(ins,pasajero,destino)
            print(pasajero,',',destino,',',len(lista),',',tiempo,',',lon,',',mod_ciudad,',',mod_prob)
            ubicacion = destino
        sys.stdout = temp
        return ' '

    temp = sys.stdout
    sys.stdout = open('ReporteDeResultados_ProblemaDidi_UberBFS.csv','a')
    print('Salida,Destino,Longuitud de camino,Tiempo de ejecucion,Dimensiones de la ciudad,Instancia de la ciudad,Distribucion de probabilidad')
    for j in range(4):
        j = j + 1
        for k in range(3):
            for i in range(30):
                i = (2 * i) + 6
                problema_UberDidi(500,i,j,k)
    sys.stdout = temp

    res = pd.read_csv('ReporteDeResultados_ProblemaDidi_uberBFS.csv')
    res = res.drop(index = res[res['Tiempo de ejecucion'] == 0].index)
    res['Dimensiones de la ciudad'] = res['Dimensiones de la ciudad'] + 1
    res['Clasificaciones'] =  res['Instancia de la ciudad'] + res['Distribucion de probabilidad'] * 4
    
    temp = sys.stdout
    sys.stdout = open('ReporteDeResultados_ProblemaDidi_UberDFS.csv','a')
    print('Salida,Destino,Longuitud de camino,Tiempo de ejecucion,Dimensiones de la ciudad,Instancia de la ciudad,Distribucion de probabilidad')

    for (s,d,n,ins,prob) in zip(res['Salida'],res['Destino'],res['Dimensiones de la ciudad'],res['Instancia de la ciudad'],res['Distribucion de probabilidad']):
        l, t = dfs(ins_BuenosAires(n,ins) ,s ,d)
        print(s,',',d,',',len(l),',',t,',',n,',',ins,',',prob)
    sys.stdout = temp

    lista_datos = []
    for i in range(30):
        i = (2 * i) + 6
        lista_datos.append(datos(i))

    temp = sys.stdout
    sys.stdout = open('ReporteDeResultados_ProblemaDidi_Uber.csv','a')
    print('Salida,Destino,Longuitud de camino,Tiempo de ejecucion,Dimensiones de la ciudad,Instancia de la ciudad,Distribucion de probabilidad')
    for (s,d,n,ins,prob) in zip(res['Salida'],res['Destino'],res['Dimensiones de la ciudad'],res['Instancia de la ciudad'],res['Distribucion de probabilidad']):
        l, t = A_est(ins_BuenosAires(n,ins) ,s ,d,lista_datos[(n-6)//2])
        print(s,',',d,',',len(l),',',t,',',n,',',ins,',',prob) 
    sys.stdout = temp 

    resd = pd.read_csv('ReporteDeResultados_ProblemaDidi_uberDFS.csv')
    resd['Dimensiones de la ciudad'] = resd['Dimensiones de la ciudad'] + 1
    resd['Clasificaciones'] =  resd['Instancia de la ciudad'] + resd['Distribucion de probabilidad'] * 4

    resa = pd.read_csv('ReporteDeResultados_ProblemaDidi_uber.csv')
    resa = resa.drop(index = resa[resa['Tiempo de ejecucion'] == 0].index)
    resa['Dimensiones de la ciudad'] = resa['Dimensiones de la ciudad'] + 1
    resa['Clasificaciones'] =  resa['Instancia de la ciudad'] + resa['Distribucion de probabilidad'] * 4


    g = se.catplot(x = 'Dimensiones de la ciudad', y = 'Tiempo de ejecucion', data = res, kind = 'box', height = 7)
    g.set_axis_labels('Dimension de la ciudad','Tiempo')
    res['Tiempo de ejecucion'].mean()

    g = se.catplot(x = 'Dimensiones de la ciudad', y = 'Tiempo de ejecucion', data = resd, kind = 'box', height = 7)
    g.set_axis_labels('Dimension de la ciudad','Tiempo')
    resd['Tiempo de ejecucion'].mean()

    g = se.catplot(x = 'Dimensiones de la ciudad', y = 'Tiempo de ejecucion', data = resa, kind = 'box', height = 7)
    g.set_axis_labels('Dimension de la ciudad','Tiempo')
    resa['Tiempo de ejecucion'].mean()

    f = se.catplot(x = 'Clasificaciones', y = 'Longuitud de camino', data = res, kind = 'box', height = 7)
    f.set_axis_labels('Dimension de la ciudad','Longuitud de camino')    

    f = se.catplot(x = 'Clasificaciones', y = 'Longuitud de camino',data = resd, kind = 'box', height = 7)
    f.set_axis_labels('Dimension de la ciudad','Longuitud de camino')

    f = se.catplot(x = 'Clasificaciones', y = 'Longuitud de camino',data = resa, kind = 'box', height = 7)
    f.set_axis_labels('Dimension de la ciudad','Longuitud de camino')

