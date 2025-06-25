def cargar_datos(path):
    f = open(path, 'r')
    lineas = f.readlines()
    f.close()

    tiempo = float(lineas[0].strip())
    valores = []
    for i in range(1, len(lineas)):
        linea = lineas[i].strip()
        valores.append(float(linea))

    cuerpos = []
    i = 0
    while i < len(valores):
        cuerpo = []
        cuerpo.append(valores[i])
        cuerpo.append(valores[i + 1])
        cuerpo.append(valores[i + 2])
        cuerpos.append(cuerpo)
        i += 3

    return tiempo, cuerpos

def error_absoluto(a, b):
    return abs(a - b)

def error_relativo(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    mayor = abs(a)
    if abs(b) > mayor:
        mayor = abs(b)
    if mayor == 0.0:
        return float('inf')
    return abs(a - b) / mayor

def contar_decimales_iguales(a, b):
    a_str = "%.15f" % a
    b_str = "%.15f" % b
    i = 0
    while i < len(a_str):
        if a_str[i] != b_str[i]:
            break
        i += 1
    punto = a_str.find('.')
    if i <= punto:
        return 0
    return i - punto - 1

def comparar_cuerpos(lista1, lista2):
    if len(lista1) != len(lista2):
        print("Diferente número de cuerpos")
        return [], [], []

    errores_abs = []
    errores_rel = []
    decimales = []

    i = 0
    while i < len(lista1):
        ea = [0.0, 0.0, 0.0]
        er = [0.0, 0.0, 0.0]
        dc = [0, 0, 0]

        j = 0
        while j < 3:
            x = lista1[i][j]
            y = lista2[i][j]
            ea[j] = error_absoluto(x, y)
            er[j] = error_relativo(x, y)
            dc[j] = contar_decimales_iguales(x, y)
            j += 1

        errores_abs.append(ea)
        errores_rel.append(er)
        decimales.append(dc)
        i += 1

    return errores_abs, errores_rel, decimales

def promedio(lista):
    total = [0.0, 0.0, 0.0]
    n = len(lista)

    i = 0
    while i < n:
        j = 0
        while j < 3:
            total[j] += lista[i][j]
            j += 1
        i += 1

    resultado = [0.0, 0.0, 0.0]
    j = 0
    while j < 3:
        resultado[j] = total[j] / n
        j += 1

    return resultado

def mayor_error(lista):
    max_error = 0.0
    indice = -1

    i = 0
    while i < len(lista):
        j = 0
        while j < 3:
            if lista[i][j] > max_error:
                max_error = lista[i][j]
                indice = i
            j += 1
        i += 1

    return indice, lista[indice]

output_secuencial = "secuencial.txt"
output_paralelo = "hibrido.txt"

tiempo1, cuerpos1 = cargar_datos(output_secuencial)
tiempo2, cuerpos2 = cargar_datos(output_paralelo)

errores_abs, errores_rel, decimales = comparar_cuerpos(cuerpos1, cuerpos2)

print("Tiempo secuencial: %.6f s" % tiempo1)
print("Tiempo paralelo:   %.6f s" % tiempo2)
print("")

prom_abs = promedio(errores_abs)
prom_rel = promedio(errores_rel)
prom_dec = promedio(decimales)

print("Error absoluto medio por coordenada (x, y, z):", prom_abs)
print("Error relativo medio por coordenada (x, y, z):", prom_rel)
print("Decimales promedio coincidentes por coordenada:", prom_dec)

indice, peor = mayor_error(errores_abs)
print("")
print("Cuerpo con mayor error absoluto: #", indice)
print("Error absoluto (x, y, z):", peor)
print("Posiciones secuencial:", cuerpos1[indice])
print("Posiciones paralelo:  ", cuerpos2[indice])
