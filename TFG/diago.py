import numpy as np

def obtener_minimo_autovalor(matriz):
    autovalores = np.linalg.eigvals(matriz)
    min_autovalor = np.min(autovalores)
    return min_autovalor

# Ejemplo de uso




def vector_a_matriz(nombre_archivo_entrada, nombre_archivo_salida):
    # Leer el vector desde el archivo de entrada
    with open(nombre_archivo_entrada, 'r') as archivo_entrada:
        vector = np.loadtxt(archivo_entrada)

    # Calcular la dimensión de la matriz
    n = int(np.sqrt(len(vector)))

    # Reshape del vector a matriz
    matriz = np.reshape(vector, (n, n))

    # Escribir la matriz en el archivo de salida
    with open(nombre_archivo_salida, 'w') as archivo_salida:
        np.savetxt(archivo_salida, matriz, fmt='%.28f')

# Ejemplo de uso
nombre_archivo_entrada = 'data.dat'
nombre_archivo_salida = 'data2.dat'

vector_a_matriz(nombre_archivo_entrada, nombre_archivo_salida)




matriz = np.loadtxt('data2.dat')
min_autovalor = obtener_minimo_autovalor(matriz)
print("El mínimo autovalor de la matriz es:", min_autovalor)

