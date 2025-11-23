#!/bin/bash

# 1. Limpiar archivos anteriores
rm -f simulador simulacion.xyz

# 2. Compilar con optimización máxima (-O3)
echo "Compilando simulador..."
g++ main.cpp -o simulador -std=c++11 -O3

# 3. Verificar si la compilación fue exitosa
if [ $? -eq 0 ]; then
    echo "Ejecutando simulación..."
    ./simulador
    echo "¡Listo! Resultados guardados en 'simulacion.xyz'"
else
    echo "Error: Falló la compilación."
fi