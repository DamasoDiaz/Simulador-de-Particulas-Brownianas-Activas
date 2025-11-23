# Simulador de Partículas Brownianas Activas en C++

**Por: Dámaso Díaz**

![Demo2](https://github.com/user-attachments/assets/9d180e44-643f-405e-b438-de4da6180074)

Este proyecto implementa un simulador en C++ para modelar la dinámica de **partículas autopropulsadas** (materia activa) en un espacio 2D. El sistema combina movimiento Browniano, autopropulsión constante, interacciones de volumen excluido y alineamiento social.

El simulador implementa **condiciones de borde periódicas** y utiliza un algoritmo de **Cell Lists (Listas de Celdas)** para optimizar la detección de colisiones, reduciendo la complejidad computacional.

## Modelo Físico

Las partículas obedecen una dinámica de Langevin sobreamortiguada que consta de dos partes:

### 1. Traslación
Cada partícula se mueve debido a tres factores: su motor interno ($V_0$), la repulsión con otras partículas y el ruido térmico.

$$\dot{\vec{r}}_i = V_0 \hat{u}_i(\theta_i) + \vec{F}_{col} + \sigma_v \vec{\xi}$$

### 2. Rotación 
El ángulo de dirección de la partícula cambia por ruido rotacional y una tendencia a alinearse con sus vecinos cercanos.

$$\dot{\theta}_i = \frac{\gamma}{N_{vec}} \sum_{j} \sin[m(\theta_j - \theta_i)] + \sigma_{\theta} \eta$$

Donde:
* **$V_0$**: Velocidad de autopropulsión constante.
* **$\vec{F}_{col}$**: Fuerza de esferas blandas (Ley de Hooke) para evitar superposiciones.
* **$\gamma$**: Fuerza de alineamiento social.
* **$m$**: Simetría del alineamiento ($m=1$ polar, $m=2$ nemático).

---

## Características Principales

* **Alta Eficiencia:** Implementación de cuadrícula de vecindad (*Cell Lists*) para simular grandes cantidades de partículas rápidamente.
* **Modelo de Esferas Blandas:** Las colisiones se manejan mediante fuerzas repulsivas proporcionales a la superposición entre partículas.
* **Condiciones de Borde Periódicas (PBC):** Las partículas que salen por un lado de la caja reingresan por el opuesto, simulando un sistema continuo.
* **Configuración Flexible:** Todos los parámetros se leen desde un archivo externo `parametros.par`.
* **Visualización:** Genera un archivo `simulacion.xyz` compatible con **OVITO**, que incluye marcadores en las esquinas para visualizar correctamente los límites de la caja.

---

## Compilación y Ejecución

El repositorio incluye un script `run.sh` que automatiza la compilación (usando optimización `-O3`) y la ejecución inmediata.

### Uso Rápido
Simplemente ejecuta el script en tu terminal:

```bash
chmod +x run.sh  # Solo la primera vez, para dar permisos
./run.sh
