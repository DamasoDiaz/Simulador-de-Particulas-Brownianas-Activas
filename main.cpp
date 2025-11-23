#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <iomanip> 

// Compilación: g++ main.cpp -o simulador -std=c++11 -O3
// Ejecución:   ./simulador

// --- Constantes Globales ---
double r = 0.1;         // Radio de las partículas
bool Cuadrilla = true;  // true = Iniciar en grilla ordenada

// --- Clase Partícula ---
class Particula {
    public:
    double radio = r;
    double x = 0.0, y = 0.0;       // Posición
    double fx = 0.0, fy = 0.0;     // Fuerzas acumuladas (Colisiones)
    double vx = 0.0, vy = 0.0;     // Velocidad
    double theta = 0.0;            // Ángulo de orientación
    double sum_theta = 0.0;        // Acumulador para alineamiento

    // Mueve la partícula y aplica condiciones de borde periódicas (PBC)
    void movimiento(double dx, double dy, double dt, double limite_borde) {
        x += dx;
        y += dy;
        
        double ancho_caja = 2.0 * limite_borde;
        
        // Si sale por un lado, entra por el opuesto
        while (x > limite_borde) x -= ancho_caja;
        while (x < -limite_borde) x += ancho_caja;
        while (y > limite_borde) y -= ancho_caja;
        while (y < -limite_borde) y += ancho_caja;
        
        vx = dx/dt;
        vy = dy/dt;
    }
};

// --- Física de Interacción ---
// Calcula fuerzas entre dos partículas (p1 y p2)
void calcular_fuerza(Particula &p1, Particula &p2, double k_elastica, double limite_borde, double r, double ancho_caja, double m=1.0) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;

    // Imagen Mínima: Buscar la distancia más corta considerando bordes
    if (dx > limite_borde) dx -= ancho_caja;
    if (dx < -limite_borde) dx += ancho_caja;
    if (dy > limite_borde) dy -= ancho_caja;
    if (dy < -limite_borde) dy += ancho_caja;

    double distancia = std::sqrt(dx * dx + dy * dy);
    
    // 1. Repulsión (Solo si se tocan)
    if (distancia < 2.0 * r) {
        double overlap = 2.0 * r - distancia;
        double fuerza_magnitud = k_elastica * overlap; // Ley de Hooke
        double fx = (fuerza_magnitud * dx) / distancia;
        double fy = (fuerza_magnitud * dy) / distancia;

        // Tercera ley de Newton (acción-reacción)
        p1.fx -= fx; p1.fy -= fy;
        p2.fx += fx; p2.fy += fy;
    }

    // 2. Alineamiento (Siempre que sean vecinas en la grilla)
    double theta_ij = sin(m * (p1.theta - p2.theta));
    p1.sum_theta -= theta_ij;
    p2.sum_theta += theta_ij;
}

int main() {
    // --- 1. Configuración ---
    int k = 0, n = 0;
    double sigma_v = 0.0, sigma_theta = 0.0, dt = 0.001, k_elastica = 0.0, limite_borde = 10.0;
    double Vo = 0.0, gamma = 0.0, m = 1.0; 
    std::string parametro;

    // Leer archivo de parámetros
    std::ifstream archivo_parametros("parametros.par");
    if (!archivo_parametros.is_open()) {
        std::cerr << "ERROR: No se encontro 'parametros.par'." << std::endl;
        return 1;
    }

    while (archivo_parametros >> parametro) {  
        // Ignorar comentarios que empiezan con #
        if (parametro[0] == '#') {
            archivo_parametros.ignore(10000, '\n'); 
            continue;
        }
        // Asignar variables
        if (parametro == "particulas") archivo_parametros >> k;
        else if (parametro == "frames") archivo_parametros >> n;
        else if (parametro == "sigma_v") archivo_parametros >> sigma_v;
        else if (parametro == "sigma_theta") archivo_parametros >> sigma_theta;
        else if (parametro == "dt") archivo_parametros >> dt;
        else if (parametro == "k_elastica") archivo_parametros >> k_elastica;
        else if (parametro == "limite_borde") archivo_parametros >> limite_borde;
        else if (parametro == "Vo") archivo_parametros >> Vo;
        else if (parametro == "gamma") archivo_parametros >> gamma;
        else if (parametro == "m") archivo_parametros >> m;
    }
    archivo_parametros.close();

    // Validar parámetros básicos
    if (k <= 0 || n <= 0) {
        std::cerr << "ERROR: Revise 'particulas' o 'frames' en parametros.par" << std::endl;
        return 1;
    }

    double sqrt_dt = std::sqrt(dt);
    double ancho_caja = 2.0 * limite_borde;

    // --- 2. Inicialización ---
    std::vector<Particula> particulas(k);

    // Generadores de azar
    std::random_device rd;
    std::mt19937 generador(rd());
    std::normal_distribution<> distribucion_v(0.0, sigma_v);         
    std::normal_distribution<> distribucion_theta(0.0, sigma_theta); 
    std::uniform_real_distribution<> distribucion_angulo_inicial(0.0, 2.0 * M_PI);

    // Ángulos aleatorios iniciales
    for (int i = 0; i < k; ++i) particulas[i].theta = distribucion_angulo_inicial(generador);

    // Posiciones iniciales (Grilla cuadrada)
    if (Cuadrilla) {
        int columnas = static_cast<int>(std::ceil(std::sqrt(static_cast<double>(k))));
        double delta = ancho_caja / columnas;
        int contador = 0;
        for (int i = 0; i < columnas; ++i) {
            for (int j = 0; j < columnas; ++j) {
                if (contador >= k) break; 
                particulas[contador].x = -limite_borde + delta * (i + 0.5);
                particulas[contador].y = -limite_borde + delta * (j + 0.5);
                contador++;
            }
        }
    }

    std::ofstream archivo_salida("simulacion.xyz");
    
    // --- 3. Sistema de Cell Lists (Optimización) ---
    int numero_celdas = static_cast<int>(std::floor(ancho_caja / (2.0*r)));
    if (numero_celdas < 1) numero_celdas = 1; 
    double L_celda = ancho_caja / numero_celdas; 
    int celdas_totales = numero_celdas * numero_celdas;
    std::vector< std::vector<int> > grilla(celdas_totales);

    std::cout << "Iniciando simulacion con " << k << " particulas..." << std::endl;

    // --- 4. Bucle Principal (Tiempo) ---
    for (int frame = 0; frame < n; ++frame) {
        // Encabezado XYZ (k particulas + 4 esquinas visuales)
        archivo_salida << (k + 4) << "\nFrame " << frame + 1 << "\n";

        // A. Resetear y llenar la grilla
        for (int i = 0; i < celdas_totales; ++i) grilla[i].clear();

        for (int i = 0; i < k; ++i) {
            int cx = static_cast<int>((particulas[i].x + limite_borde) / L_celda);
            int cy = static_cast<int>((particulas[i].y + limite_borde) / L_celda);
            
            // Asegurar índices dentro de rango
            if (cx >= numero_celdas) cx = numero_celdas - 1; if (cx < 0) cx = 0;
            if (cy >= numero_celdas) cy = numero_celdas - 1; if (cy < 0) cy = 0;

            grilla[cy * numero_celdas + cx].push_back(i);
        }

        // B. Calcular interacciones (Vecinos)
        for (int cy = 0; cy < numero_celdas; ++cy) {
            for (int cx = 0; cx < numero_celdas; ++cx) {
                int idx_actual = cy * numero_celdas + cx;
                std::vector<int>& vecinos = grilla[idx_actual];

                // 1. Misma celda
                for (size_t p1 = 0; p1 < vecinos.size(); ++p1) {
                    for (size_t p2 = p1 + 1; p2 < vecinos.size(); ++p2) {
                        calcular_fuerza(particulas[vecinos[p1]], particulas[vecinos[p2]], 
                                        k_elastica, limite_borde, r, ancho_caja, m);
                    }
                }

                // 2. Celdas vecinas (Derecha, Abajo, Diagonales)
                int desplazamientos[4][2] = {{1,0}, {0,1}, {1,1}, {-1,1}};
                for (auto& desp : desplazamientos) {
                    int nx = cx + desp[0];
                    int ny = cy + desp[1];

                    // PBC para celdas
                    if (nx >= numero_celdas) nx = 0; else if (nx < 0) nx = numero_celdas - 1;
                    if (ny >= numero_celdas) ny = 0; // ny solo crece en este loop

                    int idx_vecina = ny * numero_celdas + nx;
                    
                    // Comparar lista actual con lista vecina
                    std::vector<int>& particulas_vecinas = grilla[idx_vecina];
                    for (int i : vecinos) {
                        for (int j : particulas_vecinas) {
                            calcular_fuerza(particulas[i], particulas[j], 
                                            k_elastica, limite_borde, r, ancho_caja, m);
                        }
                    }
                }
            }
        }

        // C. Integración (Movimiento)
        for (int i = 0; i < k; ++i) {
            // Guardar datos
            archivo_salida << "1 " << particulas[i].x << " " << particulas[i].y << " "
                           << particulas[i].radio << " " << particulas[i].fx << " " << particulas[i].fy << "\n";

            // 1. Actualizar Ángulo (Ruido + Alineamiento)
            particulas[i].theta += distribucion_theta(generador) * sqrt_dt + (gamma * particulas[i].sum_theta * dt);

            // 2. Calcular Desplazamiento (Ruido + Fuerza + Nado)
            double dx = distribucion_v(generador) * sqrt_dt + 
                        (particulas[i].fx * dt) + 
                        (Vo * std::cos(particulas[i].theta) * dt);
            
            double dy = distribucion_v(generador) * sqrt_dt + 
                        (particulas[i].fy * dt) + 
                        (Vo * std::sin(particulas[i].theta) * dt);

            // 3. Mover
            particulas[i].movimiento(dx, dy, dt, limite_borde);
            
            // 4. Resetear acumuladores
            particulas[i].fx = 0.0;
            particulas[i].fy = 0.0;
            particulas[i].sum_theta = 0.0;
        }

        // D. Dibujar bordes (Visualización)
        double r_esq = 0.001; 
        archivo_salida << "2 " << -limite_borde << " " << -limite_borde << " " << r_esq << " 0 0\n";
        archivo_salida << "2 " << -limite_borde << " " << +limite_borde << " " << r_esq << " 0 0\n";
        archivo_salida << "2 " << +limite_borde << " " << -limite_borde << " " << r_esq << " 0 0\n";
        archivo_salida << "2 " << +limite_borde << " " << +limite_borde << " " << r_esq << " 0 0\n";
    }

    std::cout << "Simulacion terminada exitosamente." << std::endl;
    return 0;
}