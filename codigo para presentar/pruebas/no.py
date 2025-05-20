import numpy as np
import matplotlib.pyplot as plt

# --- CONSTANTES ---
GRAVEDAD_STD = 9.81  # m/s²

# --- FUNCIÓN PARA PEDIR DATOS ---
def pedir_datos():
    print("\n=== ENTRADA DE DATOS ===")
    datos = {
        'x0': float(input("Posición inicial x0 (m): ")),
        'y0': float(input("Posición inicial y0 (m): ")),
        'v0': float(input("Velocidad inicial (m/s): ")),
        'angulo': np.radians(float(input("Ángulo de lanzamiento (°): "))),
        'g': float(input(f"Gravedad (m/s²) [Enter para {GRAVEDAD_STD}]: ") or GRAVEDAD_STD),
        'resistencia': input("¿Considerar resistencia del aire? (s/n): ").lower() == 's'
    }
    
    if datos['resistencia']:
        datos.update({
            'masa': float(input("Masa del proyectil (kg): ")),
            'coef_arrastre': float(input("Coeficiente de arrastre: ")),
            'area': float(input("Área transversal (m²): ")),
            'densidad_aire': float(input("Densidad del aire (kg/m³) [1.225 por defecto]: ") or 1.225)
        })
    return datos

# --- FUNCIONES DE CÁLCULO ---
def descomposicion_vectorial(v0, theta):
    return v0 * np.cos(theta), v0 * np.sin(theta)

def tiempo_vuelo(v0y, g, y0=0):
    """Tiempo hasta que y(t) = 0."""
    return (v0y + np.sqrt(v0y**2 + 2 * g * y0)) / g

def alcance_maximo(v0x, t_vuelo):
    return v0x * t_vuelo

def altura_maxima(v0y, g):
    return (v0y**2) / (2 * g)

def movimiento_parabolico(x0, y0, v0x, v0y, g, t):
    x = x0 + v0x * t
    y = y0 + v0y * t - 0.5 * g * t**2
    vx = v0x
    vy = v0y - g * t
    return x, y, vx, vy

def movimiento_con_resistencia(x0, y0, v0x, v0y, g, masa, coef_arrastre, area, densidad_aire, t_max, dt=0.01):
    """Simulación numérica con resistencia del aire (Euler method)."""
    tiempos = np.arange(0, t_max, dt)
    x, y = [x0], [y0]
    vx, vy = v0x, v0y
    
    for t in tiempos[1:]:
        velocidad = np.sqrt(vx**2 + vy**2)
        fuerza_arrastre = 0.5 * densidad_aire * coef_arrastre * area * velocidad**2
        ax = - (fuerza_arrastre / masa) * (vx / velocidad)
        ay = -g - (fuerza_arrastre / masa) * (vy / velocidad)
        
        vx += ax * dt
        vy += ay * dt
        
        x.append(x[-1] + vx * dt)
        y.append(y[-1] + vy * dt)
        
        if y[-1] < 0:  # Detener si toca el suelo
            break
    
    return x, y, tiempos[:len(x)]

# --- GENERAR INFORME ---
def generar_informe(datos, resultados):
    informe = f"""
=== INFORME DE TIRO PARABÓLICO ===
- Posición inicial: ({datos['x0']:.2f}, {datos['y0']:.2f}) m
- Velocidad inicial: {datos['v0']:.2f} m/s
- Ángulo: {np.degrees(datos['angulo']):.2f}°
- Gravedad: {datos['g']:.2f} m/s²
- {'CON resistencia del aire' if datos['resistencia'] else 'SIN resistencia del aire'}

=== RESULTADOS ===
- Tiempo de vuelo: {resultados['t_vuelo']:.2f} s
- Alcance máximo: {resultados['alcance']:.2f} m
- Altura máxima: {resultados['altura_max']:.2f} m
- Velocidad final (impacto): {np.sqrt(resultados['vx_final']**2 + resultados['vy_final']**2):.2f} m/s
"""
    if datos['resistencia']:
        informe += f"- Coeficiente de arrastre: {datos['coef_arrastre']:.4f}\n"
    
    print(informe)
    return informe

# --- CONSULTA POR TIEMPO ESPECÍFICO ---
def consultar_tiempo(datos, t):
    if t < 0:
        print("¡Error: El tiempo no puede ser negativo!")
        return
    
    if datos['resistencia']:
        x, y, _ = movimiento_con_resistencia(
            datos['x0'], datos['y0'], resultados['v0x'], resultados['v0y'], 
            datos['g'], datos['masa'], datos['coef_arrastre'], datos['area'], 
            datos['densidad_aire'], t_max=t, dt=0.01
        )
        pos_x, pos_y = x[-1], y[-1]
        vx = (x[-1] - x[-2]) / 0.01 if len(x) > 1 else resultados['v0x']
        vy = (y[-1] - y[-2]) / 0.01 if len(y) > 1 else resultados['v0y']
    else:
        pos_x, pos_y, vx, vy = movimiento_parabolico(
            datos['x0'], datos['y0'], resultados['v0x'], resultados['v0y'], 
            datos['g'], t
        )
    
    # Verificar si el proyectil ya cayó
    if pos_y < 0:
        print(f"\nEn t = {t:.2f} s: El proyectil ya ha impactado en el suelo.")
    else:
        print(f"\nEn t = {t:.2f} s:")
        print(f"- Posición: ({pos_x:.2f}, {pos_y:.2f}) m")
        print(f"- Velocidad: ({vx:.2f}, {vy:.2f}) m/s")

# --- GRÁFICA ---
def graficar_trayectoria(x, y, titulo):
    plt.figure(figsize=(10, 5))
    plt.plot(x, y, 'b-', label="Trayectoria")
    plt.title(titulo)
    plt.xlabel("Distancia (m)")
    plt.ylabel("Altura (m)")
    plt.grid(True)
    plt.legend()
    plt.show()

# --- PROGRAMA PRINCIPAL ---
if __name__ == "__main__":
    datos = pedir_datos()
    resultados = {}
    
    # Descomposición vectorial
    resultados['v0x'], resultados['v0y'] = descomposicion_vectorial(datos['v0'], datos['angulo'])
    
    # Cálculos básicos
    resultados['t_vuelo'] = tiempo_vuelo(resultados['v0y'], datos['g'], datos['y0'])
    resultados['alcance'] = alcance_maximo(resultados['v0x'], resultados['t_vuelo'])
    resultados['altura_max'] = datos['y0'] + altura_maxima(resultados['v0y'], datos['g'])
    
    # Simulación
    if datos['resistencia']:
        x, y, tiempos = movimiento_con_resistencia(
            datos['x0'], datos['y0'], resultados['v0x'], resultados['v0y'], 
            datos['g'], datos['masa'], datos['coef_arrastre'], datos['area'], 
            datos['densidad_aire'], resultados['t_vuelo']
        )
        resultados['vx_final'] = (x[-1] - x[-2]) / 0.01 if len(x) > 1 else 0
        resultados['vy_final'] = (y[-1] - y[-2]) / 0.01 if len(y) > 1 else 0
    else:
        x, y, resultados['vx_final'], resultados['vy_final'] = movimiento_parabolico(
            datos['x0'], datos['y0'], resultados['v0x'], resultados['v0y'], 
            datos['g'], resultados['t_vuelo']
        )
        x = [datos['x0'] + resultados['v0x'] * t for t in np.linspace(0, resultados['t_vuelo'], 100)]
        y = [datos['y0'] + resultados['v0y'] * t - 0.5 * datos['g'] * t**2 for t in np.linspace(0, resultados['t_vuelo'], 100)]
    
    # Informe y gráfica
    informe = generar_informe(datos, resultados)
    graficar_trayectoria(x, y, "Trayectoria del Proyectil")
    
    # Consulta interactiva
    while True:
        opcion = input("\n¿Consultar posición en un tiempo específico? (s/n): ").lower()
        if opcion != 's':
            break
        t = float(input("Tiempo (s): "))
        consultar_tiempo(datos, t)