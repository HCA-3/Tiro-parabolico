import pygame
import numpy as np
import math
import tkinter as tk
from tkinter import simpledialog, messagebox

# Inicialización de Pygame
pygame.init()

# Configuración de la pantalla
WIDTH, HEIGHT = 800, 600
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Simulación de Tiro Parabólico")

# Colores
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
BLUE = (0, 0, 255)
GREEN = (0, 255, 0)
GRAY = (200, 200, 200)

# Parámetros de simulación
g = 9.81  # Gravedad (m/s^2)
time_step = 0.05  # Paso de tiempo

# Variables para almacenar datos previos
previous_data = None

# Función para obtener valores desde una ventana de entrada
def get_input():
    global previous_data
    root = tk.Tk()
    root.withdraw()
    
    if previous_data:
        use_previous = messagebox.askyesno("Opciones", "¿Desea usar los datos anteriores?")
        if use_previous:
            return previous_data
    
    def submit():
        try:
            v0 = max(0, float(v0_entry.get()))
            angle = max(0, float(angle_entry.get()))
            wind_speed = float(wind_entry.get())
            mass = max(0.1, float(mass_entry.get()))
            drag_coefficient = max(0, float(drag_entry.get()))
            
            global previous_data
            previous_data = (v0, angle, wind_speed, mass, drag_coefficient)
            input_window.destroy()
        except ValueError:
            messagebox.showerror("Error", "Ingrese valores numéricos válidos.")
    
    input_window = tk.Toplevel()
    input_window.title("Entrada de Datos")
    
    tk.Label(input_window, text="Velocidad inicial (m/s):").pack()
    v0_entry = tk.Entry(input_window)
    v0_entry.pack()
    
    tk.Label(input_window, text="Ángulo de lanzamiento (grados):").pack()
    angle_entry = tk.Entry(input_window)
    angle_entry.pack()
    
    tk.Label(input_window, text="Velocidad del viento (m/s):").pack()
    wind_entry = tk.Entry(input_window)
    wind_entry.pack()
    
    tk.Label(input_window, text="Masa del proyectil (kg):").pack()
    mass_entry = tk.Entry(input_window)
    mass_entry.pack()
    
    tk.Label(input_window, text="Coeficiente de arrastre:").pack()
    drag_entry = tk.Entry(input_window)
    drag_entry.pack()
    
    tk.Button(input_window, text="Aceptar", command=submit).pack()
    input_window.wait_window()
    
    return previous_data

# Función para simular el movimiento parabólico
def simulate_projectile(v0, angle, wind_speed, mass, drag_coefficient):
    angle_rad = math.radians(angle)
    vx = v0 * math.cos(angle_rad) + wind_speed  # Componente x con viento
    vy = v0 * math.sin(angle_rad)
    x, y = 0, 0
    trajectory = []
    max_x, max_y = 0, 0
    time = 0
    while y >= 0:
        trajectory.append((x, y, vx, vy))
        x += vx * time_step
        vy -= g * time_step
        y += vy * time_step
        max_x = max(max_x, x)
        max_y = max(max_y, y)
        time += time_step
    return trajectory, max_x, max_y, time

# Función para dibujar la cuadrícula
def draw_grid():
    for x in range(50, WIDTH - 50, 50):
        pygame.draw.line(screen, GRAY, (x, 50), (x, HEIGHT - 50), 1)
    for y in range(50, HEIGHT - 50, 50):
        pygame.draw.line(screen, GRAY, (50, y), (WIDTH - 50, y), 1)

# Función para dibujar el plano cartesiano
def draw_cartesian_plane():
    draw_grid()
    pygame.draw.line(screen, BLACK, (50, HEIGHT - 50), (WIDTH - 50, HEIGHT - 50), 2)  # Eje X
    pygame.draw.line(screen, BLACK, (50, HEIGHT - 50), (50, 50), 2)  # Eje Y
    
    font = pygame.font.Font(None, 24)
    for i in range(0, WIDTH, 100):
        pygame.draw.line(screen, BLACK, (i + 50, HEIGHT - 55), (i + 50, HEIGHT - 45), 2)
        text = font.render(f"{(i//10)} m", True, BLACK)
        screen.blit(text, (i + 40, HEIGHT - 40))
    for i in range(0, HEIGHT, 100):
        pygame.draw.line(screen, BLACK, (45, HEIGHT - 50 - i), (55, HEIGHT - 50 - i), 2)
        text = font.render(f"{(i//10)} m", True, BLACK)
        screen.blit(text, (10, HEIGHT - 55 - i))

# Loop principal
def main():
    v0, angle, wind_speed, mass, drag_coefficient = get_input()
    if not v0 or not angle:
        return  # Evita que el programa continúe sin datos válidos
    trajectory, max_x, max_y, flight_time = simulate_projectile(v0, angle, wind_speed, mass, drag_coefficient)
    scale_x = (WIDTH - 100) / max_x if max_x > 0 else 1
    scale_y = (HEIGHT - 100) / max_y if max_y > 0 else 1
    running = True
    clock = pygame.time.Clock()
    
    index = 0  # Índice para la animación progresiva
    
    while running:
        screen.fill(WHITE)
        draw_cartesian_plane()
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
        
        if index < len(trajectory):
            index += 1  # Avanza en la trayectoria lentamente
        
        for i in range(index):  # Dibuja toda la trayectoria acumulada
            point = trajectory[i]
            scaled_x = int(50 + point[0] * scale_x)
            scaled_y = int(HEIGHT - 50 - point[1] * scale_y)
            pygame.draw.circle(screen, RED, (scaled_x, scaled_y), 3)
        
        pygame.display.flip()
        clock.tick(10)  # Reduce la velocidad de la animación
    
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                return

if __name__ == "__main__":
    main()
