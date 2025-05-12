import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import json
from PIL import Image, ImageTk
import matplotlib.animation as animation
from matplotlib.widgets import Slider
import os
import time

class TiroParabolicoApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Simulador de Tiro Parabólico Avanzado")
        self.root.geometry("1100x750")
        self.style = ttk.Style()
        self.style.configure('TNotebook.Tab', font=('Helvetica', 10, 'bold'))
        self.style.configure('Accent.TButton', font=('Helvetica', 10, 'bold'), foreground='blue')
        
        # Variables de control
        self.datos = {
            'x0': tk.DoubleVar(value=0),
            'y0': tk.DoubleVar(value=0),
            'v0': tk.DoubleVar(value=10),
            'angulo': tk.DoubleVar(value=45),
            'g': tk.DoubleVar(value=9.81),
            'resistencia': tk.BooleanVar(value=False),
            'masa': tk.DoubleVar(value=1),
            'coef_arrastre': tk.DoubleVar(value=0.47),
            'area': tk.DoubleVar(value=0.1),
            'densidad_aire': tk.DoubleVar(value=1.225),
            'unidades': tk.StringVar(value="métrico"),
            'comparar': tk.BooleanVar(value=False)
        }
        
        # Variables para el segundo tiro (comparación)
        self.datos_comparacion = {
            'x0': tk.DoubleVar(value=0),
            'y0': tk.DoubleVar(value=0),
            'v0': tk.DoubleVar(value=15),
            'angulo': tk.DoubleVar(value=30),
            'g': tk.DoubleVar(value=9.81),
            'resistencia': tk.BooleanVar(value=False),
            'masa': tk.DoubleVar(value=1),
            'coef_arrastre': tk.DoubleVar(value=0.47),
            'area': tk.DoubleVar(value=0.1),
            'densidad_aire': tk.DoubleVar(value=1.225)
        }
        
        self.resultados = {}
        self.resultados_comparacion = {}
        self.figura = None
        self.canvas = None
        self.animacion = None
        self.historial = []
        self.historial_file = "historial_simulaciones.json"
        
        # Cargar historial si existe
        self.cargar_historial()
        
        self.crear_interfaz()
    
    def crear_interfaz(self):
        # Notebook (pestañas)
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True)
        
        # Pestaña de Entrada de Datos
        frame_entrada = ttk.Frame(notebook)
        notebook.add(frame_entrada, text="Datos de Entrada")
        
        # Panel de parámetros básicos
        self.crear_panel_basico(frame_entrada)
        
        # Panel de resistencia del aire
        self.crear_panel_resistencia(frame_entrada)
        
        # Panel de comparación
        self.crear_panel_comparacion(frame_entrada)
        
        # Panel de unidades
        self.crear_panel_unidades(frame_entrada)
        
        # Pestaña de Resultados
        frame_resultados = ttk.Frame(notebook)
        notebook.add(frame_resultados, text="Resultados")
        
        # Panel de resultados
        self.crear_panel_resultados(frame_resultados)
        
        # Pestaña de Gráficos
        frame_graficos = ttk.Frame(notebook)
        notebook.add(frame_graficos, text="Gráficos")
        
        # Panel de gráficos
        self.crear_panel_graficos(frame_graficos)
        
        # Pestaña de Historial
        frame_historial = ttk.Frame(notebook)
        notebook.add(frame_historial, text="Historial")
        
        # Panel de historial
        self.crear_panel_historial(frame_historial)
        
        # Barra de herramientas
        self.crear_barra_herramientas()
    
    def crear_panel_basico(self, parent):
        frame = ttk.LabelFrame(parent, text="Parámetros Básicos", padding=10)
        frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        ttk.Label(frame, text="Posición inicial x0:").grid(row=0, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.datos['x0']).grid(row=0, column=1)
        
        ttk.Label(frame, text="Posición inicial y0:").grid(row=1, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.datos['y0']).grid(row=1, column=1)
        
        ttk.Label(frame, text="Velocidad inicial:").grid(row=2, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.datos['v0']).grid(row=2, column=1)
        
        ttk.Label(frame, text="Ángulo de lanzamiento:").grid(row=3, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.datos['angulo']).grid(row=3, column=1)
        
        ttk.Label(frame, text="Gravedad:").grid(row=4, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.datos['g']).grid(row=4, column=1)
        
        ttk.Checkbutton(frame, text="Considerar resistencia del aire", 
                       variable=self.datos['resistencia'], 
                       command=self.toggle_resistencia).grid(row=5, column=0, columnspan=2, pady=5)
    
    def crear_panel_resistencia(self, parent):
        self.frame_resistencia = ttk.LabelFrame(parent, text="Parámetros de Resistencia del Aire", padding=10)
        self.frame_resistencia.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
        
        ttk.Label(self.frame_resistencia, text="Masa del proyectil:").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.frame_resistencia, textvariable=self.datos['masa']).grid(row=0, column=1)
        
        ttk.Label(self.frame_resistencia, text="Coeficiente de arrastre:").grid(row=1, column=0, sticky="w")
        ttk.Entry(self.frame_resistencia, textvariable=self.datos['coef_arrastre']).grid(row=1, column=1)
        
        ttk.Label(self.frame_resistencia, text="Área transversal:").grid(row=2, column=0, sticky="w")
        ttk.Entry(self.frame_resistencia, textvariable=self.datos['area']).grid(row=2, column=1)
        
        ttk.Label(self.frame_resistencia, text="Densidad del aire:").grid(row=3, column=0, sticky="w")
        ttk.Entry(self.frame_resistencia, textvariable=self.datos['densidad_aire']).grid(row=3, column=1)
        
        # Inicialmente oculto
        self.toggle_resistencia()
    
    def crear_panel_comparacion(self, parent):
        self.frame_comparacion = ttk.LabelFrame(parent, text="Comparación de Tiros", padding=10)
        self.frame_comparacion.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")
        
        ttk.Checkbutton(self.frame_comparacion, text="Activar comparación", 
                       variable=self.datos['comparar'], 
                       command=self.toggle_comparacion).grid(row=0, column=0, columnspan=2, pady=5)
        
        # Frame para los parámetros de comparación (inicialmente oculto)
        self.frame_comparacion_params = ttk.Frame(self.frame_comparacion)
        self.frame_comparacion_params.grid(row=1, column=0, columnspan=2, sticky="nsew")
        
        ttk.Label(self.frame_comparacion_params, text="Posición inicial x0:").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.frame_comparacion_params, textvariable=self.datos_comparacion['x0']).grid(row=0, column=1)
        
        ttk.Label(self.frame_comparacion_params, text="Posición inicial y0:").grid(row=1, column=0, sticky="w")
        ttk.Entry(self.frame_comparacion_params, textvariable=self.datos_comparacion['y0']).grid(row=1, column=1)
        
        ttk.Label(self.frame_comparacion_params, text="Velocidad inicial:").grid(row=2, column=0, sticky="w")
        ttk.Entry(self.frame_comparacion_params, textvariable=self.datos_comparacion['v0']).grid(row=2, column=1)
        
        ttk.Label(self.frame_comparacion_params, text="Ángulo de lanzamiento:").grid(row=3, column=0, sticky="w")
        ttk.Entry(self.frame_comparacion_params, textvariable=self.datos_comparacion['angulo']).grid(row=3, column=1)
        
        ttk.Checkbutton(self.frame_comparacion_params, text="Resistencia del aire", 
                       variable=self.datos_comparacion['resistencia']).grid(row=4, column=0, columnspan=2, pady=5)
        
        # Inicialmente oculto
        self.toggle_comparacion()
    
    def crear_panel_unidades(self, parent):
        frame = ttk.LabelFrame(parent, text="Configuración de Unidades", padding=10)
        frame.grid(row=3, column=0, padx=10, pady=10, sticky="nsew")
        
        ttk.Radiobutton(frame, text="Sistema Métrico (m, kg)", 
                       variable=self.datos['unidades'], value="métrico",
                       command=self.actualizar_unidades).grid(row=0, column=0, sticky="w")
        
        ttk.Radiobutton(frame, text="Sistema Imperial (pies, lb)", 
                       variable=self.datos['unidades'], value="imperial",
                       command=self.actualizar_unidades).grid(row=1, column=0, sticky="w")
    
    def crear_panel_resultados(self, parent):
        # Frame principal con scroll
        canvas = tk.Canvas(parent)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Contenido de resultados
        self.resultados_text = tk.Text(scrollable_frame, wrap=tk.WORD, height=20, width=80)
        self.resultados_text.pack(padx=10, pady=10, fill="both", expand=True)
        
        # Frame de consulta de tiempo
        frame_consulta = ttk.LabelFrame(scrollable_frame, text="Consulta por Tiempo", padding=10)
        frame_consulta.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(frame_consulta, text="Tiempo:").pack(side="left")
        self.tiempo_consulta = ttk.Entry(frame_consulta, width=10)
        self.tiempo_consulta.pack(side="left", padx=5)
        
        ttk.Button(frame_consulta, text="Consultar", command=self.consultar_tiempo).pack(side="left")
        
        self.consulta_resultado = ttk.Label(scrollable_frame, text="", wraplength=500)
        self.consulta_resultado.pack(padx=10, pady=5)
        
        # Frame de energía
        frame_energia = ttk.LabelFrame(scrollable_frame, text="Energías", padding=10)
        frame_energia.pack(fill="x", padx=10, pady=5)
        
        self.energia_text = tk.Text(frame_energia, wrap=tk.WORD, height=5, width=80)
        self.energia_text.pack(fill="both", expand=True)
    
    def crear_panel_graficos(self, parent):
        self.figura = plt.Figure(figsize=(8, 4), dpi=100)
        self.ax = self.figura.add_subplot(111)
        
        # Frame principal para el gráfico y controles
        frame_principal = ttk.Frame(parent)
        frame_principal.pack(fill="both", expand=True)
        
        # Canvas para el gráfico
        self.canvas = FigureCanvasTkAgg(self.figura, master=frame_principal)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        # Frame para controles del gráfico
        frame_controles = ttk.Frame(frame_principal)
        frame_controles.pack(fill="x", padx=10, pady=5)
        
        # Frame para la información dinámica del punto
        self.frame_info = ttk.LabelFrame(frame_principal, text="Información del Punto", padding=10)
        self.frame_info.pack(fill="x", padx=10, pady=5)
        
        # Etiquetas para mostrar la información
        self.lbl_posicion = ttk.Label(self.frame_info, text="Posición: (0.00, 0.00) m | Tiempo: 0.00 s")
        self.lbl_posicion.pack(anchor="w")
        
        self.lbl_velocidad = ttk.Label(self.frame_info, text="Velocidad: (0.00, 0.00) m/s | Magnitud: 0.00 m/s")
        self.lbl_velocidad.pack(anchor="w")
        
        self.lbl_altura_distancia = ttk.Label(self.frame_info, text="Altura: 0.00 m | Distancia: 0.00 m")
        self.lbl_altura_distancia.pack(anchor="w")
        
        self.lbl_angulo = ttk.Label(self.frame_info, text="Ángulo: 0.00°")
        self.lbl_angulo.pack(anchor="w")
        
        # Slider para navegar por la trayectoria
        self.slider_trayectoria_frame = ttk.Frame(frame_controles)
        self.slider_trayectoria_frame.pack(fill="x", pady=5)
        
        ttk.Label(self.slider_trayectoria_frame, text="Posición en trayectoria:").pack(side="left")
        self.slider_trayectoria = ttk.Scale(
            self.slider_trayectoria_frame,
            from_=0,
            to=100,
            orient="horizontal",
            command=self.actualizar_punto_trayectoria
        )
        self.slider_trayectoria.pack(side="left", expand=True, fill="x", padx=5)
        
        # Botones de control del gráfico
        ttk.Button(frame_controles, text="Exportar Gráfico", command=self.exportar_grafico).pack(side="left", padx=5)
        ttk.Button(frame_controles, text="Animación", command=self.toggle_animacion).pack(side="left", padx=5)
        
        # Slider para velocidad de animación
        self.slider_velocidad_frame = ttk.Frame(frame_controles)
        self.slider_velocidad_frame.pack(side="left", padx=10)
        
        ttk.Label(self.slider_velocidad_frame, text="Velocidad:").pack(side="left")
        self.velocidad_animacion = tk.IntVar(value=50)
        ttk.Scale(
            self.slider_velocidad_frame,
            from_=1,
            to=100,
            variable=self.velocidad_animacion,
            orient="horizontal",
            length=100
        ).pack(side="left")
    
    def actualizar_punto_trayectoria(self, valor):
        if not hasattr(self, 'resultados') or not self.resultados:
            return
        
        # Obtener el índice basado en la posición del slider
        idx = int(float(valor) / 100 * (len(self.resultados['x']) - 1))
        idx = max(0, min(idx, len(self.resultados['x']) - 1))
        
        # Obtener los datos del punto
        x = self.resultados['x'][idx]
        y = self.resultados['y'][idx]
        t = self.resultados['tiempos'][idx]
        
        # Calcular velocidad (aproximación por diferencia finita)
        if idx == 0:
            vx = self.resultados['v0x']
            vy = self.resultados['v0y']
        else:
            dt = self.resultados['tiempos'][idx] - self.resultados['tiempos'][idx-1]
            vx = (self.resultados['x'][idx] - self.resultados['x'][idx-1]) / dt
            vy = (self.resultados['y'][idx] - self.resultados['y'][idx-1]) / dt
        
        velocidad_magnitud = np.sqrt(vx**2 + vy**2)
        angulo = np.degrees(np.arctan2(vy, vx))
        
        # Actualizar etiquetas
        self.lbl_posicion.config(text=f"Posición: ({x:.2f}, {y:.2f}) m | Tiempo: {t:.2f} s")
        self.lbl_velocidad.config(text=f"Velocidad: ({vx:.2f}, {vy:.2f}) m/s | Magnitud: {velocidad_magnitud:.2f} m/s")
        self.lbl_altura_distancia.config(text=f"Altura: {y:.2f} m | Distancia: {x:.2f} m")
        self.lbl_angulo.config(text=f"Ángulo: {angulo:.2f}°")
        
        # Actualizar punto en el gráfico
        if not hasattr(self, 'punto_interactivo'):
            self.punto_interactivo, = self.ax.plot([x], [y], 'ro', markersize=8)
        else:
            self.punto_interactivo.set_data([x], [y])
        
        self.canvas.draw()
    
    def crear_panel_historial(self, parent):
        # Frame principal con scroll
        canvas = tk.Canvas(parent)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Contenido del historial
        self.historial_listbox = tk.Listbox(scrollable_frame, width=100, height=20)
        self.historial_listbox.pack(padx=10, pady=10, fill="both", expand=True)
        
        # Botones de historial
        frame_botones = ttk.Frame(scrollable_frame)
        frame_botones.pack(fill="x", padx=10, pady=5)
        
        ttk.Button(frame_botones, text="Cargar Simulación", command=self.cargar_desde_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Eliminar Selección", command=self.eliminar_del_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Limpiar Historial", command=self.limpiar_historial).pack(side="right", padx=5)
        
        # Actualizar lista de historial
        self.actualizar_lista_historial()
    
    def crear_barra_herramientas(self):
        toolbar = ttk.Frame(self.root)
        toolbar.pack(fill="x", padx=5, pady=5)
        
        ttk.Button(toolbar, text="Calcular", command=self.calcular_trayectoria, style='Accent.TButton').pack(side="left", padx=5)
        ttk.Button(toolbar, text="Exportar Informe", command=self.exportar_informe).pack(side="left", padx=5)
        ttk.Button(toolbar, text="Limpiar", command=self.limpiar).pack(side="left", padx=5)
        ttk.Button(toolbar, text="Modo Avanzado", command=self.modo_avanzado).pack(side="left", padx=5)
        ttk.Button(toolbar, text="Salir", command=self.root.quit).pack(side="right", padx=5)
    
    def toggle_resistencia(self):
        if self.datos['resistencia'].get():
            self.frame_resistencia.grid()
        else:
            self.frame_resistencia.grid_remove()
    
    def toggle_comparacion(self):
        if self.datos['comparar'].get():
            self.frame_comparacion_params.grid()
        else:
            self.frame_comparacion_params.grid_remove()
    
    def toggle_animacion(self):
        if not hasattr(self, 'animacion') or self.animacion is None:
            self.iniciar_animacion()
        else:
            self.detener_animacion()
    
    def actualizar_unidades(self):
        sistema = self.datos['unidades'].get()
        if sistema == "métrico":
            # Convertir de imperial a métrico si es necesario
            pass
        else:
            # Convertir de métrico a imperial si es necesario
            pass
    
    def validar_datos(self, datos=None, es_comparacion=False):
        if datos is None:
            datos = {k: v.get() for k, v in self.datos.items() if k not in ['unidades', 'comparar']}
        
        try:
            # Validar que los valores sean positivos donde corresponda
            if datos['v0'] <= 0:
                raise ValueError("La velocidad inicial debe ser positiva")
                
            if not (0 <= datos['angulo'] <= 90):
                raise ValueError("El ángulo debe estar entre 0° y 90°")
                
            if datos.get('resistencia', False):
                if datos['masa'] <= 0:
                    raise ValueError("La masa debe ser positiva")
                if datos['coef_arrastre'] <= 0:
                    raise ValueError("El coeficiente de arrastre debe ser positivo")
                if datos['area'] <= 0:
                    raise ValueError("El área transversal debe ser positiva")
            
            return True
        except ValueError as e:
            mensaje = f"Error en {'comparación' if es_comparacion else 'datos principales'}: {str(e)}"
            messagebox.showerror("Error de validación", mensaje)
            return False
    
    def calcular_trayectoria(self):
        # Validar datos principales
        if not self.validar_datos():
            return
        
        # Validar datos de comparación si está activado
        if self.datos['comparar'].get():
            datos_comparacion = {k: v.get() for k, v in self.datos_comparacion.items()}
            if not self.validar_datos(datos_comparacion, es_comparacion=True):
                return
        
        # Obtener datos principales
        datos = {k: v.get() for k, v in self.datos.items() if k not in ['unidades', 'comparar']}
        datos['angulo'] = np.radians(datos['angulo'])
        
        # Calcular trayectoria principal
        self.resultados = self.calcular_trayectoria_individual(datos)
        
        # Calcular trayectoria de comparación si está activado
        if self.datos['comparar'].get():
            datos_comparacion = {k: v.get() for k, v in self.datos_comparacion.items()}
            datos_comparacion['angulo'] = np.radians(datos_comparacion['angulo'])
            self.resultados_comparacion = self.calcular_trayectoria_individual(datos_comparacion)
        
        # Actualizar interfaz
        self.actualizar_resultados(datos)
        self.actualizar_grafico(datos)
        
        # Guardar en historial
        self.guardar_en_historial(datos)
    
    def calcular_trayectoria_individual(self, datos):
        resultados = {}
        
        # Descomposición vectorial
        v0x, v0y = descomposicion_vectorial(datos['v0'], datos['angulo'])
        resultados['v0x'], resultados['v0y'] = v0x, v0y
        
        # Cálculos básicos
        resultados['t_vuelo'] = tiempo_vuelo(v0y, datos['g'], datos['y0'])
        resultados['alcance'] = alcance_maximo(v0x, resultados['t_vuelo'])
        resultados['altura_max'] = datos['y0'] + altura_maxima(v0y, datos['g'])
        
        # Simulación
        if datos['resistencia']:
            x, y, tiempos = movimiento_con_resistencia(
                datos['x0'], datos['y0'], v0x, v0y, 
                datos['g'], datos['masa'], datos['coef_arrastre'], datos['area'], 
                datos['densidad_aire'], resultados['t_vuelo']
            )
            resultados['vx_final'] = (x[-1] - x[-2]) / 0.01 if len(x) > 1 else 0
            resultados['vy_final'] = (y[-1] - y[-2]) / 0.01 if len(y) > 1 else 0
        else:
            x, y, resultados['vx_final'], resultados['vy_final'] = movimiento_parabolico(
                datos['x0'], datos['y0'], v0x, v0y, 
                datos['g'], resultados['t_vuelo']
            )
            x = [datos['x0'] + v0x * t for t in np.linspace(0, resultados['t_vuelo'], 100)]
            y = [datos['y0'] + v0y * t - 0.5 * datos['g'] * t**2 for t in np.linspace(0, resultados['t_vuelo'], 100)]
        
        resultados['x'] = x
        resultados['y'] = y
        resultados['tiempos'] = np.linspace(0, resultados['t_vuelo'], len(x))
        
        # Calcular energías
        resultados['energias'] = self.calcular_energias(datos, resultados)
        
        return resultados
    
    def calcular_energias(self, datos, resultados):
        energias = {
            'cinetica': [],
            'potencial': [],
            'total': []
        }
        
        g = datos['g']
        masa = datos['masa'] if datos['resistencia'] else 1  # Masa ficticia si no hay resistencia
        
        for i in range(len(resultados['x'])):
            if i == 0:
                vx = resultados['v0x']
                vy = resultados['v0y']
            else:
                dt = resultados['tiempos'][i] - resultados['tiempos'][i-1]
                vx = (resultados['x'][i] - resultados['x'][i-1]) / dt
                vy = (resultados['y'][i] - resultados['y'][i-1]) / dt
            
            velocidad = np.sqrt(vx**2 + vy**2)
            altura = resultados['y'][i]
            
            ec = 0.5 * masa * velocidad**2
            ep = masa * g * altura
            et = ec + ep
            
            energias['cinetica'].append(ec)
            energias['potencial'].append(ep)
            energias['total'].append(et)
        
        return energias
    
    def actualizar_resultados(self, datos):
        informe = f"""
=== INFORME DE TIRO PARABÓLICO ===
- Posición inicial: ({datos['x0']:.2f}, {datos['y0']:.2f}) m
- Velocidad inicial: {datos['v0']:.2f} m/s
- Ángulo: {np.degrees(datos['angulo']):.2f}°
- Gravedad: {datos['g']:.2f} m/s²
- {'CON resistencia del aire' if datos['resistencia'] else 'SIN resistencia del aire'}

=== RESULTADOS ===
- Tiempo de vuelo: {self.resultados['t_vuelo']:.2f} s
- Alcance máximo: {self.resultados['alcance']:.2f} m
- Altura máxima: {self.resultados['altura_max']:.2f} m
- Velocidad final (impacto): {np.sqrt(self.resultados['vx_final']**2 + self.resultados['vy_final']**2):.2f} m/s
"""
        if datos['resistencia']:
            informe += f"- Coeficiente de arrastre: {datos['coef_arrastre']:.4f}\n"
        
        if self.datos['comparar'].get():
            informe += "\n=== COMPARACIÓN ===\n"
            informe += f"- Tiempo de vuelo: {self.resultados_comparacion['t_vuelo']:.2f} s (vs {self.resultados['t_vuelo']:.2f} s)\n"
            informe += f"- Alcance máximo: {self.resultados_comparacion['alcance']:.2f} m (vs {self.resultados['alcance']:.2f} m)\n"
            informe += f"- Altura máxima: {self.resultados_comparacion['altura_max']:.2f} m (vs {self.resultados['altura_max']:.2f} m)\n"
        
        self.resultados_text.delete(1.0, tk.END)
        self.resultados_text.insert(tk.END, informe)
        
        # Actualizar energías
        self.actualizar_energias(datos)
    
    def actualizar_energias(self, datos):
        if not self.resultados or 'energias' not in self.resultados:
            return
            
        energias = self.resultados['energias']
        texto_energia = "=== ENERGÍAS ===\n"
        texto_energia += f"- Energía cinética inicial: {energias['cinetica'][0]:.2f} J\n"
        texto_energia += f"- Energía potencial inicial: {energias['potencial'][0]:.2f} J\n"
        texto_energia += f"- Energía total inicial: {energias['total'][0]:.2f} J\n\n"
        
        texto_energia += f"- Energía cinética final: {energias['cinetica'][-1]:.2f} J\n"
        texto_energia += f"- Energía potencial final: {energias['potencial'][-1]:.2f} J\n"
        texto_energia += f"- Energía total final: {energias['total'][-1]:.2f} J\n"
        
        if datos['resistencia']:
            perdida = (energias['total'][0] - energias['total'][-1]) / energias['total'][0] * 100
            texto_energia += f"\n- Pérdida de energía por resistencia: {perdida:.2f}%"
        
        self.energia_text.delete(1.0, tk.END)
        self.energia_text.insert(tk.END, texto_energia)
    
    def actualizar_grafico(self, datos):
        if self.animacion:
            self.detener_animacion()
        
        self.ax.clear()
        
        # Graficar trayectoria principal
        self.ax.plot(self.resultados['x'], self.resultados['y'], 'b-', label="Trayectoria 1")
        self.ax.plot(datos['x0'], datos['y0'], 'bo', label="Punto de lanzamiento 1")
        self.ax.plot(self.resultados['alcance'], 0, 'bx', label="Punto de impacto 1")
        
        # Graficar trayectoria de comparación si existe
        if self.datos['comparar'].get() and self.resultados_comparacion:
            self.ax.plot(self.resultados_comparacion['x'], self.resultados_comparacion['y'], 'r-', label="Trayectoria 2")
            self.ax.plot(self.datos_comparacion['x0'].get(), self.datos_comparacion['y0'].get(), 'ro', label="Punto de lanzamiento 2")
            self.ax.plot(self.resultados_comparacion['alcance'], 0, 'rx', label="Punto de impacto 2")
        
        self.ax.set_title("Trayectoria del Proyectil")
        self.ax.set_xlabel("Distancia (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.grid(True)
        self.ax.legend()
        
        # Configurar el slider para la nueva trayectoria
        self.slider_trayectoria.set(0)
        self.actualizar_punto_trayectoria(0)
        
        self.canvas.draw()
    
    def iniciar_animacion(self):
        if not self.resultados:
            messagebox.showwarning("Advertencia", "No hay datos para animar. Calcule primero una trayectoria.")
            return
            
        # Configurar animación
        self.linea_animada, = self.ax.plot([], [], 'go', markersize=8)
        self.tiempo_texto = self.ax.text(0.02, 0.95, '', transform=self.ax.transAxes)
        
        # Función de inicialización
        def init():
            self.linea_animada.set_data([], [])
            self.tiempo_texto.set_text('')
            return self.linea_animada, self.tiempo_texto
        
        # Función de animación
        def animate(i):
            if i >= len(self.resultados['x']):
                i = len(self.resultados['x']) - 1
                
            self.linea_animada.set_data(self.resultados['x'][i], self.resultados['y'][i])
            self.tiempo_texto.set_text(f'Tiempo: {self.resultados["tiempos"][i]:.2f} s')
            
            # Actualizar slider y etiquetas
            self.slider_trayectoria.set(i / len(self.resultados['x']) * 100)
            self.actualizar_punto_trayectoria(i / len(self.resultados['x']) * 100)
            
            return self.linea_animada, self.tiempo_texto
        
        # Calcular intervalo basado en la velocidad del slider
        intervalo = 110 - self.velocidad_animacion.get()  # Invertir el valor (1-100 a 110-10)
        
        # Crear animación
        self.animacion = animation.FuncAnimation(
            self.figura, animate, frames=len(self.resultados['x'])*2,  # Multiplicar por 2 para más suavidad
            init_func=init, blit=True, interval=intervalo, repeat=True
        )
        
        self.canvas.draw()
    
    def detener_animacion(self):
        if hasattr(self, 'animacion') and self.animacion:
            self.animacion.event_source.stop()
            self.animacion = None
            
            # Limpiar elementos de animación
            if hasattr(self, 'linea_animada'):
                self.linea_animada.remove()
            if hasattr(self, 'tiempo_texto'):
                self.tiempo_texto.remove()
            
            self.canvas.draw()
    
    def consultar_tiempo(self):
        try:
            t = float(self.tiempo_consulta.get())
            if t < 0:
                raise ValueError("El tiempo no puede ser negativo")
                
            datos = {k: v.get() for k, v in self.datos.items() if k not in ['unidades', 'comparar']}
            datos['angulo'] = np.radians(datos['angulo'])
            
            if datos['resistencia']:
                x, y, _ = movimiento_con_resistencia(
                    datos['x0'], datos['y0'], self.resultados['v0x'], self.resultados['v0y'], 
                    datos['g'], datos['masa'], datos['coef_arrastre'], datos['area'], 
                    datos['densidad_aire'], t_max=t, dt=0.01
                )
                pos_x, pos_y = x[-1], y[-1]
                vx = (x[-1] - x[-2]) / 0.01 if len(x) > 1 else self.resultados['v0x']
                vy = (y[-1] - y[-2]) / 0.01 if len(y) > 1 else self.resultados['v0y']
            else:
                pos_x, pos_y, vx, vy = movimiento_parabolico(
                    datos['x0'], datos['y0'], self.resultados['v0x'], self.resultados['v0y'], 
                    datos['g'], t
                )
            
            if pos_y < 0:
                resultado = f"En t = {t:.2f} s: El proyectil ya ha impactado en el suelo."
            else:
                resultado = (
                    f"En t = {t:.2f} s:\n"
                    f"- Posición: ({pos_x:.2f}, {pos_y:.2f}) m\n"
                    f"- Velocidad: ({vx:.2f}, {vy:.2f}) m/s"
                )
            
            self.consulta_resultado.config(text=resultado)
            
        except ValueError as e:
            messagebox.showerror("Error", str(e))
    
    def exportar_informe(self):
        if not self.resultados:
            messagebox.showwarning("Advertencia", "No hay resultados para exportar")
            return
            
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Archivos de texto", "*.txt"), ("Archivos PDF", "*.pdf"), ("Todos los archivos", "*.*")],
            title="Guardar informe como"
        )
        
        if filepath:
            try:
                contenido = self.resultados_text.get(1.0, tk.END) + "\n" + self.energia_text.get(1.0, tk.END)
                
                if filepath.endswith('.pdf'):
                    # Exportar a PDF (requiere fpdf2)
                    try:
                        from fpdf import FPDF
                        pdf = FPDF()
                        pdf.add_page()
                        pdf.set_font("Arial", size=12)
                        
                        # Dividir el texto en líneas que quepan en la página
                        for line in contenido.split('\n'):
                            pdf.cell(200, 10, txt=line, ln=1)
                        
                        pdf.output(filepath)
                    except ImportError:
                        messagebox.showerror("Error", "Para exportar a PDF, instale fpdf2 con: pip install fpdf2")
                        return
                else:
                    # Exportar a texto plano
                    with open(filepath, 'w') as f:
                        f.write(contenido)
                
                messagebox.showinfo("Éxito", f"Informe guardado en:\n{filepath}")
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo guardar el archivo:\n{str(e)}")
    
    def exportar_grafico(self):
        if not self.figura:
            messagebox.showwarning("Advertencia", "No hay gráfico para exportar")
            return
            
        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("JPEG", "*.jpg"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("Todos los archivos", "*.*")],
            title="Guardar gráfico como"
        )
        
        if filepath:
            try:
                self.figura.savefig(filepath, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Éxito", f"Gráfico guardado en:\n{filepath}")
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo guardar el gráfico:\n{str(e)}")
    
    def limpiar(self):
        for var in self.datos.values():
            if isinstance(var, (tk.DoubleVar, tk.BooleanVar)):
                if var == self.datos['resistencia']:
                    var.set(False)
                elif var == self.datos['g']:
                    var.set(9.81)
                elif var == self.datos['densidad_aire']:
                    var.set(1.225)
                elif var == self.datos['unidades']:
                    var.set("métrico")
                elif var == self.datos['comparar']:
                    var.set(False)
                else:
                    var.set(0)
        
        # Limpiar datos de comparación
        for var in self.datos_comparacion.values():
            if isinstance(var, (tk.DoubleVar, tk.BooleanVar)):
                if var == self.datos_comparacion['resistencia']:
                    var.set(False)
                elif var == self.datos_comparacion['g']:
                    var.set(9.81)
                elif var == self.datos_comparacion['densidad_aire']:
                    var.set(1.225)
                else:
                    var.set(0)
        
        self.resultados_text.delete(1.0, tk.END)
        self.energia_text.delete(1.0, tk.END)
        self.consulta_resultado.config(text="")
        self.tiempo_consulta.delete(0, tk.END)
        
        if self.figura:
            self.ax.clear()
            self.canvas.draw()
        
        if hasattr(self, 'animacion') and self.animacion:
            self.detener_animacion()
        
        # Limpiar información del punto
        if hasattr(self, 'lbl_posicion'):
            self.lbl_posicion.config(text="Posición: (0.00, 0.00) m | Tiempo: 0.00 s")
            self.lbl_velocidad.config(text="Velocidad: (0.00, 0.00) m/s | Magnitud: 0.00 m/s")
            self.lbl_altura_distancia.config(text="Altura: 0.00 m | Distancia: 0.00 m")
            self.lbl_angulo.config(text="Ángulo: 0.00°")
        
        if hasattr(self, 'punto_interactivo'):
            self.punto_interactivo.remove()
            del self.punto_interactivo
    
    def modo_avanzado(self):
        ventana_avanzada = tk.Toplevel(self.root)
        ventana_avanzada.title("Modo Avanzado - Edición de Ecuaciones")
        ventana_avanzada.geometry("600x400")
        
        # Frame principal
        frame_principal = ttk.Frame(ventana_avanzada)
        frame_principal.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Área de texto para ecuaciones
        ttk.Label(frame_principal, text="Editar ecuaciones del modelo:").pack(anchor="w")
        self.texto_ecuaciones = tk.Text(frame_principal, wrap=tk.WORD, height=15)
        self.texto_ecuaciones.pack(fill="both", expand=True, pady=5)
        
        # Insertar ecuaciones actuales
        ecuaciones = """# Ecuaciones del modelo de tiro parabólico

def descomposicion_vectorial(v0, theta):
    return v0 * np.cos(theta), v0 * np.sin(theta)

def tiempo_vuelo(v0y, g, y0=0):
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
        
        if y[-1] < 0:
            break
    
    return x, y, tiempos[:len(x)]
"""
        self.texto_ecuaciones.insert(tk.END, ecuaciones)
        
        # Botones
        frame_botones = ttk.Frame(frame_principal)
        frame_botones.pack(fill="x", pady=5)
        
        ttk.Button(frame_botones, text="Guardar Cambios", command=self.guardar_ecuaciones).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Restaurar Predeterminado", command=self.restaurar_ecuaciones).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Cancelar", command=ventana_avanzada.destroy).pack(side="right", padx=5)
    
    def guardar_ecuaciones(self):
        # Aquí deberías implementar la lógica para guardar las ecuaciones editadas
        # Esto es complejo porque implica modificar el código en tiempo de ejecución
        messagebox.showinfo("Información", "Esta funcionalidad es compleja y requiere implementación adicional.")
    
    def restaurar_ecuaciones(self):
        # Restaurar ecuaciones predeterminadas
        pass
    
    def cargar_historial(self):
        if os.path.exists(self.historial_file):
            try:
                with open(self.historial_file, 'r') as f:
                    self.historial = json.load(f)
            except:
                self.historial = []
    
    def guardar_historial(self):
        try:
            with open(self.historial_file, 'w') as f:
                json.dump(self.historial, f)
        except:
            pass
    
    def guardar_en_historial(self, datos):
        # Crear entrada de historial
        entrada = {
            'timestamp': time.time(),
            'datos': datos,
            'resultados': {
                't_vuelo': self.resultados['t_vuelo'],
                'alcance': self.resultados['alcance'],
                'altura_max': self.resultados['altura_max']
            }
        }
        
        # Agregar al historial
        self.historial.append(entrada)
        
        # Mantener solo las últimas 50 entradas
        if len(self.historial) > 50:
            self.historial = self.historial[-50:]
        
        # Guardar en archivo
        self.guardar_historial()
        
        # Actualizar lista
        self.actualizar_lista_historial()
    
    def actualizar_lista_historial(self):
        self.historial_listbox.delete(0, tk.END)
        
        for i, entrada in enumerate(self.historial):
            timestamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(entrada['timestamp']))
            datos = entrada['datos']
            res = entrada['resultados']
            
            texto = f"{timestamp} - v0={datos['v0']:.1f} m/s, ángulo={np.degrees(datos['angulo']):.1f}° | "
            texto += f"Alcance: {res['alcance']:.1f}m, Tiempo: {res['t_vuelo']:.1f}s"
            
            self.historial_listbox.insert(tk.END, texto)
    
    def cargar_desde_historial(self):
        seleccion = self.historial_listbox.curselection()
        if not seleccion:
            messagebox.showwarning("Advertencia", "Seleccione una simulación del historial")
            return
            
        indice = seleccion[0]
        entrada = self.historial[indice]
        datos = entrada['datos']
        
        # Cargar datos principales
        self.datos['x0'].set(datos['x0'])
        self.datos['y0'].set(datos['y0'])
        self.datos['v0'].set(datos['v0'])
        self.datos['angulo'].set(np.degrees(datos['angulo']))
        self.datos['g'].set(datos['g'])
        self.datos['resistencia'].set(datos['resistencia'])
        
        if datos['resistencia']:
            self.datos['masa'].set(datos['masa'])
            self.datos['coef_arrastre'].set(datos['coef_arrastre'])
            self.datos['area'].set(datos['area'])
            self.datos['densidad_aire'].set(datos['densidad_aire'])
        
        # Calcular automáticamente
        self.calcular_trayectoria()
    
    def eliminar_del_historial(self):
        seleccion = self.historial_listbox.curselection()
        if not seleccion:
            messagebox.showwarning("Advertencia", "Seleccione una simulación del historial")
            return
            
        indice = seleccion[0]
        del self.historial[indice]
        self.guardar_historial()
        self.actualizar_lista_historial()
    
    def limpiar_historial(self):
        if messagebox.askyesno("Confirmar", "¿Está seguro que desea borrar todo el historial?"):
            self.historial = []
            self.guardar_historial()
            self.actualizar_lista_historial()

# Funciones de cálculo
def descomposicion_vectorial(v0, theta):
    return v0 * np.cos(theta), v0 * np.sin(theta)

def tiempo_vuelo(v0y, g, y0=0):
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
        
        if y[-1] < 0:
            break
    
    return x, y, tiempos[:len(x)]

# Iniciar aplicación
if __name__ == "__main__":
    root = tk.Tk()
    app = TiroParabolicoApp(root)
    root.mainloop()