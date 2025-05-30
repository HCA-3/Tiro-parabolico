import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import json
import matplotlib.animation as animation
import os
import time

class TiroParabolicoApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Simulador de Tiro Parabólico Avanzado")
        self.root.geometry("1100x750")
        
        # Configure styles
        self.style = ttk.Style()
        self.style.configure('TNotebook.Tab', font=('Helvetica', 10, 'bold'))
        self.style.configure('Accent.TButton', font=('Helvetica', 10, 'bold'), foreground='blue')
        
        # Initialize variables
        self.initialize_variables()
        self.cargar_historial()
        self.crear_interfaz()
    
    def initialize_variables(self):
        # Main trajectory variables
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
        
        # Comparison trajectory variables
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
    
    def crear_interfaz(self):
        # Create notebook (tabs)
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True)
        
        # Create tabs
        self.crear_pestana_entrada(notebook)
        self.crear_pestana_resultados(notebook)
        self.crear_pestana_graficos(notebook)
        self.crear_pestana_historial(notebook)
        
        # Create toolbar
        self.crear_barra_herramientas()
    
    def crear_pestana_entrada(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Datos de Entrada")
        
        # Create panels
        self.crear_panel_basico(frame)
        self.crear_panel_resistencia(frame)
        self.crear_panel_comparacion(frame)
        self.crear_panel_unidades(frame)
    
    def crear_panel_basico(self, parent):
        frame = ttk.LabelFrame(parent, text="Parámetros Básicos", padding=10)
        frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        entries = [
            ("Posición inicial x0:", 'x0'),
            ("Posición inicial y0:", 'y0'),
            ("Velocidad inicial:", 'v0'),
            ("Ángulo de lanzamiento:", 'angulo'),
            ("Gravedad:", 'g')
        ]
        
        for i, (text, var) in enumerate(entries):
            ttk.Label(frame, text=text).grid(row=i, column=0, sticky="w")
            ttk.Entry(frame, textvariable=self.datos[var]).grid(row=i, column=1)
        
        ttk.Checkbutton(frame, text="Considerar resistencia del aire", 
                       variable=self.datos['resistencia'], 
                       command=self.toggle_resistencia).grid(row=5, column=0, columnspan=2, pady=5)
    
    def crear_panel_resistencia(self, parent):
        self.frame_resistencia = ttk.LabelFrame(parent, text="Parámetros de Resistencia del Aire", padding=10)
        self.frame_resistencia.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
        
        entries = [
            ("Masa del proyectil:", 'masa'),
            ("Coeficiente de arrastre:", 'coef_arrastre'),
            ("Área transversal:", 'area'),
            ("Densidad del aire:", 'densidad_aire')
        ]
        
        for i, (text, var) in enumerate(entries):
            ttk.Label(self.frame_resistencia, text=text).grid(row=i, column=0, sticky="w")
            ttk.Entry(self.frame_resistencia, textvariable=self.datos[var]).grid(row=i, column=1)
        
        self.toggle_resistencia()
    
    def crear_panel_comparacion(self, parent):
        self.frame_comparacion = ttk.LabelFrame(parent, text="Comparación de Tiros", padding=10)
        self.frame_comparacion.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")
        
        ttk.Checkbutton(self.frame_comparacion, text="Activar comparación", 
                       variable=self.datos['comparar'], 
                       command=self.toggle_comparacion).grid(row=0, column=0, columnspan=2, pady=5)
        
        self.frame_comparacion_params = ttk.Frame(self.frame_comparacion)
        self.frame_comparacion_params.grid(row=1, column=0, columnspan=2, sticky="nsew")
        
        entries = [
            ("Posición inicial x0:", 'x0'),
            ("Posición inicial y0:", 'y0'),
            ("Velocidad inicial:", 'v0'),
            ("Ángulo de lanzamiento:", 'angulo')
        ]
        
        for i, (text, var) in enumerate(entries):
            ttk.Label(self.frame_comparacion_params, text=text).grid(row=i, column=0, sticky="w")
            ttk.Entry(self.frame_comparacion_params, textvariable=self.datos_comparacion[var]).grid(row=i, column=1)
        
        ttk.Checkbutton(self.frame_comparacion_params, text="Resistencia del aire", 
                       variable=self.datos_comparacion['resistencia']).grid(row=4, column=0, columnspan=2, pady=5)
        
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
    
    def crear_pestana_resultados(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Resultados")
        
        # Create scrollable frame
        canvas = tk.Canvas(frame)
        scrollbar = ttk.Scrollbar(frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Results content
        self.resultados_text = tk.Text(scrollable_frame, wrap=tk.WORD, height=20, width=80)
        self.resultados_text.pack(padx=10, pady=10, fill="both", expand=True)
        
        # Time query frame
        frame_consulta = ttk.LabelFrame(scrollable_frame, text="Consulta por Tiempo", padding=10)
        frame_consulta.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(frame_consulta, text="Tiempo:").pack(side="left")
        self.tiempo_consulta = ttk.Entry(frame_consulta, width=10)
        self.tiempo_consulta.pack(side="left", padx=5)
        
        ttk.Button(frame_consulta, text="Consultar", command=self.consultar_tiempo).pack(side="left")
        
        self.consulta_resultado = ttk.Label(scrollable_frame, text="", wraplength=500)
        self.consulta_resultado.pack(padx=10, pady=5)
        
        # Energy frame
        frame_energia = ttk.LabelFrame(scrollable_frame, text="Energías", padding=10)
        frame_energia.pack(fill="x", padx=10, pady=5)
        
        self.energia_text = tk.Text(frame_energia, wrap=tk.WORD, height=5, width=80)
        self.energia_text.pack(fill="both", expand=True)
    
    def crear_pestana_graficos(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Gráficos")
        
        # Create figure and canvas
        self.figura = plt.Figure(figsize=(8, 4), dpi=100)
        self.ax = self.figura.add_subplot(111)
        
        # Main frame
        frame_principal = ttk.Frame(frame)
        frame_principal.pack(fill="both", expand=True)
        
        # Canvas for the plot
        self.canvas = FigureCanvasTkAgg(self.figura, master=frame_principal)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        # Controls frame
        frame_controles = ttk.Frame(frame_principal)
        frame_controles.pack(fill="x", padx=10, pady=5)
        
        # Point information frame
        self.frame_info = ttk.LabelFrame(frame_principal, text="Información del Punto", padding=10)
        self.frame_info.pack(fill="x", padx=10, pady=5)
        
        # Labels for trajectory 1 (blue)
        ttk.Label(self.frame_info, text="Trayectoria 1 (Azul):").pack(anchor="w")
        self.lbl_posicion1 = ttk.Label(self.frame_info, text="Posición: (0.00, 0.00) m | Tiempo: 0.00 s")
        self.lbl_posicion1.pack(anchor="w")
        self.lbl_velocidad1 = ttk.Label(self.frame_info, text="Velocidad: (0.00, 0.00) m/s | Magnitud: 0.00 m/s")
        self.lbl_velocidad1.pack(anchor="w")
        self.lbl_aceleracion1 = ttk.Label(self.frame_info, text="Aceleración: (0.00, -9.81) m/s² | Magnitud: 9.81 m/s²")
        self.lbl_aceleracion1.pack(anchor="w")
        
        # Separator
        ttk.Separator(self.frame_info, orient='horizontal').pack(fill='x', pady=5)
        
        # Labels for trajectory 2 (red)
        self.lbl_trayectoria2 = ttk.Label(self.frame_info, text="Trayectoria 2 (Roja):")
        self.lbl_posicion2 = ttk.Label(self.frame_info, text="Posición: (0.00, 0.00) m | Tiempo: 0.00 s")
        self.lbl_velocidad2 = ttk.Label(self.frame_info, text="Velocidad: (0.00, 0.00) m/s | Magnitud: 0.00 m/s")
        self.lbl_aceleracion2 = ttk.Label(self.frame_info, text="Aceleración: (0.00, -9.81) m/s² | Magnitud: 9.81 m/s²")
        
        # Slider for trajectory 1 (blue)
        self.slider_trayectoria_frame = ttk.Frame(frame_controles)
        self.slider_trayectoria_frame.pack(fill="x", pady=5)
        
        ttk.Label(self.slider_trayectoria_frame, text="Posición trayectoria 1:").pack(side="left")
        self.slider_trayectoria = ttk.Scale(
            self.slider_trayectoria_frame,
            from_=0,
            to=100,
            orient="horizontal",
            command=self.actualizar_punto_trayectoria
        )
        self.slider_trayectoria.pack(side="left", expand=True, fill="x", padx=5)
        
        # Slider for trajectory 2 (red) - initially hidden
        self.slider_comparacion_frame = ttk.Frame(frame_controles)
        ttk.Label(self.slider_comparacion_frame, text="Posición trayectoria 2:").pack(side="left")
        self.slider_comparacion = ttk.Scale(
            self.slider_comparacion_frame,
            from_=0,
            to=100,
            orient="horizontal",
            command=self.actualizar_punto_comparacion
        )
        self.slider_comparacion.pack(side="left", expand=True, fill="x", padx=5)
        
        # Control buttons
        ttk.Button(frame_controles, text="Exportar Gráfico", command=self.exportar_grafico).pack(side="left", padx=5)
        ttk.Button(frame_controles, text="Animación", command=self.toggle_animacion).pack(side="left", padx=5)
        
        # Animation speed slider
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
    
    def crear_pestana_historial(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Historial")
        
        # Create scrollable frame
        canvas = tk.Canvas(frame)
        scrollbar = ttk.Scrollbar(frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # History content
        self.historial_listbox = tk.Listbox(scrollable_frame, width=100, height=20)
        self.historial_listbox.pack(padx=10, pady=10, fill="both", expand=True)
        
        # History buttons
        frame_botones = ttk.Frame(scrollable_frame)
        frame_botones.pack(fill="x", padx=10, pady=5)
        
        ttk.Button(frame_botones, text="Cargar Simulación", command=self.cargar_desde_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Eliminar Selección", command=self.eliminar_del_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Limpiar Historial", command=self.limpiar_historial).pack(side="right", padx=5)
        
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
            pass  # Convert from imperial to metric if needed
        else:
            pass  # Convert from metric to imperial if needed
    
    def validar_datos(self, datos=None, es_comparacion=False):
        if datos is None:
            datos = {k: v.get() for k, v in self.datos.items() if k not in ['unidades', 'comparar']}
        
        try:
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
        if not self.validar_datos():
            return
        
        if self.datos['comparar'].get():
            datos_comparacion = {k: v.get() for k, v in self.datos_comparacion.items()}
            if not self.validar_datos(datos_comparacion, es_comparacion=True):
                return
        
        # Calculate main trajectory
        datos = {k: v.get() for k, v in self.datos.items() if k not in ['unidades', 'comparar']}
        datos['angulo'] = np.radians(datos['angulo'])
        self.resultados = self.calcular_trayectoria_individual(datos)
        
        # Calculate comparison trajectory if enabled
        if self.datos['comparar'].get():
            datos_comparacion = {k: v.get() for k, v in self.datos_comparacion.items()}
            datos_comparacion['angulo'] = np.radians(datos_comparacion['angulo'])
            self.resultados_comparacion = self.calcular_trayectoria_individual(datos_comparacion)
        
        self.actualizar_resultados(datos)
        self.actualizar_grafico(datos)
        self.guardar_en_historial(datos)
    
    def calcular_trayectoria_individual(self, datos):
        resultados = {}
        
        # Vector decomposition
        v0x, v0y = descomposicion_vectorial(datos['v0'], datos['angulo'])
        resultados['v0x'], resultados['v0y'] = v0x, v0y
        
        # Basic calculations
        resultados['t_vuelo'] = tiempo_vuelo(v0y, datos['g'], datos['y0'])
        resultados['alcance'] = alcance_maximo(v0x, resultados['t_vuelo'])
        resultados['altura_max'] = datos['y0'] + altura_maxima(v0y, datos['g'])
        
        # Simulation
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
        
        # Calculate energies
        resultados['energias'] = self.calcular_energias(datos, resultados)
        
        return resultados
    
    def calcular_energias(self, datos, resultados):
        energias = {
            'cinetica': [],
            'potencial': [],
            'total': []
        }
        
        g = datos['g']
        masa = datos['masa'] if datos['resistencia'] else 1  # Fictitious mass if no resistance
        
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
        
        # Plot main trajectory (blue)
        self.ax.plot(self.resultados['x'], self.resultados['y'], 'b-', label="Trayectoria 1")
        self.ax.plot(datos['x0'], datos['y0'], 'bo', label="Punto de lanzamiento 1")
        self.ax.plot(self.resultados['alcance'], 0, 'bx', label="Punto de impacto 1")
        
        # Show/hide comparison elements
        if self.datos['comparar'].get() and self.resultados_comparacion:
            # Plot comparison trajectory (red)
            self.ax.plot(self.resultados_comparacion['x'], self.resultados_comparacion['y'], 'r-', label="Trayectoria 2")
            self.ax.plot(self.datos_comparacion['x0'].get(), self.datos_comparacion['y0'].get(), 'ro', label="Punto de lanzamiento 2")
            self.ax.plot(self.resultados_comparacion['alcance'], 0, 'rx', label="Punto de impacto 2")
            
            # Show slider and labels for trajectory 2
            self.slider_comparacion_frame.pack(fill="x", pady=5)
            self.lbl_trayectoria2.pack(anchor="w")
            self.lbl_posicion2.pack(anchor="w")
            self.lbl_velocidad2.pack(anchor="w")
            self.lbl_aceleracion2.pack(anchor="w")
        else:
            # Hide comparison elements
            self.slider_comparacion_frame.pack_forget()
            self.lbl_trayectoria2.pack_forget()
            self.lbl_posicion2.pack_forget()
            self.lbl_velocidad2.pack_forget()
            self.lbl_aceleracion2.pack_forget()
        
        # Add legend for vectors
        self.ax.plot([], [], 'r-', label='Vector Velocidad')
        self.ax.plot([], [], 'g-', label='Vector Aceleración')
        if self.datos['comparar'].get() and self.resultados_comparacion:
            self.ax.plot([], [], 'c-', label='Vector Velocidad (Comp)')
            self.ax.plot([], [], 'm-', label='Vector Aceleración (Comp)')
        
        self.ax.set_title("Trayectoria del Proyectil")
        self.ax.set_xlabel("Distancia (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.grid(True)
        self.ax.legend(loc='upper right')
        
        # Configure sliders
        self.slider_trayectoria.set(0)
        if self.datos['comparar'].get() and self.resultados_comparacion:
            self.slider_comparacion.set(0)
        
        self.actualizar_punto_trayectoria(0)
        if self.datos['comparar'].get() and self.resultados_comparacion:
            self.actualizar_punto_comparacion(0)
        
        self.canvas.draw()
    
    def actualizar_punto_trayectoria(self, valor):
        if not hasattr(self, 'resultados') or not self.resultados:
            return
        
        idx = int(float(valor) / 100 * (len(self.resultados['x']) - 1))
        idx = max(0, min(idx, len(self.resultados['x']) - 1))
        
        x = self.resultados['x'][idx]
        y = self.resultados['y'][idx]
        t = self.resultados['tiempos'][idx]
        
        if idx == 0:
            vx = self.resultados['v0x']
            vy = self.resultados['v0y']
            ax = 0
            ay = -self.datos['g'].get()
        elif idx == 1:
            dt = self.resultados['tiempos'][idx] - self.resultados['tiempos'][idx-1]
            vx_prev = self.resultados['v0x']
            vy_prev = self.resultados['v0y']
            vx = (self.resultados['x'][idx] - self.resultados['x'][idx-1]) / dt
            vy = (self.resultados['y'][idx] - self.resultados['y'][idx-1]) / dt
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        else:
            dt = self.resultados['tiempos'][idx] - self.resultados['tiempos'][idx-1]
            dt_prev = self.resultados['tiempos'][idx-1] - self.resultados['tiempos'][idx-2]
            vx = (self.resultados['x'][idx] - self.resultados['x'][idx-1]) / dt
            vy = (self.resultados['y'][idx] - self.resultados['y'][idx-1]) / dt
            vx_prev = (self.resultados['x'][idx-1] - self.resultados['x'][idx-2]) / dt_prev
            vy_prev = (self.resultados['y'][idx-1] - self.resultados['y'][idx-2]) / dt_prev
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        
        velocidad_magnitud = np.sqrt(vx**2 + vy**2)
        aceleracion_magnitud = np.sqrt(ax**2 + ay**2)
        
        self.lbl_posicion1.config(text=f"Posición: ({x:.2f}, {y:.2f}) m | Tiempo: {t:.2f} s")
        self.lbl_velocidad1.config(text=f"Velocidad: ({vx:.2f}, {vy:.2f}) m/s | Magnitud: {velocidad_magnitud:.2f} m/s")
        self.lbl_aceleracion1.config(text=f"Aceleración: ({ax:.2f}, {ay:.2f}) m/s² | Magnitud: {aceleracion_magnitud:.2f} m/s²")
        
        # Actualizar o crear el punto interactivo
        if not hasattr(self, 'punto_interactivo'):
            self.punto_interactivo, = self.ax.plot([x], [y], 'ro', markersize=8, label="Posición actual 1")
        else:
            self.punto_interactivo.set_data([x], [y])
        
        # Línea desde el punto de lanzamiento hasta la posición actual (ahora con flecha)
        x0 = self.datos['x0'].get()
        y0 = self.datos['y0'].get()
        if not hasattr(self, 'linea_trayectoria'):
            self.linea_trayectoria = self.ax.arrow(
                x0, y0, x-x0, y-y0, 
                head_width=0.3, head_length=0.5, 
                fc='g', ec='g', linestyle='--', linewidth=1, alpha=0.5
            )
        else:
            self.linea_trayectoria.remove()
            self.linea_trayectoria = self.ax.arrow(
                x0, y0, x-x0, y-y0, 
                head_width=0.3, head_length=0.5, 
                fc='g', ec='g', linestyle='--', linewidth=1, alpha=0.5
            )
        
        # Línea de dirección (vector velocidad) con punta de flecha
        scale = 0.5  # Escala para visualización del vector
        if not hasattr(self, 'linea_direccion'):
            self.linea_direccion = self.ax.arrow(
                x, y, vx*scale, vy*scale, 
                head_width=0.5, head_length=0.8, 
                fc='r', ec='r', linewidth=2, alpha=0.7
            )
        else:
            # Eliminar la flecha anterior y crear una nueva
            self.linea_direccion.remove()
            self.linea_direccion = self.ax.arrow(
                x, y, vx*scale, vy*scale, 
                head_width=0.5, head_length=0.8, 
                fc='r', ec='r', linewidth=2, alpha=0.7
            )
        
        # Vector aceleración (solo gravedad o gravedad + resistencia)
        ay_total = -self.datos['g'].get()
        if self.datos['resistencia'].get() and idx > 0:
            # Calcular aceleración por resistencia del aire
            velocidad = np.sqrt(vx**2 + vy**2)
            fuerza_arrastre = 0.5 * self.datos['densidad_aire'].get() * \
                             self.datos['coef_arrastre'].get() * \
                             self.datos['area'].get() * velocidad**2
            ax_res = - (fuerza_arrastre / self.datos['masa'].get()) * (vx / velocidad)
            ay_res = - (fuerza_arrastre / self.datos['masa'].get()) * (vy / velocidad)
            ay_total += ay_res
        
        # Mostrar vector aceleración (en verde)
        scale_a = 0.2  # Escala diferente para la aceleración
        if not hasattr(self, 'linea_aceleracion'):
            self.linea_aceleracion = self.ax.arrow(
                x, y, 0, ay_total*scale_a, 
                head_width=0.3, head_length=0.5, 
                fc='g', ec='g', linewidth=2, alpha=0.7
            )
        else:
            self.linea_aceleracion.remove()
            self.linea_aceleracion = self.ax.arrow(
                x, y, 0, ay_total*scale_a, 
                head_width=0.3, head_length=0.5, 
                fc='g', ec='g', linewidth=2, alpha=0.7
            )
        
        self.canvas.draw()
    
    def actualizar_punto_comparacion(self, valor):
        if not hasattr(self, 'resultados_comparacion') or not self.resultados_comparacion:
            return
        
        idx = int(float(valor) / 100 * (len(self.resultados_comparacion['x']) - 1))
        idx = max(0, min(idx, len(self.resultados_comparacion['x']) - 1))
        
        x = self.resultados_comparacion['x'][idx]
        y = self.resultados_comparacion['y'][idx]
        t = self.resultados_comparacion['tiempos'][idx]
        
        if idx == 0:
            vx = self.resultados_comparacion['v0x']
            vy = self.resultados_comparacion['v0y']
            ax = 0
            ay = -self.datos_comparacion['g'].get()
        elif idx == 1:
            dt = self.resultados_comparacion['tiempos'][idx] - self.resultados_comparacion['tiempos'][idx-1]
            vx_prev = self.resultados_comparacion['v0x']
            vy_prev = self.resultados_comparacion['v0y']
            vx = (self.resultados_comparacion['x'][idx] - self.resultados_comparacion['x'][idx-1]) / dt
            vy = (self.resultados_comparacion['y'][idx] - self.resultados_comparacion['y'][idx-1]) / dt
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        else:
            dt = self.resultados_comparacion['tiempos'][idx] - self.resultados_comparacion['tiempos'][idx-1]
            dt_prev = self.resultados_comparacion['tiempos'][idx-1] - self.resultados_comparacion['tiempos'][idx-2]
            vx = (self.resultados_comparacion['x'][idx] - self.resultados_comparacion['x'][idx-1]) / dt
            vy = (self.resultados_comparacion['y'][idx] - self.resultados_comparacion['y'][idx-1]) / dt
            vx_prev = (self.resultados_comparacion['x'][idx-1] - self.resultados_comparacion['x'][idx-2]) / dt_prev
            vy_prev = (self.resultados_comparacion['y'][idx-1] - self.resultados_comparacion['y'][idx-2]) / dt_prev
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        
        velocidad_magnitud = np.sqrt(vx**2 + vy**2)
        aceleracion_magnitud = np.sqrt(ax**2 + ay**2)
        
        self.lbl_posicion2.config(text=f"Posición: ({x:.2f}, {y:.2f}) m | Tiempo: {t:.2f} s")
        self.lbl_velocidad2.config(text=f"Velocidad: ({vx:.2f}, {vy:.2f}) m/s | Magnitud: {velocidad_magnitud:.2f} m/s")
        self.lbl_aceleracion2.config(text=f"Aceleración: ({ax:.2f}, {ay:.2f}) m/s² | Magnitud: {aceleracion_magnitud:.2f} m/s²")
        
        if not hasattr(self, 'punto_comparacion'):
            self.punto_comparacion, = self.ax.plot([x], [y], 'go', markersize=8, label="Posición actual 2")
        else:
            self.punto_comparacion.set_data([x], [y])
        
        # Línea desde el punto de lanzamiento hasta la posición actual (comparación)
        x0_comp = self.datos_comparacion['x0'].get()
        y0_comp = self.datos_comparacion['y0'].get()
        if not hasattr(self, 'linea_trayectoria_comp'):
            self.linea_trayectoria_comp = self.ax.arrow(
                x0_comp, y0_comp, x-x0_comp, y-y0_comp,
                head_width=0.3, head_length=0.5,
                fc='m', ec='m', linestyle='--', linewidth=1, alpha=0.5
            )
        else:
            self.linea_trayectoria_comp.remove()
            self.linea_trayectoria_comp = self.ax.arrow(
                x0_comp, y0_comp, x-x0_comp, y-y0_comp,
                head_width=0.3, head_length=0.5,
                fc='m', ec='m', linestyle='--', linewidth=1, alpha=0.5
            )
        
        # Línea de dirección (vector velocidad) con punta de flecha para comparación
        scale = 0.5
        if not hasattr(self, 'linea_direccion_comp'):
            self.linea_direccion_comp = self.ax.arrow(
                x, y, vx*scale, vy*scale, 
                head_width=0.5, head_length=0.8, 
                fc='c', ec='c', linewidth=2, alpha=0.7
            )
        else:
            self.linea_direccion_comp.remove()
            self.linea_direccion_comp = self.ax.arrow(
                x, y, vx*scale, vy*scale, 
                head_width=0.5, head_length=0.8, 
                fc='c', ec='c', linewidth=2, alpha=0.7
            )
        
        # Vector aceleración para comparación
        ay_total = -self.datos_comparacion['g'].get()
        if self.datos_comparacion['resistencia'].get() and idx > 0:
            velocidad = np.sqrt(vx**2 + vy**2)
            fuerza_arrastre = 0.5 * self.datos_comparacion['densidad_aire'].get() * \
                             self.datos_comparacion['coef_arrastre'].get() * \
                             self.datos_comparacion['area'].get() * velocidad**2
            ax_res = - (fuerza_arrastre / self.datos_comparacion['masa'].get()) * (vx / velocidad)
            ay_res = - (fuerza_arrastre / self.datos_comparacion['masa'].get()) * (vy / velocidad)
            ay_total += ay_res
        
        if not hasattr(self, 'linea_aceleracion_comp'):
            self.linea_aceleracion_comp = self.ax.arrow(
                x, y, 0, ay_total*0.2, 
                head_width=0.3, head_length=0.5, 
                fc='m', ec='m', linewidth=2, alpha=0.7
            )
        else:
            self.linea_aceleracion_comp.remove()
            self.linea_aceleracion_comp = self.ax.arrow(
                x, y, 0, ay_total*0.2, 
                head_width=0.3, head_length=0.5, 
                fc='m', ec='m', linewidth=2, alpha=0.7
            )
        
        self.canvas.draw()
    
    def iniciar_animacion(self):
        if not self.resultados:
            messagebox.showwarning("Advertencia", "No hay datos para animar. Calcule primero una trayectoria.")
            return
            
        # Animation setup
        self.linea_animada, = self.ax.plot([], [], 'ko', markersize=8)
        self.tiempo_texto = self.ax.text(0.02, 0.95, '', transform=self.ax.transAxes)
        
        def init():
            self.linea_animada.set_data([], [])
            self.tiempo_texto.set_text('')
            return self.linea_animada, self.tiempo_texto
        
        def animate(i):
            if i >= len(self.resultados['x']):
                i = len(self.resultados['x']) - 1
                
            x = self.resultados['x'][i]
            y = self.resultados['y'][i]
            t = self.resultados['tiempos'][i]
            
            self.linea_animada.set_data(x, y)
            self.tiempo_texto.set_text(f'Tiempo: {t:.2f} s')
            
            # Update slider and labels
            self.slider_trayectoria.set(i / len(self.resultados['x']) * 100)
            self.actualizar_punto_trayectoria(i / len(self.resultados['x']) * 100)
            
            return self.linea_animada, self.tiempo_texto
        
        # Calculate interval based on speed slider
        intervalo = 110 - self.velocidad_animacion.get()  # Invert value (1-100 to 110-10)
        
        # Create animation
        self.animacion = animation.FuncAnimation(
            self.figura, animate, frames=len(self.resultados['x'])*2,  # Multiply by 2 for smoother animation
            init_func=init, blit=True, interval=intervalo, repeat=True
        )
        
        self.canvas.draw()
    
    def detener_animacion(self):
        if hasattr(self, 'animacion') and self.animacion:
            self.animacion.event_source.stop()
            self.animacion = None
            
            # Clean animation elements
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
            angulo = np.degrees(datos['angulo']) if 'angulo' in datos else datos['angulo']
            
            # Calcular parametrización
            param = parametrizacion_tiro_parabolico(
                datos['x0'], datos['y0'], datos['v0'], angulo, 
                datos['g'], t
            )
            
            if param['posicion'][1] < 0:
                resultado = f"En t = {t:.2f} s: El proyectil ya ha impactado en el suelo."
            else:
                resultado = (
                    f"En t = {t:.2f} s:\n"
                    f"- Posición: ({param['posicion'][0]:.2f}, {param['posicion'][1]:.2f}) m\n"
                    f"- Velocidad: ({param['velocidad'][0]:.2f}, {param['velocidad'][1]:.2f}) m/s\n"
                    f"- Aceleración: ({param['aceleracion'][0]:.2f}, {param['aceleracion'][1]:.2f}) m/s²\n\n"
                    f"Ecuaciones paramétricas:\n"
                    f"x(t) = {param['ecuaciones']['x(t)']}\n"
                    f"y(t) = {param['ecuaciones']['y(t)']}\n"
                    f"vx(t) = {param['ecuaciones']['vx(t)']}\n"
                    f"vy(t) = {param['ecuaciones']['vy(t)']}"
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
                    try:
                        from fpdf import FPDF
                        pdf = FPDF()
                        pdf.add_page()
                        pdf.set_font("Arial", size=12)
                        
                        for line in contenido.split('\n'):
                            pdf.cell(200, 10, txt=line, ln=1)
                        
                        pdf.output(filepath)
                    except ImportError:
                        messagebox.showerror("Error", "Para exportar a PDF, instale fpdf2 con: pip install fpdf2")
                        return
                else:
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
        # Reset variables
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
        
        # Reset comparison variables
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
        
        # Clear text fields
        self.resultados_text.delete(1.0, tk.END)
        self.energia_text.delete(1.0, tk.END)
        self.consulta_resultado.config(text="")
        self.tiempo_consulta.delete(0, tk.END)
        
        # Clear plot
        if self.figura:
            self.ax.clear()
            self.canvas.draw()
        
        # Stop animation
        if hasattr(self, 'animacion') and self.animacion:
            self.detener_animacion()
        
        # Clear interactive points and lines
        if hasattr(self, 'punto_interactivo'):
            self.punto_interactivo.remove()
            del self.punto_interactivo
        
        if hasattr(self, 'punto_comparacion'):
            self.punto_comparacion.remove()
            del self.punto_comparacion
        
        if hasattr(self, 'linea_trayectoria'):
            self.linea_trayectoria.remove()
            del self.linea_trayectoria
        
        if hasattr(self, 'linea_direccion'):
            self.linea_direccion.remove()
            del self.linea_direccion
        
        if hasattr(self, 'linea_trayectoria_comp'):
            self.linea_trayectoria_comp.remove()
            del self.linea_trayectoria_comp
        
        if hasattr(self, 'linea_direccion_comp'):
            self.linea_direccion_comp.remove()
            del self.linea_direccion_comp
        
        # Limpiar vectores de aceleración si existen
        if hasattr(self, 'linea_aceleracion'):
            self.linea_aceleracion.remove()
            del self.linea_aceleracion
        
        if hasattr(self, 'linea_aceleracion_comp'):
            self.linea_aceleracion_comp.remove()
            del self.linea_aceleracion_comp
        
        # Reset labels
        self.lbl_posicion1.config(text="Posición: (0.00, 0.00) m | Tiempo: 0.00 s")
        self.lbl_velocidad1.config(text="Velocidad: (0.00, 0.00) m/s | Magnitud: 0.00 m/s")
        self.lbl_aceleracion1.config(text="Aceleración: (0.00, -9.81) m/s² | Magnitud: 9.81 m/s²")
        self.lbl_posicion2.config(text="Posición: (0.00, 0.00) m | Tiempo: 0.00 s")
        self.lbl_velocidad2.config(text="Velocidad: (0.00, 0.00) m/s | Magnitud: 0.00 m/s")
        self.lbl_aceleracion2.config(text="Aceleración: (0.00, -9.81) m/s² | Magnitud: 9.81 m/s²")
    
    def modo_avanzado(self):
        ventana_avanzada = tk.Toplevel(self.root)
        ventana_avanzada.title("Modo Avanzado - Edición de Ecuaciones")
        ventana_avanzada.geometry("600x400")
        
        frame_principal = ttk.Frame(ventana_avanzada)
        frame_principal.pack(fill="both", expand=True, padx=10, pady=10)
        
        ttk.Label(frame_principal, text="Editar ecuaciones del modelo:").pack(anchor="w")
        self.texto_ecuaciones = tk.Text(frame_principal, wrap=tk.WORD, height=15)
        self.texto_ecuaciones.pack(fill="both", expand=True, pady=5)
        
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

def parametrizacion_tiro_parabolico(x0, y0, v0, angulo, g, t):
    \"\"\"
    Devuelve la parametrización del tiro parabólico en el tiempo t
    \"\"\"
    ang_rad = np.radians(angulo)
    v0x = v0 * np.cos(ang_rad)
    v0y = v0 * np.sin(ang_rad)
    
    x = x0 + v0x * t
    y = y0 + v0y * t - 0.5 * g * t**2
    vx = v0x
    vy = v0y - g * t
    ax = 0
    ay = -g
    
    return {
        'posicion': (x, y),
        'velocidad': (vx, vy),
        'aceleracion': (ax, ay),
        'tiempo': t,
        'ecuaciones': {
            'x(t)': f"{x0:.2f} + {v0x:.2f}·t",
            'y(t)': f"{y0:.2f} + {v0y:.2f}·t - 0.5·{g:.2f}·t²",
            'vx(t)': f"{v0x:.2f}",
            'vy(t)': f"{v0y:.2f} - {g:.2f}·t"
        }
    }
"""
        self.texto_ecuaciones.insert(tk.END, ecuaciones)
        
        frame_botones = ttk.Frame(frame_principal)
        frame_botones.pack(fill="x", pady=5)
        
        ttk.Button(frame_botones, text="Guardar Cambios", command=self.guardar_ecuaciones).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Restaurar Predeterminado", command=self.restaurar_ecuaciones).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Cancelar", command=ventana_avanzada.destroy).pack(side="right", padx=5)
    
    def guardar_ecuaciones(self):
        messagebox.showinfo("Información", "Esta funcionalidad es compleja y requiere implementación adicional.")
    
    def restaurar_ecuaciones(self):
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
        entrada = {
            'timestamp': time.time(),
            'datos': datos,
            'resultados': {
                't_vuelo': self.resultados['t_vuelo'],
                'alcance': self.resultados['alcance'],
                'altura_max': self.resultados['altura_max']
            }
        }
        
        self.historial.append(entrada)
        
        if len(self.historial) > 50:
            self.historial = self.historial[-50:]
        
        self.guardar_historial()
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

# Physics calculation functions
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

def parametrizacion_tiro_parabolico(x0, y0, v0, angulo, g, t):
    """
    Devuelve la parametrización del tiro parabólico en el tiempo t
    """
    ang_rad = np.radians(angulo)
    v0x = v0 * np.cos(ang_rad)
    v0y = v0 * np.sin(ang_rad)
    
    x = x0 + v0x * t
    y = y0 + v0y * t - 0.5 * g * t**2
    vx = v0x
    vy = v0y - g * t
    ax = 0
    ay = -g
    
    return {
        'posicion': (x, y),
        'velocidad': (vx, vy),
        'aceleracion': (ax, ay),
        'tiempo': t,
        'ecuaciones': {
            'x(t)': f"{x0:.2f} + {v0x:.2f}·t",
            'y(t)': f"{y0:.2f} + {v0y:.2f}·t - 0.5·{g:.2f}·t²",
            'vx(t)': f"{v0x:.2f}",
            'vy(t)': f"{v0y:.2f} - {g:.2f}·t"
        }
    }

if __name__ == "__main__":
    root = tk.Tk()
    app = TiroParabolicoApp(root)
    root.mainloop()