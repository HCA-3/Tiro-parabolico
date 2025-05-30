import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import json
import matplotlib.animation as animation
import os
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from matplotlib.patches import FancyArrowPatch

# ==============================================
# Clases para el modelo físico (Dominio)
# ==============================================

@dataclass
class ParametrosTiro:
    """Clase para almacenar los parámetros de un tiro parabólico"""
    x0: float = 0
    y0: float = 0
    v0: float = 10
    angulo: float = 45  # en grados
    g: float = 9.81
    resistencia: bool = False
    masa: float = 1
    coef_arrastre: float = 0.47
    area: float = 0.1
    densidad_aire: float = 1.225


@dataclass
class ResultadosTiro:
    """Clase para almacenar los resultados de un cálculo de tiro parabólico"""
    x: List[float]
    y: List[float]
    tiempos: List[float]
    v0x: float
    v0y: float
    t_vuelo: float
    alcance: float
    altura_max: float
    vx_final: float
    vy_final: float
    energias: Dict[str, List[float]]


class ModeloTiroParabolico:
    """Clase que encapsula la lógica de cálculo del tiro parabólico"""
    
    @staticmethod
    def calcular_trayectoria(parametros: ParametrosTiro) -> ResultadosTiro:
        """Calcula la trayectoria basada en los parámetros dados"""
        angulo_rad = np.radians(parametros.angulo)
        v0x, v0y = ModeloTiroParabolico.descomposicion_vectorial(parametros.v0, angulo_rad)
        
        t_vuelo = ModeloTiroParabolico.tiempo_vuelo(v0y, parametros.g, parametros.y0)
        alcance = ModeloTiroParabolico.alcance_maximo(v0x, t_vuelo)
        altura_max = parametros.y0 + ModeloTiroParabolico.altura_maxima(v0y, parametros.g)
        
        if parametros.resistencia:
            x, y, tiempos = ModeloTiroParabolico.movimiento_con_resistencia(
                parametros.x0, parametros.y0, v0x, v0y, 
                parametros.g, parametros.masa, parametros.coef_arrastre, 
                parametros.area, parametros.densidad_aire, t_vuelo
            )
            vx_final = (x[-1] - x[-2]) / 0.01 if len(x) > 1 else 0
            vy_final = (y[-1] - y[-2]) / 0.01 if len(y) > 1 else 0
        else:
            x = [parametros.x0 + v0x * t for t in np.linspace(0, t_vuelo, 100)]
            y = [parametros.y0 + v0y * t - 0.5 * parametros.g * t**2 for t in np.linspace(0, t_vuelo, 100)]
            tiempos = np.linspace(0, t_vuelo, len(x))
            vx_final = v0x
            vy_final = v0y - parametros.g * t_vuelo
        
        energias = ModeloTiroParabolico.calcular_energias(parametros, x, y, tiempos)
        
        return ResultadosTiro(
            x=x, y=y, tiempos=tiempos,
            v0x=v0x, v0y=v0y,
            t_vuelo=t_vuelo, alcance=alcance, altura_max=altura_max,
            vx_final=vx_final, vy_final=vy_final,
            energias=energias
        )
    
    @staticmethod
    def descomposicion_vectorial(v0: float, theta: float) -> Tuple[float, float]:
        return v0 * np.cos(theta), v0 * np.sin(theta)
    
    @staticmethod
    def tiempo_vuelo(v0y: float, g: float, y0: float = 0) -> float:
        return (v0y + np.sqrt(v0y**2 + 2 * g * y0)) / g
    
    @staticmethod
    def alcance_maximo(v0x: float, t_vuelo: float) -> float:
        return v0x * t_vuelo
    
    @staticmethod
    def altura_maxima(v0y: float, g: float) -> float:
        return (v0y**2) / (2 * g)
    
    @staticmethod
    def movimiento_con_resistencia(x0: float, y0: float, v0x: float, v0y: float, 
                                 g: float, masa: float, coef_arrastre: float, 
                                 area: float, densidad_aire: float, t_max: float, 
                                 dt: float = 0.01) -> Tuple[List[float], List[float], List[float]]:
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
    
    @staticmethod
    def calcular_energias(parametros: ParametrosTiro, x: List[float], y: List[float], 
                         tiempos: List[float]) -> Dict[str, List[float]]:
        energias = {
            'cinetica': [],
            'potencial': [],
            'total': []
        }
        
        g = parametros.g
        masa = parametros.masa if parametros.resistencia else 1
        
        for i in range(len(x)):
            if i == 0:
                vx = (x[1] - x[0]) / (tiempos[1] - tiempos[0]) if len(x) > 1 else 0
                vy = (y[1] - y[0]) / (tiempos[1] - tiempos[0]) if len(y) > 1 else 0
            else:
                dt = tiempos[i] - tiempos[i-1]
                vx = (x[i] - x[i-1]) / dt
                vy = (y[i] - y[i-1]) / dt
            
            velocidad = np.sqrt(vx**2 + vy**2)
            altura = y[i]
            
            ec = 0.5 * masa * velocidad**2
            ep = masa * g * altura
            et = ec + ep
            
            energias['cinetica'].append(ec)
            energias['potencial'].append(ep)
            energias['total'].append(et)
        
        return energias


# ==============================================
# Clases para la interfaz de usuario
# ==============================================

class PanelBase(ttk.Frame, ABC):
    """Clase base abstracta para los paneles de la interfaz"""
    
    @abstractmethod
    def crear_widgets(self):
        pass
    
    @abstractmethod
    def obtener_datos(self) -> Dict:
        pass
    
    @abstractmethod
    def establecer_datos(self, datos: Dict):
        pass


class PanelParametrosBasicos(PanelBase):
    """Panel para los parámetros básicos del tiro parabólico"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
    
    def crear_widgets(self):
        # Frame principal con scroll
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True)
        
        # Create scrollable frame
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        frame = ttk.LabelFrame(scrollable_frame, text="Parámetros Básicos", padding=10)
        frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.vars = {
            'x0': tk.DoubleVar(value=0),
            'y0': tk.DoubleVar(value=0),
            'v0': tk.DoubleVar(value=10),
            'angulo': tk.DoubleVar(value=45),
            'g': tk.DoubleVar(value=9.81),
            'resistencia': tk.BooleanVar(value=False)
        }
        
        entries = [
            ("Posición inicial x0:", 'x0'),
            ("Posición inicial y0:", 'y0'),
            ("Velocidad inicial:", 'v0'),
            ("Ángulo de lanzamiento:", 'angulo'),
            ("Gravedad:", 'g')
        ]
        
        for i, (text, var) in enumerate(entries):
            ttk.Label(frame, text=text).grid(row=i, column=0, sticky="w")
            ttk.Entry(frame, textvariable=self.vars[var]).grid(row=i, column=1)
        
        ttk.Checkbutton(frame, text="Considerar resistencia del aire", 
                       variable=self.vars['resistencia'], 
                       command=self.controller.toggle_resistencia).grid(row=5, column=0, columnspan=2, pady=5)
    
    def obtener_datos(self) -> Dict:
        return {k: v.get() for k, v in self.vars.items()}
    
    def establecer_datos(self, datos: Dict):
        for k, v in datos.items():
            if k in self.vars:
                self.vars[k].set(v)


class PanelResistenciaAire(PanelBase):
    """Panel para los parámetros de resistencia del aire"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
    
    def crear_widgets(self):
        # Frame principal con scroll
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True)
        
        # Create scrollable frame
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        frame = ttk.LabelFrame(scrollable_frame, text="Parámetros de Resistencia del Aire", padding=10)
        frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.vars = {
            'masa': tk.DoubleVar(value=1),
            'coef_arrastre': tk.DoubleVar(value=0.47),
            'area': tk.DoubleVar(value=0.1),
            'densidad_aire': tk.DoubleVar(value=1.225)
        }
        
        entries = [
            ("Masa del proyectil:", 'masa'),
            ("Coeficiente de arrastre:", 'coef_arrastre'),
            ("Área transversal:", 'area'),
            ("Densidad del aire:", 'densidad_aire')
        ]
        
        for i, (text, var) in enumerate(entries):
            ttk.Label(frame, text=text).grid(row=i, column=0, sticky="w")
            ttk.Entry(frame, textvariable=self.vars[var]).grid(row=i, column=1)
    
    def obtener_datos(self) -> Dict:
        return {k: v.get() for k, v in self.vars.items()}
    
    def establecer_datos(self, datos: Dict):
        for k, v in datos.items():
            if k in self.vars:
                self.vars[k].set(v)


class PanelComparacion(PanelBase):
    """Panel para la comparación de tiros"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
    
    def crear_widgets(self):
        # Frame principal con scroll
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True)
        
        # Create scrollable frame
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        frame = ttk.LabelFrame(scrollable_frame, text="Comparación de Tiros", padding=10)
        frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.vars = {
            'comparar': tk.BooleanVar(value=False),
            'x0': tk.DoubleVar(value=0),
            'y0': tk.DoubleVar(value=0),
            'v0': tk.DoubleVar(value=15),
            'angulo': tk.DoubleVar(value=30),
            'resistencia': tk.BooleanVar(value=False)
        }
        
        ttk.Checkbutton(frame, text="Activar comparación", 
                       variable=self.vars['comparar'], 
                       command=self.controller.toggle_comparacion).pack(anchor="w", pady=5)
        
        self.frame_params = ttk.Frame(frame)
        self.frame_params.pack(fill="x", padx=5, pady=5)
        
        entries = [
            ("Posición inicial x0:", 'x0'),
            ("Posición inicial y0:", 'y0'),
            ("Velocidad inicial:", 'v0'),
            ("Ángulo de lanzamiento:", 'angulo')
        ]
        
        for i, (text, var) in enumerate(entries):
            ttk.Label(self.frame_params, text=text).grid(row=i, column=0, sticky="w")
            ttk.Entry(self.frame_params, textvariable=self.vars[var]).grid(row=i, column=1)
        
        ttk.Checkbutton(self.frame_params, text="Resistencia del aire", 
                       variable=self.vars['resistencia']).grid(row=4, column=0, columnspan=2, pady=5)
    
    def obtener_datos(self) -> Dict:
        return {k: v.get() for k, v in self.vars.items()}
    
    def establecer_datos(self, datos: Dict):
        for k, v in datos.items():
            if k in self.vars:
                self.vars[k].set(v)


class PanelUnidades(PanelBase):
    """Panel para la configuración de unidades"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
    
    def crear_widgets(self):
        # Frame principal con scroll
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True)
        
        # Create scrollable frame
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        frame = ttk.LabelFrame(scrollable_frame, text="Configuración de Unidades", padding=10)
        frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.var_unidades = tk.StringVar(value="métrico")
        
        ttk.Radiobutton(frame, text="Sistema Métrico (m, kg)", 
                       variable=self.var_unidades, value="métrico",
                       command=self.controller.actualizar_unidades).pack(anchor="w")
        
        ttk.Radiobutton(frame, text="Sistema Imperial (pies, lb)", 
                       variable=self.var_unidades, value="imperial",
                       command=self.controller.actualizar_unidades).pack(anchor="w")
    
    def obtener_datos(self) -> Dict:
        return {'unidades': self.var_unidades.get()}
    
    def establecer_datos(self, datos: Dict):
        if 'unidades' in datos:
            self.var_unidades.set(datos['unidades'])


class PanelResultados(PanelBase):
    """Panel para mostrar los resultados del cálculo"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
    
    def crear_widgets(self):
        # Frame principal con scroll
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True)
        
        # Create scrollable frame
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
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
        
        ttk.Button(frame_consulta, text="Consultar", command=self.controller.consultar_tiempo).pack(side="left")
        
        self.consulta_resultado = ttk.Label(scrollable_frame, text="", wraplength=500)
        self.consulta_resultado.pack(padx=10, pady=5)
        
        # Energy frame
        frame_energia = ttk.LabelFrame(scrollable_frame, text="Energías", padding=10)
        frame_energia.pack(fill="x", padx=10, pady=5)
        
        self.energia_text = tk.Text(frame_energia, wrap=tk.WORD, height=5, width=80)
        self.energia_text.pack(fill="both", expand=True)
    
    def obtener_datos(self) -> Dict:
        return {}
    
    def establecer_datos(self, datos: Dict):
        pass
    
    def actualizar_resultados(self, informe: str):
        self.resultados_text.delete(1.0, tk.END)
        self.resultados_text.insert(tk.END, informe)
    
    def actualizar_energias(self, texto_energia: str):
        self.energia_text.delete(1.0, tk.END)
        self.energia_text.insert(tk.END, texto_energia)
    
    def mostrar_consulta(self, resultado: str):
        self.consulta_resultado.config(text=resultado)


class PanelGraficos(PanelBase):
    """Panel para mostrar los gráficos de la trayectoria"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
    
    def crear_widgets(self):
        # Frame principal con scroll
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True)
        
        # Create scrollable frame
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Create figure and canvas
        self.figura = plt.Figure(figsize=(8, 4), dpi=100)
        self.ax = self.figura.add_subplot(111)
        
        # Canvas for the plot
        self.canvas = FigureCanvasTkAgg(self.figura, master=scrollable_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        # Controls frame
        frame_controles = ttk.Frame(scrollable_frame)
        frame_controles.pack(fill="x", padx=10, pady=5)
        
        # Point information frame
        self.frame_info = ttk.LabelFrame(scrollable_frame, text="Información del Punto", padding=10)
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
            command=self.controller.actualizar_punto_trayectoria
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
            command=self.controller.actualizar_punto_comparacion
        )
        self.slider_comparacion.pack(side="left", expand=True, fill="x", padx=5)
        
        # Control buttons
        ttk.Button(frame_controles, text="Exportar Gráfico", command=self.controller.exportar_grafico).pack(side="left", padx=5)
        ttk.Button(frame_controles, text="Animación", command=self.controller.toggle_animacion).pack(side="left", padx=5)
        
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
    
    def obtener_datos(self) -> Dict:
        return {}
    
    def establecer_datos(self, datos: Dict):
        pass
    
    def actualizar_grafico(self, x: List[float], y: List[float], alcance: float, 
                          x_comparacion: Optional[List[float]] = None, 
                          y_comparacion: Optional[List[float]] = None, 
                          alcance_comparacion: Optional[float] = None):
        self.ax.clear()
        
        # Plot main trajectory (blue)
        self.ax.plot(x, y, 'b-', label="Trayectoria 1")
        self.ax.plot(x[0], y[0], 'bo', label="Punto de lanzamiento 1")
        self.ax.plot(alcance, 0, 'bx', label="Punto de impacto 1")
        
        # Show/hide comparison elements
        if x_comparacion is not None and y_comparacion is not None and alcance_comparacion is not None:
            self.ax.plot(x_comparacion, y_comparacion, 'r-', label="Trayectoria 2")
            self.ax.plot(x_comparacion[0], y_comparacion[0], 'ro', label="Punto de lanzamiento 2")
            self.ax.plot(alcance_comparacion, 0, 'rx', label="Punto de impacto 2")
            
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
        
        # Eliminar elementos de vectores anteriores si existen
        if hasattr(self, 'punto_interactivo'):
            self.punto_interactivo.remove()
            del self.punto_interactivo
        if hasattr(self, 'linea_trayectoria'):
            self.linea_trayectoria.remove()
            del self.linea_trayectoria
        if hasattr(self, 'linea_direccion'):
            self.linea_direccion.remove()
            del self.linea_direccion
        if hasattr(self, 'linea_aceleracion'):
            self.linea_aceleracion.remove()
            del self.linea_aceleracion
        
        self.ax.set_title("Trayectoria del Proyectil")
        self.ax.set_xlabel("Distancia (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.grid(True)
        self.ax.legend()
        
        # Configure sliders
        self.slider_trayectoria.set(0)
        if x_comparacion is not None:
            self.slider_comparacion.set(0)
        
        self.canvas.draw()
    
    def actualizar_punto_trayectoria(self, x: float, y: float, t: float, 
                                   vx: float, vy: float, ax: float, ay: float):
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
        
        # Línea desde el punto de lanzamiento hasta la posición actual (con flecha)
        x0 = self.controller.parametros_actuales.x0
        y0 = self.controller.parametros_actuales.y0
        
        # Eliminar la línea anterior si existe
        if hasattr(self, 'linea_trayectoria'):
            self.linea_trayectoria.remove()
        
        # Crear una nueva línea con flecha
        self.linea_trayectoria = FancyArrowPatch(
            (x0, y0), (x, y),
            arrowstyle='->', mutation_scale=15,
            color='green', linestyle='--', linewidth=1, alpha=0.5
        )
        self.ax.add_patch(self.linea_trayectoria)
        
        # Línea de dirección (vector velocidad)
        scale = 0.5  # Escala para visualización del vector
        if hasattr(self, 'linea_direccion'):
            self.linea_direccion.remove()
        
        self.linea_direccion = FancyArrowPatch(
            (x, y), (x + vx*scale, y + vy*scale),
            arrowstyle='->', mutation_scale=15,
            color='red', linewidth=2, alpha=0.7
        )
        self.ax.add_patch(self.linea_direccion)
        
        # Línea de aceleración
        if hasattr(self, 'linea_aceleracion'):
            self.linea_aceleracion.remove()
        
        a_scale = 0.1  # Escala diferente para la aceleración
        self.linea_aceleracion = FancyArrowPatch(
            (x, y), (x + ax*a_scale, y + ay*a_scale),
            arrowstyle='->', mutation_scale=15,
            color='blue', linewidth=2, alpha=0.7
        )
        self.ax.add_patch(self.linea_aceleracion)
        
        self.canvas.draw()
    
    def actualizar_punto_comparacion(self, x: float, y: float, t: float, 
                                   vx: float, vy: float, ax: float, ay: float):
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
        x0_comp = self.controller.parametros_comparacion.x0
        y0_comp = self.controller.parametros_comparacion.y0
        
        if hasattr(self, 'linea_trayectoria_comp'):
            self.linea_trayectoria_comp.remove()
        
        self.linea_trayectoria_comp = FancyArrowPatch(
            (x0_comp, y0_comp), (x, y),
            arrowstyle='->', mutation_scale=15,
            color='magenta', linestyle='--', linewidth=1, alpha=0.5
        )
        self.ax.add_patch(self.linea_trayectoria_comp)
        
        # Línea de dirección (vector velocidad) para comparación
        scale = 0.5
        if hasattr(self, 'linea_direccion_comp'):
            self.linea_direccion_comp.remove()
        
        self.linea_direccion_comp = FancyArrowPatch(
            (x, y), (x + vx*scale, y + vy*scale),
            arrowstyle='->', mutation_scale=15,
            color='cyan', linewidth=2, alpha=0.7
        )
        self.ax.add_patch(self.linea_direccion_comp)
        
        self.canvas.draw()


class PanelHistorial(PanelBase):
    """Panel para mostrar el historial de simulaciones"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
    
    def crear_widgets(self):
        # Frame principal con scroll
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True)
        
        # Create scrollable frame
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
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
        
        ttk.Button(frame_botones, text="Cargar Simulación", command=self.controller.cargar_desde_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Eliminar Selección", command=self.controller.eliminar_del_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Limpiar Historial", command=self.controller.limpiar_historial).pack(side="right", padx=5)
    
    def obtener_datos(self) -> Dict:
        return {}
    
    def establecer_datos(self, datos: Dict):
        pass
    
    def actualizar_lista_historial(self, historial: List[Dict]):
        self.historial_listbox.delete(0, tk.END)
        
        for entrada in historial:
            timestamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(entrada['timestamp']))
            datos = entrada['datos']
            res = entrada['resultados']
            
            texto = f"{timestamp} - v0={datos['v0']:.1f} m/s, ángulo={datos['angulo']:.1f}° | "
            texto += f"Alcance: {res['alcance']:.1f}m, Tiempo: {res['t_vuelo']:.1f}s"
            
            self.historial_listbox.insert(tk.END, texto)


class BarraHerramientas(ttk.Frame):
    """Barra de herramientas principal de la aplicación"""
    
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        self.crear_widgets()
        
        # Configurar para que no cambie de tamaño
        self.grid_propagate(False)
        self.config(height=40)  # Altura fija
    
    def crear_widgets(self):
        self.style = ttk.Style()
        self.style.configure('Accent.TButton', font=('Helvetica', 10, 'bold'), foreground='blue')
        
        # Frame para centrar los botones
        btn_frame = ttk.Frame(self)
        btn_frame.pack(expand=True)
        
        ttk.Button(btn_frame, text="Calcular", command=self.controller.calcular_trayectoria, 
                  style='Accent.TButton').pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Exportar Informe", command=self.controller.exportar_informe).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Limpiar", command=self.controller.limpiar).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="Modo Avanzado", command=self.controller.modo_avanzado).pack(side='left', padx=5)
        
        # Añadir fecha/hora como en la imagen
        fecha_hora = time.strftime('%d/%m/%Y %H:%M')
        lbl_fecha = ttk.Label(self, text=fecha_hora)
        lbl_fecha.pack(side='right', padx=10)


class TiroParabolicoApp:
    """Clase principal de la aplicación"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Simulador de Tiro Parabólico Avanzado")
        self.root.geometry("1100x750")
        
        # Modelo
        self.modelo = ModeloTiroParabolico()
        self.parametros_actuales = ParametrosTiro()
        self.parametros_comparacion = ParametrosTiro(v0=15, angulo=30)
        self.resultados_actuales: Optional[ResultadosTiro] = None
        self.resultados_comparacion: Optional[ResultadosTiro] = None
        self.historial: List[Dict] = []
        self.historial_file = "historial_simulaciones.json"
        
        # Cargar historial
        self.cargar_historial()
        
        # Interfaz
        self.crear_interfaz()
    
    def crear_interfaz(self):
        # Frame principal con configuración de grid
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill='both', expand=True)
        
        # Configurar grid
        main_frame.grid_rowconfigure(0, weight=1)  # Fila para el notebook (expande)
        main_frame.grid_rowconfigure(1, weight=0)  # Fila para la barra de herramientas (tamaño fijo)
        main_frame.grid_columnconfigure(0, weight=1)  # Columna única
        
        # Create notebook (tabs) - en la fila 0
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=0, column=0, sticky='nsew', padx=5, pady=5)
        
        # Create tabs
        self.crear_pestana_entrada(notebook)
        self.crear_pestana_resultados(notebook)
        self.crear_pestana_graficos(notebook)
        self.crear_pestana_historial(notebook)
        
        # Create toolbar - en la fila 1
        self.barra_herramientas = BarraHerramientas(main_frame, self)
        self.barra_herramientas.grid(row=1, column=0, sticky='ew', padx=5, pady=(0,5))
    
    def crear_pestana_entrada(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Datos de Entrada")
        
        # Create panels
        self.panel_basico = PanelParametrosBasicos(frame, self)
        self.panel_basico.pack(fill='both', expand=True)
        
        self.panel_resistencia = PanelResistenciaAire(frame, self)
        self.panel_resistencia.pack(fill='both', expand=True)
        
        self.panel_comparacion = PanelComparacion(frame, self)
        self.panel_comparacion.pack(fill='both', expand=True)
        
        self.panel_unidades = PanelUnidades(frame, self)
        self.panel_unidades.pack(fill='both', expand=True)
    
    def crear_pestana_resultados(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Resultados")
        
        self.panel_resultados = PanelResultados(frame, self)
        self.panel_resultados.pack(fill='both', expand=True)
    
    def crear_pestana_graficos(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Gráficos")
        
        self.panel_graficos = PanelGraficos(frame, self)
        self.panel_graficos.pack(fill='both', expand=True)
    
    def crear_pestana_historial(self, notebook):
        frame = ttk.Frame(notebook)
        notebook.add(frame, text="Historial")
        
        self.panel_historial = PanelHistorial(frame, self)
        self.panel_historial.pack(fill='both', expand=True)
        self.panel_historial.actualizar_lista_historial(self.historial)
    
    def toggle_resistencia(self):
        if self.panel_basico.vars['resistencia'].get():
            self.panel_resistencia.pack(fill='both', expand=True)
        else:
            self.panel_resistencia.pack_forget()
    
    def toggle_comparacion(self):
        if self.panel_comparacion.vars['comparar'].get():
            self.panel_comparacion.frame_params.pack(fill='x', padx=5, pady=5)
        else:
            self.panel_comparacion.frame_params.pack_forget()
    
    def toggle_animacion(self):
        """Alternar entre iniciar y detener la animación"""
        if not hasattr(self, 'animacion') or self.animacion is None:
            self.iniciar_animacion()
        else:
            self.detener_animacion()
    
    def iniciar_animacion(self):
        """Iniciar la animación de la trayectoria"""
        if not self.resultados_actuales:
            messagebox.showwarning("Advertencia", "No hay datos para animar. Calcule primero una trayectoria.")
            return
            
        # Animation setup
        self.linea_animada, = self.panel_graficos.ax.plot([], [], 'ko', markersize=8)
        self.tiempo_texto = self.panel_graficos.ax.text(0.02, 0.95, '', transform=self.panel_graficos.ax.transAxes)
        
        def init():
            self.linea_animada.set_data([], [])
            self.tiempo_texto.set_text('')
            return self.linea_animada, self.tiempo_texto
        
        def animate(i):
            if i >= len(self.resultados_actuales.x):
                i = len(self.resultados_actuales.x) - 1
                
            x = self.resultados_actuales.x[i]
            y = self.resultados_actuales.y[i]
            t = self.resultados_actuales.tiempos[i]
            
            self.linea_animada.set_data(x, y)
            self.tiempo_texto.set_text(f'Tiempo: {t:.2f} s')
            
            # Update slider and labels
            self.panel_graficos.slider_trayectoria.set(i / len(self.resultados_actuales.x) * 100)
            self.actualizar_punto_trayectoria(i / len(self.resultados_actuales.x) * 100)
            
            return self.linea_animada, self.tiempo_texto
        
        # Calculate interval based on speed slider
        intervalo = 110 - self.panel_graficos.velocidad_animacion.get()
        
        # Create animation
        self.animacion = animation.FuncAnimation(
            self.panel_graficos.figura, animate, frames=len(self.resultados_actuales.x)*2,
            init_func=init, blit=True, interval=intervalo, repeat=True
        )
        
        self.panel_graficos.canvas.draw()
    
    def detener_animacion(self):
        """Detener la animación en curso"""
        if hasattr(self, 'animacion') and self.animacion:
            self.animacion.event_source.stop()
            self.animacion = None
            
            # Clean animation elements
            if hasattr(self, 'linea_animada'):
                self.linea_animada.remove()
                del self.linea_animada
            if hasattr(self, 'tiempo_texto'):
                self.tiempo_texto.remove()
                del self.tiempo_texto
            
            self.panel_graficos.canvas.draw()
    
    def actualizar_unidades(self):
        sistema = self.panel_unidades.var_unidades.get()
        if sistema == "métrico":
            pass  # Convert from imperial to metric if needed
        else:
            pass  # Convert from metric to imperial if needed
    
    def calcular_trayectoria(self):
        # Obtener parámetros actuales
        self.parametros_actuales = ParametrosTiro(
            **self.panel_basico.obtener_datos(),
            **self.panel_resistencia.obtener_datos()
        )
        
        # Validar parámetros
        if not self.validar_parametros(self.parametros_actuales):
            return
        
        # Calcular trayectoria principal
        self.resultados_actuales = self.modelo.calcular_trayectoria(self.parametros_actuales)
        
        # Calcular trayectoria de comparación si está activada
        if self.panel_comparacion.vars['comparar'].get():
            self.parametros_comparacion = ParametrosTiro(
                **{k: v for k, v in self.panel_comparacion.obtener_datos().items() if k != 'comparar'},
                g=self.parametros_actuales.g  # Usar misma gravedad que la trayectoria principal
            )
            
            if not self.validar_parametros(self.parametros_comparacion, es_comparacion=True):
                return
            
            self.resultados_comparacion = self.modelo.calcular_trayectoria(self.parametros_comparacion)
        
        self.actualizar_interfaz()
        self.guardar_en_historial()
    
    def validar_parametros(self, parametros: ParametrosTiro, es_comparacion: bool = False) -> bool:
        try:
            if parametros.v0 <= 0:
                raise ValueError("La velocidad inicial debe ser positiva")
                
            if not (0 <= parametros.angulo <= 90):
                raise ValueError("El ángulo debe estar entre 0° y 90°")
                
            if parametros.resistencia:
                if parametros.masa <= 0:
                    raise ValueError("La masa debe ser positiva")
                if parametros.coef_arrastre <= 0:
                    raise ValueError("El coeficiente de arrastre debe ser positivo")
                if parametros.area <= 0:
                    raise ValueError("El área transversal debe ser positiva")
            
            return True
        except ValueError as e:
            mensaje = f"Error en {'comparación' if es_comparacion else 'datos principales'}: {str(e)}"
            messagebox.showerror("Error de validación", mensaje)
            return False
    
    def actualizar_interfaz(self):
        # Actualizar resultados textuales
        informe = self.generar_informe()
        self.panel_resultados.actualizar_resultados(informe)
        
        # Actualizar energías
        texto_energia = self.generar_informe_energias()
        self.panel_resultados.actualizar_energias(texto_energia)
        
        # Actualizar gráfico
        if self.resultados_comparacion and self.panel_comparacion.vars['comparar'].get():
            self.panel_graficos.actualizar_grafico(
                self.resultados_actuales.x, self.resultados_actuales.y, self.resultados_actuales.alcance,
                self.resultados_comparacion.x, self.resultados_comparacion.y, self.resultados_comparacion.alcance
            )
        else:
            self.panel_graficos.actualizar_grafico(
                self.resultados_actuales.x, self.resultados_actuales.y, self.resultados_actuales.alcance
            )
    
    def generar_informe(self) -> str:
        informe = f"""
=== INFORME DE TIRO PARABÓLICO ===
- Posición inicial: ({self.parametros_actuales.x0:.2f}, {self.parametros_actuales.y0:.2f}) m
- Velocidad inicial: {self.parametros_actuales.v0:.2f} m/s
- Ángulo: {self.parametros_actuales.angulo:.2f}°
- Gravedad: {self.parametros_actuales.g:.2f} m/s²
- {'CON resistencia del aire' if self.parametros_actuales.resistencia else 'SIN resistencia del aire'}

=== RESULTADOS ===
- Tiempo de vuelo: {self.resultados_actuales.t_vuelo:.2f} s
- Alcance máximo: {self.resultados_actuales.alcance:.2f} m
- Altura máxima: {self.resultados_actuales.altura_max:.2f} m
- Velocidad final (impacto): {np.sqrt(self.resultados_actuales.vx_final**2 + self.resultados_actuales.vy_final**2):.2f} m/s
"""
        if self.parametros_actuales.resistencia:
            informe += f"- Coeficiente de arrastre: {self.parametros_actuales.coef_arrastre:.4f}\n"
        
        if self.resultados_comparacion and self.panel_comparacion.vars['comparar'].get():
            informe += "\n=== COMPARACIÓN ===\n"
            informe += f"- Tiempo de vuelo: {self.resultados_comparacion.t_vuelo:.2f} s (vs {self.resultados_actuales.t_vuelo:.2f} s)\n"
            informe += f"- Alcance máximo: {self.resultados_comparacion.alcance:.2f} m (vs {self.resultados_actuales.alcance:.2f} m)\n"
            informe += f"- Altura máxima: {self.resultados_comparacion.altura_max:.2f} m (vs {self.resultados_actuales.altura_max:.2f} m)\n"
        
        return informe
    
    def generar_informe_energias(self) -> str:
        if not self.resultados_actuales or not hasattr(self.resultados_actuales, 'energias'):
            return ""
            
        energias = self.resultados_actuales.energias
        texto_energia = "=== ENERGÍAS ===\n"
        texto_energia += f"- Energía cinética inicial: {energias['cinetica'][0]:.2f} J\n"
        texto_energia += f"- Energía potencial inicial: {energias['potencial'][0]:.2f} J\n"
        texto_energia += f"- Energía total inicial: {energias['total'][0]:.2f} J\n\n"
        
        texto_energia += f"- Energía cinética final: {energias['cinetica'][-1]:.2f} J\n"
        texto_energia += f"- Energía potencial final: {energias['potencial'][-1]:.2f} J\n"
        texto_energia += f"- Energía total final: {energias['total'][-1]:.2f} J\n"
        
        if self.parametros_actuales.resistencia:
            perdida = (energias['total'][0] - energias['total'][-1]) / energias['total'][0] * 100
            texto_energia += f"\n- Pérdida de energía por resistencia: {perdida:.2f}%"
        
        return texto_energia
    
    def actualizar_punto_trayectoria(self, valor):
        if not self.resultados_actuales:
            return
        
        idx = int(float(valor) / 100 * (len(self.resultados_actuales.x) - 1))
        idx = max(0, min(idx, len(self.resultados_actuales.x) - 1))
        
        x = self.resultados_actuales.x[idx]
        y = self.resultados_actuales.y[idx]
        t = self.resultados_actuales.tiempos[idx]
        
        if idx == 0:
            vx = self.resultados_actuales.v0x
            vy = self.resultados_actuales.v0y
            ax = 0
            ay = -self.parametros_actuales.g
        elif idx == 1:
            dt = self.resultados_actuales.tiempos[idx] - self.resultados_actuales.tiempos[idx-1]
            vx_prev = self.resultados_actuales.v0x
            vy_prev = self.resultados_actuales.v0y
            vx = (self.resultados_actuales.x[idx] - self.resultados_actuales.x[idx-1]) / dt
            vy = (self.resultados_actuales.y[idx] - self.resultados_actuales.y[idx-1]) / dt
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        else:
            dt = self.resultados_actuales.tiempos[idx] - self.resultados_actuales.tiempos[idx-1]
            dt_prev = self.resultados_actuales.tiempos[idx-1] - self.resultados_actuales.tiempos[idx-2]
            vx = (self.resultados_actuales.x[idx] - self.resultados_actuales.x[idx-1]) / dt
            vy = (self.resultados_actuales.y[idx] - self.resultados_actuales.y[idx-1]) / dt
            vx_prev = (self.resultados_actuales.x[idx-1] - self.resultados_actuales.x[idx-2]) / dt_prev
            vy_prev = (self.resultados_actuales.y[idx-1] - self.resultados_actuales.y[idx-2]) / dt_prev
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        
        self.panel_graficos.actualizar_punto_trayectoria(x, y, t, vx, vy, ax, ay)
    
    def actualizar_punto_comparacion(self, valor):
        if not self.resultados_comparacion:
            return
        
        idx = int(float(valor) / 100 * (len(self.resultados_comparacion.x) - 1))
        idx = max(0, min(idx, len(self.resultados_comparacion.x) - 1))
        
        x = self.resultados_comparacion.x[idx]
        y = self.resultados_comparacion.y[idx]
        t = self.resultados_comparacion.tiempos[idx]
        
        if idx == 0:
            vx = self.resultados_comparacion.v0x
            vy = self.resultados_comparacion.v0y
            ax = 0
            ay = -self.parametros_comparacion.g
        elif idx == 1:
            dt = self.resultados_comparacion.tiempos[idx] - self.resultados_comparacion.tiempos[idx-1]
            vx_prev = self.resultados_comparacion.v0x
            vy_prev = self.resultados_comparacion.v0y
            vx = (self.resultados_comparacion.x[idx] - self.resultados_comparacion.x[idx-1]) / dt
            vy = (self.resultados_comparacion.y[idx] - self.resultados_comparacion.y[idx-1]) / dt
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        else:
            dt = self.resultados_comparacion.tiempos[idx] - self.resultados_comparacion.tiempos[idx-1]
            dt_prev = self.resultados_comparacion.tiempos[idx-1] - self.resultados_comparacion.tiempos[idx-2]
            vx = (self.resultados_comparacion.x[idx] - self.resultados_comparacion.x[idx-1]) / dt
            vy = (self.resultados_comparacion.y[idx] - self.resultados_comparacion.y[idx-1]) / dt
            vx_prev = (self.resultados_comparacion.x[idx-1] - self.resultados_comparacion.x[idx-2]) / dt_prev
            vy_prev = (self.resultados_comparacion.y[idx-1] - self.resultados_comparacion.y[idx-2]) / dt_prev
            ax = (vx - vx_prev) / dt
            ay = (vy - vy_prev) / dt
        
        self.panel_graficos.actualizar_punto_comparacion(x, y, t, vx, vy, ax, ay)
    
    def consultar_tiempo(self):
        try:
            t = float(self.panel_resultados.tiempo_consulta.get())
            if t < 0:
                raise ValueError("El tiempo no puede ser negativo")
                
            param = parametrizacion_tiro_parabolico(
                self.parametros_actuales.x0, self.parametros_actuales.y0, 
                self.parametros_actuales.v0, self.parametros_actuales.angulo, 
                self.parametros_actuales.g, t
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
            
            self.panel_resultados.mostrar_consulta(resultado)
            
        except ValueError as e:
            messagebox.showerror("Error", str(e))
    
    def exportar_informe(self):
        if not self.resultados_actuales:
            messagebox.showwarning("Advertencia", "No hay resultados para exportar")
            return
            
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Archivos de texto", "*.txt"), ("Archivos PDF", "*.pdf"), ("Todos los archivos", "*.*")],
            title="Guardar informe como"
        )
        
        if filepath:
            try:
                contenido = self.generar_informe() + "\n" + self.generar_informe_energias()
                
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
        if not hasattr(self.panel_graficos, 'figura'):
            messagebox.showwarning("Advertencia", "No hay gráfico para exportar")
            return
            
        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("JPEG", "*.jpg"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("Todos los archivos", "*.*")],
            title="Guardar gráfico como"
        )
        
        if filepath:
            try:
                self.panel_graficos.figura.savefig(filepath, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Éxito", f"Gráfico guardado en:\n{filepath}")
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo guardar el gráfico:\n{str(e)}")
    
    def limpiar(self):
        # Reset panels
        self.panel_basico.establecer_datos({
            'x0': 0, 'y0': 0, 'v0': 10, 'angulo': 45, 'g': 9.81, 'resistencia': False
        })
        self.panel_resistencia.establecer_datos({
            'masa': 1, 'coef_arrastre': 0.47, 'area': 0.1, 'densidad_aire': 1.225
        })
        self.panel_comparacion.establecer_datos({
            'comparar': False, 'x0': 0, 'y0': 0, 'v0': 15, 'angulo': 30, 'resistencia': False
        })
        self.panel_unidades.establecer_datos({'unidades': "métrico"})
        
        # Clear results
        self.panel_resultados.actualizar_resultados("")
        self.panel_resultados.actualizar_energias("")
        self.panel_resultados.mostrar_consulta("")
        self.panel_resultados.tiempo_consulta.delete(0, tk.END)
        
        # Clear plot
        if hasattr(self.panel_graficos, 'ax'):
            self.panel_graficos.ax.clear()
            # Dibujar ejes vacíos para mantener los vectores visibles
            self.panel_graficos.ax.set_title("Trayectoria del Proyectil")
            self.panel_graficos.ax.set_xlabel("Distancia (m)")
            self.panel_graficos.ax.set_ylabel("Altura (m)")
            self.panel_graficos.ax.grid(True)
            
            # Eliminar cualquier elemento de vector existente
            if hasattr(self.panel_graficos, 'punto_interactivo'):
                self.panel_graficos.punto_interactivo.remove()
                del self.panel_graficos.punto_interactivo
            if hasattr(self.panel_graficos, 'linea_trayectoria'):
                self.panel_graficos.linea_trayectoria.remove()
                del self.panel_graficos.linea_trayectoria
            if hasattr(self.panel_graficos, 'linea_direccion'):
                self.panel_graficos.linea_direccion.remove()
                del self.panel_graficos.linea_direccion
            if hasattr(self.panel_graficos, 'linea_aceleracion'):
                self.panel_graficos.linea_aceleracion.remove()
                del self.panel_graficos.linea_aceleracion
            
            self.panel_graficos.canvas.draw()
        
        # Stop animation
        if hasattr(self, 'animacion') and self.animacion:
            self.detener_animacion()
        
        # Reset variables
        self.resultados_actuales = None
        self.resultados_comparacion = None
    
    def modo_avanzado(self):
        messagebox.showinfo("Información", "Modo avanzado no implementado en esta versión")
    
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
    
    def guardar_en_historial(self):
        entrada = {
            'timestamp': time.time(),
            'datos': {
                'x0': self.parametros_actuales.x0,
                'y0': self.parametros_actuales.y0,
                'v0': self.parametros_actuales.v0,
                'angulo': self.parametros_actuales.angulo,
                'g': self.parametros_actuales.g,
                'resistencia': self.parametros_actuales.resistencia,
                'masa': self.parametros_actuales.masa,
                'coef_arrastre': self.parametros_actuales.coef_arrastre,
                'area': self.parametros_actuales.area,
                'densidad_aire': self.parametros_actuales.densidad_aire
            },
            'resultados': {
                't_vuelo': self.resultados_actuales.t_vuelo,
                'alcance': self.resultados_actuales.alcance,
                'altura_max': self.resultados_actuales.altura_max
            }
        }
        
        self.historial.append(entrada)
        
        if len(self.historial) > 50:
            self.historial = self.historial[-50:]
        
        self.guardar_historial()
        self.panel_historial.actualizar_lista_historial(self.historial)
    
    def cargar_desde_historial(self):
        seleccion = self.panel_historial.historial_listbox.curselection()
        if not seleccion:
            messagebox.showwarning("Advertencia", "Seleccione una simulación del historial")
            return
            
        indice = seleccion[0]
        entrada = self.historial[indice]
        datos = entrada['datos']
        
        # Establecer datos en los paneles
        self.panel_basico.establecer_datos({
            'x0': datos['x0'],
            'y0': datos['y0'],
            'v0': datos['v0'],
            'angulo': datos['angulo'],
            'g': datos['g'],
            'resistencia': datos['resistencia']
        })
        
        if datos['resistencia']:
            self.panel_resistencia.establecer_datos({
                'masa': datos['masa'],
                'coef_arrastre': datos['coef_arrastre'],
                'area': datos['area'],
                'densidad_aire': datos['densidad_aire']
            })
        
        # Recalcular
        self.calcular_trayectoria()
    
    def eliminar_del_historial(self):
        seleccion = self.panel_historial.historial_listbox.curselection()
        if not seleccion:
            messagebox.showwarning("Advertencia", "Seleccione una simulación del historial")
            return
            
        indice = seleccion[0]
        del self.historial[indice]
        self.guardar_historial()
        self.panel_historial.actualizar_lista_historial(self.historial)
    
    def limpiar_historial(self):
        if messagebox.askyesno("Confirmar", "¿Está seguro que desea borrar todo el historial?"):
            self.historial = []
            self.guardar_historial()
            self.panel_historial.actualizar_lista_historial(self.historial)


def parametrizacion_tiro_parabolico(x0: float, y0: float, v0: float, angulo: float, 
                                   g: float, t: float) -> Dict:
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