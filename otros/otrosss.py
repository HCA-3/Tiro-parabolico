# -*- coding: utf-8 -*-
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
from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple, Optional
from matplotlib.patches import FancyArrowPatch

# ==============================================
# Clases para el modelo físico (Dominio)
# ==============================================

@dataclass
class ParametrosTiro:
    """Clase para almacenar los parámetros de un tiro parabólico"""
    x0: float = 0.0
    y0: float = 0.0
    v0: float = 10.0
    angulo: float = 45.0  # en grados
    g: float = 9.81
    resistencia: bool = False
    masa: float = 1.0
    coef_arrastre: float = 0.47
    area: float = 0.1
    densidad_aire: float = 1.225
    comparar: bool = False  # Para uso en la interfaz, no afecta los cálculos físicos


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
    def descomposicion_vectorial(v0: float, theta: float) -> Tuple[float, float]:
        """Descompone la velocidad inicial en componentes x e y."""
        return v0 * np.cos(theta), v0 * np.sin(theta)

    @staticmethod
    def tiempo_vuelo(v0y: float, g: float, y0: float = 0.0) -> float:
        """Calcula el tiempo de vuelo para un tiro sin resistencia."""
        # Resuelve y0 + v0y*t - 0.5*g*t^2 = 0 para t > 0
        discriminante = v0y**2 + 2 * g * y0
        if discriminante < 0:
            # Esto no debería ocurrir si y0 >= 0
            return 0
        t = (v0y + np.sqrt(discriminante)) / g
        return t if t > 0 else 0

    @staticmethod
    def alcance_maximo(v0x: float, t_vuelo: float, x0: float = 0.0) -> float:
        """Calcula el alcance máximo para un tiro sin resistencia."""
        return x0 + v0x * t_vuelo

    @staticmethod
    def altura_maxima(v0y: float, g: float, y0: float = 0.0) -> float:
        """Calcula la altura máxima alcanzada para un tiro sin resistencia."""
        # Tiempo para alcanzar la altura máxima (vy = 0): t_hmax = v0y / g
        # h_max = y0 + v0y * t_hmax - 0.5 * g * t_hmax**2
        return y0 + (v0y**2) / (2 * g)

    @staticmethod
    def movimiento_con_resistencia(
        x0: float, y0: float, v0x: float, v0y: float,
        g: float, masa: float, coef_arrastre: float,
        area: float, densidad_aire: float,
        dt: float = 0.01
    ) -> Tuple[List[float], List[float], List[float]]:
        """Calcula la trayectoria usando el método de Euler con resistencia del aire."""
        x, y = [x0], [y0]
        vx, vy = v0x, v0y
        t = 0.0
        tiempos = [t]

        while y[-1] >= 0 or len(y) < 2: # Continúa mientras esté sobre el suelo o al inicio
            velocidad = np.sqrt(vx**2 + vy**2)
            if velocidad == 0: # Evitar división por cero si la velocidad es nula
                ax = 0
                ay = -g
            else:
                fuerza_arrastre = 0.5 * densidad_aire * coef_arrastre * area * velocidad**2
                ax = - (fuerza_arrastre / masa) * (vx / velocidad)
                ay = -g - (fuerza_arrastre / masa) * (vy / velocidad)

            vx_new = vx + ax * dt
            vy_new = vy + ay * dt
            x_new = x[-1] + vx * dt # Usar velocidad *anterior* para posición (Euler semi-implícito)
            y_new = y[-1] + vy * dt

            # Corrección si cruza y=0
            if y_new < 0 and y[-1] >= 0:
                # Interpolación lineal para encontrar el tiempo exacto de impacto
                t_impacto = t - y[-1] * dt / (y_new - y[-1])
                dt_final = t_impacto - t
                
                # Recalcular último paso con dt_final
                vx_impacto = vx + ax * dt_final
                vy_impacto = vy + ay * dt_final
                x_impacto = x[-1] + vx * dt_final
                y_impacto = 0.0 # Por definición
                
                x.append(x_impacto)
                y.append(y_impacto)
                tiempos.append(t_impacto)
                vx = vx_impacto
                vy = vy_impacto
                break # Termina la simulación
            
            # Actualizar valores para el siguiente paso
            vx, vy = vx_new, vy_new
            x.append(x_new)
            y.append(y_new)
            t += dt
            tiempos.append(t)
            
            # Condición de parada de seguridad (si vuela demasiado tiempo)
            if t > 1000: # Límite de tiempo arbitrario
                print("Advertencia: Simulación detenida por exceder tiempo límite.")
                break
                
        return x, y, tiempos

    @staticmethod
    def calcular_energias(
        parametros: ParametrosTiro, x: List[float], y: List[float], tiempos: List[float]
    ) -> Dict[str, List[float]]:
        """Calcula las energías cinética, potencial y total en cada punto."""
        energias = {
            'cinetica': [],
            'potencial': [],
            'total': []
        }
        g = parametros.g
        masa = parametros.masa
        
        for i in range(len(x)):
            # Calcular velocidad en el punto i
            if len(x) == 1:
                vx = 0
                vy = 0
            elif i == 0:
                # Usar velocidad inicial descompuesta
                angulo_rad = np.radians(parametros.angulo)
                vx, vy = ModeloTiroParabolico.descomposicion_vectorial(parametros.v0, angulo_rad)
            else:
                dt = tiempos[i] - tiempos[i-1]
                if dt == 0: # Evitar división por cero si dt es muy pequeño
                     vx = energias['vx'][-1] if 'vx' in energias and energias['vx'] else 0
                     vy = energias['vy'][-1] if 'vy' in energias and energias['vy'] else 0
                else:
                    vx = (x[i] - x[i-1]) / dt
                    vy = (y[i] - y[i-1]) / dt
            
            velocidad_sq = vx**2 + vy**2
            altura = y[i]
            
            ec = 0.5 * masa * velocidad_sq
            # Asegurar que la altura no sea negativa para el cálculo de energía potencial
            ep = masa * g * max(altura, 0)
            et = ec + ep
            
            energias['cinetica'].append(ec)
            energias['potencial'].append(ep)
            energias['total'].append(et)
            # Guardar velocidades para referencia (opcional)
            if 'vx' not in energias: energias['vx'] = []
            if 'vy' not in energias: energias['vy'] = []
            energias['vx'].append(vx)
            energias['vy'].append(vy)
            
        return energias

    @staticmethod
    def calcular_trayectoria(parametros: ParametrosTiro) -> ResultadosTiro:
        """Calcula la trayectoria completa basada en los parámetros dados."""
        angulo_rad = np.radians(parametros.angulo)
        v0x, v0y = ModeloTiroParabolico.descomposicion_vectorial(parametros.v0, angulo_rad)
        
        if parametros.resistencia:
            x, y, tiempos = ModeloTiroParabolico.movimiento_con_resistencia(
                parametros.x0, parametros.y0, v0x, v0y,
                parametros.g, parametros.masa, parametros.coef_arrastre,
                parametros.area, parametros.densidad_aire
            )
            if not tiempos: # Si no hubo movimiento
                t_vuelo = 0.0
                alcance = parametros.x0
                altura_max = parametros.y0
                vx_final = 0.0
                vy_final = 0.0
            else:
                t_vuelo = tiempos[-1]
                alcance = x[-1]
                altura_max = max(y) if y else parametros.y0
                # Calcular velocidad final (aproximada si dt es pequeño)
                if len(tiempos) > 1:
                    dt_final = tiempos[-1] - tiempos[-2]
                    vx_final = (x[-1] - x[-2]) / dt_final if dt_final > 0 else 0
                    vy_final = (y[-1] - y[-2]) / dt_final if dt_final > 0 else 0
                else:
                    vx_final = v0x
                    vy_final = v0y
        else:
            t_vuelo = ModeloTiroParabolico.tiempo_vuelo(v0y, parametros.g, parametros.y0)
            tiempos = np.linspace(0, t_vuelo, 200) # Más puntos para suavidad
            x = parametros.x0 + v0x * tiempos
            y = parametros.y0 + v0y * tiempos - 0.5 * parametros.g * tiempos**2
            # Asegurar que el último punto esté en y=0
            x = list(x)
            y = list(y)
            tiempos = list(tiempos)
            if y[-1] != 0:
                 x.append(parametros.x0 + v0x * t_vuelo)
                 y.append(0.0)
                 tiempos.append(t_vuelo)
                 
            alcance = ModeloTiroParabolico.alcance_maximo(v0x, t_vuelo, parametros.x0)
            altura_max = ModeloTiroParabolico.altura_maxima(v0y, parametros.g, parametros.y0)
            vx_final = v0x
            vy_final = v0y - parametros.g * t_vuelo

        # Asegurar que y no sea negativo debido a errores de precisión
        y = [max(0, val) for val in y]
        
        energias = ModeloTiroParabolico.calcular_energias(parametros, x, y, tiempos)
        
        return ResultadosTiro(
            x=x, y=y, tiempos=tiempos,
            v0x=v0x, v0y=v0y,
            t_vuelo=t_vuelo, alcance=alcance, altura_max=altura_max,
            vx_final=vx_final, vy_final=vy_final,
            energias=energias
        )

    @staticmethod
    def parametrizacion_tiro_sin_resistencia(x0, y0, v0, angulo, g, t):
        """Calcula posición, velocidad y aceleración en un tiempo t (sin resistencia)."""
        angulo_rad = np.radians(angulo)
        v0x = v0 * np.cos(angulo_rad)
        v0y = v0 * np.sin(angulo_rad)
        
        xt = x0 + v0x * t
        yt = y0 + v0y * t - 0.5 * g * t**2
        
        vxt = v0x
        vyt = v0y - g * t
        
        axt = 0
        ayt = -g
        
        ecuaciones = {
            'x(t)': f"{x0:.2f} + {v0x:.2f}*t",
            'y(t)': f"{y0:.2f} + {v0y:.2f}*t - {0.5*g:.2f}*t^2",
            'vx(t)': f"{v0x:.2f}",
            'vy(t)': f"{v0y:.2f} - {g:.2f}*t"
        }
        
        return {
            'posicion': (xt, yt),
            'velocidad': (vxt, vyt),
            'aceleracion': (axt, ayt),
            'ecuaciones': ecuaciones
        }

# ==============================================
# Clases para la interfaz de usuario (Vista)
# ==============================================

class ScrollableFrame(ttk.Frame):
    """Un frame con scrollbars vertical y opcionalmente horizontal."""
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        canvas = tk.Canvas(self)
        scrollbar_y = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = ttk.Frame(canvas) # Frame donde van los widgets

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar_y.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar_y.pack(side="right", fill="y")

        # Bind mouse wheel scrolling
        self.bind_all("<MouseWheel>", lambda event: self._on_mousewheel(event, canvas), add='+')
        self.bind_all("<Button-4>", lambda event: self._on_mousewheel_linux(event, canvas, -1), add='+') # Linux scroll up
        self.bind_all("<Button-5>", lambda event: self._on_mousewheel_linux(event, canvas, 1), add='+') # Linux scroll down

    def _on_mousewheel(self, event, canvas):
        # Cross-platform mouse wheel scroll
        if event.delta:
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")

    def _on_mousewheel_linux(self, event, canvas, direction):
        # Linux specific mouse wheel scroll
        canvas.yview_scroll(direction, "units")


class PanelBase(ttk.Frame, ABC):
    """Clase base abstracta para los paneles de la interfaz"""
    def __init__(self, parent, controller, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.controller = controller
        # Configurar expansión para adaptabilidad
        self.pack(fill='both', expand=True, padx=5, pady=5)
        self.crear_widgets()

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
    def crear_widgets(self):
        scrollable = ScrollableFrame(self)
        scrollable.pack(fill="both", expand=True)
        frame = ttk.LabelFrame(scrollable.scrollable_frame, text="Parámetros Básicos", padding=10)
        frame.pack(fill="x", expand=True, padx=10, pady=10)

        self.vars = {
            'x0': tk.DoubleVar(value=0.0),
            'y0': tk.DoubleVar(value=0.0),
            'v0': tk.DoubleVar(value=10.0),
            'angulo': tk.DoubleVar(value=45.0),
            'g': tk.DoubleVar(value=9.81),
            'resistencia': tk.BooleanVar(value=False)
        }

        entries = [
            ("Posición inicial x₀ (m):", 'x0'),
            ("Posición inicial y₀ (m):", 'y0'),
            ("Velocidad inicial |v₀| (m/s):", 'v0'),
            ("Ángulo de lanzamiento θ (°):", 'angulo'),
            ("Aceleración gravedad g (m/s²):", 'g')
        ]

        # Usar grid para mejor alineación y adaptabilidad
        frame.columnconfigure(1, weight=1) # Columna de entrada se expande
        for i, (text, var) in enumerate(entries):
            lbl = ttk.Label(frame, text=text)
            lbl.grid(row=i, column=0, sticky="w", padx=5, pady=2)
            entry = ttk.Entry(frame, textvariable=self.vars[var], width=15)
            entry.grid(row=i, column=1, sticky="ew", padx=5, pady=2)
            # Añadir validación (ejemplo simple)
            entry.config(validate="key", validatecommand=(self.register(self.validate_float), '%P'))

        chk_resistencia = ttk.Checkbutton(frame, text="Considerar resistencia del aire",
                       variable=self.vars['resistencia'],
                       command=self.controller.toggle_resistencia)
        chk_resistencia.grid(row=len(entries), column=0, columnspan=2, pady=10, sticky='w')

    def validate_float(self, value_if_allowed):
        """Valida que la entrada sea un número flotante."""
        if value_if_allowed == "" or value_if_allowed == "-":
            return True
        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False

    def obtener_datos(self) -> Dict:
        return {k: v.get() for k, v in self.vars.items()}

    def establecer_datos(self, datos: Dict):
        for k, v in datos.items():
            if k in self.vars:
                try:
                    self.vars[k].set(v)
                except tk.TclError:
                    # Manejar posible error si el valor no es compatible con el tipo de variable
                    print(f"Error al establecer {k}={v}")
                    pass # O establecer un valor por defecto


class PanelResistenciaAire(PanelBase):
    """Panel para los parámetros de resistencia del aire"""
    def crear_widgets(self):
        scrollable = ScrollableFrame(self)
        scrollable.pack(fill="both", expand=True)
        self.frame = ttk.LabelFrame(scrollable.scrollable_frame, text="Parámetros de Resistencia del Aire", padding=10)
        self.frame.pack(fill="x", expand=True, padx=10, pady=10)

        self.vars = {
            'masa': tk.DoubleVar(value=1.0),
            'coef_arrastre': tk.DoubleVar(value=0.47),
            'area': tk.DoubleVar(value=0.1),
            'densidad_aire': tk.DoubleVar(value=1.225)
        }

        entries = [
            ("Masa del proyectil (kg):", 'masa'),
            ("Coeficiente de arrastre C<sub>d</sub>:", 'coef_arrastre'),
            ("Área transversal A (m²):", 'area'),
            ("Densidad del aire ρ (kg/m³):", 'densidad_aire')
        ]

        self.frame.columnconfigure(1, weight=1)
        for i, (text, var) in enumerate(entries):
            lbl = ttk.Label(self.frame, text=text)
            lbl.grid(row=i, column=0, sticky="w", padx=5, pady=2)
            entry = ttk.Entry(self.frame, textvariable=self.vars[var], width=15)
            entry.grid(row=i, column=1, sticky="ew", padx=5, pady=2)
            entry.config(validate="key", validatecommand=(self.register(self.validate_float), '%P'))

        # Inicialmente deshabilitado si la resistencia no está activa
        self.set_enabled(self.controller.parametros_actuales.resistencia)

    def validate_float(self, value_if_allowed):
        if value_if_allowed == "" or value_if_allowed == "-":
            return True
        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False

    def obtener_datos(self) -> Dict:
        return {k: v.get() for k, v in self.vars.items()}

    def establecer_datos(self, datos: Dict):
        for k, v in datos.items():
            if k in self.vars:
                 try:
                    self.vars[k].set(v)
                 except tk.TclError:
                    print(f"Error al establecer {k}={v}")
                    pass

    def set_enabled(self, enabled: bool):
        """Habilita o deshabilita los widgets del panel."""
        state = tk.NORMAL if enabled else tk.DISABLED
        for child in self.frame.winfo_children():
            try:
                child.configure(state=state)
            except tk.TclError:
                pass # Algunos widgets como Labels no tienen 'state'


class PanelComparacion(PanelBase):
    """Panel para la comparación de tiros"""
    def crear_widgets(self):
        scrollable = ScrollableFrame(self)
        scrollable.pack(fill="both", expand=True)
        frame = ttk.LabelFrame(scrollable.scrollable_frame, text="Comparación de Tiros", padding=10)
        frame.pack(fill="x", expand=True, padx=10, pady=10)

        self.vars = {
            'comparar': tk.BooleanVar(value=False),
            'x0': tk.DoubleVar(value=0.0),
            'y0': tk.DoubleVar(value=0.0),
            'v0': tk.DoubleVar(value=15.0),
            'angulo': tk.DoubleVar(value=30.0),
            'resistencia': tk.BooleanVar(value=False)
            # No incluir 'g' aquí, se asume la misma que el tiro principal
        }

        chk_comparar = ttk.Checkbutton(frame, text="Activar comparación",
                       variable=self.vars['comparar'],
                       command=self.controller.toggle_comparacion)
        chk_comparar.pack(anchor="w", pady=5)

        self.frame_params = ttk.Frame(frame, padding=5)
        self.frame_params.pack(fill="x", expand=True)

        entries = [
            ("Posición inicial x₀ (m):", 'x0'),
            ("Posición inicial y₀ (m):", 'y0'),
            ("Velocidad inicial |v₀| (m/s):", 'v0'),
            ("Ángulo de lanzamiento θ (°):", 'angulo')
        ]

        self.frame_params.columnconfigure(1, weight=1)
        for i, (text, var) in enumerate(entries):
            lbl = ttk.Label(self.frame_params, text=text)
            lbl.grid(row=i, column=0, sticky="w", padx=5, pady=2)
            entry = ttk.Entry(self.frame_params, textvariable=self.vars[var], width=15)
            entry.grid(row=i, column=1, sticky="ew", padx=5, pady=2)
            entry.config(validate="key", validatecommand=(self.register(self.validate_float), '%P'))

        chk_resistencia_comp = ttk.Checkbutton(self.frame_params, text="Resistencia del aire (Comparación)",
                       variable=self.vars['resistencia'])
        chk_resistencia_comp.grid(row=len(entries), column=0, columnspan=2, pady=5, sticky='w')

        # Inicialmente deshabilitado
        self.set_enabled(False)

    def validate_float(self, value_if_allowed):
        if value_if_allowed == "" or value_if_allowed == "-":
            return True
        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False

    def obtener_datos(self) -> Dict:
        return {k: v.get() for k, v in self.vars.items()}

    def establecer_datos(self, datos: Dict):
        for k, v in datos.items():
            if k in self.vars:
                 try:
                    self.vars[k].set(v)
                 except tk.TclError:
                    print(f"Error al establecer {k}={v}")
                    pass

    def set_enabled(self, enabled: bool):
        """Habilita o deshabilita los widgets del panel de comparación."""
        state = tk.NORMAL if enabled else tk.DISABLED
        for child in self.frame_params.winfo_children():
            try:
                child.configure(state=state)
            except tk.TclError:
                pass


class PanelUnidades(PanelBase):
    """Panel para la configuración de unidades (Funcionalidad Futura)"""
    def crear_widgets(self):
        scrollable = ScrollableFrame(self)
        scrollable.pack(fill="both", expand=True)
        frame = ttk.LabelFrame(scrollable.scrollable_frame, text="Configuración de Unidades (Futuro)", padding=10)
        frame.pack(fill="x", expand=True, padx=10, pady=10)

        self.var_unidades = tk.StringVar(value="métrico")

        rb_metrico = ttk.Radiobutton(frame, text="Sistema Métrico (m, kg, s)",
                       variable=self.var_unidades, value="métrico",
                       command=self.controller.actualizar_unidades, state=tk.DISABLED)
        rb_metrico.pack(anchor="w")

        rb_imperial = ttk.Radiobutton(frame, text="Sistema Imperial (pies, lb, s)",
                       variable=self.var_unidades, value="imperial",
                       command=self.controller.actualizar_unidades, state=tk.DISABLED)
        rb_imperial.pack(anchor="w")
        
        ttk.Label(frame, text="(Funcionalidad no implementada)").pack(anchor="w", pady=5)

    def obtener_datos(self) -> Dict:
        return {'unidades': self.var_unidades.get()}

    def establecer_datos(self, datos: Dict):
        if 'unidades' in datos:
            self.var_unidades.set(datos['unidades'])


class PanelResultados(PanelBase):
    """Panel para mostrar los resultados del cálculo y consultas."""
    def crear_widgets(self):
        scrollable = ScrollableFrame(self)
        scrollable.pack(fill="both", expand=True)
        frame = ttk.LabelFrame(scrollable.scrollable_frame, text="Resultados y Consultas", padding=10)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Sección de Resultados Principales
        frame_resultados = ttk.Frame(frame)
        frame_resultados.pack(fill="x", pady=5)
        ttk.Label(frame_resultados, text="Resultados Principales:", font=("Helvetica", 10, "bold")).pack(anchor="w")
        self.resultados_text = tk.Text(frame_resultados, height=10, width=60, wrap=tk.WORD, state=tk.DISABLED)
        self.resultados_text.pack(fill="x", expand=True, pady=2)

        # Sección de Energías
        frame_energias = ttk.Frame(frame)
        frame_energias.pack(fill="x", pady=5)
        ttk.Label(frame_energias, text="Análisis de Energías:", font=("Helvetica", 10, "bold")).pack(anchor="w")
        self.energias_text = tk.Text(frame_energias, height=8, width=60, wrap=tk.WORD, state=tk.DISABLED)
        self.energias_text.pack(fill="x", expand=True, pady=2)

        # Sección de Consulta por Tiempo (Solo sin resistencia)
        frame_consulta = ttk.Frame(frame)
        frame_consulta.pack(fill="x", pady=10)
        ttk.Label(frame_consulta, text="Consultar Estado en Tiempo t (sin resistencia):", font=("Helvetica", 10, "bold")).pack(anchor="w")
        
        subframe_consulta = ttk.Frame(frame_consulta)
        subframe_consulta.pack(fill="x")
        ttk.Label(subframe_consulta, text="Tiempo t (s):").pack(side="left", padx=5)
        self.tiempo_consulta = ttk.Entry(subframe_consulta, width=10)
        self.tiempo_consulta.pack(side="left", padx=5)
        self.tiempo_consulta.config(validate="key", validatecommand=(self.register(self.validate_float), '%P'))
        self.btn_consultar = ttk.Button(subframe_consulta, text="Consultar", command=self.controller.consultar_tiempo)
        self.btn_consultar.pack(side="left", padx=5)
        
        self.consulta_text = tk.Text(frame_consulta, height=10, width=60, wrap=tk.WORD, state=tk.DISABLED)
        self.consulta_text.pack(fill="x", expand=True, pady=5)

    def validate_float(self, value_if_allowed):
        if value_if_allowed == "" or value_if_allowed == "-":
            return True
        try:
            # Permitir solo números positivos para el tiempo
            val = float(value_if_allowed)
            return val >= 0
        except ValueError:
            return False

    def _update_text_widget(self, widget: tk.Text, content: str):
        widget.config(state=tk.NORMAL)
        widget.delete('1.0', tk.END)
        widget.insert('1.0', content)
        widget.config(state=tk.DISABLED)

    def actualizar_resultados(self, informe: str):
        self._update_text_widget(self.resultados_text, informe)

    def actualizar_energias(self, informe_energias: str):
        self._update_text_widget(self.energias_text, informe_energias)

    def mostrar_consulta(self, resultado_consulta: str):
        self._update_text_widget(self.consulta_text, resultado_consulta)
        
    def set_consulta_enabled(self, enabled: bool):
        """Habilita o deshabilita la sección de consulta."""
        state = tk.NORMAL if enabled else tk.DISABLED
        self.tiempo_consulta.config(state=state)
        self.btn_consultar.config(state=state)
        if not enabled:
            self.mostrar_consulta("(Consulta deshabilitada con resistencia del aire)")
        else:
             self.mostrar_consulta("") # Limpiar si se habilita

    def obtener_datos(self) -> Dict:
        return {}

    def establecer_datos(self, datos: Dict):
        pass


class PanelGraficos(PanelBase):
    """Panel para mostrar el gráfico de la trayectoria y controles asociados."""
    def crear_widgets(self):
        # Usar PanedWindow para permitir redimensionar entre gráfico y controles
        paned_window = tk.PanedWindow(self, orient=tk.VERTICAL, sashrelief=tk.RAISED)
        paned_window.pack(fill="both", expand=True)

        # Frame superior para el gráfico
        frame_grafico = ttk.Frame(paned_window, relief=tk.SUNKEN)
        paned_window.add(frame_grafico) # Sin weight para compatibilidad

        # Frame inferior para controles e información
        frame_controles_info = ttk.Frame(paned_window)
        paned_window.add(frame_controles_info)

        # --- Gráfico --- 
        self.figura = plt.Figure(figsize=(6, 4), dpi=100) # Tamaño inicial más pequeño
        self.ax = self.figura.add_subplot(111)
        self.ax.set_title("Trayectoria del Proyectil")
        self.ax.set_xlabel("Distancia (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.grid(True)
        self.ax.axis('equal') # Intentar mantener aspecto igual
        self.figura.tight_layout() # Ajustar layout

        self.canvas = FigureCanvasTkAgg(self.figura, master=frame_grafico)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill="both", expand=True)
        self.canvas.draw()

        # --- Controles e Información --- 
        scrollable_controls = ScrollableFrame(frame_controles_info)
        scrollable_controls.pack(fill="both", expand=True)
        inner_controls_frame = scrollable_controls.scrollable_frame

        # Frame para sliders
        frame_sliders = ttk.Frame(inner_controls_frame)
        frame_sliders.pack(fill="x", padx=10, pady=5)
        frame_sliders.columnconfigure(1, weight=1)

        # Slider trayectoria 1
        self.slider_trayectoria_frame = ttk.Frame(frame_sliders)
        self.slider_trayectoria_frame.grid(row=0, column=0, columnspan=2, sticky='ew', pady=2)
        ttk.Label(self.slider_trayectoria_frame, text="Punto Tray. 1:").pack(side="left")
        self.slider_trayectoria = ttk.Scale(
            self.slider_trayectoria_frame,
            from_=0, to=100, orient="horizontal",
            command=self.controller.actualizar_punto_trayectoria, state=tk.DISABLED
        )
        self.slider_trayectoria.pack(side="left", expand=True, fill="x", padx=5)

        # Slider trayectoria 2 (comparación)
        self.slider_comparacion_frame = ttk.Frame(frame_sliders)
        # No se empaqueta inicialmente, se hace en actualizar_grafico
        ttk.Label(self.slider_comparacion_frame, text="Punto Tray. 2:").pack(side="left")
        self.slider_comparacion = ttk.Scale(
            self.slider_comparacion_frame,
            from_=0, to=100, orient="horizontal",
            command=self.controller.actualizar_punto_comparacion, state=tk.DISABLED
        )
        self.slider_comparacion.pack(side="left", expand=True, fill="x", padx=5)

        # Frame para información del punto
        self.frame_info = ttk.LabelFrame(inner_controls_frame, text="Información del Punto", padding=10)
        self.frame_info.pack(fill="x", padx=10, pady=5)

        # Info Trayectoria 1
        ttk.Label(self.frame_info, text="Trayectoria 1 (Azul):", font=("Helvetica", 9, "bold")).pack(anchor="w")
        self.lbl_posicion1 = ttk.Label(self.frame_info, text="Pos: (N/A, N/A) m | t: N/A s")
        self.lbl_posicion1.pack(anchor="w")
        self.lbl_velocidad1 = ttk.Label(self.frame_info, text="Vel: (N/A, N/A) m/s | |v|: N/A m/s")
        self.lbl_velocidad1.pack(anchor="w")
        self.lbl_aceleracion1 = ttk.Label(self.frame_info, text="Acel: (N/A, N/A) m/s² | |a|: N/A m/s²")
        self.lbl_aceleracion1.pack(anchor="w")

        # Separador e Info Trayectoria 2
        self.sep_info = ttk.Separator(self.frame_info, orient='horizontal')
        self.lbl_trayectoria2 = ttk.Label(self.frame_info, text="Trayectoria 2 (Roja):", font=("Helvetica", 9, "bold"))
        self.lbl_posicion2 = ttk.Label(self.frame_info, text="Pos: (N/A, N/A) m | t: N/A s")
        self.lbl_velocidad2 = ttk.Label(self.frame_info, text="Vel: (N/A, N/A) m/s | |v|: N/A m/s")
        self.lbl_aceleracion2 = ttk.Label(self.frame_info, text="Acel: (N/A, N/A) m/s² | |a|: N/A m/s²")
        # Se empaquetan/desempaquetan en actualizar_grafico

        # Frame para botones de acción del gráfico
        frame_botones_graf = ttk.Frame(inner_controls_frame)
        frame_botones_graf.pack(fill="x", padx=10, pady=5)
        
        ttk.Button(frame_botones_graf, text="Exportar Gráfico", command=self.controller.exportar_grafico).pack(side="left", padx=5)
        self.btn_animacion = ttk.Button(frame_botones_graf, text="Iniciar Animación", command=self.controller.toggle_animacion)
        self.btn_animacion.pack(side="left", padx=5)

        # Controles de animación (inicialmente ocultos)
        self.frame_anim_controls = ttk.Frame(inner_controls_frame)
        # No se empaqueta inicialmente
        ttk.Label(self.frame_anim_controls, text="Velocidad Animación:").pack(side="left")
        self.velocidad_animacion = tk.IntVar(value=50) # Representa el intervalo en ms
        scale_anim = ttk.Scale(
            self.frame_anim_controls,
            from_=200, to=10, # Invertido: más a la derecha = más rápido (menor intervalo)
            variable=self.velocidad_animacion,
            orient="horizontal", length=100,
            command=lambda v: self.controller.set_velocidad_animacion(int(float(v)))
        )
        scale_anim.pack(side="left", padx=5)
        
        # Elementos dinámicos del gráfico (se crean/actualizan después)
        self.linea_trayectoria1 = None
        self.punto_lanzamiento1 = None
        self.punto_impacto1 = None
        self.punto_interactivo1 = None
        self.vector_posicion1 = None
        self.vector_velocidad1 = None
        self.vector_aceleracion1 = None
        
        self.linea_trayectoria2 = None
        self.punto_lanzamiento2 = None
        self.punto_impacto2 = None
        self.punto_interactivo2 = None
        self.vector_posicion2 = None
        self.vector_velocidad2 = None
        
        self.anim = None # Para guardar la referencia a la animación

    def obtener_datos(self) -> Dict:
        return {}

    def establecer_datos(self, datos: Dict):
        pass

    def _limpiar_elementos_dinamicos(self, ax):
        """Elimina elementos gráficos que se regeneran."""
        elementos = [
            self.linea_trayectoria1, self.punto_lanzamiento1, self.punto_impacto1,
            self.punto_interactivo1, self.vector_posicion1, self.vector_velocidad1, self.vector_aceleracion1,
            self.linea_trayectoria2, self.punto_lanzamiento2, self.punto_impacto2,
            self.punto_interactivo2, self.vector_posicion2, self.vector_velocidad2
        ]
        for elem in elementos:
            if elem:
                try:
                    # Algunos son patches (flechas), otros son Line2D
                    if isinstance(elem, FancyArrowPatch):
                        elem.remove()
                    elif hasattr(elem, 'remove'): # Line2D
                         elem.remove()
                except ValueError: # Puede que ya haya sido eliminado
                    pass
        # Resetear referencias
        self.linea_trayectoria1 = self.punto_lanzamiento1 = self.punto_impacto1 = None
        self.punto_interactivo1 = self.vector_posicion1 = self.vector_velocidad1 = self.vector_aceleracion1 = None
        self.linea_trayectoria2 = self.punto_lanzamiento2 = self.punto_impacto2 = None
        self.punto_interactivo2 = self.vector_posicion2 = self.vector_velocidad2 = None
        
        # Limpiar leyenda también si existe
        if ax.get_legend():
            ax.get_legend().remove()

    def actualizar_grafico(self, res1: Optional[ResultadosTiro], res2: Optional[ResultadosTiro]):
        self.ax.clear() # Limpia completamente el eje
        self._limpiar_elementos_dinamicos(self.ax) # Asegura limpiar referencias
        self.ax.grid(True)
        self.ax.set_xlabel("Distancia (m)")
        self.ax.set_ylabel("Altura (m)")
        self.ax.set_title("Trayectoria del Proyectil")
        
        hay_trayectoria1 = res1 and res1.x and res1.y
        hay_trayectoria2 = res2 and res2.x and res2.y
        
        # Plot trayectoria 1 (Azul)
        if hay_trayectoria1:
            self.linea_trayectoria1, = self.ax.plot(res1.x, res1.y, 'b-', label="Trayectoria 1")
            self.punto_lanzamiento1, = self.ax.plot(res1.x[0], res1.y[0], 'bo', markersize=6, label="Inicio 1")
            self.punto_impacto1, = self.ax.plot(res1.alcance, 0, 'bx', markersize=8, markeredgewidth=2, label="Impacto 1")
            self.slider_trayectoria.config(state=tk.NORMAL, to=len(res1.x)-1 if len(res1.x)>1 else 0)
            self.slider_trayectoria.set(0)
            self.actualizar_punto_info(1, 0, res1) # Mostrar info inicial
        else:
            self.slider_trayectoria.config(state=tk.DISABLED)
            self.actualizar_punto_info(1, None, None) # Limpiar info

        # Plot trayectoria 2 (Roja) - Comparación
        if hay_trayectoria2:
            self.linea_trayectoria2, = self.ax.plot(res2.x, res2.y, 'r-', label="Trayectoria 2")
            self.punto_lanzamiento2, = self.ax.plot(res2.x[0], res2.y[0], 'ro', markersize=6, label="Inicio 2")
            self.punto_impacto2, = self.ax.plot(res2.alcance, 0, 'rx', markersize=8, markeredgewidth=2, label="Impacto 2")
            # Mostrar slider y labels para trayectoria 2
            self.slider_comparacion_frame.grid(row=1, column=0, columnspan=2, sticky='ew', pady=2)
            self.sep_info.pack(fill='x', pady=3, before=self.lbl_trayectoria2)
            self.lbl_trayectoria2.pack(anchor="w")
            self.lbl_posicion2.pack(anchor="w")
            self.lbl_velocidad2.pack(anchor="w")
            self.lbl_aceleracion2.pack(anchor="w")
            self.slider_comparacion.config(state=tk.NORMAL, to=len(res2.x)-1 if len(res2.x)>1 else 0)
            self.slider_comparacion.set(0)
            self.actualizar_punto_info(2, 0, res2) # Mostrar info inicial
        else:
            # Ocultar comparación
            self.slider_comparacion_frame.grid_forget()
            self.sep_info.pack_forget()
            self.lbl_trayectoria2.pack_forget()
            self.lbl_posicion2.pack_forget()
            self.lbl_velocidad2.pack_forget()
            self.lbl_aceleracion2.pack_forget()
            self.slider_comparacion.config(state=tk.DISABLED)
            self.actualizar_punto_info(2, None, None) # Limpiar info

        # Ajustar límites y leyenda
        if hay_trayectoria1 or hay_trayectoria2:
            self.ax.legend(fontsize='small')
            # Ajustar límites para incluir ambas trayectorias
            all_x = (res1.x if hay_trayectoria1 else []) + (res2.x if hay_trayectoria2 else [])
            all_y = (res1.y if hay_trayectoria1 else []) + (res2.y if hay_trayectoria2 else [])
            if all_x and all_y:
                min_x, max_x = min(all_x), max(all_x)
                min_y, max_y = 0, max(all_y) # Y siempre desde 0
                range_x = max_x - min_x
                range_y = max_y - min_y
                self.ax.set_xlim(min_x - 0.05 * range_x, max_x + 0.05 * range_x)
                self.ax.set_ylim(min_y - 0.05 * range_y, max_y + 0.1 * range_y) # Más espacio arriba
            self.ax.axis('equal') # Reintentar aspecto igual después de ajustar límites
        else:
             self.ax.set_xlim(-1, 10)
             self.ax.set_ylim(0, 5)
             
        self.figura.tight_layout()
        self.canvas.draw()

    def _crear_o_actualizar_punto(self, punto_attr, x, y, color, marker='o', size=8):
        punto_actual = getattr(self, punto_attr)
        if punto_actual is None:
            punto_actual, = self.ax.plot([x], [y], color + marker, markersize=size, label=f"Actual {punto_attr[-1]}")
            setattr(self, punto_attr, punto_actual)
        else:
            punto_actual.set_data([x], [y])
        return punto_actual

    def _crear_o_actualizar_vector(self, vector_attr, x_base, y_base, vx, vy, color, escala=0.5, estilo='->', mutacion=15):
        vector_actual = getattr(self, vector_attr)
        if vector_actual:
            vector_actual.remove()
            
        # Evitar dibujar vector si la magnitud es cero
        if vx == 0 and vy == 0:
             setattr(self, vector_attr, None)
             return None
             
        vector_actual = FancyArrowPatch(
            (x_base, y_base), (x_base + vx * escala, y_base + vy * escala),
            arrowstyle=estilo, mutation_scale=mutacion,
            color=color, linewidth=1.5, alpha=0.8
        )
        self.ax.add_patch(vector_actual)
        setattr(self, vector_attr, vector_actual)
        return vector_actual

    def actualizar_punto_info(self, num_trayectoria: int, idx: Optional[int], resultados: Optional[ResultadosTiro]):
        if idx is None or resultados is None or idx >= len(resultados.x):
            # Limpiar información
            pos_text = "Pos: (N/A, N/A) m | t: N/A s"
            vel_text = "Vel: (N/A, N/A) m/s | |v|: N/A m/s"
            acel_text = "Acel: (N/A, N/A) m/s² | |a|: N/A m/s²"
            x, y = None, None # Para limpiar puntos/vectores
        else:
            x = resultados.x[idx]
            y = resultados.y[idx]
            t = resultados.tiempos[idx]
            
            # Calcular velocidad y aceleración en el punto idx (aproximado)
            params = self.controller.parametros_actuales if num_trayectoria == 1 else self.controller.parametros_comparacion
            g = params.g
            masa = params.masa
            Cd = params.coef_arrastre
            A = params.area
            rho = params.densidad_aire
            resistencia = params.resistencia
            
            if idx == 0:
                vx = resultados.v0x
                vy = resultados.v0y
            elif len(resultados.tiempos) > 1:
                dt = resultados.tiempos[idx] - resultados.tiempos[idx-1]
                vx = (resultados.x[idx] - resultados.x[idx-1]) / dt if dt > 0 else 0
                vy = (resultados.y[idx] - resultados.y[idx-1]) / dt if dt > 0 else 0
            else: # Solo un punto
                vx = 0
                vy = 0
                
            velocidad_mag = np.sqrt(vx**2 + vy**2)
            
            if resistencia:
                if velocidad_mag == 0:
                    ax = 0
                    ay = -g
                else:
                    f_arrastre = 0.5 * rho * Cd * A * velocidad_mag**2
                    ax = - (f_arrastre / masa) * (vx / velocidad_mag)
                    ay = -g - (f_arrastre / masa) * (vy / velocidad_mag)
            else:
                ax = 0
                ay = -g
                
            aceleracion_mag = np.sqrt(ax**2 + ay**2)
            
            pos_text = f"Pos: ({x:.2f}, {y:.2f}) m | t: {t:.2f} s"
            vel_text = f"Vel: ({vx:.2f}, {vy:.2f}) m/s | |v|: {velocidad_mag:.2f} m/s"
            acel_text = f"Acel: ({ax:.2f}, {ay:.2f}) m/s² | |a|: {aceleracion_mag:.2f} m/s²"

        # Actualizar etiquetas de texto
        if num_trayectoria == 1:
            self.lbl_posicion1.config(text=pos_text)
            self.lbl_velocidad1.config(text=vel_text)
            self.lbl_aceleracion1.config(text=acel_text)
            color = 'b'
            punto_attr = 'punto_interactivo1'
            vec_pos_attr = 'vector_posicion1'
            vec_vel_attr = 'vector_velocidad1'
            vec_acel_attr = 'vector_aceleracion1'
            x0, y0 = (resultados.x[0], resultados.y[0]) if resultados else (0,0)
        else:
            self.lbl_posicion2.config(text=pos_text)
            self.lbl_velocidad2.config(text=vel_text)
            self.lbl_aceleracion2.config(text=acel_text)
            color = 'r'
            punto_attr = 'punto_interactivo2'
            vec_pos_attr = 'vector_posicion2'
            vec_vel_attr = 'vector_velocidad2'
            vec_acel_attr = None # No mostrar aceleración para el 2do
            x0, y0 = (resultados.x[0], resultados.y[0]) if resultados else (0,0)

        # Actualizar punto interactivo y vectores en el gráfico
        if x is not None and y is not None:
            self._crear_o_actualizar_punto(punto_attr, x, y, color, marker='o', size=7)
            # Vector Posición (desde origen)
            self._crear_o_actualizar_vector(vec_pos_attr, x0, y0, x-x0, y-y0, 'g', escala=1, estilo='->', mutacion=10)
            # Vector Velocidad (desde el punto)
            self._crear_o_actualizar_vector(vec_vel_attr, x, y, vx, vy, 'm', escala=0.3, estilo='->', mutacion=12)
            # Vector Aceleración (solo para trayectoria 1)
            if vec_acel_attr:
                 self._crear_o_actualizar_vector(vec_acel_attr, x, y, ax, ay, 'c', escala=0.15, estilo='->', mutacion=12)
        else:
            # Limpiar punto y vectores si no hay datos
            punto = getattr(self, punto_attr)
            if punto: punto.remove(); setattr(self, punto_attr, None)
            vec_p = getattr(self, vec_pos_attr)
            if vec_p: vec_p.remove(); setattr(self, vec_pos_attr, None)
            vec_v = getattr(self, vec_vel_attr)
            if vec_v: vec_v.remove(); setattr(self, vec_vel_attr, None)
            if vec_acel_attr:
                vec_a = getattr(self, vec_acel_attr)
                if vec_a: vec_a.remove(); setattr(self, vec_acel_attr, None)

        self.canvas.draw_idle() # Usar draw_idle para eficiencia

    def iniciar_animacion(self, res1: ResultadosTiro, res2: Optional[ResultadosTiro]):
        if self.anim is not None:
            self.parar_animacion() # Detener animación anterior si existe
            
        self.btn_animacion.config(text="Parar Animación")
        self.frame_anim_controls.pack(side="left", padx=10, after=self.btn_animacion) # Mostrar controles
        
        frames_totales = len(res1.tiempos)
        intervalo_ms = self.velocidad_animacion.get()
        
        # Crear puntos y vectores iniciales para la animación
        p1 = self._crear_o_actualizar_punto('punto_interactivo1', res1.x[0], res1.y[0], 'b', marker='o', size=7)
        vp1 = self._crear_o_actualizar_vector('vector_posicion1', res1.x[0], res1.y[0], 0, 0, 'g', escala=1)
        vv1 = self._crear_o_actualizar_vector('vector_velocidad1', res1.x[0], res1.y[0], 0, 0, 'm', escala=0.3)
        va1 = self._crear_o_actualizar_vector('vector_aceleracion1', res1.x[0], res1.y[0], 0, 0, 'c', escala=0.15)
        elementos_anim1 = [p1, vp1, vv1, va1]
        
        elementos_anim2 = []
        if res2:
            p2 = self._crear_o_actualizar_punto('punto_interactivo2', res2.x[0], res2.y[0], 'r', marker='o', size=7)
            vp2 = self._crear_o_actualizar_vector('vector_posicion2', res2.x[0], res2.y[0], 0, 0, 'g', escala=1)
            vv2 = self._crear_o_actualizar_vector('vector_velocidad2', res2.x[0], res2.y[0], 0, 0, 'm', escala=0.3)
            elementos_anim2 = [p2, vp2, vv2]
            frames_totales = max(frames_totales, len(res2.tiempos))
            
        def update(frame_idx):
            # Actualizar trayectoria 1
            idx1 = min(frame_idx, len(res1.x) - 1)
            self.actualizar_punto_info(1, idx1, res1)
            self.slider_trayectoria.set(idx1)
            
            # Actualizar trayectoria 2 si existe
            if res2:
                idx2 = min(frame_idx, len(res2.x) - 1)
                self.actualizar_punto_info(2, idx2, res2)
                self.slider_comparacion.set(idx2)
                
            # Devolver los elementos gráficos que se actualizaron
            return elementos_anim1 + elementos_anim2

        self.anim = animation.FuncAnimation(
            self.figura, update, frames=frames_totales,
            interval=intervalo_ms, blit=True, repeat=False
        )
        self.canvas.draw()

    def parar_animacion(self):
        if self.anim is not None:
            # Detener la animación existente
            self.anim.event_source.stop()
            self.anim = None
            self.btn_animacion.config(text="Iniciar Animación")
            self.frame_anim_controls.pack_forget() # Ocultar controles
            # Restaurar puntos/vectores a la posición del slider
            self.controller.actualizar_punto_trayectoria(self.slider_trayectoria.get())
            if self.controller.parametros_comparacion.comparar:
                 self.controller.actualizar_punto_comparacion(self.slider_comparacion.get())
            print("Animación detenida.")
        else:
             print("No hay animación activa para detener.")
             self.btn_animacion.config(text="Iniciar Animación") # Asegurar texto correcto
             self.frame_anim_controls.pack_forget()

    def set_velocidad_animacion(self, intervalo_ms: int):
        if self.anim:
            self.anim.event_source.interval = intervalo_ms
            print(f"Intervalo de animación ajustado a {intervalo_ms} ms")


class PanelHistorial(PanelBase):
    """Panel para mostrar y gestionar el historial de simulaciones."""
    def crear_widgets(self):
        scrollable = ScrollableFrame(self)
        scrollable.pack(fill="both", expand=True)
        frame = ttk.LabelFrame(scrollable.scrollable_frame, text="Historial de Simulaciones", padding=10)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Listbox para mostrar el historial
        self.historial_listbox = tk.Listbox(frame, width=80, height=15)
        self.historial_listbox.pack(padx=5, pady=5, fill="both", expand=True)
        # Añadir scrollbar horizontal también si el texto es largo
        scrollbar_x = ttk.Scrollbar(frame, orient="horizontal", command=self.historial_listbox.xview)
        self.historial_listbox.configure(xscrollcommand=scrollbar_x.set)
        scrollbar_x.pack(fill="x", padx=5, before=self.historial_listbox)

        # Frame para botones de acción del historial
        frame_botones = ttk.Frame(frame)
        frame_botones.pack(fill="x", padx=5, pady=5)

        ttk.Button(frame_botones, text="Cargar Selección", command=self.controller.cargar_desde_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Eliminar Selección", command=self.controller.eliminar_del_historial).pack(side="left", padx=5)
        ttk.Button(frame_botones, text="Limpiar Historial", command=self.controller.limpiar_historial).pack(side="right", padx=5)

    def obtener_datos(self) -> Dict:
        return {}

    def establecer_datos(self, datos: Dict):
        pass

    def actualizar_lista_historial(self, historial: List[Dict]):
        self.historial_listbox.delete(0, tk.END)
        for i, entrada in enumerate(historial):
            try:
                timestamp = entrada.get('timestamp', time.time())
                fecha_hora = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timestamp))
                datos = entrada.get('datos', {})
                res = entrada.get('resultados', {})
                
                # Crear texto resumen
                resistencia_str = "R" if datos.get('resistencia', False) else "-"
                texto = (
                    f"{i+1:02d}. [{fecha_hora}] v₀={datos.get('v0', '?'):.1f}, θ={datos.get('angulo', '?'):.1f}°, {resistencia_str} | "
                    f"Alc={res.get('alcance', '?'):.1f}m, Hmax={res.get('altura_max', '?'):.1f}m, T={res.get('t_vuelo', '?'):.1f}s"
                )
                self.historial_listbox.insert(tk.END, texto)
            except Exception as e:
                # Si hay un error en una entrada, mostrarlo pero continuar
                self.historial_listbox.insert(tk.END, f"{i+1:02d}. [Error al cargar entrada: {e}]")
                print(f"Error procesando entrada de historial {i}: {e}")

    def obtener_seleccion_indice(self) -> Optional[int]:
        seleccion = self.historial_listbox.curselection()
        return seleccion[0] if seleccion else None


class BarraHerramientas(ttk.Frame):
    """Barra de herramientas principal de la aplicación."""
    def __init__(self, parent, controller):
        super().__init__(parent, height=40, relief=tk.RAISED, borderwidth=1)
        self.controller = controller
        self.pack(side=tk.TOP, fill=tk.X)
        # self.grid_propagate(False) # Evitar que los widgets internos cambien el tamaño del frame
        self.crear_widgets()

    def crear_widgets(self):
        self.style = ttk.Style()
        # Estilo para el botón principal
        self.style.configure('Accent.TButton', font=('Helvetica', 10, 'bold'), foreground='blue')
        # Estilo para botones secundarios
        self.style.configure('Tool.TButton', font=('Helvetica', 9))

        # Frame para centrar los botones principales
        btn_frame_main = ttk.Frame(self)
        btn_frame_main.pack(side=tk.LEFT, padx=10)

        ttk.Button(btn_frame_main, text="Calcular Trayectoria", command=self.controller.calcular_trayectoria,
                  style='Accent.TButton').pack(side='left', padx=5, pady=5)
        ttk.Button(btn_frame_main, text="Limpiar Todo", command=self.controller.limpiar,
                   style='Tool.TButton').pack(side='left', padx=5, pady=5)

        # Frame para botones de exportación/otros
        btn_frame_extra = ttk.Frame(self)
        btn_frame_extra.pack(side=tk.LEFT, padx=10)
        ttk.Button(btn_frame_extra, text="Exportar Informe", command=self.controller.exportar_informe,
                   style='Tool.TButton').pack(side='left', padx=5, pady=5)
        # ttk.Button(btn_frame_extra, text="Modo Avanzado", command=self.controller.modo_avanzado,
        #            style='Tool.TButton').pack(side='left', padx=5, pady=5) # Funcionalidad futura

        # Etiqueta de estado o fecha/hora a la derecha
        self.status_label = ttk.Label(self, text="Listo")
        self.status_label.pack(side='right', padx=10, pady=5)
        # Actualizar hora periódicamente
        self.actualizar_hora()

    def actualizar_hora(self):
        fecha_hora = time.strftime('%d/%m/%Y %H:%M:%S')
        self.status_label.config(text=fecha_hora)
        self.after(1000, self.actualizar_hora) # Actualizar cada segundo

    def set_status(self, message: str):
        # Podría usarse para mostrar mensajes como "Calculando..."
        # self.status_label.config(text=message)
        pass # Por ahora solo muestra la hora

# ==============================================
# Clase principal de la aplicación (Controlador)
# ==============================================

class TiroParabolicoApp:
    """Clase principal que coordina el modelo y la vista."""
    def __init__(self, root):
        self.root = root
        self.root.title("Simulador de Tiro Parabólico Interactivo")
        # Geometría inicial, pero permitir redimensionar
        self.root.geometry("1200x800")
        self.root.minsize(800, 600) # Tamaño mínimo

        # Estilo ttk
        self.style = ttk.Style()
        # Intentar usar un tema más moderno si está disponible
        available_themes = self.style.theme_names()
        print(f"Temas disponibles: {available_themes}")
        if 'clam' in available_themes:
            self.style.theme_use('clam')
        elif 'vista' in available_themes:
             self.style.theme_use('vista') # Para Windows
        elif 'aqua' in available_themes:
             self.style.theme_use('aqua') # Para macOS

        # Modelo
        self.modelo = ModeloTiroParabolico()
        self.parametros_actuales = ParametrosTiro() # Valores por defecto
        self.parametros_comparacion = ParametrosTiro(comparar=False, v0=15.0, angulo=30.0) # Valores por defecto comparación
        self.resultados_actuales: Optional[ResultadosTiro] = None
        self.resultados_comparacion: Optional[ResultadosTiro] = None
        self.historial: List[Dict] = []
        self.historial_file = "historial_simulaciones.json"
        self.animacion_activa = False
        self.velocidad_anim_ms = 50 # Intervalo en ms

        # Cargar historial
        self.cargar_historial()

        # Crear interfaz
        self.crear_interfaz()
        
        # Calcular trayectoria inicial con valores por defecto
        self.calcular_trayectoria()

    def crear_interfaz(self):
        # Barra de herramientas superior
        self.barra_herramientas = BarraHerramientas(self.root, self)

        # PanedWindow principal para dividir parámetros y resultados/gráficos
        main_paned_window = tk.PanedWindow(self.root, orient=tk.HORIZONTAL, sashrelief=tk.RAISED)
        main_paned_window.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # --- Panel Izquierdo (Parámetros e Historial) --- 
        left_pane = tk.PanedWindow(main_paned_window, orient=tk.VERTICAL, sashrelief=tk.RAISED)
        main_paned_window.add(left_pane, width=350, stretch="never") # Ancho fijo inicial

        # Notebook para organizar paneles de parámetros
        param_notebook = ttk.Notebook(left_pane)
        left_pane.add(param_notebook) # Añadir sin weight para compatibilidad

        self.panel_basico = PanelParametrosBasicos(param_notebook, self)
        param_notebook.add(self.panel_basico, text='Básicos')

        self.panel_resistencia = PanelResistenciaAire(param_notebook, self)
        param_notebook.add(self.panel_resistencia, text='Resistencia Aire')

        self.panel_comparacion = PanelComparacion(param_notebook, self)
        param_notebook.add(self.panel_comparacion, text='Comparación')

        # self.panel_unidades = PanelUnidades(param_notebook, self)
        # param_notebook.add(self.panel_unidades, text='Unidades')

        # Panel de Historial debajo de los parámetros
        self.panel_historial = PanelHistorial(left_pane, self)
        left_pane.add(self.panel_historial)

        # --- Panel Derecho (Gráficos y Resultados) --- 
        right_pane = tk.PanedWindow(main_paned_window, orient=tk.VERTICAL, sashrelief=tk.RAISED)
        main_paned_window.add(right_pane, stretch="always") # Este panel se estira

        # Panel de Gráficos en la parte superior derecha
        self.panel_graficos = PanelGraficos(right_pane, self)
        right_pane.add(self.panel_graficos) # Sin weight para compatibilidad

        # Panel de Resultados en la parte inferior derecha
        self.panel_resultados = PanelResultados(right_pane, self)
        right_pane.add(self.panel_resultados)
        
        # Configurar estado inicial de paneles dependientes
        self.toggle_resistencia() # Actualiza estado panel resistencia
        self.toggle_comparacion() # Actualiza estado panel comparación

    # --- Métodos del Controlador (Callbacks y Lógica) --- 

    def toggle_resistencia(self):
        """Actualiza estado del panel de resistencia basado en el checkbox."""
        activar_resistencia = self.panel_basico.vars['resistencia'].get()
        self.panel_resistencia.set_enabled(activar_resistencia)
        # Habilitar/deshabilitar consulta por tiempo
        self.panel_resultados.set_consulta_enabled(not activar_resistencia)
        # Actualizar parámetros actuales
        self.parametros_actuales.resistencia = activar_resistencia

    def toggle_comparacion(self):
        """Actualiza estado del panel de comparación y recalcula si es necesario."""
        activar_comparacion = self.panel_comparacion.vars['comparar'].get()
        self.panel_comparacion.set_enabled(activar_comparacion)
        # Actualizar parámetro
        self.parametros_comparacion.comparar = activar_comparacion
        # Recalcular y actualizar gráfico para mostrar/ocultar la 2da trayectoria
        if activar_comparacion:
            self.calcular_trayectoria_comparacion()
        else:
            self.resultados_comparacion = None
            self.actualizar_interfaz()

    def actualizar_unidades(self):
        """(Funcionalidad Futura) Maneja el cambio de sistema de unidades."""
        sistema = self.panel_unidades.var_unidades.get()
        messagebox.showinfo("Información", f"Cambio a sistema {sistema} no implementado aún.")
        # Aquí iría la lógica de conversión de valores en los paneles

    def validar_parametros(self, parametros: ParametrosTiro, nombre: str = "principal") -> bool:
        """Valida los parámetros de entrada antes del cálculo."""
        try:
            if parametros.v0 <= 0:
                raise ValueError("La velocidad inicial (|v₀|) debe ser positiva.")
            if not (0 <= parametros.angulo <= 90):
                raise ValueError("El ángulo de lanzamiento (θ) debe estar entre 0° y 90°.")
            if parametros.g <= 0:
                 raise ValueError("La aceleración de la gravedad (g) debe ser positiva.")
            if parametros.resistencia:
                if parametros.masa <= 0:
                    raise ValueError("La masa debe ser positiva.")
                if parametros.coef_arrastre < 0: # Puede ser 0
                    raise ValueError("El coeficiente de arrastre (C<sub>d</sub>) no puede ser negativo.")
                if parametros.area <= 0:
                    raise ValueError("El área transversal (A) debe ser positiva.")
                if parametros.densidad_aire <= 0:
                    raise ValueError("La densidad del aire (ρ) debe ser positiva.")
            return True
        except ValueError as e:
            messagebox.showerror("Error de Validación", f"Parámetros ({nombre}): {str(e)}")
            return False
        except tk.TclError as e:
             messagebox.showerror("Error de Entrada", f"Parámetros ({nombre}): Valor inválido detectado - {str(e)}")
             return False # Error al obtener valor de tk.Var

    def calcular_trayectoria_principal(self):
        """Obtiene datos, valida y calcula la trayectoria principal."""
        try:
            datos_basicos = self.panel_basico.obtener_datos()
            datos_resistencia = self.panel_resistencia.obtener_datos()
            self.parametros_actuales = ParametrosTiro(**datos_basicos, **datos_resistencia)
        except tk.TclError:
             messagebox.showerror("Error de Entrada", "Error al leer los parámetros principales. Verifique los valores.")
             return False
             
        if not self.validar_parametros(self.parametros_actuales, "principal"):
            return False
            
        print("Calculando trayectoria principal...")
        self.resultados_actuales = self.modelo.calcular_trayectoria(self.parametros_actuales)
        print("Cálculo principal completado.")
        return True

    def calcular_trayectoria_comparacion(self):
        """Obtiene datos, valida y calcula la trayectoria de comparación."""
        if not self.panel_comparacion.vars['comparar'].get():
            self.resultados_comparacion = None
            return True # No hay nada que calcular
            
        try:
            datos_comp = self.panel_comparacion.obtener_datos()
            # Usar la misma 'g' que la principal por defecto
            g_actual = self.parametros_actuales.g
            # Crear ParametrosTiro para comparación
            self.parametros_comparacion = ParametrosTiro(
                x0=datos_comp['x0'], y0=datos_comp['y0'], v0=datos_comp['v0'],
                angulo=datos_comp['angulo'], g=g_actual, # Usar g principal
                resistencia=datos_comp['resistencia'],
                # Usar params de resistencia principal si la comparación usa resistencia
                masa=self.parametros_actuales.masa if datos_comp['resistencia'] else 1.0,
                coef_arrastre=self.parametros_actuales.coef_arrastre if datos_comp['resistencia'] else 0.47,
                area=self.parametros_actuales.area if datos_comp['resistencia'] else 0.1,
                densidad_aire=self.parametros_actuales.densidad_aire if datos_comp['resistencia'] else 1.225
            )
        except tk.TclError:
             messagebox.showerror("Error de Entrada", "Error al leer los parámetros de comparación. Verifique los valores.")
             return False

        if not self.validar_parametros(self.parametros_comparacion, "comparación"):
            return False
            
        print("Calculando trayectoria de comparación...")
        self.resultados_comparacion = self.modelo.calcular_trayectoria(self.parametros_comparacion)
        print("Cálculo de comparación completado.")
        return True

    def calcular_trayectoria(self):
        """Calcula la trayectoria principal y la de comparación si está activa."""
        # Detener animación si está activa antes de recalcular
        if self.animacion_activa:
            self.toggle_animacion()
            
        if not self.calcular_trayectoria_principal():
            return # Error en principal, no continuar
            
        if not self.calcular_trayectoria_comparacion():
            # Error en comparación, pero mostrar la principal
            self.resultados_comparacion = None 
            
        self.actualizar_interfaz()
        self.guardar_en_historial()

    def actualizar_interfaz(self):
        """Actualiza todos los paneles de la interfaz con los nuevos resultados."""
        print("Actualizando interfaz...")
        # Actualizar resultados textuales
        informe_res = self.generar_informe_resultados()
        self.panel_resultados.actualizar_resultados(informe_res)
        
        # Actualizar energías
        informe_energia = self.generar_informe_energias()
        self.panel_resultados.actualizar_energias(informe_energia)
        
        # Actualizar gráfico (pasa ambos resultados, pueden ser None)
        self.panel_graficos.actualizar_grafico(self.resultados_actuales, self.resultados_comparacion)
        print("Interfaz actualizada.")

    def generar_informe_resultados(self) -> str:
        """Genera el texto para el panel de resultados."""
        if not self.resultados_actuales:
            return "No hay resultados para mostrar. Presione 'Calcular Trayectoria'."
            
        res = self.resultados_actuales
        param = self.parametros_actuales
        
        informe = f"--- TRAYECTORIA PRINCIPAL ---\n"
        informe += f" Parámetros:\n"
        informe += f"  Pos. Inicial: ({param.x0:.2f}, {param.y0:.2f}) m\n"
        informe += f"  Vel. Inicial: {param.v0:.2f} m/s, Ángulo: {param.angulo:.2f}°\n"
        informe += f"  Gravedad: {param.g:.2f} m/s²\n"
        if param.resistencia:
            informe += f"  Resistencia Aire: SÍ (m={param.masa:.2f}kg, Cd={param.coef_arrastre:.3f}, A={param.area:.3f}m², ρ={param.densidad_aire:.3f}kg/m³)\n"
        else:
            informe += f"  Resistencia Aire: NO\n"
        informe += f" Resultados:\n"
        informe += f"  Tiempo de Vuelo: {res.t_vuelo:.3f} s\n"
        informe += f"  Alcance Máximo: {res.alcance:.3f} m\n"
        informe += f"  Altura Máxima: {res.altura_max:.3f} m\n"
        v_final_mag = np.sqrt(res.vx_final**2 + res.vy_final**2)
        informe += f"  Velocidad Final (Impacto): ({res.vx_final:.2f}, {res.vy_final:.2f}) m/s, |v|={v_final_mag:.2f} m/s\n"
        
        if self.resultados_comparacion:
            res_comp = self.resultados_comparacion
            param_comp = self.parametros_comparacion
            informe += f"\n--- TRAYECTORIA COMPARACIÓN ---\n"
            informe += f" Parámetros:\n"
            informe += f"  Pos. Inicial: ({param_comp.x0:.2f}, {param_comp.y0:.2f}) m\n"
            informe += f"  Vel. Inicial: {param_comp.v0:.2f} m/s, Ángulo: {param_comp.angulo:.2f}°\n"
            if param_comp.resistencia:
                 informe += f"  Resistencia Aire: SÍ\n"
            else:
                 informe += f"  Resistencia Aire: NO\n"
            informe += f" Resultados Comparación:\n"
            informe += f"  Tiempo de Vuelo: {res_comp.t_vuelo:.3f} s\n"
            informe += f"  Alcance Máximo: {res_comp.alcance:.3f} m\n"
            informe += f"  Altura Máxima: {res_comp.altura_max:.3f} m\n"
            v_final_mag_comp = np.sqrt(res_comp.vx_final**2 + res_comp.vy_final**2)
            informe += f"  Velocidad Final (Impacto): ({res_comp.vx_final:.2f}, {res_comp.vy_final:.2f}) m/s, |v|={v_final_mag_comp:.2f} m/s\n"
            
            # Diferencias
            informe += f" Diferencias (Comp - Princ):\n"
            informe += f"  Δ Tiempo Vuelo: {res_comp.t_vuelo - res.t_vuelo:+.3f} s\n"
            informe += f"  Δ Alcance: {res_comp.alcance - res.alcance:+.3f} m\n"
            informe += f"  Δ Altura Máx: {res_comp.altura_max - res.altura_max:+.3f} m\n"
            
        return informe

    def generar_informe_energias(self) -> str:
        """Genera el texto para el panel de energías."""
        if not self.resultados_actuales or not self.resultados_actuales.energias:
            return "N/A"
            
        ener = self.resultados_actuales.energias
        param = self.parametros_actuales
        
        if not ener['cinetica'] or not ener['potencial'] or not ener['total']:
             return "Datos de energía incompletos."
             
        texto = f"--- ENERGÍAS (Trayectoria Principal) ---\n"
        texto += f" Inicial:\n"
        texto += f"  E. Cinética (Ec): {ener['cinetica'][0]:.2f} J\n"
        texto += f"  E. Potencial (Ep): {ener['potencial'][0]:.2f} J\n"
        texto += f"  E. Total (Et): {ener['total'][0]:.2f} J\n"
        texto += f" Final (Impacto):\n"
        texto += f"  E. Cinética (Ec): {ener['cinetica'][-1]:.2f} J\n"
        texto += f"  E. Potencial (Ep): {ener['potencial'][-1]:.2f} J\n"
        texto += f"  E. Total (Et): {ener['total'][-1]:.2f} J\n"
        
        if param.resistencia:
            e_inicial = ener['total'][0]
            e_final = ener['total'][-1]
            if e_inicial > 1e-6: # Evitar división por cero
                perdida_abs = e_inicial - e_final
                perdida_rel = (perdida_abs / e_inicial) * 100
                texto += f"\n Pérdida de Energía (Resistencia):\n"
                texto += f"  ΔE Total: {perdida_abs:.2f} J ({perdida_rel:.2f}%)\n"
            else:
                 texto += f"\n Pérdida de Energía: N/A (E inicial ≈ 0)\n"
        else:
            # Verificar conservación (puede haber pequeños errores numéricos)
            delta_e = abs(ener['total'][-1] - ener['total'][0])
            texto += f"\n Conservación Energía: {'Sí' if delta_e < 1e-3 else f'No (ΔE ≈ {delta_e:.3f} J)'}\n"
            
        # Podría añadirse energía para la comparación si existe
        
        return texto

    def actualizar_punto_trayectoria(self, valor_slider):
        """Callback cuando se mueve el slider de la trayectoria principal."""
        if not self.resultados_actuales or self.animacion_activa:
            return
        try:
            idx = int(float(valor_slider))
            self.panel_graficos.actualizar_punto_info(1, idx, self.resultados_actuales)
        except ValueError:
            pass # Ignorar si el valor no es válido

    def actualizar_punto_comparacion(self, valor_slider):
        """Callback cuando se mueve el slider de la trayectoria de comparación."""
        if not self.resultados_comparacion or self.animacion_activa:
            return
        try:
            idx = int(float(valor_slider))
            self.panel_graficos.actualizar_punto_info(2, idx, self.resultados_comparacion)
        except ValueError:
            pass

    def consultar_tiempo(self):
        """Calcula y muestra el estado del proyectil en un tiempo t (solo sin resistencia)."""
        if self.parametros_actuales.resistencia:
            messagebox.showwarning("Consulta no disponible", "La consulta por tiempo solo está disponible para simulaciones sin resistencia del aire.")
            return
            
        try:
            t_consulta_str = self.panel_resultados.tiempo_consulta.get()
            if not t_consulta_str:
                 raise ValueError("Ingrese un valor de tiempo.")
            t = float(t_consulta_str)
            if t < 0:
                raise ValueError("El tiempo no puede ser negativo.")
                
            # Usar parámetros actuales (sin resistencia)
            param = self.parametros_actuales
            estado = self.modelo.parametrizacion_tiro_sin_resistencia(
                param.x0, param.y0, param.v0, param.angulo, param.g, t
            )
            
            # Verificar si el tiempo es mayor al tiempo de vuelo
            t_vuelo_actual = self.resultados_actuales.t_vuelo if self.resultados_actuales else 0
            
            resultado_str = f"--- ESTADO EN t = {t:.3f} s ---\n"
            if t > t_vuelo_actual:
                 resultado_str += f"(Tiempo mayor al tiempo de vuelo total: {t_vuelo_actual:.3f} s)\n"
                 resultado_str += f"El proyectil ya impactó en el suelo.\n"
            elif estado['posicion'][1] < -1e-6: # Pequeña tolerancia por si acaso
                 resultado_str += f"El proyectil ya impactó en el suelo (y(t) < 0).\n"
            else:
                pos = estado['posicion']
                vel = estado['velocidad']
                acel = estado['aceleracion']
                vel_mag = np.sqrt(vel[0]**2 + vel[1]**2)
                acel_mag = np.sqrt(acel[0]**2 + acel[1]**2)
                
                resultado_str += f" Posición: ({pos[0]:.3f}, {pos[1]:.3f}) m\n"
                resultado_str += f" Velocidad: ({vel[0]:.3f}, {vel[1]:.3f}) m/s, |v|={vel_mag:.3f} m/s\n"
                resultado_str += f" Aceleración: ({acel[0]:.3f}, {acel[1]:.3f}) m/s², |a|={acel_mag:.3f} m/s²\n"
                
            resultado_str += f"\n--- Ecuaciones Paramétricas (sin resistencia) ---\n"
            resultado_str += f" x(t) = {estado['ecuaciones']['x(t)']}\n"
            resultado_str += f" y(t) = {estado['ecuaciones']['y(t)']}\n"
            resultado_str += f" vx(t) = {estado['ecuaciones']['vx(t)']}\n"
            resultado_str += f" vy(t) = {estado['ecuaciones']['vy(t)']}\n"
            
            self.panel_resultados.mostrar_consulta(resultado_str)
            
        except ValueError as e:
            messagebox.showerror("Error en Consulta", str(e))
        except tk.TclError:
             messagebox.showerror("Error de Entrada", "Valor de tiempo inválido.")
        except Exception as e:
             messagebox.showerror("Error Inesperado", f"Ocurrió un error al consultar: {e}")
             self.panel_resultados.mostrar_consulta("Error al realizar la consulta.")

    def exportar_informe(self):
        """Exporta los informes de resultados y energías a un archivo."""
        if not self.resultados_actuales:
            messagebox.showwarning("Sin Datos", "No hay resultados para exportar. Calcule una trayectoria primero.")
            return
            
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Archivos de texto", "*.txt"), ("Todos los archivos", "*.*")],
            title="Guardar Informe Como",
            initialfile=f"informe_tiro_{time.strftime('%Y%m%d_%H%M%S')}.txt"
        )
        
        if not filepath:
            return # Usuario canceló
            
        try:
            contenido = self.generar_informe_resultados()
            contenido += "\n" + "="*40 + "\n\n"
            contenido += self.generar_informe_energias()
            
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(contenido)
            
            messagebox.showinfo("Éxito", f"Informe guardado correctamente en:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Error al Guardar", f"No se pudo guardar el informe:\n{str(e)}")

    def exportar_grafico(self):
        """Exporta el gráfico actual a un archivo de imagen."""
        if not hasattr(self.panel_graficos, 'figura') or not self.resultados_actuales:
            messagebox.showwarning("Sin Gráfico", "No hay gráfico para exportar o no se ha calculado ninguna trayectoria.")
            return
            
        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("JPEG", "*.jpg"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("Todos los archivos", "*.*")],
            title="Guardar Gráfico Como",
            initialfile=f"grafico_tiro_{time.strftime('%Y%m%d_%H%M%S')}.png"
        )
        
        if not filepath:
            return # Usuario canceló
            
        try:
            # Asegurar que el layout esté ajustado antes de guardar
            self.panel_graficos.figura.tight_layout()
            self.panel_graficos.figura.savefig(filepath, dpi=300, bbox_inches='tight')
            messagebox.showinfo("Éxito", f"Gráfico guardado correctamente en:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Error al Guardar", f"No se pudo guardar el gráfico:\n{str(e)}")

    def limpiar(self):
        """Limpia todos los paneles, resultados y el gráfico."""
        print("Limpiando interfaz...")
        # Detener animación si está activa
        if self.animacion_activa:
            self.toggle_animacion()
            
        # Resetear parámetros a valores por defecto
        default_params = ParametrosTiro() # Crea instancia con defaults
        default_comp_params = ParametrosTiro(comparar=False, v0=15.0, angulo=30.0)
        
        self.panel_basico.establecer_datos(asdict(default_params))
        self.panel_resistencia.establecer_datos(asdict(default_params))
        self.panel_comparacion.establecer_datos(asdict(default_comp_params))
        # self.panel_unidades.establecer_datos({'unidades': "métrico"})
        
        # Actualizar estado de paneles dependientes
        self.toggle_resistencia()
        self.toggle_comparacion()
        
        # Limpiar resultados y modelo
        self.resultados_actuales = None
        self.resultados_comparacion = None
        self.parametros_actuales = default_params
        self.parametros_comparacion = default_comp_params
        
        # Limpiar paneles de texto
        self.panel_resultados.actualizar_resultados("Interfaz limpiada. Ingrese parámetros y calcule.")
        self.panel_resultados.actualizar_energias("")
        self.panel_resultados.mostrar_consulta("")
        self.panel_resultados.tiempo_consulta.delete(0, tk.END)
        
        # Limpiar gráfico
        self.panel_graficos.actualizar_grafico(None, None)
        
        print("Limpieza completada.")

    def modo_avanzado(self):
        """(Funcionalidad Futura) Cambia a un modo con más opciones."""
        messagebox.showinfo("Información", "Modo Avanzado no implementado aún.")

    def guardar_en_historial(self):
        """Guarda la simulación actual (parámetros y resultados) en el historial."""
        if not self.resultados_actuales:
            return # No guardar si no hay resultados
            
        # Crear entrada de historial
        entrada = {
            'timestamp': time.time(),
            'datos': asdict(self.parametros_actuales),
            'resultados': {
                # Guardar solo resultados clave para evitar archivos grandes
                't_vuelo': self.resultados_actuales.t_vuelo,
                'alcance': self.resultados_actuales.alcance,
                'altura_max': self.resultados_actuales.altura_max
            }
            # Podría añadirse info de comparación si existe
        }
        
        # Añadir al principio de la lista y limitar tamaño (ej: 50 últimas)
        self.historial.insert(0, entrada)
        max_historial = 50
        self.historial = self.historial[:max_historial]
        
        # Actualizar Listbox
        self.panel_historial.actualizar_lista_historial(self.historial)
        
        # Guardar en archivo JSON
        self.guardar_historial_archivo()

    def guardar_historial_archivo(self):
        """Guarda la lista de historial en un archivo JSON."""
        try:
            with open(self.historial_file, 'w', encoding='utf-8') as f:
                json.dump(self.historial, f, indent=2)
        except IOError as e:
            print(f"Error al guardar el historial en '{self.historial_file}': {e}")
        except Exception as e:
             print(f"Error inesperado al guardar historial: {e}")

    def cargar_historial(self):
        """Carga el historial desde el archivo JSON al iniciar."""
        if os.path.exists(self.historial_file):
            try:
                with open(self.historial_file, 'r', encoding='utf-8') as f:
                    self.historial = json.load(f)
                print(f"Historial cargado desde '{self.historial_file}'")
            except json.JSONDecodeError as e:
                print(f"Error al decodificar el archivo de historial '{self.historial_file}': {e}")
                self.historial = []
            except IOError as e:
                 print(f"Error al leer el archivo de historial '{self.historial_file}': {e}")
                 self.historial = []
            except Exception as e:
                 print(f"Error inesperado al cargar historial: {e}")
                 self.historial = []
        else:
            print("No se encontró archivo de historial, iniciando con historial vacío.")
            self.historial = []
            
        # Asegurar que la interfaz se actualice si se carga al inicio (después de crear widgets)
        if hasattr(self, 'panel_historial'):
             self.panel_historial.actualizar_lista_historial(self.historial)

    def cargar_desde_historial(self):
        """Carga los parámetros de una simulación seleccionada del historial."""
        indice = self.panel_historial.obtener_seleccion_indice()
        if indice is None:
            messagebox.showwarning("Sin Selección", "Seleccione una simulación del historial para cargar.")
            return
            
        if 0 <= indice < len(self.historial):
            entrada = self.historial[indice]
            datos_guardados = entrada.get('datos')
            if not datos_guardados:
                 messagebox.showerror("Error", "La entrada de historial seleccionada no contiene datos de parámetros.")
                 return
                 
            print(f"Cargando simulación {indice+1} del historial...")
            # Establecer datos en los paneles
            self.panel_basico.establecer_datos(datos_guardados)
            self.panel_resistencia.establecer_datos(datos_guardados)
            # No cargar en comparación, solo en principal
            self.panel_comparacion.establecer_datos({'comparar': False}) # Desactivar comparación al cargar
            
            # Actualizar estado de paneles dependientes
            self.toggle_resistencia()
            self.toggle_comparacion()
            
            # Recalcular trayectoria con los parámetros cargados
            self.calcular_trayectoria()
            messagebox.showinfo("Historial Cargado", f"Simulación {indice+1} cargada y recalculada.")
        else:
             messagebox.showerror("Error", "Índice de historial seleccionado inválido.")

    def eliminar_del_historial(self):
        """Elimina la simulación seleccionada del historial."""
        indice = self.panel_historial.obtener_seleccion_indice()
        if indice is None:
            messagebox.showwarning("Sin Selección", "Seleccione una simulación del historial para eliminar.")
            return
            
        if 0 <= indice < len(self.historial):
            confirmar = messagebox.askyesno("Confirmar Eliminación", f"¿Está seguro de que desea eliminar la simulación {indice+1} del historial?")
            if confirmar:
                del self.historial[indice]
                self.panel_historial.actualizar_lista_historial(self.historial)
                self.guardar_historial_archivo()
                print(f"Simulación {indice+1} eliminada del historial.")
        else:
             messagebox.showerror("Error", "Índice de historial seleccionado inválido.")

    def limpiar_historial(self):
        """Limpia todo el historial de simulaciones."""
        if not self.historial:
             messagebox.showinfo("Historial Vacío", "El historial ya está vacío.")
             return
             
        confirmar = messagebox.askyesno("Confirmar Limpieza", "¿Está seguro de que desea eliminar TODAS las entradas del historial? Esta acción no se puede deshacer.")
        if confirmar:
            self.historial = []
            self.panel_historial.actualizar_lista_historial(self.historial)
            self.guardar_historial_archivo()
            print("Historial limpiado.")
            
    def toggle_animacion(self):
        """Inicia o detiene la animación de la trayectoria."""
        if self.animacion_activa:
            # Detener animación
            self.panel_graficos.parar_animacion()
            self.animacion_activa = False
        else:
            # Iniciar animación
            if not self.resultados_actuales:
                messagebox.showwarning("Sin Datos", "Calcule una trayectoria antes de iniciar la animación.")
                return
            self.panel_graficos.iniciar_animacion(self.resultados_actuales, self.resultados_comparacion)
            self.animacion_activa = True
            
    def set_velocidad_animacion(self, intervalo_ms: int):
        """Establece la velocidad (intervalo) de la animación."""
        self.velocidad_anim_ms = intervalo_ms
        self.panel_graficos.set_velocidad_animacion(intervalo_ms)

    def on_closing(self):
        """Acciones a realizar al cerrar la ventana."""
        # Detener animación si está activa
        if self.animacion_activa:
            self.panel_graficos.parar_animacion()
        # Preguntar si guardar cambios o historial si es necesario
        # ... (lógica adicional si se requiere)
        print("Cerrando aplicación...")
        self.root.destroy()

# ==============================================
# Punto de entrada principal
# ==============================================

if __name__ == "__main__":
    root = tk.Tk()
    app = TiroParabolicoApp(root)
    # Manejar cierre de ventana
    root.protocol("WM_DELETE_WINDOW", app.on_closing)
    root.mainloop()

