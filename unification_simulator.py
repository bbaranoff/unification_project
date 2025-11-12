# -*- coding: utf-8 -*-
"""
SIMULATEUR DE CHAMP UNIFIÉ - Module principal
Basé sur "Charte des énergies : principes et éthique" - Bastien Baranoff
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

class UnifiedFieldSimulator:
    """
    Simulateur de la connexion universelle et du tenseur de courbure
    selon le formalisme géométrique unifié
    """
    
    def __init__(self, grid_size=50):
        # Constantes fondamentales
        self.G = 6.67430e-11      # Constante gravitationnelle
        self.c = 3e8              # Vitesse de la lumière
        self.hbar = 1.0545718e-34 # Constante de Planck réduite
        
        # Paramètres de simulation
        self.grid_size = grid_size
        self.time_steps = 20
        self.animation_speed = 50
        
        # Initialisation des champs
        self.setup_coordinates()
        
    def setup_coordinates(self):
        """Initialise le système de coordonnées"""
        self.x = np.linspace(-3, 3, self.grid_size)
        self.y = np.linspace(-3, 3, self.grid_size)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.time_values = np.linspace(0, 2*np.pi, self.time_steps)
        
    def spin_connection_component(self, x, y, t):
        """Composante de connexion de spin ω_μ^ab (gravitation)"""
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        
        omega = (np.exp(-r/1.5) * 
                np.cos(2*t - 2*theta) * 
                (1 + 0.3*np.sin(3*theta)))
        
        return omega
    
    def tetrad_component(self, x, y, t):
        """Composante tétrade e_μ^a (géométrie locale)"""
        e_component = (0.5 * np.sin(2*np.pi*x) * 
                      np.cos(2*np.pi*y) * 
                      np.sin(t) * 
                      np.exp(-(x**2 + y**2)/8))
        return e_component
    
    def gauge_connection_component(self, x, y, t, gauge_type='electromagnetic'):
        """Composante de connexion de jauge A_μ^I"""
        if gauge_type == 'electromagnetic':
            return (np.sin(3*x + t) * np.cos(2*y - t) * 
                   np.exp(-(x**2 + y**2)/6))
        elif gauge_type == 'weak':
            return (0.7 * np.cos(2*x - t) * np.sin(3*y + t) *
                   np.exp(-(x**2 + y**2)/5))
        elif gauge_type == 'strong':
            return (0.5 * np.sin(4*x + 0.5*t) * np.cos(3*y - 0.5*t) *
                   np.exp(-(x**2 + y**2)/4))
    
    def universal_connection_A_mu(self, x, y, t):
        """
        Connexion universelle complète :
        A_μ = ½ω_μ^ab J_ab + 1/ℓ e_μ^a P_a + A_μ^I T_I
        """
        omega_component = self.spin_connection_component(x, y, t)
        tetrad_component = self.tetrad_component(x, y, t)
        em_component = self.gauge_connection_component(x, y, t, 'electromagnetic')
        weak_component = self.gauge_connection_component(x, y, t, 'weak')
        strong_component = self.gauge_connection_component(x, y, t, 'strong')
        
        A_mu = (0.4 * omega_component + 
                0.3 * tetrad_component + 
                0.2 * em_component + 
                0.1 * weak_component + 
                0.1 * strong_component)
        
        return A_mu
    
    def curvature_tensor_F_munu(self, x, y, t):
        """Tenseur de courbure F_μν = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν]"""
        A_mu = self.universal_connection_A_mu(x, y, t)
        
        dx = self.x[1] - self.x[0]
        dy = self.y[1] - self.y[0]
        
        dA_dx = np.gradient(A_mu, dx, axis=1)
        dA_dy = np.gradient(A_mu, dy, axis=0)
        
        commutator = A_mu * np.roll(A_mu, 1, axis=0) - np.roll(A_mu, 1, axis=1) * A_mu
        
        F_munu = dA_dx - dA_dy + 0.1 * commutator
        
        return np.abs(F_munu)
    
    def energy_density(self, x, y, t):
        """Densité d'énergie associée au champ ∝ Tr(F_μν F^μν)"""
        F_munu = self.curvature_tensor_F_munu(x, y, t)
        energy = F_munu**2
        return energy
    
    def get_field_components(self, t=0):
        """Retourne tous les composants de champ pour un temps donné"""
        return {
            'A_mu': self.universal_connection_A_mu(self.X, self.Y, t),
            'F_munu': self.curvature_tensor_F_munu(self.X, self.Y, t),
            'energy': self.energy_density(self.X, self.Y, t),
            'omega': self.spin_connection_component(self.X, self.Y, t),
            'em_field': self.gauge_connection_component(self.X, self.Y, t, 'electromagnetic'),
            'X': self.X,
            'Y': self.Y,
            'time': t
        }
    
    def calculate_field_invariants(self, t=0):
        """Calcule les invariants du champ"""
        fields = self.get_field_components(t)
        
        invariants = {
            'A_mu_max': np.max(fields['A_mu']),
            'A_mu_min': np.min(fields['A_mu']),
            'A_mu_mean': np.mean(fields['A_mu']),
            'F_munu_max': np.max(fields['F_munu']),
            'F_munu_min': np.min(fields['F_munu']),
            'F_munu_mean': np.mean(fields['F_munu']),
            'energy_total': np.sum(fields['energy'])
        }
        
        return invariants
