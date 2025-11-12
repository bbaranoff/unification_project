# -*- coding: utf-8 -*-
"""
Module de visualisation pour le simulateur de champ unifié
"""

import numpy as np  # IMPORT AJOUTÉ
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class FieldVisualizer:
    """Classe pour visualiser les champs unifiés"""
    
    def __init__(self, simulator):
        self.simulator = simulator
    
    def create_static_visualization(self, t=0):
        """Crée une visualisation statique à un instant t donné"""
        fields = self.simulator.get_field_components(t)
        
        fig = plt.figure(figsize=(20, 12))
        fig.suptitle(f'CONNEXION UNIVERSELLE - Théorie Géométrique Unifiée\nInstant t = {t:.2f}', 
                    fontsize=16, fontweight='bold', y=0.95)
        
        # 1. Connexion universelle A_μ
        ax1 = fig.add_subplot(231, projection='3d')
        surf1 = ax1.plot_surface(fields['X'], fields['Y'], fields['A_mu'], 
                               cmap='viridis', alpha=0.8)
        ax1.set_title('Connexion Universelle A_μ', fontsize=10, pad=20)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('A_μ')
        
        # 2. Courbure gravitationnelle ω_μ^ab
        ax2 = fig.add_subplot(232, projection='3d')
        surf2 = ax2.plot_surface(fields['X'], fields['Y'], fields['omega'], 
                               cmap='plasma', alpha=0.8)
        ax2.set_title('Connexion de Spin ω_μ^ab', fontsize=10, pad=20)
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_zlabel('ω_μ^ab')
        
        # 3. Champ électromagnétique
        ax3 = fig.add_subplot(233, projection='3d')
        surf3 = ax3.plot_surface(fields['X'], fields['Y'], fields['em_field'], 
                               cmap='coolwarm', alpha=0.8)
        ax3.set_title('Connexion de Jauge A_μ^I', fontsize=10, pad=20)
        ax3.set_xlabel('X')
        ax3.set_ylabel('Y')
        ax3.set_zlabel('A_μ^I')
        
        # 4. Tenseur de courbure F_μν
        ax4 = fig.add_subplot(234, projection='3d')
        surf4 = ax4.plot_surface(fields['X'], fields['Y'], fields['F_munu'], 
                               cmap='hot', alpha=0.8)
        ax4.set_title('Tenseur de Courbure F_μν', fontsize=10, pad=20)
        ax4.set_xlabel('X')
        ax4.set_ylabel('Y')
        ax4.set_zlabel('|F_μν|')
        
        # 5. Densité d'énergie
        ax5 = fig.add_subplot(235, projection='3d')
        surf5 = ax5.plot_surface(fields['X'], fields['Y'], fields['energy'], 
                               cmap='inferno', alpha=0.8)
        ax5.set_title('Densité d\'Énergie', fontsize=10, pad=20)
        ax5.set_xlabel('X')
        ax5.set_ylabel('Y')
        ax5.set_zlabel('Énergie')
        
        # 6. Diagramme de phase
        ax6 = fig.add_subplot(236)
        contour = ax6.contourf(fields['X'], fields['Y'], fields['A_mu'], 
                              levels=20, cmap='RdYlBu')
        ax6.streamplot(fields['X'], fields['Y'], 
                      np.gradient(fields['A_mu'], axis=1), 
                      np.gradient(fields['A_mu'], axis=0), 
                      color='black', linewidth=0.5, density=1.5)
        ax6.set_title('Lignes de Champ - Projection 2D', fontsize=10, pad=20)
        ax6.set_xlabel('X')
        ax6.set_ylabel('Y')
        plt.colorbar(contour, ax=ax6)
        
        plt.tight_layout()
        return fig
    
    def create_interactive_plotly(self, t=0):
        """Crée une visualisation interactive avec Plotly"""
        fields = self.simulator.get_field_components(t)
        
        fig = make_subplots(
            rows=2, cols=3,
            specs=[[{'is_3d': True}, {'is_3d': True}, {'is_3d': True}],
                   [{'is_3d': True}, {'is_3d': True}, {'type': 'contour'}]],
            subplot_titles=(
                'Connexion Universelle A_μ',
                'Courbure Gravitationnelle ω_μ^ab', 
                'Champ Électromagnétique A_μ^I',
                'Tenseur de Courbure F_μν',
                'Densité d\'Énergie',
                'Lignes de Champ 2D'
            ),
            print_grid=False
        )
        
        fig.add_trace(go.Surface(x=fields['X'], y=fields['Y'], z=fields['A_mu'], 
                               colorscale='Viridis', name='A_μ'), row=1, col=1)
        fig.add_trace(go.Surface(x=fields['X'], y=fields['Y'], z=fields['omega'],
                               colorscale='Plasma', name='ω_μ^ab'), row=1, col=2)
        fig.add_trace(go.Surface(x=fields['X'], y=fields['Y'], z=fields['em_field'],
                               colorscale='RdBu', name='A_μ^I'), row=1, col=3)
        fig.add_trace(go.Surface(x=fields['X'], y=fields['Y'], z=fields['F_munu'],
                               colorscale='Hot', name='F_μν'), row=2, col=1)
        fig.add_trace(go.Surface(x=fields['X'], y=fields['Y'], z=fields['energy'],
                               colorscale='Inferno', name='Énergie'), row=2, col=2)
        fig.add_trace(go.Contour(x=self.simulator.x, y=self.simulator.y, z=fields['A_mu'],
                               colorscale='RdYlBu', name='Lignes de champ'), row=2, col=3)
        
        fig.update_layout(
            title_text="SIMULATEUR DE CHAMP UNIFIÉ - Visualisation Interactive",
            title_x=0.5,
            height=800
        )
        
        return fig
