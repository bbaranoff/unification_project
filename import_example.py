#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exemple d'import et d'utilisation du simulateur
"""

import numpy as np
from unification_simulator import UnifiedFieldSimulator
from gw_calculator import GWEnergyCalculator
from visualization import FieldVisualizer

# Exemple 1: Utilisation simple
def simple_example():
    """Exemple d'utilisation basique"""
    print("=== EXEMPLE SIMPLE ===")
    
    simulator = UnifiedFieldSimulator(grid_size=30)
    
    # Obtenir les champs à un instant donné
    fields = simulator.get_field_components(t=1.0)
    print(f"Amplitude max de A_μ: {np.max(fields['A_mu']):.4f}")
    print(f"Courbure moyenne F_μν: {np.mean(fields['F_munu']):.4f}")
    
    return fields

# Exemple 2: Calcul des paramètres GW
def gw_analysis_example():
    """Exemple d'analyse des ondes gravitationnelles"""
    print("\n=== ANALYSE ONDES GRAVITATIONNELLES ===")
    
    simulator = UnifiedFieldSimulator()
    gw_calc = GWEnergyCalculator(simulator)
    
    # Analyser sur une plage de temps
    results = gw_calc.analyze_gw_parameters(t_values=np.linspace(0, np.pi, 20))
    
    return results

# Exemple 3: Visualisation spécifique
def visualization_example():
    """Exemple de visualisation"""
    print("\n=== VISUALISATION ===")
    
    simulator = UnifiedFieldSimulator()
    visualizer = FieldVisualizer(simulator)
    
    # Créer une visualisation à un instant spécifique
    fig = visualizer.create_static_visualization(t=2.0)
    
    return fig

if __name__ == "__main__":
    # Exécuter les exemples
    fields = simple_example()
    gw_results = gw_analysis_example()
    fig = visualization_example()
    
    print("\n✅ Tous les exemples exécutés avec succès!")
