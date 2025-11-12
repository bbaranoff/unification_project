# -*- coding: utf-8 -*-
"""
Analyse quantique de la gravité : échelles de Planck, fond stochastique et formalisme spectral
Basé sur "Ondes gravitationnelles et gravité quantique" - Bastien Baranoff
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

class QuantumGravityAnalyzer:
    """
    Analyseur quantique-gravitationnel unifiant :
    - Échelles de Planck et quantification
    - Fond stochastique d'ondes GW  
    - Formalisme spectral et vide quantique
    - Analogie vide vibrant
    """
    
    def __init__(self):
        # Constantes fondamentales
        self.h = 6.62607015e-34      # Constante de Planck (J·s)
        self.hbar = 1.054571817e-34  # Constante de Planck réduite (J·s)
        self.c = 299792458           # Vitesse de la lumière (m/s)
        self.G = 6.67430e-11         # Constante gravitationnelle (m³/kg/s²)
        self.kB = 1.380649e-23       # Constante de Boltzmann (J/K)
        
        # Paramètres cosmologiques
        self.H0 = 2.2e-18            # s⁻¹ (H₀ ~ 70 km/s/Mpc)
        self.rho_c = self.critical_density()
        
        # Échelles de Planck
        self.setup_planck_scales()
        
    def setup_planck_scales(self):
        """Calcule les échelles caractéristiques de Planck"""
        # Énergie de Planck
        self.E_planck = np.sqrt(self.hbar * self.c**5 / self.G)
        self.E_planck_GeV = self.E_planck / (1.602176634e-10)  # Conversion en GeV
        
        # Longueur de Planck
        self.l_planck = np.sqrt(self.hbar * self.G / self.c**3)
        
        # Masse de Planck
        self.m_planck = np.sqrt(self.hbar * self.c / self.G)
        
        # Temps de Planck
        self.t_planck = np.sqrt(self.hbar * self.G / self.c**5)
        
        # Fréquence de Planck
        self.nu_planck = 1 / (2 * np.pi * self.t_planck)
    
    def critical_density(self):
        """Densité critique de l'univers ρ_c = 3H₀²/8πG"""
        return (3 * self.H0**2) / (8 * np.pi * self.G)
    
    def planck_einstein_relation(self, nu):
        """
        Relation de Planck-Einstein : E = hν
        Appliquée aux gravitons hypothétiques
        """
        return self.h * nu
    
    def gravitational_wave_energy_density(self, f, S_h):
        """
        Densité d'énergie des ondes gravitationnelles
        Ω_GW(f) = (2π²/3H₀²) f³ S_h(f)
        """
        return (2 * np.pi**2 / (3 * self.H0**2)) * f**3 * S_h
    
    def stochastic_background_estimation(self, f_range=None):
        """
        Estimation du fond stochastique d'ondes gravitationnelles
        """
        if f_range is None:
            f_range = np.logspace(-18, 3, 1000)  # De Hubble à kHz
            
        # Différents modèles de fond stochastique
        models = {
            'Primordial': self.primordial_background_model(f_range),
            'Binary Black Holes': self.bbh_background_model(f_range),
            'Inflationary': self.inflationary_background_model(f_range)
        }
        
        return f_range, models
    
    def primordial_background_model(self, f):
        """Modèle de fond primordial"""
        return 1e-15 * (f/1e-9)**(-1)
    
    def bbh_background_model(self, f):
        """Modèle de fond des binaires de trous noirs"""
        return 1e-9 * (f/10)**(2/3)
    
    def inflationary_background_model(self, f):
        """Modèle de fond inflationnaire"""
        return 1e-17 * np.ones_like(f)
    
    def quantum_gravity_threshold(self, h_amplitude):
        """
        Calcule le seuil de quantification quantique
        où chaque mode porterait un seul quantum d'énergie
        """
        # Pour un détecteur typique (LIGO) à 100 Hz
        f_gw = 100  # Hz
        E_quantum = self.planck_einstein_relation(f_gw)
        
        # Énergie GW typique pour une amplitude h
        E_gw = (self.c**3 / (16 * np.pi * self.G)) * (2 * np.pi * f_gw * h_amplitude)**2
        
        # Nombre de quanta
        n_quanta = E_gw / E_quantum
        
        return {
            'E_quantum': E_quantum,
            'E_gw': E_gw,
            'n_quanta': n_quanta,
            'threshold_reached': n_quanta < 1
        }
    
    def vibrating_vacuum_analogy(self, f):
        """
        Analogie avec le 'vide vibrant' :
        Les fluctuations métriques comme modes de vibration du vide
        CORRECTION : Gestion correcte des tableaux numpy
        """
        # Convertir en tableau numpy si nécessaire
        f_array = np.asarray(f)
        
        # Énergie de point zéro pour un mode de fréquence f
        E_zero_point = 0.5 * self.h * f_array
        
        # Amplitude quantique caractéristique - CORRECTION ICI
        # Éviter la division par zéro pour f=0
        f_safe = np.where(f_array == 0, 1e-100, f_array)
        
        # Calcul sécurisé de h_quantum
        numerator = 16 * np.pi * self.G * E_zero_point
        denominator = self.c**3 * f_safe**2
        h_quantum = np.sqrt(numerator / denominator)
        
        return {
            'frequencies': f_array,
            'zero_point_energy': E_zero_point,
            'quantum_amplitude': h_quantum
        }
    
    def spectral_formalism_connection(self, alpha=0.1):
        """
        Connexion avec le formalisme spectral :
        p_A = α (c²/G) H₀²
        """
        p_A = alpha * (self.c**2 / self.G) * self.H0**2
        Omega_Lambda = 8 * np.pi * alpha / 3
        
        return {
            'p_A': p_A,
            'Omega_Lambda': Omega_Lambda,
            'alpha': alpha,
            'description': 'Formalisme spectral: énergie du vide comme cohérence IR'
        }
    
    def analyze_quantum_gravity_context(self, measured_h=None, measured_f=None):
        """
        Analyse complète du contexte quantique-gravitationnel
        """
        print("\n🌌 CONTEXTE QUANTIQUE-GRAVITATIONNEL")
        print("=" * 60)
        
        # 1. Échelles de Planck
        print("1. 📏 ÉCHELLES DE PLANCK:")
        print(f"   Énergie: {self.E_planck:.2e} J = {self.E_planck_GeV:.2e} GeV")
        print(f"   Longueur: {self.l_planck:.2e} m")
        print(f"   Masse: {self.m_planck:.2e} kg")
        print(f"   Temps: {self.t_planck:.2e} s")
        print(f"   Fréquence: {self.nu_planck:.2e} Hz")
        
        # 2. Seuil de quantification
        quantum_analysis = None
        if measured_h is not None:
            quantum_analysis = self.quantum_gravity_threshold(measured_h)
            print(f"\n2. ⚛️  SEUIL DE QUANTIFICATION:")
            print(f"   Amplitude mesurée: {measured_h:.1e}")
            print(f"   Énergie quantum: {quantum_analysis['E_quantum']:.2e} J")
            print(f"   Énergie GW: {quantum_analysis['E_gw']:.2e} J")
            print(f"   Nombre de quanta: {quantum_analysis['n_quanta']:.2e}")
            
            if quantum_analysis['threshold_reached']:
                print("   → Régime quantique atteint!")
            else:
                print("   → Régime classique (n_quanta >> 1)")
        
        # 3. Fond stochastique
        print(f"\n3. 🌊 FOND STOCHASTIQUE:")
        f_range, bg_models = self.stochastic_background_estimation()
        for model_name, model_values in bg_models.items():
            typical_value = model_values[len(model_values)//2]
            print(f"   {model_name}: Ω_GW ~ {typical_value:.1e}")
        
        # 4. Formalisme spectral
        spectral = self.spectral_formalism_connection()
        print(f"\n4. 📊 FORMALISME SPECTRAL:")
        print(f"   α = {spectral['alpha']}")
        print(f"   Ω_Λ = {spectral['Omega_Lambda']:.3f}")
        print(f"   ρ_vide = {spectral['p_A']:.2e} J/m³")
        
        return {
            'planck_scales': {
                'energy': self.E_planck,
                'length': self.l_planck,
                'mass': self.m_planck,
                'time': self.t_planck,
                'frequency': self.nu_planck
            },
            'quantum_analysis': quantum_analysis,
            'stochastic_background': bg_models,
            'spectral_formalism': spectral
        }
    
    def plot_comprehensive_analysis(self, analysis_results):
        """
        Graphiques complets de l'analyse quantique-gravitationnelle
        CORRECTION : Gestion robuste des calculs
        """
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('ANALYSE QUANTIQUE-GRAVITATIONNELLE COMPLÈTE', 
                    fontsize=16, fontweight='bold')
        
        # 1. Échelles de Planck
        scales = analysis_results['planck_scales']
        scale_names = ['Énergie (J)', 'Longueur (m)', 'Masse (kg)', 'Temps (s)']
        scale_values = [scales['energy'], scales['length'], scales['mass'], scales['time']]
        
        # Conversion logarithmique sécurisée
        log_values = [np.log10(abs(v)) if v != 0 else -100 for v in scale_values]
        
        bars = axes[0,0].bar(scale_names, log_values, color=['red', 'blue', 'green', 'purple'])
        axes[0,0].set_title('Échelles de Planck (log₁₀)')
        axes[0,0].set_ylabel('log₁₀(valeur)')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Ajouter les valeurs sur les barres
        for bar, value in zip(bars, scale_values):
            height = bar.get_height()
            axes[0,0].text(bar.get_x() + bar.get_width()/2., height,
                          f'{value:.1e}', ha='center', va='bottom', rotation=0, fontsize=8)
        
        # 2. Fond stochastique
        f_range = np.logspace(-18, 3, 1000)
        bg_models = analysis_results['stochastic_background']
        
        for model_name, model_values in bg_models.items():
            axes[0,1].loglog(f_range, model_values, label=model_name, linewidth=2)
        
        axes[0,1].set_xlabel('Fréquence (Hz)')
        axes[0,1].set_ylabel('Ω_GW(f)')
        axes[0,1].set_title('Fond stochastique d\'ondes GW')
        axes[0,1].legend()
        axes[0,1].grid(True, alpha=0.3)
        
        # 3. Analogie vide vibrant - CORRECTION COMPLÈTE
        f_vac = np.logspace(-18, 10, 1000)
        
        try:
            vacuum_analogy = self.vibrating_vacuum_analogy(f_vac)
            valid_mask = ~np.isnan(vacuum_analogy['quantum_amplitude']) & ~np.isinf(vacuum_analogy['quantum_amplitude'])
            f_valid = f_vac[valid_mask]
            h_valid = vacuum_analogy['quantum_amplitude'][valid_mask]
            
            if len(f_valid) > 0:
                axes[0,2].loglog(f_valid, h_valid, color='orange', linewidth=2)
                axes[0,2].axhline(1e-21, color='red', linestyle='--', 
                                 label='Sensibilité LIGO')
                axes[0,2].set_xlabel('Fréquence (Hz)')
                axes[0,2].set_ylabel('Amplitude quantique h_quantum')
                axes[0,2].set_title('Analogique du vide vibrant')
                axes[0,2].legend()
                axes[0,2].grid(True, alpha=0.3)
            else:
                axes[0,2].text(0.5, 0.5, 'Données non disponibles', 
                              transform=axes[0,2].transAxes, ha='center')
                axes[0,2].set_title('Analogique du vide vibrant')
        except Exception as e:
            axes[0,2].text(0.5, 0.5, f'Erreur: {str(e)}', 
                          transform=axes[0,2].transAxes, ha='center')
            axes[0,2].set_title('Analogique du vide vibrant')
        
        # 4. Seuil de quantification
        if analysis_results['quantum_analysis']:
            quantum = analysis_results['quantum_analysis']
            labels = ['Énergie quantum', 'Énergie GW']
            values = [quantum['E_quantum'], quantum['E_gw']]
            
            # Conversion logarithmique sécurisée
            log_energy_values = [np.log10(abs(v)) if v != 0 else -100 for v in values]
            
            bars = axes[1,0].bar(labels, log_energy_values, color=['blue', 'red'])
            axes[1,0].set_title('Comparaison énergies (log₁₀)')
            axes[1,0].set_ylabel('log₁₀(Énergie en J)')
            
            # Ajouter les valeurs sur les barres
            for bar, value in zip(bars, values):
                height = bar.get_height()
                axes[1,0].text(bar.get_x() + bar.get_width()/2., height,
                              f'{value:.1e}', ha='center', va='bottom', rotation=0, fontsize=8)
        else:
            axes[1,0].text(0.5, 0.5, 'Aucune donnée GW\npour l\'analyse quantique', 
                          transform=axes[1,0].transAxes, ha='center', va='center')
            axes[1,0].set_title('Seuil de quantification')
        
        # 5. Formalisme spectral
        alpha_range = np.linspace(0.05, 0.15, 50)
        omega_lambda_range = 8 * np.pi * alpha_range / 3
        
        axes[1,1].plot(alpha_range, omega_lambda_range, 'purple', linewidth=2)
        axes[1,1].axvline(0.083, color='red', linestyle='--', label='α optimal')
        axes[1,1].axhline(0.7, color='green', linestyle='--', label='Ω_Λ observé')
        axes[1,1].set_xlabel('Paramètre spectral α')
        axes[1,1].set_ylabel('Ω_Λ')
        axes[1,1].set_title('Formalisme spectral: Ω_Λ(α)')
        axes[1,1].legend()
        axes[1,1].grid(True, alpha=0.3)
        
        # 6. Synthèse conceptuelle
        axes[1,2].axis('off')
        text_str = """
        SYNTHÈSE CONCEPTUELLE
        ───────────────────
        • Énergie Planck: 10¹⁹ GeV
        • Longueur Planck: 10⁻³⁵ m
        • Seuil quantification: ~10⁻³⁵ m
        • Formalisme spectral: Ω_Λ = 8πα/3
        • Vide vibrant: E = ½hν
        • Fond stochastique: trace primordiale
        ───────────────────
        Défi: Atteindre le régime
        quantique de la gravitation
        """
        axes[1,2].text(0.1, 0.9, text_str, transform=axes[1,2].transAxes, 
                      fontfamily='monospace', verticalalignment='top',
                      bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
        
        plt.tight_layout()
        return fig
