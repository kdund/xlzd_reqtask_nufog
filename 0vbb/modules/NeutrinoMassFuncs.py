#================================NeutrinoMassFuncs.py==================================#
# Created by Yanina Biondi 2019

# This file contains neutrino oscillation parameters from V41 nu-Fit http://www.nu-fit.org/?q=node/256 no SK atmospheric data

# Description of the parameters

'''
    :param m1: Mass of the m1 eigenstate
    :param s12sq: sin^2 theta12
    :param s13sq: sin^2 theta13
    :param s23sq: sin^2 theta23
    :param dm21sq: Delta m21 squared
    :param dm32sq: Delta m32 squared
    :param phi2: First Majorana phase
    :param phi3: Second Majorana phase
    :param hierarchy: Positive for normal, negative for inverted
    :returns: mbb in the same units as m1
    
    suffix 3s stands for 3 sigma, while L and H stands for lower and upper limit.
    suffix IO stands for Inverted Ordering
'''

#==============================================================================#

import numpy as np

class EffectiveMajoranaMass:
    """
    A class that represents the neutrino oscillation parameters
    ...
    Attributes
    ----------
    s12sq : double
        sin^2 theta12
    s13sq : double
        sin^2 theta13
    s23sq : double
        sin^2 theta23
    dm21sq : double
        Delta m21 squared
    dm32sq : double
        Delta m32 squared
        
    Methods
    -------
    
    calculateEffectiveMajoranaMass(m1,  phi2, phi3):
        Compute the effective Majorana mass for given parameters.
    EffectiveMajoranaMassRange(m1, step=0.001):
        Calculates the effective Majorana mass range for a given m1 lightest neutrino mass
    findEffectiveMajoranaMassRegion(region = np.logspace(-5,0,100)  
        Finds the curve of effective Majorana mass for a lightest neutrino mass in region.
        
    """    
    def __init__(self, s12sq, s23sq , s13sq , dm21sq, dm32sq, hierarchy):
        self.hierarchy = hierarchy
        self.s12sq = s12sq
        self.s13sq = s13sq
        self.s23sq = s23sq
        self.dm21sq = dm21sq
        self.dm32sq = dm32sq

    def calculateEffectiveMajoranaMass(self, m1,  phi2, phi3):
        """
        Compute the effective Majorana mass for given parameters.
        Parameters
        ----------
        m1 : double
            lightest neutrino
        phi2: double
            First Majorana phase
        phi3: double
            Second Majorana phase
            
        Returns
        -------
        mbb : double
            Effective Majorana mass
        """
        t12 = np.arcsin(np.sqrt(self.s12sq))
        t13 = np.arcsin(np.sqrt(self.s13sq))
        t23 = np.arcsin(np.sqrt(self.s23sq))

        if self.hierarchy > 0:
            m2 = np.sqrt(np.square(m1) + self.dm21sq)
            m3 = np.sqrt(np.square(m2) + self.dm32sq)
        else:
            # m3 is the lighest!
            m2 = np.sqrt(np.square(m1) + self.dm32sq)
            m3 = np.sqrt(np.square(m2) + self.dm21sq)
            m1, m3 = m3, m1

        mbb = np.abs(np.square(np.cos(t13)) * np.square(np.cos(t12)) * m1 +
                     np.square(np.cos(t13)) * self.s12sq * m2 * np.exp(1j * phi2) +
                     self.s13sq * m3 * np.exp(1j * phi3))

        return mbb
    
    def EffectiveMajoranaMassRange(self, m1, step=0.01):
        """
        Compute range of efective Majorana mass for the lightest neutrino of mass m1.
        Parameters
        ----------
        m1 : double
            lightest neutrino
        step: double, optional
            Granularity of points to scan for Majorana phases
        Returns
        -------
        (np.min(mbbs), np.max(mbbs)) : tuple
            Lowest and highest values of effective Majorana mass given the lightest neutrino mass m1
        
        """
        phis = np.arange(0, 2 * np.pi, step)

        # Create a grid of mbbs with phi2 and phi3 taking on all combinations
        mbbs = self.calculateEffectiveMajoranaMass(m1, phis, phis[:,np.newaxis])

        return (np.min(mbbs), np.max(mbbs))
        
    def findEffectiveMajoranaMassRegion(self, region = np.logspace(-5,0,100)):
        """
        Compute range of efective Majorana mass for the lightest neutrino of mass m1.
        Parameters
        ----------
        region : array of len 100, optional
            Values of lightest neutrino mass to obtain the effective Majorana mass
        Returns
        -------
        g : tuple
            Effective Majorana mass curves of dim (100, 2)
        
        """
        
        # Values of m_lightest to evaluate [eV]
        # Loop over m_lightests to find the mbb range for each
        y = []
        for m in region:
            r = self.EffectiveMajoranaMassRange(m)
            y.append(r)

        # Invert to get a tuple of (lower bounds, upper bounds)
        y = np.array(list(zip(*y)))

        g = np.zeros((2*len(region),2))
        for i in range(len(region)):
            g[i,0] = region[i]
            g[i,1] = y[0][i]
            g[len(region)+i,0] = region[-i-1]
            g[len(region)+i,1] = y[1][-i-1]

        return g
        
#===== Normal Ordering =====#

s12sqNO = 0.303
s23sqNO = 0.572
s13sqNO = 0.02203
dm21sqNO = 7.41e-5
dm32sqNO = 2.511e-3

OscParametersNO = EffectiveMajoranaMass(s12sqNO,s23sqNO,s13sqNO,dm21sqNO,dm32sqNO, 1)

s12sqNO3sL= 0.270
s23sqNO3sL= 0.406
s13sqNO3sL= 0.02029
dm21sqNO3sL= 6.82e-5
dm32sqNO3sL= 2.428e-3

OscParametersNOL = EffectiveMajoranaMass(s12sqNO3sL, s23sqNO3sL, s13sqNO3sL, dm21sqNO3sL, dm32sqNO3sL,1)

s12sqNO3sH = 0.341
s23sqNO3sH = 0.620
s13sqNO3sH = 0.02391
dm21sqNO3sH = 8.03e-5
dm32sqNO3sH = 2.597e-3

OscParametersNOH = EffectiveMajoranaMass(s12sqNO3sH, s23sqNO3sH, s13sqNO3sH, dm21sqNO3sH, dm32sqNO3sH,1)


#===== Inverted Ordering =====#

s12sqIO = 0.303
s23sqIO = 0.578
s13sqIO = 0.02219
dm21sqIO = 7.41e-5
dm32sqIO = 2.498e-3

OscParametersIO = EffectiveMajoranaMass(s12sqIO,s23sqIO,s13sqIO,dm21sqIO,dm32sqIO, -1)

s12sqIO3sL= 0.270
s23sqIO3sL= 0.412
s13sqIO3sL= 0.02047
dm21sqIO3sL= 6.82e-5
dm32sqIO3sL= 2.581e-3

OscParametersIOL = EffectiveMajoranaMass(s12sqIO3sL, s23sqIO3sL, s13sqIO3sL, dm21sqIO3sL, dm32sqIO3sL, -1)

s12sqIO3sH = 0.341
s23sqIO3sH = 0.623
s13sqIO3sH = 0.02396
dm21sqIO3sH = 8.03e-5
dm32sqIO3sH = 2.408e-3

OscParametersIOH = EffectiveMajoranaMass(s12sqIO3sH, s23sqIO3sH, s13sqIO3sH, dm21sqIO3sH, dm32sqIO3sH, -1)