regex_parsing = {
    "orca": {
        "B": "Rotational constants in cm-1",
        "m": "Total Dipole Moment",
        "E": "FINAL SINGLE POINT ENERGY",
        "start_spec": "SPECTRA",
        "end_spec": "***",
        "s_UV": """ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy    Wavelength   fosc        P2         PX        PY        PZ  
        (cm-1)      (nm)                 (au**2)     (au)      (au)      (au) 
-----------------------------------------------------------------------------""",
        "s_ECD": """CD SPECTRUM
-------------------------------------------------------------------
State  Energy     Wavelength     R         MX        MY        MZ   
       (cm-1)       (nm)     (1e40*cgs)   (au)      (au)      (au)  
-------------------------------------------------------------------""",
        "break_spectra": "\n\n",
        "idx_en_UV": 2, # index for energy in the UV table
        "idx_imp_UV": 3, # index for oscillator strength in the UV table
        "idx_en_ECD": 2, # index for energy in the ECD table
        "idx_imp_ECD": 3, # index for rotatory strength in the ECD table
        "s_freq": "VIBRATIONAL FREQUENCIES",
        "e_freq": "------------",
        "idx_freq": 1, # index for frequency in frequency table
        "opt_done" : "                    ***        THE OPTIMIZATION HAS CONVERGED     ***",
        
    }
}
