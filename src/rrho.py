
import numpy as np
from scipy.constants import R, h, c, Boltzmann, N_A
from scipy.constants import physical_constants


c = c*100 # convert speed of light in cm/s
J_TO_H = physical_constants['joule-hartree relationship'][0]
AMU_TO_KG = physical_constants['atomic mass constant'][0]


# R = 8.314462618 J/(mol K)
# h = 6.62607015e-34 J*s
# c = 2.9979245800E+10 cm/s
# Boltzmann = 1.380649e-23 J/K
# J_TO_H = 2.2937122783963e+17 Eh/J
# AMU_TO_KG = 1.6605390666e-27 kg*mol/g

def calc_damp(frequency: np.ndarray, cut_off: float, alpha: int) -> np.ndarray:
    """
    Damping factor proportionate to frequency
    1/(1+(cut_off/ν)^α)

    Damping factory NO measure unit
    """
    return 1/(1+(cut_off/frequency)**alpha)



def calc_zpe(frequency: np.ndarray = np.array([0])) -> float: 
    """
    Calculate the Zero Point Energy

    ZPE = Σ(1/2 * (hνc) )
    Return in Eh
    """

    return np.sum( (h*frequency*c)/(2) ) * J_TO_H



def calc_translational_energy(T:float) -> float:
    """
    Translational energy 
    U_trans = 3/2 KbT
    Return in Eh
    """

    return 1.5*Boltzmann*T * J_TO_H

def calc_rotational_energy(T:float, linear=False) -> float:
    """
    Translational energy 
    U_rot = 3/2 KbT or KbT if linear
    Return in Eh
    """
    if linear: return Boltzmann*T * J_TO_H
    return 1.5*Boltzmann*T * J_TO_H


def calc_qRRHO_energy(freq: np.ndarray, T: float) -> np.ndarray:
    """
    quasi-Rigid Rotor Harmonic Oscillator energy

    U = hvc*e^(-hvc/kbT)/(1-e^(-hvc/kbT))
    Return vibrational energy for each vibrational mode in Joule
    """
    f = h*freq*c/(Boltzmann*T)
    return (h*freq*c * np.exp(-f)/(1-np.exp(-f)))




def calc_vibrational_energy(freq : np.ndarray, T : float, cut_off : float, alpha : int) -> float:
    """
    Damp the energy in proportion to the frequency
    
    freq | array : frequency
    T | float : temperature
    alpha | int : damping factor. Default and unchangeable value is 4.
    cut_off | float : damping frequency. Default 100 cm-1

    Return the energy in Eh
    """
    h_damp = calc_damp(freq, cut_off=cut_off, alpha=alpha) 
    return np.sum(h_damp*calc_qRRHO_energy(freq, T) + (1-h_damp) * Boltzmann * T * 0.5) * J_TO_H

def calc_translational_entropy(MW: float, T: float, P: float) -> float:
    """
    Translational entropy

    MW float : molecular weight
    T float  : temperature
    P float  : pressure [Pa]
    solv str : solvent string
    
    S_trans = kB(5/2 + ln(sqrt((2π*MW*kBT)/(N_A*h^2))^3*kB*T/p))
    Return in Eh
    """


    lambda_ = np.sqrt((2* np.pi * MW * Boltzmann * T)/(1000 * N_A * h**2))
    V = (Boltzmann*T)/(P*1000)

    return Boltzmann*(5/2 + np.log(lambda_**3 * V)) * J_TO_H


def calc_rotational_entropy(B, T, zpve, symno: int = 1, linear: bool=False) -> float:
    """
    Rotational entropy

    B: rotational constant [cm-1]
    zpve: to ensure that is not an atom
    linear: if molecule is linear
    symno (σ): numer of symmetry, in relation of the Point Group of the molecule

    θ_R = hcB/kB
    qrot = sqrt(πT^3/(θ_Rx*θ_Ry*θ_Rz))
    S_R = kB * (ln(qrot/σ) + 1.5)

    Return in Eh
    """
    rot_temperature = h*c*B/Boltzmann

    if zpve == 0: 
        return 0
    if linear: 
        qrot = T / rot_temperature
    else:
        qrot = np.sqrt(np.pi * T**3 / np.prod(rot_temperature))

    return  Boltzmann * (np.log(qrot/symno) + 1 + (0 if linear else .5)) * J_TO_H



def calc_S_V_grimme(freq:np.array, T) -> np.array:
    """
    V factor used for the damping of the frequency

    freq: frequencies [cm-1]
    T: temperature
    
    Return in J
    """
    f = h*freq*c/(Boltzmann*T) 
    return (f*Boltzmann)/(np.exp(f)-1) - Boltzmann * np.log(1-np.exp(-f))


def calc_S_R_grimme(freq:np.array, T: float, B:np.array) -> np.array:

    """
    R factor used for the damping of the frequency

    freq: frequencies [cm-1]
    T: temperature
    B: rotatory constant [cm-1]
    
    μ = h/8π^2*ν*c
    (1+ln(8*π^3*(μB)/(µ+B)*kB*T)/h^2) * kB * 0.5

    Return in J
    """

    B = (np.sum(B*c)/len(B))**-1 * h
    mu = h/(8*np.pi**2*freq*c)    
    f = 8*np.pi**3*(mu*B/(mu+B))*Boltzmann*T / h**2

    return (0.5+np.log(f**0.5)) * Boltzmann


def calc_vibrational_entropy(freq, T, B, cut_off=100, alpha=4) -> float:
    """
    Vibrational entropy 

    freq | list : frequencies [cm-1]
    T | float : temperature [K]
    B | np.array : 
    cut_off | float : cut off for the damping of the frequency
    alpha | float : damping factor

    return | float : vibrational entropy [Eh]
    """

    s_damp = calc_damp(freq, cut_off, alpha)
    return np.sum(calc_S_V_grimme(freq, T) * s_damp + (1-s_damp) * calc_S_R_grimme(freq, T, B))*J_TO_H


def calc_electronic_entropy(m) -> float:
    """
    Electronic entropy

    m: electronic multiplicity

    kB*ln(m)

    Return in Eh
    """
    return Boltzmann * np.log(m) * J_TO_H




def free_gibbs_energy(
        SCF : float, T : float, freq : np.ndarray, mw: float, B: np.ndarray, m:int,
        
        # defaults 
        linear : bool =False, cut_off=100, alpha=4,  P: float = 101.325) -> float:
    """
    SCF: self consistent field energy [Eh] + dispersions 
    T: temperature [K]
    P: pressure [kPa]
    B: rotational constant [cm-1]
    m: spin multiplicity

    linear | bool : if molecule is linear
    cut_off | float : frequency cut_off
    alpha | int : frequency damping factor
    P | float : pressure [kPa]
    """
    freq = freq[freq>0]

    zpve = calc_zpe(freq)

    U_trans = calc_translational_energy(T)
    U_rot = calc_rotational_energy(T, linear) if zpve > 0 else 0
    U_vib = calc_vibrational_energy(freq, T, cut_off, alpha)

    H = SCF + zpve + U_trans + U_rot + U_vib + Boltzmann*T*J_TO_H

    S_elec = calc_electronic_entropy(m)
    S_vib = calc_vibrational_entropy(freq, T, B, cut_off, alpha)
    S_rot = calc_rotational_entropy(B, T, zpve)
    S_trans = calc_translational_entropy(mw, T, P)

    S = S_trans + S_rot + S_vib + S_elec

    return H-T*S



