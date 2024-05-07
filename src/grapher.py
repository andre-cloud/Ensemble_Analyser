import numpy as np
import os
import scipy.optimize as opt
from scipy.signal import argrelextrema
from scipy.constants import c, h, electron_volt, R
import matplotlib.pyplot as plt
from typing import Union

plt.set_loglevel("error")

try:
    from src.regex_parsing import regex_parsing
except ImportError as e:  # pragma: no cover
    print(e)
    from regex_parsing import regex_parsing

FACTOR_EV_NM = h * c / (10**-9 * electron_volt)


class Graph:
    """
    Generate electronic graphs (UV and ECD) and, if present "uv_ref.dat" or "ecd_ref.dat" file as mere xy file, it will automatically align the calculated graph with the reference.
    """

    def __init__(
        self,
        confs,
        protocol,
        log,
        T,
        final_lambda=800.0,
        definition=4,
        FWHM: Union[None, float] = None,
        shift: Union[None, float] = None,
        invert: bool = False,
        regraph: bool = False,
    ):
        self.confs = [i for i in confs if i.active]
        self.protocol = protocol
        self.log = log
        self.pop = self.calc_pop(T, regraph=regraph)

        self.log.debug(self.pop)

        self.x = np.linspace(
            FACTOR_EV_NM / 150, FACTOR_EV_NM / (final_lambda), 10**definition
        )  # eV x axis

        self.spectra = []
        self.filter_outputs()

        self.uv_impulses = [self.get_uv(i, pop) for i, pop in zip(self.spectra, self.pop)]
        self.ecd_impulses = [self.get_ecd(i, pop) for i, pop in zip(self.spectra, self.pop)]

        # creating the ECD and UV graph. If uv_ref.dat and/or ecd_ref.dat auto-convolution of the reference is performed
        if os.path.exists(os.path.join(os.getcwd(), "ecd_ref.dat")):
            ecd = self.auto_convolution(
                os.path.join(os.getcwd(), "ecd_ref.dat"),
                fname_ref_damp=os.path.join(os.getcwd(), "ecd_ref_norm_eV.dat"),
                impulses=self.ecd_impulses,
                fname=f"ecd_protocol_{self.protocol.number}_auto_conv.dat",
                user_sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else None,
                user_shift=shift,
                invert=invert,
            )
        else:
            ecd = self.calc_graph(
                impulses=self.ecd_impulses,
                sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else 1 / 3,
                fname=f"ecd_protocol_{self.protocol.number}.dat",
                save=True,
            )

        if os.path.exists(os.path.join(os.getcwd(), "uv_ref.dat")):
            uv = self.auto_convolution(
                os.path.join(os.getcwd(), "uv_ref.dat"),
                fname_ref_damp=os.path.join(os.getcwd(), "uv_ref_norm_eV.dat"),
                impulses=self.uv_impulses,
                fname=f"uv_protocol_{self.protocol.number}_auto_conv.dat",
                user_sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else None,
                user_shift=shift,
            )
        else:
            uv = self.calc_graph(
                impulses=self.uv_impulses,
                sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else 1 / 3,
                fname=f"uv_protocol_{self.protocol.number}.dat",
                save=True,
            )

        Graph.damp_graph(f"ecd_protocol_{self.protocol.number}.dat", self.x, ecd)
        Graph.damp_graph(f"uv_protocol_{self.protocol.number}.dat", self.x, uv)

    def filter_outputs(self) -> None:
        """
        Read once all conformers output, keeping only graph information

        :rtype: None
        """

        st, en = (
            regex_parsing[self.protocol.calculator]["start_spec"],
            regex_parsing[self.protocol.calculator]["end_spec"],
        )

        for i in self.confs:
            with open(
                os.path.join(
                    os.getcwd(), i.folder, f"protocol_{self.protocol.number}.out"
                )
            ) as f:
                fl = f.read()

            sp = fl.split(st)[-1].split(en)[0]
            self.spectra.append(sp)

        return None

    def get_uv(self, spectra, pop):
        """
        Get the impulses for the UV spectra calculation

        :param spectra: parsed output file
        :type spectra: str
        :param pop: bolzmann population
        :type pop: float

        :return: energy and impulse tuple
        :rtype: tuple(float, float)
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]["s_UV"])[
            -1
        ].split(regex_parsing[self.protocol.calculator]["break"])[0]

        return [
            (
                FACTOR_EV_NM
                / float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_en_UV"]
                    ]
                ),
                float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_imp_UV"]
                    ]
                ) * pop,
            )
            for i in graph.splitlines()
            if i
        ]

    def get_ecd(self, spectra, pop):
        """
        Get the impulses for the ECD spectra calculation

        :param spectra: parse output file
        :type spectra: str
        :param pop: bolzmann population
        :type pop: float


        :return: energy and impulse tuple
        :rtype: tuple(float, float)
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]["s_ECD"])[
            -1
        ].split(regex_parsing[self.protocol.calculator]["break"])[0]

        return [
            (
                FACTOR_EV_NM
                / float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_en_ECD"]
                    ]
                ),
                float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_imp_ECD"]
                    ]
                ) * pop,
            )
            for i in graph.splitlines()
            if i
        ]

    @staticmethod
    def gaussian(x, ev, I, sigma) -> np.ndarray:
        r"""
        Create a gaussian convolution for each impulse

        :param ev: energy of the impulse
        :type ev: float
        :param I: intensity of the impulse (Fosc for UV, R(vel) for ECD)
        :type I: float
        :param sigma: sigma of the gaussian distribution
        :type sigma: float

        .. math::
            f(x) = \frac {I}{σ*\sqrt{2π}} e^{-\frac 12(\frac {x-eV}{σ})^2}

        :return: Gaussian of the impulse
        :rtype: np.array
        """
        return I / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - ev) / sigma) ** 2)

    def calc_pop(self, T, regraph) -> np.ndarray:
        """
        Boltzmann population, if necessary with correction

        :param T: temperature [K]
        :type T: float
        :param regraph: the routine ragraph is running
        :type regraph: bool

        :return: population distribution
        :rtype: np.array (1D array)
        """

        if not regraph:
            n, n_1 = self.protocol.number, str(int(self.protocol.number) - 1)

            # CONF do have frequency calculation before
            if self.protocol.freq:
                pass
            elif self.confs[0].energies[n_1]["G"]:
                # So energy is corrected: if functionals are the same, nothing change; else energy of the new function is corrected with lower frequency correction
                for i in self.confs:
                    i.energies[n]["G"] = i.energies[n]["E"] + (
                        i.energies[n_1]["G"] - i.energies[n_1]["E"]
                    )
            else:
                raise IOError("No frequency calculation")
        

            ens = np.array([i.get_energy for i in self.confs])
            ens_rel = ens - min(ens)
            bolz = np.exp((-ens_rel * 4186) / (R * T))
            pop = bolz / np.sum(bolz)
            for idx, i in enumerate(list(ens)):
                self.confs[idx]._last_energy["G"] = ens[idx]
                self.confs[idx]._last_energy["Erel"] = i
                self.confs[idx]._last_energy["Pop"] = pop[idx] * 100

        else:
            pop = np.array([i.energies[self.protocol.number]["Pop"] for i in self.confs])


        return pop

    def auto_convolution(
        self,
        fname_ref,
        fname_ref_damp,
        impulses,
        fname,
        norm=1,
        user_sigma=None,
        user_shift=None,
        invert=False,
    ) -> np.ndarray:
        """
        Optimization to find the best fitting values for the Gaussian convolution.
        Optimization "Fitness Function" is the sum of the absolute value of the differences between the computed and the experimental graph that lay above the threshold.

        :param fname_ref: filename of the reference file
        :type fname_ref: str
        :param fname_ref_damp: filename to save the autoconvoluted graph
        :type fname_ref_damp: str
        :param impulses: list of single excitation [eV, I] where I can be a UV or ECD excitation
        :type impulses: np.array((eV, I))
        :param fname: filename to save the final convoluted graph
        :type fname: str
        :param user_sigma: user defined sigma for gaussian convolution
        :type user_sigma: float
        :param user_shift: user defined shift after gaussian convolution
        :type user_shift: float
        :param invert: invert the ECD graph when True
        :type invert: bool

        :return: normalized graph
        :rtype: np.array (1D array)
        """

        ref = Ref_graph(fname_ref, None, invert=invert)
        x_min, x_max = ref.x_min, ref.x_max
        ref.y = Graph.normalise(ref.y, norm=norm)

        # resampling the experimental data, in order to fetch the x_exp.size
        X = self.x.copy()
        Y_exp_interp = np.interp(X, ref.x, ref.y, left=0, right=0)
        
        max_exp, limit = ref.get_maximum()


        if set(list(Y_exp_interp)) == set([0. ]):
            Y_exp_interp = np.interp(X, ref.x[::-1], ref.y[::-1], left=0, right=0)

        Graph.damp_graph(fname_ref_damp, X, Y_exp_interp)

        def optimiser(variables):
            """
            Callback for the scipy.optimize.minimize
            """
            sigma, shift, = variables

            # Here the normalization! 
            Y_comp = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=shift, sigma=sigma, save=False
                ),
                norm=norm, x_min=idx_x_min, x_max=idx_x_max
            )

            # RMSD calculation
            rmsd = np.sqrt(
                np.mean(
                    (
                        Y_comp[np.where((X >= x_min) & (X <= x_max))]
                        - Y_exp_interp[np.where((X >= x_min) & (X <= x_max))]
                    )
                    ** 2
                )
            )
            return rmsd

        initial_guess = [0.1415, -0.000001]
        shift = (user_shift, user_shift) if user_shift else (-1, 1)
        sigma = (user_sigma, user_sigma) if user_sigma else (0.08, 0.21)
        default_guess = [0.1415, 0]  # the σ correspond to a FWHM of 0.33 eV
        result = opt.minimize(
            optimiser,
            initial_guess,
            bounds=[sigma, shift],  # FWHM is between 0.2 and 0.5 eV
            options={"maxiter": 10000},
            method="Powell",
        )
        if result.success:
            sigma, shift = result.x
            self.log.info(
                f"Convergence of parameters succeeded in {result.nfev} steps.\n"
                f"Confidence level: {(1-result.fun/2)*100:.2f}%. Parameters obtained\n"
                f"\t- σ = {sigma:.4f} eV (that correspond to a FWHM = {(sigma*np.sqrt(2*np.log(2))*2):.4f} eV)\n"
                f"\t- Δ = {shift:.4f} eV (in this case, a negative shift corresponds to a RED-shift)"
            )
            Y_COMP = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=shift, sigma=sigma, save=True, fname=fname
                ),
                norm=norm,
            )

        else:
            self.log.info(
                f"Convergence of parameters NOT succeeded.\n"
                "Parameters used to convolute the saved graph\n"
                f"\t- σ = {default_guess[0]:.4f} eV (that correspond to a FWHM = {(default_guess[0]*np.sqrt(2*np.log(2))*2):.4f} eV)\n"
                f"\t- Δ = 0.0000 eV"
            )
            Y_COMP = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=0, sigma=1 / 3, save=True, fname=fname
                ),
                norm=norm,
            )

        return Y_COMP

    def calc_graph(self, impulses, sigma, shift=0, fname="", save=False) -> np.ndarray:
        """
        Build the Spectra

        :param impulses: list of the (eV, I)s for each root of each conformer
        :type impulses: list
        :param sigma: dispersion for the gaussian convolution [eV]
        :type sigma: float
        :param fname: the name of the file to store the graph
        :type fname: str
        :param save: save the .dat file for each conformer
        :type save: bool

        :return: Convoluted and weighted spectra
        :rtype: np.array (1D array)
        """

        x = self.x.copy()
        y = np.zeros(x.shape)

        for idx in range(len(self.confs)):
            y_ = np.zeros(x.shape)
            for ev, I in impulses[idx]:
                y_ += Graph.gaussian(x + shift, ev, I, sigma)

            if save:
                Graph.damp_graph(
                    fname=os.path.join(os.getcwd(), self.confs[idx].folder, fname),
                    x=x,
                    y=y_,
                )

            y += y_

        return y

    @staticmethod
    def damp_graph(fname: str, x: np.ndarray, y: np.ndarray) -> None:
        """
        Damp an xy graph into file

        :param fname: filename to store the graph
        :type fname: str
        :param x: x axis
        :type x: np.array (1D array)
        :param y: y axis
        :type y: np.array (1D array)

        :return: None
        """
        data = np.array([(xi, yi) for xi, yi in zip(x, y)])
        np.savetxt(fname, data)
        return

    @staticmethod
    def load_graph(fname, is_ev=False):
        """
        Load an already saved graph from a file XY

        :param fname: filename
        :type fname: str
        :param is_ev: if the Y is already converted in eV, defaults to False
        :type is_ev: bool, optional
        :return: the graph it self divided in X and Y
        :rtype: tuple (X, Y)
        """
        arr = np.loadtxt(fname)
        if not is_ev:
            return FACTOR_EV_NM / arr[:, 0], arr[:, 1]
        return arr[:, 0], arr[:, 1]

    @staticmethod
    def normalise(y: np.ndarray, norm=1, x_min=0, x_max=-1) -> np.array:
        """
        Normalize an ensemble between 1 and -1, if not set otherwise.

        :param y: 1D array
        :type y: np.array
        :param norm: max value to normalize at
        :type norm: float
        :param x_min: index of the minimum value of the experimental value
        :type x_min: int
        :param x_max: index of the maximum value of the experimental value
        :type x_max: int


        :return: 1D normalized array
        :rtype: np.array
        """
        return (
            y / (np.max([np.max(y[x_min:x_max]), np.min(y[x_min:x_max]) * (-1 if np.min(y[x_min:x_max]) < 0 else 1)])) * norm
        )


class Ref_graph:
    """
    Load the reference graph in order to shift and convolute properly the calculated one.
    """

    def __init__(self, fname: str, log, is_ev: bool = False, invert: bool = False):
        data = np.loadtxt(fname, dtype=float)
        self.data = data[np.argsort(data[:, 0])]
        self.x = data[:, 0] if is_ev else FACTOR_EV_NM / data[:, 0]
        self.y = data[:, 1] * (1 if not invert else -1)

        self.log = log

    @property
    def x_min(self):
        return float(min(self.x))

    @property
    def x_max(self):
        return float(max(self.x))
    

    def get_maximum(self):
        """Get the maximum between the maxima. 

        :return: X,Y of the maximum of the maxima. This point must lay between two minima to be considered
        :rtype: tuple
        """


        max_indices = argrelextrema(self.y, np.greater, order=100)[0]
        min_indices = argrelextrema(self.y, np.less, order=100)[0]
        
        first_nonzero_index = np.argmax(self.y != 0)
        last_nonzero_index = len(self.y) - np.argmax(self.y[::-1] != 0) - 1
        
        if self.y[first_nonzero_index] >= self.y[first_nonzero_index+1]:
            max_indices = np.concatenate((max_indices, [first_nonzero_index])) if len(max_indices)!=0 else np.array([self.y[0]])
        elif self.y[first_nonzero_index] < self.y[first_nonzero_index+1]:
            min_indices = np.concatenate((min_indices, [first_nonzero_index])) if len(min_indices)!=0 else np.array([self.y[0]])

        if self.y[last_nonzero_index] >= self.y[last_nonzero_index-1]:
            max_indices = np.concatenate((max_indices, [last_nonzero_index])) if len(max_indices)!=0 else np.array([last_nonzero_index])
        elif self.y[last_nonzero_index] < self.y[last_nonzero_index-1]:
            min_indices = np.concatenate((min_indices, [last_nonzero_index])) if len(min_indices)!=0 else np.array([last_nonzero_index])
            
        massimi = np.array([(self.x[i], self.y[i]) for i in max_indices])

        minimi = np.array([(self.x[i], self.y[i]) for i in min_indices])
        print(minimi)


        max_with_minima = []
        for i in max_indices:
            left_min_index = np.argmax(min_indices < i)
            right_min_index = np.argmax(min_indices > i)
            if left_min_index >= 0 and right_min_index <= len(min_indices) - 1:
                max_with_minima.append((self.x[i], self.y[i]))
        
        max_with_minima = np.array(max_with_minima)
        tmp = np.argmax(max_with_minima[:, 1])
        max_with_minima = list(max_with_minima[tmp])
        
        limits = [minimi[np.where(minimi[:, 0]<max_with_minima[0])[0][-1]], minimi[np.where(minimi[:, 0]>max_with_minima[0])[0][0]] ]

        return max_with_minima, limits
