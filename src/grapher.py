import numpy as np
import os
import scipy.optimize as opt
from scipy.constants import c, h, electron_volt, R
import matplotlib.pyplot as plt

try:
    from src.regex_parsing import regex_parsing
except ImportError as e:
    print(e)
    from regex_parsing import regex_parsing

FACTOR_EV_NM = h * c / (10**-9 * electron_volt)


class Graph:
    def __init__(self, confs, protocol, log, T, final_lambda=800.0, definition=4):
        """

        confs | list : whole ensemble list
        protocol | Protocol : protocol that calculates the electronic spectra
        log : logger instance
        T | float : temperature [K]
        final_lambda | float : last wavelength to convolute the spectra
        definition | int : number of point for the wavelength interval
        """

        self.confs = [i for i in confs if i.active]
        self.protocol = protocol
        self.log = log
        self.pop = self.calc_pop(T)
        self.log.debug(self.pop)

        self.x = np.linspace(
            FACTOR_EV_NM / 100, FACTOR_EV_NM / (final_lambda), 10**definition
        )  # eV x axis

        self.spectra = []
        self.filter_outputs()

        self.uv_impulses = [self.get_uv(i) for i in self.spectra]
        self.ecd_impulses = [self.get_ecd(i) for i in self.spectra]

        # creating the ECD and UV graph. If uv_ref.dat and/or ecd_ref.dat autoconvolution of the reference is performed
        if os.path.exists(os.path.join(os.getcwd(), "ecd_ref.dat")):
            ecd = self.auto_convolution(
                os.path.join(os.getcwd(), "ecd_ref.dat"),
                fname_ref_damp=os.path.join(os.getcwd(), "ecd_ref_norm_eV.dat"),
                impulses=self.ecd_impulses,
                fname=f"ecd_protocol_{self.protocol.number}_auto_conv.dat",
            )
        else:
            ecd = self.calc_graph(
                impulses=self.ecd_impulses,
                sigma=1 / 3,
                fname=f"ecd_protocol_{self.protocol.number}.dat",
                save=True,
            )

        if os.path.exists(os.path.join(os.getcwd(), "uv_ref.dat")):
            uv = self.auto_convolution(
                os.path.join(os.getcwd(), "uv_ref.dat"),
                fname_ref_damp=os.path.join(os.getcwd(), "uv_ref_norm_eV.dat"),
                impulses=self.ecd_impulses,
                fname=f"uv_protocol_{self.protocol.number}_auto_conv.dat",
            )
        else:
            uv = self.calc_graph(
                impulses=self.ecd_impulses,
                sigma=1 / 3,
                fname=f"uv_protocol_{self.protocol.number}.dat",
                save=True,
            )

        Graph.damp_graph(f"ecd_protocol_{self.protocol.number}.dat", self.x, ecd)
        Graph.damp_graph(f"uv_protocol_{self.protocol.number}.dat", self.x, uv)

    def filter_outputs(self) -> None:
        """
        Read once all conformers output, keeping only graph information

        return None
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

    def get_uv(self, spectra):
        """
        Get the impulses for the UV spectra calculation

        spectra | str : parse output file

        return | tuple(float, float) : energy and impulse tuple
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]["s_UV"])[
            -1
        ].split(regex_parsing[self.protocol.calculator]["break_spectra"])[0]

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
                ),
            )
            for i in graph.splitlines()
            if i
        ]

    def get_ecd(self, spectra):
        """
        Get the impulses for the ECD spectra calculation

        spectra | str : parse output file

        return | tuple(float, float) : energy and impulse tuple
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]["s_ECD"])[
            -1
        ].split(regex_parsing[self.protocol.calculator]["break_spectra"])[0]

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
                ),
            )
            for i in graph.splitlines()
            if i
        ]

    @staticmethod
    def gaussian(x, ev, I, sigma) -> np.ndarray:
        """
        Create a gaussian convolution for each impulse

        ev | float : energy of the impulse
        I | float : intensity of the impulse (Fosc for UV, R(vel) for ECD)
        sigma | float : sigma of the gaussian distribution

        I/(σ*sqrt(2π)) * exp(-1/2*((x-ev)/σ)^2)

        return | np.array : Gaussian of the impulse
        """
        return I / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - ev) / sigma) ** 2)

    def calc_pop(self, T) -> np.ndarray:
        """
        Boltzmann population, if necessary with correction

        T | float : temperature [K]

        return | np.array (1D) : population distribution
        """

        # if not self.protocol.take_from:
        n, n_1 = self.protocol.number, str(int(self.protocol.number) - 1)
        # else:
        #     n , n_1 = self.protocol.number, str(self.protocol.take_from)

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

        return pop

    def auto_convolution(
        self, fname_ref, fname_ref_damp, impulses, fname, norm=1
    ) -> np.ndarray:
        """
        Optimization to find the best fitting values for the Gaussian convolution.
        Optimization "Fitness Function" is the sum of the absolute value of the differences between the computed and the experimental graph that lay above the threshold.

        fname_ref | str : filename of the reference file
        fname_ref_damp | str : filename to save the autoconvoluted graph
        impulses | np.array((eV, I)) : list of single excitation [eV, I] where I can be a UV or ECD excitation
        fname | str : filename to save the final convoluted graph

        return | np.array : normalized graph
        """

        ref = Ref_graph(fname_ref, None)
        x_min, x_max = ref.x_min, ref.x_max
        ref.y = Graph.normalise(ref.y, norm=norm)
        # area_ref = trapezoid(a.y, a.x)

        # resampling the experimental data, in order to fetch the x_exp.size
        X = self.x.copy()
        Y_exp_interp = np.interp(X, ref.x, ref.y, left=0, right=0)
        Graph.damp_graph(fname_ref_damp, X, Y_exp_interp)

        def optimiser(variables):
            """
            Callback for the scipy.optimize.minimize
            """
            sigma, shift, threshold = variables
            Y_comp = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=shift, sigma=sigma, save=False
                ),
                norm=norm,
            )
            # graphs[len(graphs)] = Y_comp

            # different with threshold
            # y = np.abs(Y_comp - Y_exp_interp)
            # y = y[np.argwhere((X>=x_min) & (X<=x_max))]
            # diff = np.sum(y[np.where((y>0) & (y>threshold))])

            # RMSD calculation
            rmsd = np.sqrt(np.mean((Y_comp - Y_exp_interp) ** 2))

            # exp = np.array([[x,y] for x, y in zip(X, Y_exp_interp)])
            # comp = np.array([[x,y] for x, y in zip(X, Y_comp)])

            # diff_area_pos = trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1]) - trapezoid(comp[comp[:, 1]>0][:, 0], comp[comp[:, 1]>0][:, 1])
            # diff_area_neg = trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1]) - trapezoid(comp[comp[:, 1]<0][:, 0], comp[comp[:, 1]<0][:, 1])
            # diff = (diff_area_pos*trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1])) + diff_area_neg*trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1])/ (trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1]) + trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1]) )

            # print(sigma, shift, threshold, diff)
            return rmsd

        confidence = 0.01
        initial_guess = [0.2, -0.000001, confidence]
        self.log.debug(initial_guess)
        default_guess = [0.214, 0]
        result = opt.minimize(
            optimiser,
            initial_guess,
            bounds=[(0.16, 1 / 5), (-1, 1), (0.01, 0.01)],
            options={"maxiter": 10000},
            method="Powell",
        )
        if result.success:
            sigma, shift, thr = result.x
            self.log.info(
                f"Convergence of parameters succeeded within a threshold of {thr:.2f}u.a. for the ∆ε. Confidence level: {(1-result.fun/(2*X[np.where((X>=x_min) & (X<=x_max))].size))*100:.2f}%. Parameters obtained\n\t- σ = {sigma:.4f} eV (that correspond to a FWHM = {(sigma*np.sqrt(2*np.log(2))*2):.4f} eV)\n\t- Δ = {shift:.4f} eV (in this case, a negative shift corresponds to a RED-shift)"
            )
            Y_COMP = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=shift, sigma=sigma, save=True, fname=fname
                ),
                norm=norm,
            )
        else:
            self.log.info(
                f"Convergence of parameters NOT succeeded within a threshold of {confidence:.2f}u.a. for the ∆ε. Parameters used to convolute the saved graph\n\t- σ = {default_guess[0]:.4f} eV (that correspond to a FWHM = {(default_guess[0]*np.sqrt(2*np.log(2))*2):.4f} eV)\n\t- Δ = 0.0000 eV"
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

        impulses | list(tuple) : list of the (eV, I)s for each root of each conformer
        sigma | float : dispersion for the gaussian convolution [eV]
        fname | str : the name of the file to store the graph
        save | bool : save the .dat file for each conformer

        return | np.array (1D array) : Convoluted and weighted spectra
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

            y += self.pop[idx] * y_

        return y

    @staticmethod
    def damp_graph(fname: str, x: np.ndarray, y: np.ndarray) -> None:
        """
        Damp an xy graph into file

        fname | str : filename to store the graph
        x | np.array (1D) : x axis
        y | np.array (1D) : y axis

        return None
        """
        data = np.array([(xi, yi) for xi, yi in zip(x, y)])
        np.savetxt(fname, data)
        return

    @staticmethod
    def load_graph(fname, is_ev=False):
        arr = np.loadtxt(fname)
        if not is_ev:
            return FACTOR_EV_NM / arr[:, 0], arr[:, 1]
        return arr[:, 0], arr[:, 1]

    @staticmethod
    def normalise(y: np.ndarray, norm=1) -> np.array:
        """
        Normalize an ensemble between 1 and -1, if not set otherwise.

        y | np.array : 1D array
        norm | float : max value to normalize at

        return | np.array : 1D normalized array
        """
        return (
            y / (np.max([np.max(y), np.min(y) * (-1 if np.min(y) < 0 else 1)])) * norm
        )


class Ref_graph:
    def __init__(self, fname: str, log, is_ev: bool = False):
        data = np.loadtxt(fname, dtype=float)
        self.x = data[:, 0] if is_ev else FACTOR_EV_NM / data[:, 0]
        self.y = data[:, 1]

        self.log = log

    # def integral(self) -> float:
    #     """
    #     Calculate the area subtended by the curve

    #     return | float
    #     """

    #     return trapezoid(self.y, self.x)

    @property
    def x_min(self):
        return float(min(self.x))

    @property
    def x_max(self):
        return float(max(self.x))


class Test_Graph:
    def __init__(self, fname):
        data = np.loadtxt(fname, dtype=float)
        self.x = np.linspace(
            FACTOR_EV_NM / 100, FACTOR_EV_NM / (800), 10**5
        )  # eV x axis
        self.eV = FACTOR_EV_NM / data[:, 0]
        self.imp = data[:, 1]
        # self.y = Graph.normalise(self.y)

    def calc_graph(self, shift, sigma):
        x = self.x.copy()
        y = np.zeros(x.shape)

        y_ = np.zeros(x.shape)
        for ev, I in zip(self.eV, self.imp):
            y_ += Graph.gaussian(x, ev + shift, I, sigma)
        y += y_

        return y


if __name__ == "__main__":
    import scipy.optimize as opt

    os.chdir("src")

    a = Ref_graph("../files/ecd_ref.txt", None)
    x_min, x_max = a.x_min, a.x_max
    a.y = Graph.normalise(a.y)
    # area_ref = trapezoid(a.y, a.x)

    # resampling the experimental data, in order to fetch the x_exp.size
    X = np.linspace(FACTOR_EV_NM / 100, FACTOR_EV_NM / (800), 10**5)

    Y_exp_interp = np.interp(X, a.x, a.y, left=0, right=0)

    computed = Test_Graph("../files/impulse.dat")

    graphs = {}
    counter = 0

    def optimiser(variables):
        global counter
        sigma, shift, threshold = variables
        Y_comp = Graph.normalise(computed.calc_graph(shift, sigma))
        graphs[len(graphs)] = Y_comp
        # graphs[len(graphs)] = Y_comp

        # different with threshold
        # y = np.abs(Y_comp - Y_exp_interp)
        # y = y[np.argwhere((X>=x_min) & (X<=x_max))]
        # diff = np.sum(y[np.where((y>0) & (y>threshold))])

        # RMSD calculation
        rmsd = np.sqrt(np.mean((Y_comp - Y_exp_interp) ** 2))

        # exp = np.array([[x,y] for x, y in zip(X, Y_exp_interp)])
        # comp = np.array([[x,y] for x, y in zip(X, Y_comp)])

        # diff_area_pos = trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1]) - trapezoid(comp[comp[:, 1]>0][:, 0], comp[comp[:, 1]>0][:, 1])
        # diff_area_neg = trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1]) - trapezoid(comp[comp[:, 1]<0][:, 0], comp[comp[:, 1]<0][:, 1])
        # diff = (diff_area_pos*trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1])) + diff_area_neg*trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1])/ (trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1]) + trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1]) )

        # print(sigma, shift, threshold, diff)
        counter += 1
        return rmsd
        sigma, shift, threshold = variables
        Y_comp = Graph.normalise(computed.calc_graph(shift, sigma))

        # different with threshold
        y = np.abs(Y_comp - Y_exp_interp)
        diff = np.sum(y[np.where((y > 0) & (y > threshold))])

        exp = np.array([[x, y] for x, y in zip(X, Y_exp_interp)])
        comp = np.array([[x, y] for x, y in zip(X, Y_comp)])

        # diff_area_pos = trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1]) - trapezoid(comp[comp[:, 1]>0][:, 0], comp[comp[:, 1]>0][:, 1])
        # diff_area_neg = trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1]) - trapezoid(comp[comp[:, 1]<0][:, 0], comp[comp[:, 1]<0][:, 1])
        # diff = (diff_area_pos*trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1])) + diff_area_neg*trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1])/ (trapezoid(exp[exp[:, 1]>0][:, 0], exp[exp[:, 1]>0][:, 1]) + trapezoid(exp[exp[:, 1]<0][:, 0], exp[exp[:, 1]<0][:, 1]) )

        print(sigma, shift, threshold, diff)
        return diff

    confidence = 0.01
    initial_guess = [0.2, -0.000001, confidence]
    result = opt.minimize(
        optimiser,
        initial_guess,
        bounds=[(0.16, 1 / 5), (-1, 1), (0.01, 0.01)],
        options={"maxiter": 10000},
        method="Powell",
    )
    if result.success:
        print((result.x), result.fun, (1 - result.fun / (2 * X.size)) * 100, counter)
    else:
        print("NO")

    sigma, shift, thr = result.x

    plt.plot(X, Y_exp_interp)
    plt.fill_between(X, Y_exp_interp - confidence, Y_exp_interp + confidence, alpha=0.2)
    for i in graphs:
        plt.plot(X, graphs[i], "--", alpha=0.4)

    plt.plot(X, Graph.normalise(computed.calc_graph(shift, sigma)))
    plt.show()
