import os, glob, sys

# Set up matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.table import Table
import numpy as np
import os
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
os.chdir("/home/clem/WSL_FERMI")


def last_divider(x, d):
    force_stop = 0
    n = d
    find = False
    while find != True and force_stop < 10000000:
        force_stop += 1
        print(force_stop)
        if x % n == 0:
            find = True
            return n
        else:
            n = n - 1


def make_divide_list(x, d):
    if x % d != 0:
        rest = x % d
        res = [d for i in range(x // d)]
        res += [x % d]
    else:
        res = [x // d for i in range(d)]
    return res


def divisability(i, output=False):
    """
    test la divisabilité de manière conne d'un nombre et retourne les diviseurs
    :param i:
    :return:
    """
    notprimer = False
    diviseur = []
    division_result = []
    if output == True:
        for n in tqdm(range(2, i // 2)):
            if i % n == 0:
                notprimer = True
                print(f"{i:_} est divisible par {n:_}")
                diviseur += [n]
                division_result += [i // n]
            else:
                pass
        if notprimer == False:
            print(f"{i:_} est peut être premier...")

        return np.array(diviseur), np.array(division_result)

    else:
        for n in tqdm(range(2, i // 2)):
            if i % n == 0:
                notprimer = True
                diviseur += [n]
                division_result += [i // n]
            else:
                pass
        return np.array(diviseur), np.array(division_result)


class GammaFits:
    def __init__(self, fitsfile):
        """Initialise un objet de la classe GammaFits à partir d'un fichier FITS.

        Cette méthode lit les données des photons gamma détectés par le télescope spatial Fermi
        et les stocke dans des attributs sous forme de tableaux numpy.

        Args:
            fitsfile (str): le nom du fichier FITS contenant les données des photons gamma.

        Attributes:
            fitsfile (str): le chemin du fichier FITS.
            f1 (HDUList): un objet HDUList ouvert avec la méthode fits.open().
            events (Table): un objet Table contenant les données des événements du télescope Fermi.
            energy (ndarray): un tableau numpy contenant les énergies des photons gamma en MeV.
            L (ndarray): un tableau numpy contenant les longitudes galactiques des photons gamma en degrés.
            B (ndarray): un tableau numpy contenant les latitudes galactiques des photons gamma en degrés.
            zenith_angle (ndarray): un tableau numpy contenant les angles zénithaux des photons gamma en degrés.
        """

        self.fitsfile = fitsfile
        self.f1 = fits.open(self.fitsfile)
        self.events = Table.read(self.fitsfile, hdu=1)
        self.energy = self.events['ENERGY']
        self.L = self.events['L']
        self.B = self.events['B']
        self.zenith_angle = self.events['ZENITH_ANGLE']
        self.TIME = self.events["TIME"]
        self.events_number = len(self.TIME)
        self.slicing_posibility, self.events_number_per_time_slice = divisability(self.events_number)
        self.events_copy = Table.copy(self.events)
        self.keep_bin_edge()

    def reset(self):
        """
        Reset l'event actuelle à la semaine actuelle
        :return:
        """
        self.fitsfile = self.fitsfile
        self.f1 = fits.open(self.fitsfile)
        self.events = Table.read(self.fitsfile, hdu=1)
        self.energy = self.events['ENERGY']
        self.L = self.events['L']
        self.B = self.events['B']
        self.zenith_angle = self.events['ZENITH_ANGLE']
        self.TIME = self.events["TIME"]
        self.events_number = len(self.TIME)
        self.slicing_posibility, self.events_number_per_time_slice = divisability(self.events_number)

    def set_curent_slice_as_event(self):
        """
        Sets the currently selected time slice as the main event for the object, so that it can be used for plotting and analysis.

        This function assumes that the currently selected time slice is stored in the `current_slice` attribute of the object.

        Upon execution, the function first sets the `events` attribute of the object to the currently selected slice.

        It then sets the following attributes of the object based on the data in the selected slice:
        - `energy`: the energy of each event in the slice
        - `L`: the longitude of each event in the slice
        - `B`: the latitude of each event in the slice
        - `zenith_angle`: the zenith angle of each event in the slice
        - `TIME`: the time of each event in the slice
        - `events_number`: the total number of events in the slice

        These attributes are used for subsequent plotting and analysis operations on the selected time slice.
        """
        self.events = self.curent_slice
        self.energy = self.events['ENERGY']
        self.L = self.events['L']
        self.B = self.events['B']
        self.zenith_angle = self.events['ZENITH_ANGLE']
        self.TIME = self.events["TIME"]
        self.events_number = len(self.TIME)

    def divby(self):
        divisability(self.events_number, output=True)

    def show_info(self):
        """affiche les information du fits
        """
        fits.info(self.fitsfile)
        print(self.events.columns)
        print(self.events['ENERGY'].unit)
        print(self.events['ENERGY'][0])

    def zenith_angle_corection(self):
        """Corrige les mesures des photons gamma en fonction de l'angle zénithal.

        Les photons gamma détectés par le télescope spatial Fermi sont filtrés
        pour ne garder que ceux dont l'angle zénithal est inférieur à 80 degrés.
        L'angle zénithal est l'angle entre le photon et la verticale.
        Les valeurs corrigées sont stockées dans les attributs L_cut et B_cut.

        Args:
            self: objet de la classe contenant les attributs L, B et zenith_angle.

        Returns:
            None
    """
        zenith_angle_mask = self.zenith_angle < 80
        self.L_cut = self.L[zenith_angle_mask]
        self.B_cut = self.B[zenith_angle_mask]

    def time_slice(self,
                   nb_total_slice,
                   slice_selcte=0,
                   force_slice_number=False
                   ):
        """
        fonction de découpe d'une table en slice.
        :param nb_total_slice: nombre de slice désiré
        :param slice_selcte: slice actuelle
        :return:
        """
        if force_slice_number == False:
            if nb_total_slice not in self.slicing_posibility:
                nb_total_slice = np.sort(np.abs(self.slicing_posibility - nb_total_slice))[0] + nb_total_slice
                print(f"Le nombre de slice total désiré ne divise pas le nombre de données.\n"
                      f" Le nombre le plus proche trouvé est : {nb_total_slice}")

            self._number_per_slice = self.events_number // nb_total_slice

            print(f"::\t le nombre d'event par slice est de {self._number_per_slice}")
            self.curent_slice = self.events[
                                self._number_per_slice * slice_selcte:self._number_per_slice * (slice_selcte + 1)]
            self.divide_table = []
            self._nb_total_slice = nb_total_slice
        if force_slice_number:
            self._number_per_slice = self.events_number // nb_total_slice

    def slice_gen(self, i, force_slice_number=False):
        """génère la slice i"""

        if force_slice_number:
            self #TODO Ajouter le cas ou on ne prend pas une division entière
            self.curent_slice = self.events_copy[self._number_per_slice * i:self._number_per_slice * (i + 1)]
        else:
            self.curent_slice = self.events_copy[self._number_per_slice * i:self._number_per_slice * (i + 1)]

    def make_all_slice(self, nb_slice):
        """génère toute les slice et le stock dans une lsite"""
        self.time_slice(nb_slice)
        self.all_slice = []
        self.name_slice = []
        print("Création des slice en cours")
        for i in tqdm(range(self._nb_total_slice)):
            self.name_slice += [str(i)]
            self.slice_gen(i)
            self.all_slice += [self.curent_slice]

    def set_slice_as_event(self, slice):
        """Définie la slice renseigné en paramètre comme jeux de donné principal"""
        self.curent_slice = slice
        self.set_curent_slice_as_event()

    def plot_all_slice_hist(self, norm='linear', vmin=0, vmax=200, show=True):
        """plot all hist"""

        i = 1
        print("Calcule de la statistique de toute les slice en cours ...")
        self.set_slice_as_event(self.all_slice[0])
        self.plot_hist2d(norm='linear', bins=(self.xedges_tot, self.yedges_tot), vmin=vmin, vmax=vmax,
                         title=self.name_slice[0], show=show)
        self.all_counts = self.counts
        self.all_xedges = self.xedges
        self.all_yedges = self.yedges

        for slice in tqdm(self.all_slice[1:]):
            self.set_slice_as_event(slice)
            self.plot_hist2d(norm='linear', bins=(self.xedges_tot, self.yedges_tot), vmin=vmin, vmax=vmax,
                             title=self.name_slice[i],show=show)
            i += 1
            self.all_counts = np.dstack((self.all_counts, self.counts))
            self.all_xedges = np.dstack((self.all_xedges, self.xedges))
            self.all_yedges = np.dstack((self.all_yedges, self.yedges))
            # self.all_im = np.dstack((self.all_im, self.im))

            # if show:
            #     plt.show()
            # else:
            #     plt.clf()
            #     plt.close()

    def make_pcolormesh_param(self, slice_number):
        """
        Donne les paramètre du pcolormesh corespondant pour replot le hist2d en pcolor mesh
        :param slice_number: la slice selectionné
        :return: X,Y,Z : X contient les bort en X, Y les bord en Y, et Z le nombre de couts compté
        """
        X, Y = np.meshgrid(self.xedges_tot, self.yedges_tot)
        Z = w15.all_counts[:, :, slice_number].T
        return X, Y, Z

    def plot_hist2d(self, norm=LogNorm(), bins=[360, 180], vmin=0, vmax=200, title='',show=True):
        """plot le hist 2D de la carte du ciel, avec et sans corection du zenith angle
        """

        self.zenith_angle_corection()

        if show:
            fig, ax = plt.subplots()
            ax.set_xlabel('L (deg.)', fontsize=10)
            ax.set_ylabel('B (deg.)', fontsize=10)

            self.counts, self.xedges, self.yedges, self.im = ax.hist2d(self.L_cut, self.B_cut, bins=bins,
                                                                       norm='linear', vmin=vmin, vmax=vmax)
            plt.title(title)

            fig.colorbar(self.im, ax=ax)
            plt.show()
        if not show:
            self.counts, self.xedges, self.yedges = np.histogram2d(self.L_cut, self.B_cut, bins=bins)
            # self.im =

    def keep_bin_edge(self):
        """plot le hist 2D de la carte du ciel, avec et sans corection du zenith angle
        """
        fig, ax = plt.subplots()
        self.zenith_angle_corection()
        self.all_counts, self.xedges_tot, self.yedges_tot, im = ax.hist2d(self.L_cut, self.B_cut, bins=[360, 180],
                                                                          norm='linear')
        plt.clf()

    def plot_hist2d_compar(self):
        """plot le hist 2D de la carte du ciel, avec et sans corection du zenith angle

        """
        fig, axs = plt.subplots(2, 1, figsize=(8, 8))

        emin_mask = self.energy > 50000.00
        axs[0].set_xlabel('L (deg.)', fontsize=10)
        axs[0].set_ylabel('B (deg.)', fontsize=10)
        counts0, xedges0, yedges0, im0 = axs[0].hist2d(self.L, self.B, bins=[360, 180], norm=LogNorm())
        fig.colorbar(im0, ax=axs[0])
        self.zenith_angle_corection()
        axs[1].set_xlabel('L (deg.)', fontsize=10)
        axs[1].set_ylabel('B (deg.)', fontsize=10)
        counts1, xedges1, yedges1, im1 = axs[1].hist2d(self.L_cut, self.B_cut, bins=[360, 180], norm=LogNorm())
        fig.colorbar(im1, ax=axs[1])
        plt.show()


class TwoWeek_Gamma_Fits:
    def __init__(self, week1, week2):
        """
        permet la comparaison de deux semaine de donné
        Args:
            week1 (_type_): _description_
            week2 (_type_): _description_
        """
        self.W1 = week1
        self.W2 = week2

    def hist2D_compar(self):
        """
        Plot les deux histogram superpose
        :return:
        """
        fig, axs = plt.subplots(2, 1, figsize=(8, 8))
        axs[0].set_xlabel('L (deg.)', fontsize=10)
        axs[0].set_ylabel('B (deg.)', fontsize=10)
        self.W1.zenith_angle_corection()
        self.counts1, self.xedges1, self.yedges1, self.im1 = axs[0].hist2d(self.W1.L_cut, self.W1.B_cut,
                                                                           bins=[360, 180], norm=LogNorm())
        fig.colorbar(self.im1, ax=axs[0])
        self.W2.zenith_angle_corection()
        axs[1].set_xlabel('L (deg.)', fontsize=10)
        axs[1].set_ylabel('B (deg.)', fontsize=10)
        self.counts2, self.xedges2, self.yedges2, self.im2 = axs[1].hist2d(self.W2.L_cut, self.W2.B_cut,
                                                                           bins=[360, 180], norm=LogNorm())
        axs[0].set_title("hist_compar\n Week A")
        axs[1].set_title("week B")
        fig.colorbar(self.im2, ax=axs[1])
        plt.tight_layout()
        plt.show()


# gf = GammaFits(fitsfile)
# gf.show_info()
# gf.plot_hist2d()

if __name__ == "__main__":
    # w15 = GammaFits("HauteEnergie/Data/source_observation/lat_photon_weekly_w015_p305_v001.fits")
    # w15.make_all_slice(23)

    # w15.plot_all_slice_hist(vmin=0, vmax=150, show=False)
    w15 = GammaFits("HauteEnergie/Data/source_observation/lat_photon_weekly_w015_p305_v001.fits")
    w15.make_all_slice(10000)
    w15.plot_all_slice_hist(vmin=0, vmax=5,show=False)
    print("DONE")
#     w15 = GammaFits("HauteEnergie/Data/source_observation/lat_photon_weekly_w015_p305_v001.fits")
#     w16 = GammaFits("HauteEnergie/Data/source_observation/lat_photon_weekly_w016_p305_v001.fits")
#     e15 = w15.events
#     e16 = w16.events
#     w15.plot_hist2d()
#     w15_16 = TwoWeek_Gamma_Fits(w15, w16)
#     w15_16.hist2D_compar()
#     plt.imshow(((np.transpose(np.abs(np.abs(w15_16.counts2) - np.abs(w15_16.counts1)))) * 100) > 10000)
#
#
