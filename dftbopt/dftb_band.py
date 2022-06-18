from ase.calculators.dftb import Dftb
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Kpoints
from ase.io import write, read
from ase.build import bulk
from dptools.bandout import BandOut
from vaspvis import Band
import matplotlib.pyplot as plt
import numpy as np
import hsd
import shutil
import copy
import numpy as np
import os


class DFTBBand:
    def __init__(
        self,
        atoms,
        maximum_angular_momenta,
        slako_dir='slako',
        dftb_dir='dftb',
        dir_index='1',
    ) -> None:
        self.atoms = atoms
        self.maximum_angular_momenta = {'': ''} | maximum_angular_momenta
        self.slako_dir = os.path.join(os.getcwd(), slako_dir, dir_index)
        if self.slako_dir[-1] != '/':
            self.slako_dir += '/'

        self.dftb_dir = dftb_dir
        self.dir_index = dir_index

        self.efermi = None
        self.dft_eigenvalues = None
        self.dftb_eigenvalues = None

    def __getstate__(self):
        print('DFTBand')
        d = self.__dict__.copy()
        for key in self.__dict__:
            print(key, type(d[key]))

    def __setstate__(self, state):
        self.__dict__.update(state)

    def _interpolate(self, point1, point2, n):
        xs = np.linspace(point1[0], point2[0], n)
        ys = np.linspace(point1[1], point2[1], n)
        zs = np.linspace(point1[2], point2[2], n)

        return np.c_[xs, ys, zs]

    def get_high_symm_points(self, coords, n):
        st = AseAtomsAdaptor().get_structure(self.atoms)
        st = SpacegroupAnalyzer(st).get_primitive_standard_structure()
        kp = HighSymmKpath(st)
        kpoints = Kpoints.automatic_linemode(n, kp)
        kpoints.labels = ';'.join(kpoints.labels).replace('\\Gamma', 'G').split(';')

        coords = coords.upper()
        if all(c in kpoints.labels for c in coords):
            kpath = [c for c in coords]

            for i in reversed(range(len(kpath))):
                if i < len(kpath) - 1 and i > 0:
                    kpath.insert(i+1, kpath[i])

            info_dict = {l: c for l, c in zip(kpoints.labels, kpoints.kpts)}

            kpoints.labels = kpath
            kpoints.kpts = [info_dict[l] for l in kpath]

        return kpoints

    def run_chg(self, kpoint_density=(7,7,7)):
        atoms = copy.deepcopy(self.atoms)

        kwargs = {}
        for k, v in self.maximum_angular_momenta.items():
            kwargs[f'Hamiltonian_MaxAngularMomentum_{k}'] = v

        kwargs['Hamiltonian_Scc'] = 'Yes'
        kwargs['Hamiltonian_SccTolerence'] = 1e-5
        kwargs['Hamiltonian_PolynomialRepulsive'] = 'SetForAll { Yes }'

        current_dir = os.getcwd()
        dftb_path = os.path.join(current_dir, self.dftb_dir)
        index_path = os.path.join(dftb_path, self.dir_index)
        chg_path = os.path.join(index_path, 'chg')

        if not os.path.isdir(dftb_path):
            os.mkdir(dftb_path)

        if not os.path.isdir(index_path):
            os.mkdir(index_path)

        if not os.path.isdir(chg_path):
            os.mkdir(chg_path)

        os.chdir(chg_path)

        calc = Dftb(
            atoms=atoms,
            slako_dir=self.slako_dir,
            kpts=kpoint_density,
            **kwargs,
        )
        atoms.set_calculator(calc)
        calc.calculate(atoms)

        os.chdir(current_dir)


    def run_band(self, coords, n=50):
        atoms = copy.deepcopy(self.atoms)

        kwargs = {}

        for k, v in self.maximum_angular_momenta.items():
            kwargs[f'Hamiltonian_MaxAngularMomentum_{k}'] = v

        kwargs['Hamiltonian_Scc'] = 'Yes'
        kwargs['Hamiltonian_SccTolerence'] = 1e-5
        kwargs['Hamiltonian_ReadInitialCharges'] = 'Yes'
        kwargs['Hamiltonian_PolynomialRepulsive'] = 'SetForAll { Yes }'


        kpoints = self.get_high_symm_points(coords=coords, n=n)
        num_kpts = kpoints.num_kpts
        kpts = kpoints.kpts

        all_kpoints = np.vstack([
            self._interpolate(kpts[i], kpts[i+1], num_kpts) for i in range(len(kpts) - 1) if not i % 2
        ])

        current_dir = os.getcwd()
        dftb_path = os.path.join(current_dir, self.dftb_dir)
        index_path = os.path.join(dftb_path, self.dir_index)
        band_path = os.path.join(index_path, 'band')

        if not os.path.isdir(dftb_path):
            os.mkdir(dftb_path)

        if not os.path.isdir(index_path):
            os.mkdir(index_path)

        if not os.path.isdir(band_path):
            os.mkdir(band_path)

        os.chdir(band_path)

        shutil.copyfile(
            os.path.join(index_path, 'chg', 'charges.bin'),
            os.path.join(band_path, 'charges.bin'),
        )

        print(os.listdir(band_path))

        os.chdir(band_path)

        calc = Dftb(
            atoms=atoms,
            slako_dir=self.slako_dir,
            kpts=all_kpoints,
            **kwargs,
        )

        atoms.set_calculator(calc)
        calc.calculate(atoms)

        os.chdir(current_dir)

        self.efermi = calc.get_fermi_level()

        band_file = os.path.join(band_path, 'band.out')
        bandout = BandOut.fromfile(band_file)
        dftb_eigenvalues = bandout.eigvalarray.transpose(0,2,1,3)[0,:,:,0] - self.efermi
        self.dftb_eigenvalues = dftb_eigenvalues


    def add_dft_bands(self, vasp_folder):
        dft_band = Band(folder=vasp_folder)
        dft_eigenvalues = dft_band.eigenvalues
        self.dft_eigenvalues = dft_eigenvalues

    def _locate_and_shift_bands(self, eigenvalues, band_range):
        band_mean = eigenvalues.mean(axis=1)

        below_index = np.where(band_mean < 0)[0]
        above_index = np.where(band_mean >= 0)[0]

        vbm = np.max(eigenvalues[below_index])
        cbm = np.min(eigenvalues[above_index])
        bg = cbm - vbm

        if cbm < vbm:
            vbm = 0.0
            cbm = 0.0
            bg = 0

        valence_bands = eigenvalues[below_index[-band_range[0]:]]
        conduction_bands = eigenvalues[above_index[:band_range[1]]]

        valence_bands -= vbm
        conduction_bands -= cbm

        shifted_bands = np.r_[conduction_bands, valence_bands]

        return shifted_bands, bg


    def get_delta_band(self, band_range=[5,5]):
        dft_eigenvalues = self.dft_eigenvalues
        dftb_eigenvalues = self.dftb_eigenvalues

        shifted_dft_eigenvalues, dft_bg = self._locate_and_shift_bands(dft_eigenvalues, band_range)
        shifted_dftb_eigenvalues, dftb_bg = self._locate_and_shift_bands(dftb_eigenvalues, band_range)

        n = shifted_dft_eigenvalues.shape[0] * shifted_dft_eigenvalues.shape[1]
        delta_band = np.sum((1/n)*np.sum((shifted_dft_eigenvalues - shifted_dftb_eigenvalues)**2))**(1/2)
        delta_bg = np.abs(dft_bg - dftb_bg)

        return delta_band, delta_bg


    def plot_band(self):
        dftb_eigenvalues = self.dftb_eigenvalues
        dft_eigenvalues = self.dft_eigenvalues
        # nbands, nkpts = dftb_eigenvalues.shape

        # pin = hsd.load('./dftb_pin.hsd')
        # kpoints = np.array(pin['Hamiltonian']['DFTB']['KPointsAndWeights'])[:,:3]

        fig, ax = plt.subplots(figsize=(4,3), dpi=400)

        for j in range(dft_eigenvalues.shape[0]):
            ax.plot(
                np.arange(len(dft_eigenvalues[j])),
                dft_eigenvalues[j],
                color='black',
            )

        for j in range(dftb_eigenvalues.shape[0]):
            ax.plot(
                np.arange(len(dftb_eigenvalues[j])),
                dftb_eigenvalues[j],
                color='red',
            )


        ax.set_ylim(-10,10)
        ax.set_xlim(0,dftb_eigenvalues.shape[1])
        fig.tight_layout(pad=0.4)
        fig.savefig('bs_test.png')

        




if __name__ == '__main__':
    atoms = read('./POSCAR')
    dftb_band = DFTBBand(
        atoms=atoms,
        maximum_angular_momenta={'In': 'd', 'As': 'd'},
        slako_dir='slako_test',
    )
    dftb_band.run_chg()
    dftb_band.run_band(coords='GXWLGK')
    dftb_band.add_dft_bands(vasp_folder='./dft/band')
    dftb_band.plot_band()

    print(dftb_band.get_delta_band())
