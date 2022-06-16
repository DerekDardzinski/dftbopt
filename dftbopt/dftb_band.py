from ase.calculators.dftb import Dftb
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Kpoints
from ase.io import write, read
from ase.build import bulk
from dptools.bandout import BandOut
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
        slako_dir,
    ) -> None:
        self.atoms = atoms
        self.maximum_angular_momenta = {'': ''} | maximum_angular_momenta
        self.slako_dir = os.path.join(os.getcwd(), slako_dir)
        if self.slako_dir[-1] != '/':
            self.slako_dir += '/'

        self.efermi = None

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
        folder_path = os.path.join(current_dir, 'chg')
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)

        os.chdir(folder_path)

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
        folder_path = os.path.join(current_dir, 'band')
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)

        shutil.copyfile(
            os.path.join(current_dir, 'chg', 'charges.bin'),
            os.path.join(folder_path, 'charges.bin'),
        )
        print(os.listdir(folder_path))

        os.chdir(folder_path)

        calc = Dftb(
            atoms=atoms,
            slako_dir=self.slako_dir,
            kpts=all_kpoints,
            **kwargs,
        )

        atoms.set_calculator(calc)
        calc.calculate(atoms)
        self.efermi = calc.get_fermi_level()

        os.chdir(current_dir)

    def plot_band(self):
        band_dir = os.path.join(os.getcwd(), 'band', 'band.out')
        bandout = BandOut.fromfile(band_dir)
        eigenvalues = bandout.eigvalarray.transpose(0,2,1,3)[:,:,:,0] - self.efermi
        nspin, nbands, nkpts = eigenvalues.shape

        # pin = hsd.load('./dftb_pin.hsd')
        # kpoints = np.array(pin['Hamiltonian']['DFTB']['KPointsAndWeights'])[:,:3]

        fig, ax = plt.subplots(figsize=(4,3), dpi=400)
        for i in range(nspin):
            for j in range(nbands):
                ax.plot(
                    np.arange(nkpts),
                    eigenvalues[i,j],
                )

        ax.set_ylim(-10,10)
        ax.set_xlim(0,nkpts)
        fig.tight_layout(pad=0.4)
        fig.savefig('bs.png')

        




if __name__ == '__main__':
    atoms = read('./POSCAR')
    dftb_band = DFTBBand(
        atoms=atoms,
        maximum_angular_momenta={'In': 'd', 'As': 'd'},
        slako_dir='slako',
    )
    dftb_band.run_chg()
    dftb_band.run_band(coords='GXWLGK')
    dftb_band.plot_band()
