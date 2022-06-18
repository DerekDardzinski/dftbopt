from multiprocessing import Pool, cpu_count
from atomic_dft import get_hubbard_values, write_slako, get_valences, get_configurations
from pymatgen.core.periodic_table import Element
from dftb_band import DFTBBand
from ase.data import atomic_numbers, covalent_radii
from ase.units import Bohr
import pyswarms as ps
import numpy as np
import os

class Optimizer:
    def __init__(self, atoms):
        self.atoms = atoms
        self.elements = np.unique(atoms.get_chemical_symbols())
        self.hubbardvalues = get_hubbard_values(elements=self.elements)
        self.valences = get_valences(self.elements)
        self.init_guess = self._get_init_guess()


    def _get_init_guess(self):
        params = {}
        for e in self.elements:
            params[f'{e}_r0'] = 2 * (covalent_radii[atomic_numbers[e]] / Bohr)
            params[f'{e}_sigma'] = 2

        return params
    
    def _get_bounds(self):
        lower_bound = []
        upper_bound = []
        for key in self.init_guess:
            if 'r0' in key:
                upper_bound.append(self.init_guess[key] + 5)
                if self.init_guess[key] - 2 < 3:
                    lower_bound.append(3)
                else:
                    lower_bound.append(self.init_guess[key] - 2)

            if 'sigma' in key:
                lower_bound.append(2)
                upper_bound.append(6)

        return np.array(lower_bound), np.array(upper_bound)

    def _f(self, params, index):
        params_dict = dict(zip(self.init_guess.keys(), params))
        print(f'Params {index} =', params_dict)
        try:
            write_slako(
                elements=self.elements,
                params=params_dict,
                hubbardvalues=self.hubbardvalues,
                slako_dir='slako',
                dir_index=str(index),
            )
        except:
            os.rmdir(os.path.join('slako', str(index)))


    def f(self, params):
        inputs = zip(params, range(len(params)))
        with Pool(cpu_count()) as p:
            p.starmap(self._f, inputs)

        l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        rl_map = dict(zip(l_map.values(), l_map.keys()))

        maximum_angular_momenta = {
            e: rl_map[max([l_map[v[1]] for v in self.valences[e]])] for e in self.elements
        }

        delta_bands = []
        for i in range(len(params)):
            try:
                dftb_band = DFTBBand(
                    atoms=self.atoms,
                    maximum_angular_momenta=maximum_angular_momenta,
                    dir_index=str(i),
                )

                dftb_band.run_chg(kpoint_density=(11,11,11))
                dftb_band.run_band(coords='GXWLGK')
                dftb_band.add_dft_bands(vasp_folder='./dft/band')
                delta_band, _ = dftb_band.get_delta_band()
                delta_bands.append(delta_band)
            except:
                print('DD ERROR')
                delta_bands.append(np.nan)

        print(delta_bands)
        return delta_bands

    def run(self, n_particles=6, c1=0.5, c2=0.3, w=0.9, n_processes=cpu_count(), iters=4):
        init_pos = np.array(list(self.init_guess.values()))
        lower_bound, upper_bound = self._get_bounds()
        bounds = (lower_bound, upper_bound)
        print('Init Guess = ', init_pos)
        print('Lower Bound =', lower_bound)
        print('Upper Bound =', upper_bound)
        opt = ps.single.GlobalBestPSO(
            n_particles=n_particles,
            dimensions=len(init_pos),
            options={'c1': c1, 'c2': c2, 'w': w},
            # init_pos=init_pos,
            bounds=bounds,
        )

        cost, pos = opt.optimize(
            self.f,
            n_processes=None,
            iters=iters,
        )

        print(opt.cost_history)


if __name__ == '__main__':
    # set_start_method('fork')
    from ase.io import read
    atoms = read('./POSCAR')
    opt = Optimizer(atoms)
    opt.run()
