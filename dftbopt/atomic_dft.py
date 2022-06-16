from hotcent.atomic_dft import AtomicDFT
from ase.units import Bohr
from ase.build import bulk
from ase.data import atomic_numbers, covalent_radii
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from pymatgen.core.periodic_table import Element


class DFT:
    def __init__(self, elements, xc='GGA_X_PBE+GGA_C_PBE'):
        self.elements = elements
        self.valences = {
            e: Element(e).electronic_structure.split('.')[1:] for e in elements
        }
        self.xc = xc
        self.atoms = self.init_atoms()

    def init_atoms(self):
        atoms = {}
        for element in self.elements:
            e = Element(element)
            atom = AtomicDFT(
                element,
                xc=self.xc,
                configuration=e.electronic_structure.replace('.', ' '),
                valence=[v[:2] for v in self.valences[element]],
                scalarrel=False,
                mix=0.05,
                maxiter=500,
                confinement=PowerConfinement(r0=40., s=4),
                txt='./dft_log.out',
                verbose=True,
            )
            atoms[element] = atom

        return atoms

    def run_init_atoms(self):
        for element in self.elements:
            self.atoms[element].run()
            self.atoms[element].info = {}
            self.atoms[element].info['eigenvalues'] = {
                nl: self.atoms[element].get_eigenvalue(nl) for nl in self.atoms[element].valence
            }

            U_tmp = {}
            for v in self.valences[element]:
                U = self.atoms[element].get_hubbard_value(
                    v[:2],
                    scheme='central',
                    maxstep=1,
                )
                U_tmp[v[1]] = U

            self.atoms[element].info['hubbardvalues'] = U_tmp
            self.atoms[element].info['occupations'] = {
                v[:2]: int(v[2:]) for v in self.valences[element]
            }


if __name__ == '__main__':
    dft = DFT(elements=['In', 'As'])
    dft.run_init_atoms()
    print(dft.atoms['In'].info)
    print(dft.atoms['As'].info)
