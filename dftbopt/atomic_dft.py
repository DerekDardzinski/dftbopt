from hotcent.atomic_dft import AtomicDFT
from ase.units import Bohr
from ase.build import bulk
from ase.data import atomic_numbers, covalent_radii
from itertools import combinations_with_replacement, repeat
from dd_hotcent.hotcent.atomic_dft import AtomicDFT
from dd_hotcent.hotcent.confinement import PowerConfinement as PC
from dd_hotcent.hotcent.slako import SlaterKosterTable
from pymatgen.core.periodic_table import Element
# from multiprocessing import Pool, cpu_count
import os

def get_hubbard_values(elements, xc='GGA_X_PBE+GGA_C_PBE'):
    valences = {
        e: Element(e).electronic_structure.split('.')[1:] for e in elements
    }
    configuration = {
        e: Element(e).electronic_structure.replace('.', ' ') for e in elements
    }
    U_values = {}
    for e in elements:
        atom = AtomicDFT(
            e,
            xc=xc,
            configuration=configuration[e],
            valence=[v[:2] for v in valences[e]],
            scalarrel=False,
            mix=0.05,
            maxiter=500,
            confinement=PC(r0=40., s=4),
            txt=None,
            verbose=False,
        )
        atom.run()

        U_tmp = {}
        for v in valences[e]:
            U = atom.get_hubbard_value(
                v[:2],
                scheme='central',
                maxstep=1,
            )
            U_tmp[v[1]] = U

        U_values[e] = U_tmp

    return U_values


def write_slako(elements, params, hubbardvalues, slako_dir='slako', dir_index='1', xc='GGA_X_PBE+GGA_C_PBE'):
    valences = {
        e: Element(e).electronic_structure.split('.')[1:] for e in elements
    }
    configuration = {
        e: Element(e).electronic_structure.replace('.', ' ') for e in elements
    }

    atoms = {}
    eigenvalues = {}
    occupations = {}
    for e in elements:
        conf = PC(r0=params[f'{e}_n_r0'], s=2)
        wf_conf = {v: PC(r0=params[f'{e}_{v[:2]}_r0'], s=2) for v in valences[e]}
        atom = AtomicDFT(
            e,
            xc=xc,
            configuration=configuration[e],
            valence=[v[:2] for v in valences[e]],
            scalarrel=False,
            confinement=conf,
            wf_confinement=wf_conf,
            verbose=False,
            mix=0.05,
            maxiter=500,
            txt=None,
        )
        atom.run()
        atoms[e] = atom
        eigenvalues[e] = {nl: atom.get_eigenvalue(nl) for nl in atom.valence}
        occupations[e] = {v[:2]: int(v[2:]) for v in valences[e]}

    atom_combos = combinations_with_replacement(elements, 2)


    if not os.path.isdir(slako_dir):
        os.mkdir(slako_dir)

    if not os.path.isdir(os.path.join(slako_dir, dir_index)):
        os.mkdir(os.path.join(slako_dir, dir_index))

    rmin, dr, N = 0.4, 0.02, 600

    for atom_combo in atom_combos:
        sk = SlaterKosterTable(atoms[atom_combo[0]], atoms[atom_combo[1]], txt=None)
        sk.run(rmin, dr, N, superposition='density', xc=xc, wflimit=1e-4)

        if atom_combo[0] == atom_combo[1]:
            sk.write(
                os.path.join(slako_dir, dir_index, f'{atom_combo[0]}-{atom_combo[1]}.skf'),
                spe=0.,
                eigenvalues=eigenvalues[atom_combo[0]],
                hubbardvalues=hubbardvalues[atom_combo[0]],
                occupations=occupations[atom_combo[0]],
            )
        else:
            sk.write(
                os.path.join(slako_dir, dir_index, f'{atom_combo[0]}-{atom_combo[1]}.skf'),
                pair=(atom_combo[0], atom_combo[1]),
                spe=0.
            )
            sk.write(
                os.path.join(slako_dir, dir_index, f'{atom_combo[1]}-{atom_combo[0]}.skf'),
                pair=(atom_combo[1], atom_combo[0]),
                spe=0.
            )


# class DFT:
#     def __init__(self, elements, xc='GGA_X_PBE+GGA_C_PBE'):
#         self.elements = elements
#         self.valences = {
#             e: Element(e).electronic_structure.split('.')[1:] for e in elements
#         }
#         self.configuration = {e: Element(e).electronic_structure.replace('.', ' ') for e in elements}
#         # self.valences = {'In': ['5s2', '5p1'], 'As': ['4s2', '4p3']}
#         # self.configuration = {'In': '[Kr] 4d10 5s2 5p1', 'As': '[Ar] 3d10 4s2 4p3'}
#         # print(self.configuration)
#         # print(self.valences)
#         self.xc = xc
#         self.atoms = self.init_atoms()


#     def __getstate__(self):
#         d = self.__dict__.copy()
#         print('DFT')
#         for key in self.__dict__:
#             print(key, type(d[key]))

#     def __setstate__(self, state):
#         self.__dict__.update(state)

#     def get_hubbard_values(self):
#         U_values = {}
#         for e in self.elements:
#             atom = AtomicDFT(
#                 e,
#                 xc=self.xc,
#                 configuration=self.configuration[e],
#                 valence=[v[:2] for v in self.valences[e]],
#                 scalarrel=False,
#                 mix=0.05,
#                 maxiter=500,
#                 confinement=PC(r0=40., s=4),
#                 txt=None,
#                 verbose=False,
#             )
#             atom.run()

#             U_tmp = {}
#             for v in self.valences[e]:
#                 U = self.atoms[e].get_hubbard_value(
#                     v[:2],
#                     scheme='central',
#                     maxstep=1,
#                 )
#                 U_tmp[v[1]] = U

#             U_values[e] = U_tmp

#         return U_values


#     def init_atoms(self):
#         atoms = {}
#         for e in self.elements:
#             atom = AtomicDFT(
#                 e,
#                 xc=self.xc,
#                 configuration=self.configuration[e],
#                 valence=[v[:2] for v in self.valences[e]],
#                 scalarrel=False,
#                 mix=0.05,
#                 maxiter=500,
#                 confinement=PC(r0=40., s=4),
#                 txt=None,
#                 verbose=False,
#             )
#             atoms[e] = atom

#         return atoms

#     def run_init_atoms(self):
#         for e in self.elements:
#             self.atoms[e].run()
#             self.atoms[e].info = {}
#             self.atoms[e].info['eigenvalues'] = {
#                 nl: self.atoms[e].get_eigenvalue(nl) for nl in self.atoms[e].valence
#             }

#             U_tmp = {}
#             for v in self.valences[e]:
#                 U = self.atoms[e].get_hubbard_value(
#                     v[:2],
#                     scheme='central',
#                     maxstep=1,
#                 )
#                 U_tmp[v[1]] = U

#             self.atoms[e].info['hubbardvalues'] = U_tmp
#             self.atoms[e].info['occupations'] = {
#                 v[:2]: int(v[2:]) for v in self.valences[e]
#             }


#     def generate_slako(self, params, U_values, slako_dir='./slako'):
#         for e in self.elements:
#             atom = AtomicDFT(
#                 element,
#                 xc=xc,
#                 configuration=Element(element).electronic_structure.replace('.', ' '),
#                 valence=[v[:2] for v in valences[element]],
#                 scalarrel=False,
#                 confinement=conf,
#                 wf_confinement=wf_conf,
#             )
#             atom.run()
#             atoms[element] = atom

#             conf = PC(r0=params[f'{e}_n_r0'], s=2)
#             wf_conf = {v: PC(r0=params[f'{e}_{v[:2]}_r0'], s=2) for v in self.valences[e]}
#             self.atoms[e].set_confinement(conf)
#             self.atoms[e].set_wf_confinement(wf_conf)
#             self.atoms[e].run()
#             self.atoms[e].info['eigenvalues'] = {
#                 nl: self.atoms[e].get_eigenvalue(nl) for nl in self.atoms[e].valence
#             }

#         atom_combos = combinations_with_replacement(self.elements, 2)

#         if not os.path.isdir(slako_dir):
#             os.mkdir(slako_dir)

#         rmin, dr, N = 0.4, 0.02, 600

#         for atom_combo in atom_combos:
#             sk = SlaterKosterTable(self.atoms[atom_combo[0]], self.atoms[atom_combo[1]])
#             sk.run(rmin, dr, N, superposition='density', xc=self.xc, wflimit=1e-4)

#             if atom_combo[0] == atom_combo[1]:
#                 sk.write(
#                     os.path.join(slako_dir, f'{atom_combo[0]}-{atom_combo[1]}.skf'),
#                     spe=0.,
#                     **self.atoms[atom_combo[0]].info
#                 )
#             else:
#                 sk.write(
#                     os.path.join(slako_dir, f'{atom_combo[0]}-{atom_combo[1]}.skf'),
#                     pair=(atom_combo[0], atom_combo[1]),
#                     spe=0.
#                 )
#                 sk.write(
#                     os.path.join(slako_dir, f'{atom_combo[1]}-{atom_combo[0]}.skf'),
#                     pair=(atom_combo[1], atom_combo[0]),
#                     spe=0.
#                 )

#     def generate_slako_old(self, params, slako_dir='./slako'):
#         print(self.__dict__)
#         for e in self.elements:
#             conf = PC(r0=params[f'{e}_n_r0'], s=2)
#             wf_conf = {v: PC(r0=params[f'{e}_{v[:2]}_r0'], s=2) for v in self.valences[e]}
#             self.atoms[e].set_confinement(conf)
#             self.atoms[e].set_wf_confinement(wf_conf)
#             self.atoms[e].run()
#             self.atoms[e].info['eigenvalues'] = {
#                 nl: self.atoms[e].get_eigenvalue(nl) for nl in self.atoms[e].valence
#             }

#         atom_combos = combinations_with_replacement(self.elements, 2)

#         if not os.path.isdir(slako_dir):
#             os.mkdir(slako_dir)

#         rmin, dr, N = 0.4, 0.02, 600

#         for atom_combo in atom_combos:
#             sk = SlaterKosterTable(self.atoms[atom_combo[0]], self.atoms[atom_combo[1]])
#             sk.run(rmin, dr, N, superposition='density', xc=self.xc, wflimit=1e-4)

#             if atom_combo[0] == atom_combo[1]:
#                 sk.write(
#                     os.path.join(slako_dir, f'{atom_combo[0]}-{atom_combo[1]}.skf'),
#                     spe=0.,
#                     **self.atoms[atom_combo[0]].info
#                 )
#             else:
#                 sk.write(
#                     os.path.join(slako_dir, f'{atom_combo[0]}-{atom_combo[1]}.skf'),
#                     pair=(atom_combo[0], atom_combo[1]),
#                     spe=0.
#                 )
#                 sk.write(
#                     os.path.join(slako_dir, f'{atom_combo[1]}-{atom_combo[0]}.skf'),
#                     pair=(atom_combo[1], atom_combo[0]),
#                     spe=0.
#                 )


if __name__ == '__main__':
    elements = ['In', 'As']
    hubbardvalues = get_hubbard_values(elements=elements)
    # print('Init')
    # dft = DFT(elements=elements)
    # print('Starting Hubbard Value Calc')
    # dft.run_init_atoms()
    # print('Finished Hubbard Value Calc')
    valences = {'In': ['4d10', '5s2', '5p1'], 'As': ['3d10', '4s2', '4p3']}
    rcovs = {e: covalent_radii[atomic_numbers[e]] / Bohr for e in elements}
    params = {}

    for e in elements:
        params[f'{e}_n_r0'] = 3 * rcovs[e]
        for v in valences[e]:
            params[f'{e}_{v[:2]}_r0'] = 2 * rcovs[e]

    print('Generate Slako Files')
    write_slako(elements=elements, params=params, hubbardvalues=hubbardvalues)
    print('Done!')
