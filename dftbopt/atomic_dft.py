from ase.units import Bohr
from ase.build import bulk
from ase.data import atomic_numbers, covalent_radii
from itertools import combinations_with_replacement, repeat
from hotcent.atomic_dft import AtomicDFT
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement as PC
from hotcent.slako import SlaterKosterTable
from pymatgen.core.periodic_table import Element
# from multiprocessing import Pool, cpu_count
import os

ORB_ORDER = [
    '1s',
    '2s',
    '2p',
    '3s',
    '3p',
    '4s',
    '3d',
    '4p',
    '5s',
    '4d',
    '5p',
    '6s',
    '4f',
    '5d',
    '6p',
    '7s',
    '5f',
    '6d',
    '7p',
    '6f',
    '7d',
    '7f',
]

def get_valences(elements):
    valences = {}
    for e in elements:
        valence = Element(e).electronic_structure.split('.')[1:]
        ind = ORB_ORDER.index(valence[-1][:2]) + 1

        found_next_empty = False
        while not found_next_empty:
            if 's' in ORB_ORDER[ind] or 'f' in ORB_ORDER[ind]:
                ind += 1
            else:
                valence.append(ORB_ORDER[ind] + '0')
                found_next_empty = True

        if len(valence) > 3:
            valence = valence[1:]

        valences[e] = valence

    return valences


def get_configurations(elements):
    configs = {}
    for e in elements:
        valence = Element(e).electronic_structure.split('.')
        ind = ORB_ORDER.index(valence[-1][:2]) + 1

        found_next_empty = False
        while not found_next_empty:
            if 's' in ORB_ORDER[ind] or 'f' in ORB_ORDER[ind]:
                ind += 1
            else:
                valence.append(ORB_ORDER[ind] + '0')
                found_next_empty = True

        configs[e] = ' '.join(valence)

    return configs

def get_hubbard_values(elements, xc='GGA_X_PBE+GGA_C_PBE'):
    valences = get_valences(elements)
    configurations = get_configurations(elements)
    U_values = {}
    for e in elements:
        atom = AtomicDFT(
            e,
            xc=xc,
            configuration=configurations[e],
            valence=[v[:2] for v in valences[e]],
            scalarrel=False,
            mix=0.05,
            maxiter=1000,
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
    valences = get_valences(elements)
    configurations = get_configurations(elements)

    atoms = {}
    occupations = {}
    eigenvalues = {}
    for e in elements:
        conf = PC(r0=1.5 * params[f'{e}_r0'], s=params[f'{e}_sigma'])
        wf_conf = {v[:2]: PC(r0=params[f'{e}_r0'], s=params[f'{e}_sigma']) for v in valences[e]}
        atom = AtomicDFT(
            e,
            xc=xc,
            configuration=configurations[e],
            valence=[v[:2] for v in valences[e]],
            scalarrel=False,
            confinement=conf,
            wf_confinement=wf_conf,
            verbose=False,
            txt='Slako.txt',
            mix=0.05,
            maxiter=1000,
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



if __name__ == '__main__':
    # elements = ['In', 'As']
    elements = ['As', 'In']
    hubbardvalues = get_hubbard_values(elements=elements)
    # eigenvalues = {
    #     'In': {'5d': 0.135383, '5p': -0.092539, '5s': -0.30165},
    #     'As': {'4d': 0.135383, '4p': -0.092539, '4s': -0.30165}
    # }
    # hubbardvalues = {
    #     'In': {'d': 0.156519, 'p': 0.189913, 's': 0.257192},
    #     'As': {'d': 0.127856, 'p': 0.271613, 's': 0.330013},
    # }
    # valences = {'In': ['4d10', '5s2', '5p1'], 'As': ['3d10', '4s2', '4p3']}
    # valences = _get_valences(elements)
    # rcovs = {e: covalent_radii[atomic_numbers[e]] / Bohr for e in elements}
    rcovs = {'In': 4.8, 'As': 4.4}
    sigmas = {'In': 13.2, 'As': 5.6}
    params = {}

    for e in elements:
        params[f'{e}_r0'] = 2 * rcovs[e]
        params[f'{e}_sigma'] = sigmas[e]

    opt = [6.77721995, 3.6238089, 4.91739859, 3.03960597]
    params = dict(zip(params.keys(), opt))

    # params = {'As_r0': 6.787252132356805, 'As_sigma': 3.151596155752452, 'In_r0': 9.661769132437286, 'In_sigma': 14.889611814185765}

    print(f'{hubbardvalues = }')
    print(f'{params =}')
    print('Generate Slako Files')
    write_slako(
        elements=elements,
        params=params,
        hubbardvalues=hubbardvalues,
        slako_dir='slako_test',
    )
    print('Done!')
