from ase.units import Bohr
from ase.data import atomic_numbers, covalent_radii
from itertools import combinations_with_replacement
from pymatgen.core.periodic_table import Element
from dd_hotcent.hotcent.atomic_dft import AtomicDFT
from dd_hotcent.hotcent.slako import SlaterKosterTable
from dd_hotcent.hotcent.confinement import PowerConfinement
import os

# Define standard, rule-of-thumb confinement potentials
xc = 'GGA_X_PBE+GGA_C_PBE'
elements = ['In', 'As']
valences = {
    e: Element(e).electronic_structure.split('.')[1:] for e in elements
}
rcovs = [covalent_radii[atomic_numbers[element]] / Bohr for element in elements]
confs = [PowerConfinement(r0=3 * r, s=2) for r in rcovs]
wf_confs = [
    {v: PowerConfinement(r0=2 * r, s=2) for v in valences[e]}
    for e, r in zip(elements, rcovs)
]

atoms = {}
for element, conf, wf_conf in zip(elements, confs, wf_confs):
    atom = AtomicDFT(
        element,
        xc=xc,
        configuration=Element(element).electronic_structure.replace('.', ' '),
        valence=[v[:2] for v in valences[element]],
        scalarrel=False,
        confinement=conf,
        wf_confinement=wf_conf,
    )
    atom.run()
    atoms[element] = atom


U_values = {}
for element, conf, wf_conf in zip(elements, confs, wf_confs):
    atom_hubbard = AtomicDFT(
        element,
        xc=xc,
        configuration=Element(element).electronic_structure.replace('.', ' '),
        valence=[v[:2] for v in valences[element]],
        scalarrel=False,
        # Add a very weak confinement potential to aid anion convergence:
        confinement=PowerConfinement(r0=40., s=4),
    )
    U_tmp = {}
    for v in valences[element]:
        U = atom_hubbard.get_hubbard_value(
            v[:2],
            scheme='central',
            maxstep=1,
        )
        U_tmp[v[1]] = U

    U_values[element] = U_tmp

rmin, dr, N = 0.4, 0.02, 600
atom_combos = combinations_with_replacement(elements, 2)

for atom_combo in atom_combos:
    sk = SlaterKosterTable(atoms[atom_combo[0]], atoms[atom_combo[1]])
    sk.run(rmin, dr, N, superposition='density', xc=xc, wflimit=1e-4)

    # Write the Slater-Koster tables to file (without two-body repulsion at this point).
    # This file also stores the eigenvalues, Hubbardvalues, occupations, as well as the
    # so-called spin-polarization error (the magnetization energy of the atom, which we
    # don't need to consider here).
    if atom_combo[0] == atom_combo[1]:
        hub = U_values[atom_combo[0]]
        occu = {v[:2]: int(v[2:]) for v in valences[atom_combo[0]]}
        eig = dict(atoms[atom_combo[0]].get_valence_energies())
        print(hub)
        print(occu)
        print(eig)
        sk.write(
            os.path.join(
                './slako', f'{atom_combo[0]}-{atom_combo[1]}_no_repulsion.skf'),
            eigenvalues=eig,
            hubbardvalues=hub,
            occupations=occu,
            spe=0.
        )
    else:
        sk.write(
            os.path.join(
                './slako', f'{atom_combo[0]}-{atom_combo[1]}_no_repulsion.skf'),
            pair=(atom_combo[0], atom_combo[1]),
            spe=0.
        )
        sk.write(
            os.path.join(
                './slako', f'{atom_combo[1]}-{atom_combo[0]}_no_repulsion.skf'),
            pair=(atom_combo[1], atom_combo[0]),
            spe=0.
        )
