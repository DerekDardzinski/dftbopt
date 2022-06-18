""" A test with relatively steep confinements. """
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement

conf_dens = PowerConfinement(r0=3., s=12.)
conf_wf = {'1s': PowerConfinement(r0=2., s=12.)}

atom = AtomicDFT('H',
                 xc='LDA',
                 configuration='1s1',
                 valence=['1s'],
                 scalarrel=False,
                 confinement=conf_dens,
                 wf_confinement=conf_wf,
                 )
atom.run()
