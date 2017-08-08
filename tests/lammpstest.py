# Richard Darst, November 2009

# Richard Darst, November 2009

from saiga12.lammps import Lammps


def make_system():
    S = Lammps()
    
    import string
    script = string.Template("""
    units           lj
    atom_style      atomic
    lattice fcc     $density
    region          box block 0 $L 0 $L 0 $L
    create_box      1 box
    create_atoms    1 box
    mass            1 1.0
    velocity        all create $temp 87287 loop geom
    
    pair_style      lj/cut 2.5
    pair_coeff      1 1 1.0 1.0 2.5  # typeI typeJ epsilon sigma cutoff
    
    #neighbor        0.3 bin
    #neigh_modify    delay 0 every 20 check no
    fix             1 all nve
    
    """)
    density = .84
    temp = 1.44
    L = 10
    script = script.substitute(locals())
    S.commands(script.split('\n'))
    return S
    print S.density


S = make_system()
