import re
import os
import numpy as np

def convertatomicnumber(X):
    if X == 1:
        return 'H'
    elif X == 6:
        return 'C'
    elif X == 7:
        return 'N'
    elif X == 8:
        return 'O'

# find_normal_modes
#               parses a guassian output and returns the normal modes
def find_normal_modes(file_string):
    # We first extract the general frequency information into a
    # vector of up to three modes per element
    p_float = re.compile("-?\\d+\\.\\d+")
    p_freq_block = re.compile("(Frequencies(.|\n)*?(?=(Thermo|Frequ)))")
    matches = re.findall(p_freq_block, file_string)

    # From each of these elements, we extract the force constants
    p_frc_const = re.compile("Frc.*")
    force_constants = []
    for i in matches:
        m = re.search(p_frc_const, i[0])
        numbers = re.findall(p_float, m.group(0))
        for n in numbers:
            force_constants.append(float(n))

    # From each of these elements, we extract the force constants
    p_frqs = re.compile("Frequencies.*")
    freqs_constants = []
    for i in matches:
        m = re.search(p_frqs, i[0])
        numbers = re.findall(p_float, m.group(0))
        for n in numbers:
            freqs_constants.append(float(n))

    for i,j in enumerate(freqs_constants):
        if j < 400.0:
            force_constants[i] = 0.1

    print(force_constants)

    # To get the normal modes, we recognize that the rows we want from the output
    # are the only ones that begin with an integer and end with a double
    p_atom = re.compile("(\n\\s+\\d+\\s+.*-?\\d+\\.\\d+)")
    atoms = re.findall(p_atom, matches[0][0])
    number_atoms = len(atoms)
    normals_of_atoms = [[] for _ in range(number_atoms)]
    for i in matches:
        atoms = re.findall(p_atom, i[0])
        for atom_number, values in enumerate(atoms):
            modes = re.findall(p_float, values)
            mode_values = []
            for i in modes:
                mode_values.append(float(i))
            normals_of_atoms[atom_number].extend(mode_values)

    # Right now we have two data object. The first one is a simple list
    # list of the force constants per frequency. The second is a list
    # of lists representing [atoms][normal values]. We now want to merge
    # the two such to get a list of of lists, in which the inner list is
    # is composed of force_constants, x_atom1, y_atom1, z_atom1, x_atom2.
    # and the outer represents each frequency.
    number_force_constants = len(force_constants)
    normal_mods = [[] for _ in range(number_force_constants)]
    for k, f_const in enumerate(force_constants):
        normal_mods[k].append(f_const)
    for i in range(number_atoms):
        x = []
        y = []
        z = []
        for j in range(len(normals_of_atoms[i])):
            if (j % 3) == 0:
                x.append(normals_of_atoms[i][j])
            if (j % 3) == 1:
                y.append(normals_of_atoms[i][j])
            if (j % 3) == 2:
                z.append(normals_of_atoms[i][j])
        for j in range(len(x)):
            normal_mods[j].extend([x[j], y[j], z[j]])
    return normal_mods, freqs_constants

wkdir = '/home/jujuman/Research/GDB-11-wB97X-6-31gd/dnnts_dissociation/scans_cc_bonds_dft/double/'

fpf = 'ethene_CC_disso' #Filename prefix

TSS = 10
LOT='WB97X/6-31g*' # High level of theory
rdm='uniform' #Random dist
type='nmrandom'
Temp='1500.0'
SCF='Tight'

files = os.listdir(wkdir)
files.sort()

xyz = open(wkdir + 'structures.xyz', 'w')

freq1 = []

Nc = 0
Nd = 0
for k in files:
    if 'structure' in k and '.out' in k:
        file = wkdir + k

        f = open(wkdir + 'inputs/' + fpf + '-' + str(Nc) + '.ipt', 'w')

        print(k)

        fi = open(file,'r').read()

        br = re.compile("(?:Input orientation:).*\n.*\n.*\n.*\n.*\n([\s\S]+?)(?:-----)")
        cr = re.compile("\s+?\d+?\s+?(\d+?)\s+?\d+?\s+?(\S+?)\s+?(\S+?)\s+?(\S+?)\s+?(?:\s|\n)")

        m1 = re.findall(br,fi)

        Na = len(re.findall(cr, m1[0]))

        #---------- Write Input Variables ------------
        dfname=fpf + '-' + str(Nc) + '_train.dat'
        vdfname=fpf + '-' + str(Nc) + '_valid.dat'
        edfname=fpf + '-' + str(Nc) + '_test.dat'

        V = 6
        if m1 is 2:
            V = 5

        DOF = (3 * Na - V)
        Nd += TSS * DOF

        f.write ('TSS=' + str(int(TSS * DOF)) + ' \n')
        f.write ('VSS=' + str(int((TSS * DOF)/10.0)) + ' \n')
        f.write ('ESS=' + str(int((TSS * DOF)/10.0)) + ' \n')
        #f.write ('TSS=' + str(0) + ' \n')
        #f.write ('VSS=' + str(0) + ' \n')
        #f.write ('ESS=' + str(20) + ' \n')
        f.write ('LOT=' + LOT + ' \n')
        f.write ('rdm=' + rdm + '\n')
        f.write ('type=' + type + '\n')
        f.write ('Temp=' + Temp + '\n')
        f.write ('mem=' + '1024' + '\n')
        f.write ('SCF=' + SCF + '\n')
        f.write ('dfname=' + dfname + ' \n')
        f.write ('vdfname=' + vdfname + ' \n')
        f.write ('edfname=' + edfname + ' \n')
        f.write ('optimize=0 \n')
        f.write ('frequency=0 \n')

        f.write ('\n\n')
        f.write ('$coordinates\n')
        for i,b in enumerate(m1):
            m2 = re.findall(cr,b)
            xyz.write('\n'+str(len(m2))+'\n')
            for j, c in enumerate(m2):
                at = convertatomicnumber(int(c[0]))
                xyz.write(at + ' ' + c[1] + ' ' + c[2] + ' ' + c[3] + '\n')
                f.write(' ' + at + ' ' + at + ' ' + c[1] + ' ' + c[2] + ' ' + c[3] + '\n')
        f.write ('&\n\n')

        f.write ('$connectivity\n')
        f.write (' NONE\n')
        f.write ('&\n\n')

        nms,frq = find_normal_modes(fi)
        nms = np.array(nms)
        freq1.append(float(frq[4]))
        f.write ('$normalmodes\n')
        for m in nms:
            f.write('FRCCNST='+ "{:.4e}".format(float(m[0])) +' {\n')
            modes = m[1:].flatten().reshape(Na,3)
            for x in modes:
                f.write(' ' + "{:.4e}".format(float(x[0])) + ' '
                            + "{:.4e}".format(float(x[1])) + ' '
                            + "{:.4e}".format(float(x[2])) + '\n')
            f.write('}\n')
        f.write ('&\n\n')

        f.close()
        Nc = Nc + 1

#import matplotlib.pyplot as plt

#plt.plot(freq1)

#plt.xlabel('step')
#plt.ylabel('freq')

#plt.show()
