"""
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

Program: generate.py

Generates a simple initial ssDNA or dsDNA configuration from a given sequence.
For dsDNA the sequence should be preceded by the 'DOUBLE' keyword.

Usage: 
$$ python generate.py box_offset box_length sequence_file [topology (strand OR ring)] [twist ((p_0-p)/p) for strand OR Delta_Lk for ring] 
"""

# for python2/3 compatibility
from __future__ import print_function
#!/usr/bin/env python

# number of ACGT base type groups
# default value 16 protects against asymmetric base pairing over 1.5 turns on the helix
N_BASE_TYPE_GROUPS = 16

"""
Import basic modules
"""
import sys, os, timeit

from timeit import default_timer as timer
start_time = timer()
"""
Try to import numpy; if failed, import a local version mynumpy
which needs to be provided
"""
try:
    import numpy as np
except:
    print("numpy not found. Exiting.", file=sys.stderr)
    sys.exit(1)

"""
Check that the required arguments (box offset and size in simulation units
and the sequence file were provided
"""
try:
    box_offset = float(sys.argv[1])
    box_length = float(sys.argv[2])
    infile = sys.argv[3]
    # default is strand with no twist, i.e. pitch = 10.5 bp
    if len(sys.argv) == 4:
        topo = 'strand'
        twist = 0
    # default is no twist, i.e. pitch = 10.5 bp
    elif len(sys.argv) == 5:
        topo = sys.argv[4]
        twist = 0
    elif len(sys.argv) == 6:
        topo = sys.argv[4]
        twist = float(sys.argv[5])

except:
    print("Usage: %s <%s> <%s> <%s> <%s> <%s> " % (sys.argv[0], \
        "box offset", "box length", "file with sequences", "[topology (strand OR ring)]", "[twist ((p_0-p)/p) for strand OR Delta_Lk for ring]"), file=sys.stderr)
    sys.exit(1)
box = np.array ([box_length, box_length, box_length])

"""
Try to open the file and fail gracefully if file cannot be opened
"""
try:
    inp = open (infile, 'r')
    inp.close()
except:
    print( "Could not open file '%s' for reading. Aborting." % infile, file=sys.stderr)
    sys.exit(2)

# return parts of a string
def partition(s, d):
    if d in s:
        sp = s.split(d, 1)
        return sp[0], d, sp[1]
    else:
        return s, "", ""

"""
Define the model constants
"""
# set model constants
PI = np.pi
POS_BASE =  0.4
POS_BACK = -0.4
EXCL_RC1 =  0.711879214356
EXCL_RC2 =  0.335388426126
EXCL_RC3 =  0.52329943261

"""
Define auxiliary variables for the construction of a helix
"""
# center of the double strand
CM_CENTER_DS = POS_BASE + 0.2

# ideal distance between base sites of two nucleotides
# which are to be base paired in a duplex
BASE_BASE = 0.3897628551303122

# cutoff distance for overlap check
RC2 = 16

# squares of the excluded volume distances for overlap check
RC2_BACK = EXCL_RC1**2
RC2_BASE = EXCL_RC2**2
RC2_BACK_BASE = EXCL_RC3**2

# enumeration to translate from letters to numbers and vice versa
number_to_base = {1 : 'A', 2 : 'C', 3 : 'G', 4 : 'T'}
base_to_number = {'A' : 1, 'a' : 1, 'C' : 2, 'c' : 2,
                  'G' : 3, 'g' : 3, 'T' : 4, 't' : 4}

# auxiliary arrays
positions = []
a1s = []
a3s = []
quaternions = []

newpositions = []
newa1s = []
newa3s = []

basetype  = []
strandnum = []

bonds = []

"""
Convert local body frame to quaternion DOF
"""
def exyz_to_quat (mya1, mya3):

    mya2 = np.cross(mya3, mya1)
    myquat = [1,0,0,0]

    q0sq = 0.25 * (mya1[0] + mya2[1] + mya3[2] + 1.0)
    q1sq = q0sq - 0.5 * (mya2[1] + mya3[2])
    q2sq = q0sq - 0.5 * (mya1[0] + mya3[2])
    q3sq = q0sq - 0.5 * (mya1[0] + mya2[1])

    # some component must be greater than 1/4 since they sum to 1
    # compute other components from it

    if q0sq >= 0.25:
        myquat[0] = np.sqrt(q0sq)
        myquat[1] = (mya2[2] - mya3[1]) / (4.0*myquat[0])
        myquat[2] = (mya3[0] - mya1[2]) / (4.0*myquat[0])
        myquat[3] = (mya1[1] - mya2[0]) / (4.0*myquat[0])
    elif q1sq >= 0.25:
        myquat[1] = np.sqrt(q1sq)
        myquat[0] = (mya2[2] - mya3[1]) / (4.0*myquat[1])
        myquat[2] = (mya2[0] + mya1[1]) / (4.0*myquat[1])
        myquat[3] = (mya1[2] + mya3[0]) / (4.0*myquat[1])
    elif q2sq >= 0.25:
        myquat[2] = np.sqrt(q2sq)
        myquat[0] = (mya3[0] - mya1[2]) / (4.0*myquat[2])
        myquat[1] = (mya2[0] + mya1[1]) / (4.0*myquat[2])
        myquat[3] = (mya3[1] + mya2[2]) / (4.0*myquat[2])
    elif q3sq >= 0.25:
        myquat[3] = np.sqrt(q3sq)
        myquat[0] = (mya1[1] - mya2[0]) / (4.0*myquat[3])
        myquat[1] = (mya3[0] + mya1[2]) / (4.0*myquat[3])
        myquat[2] = (mya3[1] + mya2[2]) / (4.0*myquat[3])

    norm = 1.0/np.sqrt(myquat[0]*myquat[0] + myquat[1]*myquat[1] + \
			  myquat[2]*myquat[2] + myquat[3]*myquat[3])
    myquat[0] *= norm
    myquat[1] *= norm
    myquat[2] *= norm
    myquat[3] *= norm

    return np.array([myquat[0],myquat[1],myquat[2],myquat[3]])

"""
Adds a strand to the system by appending it to the array of previous strands
"""
def add_strands (mynewpositions, mynewa1s, mynewa3s):
    overlap = False

    # This is a simple check for each of the particles where for previously
    # placed particles i we check whether it overlaps with any of the
    # newly created particles j

    print( "## Checking for overlaps", file=sys.stdout)

    for i in range(len(positions)):

        p = positions[i]
        pa1 = a1s[i]

        for j in range (len(mynewpositions)):

            q = mynewpositions[j]
            qa1 = mynewa1s[j]

            # skip particles that are anyway too far away
            dr = p - q
            dr -= box * np.rint(dr / box)
            if np.dot(dr, dr) > RC2:
                continue

            # base site and backbone site of the two particles
            p_pos_back = p + pa1 * POS_BACK
            p_pos_base = p + pa1 * POS_BASE
            q_pos_back = q + qa1 * POS_BACK
            q_pos_base = q + qa1 * POS_BASE

            # check for no overlap between the two backbone sites
            dr = p_pos_back - q_pos_back
            dr -= box * np.rint(dr / box)
            if np.dot(dr, dr) < RC2_BACK:
                overlap = True

            # check for no overlap between the two base sites
            dr = p_pos_base -  q_pos_base
            dr -= box * np.rint(dr / box)
            if np.dot(dr, dr) < RC2_BASE:
                overlap = True

            # check for no overlap between backbone site of particle p
            # with base site of particle q
            dr = p_pos_back - q_pos_base
            dr -= box * np.rint (dr / box)
            if np.dot(dr, dr) < RC2_BACK_BASE:
                overlap = True

            # check for no overlap between base site of particle p and
            # backbone site of particle q
            dr = p_pos_base - q_pos_back
            dr -= box * np.rint (dr / box)
            if np.dot(dr, dr) < RC2_BACK_BASE:
                overlap = True

            # exit if there is an overlap
            if overlap:
                return False

    # append to the existing list if no overlap is found
    if not overlap:

        for p in mynewpositions:
            positions.append(p)
        for p in mynewa1s:
            a1s.append (p)
        for p in mynewa3s:
            a3s.append (p)
        # calculate quaternion from local body frame and append
        for ia in range(len(mynewpositions)):
            mynewquaternions = exyz_to_quat(mynewa1s[ia],mynewa3s[ia])
            quaternions.append(mynewquaternions)

    return True

"""
Calculate angle of rotation site to site
"""
def get_angle(bp):
    #n, default number of bases per turn
    n = 10.55
    np = n
    nm = n
    found = False
    angle = 0
    while found == False:

        turnsp = bp/np
        turnsm = bp/nm
    
        diffp = abs(turnsp - round(turnsp))
        diffm = abs(turnsm - round(turnsm))

        if diffp < 0.03:
            found = True
            turnsp = round(turnsp)+int(twist)
            angle = (360*turnsp)/bp
            angle = round(angle,2)
            break

        elif diffm < 0.03:
            found = True
            turnsm = round(turnsm)+int(twist)
            angle = (360*turnsm)/bp
            angle = round(angle,2)
            break

        else:
            np += 0.001
            nm -= 0.001

    return angle

"""
Returns the rotation matrix defined by an axis and angle
"""
def get_rotation_matrix(axis, anglest, nbp=0):
    # The argument anglest can be either an angle in radiants
    # (accepted types are float, int or np.float64 or np.float64)
    # or a tuple [angle, units] where angle is a number and
    # units is a string. It tells the routine whether to use degrees,
    # radiants (the default) or base pairs turns.
    if not isinstance (anglest, (np.float64, np.float32, float, int)):
        if len(anglest) > 1:
            if anglest[1] in ["degrees", "deg", "o"]:
                angle = (np.pi / 180.) * (anglest[0])
            elif anglest[1] in ["bp"]:

                # strand: angle is 360/pitch
                # pitch is p = 10.5 * (1-twist)
                if nbp == 0:
                    ang = 360./(10.5*(1-twist))

                # ring: angle depends on Delta_Lk = Lk - Lk_0
                else:
                    ang = get_angle(nbp)

                angle = int(anglest[0]) * (np.pi / 180.) * (ang)
            else:
                angle = float(anglest[0])
        else:
            angle = float(anglest[0])
    else:
        angle = float(anglest) # in degrees (?)

    axis = np.array(axis)
    axis /= np.sqrt(np.dot(axis, axis))

    ct = np.cos(angle)
    st = np.sin(angle)
    olc = 1. - ct
    x, y, z = axis

    return np.array([[olc*x*x+ct, olc*x*y-st*z, olc*x*z+st*y],
                    [olc*x*y+st*z, olc*y*y+ct, olc*y*z-st*x],
                    [olc*x*z-st*y, olc*y*z+st*x, olc*z*z+ct]])

"""
Generates the position and orientation vectors of a
(single or double) strand from a sequence string
"""
def generate_strand(bp, sequence=None, start_pos=np.array([0, 0, 0]), \
	  dir=np.array([0, 0, 1]), perp=False, double=True, rot=0.):
    # generate empty arrays
    mynewpositions, mynewa1s, mynewa3s = [], [], []

    # cast the provided start_pos array into a numpy array
    start_pos = np.array(start_pos, dtype=float)

    # overall direction of the helix
    dir = np.array(dir, dtype=float)
    if sequence == None:
        sequence = np.random.randint(1, 5, bp)

    # the elseif here is most likely redundant
    elif len(sequence) != bp:
        n = bp - len(sequence)
        sequence += np.random.randint(1, 5, n)
        print( "sequence is too short, adding %d random bases" % n, file=sys.stderr)

    # normalize direction
    dir_norm = np.sqrt(np.dot(dir,dir))
    if dir_norm < 1e-10:
        print( "direction must be a valid vector, defaulting to (0, 0, 1)", file=sys.stderr)
        dir = np.array([0, 0, 1])
    else: dir /= dir_norm

    # find a vector orthogonal to dir to act as helix direction,
    # if not provided switch off random orientation
    if perp is None or perp is False:
        v1 = np.random.random_sample(3)
        # comment in to suppress randomised base vector
        v1 = [1,0,0]
        v1 -= dir * (np.dot(dir, v1))
        v1 /= np.sqrt(sum(v1*v1))
    else:
        v1 = perp;

    # generate rotational matrix representing the overall rotation of the helix
    R0 = get_rotation_matrix(dir, rot)

    # rotation matrix corresponding to one step along the helix
    R = get_rotation_matrix(dir, [1, "bp"])

    # set the vector a1 (backbone to base) to v1
    a1 = v1

    # apply the global rotation to a1
    a1 = np.dot(R0, a1)

    # set the position of the fist backbone site to start_pos
    rb = np.array(start_pos)

    # set a3 to the direction of the helix
    a3 = dir
    for i in range(bp):
    # work out the position of the centre of mass of the nucleotide
        rcom = rb - CM_CENTER_DS * a1

        # append to newpositions
        mynewpositions.append(rcom)
        mynewa1s.append(a1)
        mynewa3s.append(a3)

        # if we are not at the end of the helix, we work out a1 and rb for the
        # next nucleotide along the helix
        if i != bp - 1:
            a1 = np.dot(R, a1)
            rb += a3 * BASE_BASE

    # if we are working on a double strand, we do a cycle similar
    # to the previous one but backwards
    if double == True:
        a1 = -a1
        a3 = -dir
        R = R.transpose()
        for i in range(bp):
            rcom = rb - CM_CENTER_DS * a1
            mynewpositions.append (rcom)
            mynewa1s.append (a1)
            mynewa3s.append (a3)
            a1 = np.dot(R, a1)
            rb += a3 * BASE_BASE

    assert (len (mynewpositions) > 0)

    return [mynewpositions, mynewa1s, mynewa3s]

"""
Generates the position and orientation vectors of a 
(single or double) ring from a sequence string
"""
def generate_ring(bp, sequence=None, start_pos=np.array([0, 0, 0]), \
          dir=np.array([0, 0, 1]), perp=False, double=True, rot=0.):
    # generate empty arrays
    mynewpositions, mynewa1s, mynewa3s = [], [], []

    # cast the provided start_pos array into a numpy array
    start_pos = np.array(start_pos, dtype=float)

    # overall direction of the helix
    dir = np.array(dir, dtype=float)
    if sequence == None:
        sequence = np.random.randint(1, 5, bp)

    # the elseif here is most likely redundant 
    elif len(sequence) != bp:
        n = bp - len(sequence)
        sequence += np.random.randint(1, 5, n)
        print("sequence is too short, adding %d random bases" % n, file=sys.stderr)

    # normalize direction
    dir_norm = np.sqrt(np.dot(dir,dir))
    if dir_norm < 1e-10:
        print("direction must be a valid vector,\
                defaulting to (0, 0, 1)", file=sys.stderr)
        dir = np.array([0, 0, 1])
    else: dir /= dir_norm

    # find a vector orthogonal to dir to act as helix direction,
    # if not provided switch off random orientation
    if perp is None or perp is False:
        v1 = np.random.random_sample(3)
        # comment in to suppress randomised base vector
        v1 = [1,0,0]
        v1 -= dir * (np.dot(dir, v1))
        v1 /= np.sqrt(sum(v1*v1))
    else:
        v1 = perp;


    #work out v2 (lab frame) (Orthogonal vector to dir and v1) - probably need to translate by R
    v2 = np.cross(dir,v1)
    v2 /= np.sqrt(np.dot(v2,v2))

    # generate rotational matrix representing the overall rotation of the helix
    R0 = get_rotation_matrix(dir, rot)

    # rotation matrix corresponding to one step along the helix
    R = get_rotation_matrix(dir, [1, "bp"],bp)

    #rotation matrix for ring creation, rotation around lab v2 by 2 pi/ len(bp)

    ang =( 2*PI)/(bp)

    Rv2 = get_rotation_matrix(v2,ang)

    # set the vector a1 (backbone to base) to v1 
    a1 = v1

    # apply the global rotation to a1 
    a1 = np.dot(R0, a1)

    # set the position of the fist backbone site to start_pos
    rb = np.array(start_pos)

    # set a3 to the direction of the helix
    a3 = dir
    for i in range(bp):

        R = get_rotation_matrix(a3,[1,"bp"],bp)

        rcom = rb - CM_CENTER_DS * a1

        # append to newpositions
        mynewpositions.append(rcom)
        mynewa1s.append(a1)
        mynewa3s.append(a3)
        if i != bp -1 :
            rb += a3*BASE_BASE
            #calculate a3p, a1p by rotating a3, a1 around v2
            a1 = np.dot(R0,a1)
            a3p = np.dot(Rv2,a3)
            a1p = np.dot(Rv2,a1)
            a1pp = np.dot(R,a1p)
            a3 = a3p
            a1 = a1pp

    if double == True:
        cor=(bp+1)/bp
        cor=1
        ang =(-2*PI)/(bp)*cor

        a1 = -a1
        a3 = -a3
        R = R.transpose()
        Rv2 = get_rotation_matrix(v2,ang)

        for i in range(bp):

            rcom = rb - CM_CENTER_DS * a1

            mynewpositions.append (rcom)
            mynewa1s.append (a1)
            mynewa3s.append(a3)

            a3p = np.dot(Rv2,a3)
            a1p = np.dot(Rv2,a1)
            a1pp = np.dot(R,a1p)
            a3 = a3p
            a1 = a1pp
            rb += a3 * BASE_BASE
            R = get_rotation_matrix(a3,[1,"bp"],bp)

    assert (len (mynewpositions) > 0)

    return [mynewpositions, mynewa1s, mynewa3s]

"""
Main function for this script.
Reads a text file with the following format:
- Each line contains the sequence for a single strand (A,C,G,T)
- Lines beginning with the keyword 'DOUBLE' produce double-stranded DNA

Ex: Two ssDNA (single stranded DNA)
ATATATA
GCGCGCG

Ex: Two strands, one double stranded, the other single stranded.
DOUBLE AGGGCT
CCTGTA

"""

def read_strands(filename):
    try:
        infile = open (filename)
    except:
        print( "Could not open file '%s'. Aborting." % filename, file=sys.stderr )
        sys.exit(2)

    # This block works out the number of nucleotides and strands by reading
    # the number of non-empty lines in the input file and the number of letters,
    # taking the possible DOUBLE keyword into account.
    nstrands, nnucl, nbonds = 0, 0, 0
    lines = infile.readlines()
    for line in lines:
        line = line.upper().strip()
        if len(line) == 0:
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1]
            length = len(line)
            print( "## Found duplex of %i base pairs" % length, file=sys.stdout)
            nnucl += 2*length
            nstrands += 2
            if topo == 'ring':
                nbonds+= 2*length
            else:
                nbonds += (2*length-2)
        else:
            line = line.split()[0]
            length = len(line)
            print( "## Found single strand of %i bases" % length, file=sys.stdout)
            nnucl += length
            nstrands += 1
            if topo == 'ring':
                nbonds =+ length
            else:
                nbonds += length-1

    # rewind the sequence input file
    infile.seek(0)

    print( "## nstrands, nnucl = ", nstrands, nnucl, file=sys.stdout)

    # generate the data file in LAMMPS format
    try:
        out = open ("data.oxdna", "w")
    except:
        print( "Could not open data file for writing. Aborting.", file=sys.stderr)
        sys.exit(2)

    lines = infile.readlines()
    nlines = len(lines)
    i = 1
    myns = 0
    noffset = 1

    for line in lines:
        line = line.upper().strip()

        # skip empty lines
        if len(line) == 0:
            i += 1
            continue

        # block for duplexes: last argument of the generate function
        # is set to 'True'
        if line[:6] == 'DOUBLE':
            line = line.split()[1]
            length = len(line)
            seq = [(base_to_number[x]) for x in line]

            # modify default type for unique base pairing
            group = 0
            for s in range(len(seq)):
                seq[s] += group*4
                group += 1
                # reset according to maximal number of base type groups
                if group == N_BASE_TYPE_GROUPS:
                    group = 0

            myns += 1

            for b in range(length):
                basetype.append(seq[b])
                strandnum.append(myns)

            for b in range(length-1):
                bondpair = [noffset + b, noffset + b + 1]
                bonds.append(bondpair)

            if topo== 'ring':
                bondpair = [noffset,noffset+length-1]
                bonds.append(bondpair)

            noffset += length

            # create the complementary sequence of the second strand
            
            seq2 = seq
            for s in range(len(seq2)):
                if seq2[s]%4 == 1:
                    seq2[s] += 3
                elif seq2[s]%4 == 2:
                    seq2[s] += 1
                elif seq2[s]%4 == 3:
                    seq2[s] -= 1
                elif seq2[s]%4 == 0:
                    seq2[s] -= 3

            seq2 = seq2[::-1]
           
            myns += 1

            for b in range(length):
                basetype.append(seq2[b])
                strandnum.append(myns)

            for b in range(length-1):
                bondpair = [noffset + b, noffset + b + 1]
                bonds.append(bondpair)

            if topo == 'ring':
                bondpair = [noffset, noffset + length-1]
                bonds.append(bondpair)

            noffset += length

            print( "## Created duplex of %i bases" % (2*length), file=sys.stdout)

            # generate random position of the first nucleotide
            com = box_offset + np.random.random_sample(3) * box
            # comment out to randomise
            com = [0,0,0]

            # generate the random direction of the helix
            axis = np.random.random_sample(3)
            # comment out to randomise
            axis = [0,0,1]
            axis /= np.sqrt(np.dot(axis, axis))

            # use the generate function defined above to create 
            # the position and orientation vector of the strand
            if topo == 'ring':
                newpositions, newa1s, newa3s = generate_ring(len(line), \
                    sequence=seq, dir=axis, start_pos=com, double=True)
            else:
                newpositions, newa1s, newa3s = generate_strand(len(line), \
                    sequence=seq, dir=axis, start_pos=com, double=True)

            # generate a new position for the strand until it does not overlap
            # with anything already present
            start = timer()
            while not add_strands(newpositions, newa1s, newa3s):
                com = box_offset + np.random.random_sample(3) * box
                axis = np.random.random_sample(3)
                axis /= np.sqrt(np.dot(axis, axis))
                if topo == 'ring':
                    newpositions, newa1s, newa3s = generate_ring(len(line), \
                        sequence=seq, dir=axis, start_pos=com, double=True)
                else:
                    newpositions, newa1s, newa3s = generate_strand(len(line), \
                        sequence=seq, dir=axis, start_pos=com, double=True)
                print( "## Trying %i" % i, file=sys.stdout)
            end = timer()
            print( "## Added duplex of %i bases (line %i/%i) in %.2fs, now at %i/%i" % \
				      (2*length, i, nlines, end-start, len(positions), nnucl), file=sys.stdout)

        # block for single strands: last argument of the generate function
        # is set to 'False'
        else:
            length = len(line)
            seq = [(base_to_number[x]) for x in line]

            myns += 1
            for b in range(length):
                basetype.append(seq[b])
                strandnum.append(myns)

            for b in range(length-1):
                bondpair = [noffset + b, noffset + b + 1]
                bonds.append(bondpair)

            if topo == 'ring':
                bondpair = [noffset, noffset + length-1]
                bonds.append(bondpair)

            noffset += length

            # generate random position of the first nucleotide
            com = box_offset + np.random.random_sample(3) * box

            # generate the random direction of the helix
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))

            print("## Created single strand of %i bases" % length, file=sys.stdout)
            if topo == 'ring':
                newpositions, newa1s, newa3s = generate_ring(length, \
                              sequence=seq, dir=axis, start_pos=com, double=False)
            else:
                newpositions, newa1s, newa3s = generate_strand(length, \
                               sequence=seq, dir=axis, start_pos=com, double=False)
            start = timer()
            while not add_strands(newpositions, newa1s, newa3s):
                com = box_offset + np.random.random_sample(3) * box
                # comment out to randomise
                com = [0,0,0]
                axis = np.random.random_sample(3)
                # comment out to randomise
                axis = [0,0,1]
                axis /= np.sqrt(np.dot(axis, axis))
                if topo == 'ring':
                    newpositions, newa1s, newa3s = generate_ring(length, \
                          sequence=seq, dir=axis, start_pos=com, double=False)

                else:
                    newpositions, newa1s, newa3s = generate_strand(length, \
                          sequence=seq, dir=axis, start_pos=com, double=False)
                print >> sys.stdout, "## Trying  %i" % (i)
            end = timer()
            print( "## Added single strand of %i bases (line %i/%i) in %.2fs, now at %i/%i" % \
				      (length, i, nlines, end-start,len(positions), nnucl), file=sys.stdout)

        i += 1

    # sanity check
    if not len(positions) == nnucl:
        print( len(positions), nnucl )
        raise AssertionError

    out.write('# LAMMPS data file\n')
    out.write('%d atoms\n' % nnucl)
    out.write('%d ellipsoids\n' % nnucl)
    out.write('%d bonds\n' % nbonds)
    out.write('\n')
    out.write('%d atom types\n' % (N_BASE_TYPE_GROUPS*4))
    out.write('1 bond types\n')
    out.write('\n')
    out.write('# System size\n')
    out.write('%f %f xlo xhi\n' % (box_offset,box_offset+box_length))
    out.write('%f %f ylo yhi\n' % (box_offset,box_offset+box_length))
    out.write('%f %f zlo zhi\n' % (box_offset,box_offset+box_length))

    out.write('\n')
    out.write('Masses\n')
    out.write('\n')
    for i in range(N_BASE_TYPE_GROUPS*4):
        out.write('%d 3.1575\n' % (i+1))

    # for each nucleotide print a line under the headers
    # Atoms, Velocities, Ellipsoids and Bonds
    out.write('\n')
    out.write(\
      '# Atom-ID, type, position, molecule-ID, ellipsoid flag, density\n')
    out.write('Atoms\n')
    out.write('\n')

    for i in range(nnucl):
        out.write('%d %d %22.15le %22.15le %22.15le %d 1 1\n' \
            % (i+1, basetype[i], positions[i][0], positions[i][1], positions[i][2], strandnum[i]))

    out.write('\n')
    out.write('# Atom-ID, translational velocity, angular momentum\n')
    out.write('Velocities\n')
    out.write('\n')

    for i in range(nnucl):
        out.write("%d %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le\n" \
            % (i+1,0.0,0.0,0.0,0.0,0.0,0.0))

    out.write('\n')
    out.write('# Atom-ID, shape, quaternion\n')
    out.write('Ellipsoids\n')
    out.write('\n')

    for i in range(nnucl):
        out.write("%d %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le %22.15le\n"  \
            % (i+1,1.1739845031423408,1.1739845031423408,1.1739845031423408, \
            quaternions[i][0],quaternions[i][1], quaternions[i][2],quaternions[i][3]))

    out.write('\n')
    out.write('# Bond topology\n')
    out.write('Bonds\n')
    out.write('\n')

    for i in range(nbonds):
        out.write("%d  %d  %d  %d\n" % (i+1,1,bonds[i][0],bonds[i][1]))

    out.close()

    print("## Wrote data to 'data.oxdna'", file=sys.stdout)
    print("## DONE", file=sys.stdout)

# call the above main() function, which executes the program
read_strands (infile)

end_time=timer()
runtime = end_time-start_time
hours = runtime/3600
minutes = (runtime-np.rint(hours)*3600)/60
seconds = (runtime-np.rint(hours)*3600-np.rint(minutes)*60)%60
print( "## Total runtime %ih:%im:%.2fs" % (hours,minutes,seconds), file=sys.stdout)


