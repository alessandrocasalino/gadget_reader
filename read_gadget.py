# Alessandro Casalino - University of Bologna
#
# Script to import header and particles position/velocity from a Gadget format = 1 or format = 2 file
#
# In the script we consider the data structure of Gadget2 format = 1 and format = 2 snapshots:
# 1. DATA STRUCTURE while using the pointer
#    blocksize (4 bytes - int32) --- data (data_type * number of data) --- blocksize (4 bytes - int32)
#    Therefore when we want to go to next block with the pointer, we also need to read the final blocksize
# 2. format = 2 snapshots also have header block as the header of the main data block
#    This header is 4 chars and its size is 8 bytes

import os
import numpy as np

# ------------ SETTINGS ------------

# Full path to the snapshot (without any .x suffix if snapshots splitted in more files)
filename = "/scratch/extra/marco.baldi5/Alessandro/snapdir_005/snap_005"
#filename = "/scratch/extra/marco.baldi5/Alessandro/C-Gadget_examples/output_lcdm_128/snapdir_005/snap_005"


# ------------ FILE IMPORT ------------
# Open the file if available
def open_snapshot(filename, nfile = 0):

    import os

    # Check the existance of the file
    if os.path.exists(filename):
          curfilename = filename
    # If the snapshot is divided in several files, I save the full path with ".nfile" instead
    elif os.path.exists(filename+".0"):
        if os.path.exists(filename+"."+str(nfile)):
            curfilename = filename+"."+str(nfile)
        else:
            print("Error: File not found:", filename+"."+str(nfile))
    # If nothing is found, throw an error
    else:
      print("Error: File not found:", filename)
      return None

    # Open the file in read mode
    snap = open(curfilename, 'rb')

    return snap



# ------------ FORMAT CHECK ------------
# Check if the file has the correct initial blocksize, the format and the endian
def check_blocksize(snap):

    # Check the block size
    blocksize = np.fromfile(snap,dtype=np.int32,count=1)

    # We must check if the format is small (most significant figure stored at small address) or big-endian
    # Firstly we check the size of the blocksize
    if blocksize[0] == 256:
        gformat = 1
        swap = 0
    elif blocksize[0] == 8:
        gformat = 2
        swap = 0
    # But if we can't find anything we try to switch to the other endian, swapping the bytes, with byteswap()
    else:
        blocksize.byteswap(True)
        if blocksize[0] == 256:
            gformat = 1
            swap = 1
        elif blocksize[0] == 8:
            gformat = 2
            swap = 1
        else:
            print("Error: The format of the selected snapshot seems not to be Gadget-1/2 compatible: ", filename)
            return -1, -1

    # If the format is 2, we have to skip the header of the block ("POS ", "VEL ", etc..) and the final blocksize
    if gformat == 2:
        snap.seek(16, os.SEEK_CUR)

    return gformat, swap



# ------------ HEADER ------------
# Read the header informations
# For informations about the structure check the Gadget-2 user guide
# https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf (page 33, table 4)
# - snap               : snapshot file
# - print_informations : print header informations
def read_header(filename, print_informations = True):

    # Open snapshot
    snap = open_snapshot(filename)
    if not snap:
        return None

    # Check if the blocksize of the header is 256, and the endian
    gformat, swap = check_blocksize(snap)
    if swap == -1:
        return None

    # Number of particles (array by type) IN THIS FILE
    # 0 - Gas particles
    # 1 - Halo particles
    # 2 - Disk particles
    # 3 - Bulge particles
    # 4 - Star particles
    # 5 - Bndry particles
    npart_file = np.fromfile(snap,dtype=np.int32,count=6)

    # Mass of the particles (Gadget units, default = 10e10 Msun/h)
    masspart = np.fromfile(snap,dtype=np.float64,count=6)

    # Time
    # Take first element of the np.array with [0]
    time = np.fromfile(snap,dtype=np.float64,count=1)[0]

    # Redshift
    redshift = np.fromfile(snap,dtype=np.float64,count=1)[0]

    # Star Formation flag
    flag_sfr = np.fromfile(snap,dtype=np.int32,count=1)[0]

    # Feedback flag
    flag_feedback = np.fromfile(snap,dtype=np.int32,count=1)[0]

    # Total number of particles IN THE SIMULATION
    # This differs from npart_file if there are more files for the snapshot
    npart_sim = np.fromfile(snap,dtype=np.int32,count=6)

    # Cooling flag
    flag_cooling = np.fromfile(snap,dtype=np.int32,count=1)[0]

    # Number of files in each snapshot
    nfile = np.fromfile(snap,dtype=np.int32,count=1)[0]

    # Boxsize (in Gadget units, default = 1e-3 MPc/h)
    boxsize = np.fromfile(snap,dtype=np.float64,count=1)[0]

    # Omega_cdm
    Omega_cdm = np.fromfile(snap,dtype=np.float64,count=1)[0]

    # Omega_Lambda
    Omega_Lambda = np.fromfile(snap,dtype=np.float64,count=1)[0]

    # Hubble parameter h
    h = np.fromfile(snap,dtype=np.float64,count=1)[0]

    # Swap for correct small/big-endian
    if swap == 1:
        npart_file.byteswap(True)
        masspart.byteswap(True)
        time = time.byteswap(True)
        redshift = redshift.byteswap(True)
        flag_sfr = flag_sfr.byteswap(True)
        flag_feedback = flag_feedbackback.byteswap(True)
        npart_sim = npart_sim.byteswap(True)
        flag_cooling = flag_cooling.byteswap(True)
        nfile.byteswap(True)
        boxsize = boxsize.byteswap(True)
        Omega_cdm = Omega_cdm.byteswap(True)
        Omega_Lambda = Omega_Lambda.byteswap(True)
        h = h.byteswap(True)

    snap.close()

    Nparticles_file = np.sum(npart_file)
    Nparticles = np.sum(npart_sim)

    # Print informations about the simulation
    if print_informations:
        print("Simulation informations")
        print("- Snapshot format:     ", gformat)
        print("- Number of particles: ", Nparticles)
        print("- Redshift:            ", redshift)
        print("- Boxsize:             ", boxsize)
        print("- h:                   ", h)
        print("- Omega_cdm:           ", Omega_cdm)
        print("- Omega_Lambda:        ", Omega_Lambda)

    return nfile, Nparticles, Nparticles_file, gformat, swap



# ------------ DATA ------------
# Read the block related to position, velocity, etc..
# - filename : snapshot file name
# - block    : name of the block to extract, i.e. POS or VEL
# - debug    : make the function more talkative
# NOTE: the particles are stored in order of type, so in principle to extract only one type you might need
#       n_part_file on each snapshot file
def extract_block(filename, block, verbose = 1):

    nfile, Nparticles, Nparticles_file, gformat, swap = read_header(filename, print_informations = False)

    if verbose > 0: print(" " + "-" * 45)

    # Block to extract:  1 for POS, 2 for VEL
    # We also set the data type of the data to extract
    # For the whole list of available blocks, check:
    # https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf (page 32)
    if block == "POS ":
        block_wanted = 1
        data_type = np.dtype((np.float32,3))
    elif block == "VEL ":
        block_wanted = 2
        data_type = np.dtype((np.float32,3))
    elif block == "MASS":
        block_wanted = 4
        data_type = np.float32
    elif block == "ID  ":
        block_wanted = 3
        data_type = np.int32
    else:
        print("Error: No Gadget compatible block selected!")
        return None

    data = np.zeros(0, dtype = data_type)

    # We loop the different files related to the same snapshot
    for i in range(nfile):

        snap = open_snapshot(filename, i)
        if not snap: return # Interrupt if file not found

        # Find size of the file
        snap.seek(0, os.SEEK_END)   # Go to the end
        filesize = snap.tell()      # Current file position
        snap.seek(0, os.SEEK_SET)   # Go to the beginning

        block_num = 0
        found = False

        # We loop until we reach the end of the file
        while(snap.tell()<filesize):

            # Gadget format = 1 condition : simply select the block number
            if gformat==1:
                found = gformat == 1 and block_num == block_wanted
            # Gadget format = 2 condition : check the header of the block
            else:
                # Advance the pointer through the blocksize (int32, 4 bytes)
                snap.seek(4,os.SEEK_CUR)
                # Read and check if the block is the one we are looking for
                # Note: reading advances the pointer by the specified amount of bytes
                found = snap.read(4).decode() == block
                # After we read, we have to advance the pointer of the header by:
                # remaining data size (4 bytes) + blocksize size (4 bytes)
                snap.seek(8,os.SEEK_CUR)

            # Check the blocksize (int32, 4 bytes)
            blocksize = np.fromfile(snap, dtype = np.int32, count = 1)[0]

            # If we found the block, read it, otherwise keep searching (else)
            if found:

                if verbose > 0: print("Find block '" + block + "' in snapshot", i ,". Reading...")

                # Fix "ID  " if long int is used
                if i == 0 and block == "ID  " and blocksize == Nparticles_file * np.dtype(data_type).itemsize * 2:
                    data_type = np.int64

                # We append the whole block of data (we split them in 3d arrays later if POS and VEL)
                # Note: This should be faster than splitting data here with a loop
                data = np.append(data, np.fromfile(snap, dtype = data_type, count = int(blocksize/np.dtype(data_type).itemsize)))

                if verbose > 1: print("Reading complete.")

                break

            else:

                # If this is not the data block we have to advance the pointer to the blocksize
                snap.seek(blocksize, os.SEEK_CUR)
                # But we also have to advance it for the final blocksize int
                # Therefore this check (or analogous pointer advance) is compulsory
                check = np.fromfile(snap, dtype = np.int32, count = 1)[0]
                # The blocksize found before should be the same as the final one. If not, something went wrong!
                if check != blocksize:
                    print("\nError: The initial and final blocksize are different.")

            block_num = block_num + 1

        if not found:
            print("\nWarning: Can not find block '" + block + "'.")
            break

        snap.close()

    # Check for little/big endian
    if swap == 1:
        data.byteswap(True)

    # Finally we split the data to have 3d arrays if POS or VEL
    if block in ["POS ", "VEL "]:
        if len(data) % 3 == 0:
            data=data.reshape(int(len(data)/3),3)
        else:
            print("\nError: Data is not divisible by 3, can not split data")
            return None

    # Check if the length of the fata is equal to the number of particles
    if block in ["POS ", "VEL ", "ID  "] and len(data) != Nparticles:
        print("\nError: Data seems not correct. Shape:", len(data), " - Number particles:", Nparticles)

    return data



# ------------ MAIN ------------

nfile, Nparticles, Nparticles_file, gformat, swap = read_header(filename)

if gformat in [1,2]:

    # Write additional messages during snapshot read
    verbose = 0

    # The block name should be for chars, add spaces to complete it, e.g. "POS " and not "POS"
    pos = extract_block(filename, "POS ", verbose = verbose)
    vel = extract_block(filename, "VEL ", verbose = verbose)
    mass = extract_block(filename, "MASS", verbose = verbose)
    pid = extract_block(filename, "ID  ", verbose = verbose)

    print(" " + "-" * 45)

    particle = int(100)
    print("Example: ")
    print("   - Particle", particle, "position :", pos[particle])
    print("   - Particle", particle, "velocity :", vel[particle])
    #print("   - Particle", particle, "mass:", mass[particle])
    print("   - Particle", particle, "ID       :", pid[particle])
