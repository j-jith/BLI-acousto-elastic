import sys
import os
import struct
import numpy as np

SIZEOF_DOUBLE = 8


def read(filename):
    """Function to read a deal.ii vector written using the Vector<>::block_write() function.

    Args:
        filename (str): Name of the deal.ii vector file.
    
    Returns:
        numpy array containing the vector.

    """

    # check if file exists
    if not os.path.isfile(filename):
        print('File not found!')
        return
    
    with open(filename, 'rb') as f:
        n = int(f.readline().strip()) # length of the vector

        # skip file till start of vector
        while f.read(1) != b'[':
            pass

        # read vector from file
        v = np.fromfile(f, dtype=float, count=n)

        # check end of file
        eof = f.read(1)

    if eof != b']':
        print('Something has gone wrong! EOF does not match "]".')
        return
    else:
        return v


def read_dof(filename, ndof=1):
    """Function to read a specific degree of freedom(dof) from a deal.ii vector written using the
    Vector<>::block_write() function.

    Args:
        filename (str): Name of the deal.ii vector file.
        ndof (int): Index of dof (starting from 1) [Default is 1].
    
    Returns:
        Value of the dof.

    """

    # check if file exists
    if not os.path.isfile(filename):
        print('File not found!')
        return
    
    with open(filename, 'rb') as f:
        n = int(f.readline().strip()) # length of the vector

        # check if ndof is greater than vector length
        if ndof > n:
            print('ndof greater than vector length!')
            return
        
        # skip file till start of vector
        while f.read(1) != b'[':
            pass

        # skip till desired dof
        f.read(SIZEOF_DOUBLE * (ndof-1))

        # read dof from file
        v, = struct.unpack('d', f.read(SIZEOF_DOUBLE))

        # skip to (end-1) of file
        f.seek(-1, os.SEEK_END)

        # check end of file
        eof = f.read(1)

    if eof != b']':
        print('Something has gone wrong! EOF does not match "]".')
        return
    else:
        return v
