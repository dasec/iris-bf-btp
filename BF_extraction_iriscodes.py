'''
Implementation of the Bloom filter based biometric template protection for iris images. More details in:

[IET14] C. Rathgeb, F. Breitinger, C. Busch, H. Baier, "On application of bloom filters to iris biometrics",
        in IET Biometrics, vol. 3, no. 4, pp. 207 - 218, 2014

[IF18] M. Gomez-Barrero, C. Rathgeb, G. Li, R. Raghavendra, J. Galbally and C. Busch
        "Multi-Biometric Template Protection Based on Bloom Filters", in Information Fusion, vol. 42, pp. 37-50, 2018.

Please remember to reference articles [IET14] and [IF18] on any work made public, whatever the form,
based directly or indirectly on these metrics.
'''

__author__ = "Marta Gomez-Barrero"
__copyright__ = "Copyright (C) 2017 Hochschule Darmstadt"
__license__ = "License Agreement provided by Hochschule Darmstadt (https://share.nbl.nislab.no/g03-06-btp/iris-bf-btp/blob/master/hda-license.pdf)"
__version__ = "1.0"

import numpy
import math
import os
import argparse
from random import sample

######################################################################
### Parameter and arguments definition

parser = argparse.ArgumentParser(description='Extract unprotected LGBPHS and protected Bloom filter templates from the FERET DB.')

# location of source images, final templates and intermediate steps (the latter for debugging purposes)
parser.add_argument('DBdir', help='directory where the compressed face DB is stored', type=str)
parser.add_argument('--DB_BFtemplates', help='directory where the unprotected face templates will be stored', type=str, nargs='?', default = './BFTemplates/')
parser.add_argument('--nXORKeys', help='number of keys for the XOR operation for the feature level fusion. for a unimodal system it should be the default: 0', type=int, nargs='?', default = 0)

args = parser.parse_args()
DBdir = args.DBdir
DB_BFtemplates = args.DB_BFtemplates
nXORKeys = args.nXORKeys

if not os.path.exists(DBdir):
    os.mkdir(DBdir)
if not os.path.exists(DB_BFtemplates):
    os.mkdir(DB_BFtemplates)

# Parameters of Bloom filter extraction
N_BITS_BF = 10
N_WORDS_BF = 32
N_BF_X = 512//N_WORDS_BF
N_BF_Y = 20//N_BITS_BF
N_BLOCKS = N_BF_X * N_BF_Y
BF_SIZE = int(math.pow(2, N_BITS_BF))


####################################################################
### Some auxiliary functions

def add_unlinkability(features, keyPERM):
    '''Permutes rows within regions of an iris-code to achieve unlinkability'''
    perm_feat = numpy.zeros(shape=feat.shape, dtype=int)

    # divide iris-code in four regions, and reshape each region to a size [N_BITS_BF * N_BF_X / 2, N_WORDS_BF]
    featsAux1 = numpy.reshape(features[0 : N_BITS_BF, 0 : (N_BF_X / 2) * N_WORDS_BF], [N_BITS_BF * N_BF_X / 2, N_WORDS_BF])
    featsAux2 = numpy.reshape(features[0 : N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF : N_BF_X * N_WORDS_BF], [N_BITS_BF * N_BF_X / 2, N_WORDS_BF])
    featsAux3 = numpy.reshape(features[N_BITS_BF : 2*N_BITS_BF, 0 : (N_BF_X / 2) * N_WORDS_BF], [N_BITS_BF * N_BF_X / 2, N_WORDS_BF])
    featsAux4 = numpy.reshape(features[N_BITS_BF : 2*N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF : N_BF_X * N_WORDS_BF], [N_BITS_BF * N_BF_X / 2, N_WORDS_BF])

    # permute rows within each region
    perm_feat[0: N_BITS_BF, 0: (N_BF_X / 2) * N_WORDS_BF] = numpy.reshape(featsAux1[keyPERM[0, :], :], [N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF])
    perm_feat[0: N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF: N_BF_X * N_WORDS_BF] = numpy.reshape(featsAux2[keyPERM[2, :], :], [N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF])
    perm_feat[N_BITS_BF: 2 * N_BITS_BF, 0: (N_BF_X / 2) * N_WORDS_BF] = numpy.reshape(featsAux3[keyPERM[2, :], :], [N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF])
    perm_feat[N_BITS_BF: 2 * N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF: N_BF_X * N_WORDS_BF] = numpy.reshape(featsAux4[keyPERM[3, :], :], [N_BITS_BF, (N_BF_X / 2) * N_WORDS_BF])

    return perm_feat


def extract_BFs_from_iriscode(feat):
    '''Extracts the BF template from the iris-code'''
    template = numpy.zeros(shape=[N_BLOCKS, BF_SIZE], dtype=int)

    index = 0
    for x in range(N_BF_X):
        for y in range(N_BF_Y):
            bf = numpy.zeros(shape=[BF_SIZE])

            ini_x = x * N_WORDS_BF
            fin_x = (x + 1) * N_WORDS_BF
            ini_y = y * N_BITS_BF
            fin_y = (y + 1) * N_BITS_BF
            new_hists = feat[ini_y: fin_y, ini_x: fin_x]

            for w in range(N_WORDS_BF):
                hist = new_hists[:, w]
                location = int('0b' + ''.join([str(a) for a in hist]), 2)
                bf[location] = int(1)

            template[index, :] = bf
            index += 1

    return template

def extract_BFs_from_iriscode_XOR(feat, keys):
    '''Extracts the BF template from the iris-code, using multiple XOR keys (for feature-level fusion)'''
    template = numpy.zeros(shape=[N_BLOCKS, BF_SIZE], dtype=int)

    index = 0
    for x in range(N_BF_X):
        for y in range(N_BF_Y):
            bf = numpy.zeros(shape=[BF_SIZE])

            ini_x = x * N_WORDS_BF
            fin_x = (x + 1) * N_WORDS_BF
            ini_y = y * N_BITS_BF
            fin_y = (y + 1) * N_BITS_BF
            new_hists = feat[ini_y: fin_y, ini_x: fin_x]

            for w in range(N_WORDS_BF):
                hist = new_hists[:, w]
                location = int('0b' + ''.join([str(a) for a in hist]), 2)
                for k in keys:
                    bf[location ^ k] = int(1)

            template[index, :] = bf
            index += 1

    return template

####################################################################
### Template extraction

# Define permutation key to provide unlinkability
key = numpy.zeros(shape = [4, N_BITS_BF * N_BF_X / 2], dtype = int)
for j in range(4):
    key[j, :] = numpy.random.permutation(N_BITS_BF * N_BF_X / 2)

print("Extracting BF templates for the DB...")

if not(nXORKeys == 0):
    xorKeys = sample(range(BF_SIZE), nXORKeys)
    for filename in os.listdir(DBdir):
        print(filename[0:-4])

        # load iriscode
        f = open(DBdir + '/' + filename, 'r')
        feat = numpy.asarray([list(map(int, list(line.rstrip()))) for line in f.readlines()])

        # permute features to provide unlinkability
        featUnlink = add_unlinkability(feat, key)

        # extract BFs
        bfs = extract_BFs_from_iriscode_XOR(featUnlink, xorKeys)
        numpy.savetxt(DB_BFtemplates + filename[0:-4] + '_BFtemplate.txt', bfs, fmt='%d')

else:
    for filename in os.listdir(DBdir):
        print(filename[0:-4])

        # load iriscode
        f = open(DBdir + '/' + filename, 'r')
        feat = numpy.asarray([list(map(int, list(line.rstrip()))) for line in f.readlines()])

        # permute features to provide unlinkability
        featUnlink = add_unlinkability(feat, key)

        # extract BFs
        bfs = extract_BFs_from_iriscode(featUnlink)
        numpy.savetxt(DB_BFtemplates + filename[0:-4] + '_BFtemplate.txt', bfs, fmt='%d')
