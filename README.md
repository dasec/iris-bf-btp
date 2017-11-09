# BTP based on Bloom filters for iris-codes

Biometric Template Protection based on Bloom filters for iris-codes. Bloom filter template extraction based on the method proposed in [[IET14]](http://ieeexplore.ieee.org/abstract/document/6985870/) and [[IF18]](http://www.sciencedirect.com/science/article/pii/S1566253516301233).

## License
This work is licensed under license agreement provided by Hochschule Darmstadt ([h_da-License](/hda-license.pdf)).

## Instructions

### Dependencies
* numpy
* math
* argparser
* os
* random.sample

### Usage

1. Run BF_extraction_iriscodes.py to extract the protected templates from iris-codes provided.

    ```python
	usage: BF_extraction_iriscodes.py [-h] [--DB_BFtemplates [DB_BFTEMPLATES]]
                                      DBdir
    
    Extract unprotected LGBPHS and protected Bloom filter templates from the FERET
    DB.
    
    positional arguments:
      DBdir                 directory where the compressed face DB is stored

    optional arguments:
      -h, --help            show this help message and exit
      --DB_BFtemplates [DB_BFTEMPLATES]
                            directory where the unprotected face templates will be
                            stored
      --nXORKeys [NXORKEYS]
                            number of keys for the XOR operation for the feature
                            level fusion. for a unimodal system it should be the
                            default: 0
    ```
	1. Input: folder containing the iris-codes to be protected with the Bloom filter based transformation. Each iris-code should be stored as a text file, with one row per line.
	2. Output: extracted templates, stored as text files in DB_BFtemplates.
	3. Other parameters for the Bloom filter template extraction might be changed at the top of the script. The values used in [IET14] and [IF18] are included as default.
2. Run computeScores.py to compute the mated and non-mated scores

    ```python
    usage: computeScores.py [-h] [--scoresDir [SCORESDIR]]
                            [--matedScoresFile [MATEDSCORESFILE]]
                            [--nonMatedScoresFile [NONMATEDSCORESFILE]]
                            DB_BFtemplates matedComparisonsFile
                            nonMatedComparisonsFile
    
    Compute protected Bloom filter scores from a given DB and protocol.
    
    positional arguments:
      DB_BFtemplates        directory where the protected BF templates are stored
      matedComparisonsFile  file comprising the mated comparisons to be carried
                            out
      nonMatedComparisonsFile
                            file comprising the non-mated comparisons to be
                            carried out
    
    optional arguments:
      -h, --help            show this help message and exit
      --scoresDir [SCORESDIR]
                            directory where unprotected and protected scores will
                            be stored
      --matedScoresFile [MATEDSCORESFILE]
                            file comprising the mated scores computed
      --nonMatedScoresFile [NONMATEDSCORESFILE]
                            file comprising the non-mated scores computed
    ```
	1. Input: folder with the protected templates and the files comprising the mated and non-mated comparisons. Such files should comprise one comparison per line, with the two files to be compared separated by a blank space. The directory should *not* be included in the filename.
	2. Output: mated and non-mated scores, stored in text files with one score per row.

## References

More details in:

- [[IET14]](http://ieeexplore.ieee.org/abstract/document/6985870/) C. Rathgeb, F. Breitinger, C. Busch, H. Baier, "On application of bloom filters to iris biometrics", in IET Biometrics, vol. 3, no. 4, pp. 207 - 218, 2014

- [[IF18]](http://www.sciencedirect.com/science/article/pii/S1566253516301233) M. Gomez-Barrero, C. Rathgeb, G. Li, R. Raghavendra, J. Galbally and C. Busch, "Multi-Biometric Template Protection 
Based on Bloom Filters", in Information Fusion, vol. 42, pp. 37-50, 2018.

Please remember to reference articles [IET14] and [IF18] on any work made public, whatever the form,
based directly or indirectly on these scripts.