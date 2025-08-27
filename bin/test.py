#!/usr/bin/env python3
err = False

try:
    import Bio
except:
    print('You need to install biopython >= 1.79')
    err = True

try:
    import pandas
except:
    print('You need to install pandas >= 2.2.2')
    err = True

try:
    import tensorflow
except:
    print('You need to install tensorflow >= 2.10.0')
    err = True

try:
    import RNA
except:
    print('You need to install viennarna >= 2.5.1 with support for your python version')
    err = True

if err:
    print('requirements not met')
    exit(1)
else:
    exit()
