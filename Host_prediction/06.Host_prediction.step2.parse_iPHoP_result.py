#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call    
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
# Aim: Parse iPHoP result to get the final host prediction result
# The following rules were used:
# (1) If "blast" or "CRISPR" results were obtained for one virus genome, 
#     the result with the highest confidence score was assigned as the final result. 
# (2) If only "iPHoP-RF" results were obtained for one virus genome, 
#     the result with the highest confidence score was assigned as the final result.

 

