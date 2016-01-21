import re
import numpy as np


m = re.search('[A-Za-z]+\s+(\d+)\s+\d+\.\d+\s+(\d+)\s+(\d+\.\d+)\s','H   1 1.089000  2  109.4710  3  120.0000')

if (m):
    print(m.group(0))