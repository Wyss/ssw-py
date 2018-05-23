import sys
from os.path import dirname, abspath
SSW_PATH = dirname(dirname(dirname(abspath(__file__))))
sys.path.append(SSW_PATH)