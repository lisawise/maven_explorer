For creating filtered_data

import pandas as pd
import numpy as np
from pyteomics import mgf, mzxml, mass
import os
import matplotlib.pyplot as plt

For aggregation

import pandas as pd
import numpy as np
from pyteomics import mgf, mzxml, mass
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy import cluster
from scipy.spatial.distance import pdist
from scipy.stats import ttest_ind