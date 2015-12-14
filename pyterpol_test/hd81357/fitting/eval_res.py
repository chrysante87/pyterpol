import numpy as np
import pyterpol

itf = pyterpol.Interface()
itf.load('hd81357.sav')
itf.populate_comparisons()
print itf.get_degrees_of_freedom()
print itf.compute_chi2_treshold()





