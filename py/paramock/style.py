import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['lines.linewidth'] = 3

mpl.rcParams['font.size']   = 14
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif']  = 'mathpazo'

mpl.rcParams['text.usetex']    = True
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.grid']      = True

mpl.rcParams['grid.linestyle'] = 'dotted'

mpl.rcParams['xtick.labelsize']     = 18
mpl.rcParams['xtick.top']           = True
mpl.rcParams['xtick.minor.visible'] = True

mpl.rcParams['ytick.labelsize']     = 18
mpl.rcParams['ytick.right']         = True
mpl.rcParams['ytick.minor.visible'] = True

mpl.rcParams['savefig.bbox'] = 'tight'

mpl.rcParams['figure.figsize'] = (12,8)
