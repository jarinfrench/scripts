import numpy as np
import proplot as plt
plt.rc['figure.facecolor'] = 'white'
plt.rc['fontsize'] = 14
plt.rc['grid.color'] = 'white'
plt.rc.cycle = 'Dark2'

data = [np.loadtxt(f,
        dtype={'names': ('time', 'grain1', 'grain2'),
               'formats': ('float', 'float', 'float')},
        skiprows = 1) for f in ["/media/jarinf/Research2/tmp/Cooper_same_initial_Basak/home/jarinf/projects/uo2/grain_growth/cylindrical/111/45degree/Cooper/T2900/large_r/all_area_data_avg.txt",
                                "/media/jarinf/Research2/tmp/home/jarinf/projects/uo2/grain_growth/cylindrical/111/45degree/Cooper/T2900/large_r/all_area_data_avg.txt"]]

fig, ax = plt.subplots(figwidth = 4, figheight = 2.5, xlabel = 'Time (ps)', ylabel = 'Area ($\AA^2$)')
ax.scatter(data[0]['time'], data[0]['grain2'], s = 10, label = "Simple Anneal")
ax.scatter(data[1]['time'], data[1]['grain2'], s = 10, label = "Complex Anneal")
ax.legend(loc = 'best', ncols = 1, frame = False)
plt.show()
