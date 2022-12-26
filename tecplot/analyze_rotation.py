# Run this file by using 'tec360-env -- python analyze_rotation.py'
import tecplot as tp
from tecplot.constant import *

import sys, os, glob
from natsort import natsorted

if '-c' in sys.argv:
    tp.session.connect()

root_dir = "/media/jarinf/Research1/U/grain_growth/gamma"
image_save_dir = "/media/jarinf/Research1/U/grain_growth/gamma/tecplot_rotation"

for axis in ['100', '110', '111']:
    for pot in ['adp', 'meam', 'ternary_eam']:
        for angle in ['20degree', '30degree', '45degree', 'sigma7']:
            for T in ['T900', 'T950', 'T1000']: # for T in ['T1050', 'T1100', 'T1150', 'T1200', 'T1250', 'T1300', 'T1350', 'T1400']:
                for dirr in ['dir_1', 'dir_2', 'dir_3', 'dir_4', 'dir_5']:
                    path = '/'.join([root_dir, axis, pot, angle, T, 'large_r', dirr])
                    if os.path.exists(path):
                        if pot == "adp":
                            potential = "ADP"
                        elif pot == "meam":
                            potential = "MEAM"
                        else:
                            potential = "Ternary_EAM"
                        file = natsorted(glob.glob(os.path.join(root_dir, axis, pot, angle, T, "large_r", dirr, "*_tracked.dat")))[-1]
                        tp.new_layout()
                        dataset = tp.data.load_tecplot(file,
                            read_data_option = ReadDataOption.Replace,
                            reset_style = True,
                            initial_plot_type = PlotType.Cartesian3D
                            )

                        frame = tp.active_frame()
                        plot = frame.plot()
                        plot.activate()
                        plot.show_scatter = True
                        plot.show_shade = False

                        # set axes of 3D Cartesian plot to variables V3, V4, V5 in the dataset
                        plot.axes.x_axis.variable = dataset.variable('V3')
                        plot.axes.y_axis.variable = dataset.variable('V4')
                        plot.axes.z_axis.variable = dataset.variable('V5')

                        plot.view.theta = 0
                        plot.view.psi = 0
                        plot.view.fit()

                        scatter = plot.fieldmap(dataset.zone("ZONE 001")).scatter
                        scatter.symbol_type = SymbolType.Geometry
                        scatter.symbol().shape = GeomShape.Sphere
                        scatter.size = 0.5

                        plot.view.theta = 0
                        plot.view.psi = 0
                        plot.view.fit()

                        tp.export.save_png(f"{image_save_dir}/{potential}_{axis}_{angle}_{T}_{dirr.replace('_', '')}_rotation.png")
