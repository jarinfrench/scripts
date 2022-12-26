# Run this file by using 'tec360-env -- python rotation_gif_maker.py'
import tecplot as tp
from tecplot.constant import *

import sys, os, glob
from natsort import natsorted

if '-c' in sys.argv:
    tp.session.connect()

root_dir = "/media/jarinf/Research3/Research1_backup/U/grain_growth/gamma"
save_dir = "/media/jarinf/Research2/Pictures/U_rotation_images"

# axis='100'
# pot='ternary_eam'
# angle='20degree'
# T='T1100'
# dirr='dir_3'
for axis in ['100', '110', '111']:
    for pot in ['adp', 'meam', 'ternary_eam']:
        for angle in ['20degree', '30degree', '45degree', 'sigma7']:
            for T in ['T900', 'T950', 'T1000', 'T1050', 'T1100', 'T1150', 'T1200', 'T1250', 'T1300', 'T1350', 'T1400']:
                for dirr in ['dir_1', 'dir_2', 'dir_3', 'dir_4', 'dir_5', 'dir_6', 'dir_7', 'dir_8', 'dir_9', 'dir_10']:
                    path = "/".join([root_dir, axis, pot, angle, T, 'large_r', dirr])
                    if os.path.exists(path):
                        if pot == "adp":
                            potential = "ADP"
                        elif pot == "meam":
                            potential = "MEAM"
                        else:
                            potential = "Ternary_EAM"
                        files = natsorted(glob.glob(os.path.join(path, "*_tracked.dat")))
                        if not len(files) == 3:
                            print(f"WARNING: {path} does not have the correct number of tracked files ({len(files)})")
                            continue
                        tp.new_layout()
                        dataset = tp.data.load_tecplot(files,
                            read_data_option = ReadDataOption.Replace,
                            reset_style = True,
                            initial_plot_type = PlotType.Cartesian3D,
                            add_zones_to_existing_strands = True
                            )

                        frame = tp.active_frame()
                        for time, z in enumerate(frame.dataset.zones()):
                            z.strand = 1
                            z.solution_time = time

                        plot = frame.plot()
                        plot.activate()
                        plot.active_fieldmaps += [0]
                        plot.show_scatter = True
                        plot.show_shade = False

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

                        plot.solution_time = 2.0

                        plot.view.theta = 0
                        plot.view.psi = 0
                        plot.view.fit()

                        tp.export.save_time_animation_jpeg(
                        f"{save_dir}/{potential}_{axis}_{angle}_{T}_{dirr.replace('_','')}_rotation.jpeg"
                        )
