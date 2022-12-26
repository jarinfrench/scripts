# Run this file by using 'tec360-env -- python vector_plot_rotation.py'
import tecplot as tp
from tecplot.constant import *

import sys, os, glob
import numpy as np

if '-c' in sys.argv:
    tp.session.connect()

root_dir = "/media/jarinf/Research2/working"
save_dir = "/media/jarinf/Research2/working"

files = glob.glob(os.path.join(root_dir, "*_displacement_data_circle.dat"))
# file = os.path.join(root_dir,"adp_100_20degree_T1000_dir_1_displacement_data.dat")
for file in files:
    tp.new_layout()
    dataset = tp.data.load_tecplot(file,
        read_data_option = ReadDataOption.Replace,
        reset_style = True,
        initial_plot_type = PlotType.Cartesian3D
        )

    frame = tp.active_frame()
    plot = frame.plot()
    plot.activate()
    plot.show_shade = False
    plot.show_contour = False
    plot.show_streamtraces = False
    plot.show_edge = False

    plot.axes.x_axis.variable = dataset.variable('Xu')
    plot.axes.y_axis.variable = dataset.variable('Yu')
    plot.axes.z_axis.variable = dataset.variable('Zu')

    plot.value_blanking.active = True
    plot.value_blanking.cell_mode = ValueBlankCellMode.AnyCorner
    constraint = plot.value_blanking.constraint(0)
    constraint.active = True
    constraint.compare_by = ConstraintOp2Mode.UseConstant
    constraint.comparison_operator = RelOp.GreaterThanOrEqual
    constraint.comparison_value = 3.0
    constraint.variable = dataset.variable('Magnitude')

    vect = plot.vector
    vect.u_variable = dataset.variable('X(K)')
    vect.v_variable = dataset.variable('Y(K)')
    vect.w_variable = dataset.variable('Z(K)')
    vect.use_relative = False
    vect.length = 4
    vect.size_arrowhead_by_fraction = False
    vect.arrowhead_size = 1.5
    vect.arrowhead_angle = 10
    vect.use_even_spacing = True
    vect.even_spacing = (6.0, 6.0, 22.0)

    cont = plot.contour(0)
    cont.variable = dataset.variable('Magnitude')
    cont.colormap_name = 'Diverging - Orange/Purple'
    plot.fieldmap(0).vector.color = cont
    cont.levels.reset_levels(np.linspace(0,3,15))
    plot.show_vector = True

    plot.view.theta = 0
    plot.view.psi = 0
    plot.view.fit()

    plot.view.theta = 0
    plot.view.psi = 0
    plot.view.fit()
    plot.axes.orientation_axis.show = False

    tp.export.save_jpeg(f"{save_dir}/{os.path.basename(file).split('_disp')[0]}_circle_vector_plot.jpeg")
