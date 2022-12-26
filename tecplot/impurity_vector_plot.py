# Run this file by using 'tec360-env -- python impurity_vector_plot.py'
import tecplot as tp
from tecplot.constant import *

import argparse, os

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file [options]', description = "Generate a vector plot of the given displacement file")
parser.add_argument('file', help = "The file to process")
parser.add_argument('-c', action = 'store_true', help = "Connect to an active TecPlot session")
args = parser.parse_args()

if args.c:
    tp.session.connect()

cwd = os.getcwd().split('/')
imp_type = cwd[6].replace('_effect','')
axis = cwd[9]
T = cwd[-3]
mis = cwd[-4]
d = cwd[-1].replace('_final_','')
if imp_type == "moly":
    imp = cwd[11].replace('at%','') + "Mo"
elif imp_type == "xenon":
    imp = cwd[11].replace('at%','') + "Xe"
else:
    imp = cwd[11].replace('at%','') + cwd[12].replace('at%','')
tp.new_layout()
dataset = tp.data.load_tecplot(args.file,
    read_data_option = ReadDataOption.Replace,
    reset_style = True,
    initial_plot_type = PlotType.Cartesian3D)

tp.data.operate.execute_equation(equation='{XY(K)} = sqrt({X(K)}**2  + {Y(K)}**2)')

tp.active_frame().load_stylesheet('/media/jarinf/Research1/Tecplot/styles/vector_plot_of_impurity_motion.sty')
tp.export.save_png(f"U{imp}_{axis}_{mis}_{T}_{d}_impurity_vector_plot.png")
