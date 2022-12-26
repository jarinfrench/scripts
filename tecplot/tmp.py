import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *

# Uncomment the following line to connect to a running instance of Tecplot 360:
# tp.session.connect()

tp.data.operate.execute_equation(equation='{XY(K)} = sqrt({X(K)}**2  + {Y(K)}**2)')
# End Macro.

