import distutils.core
import Cython.Build

distutils.core.setup(ext_modules = Cython.Build.cythonize("convertData.pyx", language_level = 3))
