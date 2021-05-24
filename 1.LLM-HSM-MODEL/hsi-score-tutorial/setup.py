from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(["hsiscore.pyx", "LLM_model_cpy.pyx"])
)
