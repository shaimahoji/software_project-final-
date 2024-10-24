from setuptools import Extension, setup
import numpy

module = Extension("mysymnmf",
                  sources=[
                    'symnmf.c',
                    'symnmfmodule.c'
                  ],
                  include_dirs=[numpy.get_include()] 
                  )
setup(name='mysymnmf',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])