from setuptools import Extension, setup

# Defining the C extension module, specifying the module name and source files
module = Extension("mysymnmf",
                  sources = [
                    'symnmf.c',       # C file containing core implementation
                    'symnmfmodule.c'  # C file for the Python module wrapper
                  ])

# Configuring the setup for the Python package
setup(name ='mysymnmf',
     version ='1.0',
     description ='Python wrapper for custom C extension',
     ext_modules = [module] # List of C extension modules to build
    )