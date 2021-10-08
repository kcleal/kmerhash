from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import os


extras = ["-Wno-sign-compare", "-Wno-unused-function", "-Wno-unused-result", '-Wno-ignored-qualifiers',
          "-Wno-deprecated-declarations"]

ext_modules = list()

root = os.path.abspath(os.path.dirname(__file__))
include_dirs = [numpy.get_include(), root]

print("Include dirs", include_dirs)
print("Extras compiler args", extras)


ext_modules.append(Extension(f"kmerhash.hasher",
                             [f"kmerhash/hasher.pyx"],
                             include_dirs=include_dirs,
                             extra_compile_args=extras,
                             define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                             language="c"))


print("Found packages", find_packages(where="."))
setup(
    name="kmerhash",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    url="https://github.com/kcleal/dysgu",
    description="Generate hash values from a sequence corresponding to kmers in the sequence",
    license="MIT",
    version='0.2',
    python_requires='>=3.7',
    install_requires=[  # runtime requires
            'cython',
        ],
    setup_requires=[
            'cython',
        ],
    packages=["kmerhash", "kmerhash.tests"],
    ext_modules=cythonize(ext_modules),
    include_package_data=True,
    zip_safe=True,

)
