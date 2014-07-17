import sys
import os
from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
      name = "ParSNP",
      version = "1.0",
      author = "Todd J. Treangen, Brian Ondov & Adam Phillippy",
      author_email = "treangen+ParSNP@gmail.com",
      description = ("Rapid bacterial core genome alignment and WGST"),
      license = "BSD",
      keywords = "multialignment snps bionformatics",
      url = "http://github.com/marbl/parsnp",
      long_description=read('README'),
      classifiers=[
                   "Development Status :: 1 - Release"
                   ],
      scripts=['src/parsnp.py'],
      #    ext_modules = ext_modules,
      cmdclass = {'build_ext': build_ext}
      )