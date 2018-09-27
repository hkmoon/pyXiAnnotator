#!/usr/bin/env python3
from setuptools import setup
import os

version_path = os.path.join(
    os.path.dirname(__file__),
    'pyxiannotator',
    'version.txt'
)
with open(version_path, 'r') as version_file:
    pyXiAnnotator_version = version_file.read().strip()

setup(
    name             = 'pyxiannotator',
    version          = pyXiAnnotator_version,
    packages         = ['pyxiannotator'],
    package_dir      = {'pyxiannotator': 'pyxiannotator'},
    package_data     = {
        'pyxiannotator': [
            'version.txt',
            'jar/*.jar'
        ]
    },
    setup_requires   = ['Cython'],
    install_requires = ['jnius', 'requests', 'urllib3'],
    description      = 'MS spectra annotation and analysis',
    long_description = 'pyXiAnnotator - mass spectrometry spectra annotation and analysis',
    author           = 'L. Kolbowski',
    author_email     = 'lars.kolbowski@campus.tu-berlin.de',
    url              = 'https://github.com/Rappsilber-Laboratory/pyXiAnnotator',
    license          = 'Apache License 2.0',
    classifiers      = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: Unix',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: POSIX :: SunOS/Solaris',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Mass Spectrometry',
    ]
)
