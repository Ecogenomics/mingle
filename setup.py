from distutils.core import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'mingle', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='mingle',
    version=version(),
    author='Joel Boyd, Ben Woodcroft, Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['mingle'],
    scripts=['bin/mingle'],
    package_data={'mingle': ['VERSION']},
    url='http://pypi.python.org/pypi/mingle/',
    license='GPL3',
    description='Infer gene trees with HMMs or BLAST homology search.',
    long_description=open('README.md').read(),
    install_requires=[
        "biopython >= 1.63",
        "graftm >= 0.0.1"],
)
