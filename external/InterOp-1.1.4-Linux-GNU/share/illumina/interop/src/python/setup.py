from setuptools import setup, Distribution
import os

class BinaryDistribution(Distribution):

    def has_ext_modules(self):
        return True

    def is_pure(self):
        return False

setup(
    name='interop',
    version='1.1.4',
    description="""The Illumina InterOp libraries are a set of common routines used for reading InterOp metric files
                   produced by Illumina sequencers. These libraries are backwards compatible and capable of supporting
                   prior releases of the software, with one exception: GA systems have been excluded.""",
    maintainer='Illumina Inc.',
    url='https://github.com/Illumina/interop',
    license='GPL',
    download_url='https://github.com/Illumina/interop/releases/latest',
    packages=['interop'],
    include_package_data=True,
    package_data={
        'interop': [ '*.so', '*.pyd', '*.dylib' ],
    },
    classifiers=[
    'Development Status :: 5 - Production/Stable',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords="Illumina sequencer HiSeqX HiSeq NextSeq MiniSeq NovaSeq MiSeq SBS genome",
    install_requires=['numpy>=1.13'],
    distclass=BinaryDistribution
)

