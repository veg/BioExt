#!/usr/bin/env python

from __future__ import division, print_function

import os.path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


sourcefiles = [os.path.join("BioExt", "tn93", "_tn93.pyx"),
               os.path.join("BioExt", "tn93", "tn93.c")]


class BioCythonNumpyBuildExt(build_ext):
    """
    This waits until the ``setup_requires`` packages are installed, and only then does the necessary configuration
    changes. Cf. https://stackoverflow.com/a/21621689.
    is interested in that.
    """
    def finalize_options(self):
        self.add_extensions()
        super().finalize_options()
        self.include_numpy()
        self.add_package_data()

    def include_numpy(self):
        # include numpy dirs, this needs numpy installed
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

    def add_package_data(self):
        # add package data, this needs biopython installed
        __builtins__.__NUMPY_SETUP__ = False
        from BioExt.references._factory import _installrefdirs
        self.package_data = {
            'BioExt': [
                'data/fonts/ttf/*.ttf',
                'data/scorematrices/*.txt'
            ] + _installrefdirs
        }

    def add_extensions(self):
        # add tn93 extension, this needs cython installed
        from Cython.Build import cythonize

        self.distribution.ext_modules.extend(cythonize([
            Extension(
                "BioExt.tn93._tn93",
                include_dirs=[os.path.join("BioExt", "tn93")],
                sources=sourcefiles
            )
        ]))


ext_modules = [
    Extension(
        'BioExt.align._align',
        sources=[
            os.path.join('BioExt', 'align', '_align.c'),
            os.path.join('BioExt', 'align', 'alignment.c')
            ],
        libraries=['m'],
        extra_compile_args=['-O3', '-I.']
        ),
    Extension(
        'BioExt.merge._merge',
        sources=[
            os.path.join('BioExt', 'merge', '_merge.c'),
            os.path.join('BioExt', 'merge', 'merge.cpp')
            ],
        extra_compile_args=['-O3', '-I.']
        ),
    Extension(
        'BioExt.rateclass._rateclass',
        sources=[
            os.path.join('BioExt', 'rateclass', '_rateclass.cpp'),
            os.path.join('BioExt', 'rateclass', 'rateclass.cpp')
            ],
        extra_compile_args=['-O3', '-I.']
        )
    ]


setup(
    name='bioext',
    version='0.20.6',
    description='Misc utilities and definitions not included or hidden in BioPython',
    author='N Lance Hepler',
    author_email='nlhepler@gmail.com',
    url='http://github.com/veg/bioext',
    license='GNU GPL version 3',
    packages=[
        'BioExt',
        'BioExt.align',
        'BioExt.args',
        'BioExt.collections',
        'BioExt.errorize',
        'BioExt.freetype',
        'BioExt.freetype.ft_enums',
        'BioExt.io',
        'BioExt.io.BamIO',
        'BioExt.io.LazyAlignIO',
        'BioExt.io.SamIO',
        'BioExt.merge',
        'BioExt.misc',
        'BioExt.ndarray',
        'BioExt.optimize',
        'BioExt.orflist',
        'BioExt.phylo',
        'BioExt.quiver',
        'BioExt.rateclass',
        'BioExt.references',
        'BioExt.scorematrices',
        'BioExt.stats',
        'BioExt.tn93',
        'BioExt.uds',
        'BioExt.untranslate'
        ],
    package_dir={
        'BioExt': 'BioExt',
        'BioExt.align': 'BioExt/align',
        'BioExt.args': 'BioExt/args',
        'BioExt.collections': 'BioExt/collections',
        'BioExt.errorize': 'BioExt/errorize',
        'BioExt.freetype': 'BioExt/freetype',
        'BioExt.freetype.ft_enums': 'BioExt/freetype/ft_enums',
        'BioExt.io': 'BioExt/io',
        'BioExt.io.BamIO': 'BioExt/io/BamIO',
        'BioExt.io.LazyAlignIO': 'BioExt/io/LazyAlignIO',
        'BioExt.io.SamIO': 'BioExt/io/SamIO',
        'BioExt.merge': 'BioExt/merge',
        'BioExt.misc': 'BioExt/misc',
        'BioExt.ndarray': 'BioExt/ndarray',
        'BioExt.optimize': 'BioExt/optimize',
        'BioExt.orflist': 'BioExt/orflist',
        'BioExt.phylo': 'BioExt/phylo',
        'BioExt.quiver': 'BioExt/quiver',
        'BioExt.rateclass': 'BioExt/rateclass',
        'BioExt.references': 'BioExt/references',
        'BioExt.scorematrices': 'BioExt/scorematrices',
        'BioExt.stats': 'BioExt/stats',
        'BioExt.tn93': 'BioExt/tn93',
        'BioExt.uds': 'BioExt/uds',
        'BioExt.untranslate': 'BioExt/untranslate'
        },
    scripts=[
        'scripts/bam2fna',
        'scripts/bam2msa',
        'scripts/bamclip',
        'scripts/bealign',
        'scripts/clipedge',
        'scripts/consensus',
        'scripts/msa2bam',
        'scripts/seqmerge',
        'scripts/translate'
        # 'scripts/variants'
        ],
    ext_modules=ext_modules,
    cmdclass={'build_ext': BioCythonNumpyBuildExt},
    setup_requires=['biopython >=1.78', 'cython', 'numpy >= 1.22, < 1.23'],
    install_requires=[
        'biopython >=1.78',
        'numpy >= 1.22, < 1.23',
        'scipy >=0.15',
        'pysam >=0.17',
        'joblib',
        'six'
        ]
    )
