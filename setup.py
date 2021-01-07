#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='fluprodia',
    version='1.3',
    license='MIT',
    description='Creating Fluid Proprety Diagrams using CoolProp',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    author='Francesco Witte',
    author_email='fluprodia@witte.sh',
    url='https://github.com/fwitte/fluprodia',

    packages=['fluprodia'] + ['fluprodia.' + p for p in find_packages('src/fluprodia')],
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Utilities',
    ],
    project_urls={
        'Documentation': 'https://fluprodia.readthedocs.io/',
        'Changelog': 'https://fluprodia.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/fwitte/fluprodia/issues',
    },
    keywords=[
        'Fluid Property Diagrams', 'CoolProp', 'TESPy',
    ],
    python_requires='>=3.6.*,<3.9',
    install_requires=[
        'CoolProp>=6.4,<7',
        'matplotlib>=3.2,<4',
        'numpy>=1.13.3,<2'
    ]
)
