#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

requirements = [
    "kipoi>=0.5.5",
    "pyfaidx",
    "numpy",
    "pandas",
    "tqdm",
    # "colorlog",
    # "related>=0.6.0",
    # sometimes required
    # "h5py",
    "gffutils",
    "kipoi-utils>=0.1.1",
    "kipoi-conda>=0.1.0",
    "pyranges"
]

test_requirements = [
    "bumpversion",
    "wheel",
    "epc",
    "jedi",
    "pytest>=3.3.1",
    "pytest-xdist",  # running tests in parallel
    "pytest-pep8",  # see https://github.com/kipoi/kipoi/issues/91
    "pytest-mock",
    "pytest-cov",
    "coveralls",
    "scikit-learn",
    "cython",
    "cyvcf2",
    "pyranges>=0.0.71",
    # "keras",
    # "tensorflow",
    "pybedtools",
    # "concise"
]

setup(
    name='kipoiseq',
    version='0.7.1',
    description="kipoiseq: sequence-based data-loaders for Kipoi",
    author="Kipoi team",
    author_email='avsec@in.tum.de',
    url='https://github.com/kipoi/kipoiseq',
    long_description="kipoiseq: sequence-based data-loaders for Kipoi",
    packages=find_packages(),
    install_requires=requirements,
    extras_require={
        "develop": test_requirements,
    },
    license="MIT license",
    zip_safe=False,
    keywords=["model zoo", "deep learning",
              "computational biology", "bioinformatics", "genomics"],
    test_suite='tests',
    include_package_data=False,
    tests_require=test_requirements,
    python_requires='>=3.6'
)
