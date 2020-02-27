#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

requirements = [
    "kipoi>=0.5.5",
    # "genomelake",
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
    "kipoi-conda>=0.1.0"
]

test_requirements = [
    "bumpversion",
    "wheel",
    "epc",
    "jedi",
    "pytest>=3.3.1",
    "pytest-xdist",  # running tests in parallel
    "pytest-pep8",  # see https://github.com/kipoi/kipoi/issues/91
    "pytest-cov",
    "coveralls",
    "scikit-learn",
    "cython",
    "cyvcf2",
    "pyranges>=0.0.71",
    # "genomelake",
    "keras",
    "tensorflow",
    "pybedtools"
]

setup(
    name='kipoiseq',
    version='0.3.3',
    description="kipoiseq: sequence-based data-laoders for Kipoi",
    author="Kipoi team",
    author_email='avsec@in.tum.de',
    url='https://github.com/kipoi/kipoiseq',
    long_description="kipoiseq: sequence-based data-laoders for Kipoi",
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
    tests_require=test_requirements
)
