#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

# with open('README.rst') as readme_file:
#     readme = readme_file.read()
#
# with open('HISTORY.rst') as history_file:
#     history = history_file.read()

# requirements = [
#     'Click>=6.0',
#     'matplotlib',
#     'numpy',
#     'networkx',
# ]

# test_requirements = [
#     # TODO: put package test requirements here
#     'hypothesis',
#     'pytest'
# ]

setup(
    name='mbtools',
    version='2016.8.1',
    description="Graph Visualization Package",
    # long_description=readme + '\n\n' + history,
    author="Eric J. Ma",
    author_email='ericmajinglong@gmail.com',
    url='https://github.com/ericmjl/nxviz',
    packages=[
        'mbtools',
    ],
    # package_dir={'src':
    #              'mbtools'},
    # entry_points={
    #     'console_scripts': [
    #         'nxviz=nxviz.cli:main'
    #     ]
    # },
    # include_package_data=True,
    # install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    # keywords='nxviz',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    # test_suite='tests',
    # tests_require=test_requirements
)
