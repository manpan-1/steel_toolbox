#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0', 'numpy', 'numpy-stl', 'matplotlib', 'scipy',
]

setup_requirements = [
    'Click>=6.0', 'numpy', 'numpy-stl', 'matplotlib', 'scipy',
]

test_requirements = [
    'Click>=6.0', 'numpy', 'numpy-stl', 'matplotlib', 'scipy',
]

setup(
    name='steeltoolbox',
    version='0.3.0',
    description="A collection of classes and tools used on structural steel research and design.",
    long_description=readme + '\n\n' + history,
    author="Panagiotis Manoleas",
    author_email='manpan@ltu.se',
    url='https://github.com/manpan-1/steel_toolbox',
    packages=find_packages(include=['steel_toolbox']),
    entry_points={
        'console_scripts': [
            'steel_toolbox=steel_toolbox.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='steel_toolbox',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
