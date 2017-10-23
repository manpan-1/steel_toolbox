#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0', 'numpy',
    # TODO: put package requirements here
]

setup_requirements = [
    # TODO(manpan-1): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='steel_toolbox',
    version='0.1.0',
    description="A collection of classes and tools used on structural steel design.",
    long_description=readme + '\n\n' + history,
    author="Panagiotis Manoleas",
    author_email='manpan@ltu.se',
    url='https://github.com/manpan-1/steel_toolbox',
    packages=find_packages(include=['steel_toolbox'], exclude=['docs', 'tests*']),
    entry_points={
        'console_scripts': [
            'steel_toolbox=steel_toolbox.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='structural steel eurocode',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
