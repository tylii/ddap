# -*- coding: utf-8 -*-

from setuptools import setup


with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")


setup(
    name = "ddap",
    version = '1.0.0',
    packages = ["ddap","."],
    package_data = {'ddap': ['data/*']},
    entry_points = {
        "console_scripts": ['ddap = run_ddap:main']
        },
    description = "Docking domain affinity prediction and pathway prediction tool for type I PKS.",
    long_description = long_descr,
    include_package_data=True,
    author = "Tingyang Li",
    install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
        'cairosvg',
        'scikit-learn>=0.20.2',
        'biopython'
    ],
    author_email = "tyli@umich.edu",
    url = "https://github.com/tylii/ddap/",
    )
