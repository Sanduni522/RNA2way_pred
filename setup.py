#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


with open("README.md", "r", encoding="UTF-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="UTF-8") as f:
    requirements = f.read().splitlines()

setup(
    name="RNA2way_pred",
    version="0.1.0",
    description="A script to predict the DMS chemical reactivity of A and C in two-way junctions using the structural parameters of these nucleotides."
                "This script will run rosetta-farfar to model pdb structures and calculate the structural parameters of the nucleotides of the models."
                "Finally, using these structural parameters as inputs, using Deep neural networks, the DMS reactivity will be predicted."
                "It will also provide the correlation between the predicted DMS reactivity and the Experimental DMS reactivity for a particualr nucleotide.",
    long_description=readme,
    long_description_content_type="test/markdown",
    author="Sanduni Deenalattha",
    author_email="sdeenalattha2@huskers.unl.edu",
    url="https://github.com/jyesselm/RNA2way_pred",
    packages=[
        "RNA2way_pred",
    ],
    package_dir={"RNA2way_pred": "RNA2way_pred"},
    py_modules=["RNA2way_pred/cli"],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords="RNA2way_pred",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={"console_scripts": ["RNA2way_pred = RNA2way_pred.cli:cli"]},
)
