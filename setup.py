#!/usr/bin/env python3
"""Setup downloader-cli"""

import io
from setuptools import setup

with io.open("README.md", encoding='utf-8') as fh:
    long_description = fh.read()

# Set your packages requirements here.
requirements = [
    'urllib3>=1.25.6'
]

exec(open("downloader-cli/__version__.py").read())


setup(
    name="downloader-cli",
    version=__version__,
    author="YOUR NAME ",
    author_email="YOUR_EMAIL_HERE@vhio.net",
    description="A simple sample code written in Python  as demo.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.dev/radiomicsvhio/repo_template",
    packages=["downloader-cli"],
    classifiers=(
        [
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ]
    ),
    entry_points={
        'console_scripts': [
            "dw = downloader-cli.download:main"
        ]
    },
    install_requires=requirements,
)