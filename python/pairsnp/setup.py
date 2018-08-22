import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pairsnp",
    version="0.0.3",
    author="Gerry Tonkin-Hill",
    author_email="g.tonkinhill@gmail.com",
    description="A simple package for calculating pairwise SNP distances",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gtonkinhill/pairsnp/",
    install_requires=[
          'numpy',
          'scipy'
      ],
    entry_points = {
        'console_scripts': ['pairsnp=pairsnp.pairsnp:main'],
    },
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)