from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(name="ImagingCytometryTools",
      version="1.0.4",
      description="Analysis of complex imaging data",
      long_description=long_description,
      url="https://github.com/NiklasVesper/ImagingCytometryTools",
      author="Niklas Vesper",
      author_email="niklas.vesper@uniklinik-freiburg.de",
      packages=find_packages(),
      include_package_data=True,
      license="LICENSE.txt",
      install_requires=open("requirements.txt").read(),
      classifiers=["Programming Language :: Python :: 3.8", "License :: MIT License", "Topic :: Scientific/Engineering :: Bio-Informatics", "Topic :: Scientific/Engineering :: Visualization",], 
      python_requires=">=3.8")
