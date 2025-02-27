from setuptools import setup, find_packages

setup(
    name="gene_nomenclature",
    version="0.1",
    packages=find_packages(),
    install_requires=["pandas", "numpy"],
    package_data={"gene_nomenclature": ["hugo_export.txt"]},
    description="A package to look up gene nomenclature information.",
    author="Miles Weatherseed",
    author_email="miles.weatherseed@gmail.com",
)
