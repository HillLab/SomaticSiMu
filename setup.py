import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SomaticSiMu", 
    version="2.0.0",
    author="David Chen",
    author_email="dchen362@uwo.ca",
    description="Mutational Signature Simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=["pandas", "numpy"],
    url="https://github.com/HillLab/SomaticSiMu",
    packages=setuptools.find_packages(),
    classifiers=[
	"Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
Public License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.8',
)