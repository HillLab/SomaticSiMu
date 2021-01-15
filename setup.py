import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SomaticSiMu", # Replace with your own username
    version="0.0.1",
    author="David Chen",
    author_email="dchen362@uwo.ca",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=["pandas", "numpy"],
    url="https://github.com/HillLab/SomaticSiMu",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.6',
)