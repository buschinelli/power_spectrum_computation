import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Bicar",
    version="0.0.1",
    author="Joao Buschinelli",
    author_email="joao.correabuschinelli@mail.mcgill.ca",
    description="Bi-cartesian power spectra computation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/buschinelli/power_spectrum_computation/packing/Bicar/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)



