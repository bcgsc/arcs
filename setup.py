import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="arcs",
    version="1.2.7",
    description="Scaffolding genome sequence assemblies using linked or long read sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bcgsc.github.io/arcs/",
    license="GPLv3",
    python_requires=">=3",
    scripts=[
        "Examples/makeTSVfile.py"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
