from setuptools import setup, find_packages
import pathlib
import pkg_resources

version = "0.0.1"

with open("README.md") as readme_file:
    readme = readme_file.read()

with pathlib.Path("requirements.txt").open() as requirement_file:
    requirements = [ 
        str(req) for req in pkg_resources.parse_requirements(requirement_file)
    ]

setup(name = "LangevinDynamics", 
        version = version, 
        description=(
            """
            openmm implementation of simple Langevin Dynamics designed to allow us to perform various operations such as VES
            """
        ),
        long_description=readme, 
        author="Yusheng Cai", 
        author_email="cys9741@seas.upenn.edu", 
        python_requires=">=3.9", 
        install_requires=requirements, 
        packages=find_packages()
)

