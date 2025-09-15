from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()
    
setup(
    name = "spqrp",
    version = "0.1.0",
    packages=find_packages(),
    install_requires=requirements,
    package_data={
        'spqrp': ['data/*.csv'],
    },
    include_package_data=True,
)
