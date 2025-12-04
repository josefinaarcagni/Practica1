from setuptools import setup, find_packages

setup(
    name="persephone",
    version="2025.1.0",
    packages=find_packages(), 
    author="Telmo Blasco, ...",
    author_email="tblasco@unav.es, ...",
    description="Package to ...",
    license='',
    url='',
    install_requires=[
        "pandas",
        "pyarrow>=12.0.1",
        "numpy",
        "cobra",
        "scipy>=1.7, <1.11",
        "tqdm"
    ],
    include_package_data=True,
    package_data={
        'persephone': ['data/**/*'],
    },
)
