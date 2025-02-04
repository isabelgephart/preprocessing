import setuptools

setuptools.setup(
    name="preproc_functions",
    version="0.1.0",
    author="Isabel Gephart",
    author_email="igephart@uchicago.edu",
    description="preprocessing software",
    url="https://github.com/isabelgephart/preprocessing.git",
    packages=setuptools.find_packages(),
    install_requires=[
        'pandas',
        'matplotlib',
        'numpy',
        'scipy',
    ],
)
