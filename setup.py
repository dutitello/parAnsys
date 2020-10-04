import setuptools

try:
    import numpy
except:
    raise Exception('Please install NumPy first.')

try:
    import scipy
except:
    raise Exception('Please install SciPy first.')


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="paransys",
    version="0.1",
    author="Eduardo Pagnussat Titello",
    author_email="du.titello+paransys@gmail.com",
    description="PARANSYS: Python pArametric Reliability Analysis on ANSYS",
    url="https://github.com/dutitello/parAnsys",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['paransys'],
    python_requires='>=3.6',
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3",
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Microsoft :: Windows',
    ],
    install_requires=['numpy>=1.18.0',
                      'scipy>=1.17']
)
