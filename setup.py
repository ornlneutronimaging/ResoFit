#!/usr/bin/env python
from setuptools import setup, find_packages
import codecs


def read_file(filename):
    """
    Read a utf8 encoded text file and return its contents.
    """
    with codecs.open(filename, 'r', 'utf8') as f:
        return f.read()


setup(
    name="ResoFit",
    version="0.0.1",
    author="Yuxuan Zhang, Jean C. Bilheux",
    author_email="zhangy6@ornl.gov, bilheuxjm@ornl.gov",
    packages=find_packages(exclude=['tests', 'notebooks']),
    test_suite='tests',

    install_requires=[
        'numpy',
        'pandas',
        'ImagingReso',
        'pprint',
        'scipy',
        'lmfit',
        'peakutils',
    ],
    dependency_links=[
    ],
    description="Fitting tool for neutron resonance imaging",
    long_description=read_file('README.rst'),
    license='BSD',
    keywords=['neutron', 'resonance', 'fitting', 'imaging'],
    url="https://github.com/ornlneutronimaging/ResoFit.git",
    classifiers=['Development Status :: 3 - Alpha',
                 'License :: OSI Approved :: BSD License',
                 'Topic :: Scientific/Engineering :: Physics',
                 'Intended Audience :: Developers',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3.5',
                 'Programming Language :: Python :: 3.6',
                 'Natural Language :: English'],
)


# End of file
