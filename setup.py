#from distutils.core import setup
from setuptools import setup

setup(
    name='scrapePDB',
    description='Scrape the DPB data base to create cutom dataset',
    version='0.0',
    url='https://github.com/DeepRank/scrapePDB',
    packages=['scrapePDB'],


    install_requires=[
        'numpy >= 1.13'],
        #'tarfiles',
        #'pickle'],

    extras_require= {
        'test': ['nose', 'coverage', 'pytest',
                 'pytest-cov','codacy-coverage','coveralls'],
    }
)
