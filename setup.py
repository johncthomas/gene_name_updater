from os import path
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

def main():
    setup(name='HGNC_converter',
          version='0.1.1',
          description='Get HGNC symbols, id, etc from obsolete symbols',
          packages=find_packages(),
          install_requires=['pandas', 'biopython'],
          include_package_data=True,
          package_data={'HGNC_converter': ['complete_HGNC.tsv']})


if __name__ == '__main__':
    main()