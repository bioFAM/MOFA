from setuptools import setup
from setuptools import find_packages

if __name__ == '__main__':
    setup(name='MOFA',
          version='0.1',
          description='Multi-Omics Factor Analysis',
          #long_description=read('README.rst'),
          url='http://github.com/rargelaguet/mofa',
          author='Ricard Argelaguet, Damien Arnol and Britta Velten',
          author_email='ricard.argelaguet@gmail.com',
          license='MIT',
          packages=find_packages(),
		      install_requires=[
            'pandas>=0.20',
            'scipy>=0.19',
            'scikit-learn>=0.18',
            'argparse>=1.4',
            'h5py >=2.7.0',
          ],
        )
