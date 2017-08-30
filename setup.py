from setuptools import setup
from setuptools import find_packages

def setup_package():
  install_requires = ['pandas', 'scipy', 'numpy', 'sklearn', 'argparse', 'h5py']
  console_scripts = [ 'mofa=mofa.run.template_run:entry_point'],
  metadata = dict(
      name = 'MOFA',
      version = '0.1',
      description = 'Multi-Omics Factor Analysis',
      #long_description=read('README.rst'),
      url = 'http://github.com/rargelaguet/mofa',
      author = 'Ricard Argelaguet, Damien Arnol and Britta Velten',
      author_email = 'ricard.argelaguet@gmail.com',
      license = 'MIT',
      packages = find_packages(),
      install_requires = install_requires,
      entry_points = {'console_scripts': console_scripts}
    )

  setup(**metadata)

if __name__ == '__main__':
  setup_package()