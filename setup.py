from setuptools import setup, find_packages

__author__ = 'shafferm'
__version__ = '0.0.0'

setup(
      name="microMetabPred",
      version=__version__,
      setup_requires=['pytest-runner'],
      test_require=['pytest'],
      install_requires=['scipy', 'biom-format', 'pandas', 'matplotlib', 'statsmodels', 'numpy', 'aiohttp', 'seaborn',
                        'matplotlib-venn'],
      scripts=['scripts/MicroMetabPred.py', 'scripts/extract_ko_genome_from_organism.py'],
      packages=find_packages(),
      description="A tool for predicting putative metabolite origins for microbes or between microbes and host with or"
                  " without metabolomics data",
      author="Michael Shaffer",
      author_email='michael.shaffer@ucdenver.edu',
      url="https://github.com/shafferm/microMetabPred/",
      download_url="https://github.com/shafferm/microMetabPred/tarball/%s" % __version__
)
