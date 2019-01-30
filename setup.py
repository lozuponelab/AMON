from setuptools import setup, find_packages
from AMON import __version__ as version

__author__ = 'shafferm'
__version__ = version

setup(
      name="AMON-bio",
      version=__version__,
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      install_requires=['scipy', 'biom-format', 'pandas', 'matplotlib', 'statsmodels', 'numpy', 'aiohttp', 'seaborn',
                        'matplotlib-venn', 'KEGG-parser'],
      scripts=['scripts/amon.py', 'scripts/extract_ko_genome_from_organism.py'],
      packages=find_packages(),
      description="Annotation of Metabolite Origin via Networks: A tool for predicting putative metabolite origins for"
                  "microbes or between microbes and host with or without metabolomics data",
      author="Michael Shaffer",
      author_email='michael.shaffer@ucdenver.edu',
      url="https://github.com/shafferm/AMON/",
      download_url="https://github.com/shafferm/AMON/tarball/%s" % __version__
)
