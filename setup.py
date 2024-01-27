from setuptools import setup, find_packages
from AMON import __version__ as version
from pathlib import Path

__author__ = 'lozuponelab'
__version__ = version

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

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
      long_description=long_description,
      long_description_content_type='text/markdown',
      author="Michael Shaffer, Kumar Thurimella",
      author_email='lozuponelab.dev@olucdenver.onmicrosoft.com',
      url="https://github.com/lozuponelab/AMON/",
      download_url="https://github.com/lozuponelab/AMON/tarball/%s" % __version__
)
