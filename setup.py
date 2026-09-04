from setuptools import setup
from pathlib import Path

setup(
    name='codonslider',
    version='0.1.0',
    packages=['codonslider'],
    url='https://github.com/bsimeonovb/codonslider',
    license='MIT',
    author='Boris Simeonov',
    author_email='bsimeonovb@gmail.com',
    description='Bioinformatic tool for codon stretch detection in CCDS sequences',
    long_description=Path('README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
    python_requires='>=3.9',
    install_requires=[
        'biopython',
        'pandas',
    ],
)
