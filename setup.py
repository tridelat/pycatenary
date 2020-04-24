import setuptools
from distutils.core import setup

with open('README.md') as file:
    long_description = file.read()

setup(name='pycatenary',
      packages=['pycatenary'],
      version='0.1',
      description='A Python library for solving catenary equations',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Tristan de Lataillade',
      author_email='delataillade.tristan@gmail.com',
      #url='https://github.com/tridelat/pycatenary',
      #download_url='https://github.com/tridelat/pycatenary/...',
      keywords=['catenary', 'mooring', 'python'],
      classifiers=["Programming Language :: Python :: 3",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: OS Independent",
      ],
)
