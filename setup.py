from setuptools import setup

with open('README.md') as fp:
    LONG_DESCRIPTION = fp.read()

setup(
      name='DissNet',
      maintainer = 'Davide Momi',
      maintainer_email = 'momi.davide89@gmail.com',
      version='1.0-alpha',
      description='brain network stimulation',
      long_description = LONG_DESCRIPTION,
      author='Daive Momi',
      author_email='momi.davide89@gmail.com',
      twitter='https://twitter.com/DaveMomi',
      website='https://davi1990.github.io/',
      license='MIT',
      zip_safe=False,
      packages=['DissNet'],
      url = 'https://github.com/Davi1990/DissNet',
      keywords = ['brain network','dwi','imaging','fMRI','diffusion','structual'],
      install_requires=[
          'numpy',
          'mayavi',
          'pysurfer',
          'pandas',
          'matplotlib',
          'scipy',
          'mne',
          'networkx',
          'bct'
      ],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'Topic :: Scientific/Engineering',
          'Operating System :: MacOS',
          'Operating System :: Unix',
          'Topic :: Software Development :: Build Tools',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7'

      ],
      )
