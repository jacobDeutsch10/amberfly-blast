from setuptools import setup, find_packages

setup(name='amberfly-blast',
      version='0.1.0',
      description='Blast analysis of data from PATMOS',
      url='https://github.com/jacobDeutsch10/amberfly-blast',
      author='Jacob Deutsch',
      license='MIT',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      install_requires=[
          'numpy>=1.16.1',
          'mathutils~=2.81.2'
      ],
      zip_safe=False)
