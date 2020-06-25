from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='fatiguepy',
      version='1.8.1',
      description='Package to estimate life of random fatigue history with frequency domain methods',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/thiagonds65/fatiguepy',
      author='Thiago Santos',
      author_email='thiagonds@id.uff.br',
      license='MIT',
      packages=['fatiguepy'])

