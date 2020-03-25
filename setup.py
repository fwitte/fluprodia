from setuptools import setup
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='fluprodia',
      version='0.0.1 dev',
      description='TESPy fluid property plotting',
      url='http://github.com/fwitte/tespy.plotting',
      author='Francesco Witte',
      author_email='francesco.witte@web.de',
      long_description=read('README.rst'),
      license='MIT',
      packages=[],
      python_requires='>=3', '<3.8',
      install_requires=['CoolProp>=6,<7',
                        'numpy>=1.13.3,<2',
                        'matplotlib>=3.2,<4'])
