import setuptools
from io import open
import sys
from subprocess import call
from os import path
from os import mkdir


if sys.version_info[0] < 3:
    sys.stdout.write('metaclock requires Python 3 or higher. Please update you Python installation')


install_reqs = ["biopython", "matplotlib", "numpy", "pandas", "seaborn"]

setuptools.setup(name='metaclock',
                 version='1.0.2',
                 author='Kun D. Huang',
                 author_email='kun.huang@unitn.it',
                 # url='',
                 license='A-GPL 3.0',
                 # scripts=['phylophlan/phylophlan_write_default_configs.sh'],
                 packages=setuptools.find_packages(),
                 package_data={
                     'metaclock': [
                         'metaclock/utils/*',
                         'metaclock/metaclock_configs/*'
                 ]},
                 entry_points={
                     'console_scripts': [
                         'metaclock_mac = metaclock.metaclock_mac:main',
                         'metaclock_visualizer = metaclock.metaclock_visualizer:main',
                         'metaclock_tailor = metaclock.metaclock_tailor:main',
                         'metaclock_combiner = metaclock.metaclock_combiner:main',
                         'bo6_screen.py = metaclock.utils.bo6_screen:main',
                         'filter.py = metaclock.utils.filter:main',
                         'metaclock_mac_template_configs = metaclock.utils.metaclock_mac_template_configs:main'
                 ]},
                 description='metaclock package test',
                 include_package_data=True,
                 long_description=open('README.md').read(),
                 long_description_content_type='text/markdown',
                 install_requires=install_reqs,
                 zip_safe=False)
