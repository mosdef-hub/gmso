from setuptools import setup, find_packages
from pathlib import Path

#####################################
VERSION = "0.2.0"
ISRELEASED = False
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + '.dev0'
#####################################


with open(Path(__file__).parent / 'requirements.txt') as req_file:
    requirements = list(filter(None, req_file.read().split('\n')))


setup(
    name='gmso',
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=True,
    author='Matthew W Thompson, Justin Gilmer',
    author_email='matt.thompson@vanderbilt.edu, justin.b.gilmer@vanderbilt.edu',
    url='https://github.com/mosdef-hub/gmso',
    download_url='https://github.com/mosdef-hub/gmso/tarball/{}'.format(__version__),
    license="MIT",
    keywords='gmso',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
    install_requires=requirements,
    python_requires='>=3.6, <4',
)
