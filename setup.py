from setuptools import setup

#####################################
VERSION = "0.0.0"
ISRELEASED = False 
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + '.dev0'
#####################################


setup(
    name='topology',
    version=__version__,
    author='Matthew W Thompson, Justin Gilmer',
    author_email='matt.thompson@vanderbilt.edu, justin.b.gilmer@vanderbilt.edu',
    url='https://github.com/mattwthompson/topology',
    download_url='https://github.com/mattwthompson/topology/tarball/{}'.format(__version__),
    package_dir={'topology': 'topology'},
    license="MIT",
    zip_safe=False,
    keywords='topology',
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
)

