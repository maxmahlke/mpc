from setuptools import setup

setup(
    name='mpc observations parser',
    version='0.1',
    py_modules=['mpc'],
    install_requires=[
        'astropy',
        'click',
        'numpy',
        'pandas',
        'tqdm'
    ],
    entry_points='''
        [console_scripts]
        mpc=mpc:cli
    ''',
)
