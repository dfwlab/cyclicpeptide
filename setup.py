from setuptools import setup, find_packages

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()
setup(
    name='cyclicpeptide',
    version='1.3.0',
    packages=find_packages(),
    description='A python package for cyclic peptides drug design',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Willow0316/cyclicpeptide/tree/master',
    author='Willow',
    author_email='willow.yang.b@gmail.com',
    python_requires='>=3.8',  # Python 版本要求
    license='MIT',
    install_requires=[
        'matplotlib==3.8',
        'networkx==3.1',
        'numpy==1.24.4',
        'pandas==2.0.3',
        'rdkit==2023.9.5',
        'scipy==1.13.0'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
    package_data={'': ['*.csv', '*.txt','*.tsv']}, #这个很重要
    include_package_data=True #也选上
)