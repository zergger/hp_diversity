from setuptools import setup, find_packages

setup(
    name="hp_diversity",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'scipy'
    ],
    author="Liu, K.",
    author_email="xqa9038@dingtalk.com",
    description="A package for calculating mtDNA heteroplasmy alpha, beta, and gamma diversity.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/zergger/hp_diversity",
)
