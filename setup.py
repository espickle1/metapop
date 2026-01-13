import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="metapop",
    version="1.0.0",
    author="Kenji Gerhardt",
    author_email="kenji.gerhardt@gmail.com",
    description="A metagenomic population statistics pipeline.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    py_modules=["metapop"],
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "pysam>=0.16.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "scikit-learn>=0.24.0",
    ],
    extras_require={
        "colab": [
            "plotly>=5.0.0",
            "ipywidgets>=7.6.0",
        ],
        "full": [
            "scikit-bio>=0.5.6",
            "plotly>=5.0.0",
            "ipywidgets>=7.6.0",
            "adjustText>=0.7.3",
        ],
    },
    entry_points={
        "console_scripts": [
            "metapop=metapop.metapop_main:main",
        ]
    }
)