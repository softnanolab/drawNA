import setuptools

with open("README.md", "r", encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="drawNA",
    version="1.0.2",
    author="Shanil Panara & Debesh Mandal",
    description="Package for creating origami",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/softnanolab/drawNA",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6"
)
