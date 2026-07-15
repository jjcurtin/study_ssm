# Starting point
FROM --platform=linux/amd64 python:3.10.14

# Install R packages
RUN pip install --upgrade pip
COPY py_packages.txt .
RUN pip install -r py_packages.txt