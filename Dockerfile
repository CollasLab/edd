FROM python:2

RUN apt-get update
# Needed for scipy
RUN apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran bedtools

# This needs to run before statsmodels (and before scipy?)
# (Modern versions work as of Aug 2022, but
# freezing here just in case)
RUN pip install Cython==0.29.32 numpy==1.16.6

# Needs to be done before statsmodels
# (version doesn't seem to matter; statsmodels requires >=0.18)
RUN pip install scipy==1.2.3

# 0.11 incompatible with python 2
RUN pip install statsmodels==0.10.2

# 0.9 seems incompatible too
RUN pip install pybedtools==0.8.2

# readme states this is required
RUN pip install pysam==0.11.2.2

RUN pip install edd

CMD ["edd"]