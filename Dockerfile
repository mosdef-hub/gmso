ARG PY_VERSION=3.7
FROM continuumio/miniconda3:4.8.2-alpine AS builder
ARG PY_VERSION

EXPOSE 8888

LABEL maintainer.name="mosdef-hub"\
      maintainer.url="https://mosdef.org"

ENV PATH /opt/conda/bin:$PATH

USER root

ADD . /gmso

WORKDIR /gmso

RUN conda update conda -yq && \
	conda config --set always_yes yes --set changeps1 no && \
	. /opt/conda/etc/profile.d/conda.sh && \
	sed -i -E "s/python.*$/python="$PY_VERSION"/" environment-dev.yml && \
	conda env create nomkl -f environment-dev.yml && \
	conda activate gmso-dev && \
        python setup.py install && \
	echo "source activate gmso-dev" >> \
	/home/anaconda/.profile && \
	conda clean -afy && \
	mkdir /home/anaconda/gmso-notebooks && \
	chown -R anaconda:anaconda /gmso && \
	chown -R anaconda:anaconda /opt && \
	chown -R anaconda:anaconda /home/anaconda

WORKDIR /home/anaconda

CMD /bin/su anaconda -s /bin/sh -l
