FROM continuumio/miniconda3:4.7.12

RUN mkdir -p /usr/share/man/man1/
WORKDIR /app
ADD environment.yaml environment.yaml

RUN apt-get update && \
    apt-get install -y \
	bash make gcc zlib1g-dev libbz2-dev libcurl4-openssl-dev \
	libssl-dev && \
    apt-get autoremove && apt-get clean \
    && /opt/conda/bin/conda env create -f environment.yaml python=3.6 \
    && /opt/conda/bin/conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
	&& apt-get purge -y --auto-remove gcc make \
	zlib1g-dev libbz2-dev libcurl4-openssl-dev \
	libssl-dev

ENV PATH /opt/conda/envs/pedia/bin:$PATH
ADD . ./
VOLUME ./data

RUN echo "source activate pedia" > ~/.bashrc
ENTRYPOINT ["python", "pedia.py"]
CMD ["-h"]

