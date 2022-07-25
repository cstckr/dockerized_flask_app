FROM python:3.9.12-bullseye
WORKDIR /23_deployment
ADD . /23_deployment

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version


RUN pip install -r requirements.txt
RUN conda install -c conda-forge rdkit==2020.09.3
CMD ["flask", "run", "-h", "0.0.0.0", "-p", "5000"]