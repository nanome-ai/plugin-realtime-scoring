FROM continuumio/miniconda3

ENV PLUGIN_SERVER=plugins.nanome.ai

COPY . /app
WORKDIR /app

RUN conda install -y -c openbabel openbabel
RUN pip install nanome

CMD python -m nanome_realtime_scoring.RealtimeScoring -a ${PLUGIN_SERVER}