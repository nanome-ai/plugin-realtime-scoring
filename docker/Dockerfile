FROM nanome/plugin-env

ENV ARGS=''

ARG CACHEBUST
RUN pip install git+https://github.com/nanome-ai/nanome-lib#nanome

COPY . /app
WORKDIR /app

CMD python3 run.py ${ARGS}
