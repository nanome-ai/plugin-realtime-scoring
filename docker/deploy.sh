#!/bin/bash

if [ "$(docker ps -aq -f name=realtime-scoring)" != "" ]; then
    echo "removing exited container"
    docker rm -f realtime-scoring
fi

docker run -d \
--name realtime-scoring \
--restart unless-stopped \
-e ARGS="$*" \
realtime-scoring
