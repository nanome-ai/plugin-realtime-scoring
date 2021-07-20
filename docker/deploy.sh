#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

existing=$(docker ps -aq -f name=realtime-scoring)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run \
--name realtime-scoring \
--restart unless-stopped \
-e ARGS="$*" \
--mount type=bind,source=/home/mike/workspace/nanome-lib/nanome,target=/opt/conda/lib/python3.7/site-packages/nanome \
--mount type=bind,source=/home/mike/workspace/plugin-realtime-scoring,target=/app \
realtime-scoring
