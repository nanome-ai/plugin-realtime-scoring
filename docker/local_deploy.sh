if [ "$(docker ps -aq -f name=realtime-scoring)" != "" ]; then
    # cleanup
    echo "removing exited container"
    docker rm -f realtime-scoring
fi

docker run -d \
--restart always \
-e ARGS="$*" \
--name realtime-scoring realtime-scoring
