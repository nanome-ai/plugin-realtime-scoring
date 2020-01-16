if [ "$(docker ps -q -f name=realtime_scoring)" != "" ]; then
    if [ "$(docker ps -aq -f status=exited -f name=realtime_scoring)" != "" ]; then
        # cleanup
        docker rm realtime_scoring
    fi
    docker rm -f realtime_scoring
fi

if [ "$1" != "" ]; then
    echo "Using specified plugin server: $1"
    docker run -d -p 8888:8888 -e PLUGIN_SERVER=$1 --name realtime_scoring realtime-scoring
else
    echo "Using default plugin server: plugins.nanome.ai"
    docker run -d -p 8888:8888 --name realtime_scoring realtime-scoring
fi