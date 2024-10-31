#!/bin/bash

# 设置调试模式标志
DEBUG=false

for arg in "$@"; do
  if [ "$arg" == "--debug" ]; then
    DEBUG=true
    break
  fi
done


# 在调试模式下输出详细信息
if [ "$DEBUG" == true ]; then
  echo "Running in debug mode...log file will be saved to app.log"
  nohup gunicorn --bind 0.0.0.0:5050 --workers 1 --worker-class gevent --timeout 200 --preload app:app > app.log 2>&1 &
else
  echo "Running in production mode...will not log to file"
  nohup gunicorn --bind 0.0.0.0:5050 --workers 1 --worker-class gevent --timeout 200 --preload app:app > /dev/null 2>&1 &
fi
