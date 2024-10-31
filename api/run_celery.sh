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
  echo "Running in debug mode...log file will be saved to celery.log"
  nohup celery -A app.celery worker -P gevent -c 1 --loglevel DEBUG -Q molecular_docking,global_docking > celery.log 2>&1 &
else
  echo "Running in production mode...will not log to file"
  nohup celery -A app.celery worker -P gevent -c 1 --loglevel DEBUG -Q molecular_docking,global_docking > /dev/null 2>&1 &
fi
