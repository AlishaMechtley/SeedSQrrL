#!/bin/bash

touch assembly.log resources.log

echo ">> Beginning automated assembly process." >> assembly.log

python automate_mitobim.py >> assembly.log &
AUTO_PID=$!

while kill -0 ${AUTO_PID};
do
  ps --forest -o pid,ppid,c,rss,time,cmd -g $(ps -o sid= -p ${AUTO_PID}) >> resources.log;
  sleep 20;
done;

echo ">> Automated assembly process complete." >> assembly.log