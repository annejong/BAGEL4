#!/bin/bash

# remove files older then a 60 days, using www-data user

# cd /tmpdrive
# sudo -u www-data /data/bagel4/00.clean_tmp_folder.sh
cd /tmpdrive
find /tmpdrive/bagel4wrapper/ -name "*" -mtime +60 -exec rm -rf {} \;

# remove empty folders (rmdir will only remove empty folders)
#find /tmpdrive/bagel4wrapper -name "*" -mtime +60 -type d -empty -exec rmdir {} \;



# check how many files are removed
df -h
ls -f1 /tmpdrive/bagel4wrapper | wc -l


