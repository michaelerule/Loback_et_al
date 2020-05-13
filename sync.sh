#!/usr/bin/env bash

echo "Staging only changes in tracked files, to avoid commiting cache files and outputs which are too large to store on github..."
git add -u

# Tidy up
# find . -iname "*.pyc" -exec rm {} \;
# find . -path "*.git*" -prune -o -path "./Version 2/Rule/PPC_cache" -prune -o -type f -size -5M | xargs -d "\n" git add

# Commit using the message specified as first argument to this script
echo "Git commit"


MSG=${1:-no_message}  
git commit -m "$MSG"

# Synchronize with master on github
echo "git pull"
git pull

echo "git push"
git push origin master


