#!/usr/bin/env bash

# Add any new files, add all updates to all files
echo "Adding all changes"
git add --all . 
git add -u :/

# Commit using the message specified as first argument to this script
echo "Git commit"
git commit -m "$1"

# Synchronize with master on github
echo "git pull"
git pull

echo "git push"
git push origin master
