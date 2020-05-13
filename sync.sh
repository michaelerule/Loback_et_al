#!/usr/bin/env bash


sdfgdsfgdsgfdgds

# Tidy up
find . -iname "*.pyc" -exec rm {} \;

# Add any new files, add all updates to all files
echo "Adding all changes"
git add -A . 
#git add -u :/

# un-add files that are too big to store on github

echo "Not committing these files because they are too large:"
IFS='
'
for i in $(find . -type f -size +5M); do
    git reset "$i"
done

# Commit using the message specified as first argument to this script
echo "Git commit"
git commit -m "..."

# Synchronize with master on github
echo "git pull"
git pull

echo "git push"
git push origin master

# We need a different strategy for git add

find . -path "*.git*" -prune -o -path "./Version 2/Rule/PPC_cache" -prune -o -type f -size -5M | xargs -d "\n" git add

for i in $(find . -path ./.git -prune -o -path "./Version 2/Rule/PPC_cache" -prune -o -type f -size -5M); do
    git add "$i"
done
git ls-files --deleted | xargs -d "\n" git add


find . -path "./.git" -prune -o -type f -size +5M | xargs -d "\n" git rm --cached
