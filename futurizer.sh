# Iterate over over futurizer suite ("--stage1" or "--stage2")
# Atomically commit result of each fixer
#
# Attempt to do the same thing for doctests

for fix in $(futurize $1 --list-fixes | tail -n +2)
do
    echo "futurize --fix $fix"
	futurize --fix $fix  --write --nobackups --no-diffs $2
    git add --all
    git commit -m "futurize --fix $fix"
    
    echo "futurize_doctests $fix"
	python futurize_doctests.py $fix  --write --nobackups --no-diffs $2
    git add --all
    git commit -m "futurize_doctests $fix"
done
