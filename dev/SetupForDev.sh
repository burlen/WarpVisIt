#!/bin/bash

nerr=0
echo "setting you up for WarpVisIt development..."

# add hook for white space checks
echo -n "installing pre-commit hook..."
cp dev/pre-commit .git/hooks/pre-commit
chmod 755 .git/hooks/pre-commit
if [[ -x .git/hooks/pre-commit ]]
then
    echo "Success!"
else
    echo "Error."
    let nerr=nerr+1
fi

# tweak what the hook checks for
# we don't want tabs or trailing white space
echo -n "define whitespace errors..."
git config core.whitespace trailing-space,space-before-tab,indent-with-non-tab
nok=$?
if [[ "$nok" == "0" ]]
then
    echo "Success!"
else
    echo "Error."
    let nerr=nerr+1
fi

# done, status
echo "setup complted with $nerr errors."
