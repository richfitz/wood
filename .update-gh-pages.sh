#!/bin/bash
# Based on
# http://sleepycoders.blogspot.com.au/2013/03/sharing-travis-ci-generated-files.html
# https://help.github.com/articles/creating-project-pages-manually
#
# This differs by deleting the history of the gh pages branch by doing
# a force push off a brand new branch.  This helps keep the size of
# the repo under control.
echo -e "Preparing to copy generated files to gh-pages branch"
if [ "$TRAVIS_PULL_REQUEST" == "false" ]; then
    if [ "$USE_PACKRAT" == "1" ]; then
        echo -e "This is a the packrat version, not copying files"
    else
        echo -e "Starting to update gh-pages\n"

        mkdir -p $HOME/keep
        cp -R index.html wood.html stylesheet.css figure doc/wood-ms.pdf doc/wood-ms-supporting.pdf $HOME/keep

        #go to home and setup git
        cd $HOME
        git config --global user.email "rich.fitzjohn@gmail.com"
        git config --global user.name "Rich FitzJohn"

        echo -e "Recloning project"
        # Reclone the project, using the secret token.  Uses /dev/null to avoid leaking decrypted key
        git clone --quiet --branch=gh-pages --single-branch https://${GH_TOKEN}@github.com/richfitz/wood.git gh-pages > /dev/null

        cd gh-pages

        # Move the old branch out of the way and create a new one:
        git branch -m gh-pages-old
        git checkout --orphan gh-pages

        # Delete all the files and replace with our good set
        git rm -rf .
        cp -Rf $HOME/keep/* .

        # add, commit and push files
        git add -f .
        git commit -m "Travis build $TRAVIS_BUILD_NUMBER pushed to gh-pages"
        echo -e "Pushing to origin/gh-pages"
        git push -fq origin gh-pages > /dev/null

        echo -e "Uploaded generated files to gh-pages\n"
    fi
else
    echo -e "This is a pull request, not copying files"
fi

echo -e "Done"
