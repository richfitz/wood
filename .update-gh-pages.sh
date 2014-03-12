#!/bin/bash
# http://sleepycoders.blogspot.com.au/2013/03/sharing-travis-ci-generated-files.html
echo -e "Preparing to copy generated files to gh-pages branch"
if [ "$TRAVIS_PULL_REQUEST" == "false" ]; then
  echo -e "Starting to update gh-pages\n"

  mkdir -p $HOME/keep
  cp -R wood.html figure doc/wood-ms.pdf $HOME/keep

  #go to home and setup git
  cd $HOME
  git config --global user.email "rich.fitzjohn@gmail.com"
  git config --global user.name "Rich FitzJohn"

  #using token clone gh-pages branch
  git clone --quiet --branch=gh-pages https://${GH_TOKEN}@github.com/richfitz/wood.git gh-pages > /dev/null

  #go into diractory and copy data we're interested in to that directory
  cd gh-pages
  cp -Rf $HOME/keep/* .

  #add, commit and push files
  git add -f .
  git commit -m "Travis build $TRAVIS_BUILD_NUMBER pushed to gh-pages"
  git push -fq origin gh-pages > /dev/null

  echo -e "Uploaded generated files to gh-pages\n"
else
  echo -e "This is a pull request, not copying files"
fi

echo -e "Done"
