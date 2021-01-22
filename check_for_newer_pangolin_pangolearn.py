#!/usr/bin/env python

import subprocess
from urllib import request
import json

def get_installed_pangolin_version():
    """
    Get the version information for the currently installed version of
    pangolin and pangolearn
    """
    # get version information for pangolin, pangolearn model
    pangolin_version = subprocess.run("pangolin --version".split(),
                                      stdout=subprocess.PIPE)
    pangolin_version = pangolin_version.stdout.decode('utf-8').strip()
    # add v to match release tag
    pangolin_version = pangolin_version.replace(' ', ' v')

    pangolearn_version = subprocess.run("pangolin --pangoLEARN-version".split(),
                                        stdout=subprocess.PIPE)
    pangolearn_version = pangolearn_version.stdout.decode('utf-8').strip()
    # add data release to match release tag
    pangolearn_version = pangolearn_version.replace(' ', ' data release ')

    return pangolin_version, pangolearn_version


def get_latest_pangolin_release_version():
    """
    Check for the latest release of pangolin and pangolearn on github
    """

    pangolin_release = request.urlopen('https://api.github.com/repos/cov-lineages/pangolin/releases')
    pangolin_release = json.load(pangolin_release)
    pangolin_release = pangolin_release[0]['name']

    pangolearn_release = request.urlopen('https://api.github.com/repos/cov-lineages/pangoLEARN/releases')
    pangolearn_release = json.load(pangolearn_release)
    pangolearn_release = pangolearn_release[0]['name']

    return pangolin_release, pangolearn_release


if __name__ == '__main__':

    # get currently installed versions
    pangolin_version, pangolearn_version = get_installed_pangolin_version()

    # get latest release version
    pangolin_release, pangolearn_release = get_latest_pangolin_release_version()

    if pangolin_version != pangolin_release:
        print(f"Newer pangolin version available: {pangolin_release}")
        print("Consider running: `pip install git+https://github.com/cov-lineages/pangolin.git --upgrade`")
    else:
        print(f"Latest pangolin already installed: {pangolin_release}")

    if pangolearn_version != pangolearn_release:
        print(f"Newer pangoLEARN model available: {pangolearn_release}")
        print("Consider running: `pip install git+https://github.com/cov-lineages/pangoLEARN.git --upgrade`")
    else:
        print(f"Latest pangoLEARN model already installed: {pangolearn_version}")

