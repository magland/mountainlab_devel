#!/bin/bash

# Install Qt5
echo "Installing Qt5"
sudo apt-add-repository ppa:ubuntu-sdk-team/ppa
sudo apt-get update
sudo apt-get install qtdeclarative5-dev
sudo apt-get install qtbase5-dev qtscript5-dev make g++