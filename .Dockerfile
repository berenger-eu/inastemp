FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y wget
RUN apt-get install -y g++
RUN apt-get install -y cmake
RUN apt-get install -y git
RUN apt-get install -y clang
RUN apt-get install -y g++-4.8
COPY sde-external-8.35.0-2019-03-11-lin.tar.bz2 ./
RUN tar -xvf sde-external-8.35.0-2019-03-11-lin.tar.bz2
RUN ln -s $(pwd)/sde-external-8.35.0-2019-03-11-lin/sde64 /usr/bin/sde64
