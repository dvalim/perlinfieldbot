FROM ubuntu:20.04

RUN apt-get update && apt-get install -y gcc g++ cmake libsfml-dev
RUN apt-get install -y xutils
COPY . /perlinfieldbot

WORKDIR /build

