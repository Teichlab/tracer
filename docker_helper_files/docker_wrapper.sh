#!/bin/bash
set -e

#the environmental variable thing apparently can't be set when making the container
export IGDATA=/ncbi-igblast-1.7.0/bin
tracer ${@:1}
