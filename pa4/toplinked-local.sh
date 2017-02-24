#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <xml>"
    exit
fi

sbt "run-main TopLinked local $1"
