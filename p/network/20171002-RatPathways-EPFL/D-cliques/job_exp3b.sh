#!/bin/bash

source ~/.bashrc

F=$(ls exp3b*.py); python3 $F 1> $F-out.log 2> $F-err.log
