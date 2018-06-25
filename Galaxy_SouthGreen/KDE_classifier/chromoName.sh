#!/bin/bash
grep -v "#" $1 | cut -f 1 | uniq
