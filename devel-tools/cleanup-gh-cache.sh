#!/bin/bash

key=${1:-key Linux-build-}

if [ -z "${key}" ] ; then
    echo "Error: please specify a cache key as argument." >& 2
    exit 1
fi

gh actions-cache --help >& /dev/null || gh extension install actions/gh-actions-cache

# Keep at least one cache for each branch
branches=$(gh actions-cache list --key ${key} | cut -f 3 | uniq)
for branch in ${branches} ; do
    caches=($(gh actions-cache list --key ${key} --branch ${branch} | cut -f 1))
    i=${#caches[@]}
    while [ ${i} -ge 1 ] ; do
        if [ -n "${caches[${i}]}" ] ; then
            gh actions-cache delete --confirm "${caches[${i}]}" || true
        fi
        i=$((${i} - 1))
    done
done
