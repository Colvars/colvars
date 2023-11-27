#!/bin/bash

key=${1:-key Linux-build-}

if [ -z "${key}" ] ; then
    echo "Error: please specify a cache key as argument." >& 2
    exit 1
fi

gh actions-cache --help > /dev/null || gh extension install actions/gh-actions-cache

# Delete all caches except those generated from master
for cache in $(gh actions-cache list --key ${key} | cut -f 1); do
    branch=$(gh actions-cache list --key ${cache} | cut -f 3)
    if [ ${branch} != 'refs/heads/master' ] ; then
        # Ignore error if another job just deleted the same cache
        gh actions-cache delete --confirm "${cache}" || true
    fi
done

# Now keep one cache for master
caches=($(gh actions-cache list --key ${key} | cut -f 1))
i=${#caches[@]}
while [ ${i} -ge 1 ] ; do
    if [ -n "${caches[${i}]}" ] ; then
        gh actions-cache delete --confirm "${caches[${i}]}" || true
    fi
    i=$((${i} - 1))
done
