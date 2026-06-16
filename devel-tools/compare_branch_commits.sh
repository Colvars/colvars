#!/usr/bin/env bash

# Iterates over the commits in the first argument (ideally: "master...branch1") and prints out
# those whose subject lines are not found in the second argument (ideally: "master...branch2")

source_branch="$1"
target_branch="$2"

mapfile -t target_branch_messages < <(git log --format=%s "${target_branch}")
mapfile -t target_branch_hashes < <(git log --format=%h "${target_branch}")


declare -A target_branch_commits
for i in "${!target_branch_hashes[@]}"; do
    target_branch_commits["${target_branch_hashes[$i]}"]="${target_branch_messages[$i]}"
done

while IFS= read -r msg; do
    if [ -z "$msg" ] ; then
        continue
    fi
    found=0
    for target_branch_hash in "${!target_branch_commits[@]}"; do
        if [ "${target_branch_commits[${target_branch_hash}]}" == "${msg}" ] ; then
            found=1
            break
        fi
    done
    if [ ${found} == 0 ] ; then
        git log --format="%h %s" "${source_branch}" | grep "${msg}"
    fi
done < <(git log --format=%s "${source_branch}")
