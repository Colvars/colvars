#!/usr/bin/env python3

import sys
import json
import datetime
from dateutil import parser as date_parser
import subprocess



def get_pr_list(state='merged'):
    # 10,000 sounds like a reasonable limit for the Colvars repo
    cmd = "gh pr list --state %s --limit 10000 --json number,url,mergedAt,title,author,labels" % state
    try:
        txt = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT).decode('UTF-8')
    except subprocess.CalledProcessError as e:
        print(e.output.decode('UTF-8'))
        raise Exception("Error calling gh command.")
    return json.loads(txt)


def get_pr_commits(number):
    cmd = "gh pr view %s --json commits" % number
    try:
        txt = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT).decode('UTF-8')
    except subprocess.CalledProcessError as e:
        print(e.output.decode('UTF-8'))
        raise Exception("Error calling gh command.")
    return json.loads(txt)['commits']


def get_commits_authors(commits):
    unique_authors = []
    for commit in commits:
        for author in commit['authors']:
            if not ("@"+author['login']) in unique_authors:
                unique_authors += ["@"+author['login']]
    return sorted(unique_authors, key=str.casefold)


def get_pr_authors(pr):
    """
    Return a PR's author first, along with others in alphabetical order
    """
    pr_authors = ["@"+pr['author']['login']]
    pr_commits = get_pr_commits(pr['number'])
    commit_authors = get_commits_authors(pr_commits)
    for commit_author in commit_authors:
        if not commit_author in pr_authors:
            pr_authors += [commit_author]
    return pr_authors


def print_pr_report(ref_date=None):

    if not ref_date is None:
        print()
        print("The following is a list of pull requests since",
              ref_date+":")
        ref_date_ts = date_parser.parse(ref_date).timestamp()
    else:
        ref_date_ts = 0
        print()
        print("The following is a list of all pull requests:")

    pr_db = get_pr_list()
    for pr in pr_db:
        pr['mergedAt'] = date_parser.parse(pr['mergedAt']).timestamp()
        if pr['mergedAt'] > ref_date_ts:
            pr_authors = "("+", ".join(get_pr_authors(pr))+")"
            print()
            print("-", pr['number'], pr['title'])
            print(" ", pr['url'], pr_authors)


if __name__ == '__main__':

    ref_date = None
    if len(sys.argv) > 1:
        ref_date = sys.argv[1]

    print_pr_report(ref_date)
