#!/usr/bin/env python3

import sys
import json
import datetime
from dateutil import parser as date_parser
import subprocess


backends = set(['GROMACS', 'LAMMPS', 'NAMD', 'TinkerHP', 'VMD'])


def affects_backend(labels, backend=None):
    if backend is None:
        # Without a backend specified, this PR affects everything
        return True
    if not backend in backends:
        raise Exception("Invalid backend:", backend)
    if 'testing' in labels:
        return False
    has_backend_labels = False
    for label in labels:
        if label in backends:
            has_backend_labels = True
    if (not has_backend_labels) or (backend in labels):
        return True
    return False


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
    authors = []
    for commit in commits:
        authors += [("@"+author['login']) for author in commit['authors']]
    return sorted(list(set(authors)), key=str.casefold)


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


def print_pr_report(kwargs):

    ref_date = kwargs.get('since')

    if not ref_date is None:
        print()
        print("The following is a list of pull requests since",
              ref_date+":")
        ref_date_ts = date_parser.parse(ref_date).timestamp()
    else:
        ref_date_ts = 0
        print()
        print("The following is a list of all pull requests:")

    pr_db = get_pr_list(kwargs.get('state'))
    all_authors = []
    for pr in pr_db:
        pr['mergedAt'] = date_parser.parse(pr['mergedAt']).timestamp()
        pr_labels = [label['name'] for label in pr['labels']]
        if pr['mergedAt'] > ref_date_ts and affects_backend(
                pr_labels, kwargs.get('backend')):
            pr_authors = get_pr_authors(pr)
            all_authors += pr_authors
            print()
            print("-", pr['number'], pr['title'])
            print(" ", pr['url'], "("+", ".join(pr_authors)+")")

    print()
    print("Authors:", ", ".join(sorted(list(set(all_authors)),
                                       key=str.casefold)))


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(
        description='''List of Colvars PRs with the given constraints''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--since',
                        type=str,
                        help="List PRs merged this date; default is all")
    parser.add_argument('--backend',
                        type=str,
                        choices=backends,
                        help="List PRs that affect only this backend; "
                        "default is all")
    parser.add_argument('--state',
                        type=str,
                        default='merged',
                        choices=['open', 'closed', 'merged', 'all'],
                        help="List PRs in this state")
    kwargs = vars(parser.parse_args())

    print_pr_report(kwargs)
