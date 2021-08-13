# Extract BibTeX entries

import re


def gen_cplusplus_entry(entry, key, url):
    key_year_sep = re.search('\d', key).start()
    key_author = key[0:key_year_sep]
    key_year = key[key_year_sep:]
    result = """
  paper_count_[std::string("%s")] = 0;
  paper_url_[std::string("%s")] = "%s";""" % (key, key, url)
    result += """
  paper_bibtex_[std::string("%s")] =
    \"\\n\"""" % key
    for line in entry.splitlines():
        result += """
    \"%s\\n\"""" % line
    result += ';'
    return result


features = {}


def print_cplusplus_features(features):
    for key in features.keys():
        if features[key] is None:
            features[key] = 'n/a'
        print("""
  feature_count_[std::string("%s")] = 0;
  feature_paper_map_[std::string("%s")] = "%s";""" % (key, key, features[key]))


def assign_undef_features(entry):
    for key in features.keys():
        if features[key] is None:
            features[key] = entry


with open('colvars-code-refs.bib') as bib:
    entry = ""
    for line in bib.readlines():
        first_char = None
        if len(line.lstrip()) > 0:
            first_char = line.lstrip()[0]
        if first_char == '%':
            if not line.lstrip()[2:5] == '---':
                feature = line.lstrip('%').lstrip().rstrip()
                features[feature] = None
            continue
        if first_char == None:
            continue
        if first_char == '@':
            words = re.split('@|{|}',
                             line.lstrip(first_char).rstrip('\n').rstrip(','))
            key = words[1]
        else:
            if line.lstrip().lower()[0:3] == 'url':
                words = re.split(' |=|{|}', line.lstrip())
                for word in words:
                    if len(word) > 5:
                        if word[0:5] == 'https':
                            url = word
        entry += line.replace('\\', '\\\\')
        if first_char == '}':
            # Closing an entry
            print(gen_cplusplus_entry(entry, key, url))
            assign_undef_features(key)
            entry = ""
            url = ""

    print("""
  paper_url_[std::string("n/a")] = "";
  paper_bibtex_[std::string("n/a")] = "";""")

    print_cplusplus_features(features)
