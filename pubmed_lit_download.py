#!/usr/bin/env python
'''Description'''
import sys

try:
    import inflect
except ImportError:
    print('Install Python Inflect module by running the following command:')
    print('pip install -e git+https://github.com/benthor/inflect.py#egg=inflect')
    raise SystemExit

import time
from Bio import Entrez, Medline

p = inflect.engine()
MAX_ARTICLES = 10

def get_search_count(terms):
    handle = Entrez.egquery(term=terms)
    record = Entrez.read(handle)
    for row in record['eGQueryResult']:
        if row['DbName'] == 'pubmed':
            count = int(row['Count'])
            print('Total items found %s ' % count)
    return count

def check_plural(terms):
    new_terms = terms
    for word in terms.replace(')','').replace('(','').split(' '):
        # ignore these words
        if (word == 'and' or word == 'or' or word == 'not'):
            continue

        if (not word.endswith('s') and p.plural(word) != word + 's'):
            print('Warning: (correct) plural form of %s is %s' %
                                                (word, p.plural(word)))
            while True:
                answer = raw_input('Include/Replace/Ignore [i/r/g]?: ')
                if answer == 'r':
                    new_terms = new_terms.replace(word, p.plural(word))
                    break
                elif answer == 'i':
                    new_terms = new_terms.replace(word, '(%s or %s)' %
                                                    (word, p.plural(word)))
                    break
                elif answer == 'g':
                    break
                else:
                    print('Unrecognized answer!')
                    continue

    return new_terms


def add_search_type(terms, stype):
    new_terms = terms
    for word in terms.replace(')','').replace('(','').split(' '):
        # ignore these words
        if (word == 'and' or word == 'or' or word == 'not'):
            continue

        new_terms = new_terms.replace(word, word + stype, 1)
    return new_terms


def get_html_head(terms):
    html = '<!DOCTYPE html>'
    html += '\n<html>'
    html += '\n<head>'
    html += '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">'
    html += '<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/js/bootstrap.min.js"></script>'
    html += '\n</head>'
    html += '\n<title>Search Results</title>\n<body>'
    html += '\n<h1>Search Result</h1>'
    html += '\n<table class="table table-hover">'
    html += '\n<tr><td>%s</td><td>%s</td></tr>' % ('Search terms', terms)
    return html


def get_html_tail():
    html = '\n</table></body></html>'
    return html

def to_html(data, terms):
    '''Converts data to HTML format'''

    html = ''
    records = Medline.parse(data)
    font_format = '<font color="red"><strong>%s</strong></font>'
    for rec in records:
        for word in terms.split(' '):

            # ignore these words
            if (word == 'and' or word == 'or' or word == 'not'):
                continue

            rec['TI'] = rec['TI'].replace(word.lower(), font_format % word)
            rec['TI'] = rec['TI'].replace(word.title(),
                    font_format % word.title())
            rec['TI'] = rec['TI'].replace(word.upper(),
                    font_format % word.upper())

            if 'AB' in rec.keys():
                rec['AB'] = rec['AB'].replace(word.lower(), font_format % word)
                rec['AB'] = rec['AB'].replace(word.title(),
                        font_format % word.title())
                rec['AB'] = rec['AB'].replace(word.upper(),
                        font_format % word.upper())

        for key, value in rec.iteritems():
            if key in ('TI', 'AU', 'AB', 'PMID', 'DP'):
                if key == 'AB':
                    html += \
                        '\n<tr class="success"><td>%s</td>' \
                        '<td>%s</td></tr>' % (key, value)
                elif key == 'TI':
                    html += \
                        '\n<tr class="active"><td>%s</td>' \
                        '<td><strong>%s</strong></td></tr>' % (key, value)
                elif key == 'PMID':
                    html += \
                        '\n<tr><td>%s</td>' \
                        '<td><a href="http://www.ncbi.nlm.nih.gov/pubmed?' \
                        'term=%s%%5BPMID%%5D" target="_blank">' \
                        '%s</a></td></tr>' % (key, value, value)
                elif key == 'AU':
                    html += \
                        '\n<tr><td>%s</td>' \
                        '<td>%s</td></tr>' % (key, ', '.join(value))
                else:
                    html += '\n<tr><td>%s</td><td>%s</td></tr>' % (key, value)
    return html


def main():
    '''Main function'''
    Entrez.email = raw_input('Enter your email: ')
    while True:
        search_term = raw_input('Enter search term(s): ')
        search_term = check_plural(search_term)
        search_term_with_types = add_search_type(search_term, '[Title/Abstract]')
        print('Search term(s) = %s' % search_term_with_types)
        count = get_search_count(search_term_with_types)
        if count > 0:
            answer = raw_input('Download all articles? [Y/n]: ')
        else:
            continue

        if answer.lower() == 'y' or answer.lower() == '':
            # use history feature for searching
            search_handle = Entrez.esearch(db='pubmed',
                    term=search_term_with_types, retmax=count, usehistory='y')
            search_results = Entrez.read(search_handle)
            search_handle.close()

            idlist = search_results['IdList']
            webenv = search_results['WebEnv']
            query_key = search_results['QueryKey']

            batch_size = 5
            out_html_handle = open('results.html', 'w')

            if count > MAX_ARTICLES:
                count = MAX_ARTICLES  # limit the number of articles

            html = get_html_head(search_term_with_types)
            for start in range(0, count, batch_size):
                end = min(count, start + batch_size)
                print('Downloading record %i to %i' %
                        (start + 1, end))
                fetch_handle = Entrez.efetch(db='pubmed', rettype='medline',
                                            retmod='text', retstart=start,
                                            retmax=batch_size, webenv=webenv,
                                            query_key=query_key)

                html += to_html(fetch_handle, search_term)
                print('Sleeping...')
                time.sleep(5)

            html += get_html_tail()
            out_html_handle.write(html)
            out_html_handle.close()
            quit = raw_input('Quit? [y/N]: ')
            if quit.lower() == 'y':
                break

        else:
            continue


if __name__=='__main__':
    main()
