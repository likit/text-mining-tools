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

def to_html(data, terms):
    '''Converts data to HTML format'''
    html = '''<html><title>Search Results</title><body><table>'''
    records = Medline.parse(data)
    font_format = '<font color="red"><b>%s</b></font>'
    for rec in records:
        for word in terms.split(' '):

            # ignore these words
            if (word == 'and' or word == 'or' or word == 'not'):
                continue

            rec['TI'] = rec['TI'].replace(word, font_format % word)
            rec['TI'] = rec['TI'].replace(word.title(),
                    font_format % word.title())
            rec['TI'] = rec['TI'].replace(word.upper(),
                    font_format % word.upper())

            if 'AB' in rec.keys():
                rec['AB'] = rec['AB'].replace(word, font_format % word)
                rec['AB'] = rec['AB'].replace(word.title(),
                        font_format % word.title())
                rec['AB'] = rec['AB'].replace(word.upper(),
                        font_format % word.upper())

        for key, value in rec.iteritems():
            if key in ('TI', 'AU', 'AB', 'PMID', 'DP'):
                html += '<tr><td>%s</td><td>%s</td></tr>' % (key, value)
    html += '''</table></body></html>'''
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
        answer = raw_input('Download all articles? [Y/n]: ')

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
            # out_handle = open('%s.txt' % '_'.join(terms.split(' AND ')), 'w')
            out_html_handle = open('results.html', 'w')

            if count > MAX_ARTICLES:
                count = MAX_ARTICLES  # limit the number of articles

            for start in range(0, count, batch_size):
                end = min(count, start + batch_size)
                print('Downloading record %i to %i' %
                        (start + 1, end))
                fetch_handle = Entrez.efetch(db='pubmed', rettype='medline',
                                            retmod='text', retstart=start,
                                            retmax=batch_size, webenv=webenv,
                                            query_key=query_key)

                html = to_html(fetch_handle, search_term)
                # out_handle.write(data)
                out_html_handle.write(html)
                print('Sleeping...')
                time.sleep(5)

            # out_handle.close()
            quit = raw_input('Quit? [y/N]: ')
            if quit.lower() == 'y': break
        else:
            continue


if __name__=='__main__':
    main()
