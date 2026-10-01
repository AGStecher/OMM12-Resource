#!/usr/bin/env python3
"""Build data/omm12_publications.tsv for the About page of the OMM12 app.

Same rules as OMM12_Publications_2015-2026.ipynb (statistics_OMM12 folder):
  * all Dimensions-Publication-*.csv exports in SOURCE_DIR are merged,
    de-duplicated by Publication ID (newest export wins), years 2015-2026
  * bioRxiv-only records and off-topic hits are dropped
  * only research articles, reviews and book chapters are kept
  * 'uses OMM12' = usage category starting with "used" (manual_decision
    overrides omm12_usage)
  * citations from omm12_citations.csv (Scopus, Web of Science, Dimensions)
  * countries from the authors' raw affiliations (country-matching code
    copied from the notebook)

Output (one row per publication):
  data/omm12_publications.tsv
  data/omm12_publications_info.tsv   (export date, citation date, counts)

Usage (from the OMM12_website folder):
  python build_publication_stats/build_publication_stats.py
  python build_publication_stats/build_publication_stats.py <folder with the exports>
"""
import glob, re, sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ROOT       = Path(__file__).resolve().parent.parent
SOURCE_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'data' / 'statistics_OMM12' / 'figures_v7' / 'OMM12_worldmap'
OUT_DIR    = ROOT / 'data'
YEARS      = list(range(2015, 2027))
EXCLUDE_SOURCES   = ['bioRxiv']
USAGE_FILE        = SOURCE_DIR / 'omm12_usage_check.csv'
CITATIONS_FILE    = SOURCE_DIR / 'omm12_citations.csv'
CITE_DATE         = '2026-09-28'   # date the Scopus/WoS counts were retrieved -- update when re-fetched
OFFTOPIC_CATEGORY = 'off-topic (different meaning of OMM / unrelated field)'
INCLUDE_TYPES     = ['Research article', 'Review', 'Book chapter']

# ISO3 codes for the interactive world map
ISO3 = {
    'Germany': 'DEU', 'United States': 'USA', 'United Kingdom': 'GBR', 'Switzerland': 'CHE', 'France': 'FRA',
    'China': 'CHN', 'Hong Kong': 'HKG', 'Taiwan': 'TWN', 'Canada': 'CAN', 'Japan': 'JPN', 'Italy': 'ITA',
    'South Korea': 'KOR', 'Denmark': 'DNK', 'Netherlands': 'NLD', 'Belgium': 'BEL', 'Luxembourg': 'LUX',
    'Austria': 'AUT', 'Czechia': 'CZE', 'Slovakia': 'SVK', 'Poland': 'POL', 'Hungary': 'HUN',
    'Slovenia': 'SVN', 'Croatia': 'HRV', 'Serbia': 'SRB', 'Romania': 'ROU', 'Bulgaria': 'BGR',
    'Greece': 'GRC', 'Cyprus': 'CYP', 'Spain': 'ESP', 'Portugal': 'PRT', 'Ireland': 'IRL', 'Sweden': 'SWE',
    'Norway': 'NOR', 'Finland': 'FIN', 'Iceland': 'ISL', 'Estonia': 'EST', 'Latvia': 'LVA',
    'Lithuania': 'LTU', 'Ukraine': 'UKR', 'Belarus': 'BLR', 'Russia': 'RUS', 'Turkey': 'TUR',
    'Israel': 'ISR', 'Iran': 'IRN', 'Iraq': 'IRQ', 'Jordan': 'JOR', 'Lebanon': 'LBN', 'Saudi Arabia': 'SAU',
    'United Arab Emirates': 'ARE', 'Qatar': 'QAT', 'Kazakhstan': 'KAZ', 'Uzbekistan': 'UZB',
    'Tajikistan': 'TJK', 'Egypt': 'EGY', 'Morocco': 'MAR', 'Algeria': 'DZA', 'Tunisia': 'TUN',
    'Nigeria': 'NGA', 'Ghana': 'GHA', 'Cameroon': 'CMR', 'Kenya': 'KEN', 'Uganda': 'UGA', 'Tanzania': 'TZA',
    'Ethiopia': 'ETH', 'South Africa': 'ZAF', 'India': 'IND', 'Pakistan': 'PAK', 'Bangladesh': 'BGD',
    'Sri Lanka': 'LKA', 'Nepal': 'NPL', 'Thailand': 'THA', 'Vietnam': 'VNM', 'Malaysia': 'MYS',
    'Singapore': 'SGP', 'Indonesia': 'IDN', 'Philippines': 'PHL', 'Australia': 'AUS', 'New Zealand': 'NZL',
    'Mexico': 'MEX', 'Costa Rica': 'CRI', 'Panama': 'PAN', 'Cuba': 'CUB', 'Brazil': 'BRA',
    'Argentina': 'ARG', 'Chile': 'CHL', 'Colombia': 'COL', 'Peru': 'PER', 'Ecuador': 'ECU',
    'Uruguay': 'URY', 'Venezuela': 'VEN',
}


def read_dimensions(path):
    with open(path, encoding='utf-8-sig') as f:
        first = f.readline()
    df = pd.read_csv(path, skiprows=1 if 'About the data' in first else 0, encoding='utf-8-sig', dtype={'PMID': str})
    df['export_file'] = Path(path).name
    return df


files = sorted(glob.glob(str(SOURCE_DIR / 'Dimensions-Publication-*.csv')), reverse=True)
if not files:
    sys.exit(f'No Dimensions-Publication-*.csv files in {SOURCE_DIR}')
raw = pd.concat([read_dimensions(f) for f in files], ignore_index=True)
n_rows = len(raw)
pubs = (raw.drop_duplicates('Publication ID', keep='first')
           .query('PubYear >= @YEARS[0] and PubYear <= @YEARS[-1]').reset_index(drop=True))
pubs['PubYear'] = pubs['PubYear'].astype(int)
n_unique = len(pubs)
pubs = pubs[~pubs['Source title'].fillna('').str.strip().str.casefold()
              .isin([s.casefold() for s in EXCLUDE_SOURCES])].reset_index(drop=True)

usage = pd.read_csv(USAGE_FILE, encoding='utf-8-sig', dtype=str)
usage['category'] = usage['manual_decision'].where(usage['manual_decision'].fillna('').str.strip() != '',
                                                   usage['omm12_usage'])
usage = usage.drop_duplicates('Publication ID').set_index('Publication ID')
off = set(usage.index[usage['category'] == OFFTOPIC_CATEGORY])
pubs = pubs[~pubs['Publication ID'].isin(off)].reset_index(drop=True)

fallback_type = pubs['Publication Type'].map({'Chapter': 'Book chapter', 'Edited Book': 'Book',
                                              'Proceeding': 'Conference abstract',
                                              'Preprint': 'Preprint'}).fillna('Research article')
pubs['article_type'] = pubs['Publication ID'].map(usage['article_type']).fillna(fallback_type)
pubs = pubs[pubs['article_type'].isin(INCLUDE_TYPES)].reset_index(drop=True)
pubs['usage_category'] = pubs['Publication ID'].map(usage['category']).fillna('not classified')
pubs['uses_omm12'] = pubs['usage_category'].str.startswith('used')

cites = pd.read_csv(CITATIONS_FILE)
cites['doi'] = cites['doi'].str.strip().str.lower()
cites = cites.drop_duplicates('doi').set_index('doi')
doi_key = pubs['DOI'].fillna('').str.strip().str.lower()
for col in ['scopus_cited_by', 'wos_times_cited', 'dimensions_times_cited']:
    pubs[col] = doi_key.map(cites[col]) if col in cites else pd.NA

# journal display name: most common spelling per journal (as in the notebook)
key = pubs['Source title'].fillna('').str.strip().str.casefold()
display = pubs.assign(k=key).groupby('k')['Source title'].agg(
    lambda s: s.dropna().str.strip().value_counts().index[0] if s.notna().any() else '')
pubs['journal'] = key.map(display)

# ---- country matching: copied from OMM12_Publications_2015-2026.ipynb ----
COUNTRY_PATTERNS = {
    'Germany': [r'Germany', r'Deutschland'],
    'United States': [r'USA', r'U\.S\.A', r'United States', r'\bUS(?=\s*$)'],
    'United Kingdom': [r'United Kingdom', r'\bUK\b', r'England', r'Scotland', r'(?<!New South )Wales',
                       r'Northern Ireland', r'Great Britain'],
    'Switzerland': [r'Switzerland', r'Schweiz', r'Suisse'],
    'France': [r'France'], 'China': [r'China', r'P\.\s?R\.\s?C', r'中国', r'广州', r'北京', r'上海'],
    'Hong Kong': [r'Hong Kong'], 'Taiwan': [r'Taiwan'], 'Canada': [r'Canada'], 'Japan': [r'Japan'],
    'Italy': [r'Italy', r'Italia'], 'South Korea': [r'South Korea', r'Republic of Korea', r'\bKorea\b'],
    'Denmark': [r'Denmark'], 'Netherlands': [r'Netherlands', r'Nederland'], 'Belgium': [r'Belgium'],
    'Luxembourg': [r'Luxembourg'], 'Austria': [r'Austria', r'Österreich'],
    'Czechia': [r'Czech Republic', r'Czechia'], 'Slovakia': [r'Slovakia'], 'Poland': [r'Poland'],
    'Hungary': [r'Hungary'], 'Slovenia': [r'Slovenia'], 'Croatia': [r'Croatia'], 'Serbia': [r'Serbia'],
    'Romania': [r'Romania'], 'Bulgaria': [r'Bulgaria'], 'Greece': [r'Greece'], 'Cyprus': [r'Cyprus'],
    'Spain': [r'Spain', r'España'], 'Portugal': [r'Portugal'], 'Ireland': [r'(?<!Northern )Ireland'],
    'Sweden': [r'Sweden'], 'Norway': [r'Norway'], 'Finland': [r'Finland'], 'Iceland': [r'Iceland'],
    'Estonia': [r'Estonia'], 'Latvia': [r'Latvia'], 'Lithuania': [r'Lithuania'], 'Ukraine': [r'Ukraine'],
    'Belarus': [r'Belarus'], 'Russia': [r'Russia'], 'Turkey': [r'Turkey', r'Türkiye', r'Turkiye'],
    'Israel': [r'Israel'], 'Iran': [r'\bIran\b'], 'Iraq': [r'\bIraq\b'], 'Jordan': [r'\bJordan\b'],
    'Lebanon': [r'Lebanon'], 'Saudi Arabia': [r'Saudi Arabia'],
    'United Arab Emirates': [r'United Arab Emirates', r'\bUAE\b'], 'Qatar': [r'Qatar'],
    'Kazakhstan': [r'Kazakhstan'], 'Uzbekistan': [r'Uzbekistan'], 'Tajikistan': [r'Tajikistan'],
    'Egypt': [r'Egypt'], 'Morocco': [r'Morocco'], 'Algeria': [r'Algeria'], 'Tunisia': [r'Tunisia'],
    'Nigeria': [r'Nigeria'], 'Ghana': [r'Ghana'], 'Cameroon': [r'Cameroon'], 'Kenya': [r'Kenya'],
    'Uganda': [r'Uganda'], 'Tanzania': [r'Tanzania'], 'Ethiopia': [r'Ethiopia'],
    'South Africa': [r'South Africa'], 'India': [r'\bIndia\b'], 'Pakistan': [r'Pakistan'],
    'Bangladesh': [r'Bangladesh'], 'Sri Lanka': [r'Sri Lanka'], 'Nepal': [r'Nepal'],
    'Thailand': [r'Thailand'], 'Vietnam': [r'Vietnam', r'Viet Nam'], 'Malaysia': [r'Malaysia'],
    'Singapore': [r'Singapore'], 'Indonesia': [r'Indonesia'], 'Philippines': [r'Philippines'],
    'Australia': [r'Australia'], 'New Zealand': [r'New Zealand'],
    'Mexico': [r'(?<!New )Mexico', r'México'], 'Costa Rica': [r'Costa Rica'], 'Panama': [r'Panama'],
    'Cuba': [r'\bCuba\b'], 'Brazil': [r'Brazil', r'Brasil'], 'Argentina': [r'Argentina'],
    'Chile': [r'\bChile\b'], 'Colombia': [r'Colombia'], 'Peru': [r'\bPeru\b'], 'Ecuador': [r'Ecuador'],
    'Uruguay': [r'Uruguay'], 'Venezuela': [r'Venezuela'],
}
# Last-resort city / institution hints for affiliations that give no country
CITY_HINTS = {
    'Germany': ['Munich', 'München', 'Berlin', 'Hamburg', 'Heidelberg', 'Aachen', 'Freising', 'Tübingen',
                'Hannover', 'Braunschweig', 'Kiel', 'Cologne', 'Köln', 'Frankfurt', 'Mainz', 'Würzburg',
                'Göttingen', 'Jena', 'Leipzig', 'Dresden', 'Freiburg', 'Ulm', 'Regensburg', 'Bonn', 'LMU',
                'Helmholtz', 'Max Planck', 'DKTK', 'DZIF'],
    'Switzerland': ['Zurich', 'Zürich', 'Basel', 'Bern', 'Lausanne', 'Geneva', 'Genève', 'ETH', 'EPFL'],
    'France': ['Paris', 'Lyon', 'Marseille', 'Toulouse', 'Strasbourg', 'Montpellier', 'Lille', 'INRAE',
               'INSERM', 'CNRS', 'Institut Pasteur'],
    'Austria': ['Vienna', 'Wien', 'Graz', 'Innsbruck', 'Salzburg'],
    'United Kingdom': ['London', 'Manchester', 'Oxford', 'Cambridge', 'Edinburgh', 'Glasgow', 'Birmingham',
                       'Norwich', 'Quadram', 'Wellcome Sanger'],
    'Spain': ['Madrid', 'Barcelona', 'Valencia', 'Sevilla', 'Seville', 'CSIC', 'Complutense', 'Pompeu Fabra',
              'ICREA'],
    'Italy': ['Rome', 'Roma', 'Milan', 'Milano', 'Padova', 'Padua', 'Torino', 'Turin', 'Bologna', 'Naples',
              'Napoli', 'Catania', 'Florence', 'Firenze', 'IRCCS'],
    'Canada': ['Calgary', 'Toronto', 'Montreal', 'Montréal', 'Vancouver', 'Edmonton', 'Ottawa', 'Hamilton',
               'McMaster', 'Alberta', 'Ontario', 'Quebec', 'Québec', 'British Columbia'],
    'China': ['Guangzhou', 'Shenzhen', 'Beijing', 'Shanghai', 'Wuhan', 'Hangzhou', 'Nanjing', 'Chengdu',
              'Jinan', 'Anhui', 'Guangdong', 'Zhejiang', 'Chinese Academy'],
    'Sweden': ['Umeå', 'Stockholm', 'Gothenburg', 'Göteborg', 'Uppsala', 'Lund', 'Karolinska'],
    'Denmark': ['Copenhagen', 'Aarhus'],
    'Netherlands': ['Amsterdam', 'Utrecht', 'Leiden', 'Rotterdam', 'Groningen', 'Wageningen'],
    'Belgium': ['Leuven', 'Ghent', 'Brussels'],
    'Czechia': ['Prague', 'Praha', 'Brno', 'Nový Hrádek'],
    'Ethiopia': ['Addis Ababa'],
    'Japan': ['Tokyo', 'Osaka', 'Kyoto', 'RIKEN'],
    'United States': ['Harvard', 'Stanford', 'Yale', 'MIT', 'NIH', 'Bethesda', 'Boston'],
}
US_STATES = ['Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut',
             'Delaware', 'Florida', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas',
             'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota',
             'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey',
             'New Mexico', 'New York', 'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon',
             'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas',
             'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming']
US_ZIP = (r'\b(?:A[KLRZ]|C[AOT]|D[CE]|FL|GA|HI|I[ADLN]|K[SY]|LA|M[ADEINOST]|N[CDEHJMVY]|O[HKR]|PA|RI|'
          r'S[CD]|T[NX]|UT|V[AT]|W[AIVY])\s+\d{5}\b')

def _rx(p):
    return p if p.startswith(('\\b', '(?')) or not p.isascii() else rf'\b(?:{p})\b'
COUNTRY_RX = {c: re.compile('|'.join(_rx(p) for p in pats)) for c, pats in COUNTRY_PATTERNS.items()}
CITY_RX    = {c: re.compile('|'.join(_rx(re.escape(p)) for p in pats)) for c, pats in CITY_HINTS.items()}
US_RX      = re.compile(r'\b(?:' + '|'.join(US_STATES) + r')\b|' + US_ZIP)

def split_affiliations(raw_text):
    # 'Name (aff1; aff2); Name2 (aff3)' -> ['aff1', 'aff2', 'aff3']  (handles nested brackets)
    chunks, depth, buf = [], 0, []
    for ch in raw_text:
        if ch == '(':
            if depth: buf.append(ch)
            depth += 1
        elif ch == ')':
            depth = max(depth - 1, 0)
            if depth: buf.append(ch)
            else: chunks.append(''.join(buf)); buf = []
        elif depth:
            buf.append(ch)
    return [a.strip() for c in chunks for a in c.split(';') if a.strip() and '@' not in a]

def match_country(aff):
    found = {c for c, rx in COUNTRY_RX.items() if rx.search(aff)}
    if 'Hong Kong' in found:
        found.discard('China')
    if not found and US_RX.search(aff):
        found = {'United States'}
    return found

# pass 1: explicit country names
pub_affs = pubs['Authors (Raw Affiliation)'].fillna('').map(split_affiliations)
matched  = pub_affs.map(lambda affs: [match_country(a) for a in affs])

# pass 2: learn institution segments -> country from affiliations that did resolve
seg_votes = defaultdict(Counter)
for affs, ms in zip(pub_affs, matched):
    for a, m in zip(affs, ms):
        if len(m) == 1:
            for seg in (s.strip() for s in a.split(',')):
                if len(seg) > 6:
                    seg_votes[seg][next(iter(m))] += 1

def fallback(aff):
    votes = Counter()
    for seg in (s.strip() for s in aff.split(',')):
        votes.update(seg_votes.get(seg, {}))
    if votes:
        return {votes.most_common(1)[0][0]}
    return {c for c, rx in CITY_RX.items() if rx.search(aff)}


pubs['countries'] = [
    sorted(set().union(*[m if m else fallback(a) for a, m in zip(affs, ms)])) if affs else []
    for affs, ms in zip(pub_affs, matched)]
# ---- end of copied code ----

unknown = sorted({c for cs in pubs['countries'] for c in cs if c not in ISO3})
if unknown:
    print('WARNING: no ISO3 code for', unknown, '-- add them to ISO3; they will be missing from the map')

out = pd.DataFrame({
    'publication_id': pubs['Publication ID'],
    'year': pubs['PubYear'],
    'title': pubs['Title'],
    'journal': pubs['journal'],
    'publication_type': pubs['Publication Type'],
    'article_type': pubs['article_type'],
    'usage_category': pubs['usage_category'],
    'uses_omm12': pubs['uses_omm12'].map({True: 'TRUE', False: 'FALSE'}),
    'doi': pubs['DOI'],
    'pmid': pubs['PMID'],
    'dimensions_url': pubs['Dimensions URL'],
    'open_access': pubs['Open Access'],
    'fcr': pubs['FCR'],
    'citations_scopus': pubs['scopus_cited_by'],
    'citations_wos': pubs['wos_times_cited'],
    'citations_dimensions': pubs['dimensions_times_cited'],
    'countries': pubs['countries'].map(lambda cs: '; '.join(cs)),
    'countries_iso3': pubs['countries'].map(lambda cs: '; '.join(ISO3[c] for c in cs if c in ISO3)),
}).sort_values(['year', 'title'], ascending=[False, True])
for c in ['title', 'journal']:
    out[c] = out[c].fillna('').str.replace(r'[\t\r\n]+', ' ', regex=True)
out.to_csv(OUT_DIR / 'omm12_publications.tsv', sep='\t', index=False)

export_date = re.search(r'(\d{4}-\d{2}-\d{2})', Path(files[0]).name).group(1)
info = pd.DataFrame([
    ('dimensions_export_date', export_date),
    ('dimensions_exports', ', '.join(Path(f).name for f in files)),
    ('citations_date', CITE_DATE),
    ('rows_in_exports', n_rows),
    ('unique_publications', n_unique),
    ('included_publications', len(out)),
    ('uses_omm12', int(pubs['uses_omm12'].sum())),
], columns=['key', 'value'])
info.to_csv(OUT_DIR / 'omm12_publications_info.tsv', sep='\t', index=False)
print(f'{len(files)} export(s), {n_rows} rows -> {n_unique} unique -> {len(out)} included, '
      f'{int(pubs["uses_omm12"].sum())} use OMM12')
print(f'wrote {OUT_DIR / "omm12_publications.tsv"} and omm12_publications_info.tsv')

# Aggregated tables for the app (the per-publication list above stays local)
import runpy
runpy.run_path(str(Path(__file__).with_name('make_publication_summary.py')), run_name='__main__')
