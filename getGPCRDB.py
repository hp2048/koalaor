import requests
import re
import os

def clean_family_name(text):
	cleanr = re.compile('<.*?>')
	cleantext = re.sub(cleanr,'',text)
	cleantext = re.sub(' ','_',cleantext)
  cleantext = re.sub(';','',cleantext)
	cleantext = re.sub('-','_',cleantext)
	cleantext = re.sub('&', '', cleantext)
	cleantext = re.sub('/','_', cleantext)
	return cleantext

# fetch children for Class A (Rhodopsin) family, slug 001
url = 'http://gpcrdb.org/services/proteinfamily/descendants/001/'
response = requests.get(url)
children_data = response.json()

for child in children_data:
	if child['slug'].count('_') == 3:
	url = 'http://gpcrdb.org/services/proteinfamily/proteins/' + child['slug'] + '/'
	response = requests.get(url)
	proteins = response.json()
	printed = 0
	family_name = clean_family_name(child['name'])
	output = open(family_name+'.fasta', 'w')
	for protein in proteins:
		if protein['source'] == 'SWISSPROT':
		printed += 1
		output.write('>'+protein['entry_name']+'\n')
		output.write(protein['sequence']+'\n')
	output.close()
	if printed < 2:
		os.remove(family_name+'.fasta')
