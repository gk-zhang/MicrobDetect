'''
@author: guokun zhang
'''

import sys
import os
import re
import argparse
from collections import defaultdict

# function: store the uniq tax_id (key) and all the possible species names (value as list) in a hash table
def store_names(names_dmp):
  l=lambda:defaultdict(l)
  table=l()
  
  with open(names_dmp) as f:
    for line in f.readlines():
      list_line = line.split("|")
      tax_id = list_line[0].strip()
      name = list_line[1].strip()
      name_class = list_line[3].strip()
      if table[tax_id][name_class]:
        table[tax_id][name_class].append(name) # if exits already
      else:
        table[tax_id][name_class] = []
        table[tax_id][name_class].append(name)
        
  # remove duplicates
  #for tax_id in table:
  #  table[tax_id] = list(set(table[tax_id]))
    
  return table

# function: store the uniq patrix FILE_NAME (key) and the NCBI_TAX_ID (value) in the patric genomes file
def store_genomes(genomes):
  l=lambda:defaultdict(l)
  table=l()
  
  with open(genomes) as f:
    next(f)
    for line in f:
      list_line = line.split("\t")
      ncbi_tax_id = list_line[3].strip()
      file_name = list_line[2].strip()
      genome_name = list_line[1].strip()
      table[file_name]["ncbi_tax_id"] = ncbi_tax_id
      table[file_name]["genome_name"] = genome_name
  
  return table
  

# function: store the tax_id (key) and its parent tax_id + its rank(value as list) in a hash table, suppose each tax_id has only one parent_tax_id
def store_nodes(nodes_dmp):
  l=lambda:defaultdict(l)
  table=l()
  
  with open(nodes_dmp) as f:
    for line in f.readlines():
      list_line = line.split("|")
      tax_id = list_line[0].strip()
      parent_tax_id = list_line[1].strip()
      rank = list_line[2].strip()
      table[tax_id]["parent"] = parent_tax_id
      table[tax_id]["rank"] = rank
   
  return table

# function: store all the tax_id in NCBI *.dmp files in a tree
#def store_tree():

# function: find the species name (scientific name) for the given tax_id, in the table of "names.dmp" file
def find_sci_name(table_names, tax_id):
  return table_names[tax_id]["scientific name"][0]

# function: find the tax_id for the given ncbi species name, in the table of "names.dmp" file
def find_tax_id(table_names, genome_name):
  for tax_id in table_names:
    for name_class in table_names[tax_id]:
      for name in table_names[tax_id][name_class]:
        if name == genome_name:
          return tax_id  

# function: find the ncbi_tax_id for the given patric file name, in the table of "genomes" file
def find_ncbi_tax_id(table_genomes, file_name):
  return table_genomes[file_name]["ncbi_tax_id"]
  
# function: find the genome_name for the given patric file name, in the table of "genomes" file
def find_ncbi_genome_name(table_genomes, file_name):
  return table_genomes[file_name]["genome_name"]

# function: find the rank(superkingdom, kingdom,...) and the parent tax_id for the given ncbi tax_id, in the table of "nodes.dmp" file
def find_parent_tax_id(table_nodes, tax_id):
  return table_nodes[tax_id]["rank"], table_nodes[tax_id]["parent"]
        
if __name__ == '__main__':
  
  #Parse Command Line
  parser = argparse.ArgumentParser()

  parser.add_argument('-patric', '--patric_list', help='a file containing all the patric bacteria entries')
  
  parser.add_argument('-names', '--names_ncbi_taxonomy', help='file names.dmp from ncbi taxonomy')

  parser.add_argument('-nodes', '--nodes_ncbi_taxonomy', help='file nodes.dmp from ncbi taxonomy')
  
  parser.add_argument('-genomes', '--genomes_patric', help='file genomes from patric relase note')
  
  args = parser.parse_args()

  list_patric = args.patric_list
  
  names_dmp = args.names_ncbi_taxonomy

  nodes_dmp = args.nodes_ncbi_taxonomy
  
  genomes = args.genomes_patric
  
  table_names = store_names(names_dmp)
  
  table_nodes = store_nodes(nodes_dmp)
  
  table_genomes = store_genomes(genomes)
  
  # for each bacteria name in patric list, find its ncbi taxonomy annotation
  with open(list_patric) as f:
    for line in f.readlines():
    
      file_name = line.strip()
      
      tax_id = find_ncbi_tax_id(table_genomes, file_name)
      
      rank, parent_tax_id = find_parent_tax_id(table_nodes, tax_id)
      
      if parent_tax_id:
        pass
      else:
        genome_name = find_ncbi_genome_name(table_genomes, file_name)
        tax_id = find_tax_id(table_names, genome_name)
        rank, parent_tax_id = find_parent_tax_id(table_nodes, tax_id)
      
      tax_annotation = rank + ":" + find_sci_name(table_names, tax_id)
      
      while tax_id != "1":
        tax_id = parent_tax_id
        rank, parent_tax_id = find_parent_tax_id(table_nodes, tax_id)
        tax_annotation = rank + ":" + find_sci_name(table_names, tax_id) + "\t" + tax_annotation
        
      print file_name + "\t" + tax_annotation
  
