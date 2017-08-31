'''
@author: guokun zhang
'''

import sys
import os
import re
import argparse
from collections import defaultdict

# function: return the number of lines of the given file
def read_file(filename):
    len=0
    with open(filename) as f:
      for line in f.readlines():
        len=len+1
    return len
    
def merge_intervals(sorted_by_lower_bound):
    """
    A simple algorithm can be used:
    #1. Sort the intervals in increasing order
    2. Push the first interval on the stack
    3. Iterate through intervals and for each one compare current interval
       with the top of the stack and:
       A. If current interval does not overlap, push on to stack
       B. If current interval does overlap, merge both intervals in to one
          and push on to stack
    4. At the end return stack
    """
    #sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
      if not merged:
        merged.append(higher)
      else:
        lower = merged[-1]
        # test for intersection between lower and higher:
        # we know via sorting that lower[0] <= higher[0]
        # if we also want to merge the ajacent intervals, e.g. (1,3),(4,7)->(1,7)
        #   use: if higher[0] <= (lower[1] + 1)
        # if we do not want to merge the ajacent intevals,
        #   use: if higher[0] <= lower[1]
        if higher[0] <= (lower[1] + 1):
          upper_bound = max(lower[1], higher[1])
          merged[-1] = (lower[0], upper_bound)  # replace by merged interval
        else:
          merged.append(higher)
    return merged

def mapped_intervals(pos, cigar, md):
    """
    The algorithm is:
    if pos = 0, means not mapped to reference sequence
      return empty list
    if pos !=0, means mapped to reference sequence
      if MD tag is not empty, use MD tag to get matched intervals
      otherwise use CIGAR tag to get matched intervals (mismatches are ignored, therefore the intervals are larger than true)
    """
    if pos == 0:
      # not mapped to reference sequence
      return [];
    else:
      list_mapped_intervals = []
      current_pos = int(pos)
      
      if md != "":
        # if MD tag is not empty, use md tag only, this is the exact way
        p = re.compile('\d+\D*')
        list_md = p.findall(md)
        
        for x in xrange(0, len(list_md)):
          p_num = re.compile('\d+')
          dist_match = int(p_num.findall(list_md[x])[0])
          if "^" in list_md[x]:
            # in the case of a deletion event 
            dist_skip = len(list_md[x]) - len(str(dist_match)) - 1
          else:
            dist_skip = len(list_md[x]) - len(str(dist_match))
          
          # add into the interval
          if dist_match != 0:
            list_mapped_intervals.append((current_pos, current_pos + dist_match -1))
          current_pos = current_pos + dist_match + dist_skip
            
      else:
        # if MD tag is empty, use CIGAR tag only, in this case, the mismatches within the "M" segments are ignored
        p = re.compile('\d+\D')
        list_cigar_ope = p.findall(cigar)
    
        for x in xrange(0, len(list_cigar_ope)):
          distance = int(list_cigar_ope[x][:-1])
          operation = list_cigar_ope[x][-1]
          if x == 0:
            # for the first operation
            if operation == "M":
              list_mapped_intervals.append((current_pos, current_pos + distance -1))
              current_pos = current_pos + distance
            elif operation == "I":
              pass 
            elif operation == "D":
              pass
            elif operation == "N":
              pass
            elif operation == "S":
              pass            
            elif operation == "H":
              pass
            elif operation == "P":
              pass            
            elif operation == "=":
              pass
            elif operation == "X":         
              pass
          else:
            if operation == "M":
              list_mapped_intervals.append((current_pos, current_pos + distance -1))
              current_pos = current_pos + distance
            elif operation == "I":
              pass
            elif operation == "D":
              current_pos = current_pos + distance
            elif operation == "N":
              current_pos = current_pos + distance
            elif operation == "S":
              pass
            elif operation == "H":
              pass
            elif operation == "P":
              pass 
            elif operation == "=":
              pass
            elif operation == "X": 
              pass
      return list_mapped_intervals
    
# function: store the mapping information into the hash/tree, for the patric all mode
def store_tree_all(filename):
    l=lambda:defaultdict(l)
    table=l()
    
    # a list to save the length of each reference sequence
    list_len=l()
    
    with open(filename) as f:
      for line in f.readlines():
        if line.startswith('@'):
          # if it is the head information line
          # create all the hash tables
          if line.startswith('@SQ'):
            info = line.split("\t")[1].split("|")
            spe_name = info[len(info)-1]
            contig_name = info[len(info)-2]
            table[spe_name][contig_name] = []
            list_len[spe_name][contig_name] = int(line.split("\t")[2].strip()[3:])
        else:
          # if it is the main content
          # fill in the hash tables
          
          # do it only if mapped
          if line.split("\t")[2].strip() != "*":
            info = line.split("\t")[2].split("|")
            spe_name = info[len(info)-1]
            contig_name = info[len(info)-2]          
          
            # the leftmost mapping position. if set as "0", means unmapped
            pos = line.split("\t")[3]
          
            # the cigar field
            cigar = line.split("\t")[5]
          
            # the MD tag list if available, otherwise leave empty
            md = ""
            for x in xrange(11, len(line.split("\t"))):
              p = re.compile('MD:Z.*')
              if len(p.findall(line.split("\t")[x].strip())) >= 1:
                md = line.split("\t")[x].strip()[5:]
  
            # based on the information of mapping postion and cigar tag, return a list of intervals which truelly matched back to the reference sequence
            list_mapped_intervals = mapped_intervals(pos, cigar, md)
          
            # add each interval in the list into the hash table
            for each_interval in list_mapped_intervals:
              table[spe_name][contig_name].append(each_interval)
            
    return table, list_len
    # for spe in table:
      # for contig in table[spe]:
        # # sort the list
        # # table[spe][contig].sort()
        # # remove duplicates from the list
        # #print spe,contig,sorted(list(set(table[spe][contig])))
        # print spe, contig, merge_intervals(table[spe][contig])

# functiong:store the mapping information into the hash/tree, for the single reference mode         
def store_tree_single(filename):        
    l=lambda:defaultdict(l)
    table=l()
    
    # a list to save the length of each reference sequence
    list_len=l()
    
    with open(filename) as f:
      for line in f.readlines():
        if line.startswith('@'):
          # if it is the head information line
          # create all the hash tables
          if line.startswith('@SQ'):
            contig_name = line.split("\t")[1].strip()[3:]
            table[contig_name] = []
            list_len[contig_name] = int(line.split("\t")[2].strip()[3:])
        else:
          # if it is the main content
          # fill in the hash tables
          
          # do it only if mapped
          if line.split("\t")[2].strip() != "*":
            contig_name = line.split("\t")[2].strip()          
          
            # the leftmost mapping position. if set as "0", means unmapped
            pos = line.split("\t")[3]
          
            # the cigar field
            cigar = line.split("\t")[5]
          
            # the MD tag list if available, otherwise leave empty
            md = ""
            for x in xrange(11, len(line.split("\t"))):
              p = re.compile('MD:Z.*')
              if len(p.findall(line.split("\t")[x].strip())) >= 1:
                md = line.split("\t")[x].strip()[5:]
  
            # based on the information of mapping postion and cigar tag, return a list of intervals which truelly matched back to the reference sequence
            list_mapped_intervals = mapped_intervals(pos, cigar, md)
          
            # add each interval in the list into the hash table
            for each_interval in list_mapped_intervals:
              table[contig_name].append(each_interval)
            
    return table, list_len
    
    
# function: calculate for each species / contigs the covered length/percentage of intervals, for patric all mode
def calculate_all(table, list_len, output): 
    
    output_all = open(output+"_all.txt", 'w')
    output_merged = open(output+"_merged.txt", 'w')
    output_contig = open(output+"_contig.txt", 'w')
    output_spe = open(output+"_spe.txt", 'w')
   
    for spe in table:
      len_covered_spe = 0
      len_ref_spe = 0
      
      for contig in table[spe]:
        len_covered_contig = 0
        sorted_by_lower_bound = sorted(table[spe][contig], key=lambda tup: tup[0])
        for each_inter in sorted_by_lower_bound:
          output_all.write(spe)
          output_all.write("\t")
          output_all.write(contig)
          output_all.write("\t")
          output_all.write(str(each_inter[0]))
          output_all.write("\t")
          output_all.write(str(each_inter[1]))
          output_all.write("\n")
 
        intervals = merge_intervals(sorted_by_lower_bound)
        
        for interval in intervals:
          len_covered_contig = len_covered_contig + abs(interval[1] - interval[0] + 1)
          output_merged.write(spe)
          output_merged.write("\t")
          output_merged.write(contig)
          output_merged.write("\t")
          output_merged.write(str(list_len[spe][contig]))
          output_merged.write("\t")
          output_merged.write(str(interval[0]))
          output_merged.write("\t")
          output_merged.write(str(interval[1]))
          output_merged.write("\t")
          output_merged.write(str(abs(interval[1]-interval[0]+1)))
          output_merged.write("\t")
          output_merged.write("%.10f" % (float(abs(interval[1]-interval[0]+1))/list_len[spe][contig]))
          output_merged.write("\n")
             
        len_covered_spe = len_covered_spe + len_covered_contig
        len_ref_spe = len_ref_spe + list_len[spe][contig]
        
        output_contig.write(spe)
        output_contig.write("\t")
        output_contig.write(contig)
        output_contig.write("\t")
        output_contig.write(str(list_len[spe][contig]))
        output_contig.write("\t")
        output_contig.write(str(len_covered_contig))
        output_contig.write("\t")
        output_contig.write("%.10f" % (float(len_covered_contig)/list_len[spe][contig]))
        output_contig.write("\n") 

      #covered = float(len_covered_spe)/len_ref_spe
      #print "%s\t%.3f" % (spe, covered)
      output_spe.write(spe)
      output_spe.write("\t")
      output_spe.write(str(len_ref_spe))
      output_spe.write("\t")
      output_spe.write(str(len_covered_spe))
      output_spe.write("\t")
      output_spe.write("%.10f" % (float(len_covered_spe)/len_ref_spe))
      output_spe.write("\n")
    
    output_all.close()
    output_merged.close()
    output_contig.close()
    output_spe.close()  

# function: calculate for each species / contigs the covered length/percentage of intervals, for single reference mode
def calculate_single(table, list_len, output):        
      output_all = open(output+"_all.txt", 'w')
      output_merged = open(output+"_merged.txt", 'w')
      output_contig = open(output+"_contig.txt", 'w')

      len_covered_spe = 0
      len_ref_spe = 0
      
      for contig in table:
        len_covered_contig = 0
        sorted_by_lower_bound = sorted(table[contig], key=lambda tup: tup[0])
        
        for each_inter in sorted_by_lower_bound:
          output_all.write(contig)
          output_all.write("\t")
          output_all.write(str(each_inter[0]))
          output_all.write("\t")
          output_all.write(str(each_inter[1]))
          output_all.write("\n")

        intervals = merge_intervals(sorted_by_lower_bound)
        
        for interval in intervals:
          len_covered_contig = len_covered_contig + abs(interval[1] - interval[0] + 1)
          output_merged.write(contig)
          output_merged.write("\t")
          output_merged.write(str(list_len[contig]))
          output_merged.write("\t")
          output_merged.write(str(interval[0]))
          output_merged.write("\t")
          output_merged.write(str(interval[1]))
          output_merged.write("\t")
          output_merged.write(str(abs(interval[1]-interval[0]+1)))
          output_merged.write("\t")
          output_merged.write("%.10f" % (float(abs(interval[1]-interval[0]+1))/list_len[contig]))
          output_merged.write("\n") 
        
        #covered_contig = float(len_covered_contig)/list_len[contig]
        
        #len_covered_spe = len_covered_spe + len_covered_contig
        #len_ref_spe = len_ref_spe + list_len[contig]
        
        #print "%s\t%.3f" % (contig, covered_contig)
        
        output_contig.write(contig)
        output_contig.write("\t")
        output_contig.write(str(list_len[contig]))
        output_contig.write("\t")
        output_contig.write(str(len_covered_contig))
        output_contig.write("\t")
        output_contig.write("%.10f" % (float(len_covered_contig)/list_len[contig]))
        output_contig.write("\n")
       
      #covered = float(len_covered_spe)/len_ref_spe
      
      #print "%.3f" % covered

      output_all.close()
      output_merged.close()
      output_contig.close()
        
if __name__ == '__main__':
  
  #Parse Command Line
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', help='input file')
  parser.add_argument("-v", "--verbosity", type=str, choices=['all', 'single'],
                    help="increase output verbosity, all: use patric all reference mode, single: use single reference mode")
  parser.add_argument('-o', '--output', help='output file')
  args = parser.parse_args()

  input = args.input
  output = args.output
    
  #sam_file = str(sys.argv[1]) 
  sam_file = input
  
  if args.verbosity == 'all':
    table, list_len = store_tree_all(sam_file)
    calculate_all(table, list_len, output)
  elif args.verbosity == 'single':
    table, list_len = store_tree_single(sam_file) 
    calculate_single(table, list_len, output)
  
