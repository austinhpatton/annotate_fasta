#!/Users/jeremias/.pyenv/versions/anaconda3-5.0.0/bin/python
# Given a fasta file and the corresponding blast output this script annotates the fasta with the hit resulting in the longest alignment and covering the beginning of the query
#     Copyright (C) 2018  Jeremias N. Brand

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
from Bio import SeqIO
from collections import defaultdict


class BlastRes():
    def __init__(self, qseqid, sseqid, pident, qlen, length, mismatch, 
		gapopen, qstart, qend, sstart, send, evalue, bitscore):
        """
		This class stores the results of a blast search in 
		-outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"
		"""

        self.qseqid   = qseqid 
        self.sseqid   = sseqid 
        self.pident   = pident 
        self.qlen     = qlen 
        self.length   = length 
        self.mismatch = mismatch
        self.gapopen  = gapopen 
        self.qstart   = qstart 
        self.qend     = qend  
        self.sstart   = sstart 
        self.send     = send  
        self.evalue   = evalue 
        self.bitscore = bitscore
        
        
def check_hit_position(BlastRes_obj, current_coverage, max_query_start = 1, max_subject_start = 3, min_coverage = 0.95):
    """
    Checks if a hit stored in a BlastRes class object meets the hit criteria we give
    returns a boolean, False if the conditions are not met

    """
    res = False
    if (int(BlastRes_obj.qstart) <= max_query_start and  int(BlastRes_obj.sstart) <= max_subject_start):
    	hit_cov = int(BlastRes_obj.length)/int(BlastRes_obj.qlen)
    	if (hit_cov > min_coverage):
        	res = True
    return res


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("input_fasta", type=str, help="fasta file to annotate")
	parser.add_argument("blast_output", type=str, help="blast output of input_fasta")
	parser.add_argument("output_file", type=str, help="name of the outfile")
	parser.add_argument("max_query_start", type=int, help="maximum start position of blast alignment in query")
	parser.add_argument("max_subject_start", type=int, help="maximum start position of blast alignment in subject")
	parser.add_argument("min_coverage", type=float, help="minimum fraction of the query hit by the alignment needed to annotate")
	args = parser.parse_args()

	seq_handle = args.input_fasta
	blast_output = args.blast_output
	output_file = args.output_file
	max_query_start = args.max_query_start
	max_subject_start = args.max_subject_start
	min_coverage = args.min_coverage
	# read all sequences into the file
	record_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

	# the defaultdict allows to append items
	blast_dict = defaultdict(list)  
	# populate the dict
	with open(blast_output, "r") as blast_f:
	    for line in blast_f:
	        blast_res = line.strip().split("\t")
	        name = blast_res[0]
	        blast_dict[name].append(BlastRes(*blast_res))
	print(blast_dict) 


	with open(output_file, "w") as output_handle:
	    for record in record_dict:
	        # print(record, record_dict[record],blast_dict[record_dict[record].id] )
	        v = blast_dict[record_dict[record].id]
	        coverage = 0
	        annotation = "No hit"
	        if len(v) > 1:
	            for res in v:
	                hit_cov = int(res.length)/int(res.qlen)
	                if check_hit_position(res, coverage,
	                 max_query_start = max_query_start, max_subject_start = max_subject_start, min_coverage = min_coverage):
	                    coverage = hit_cov 
	                    annotation = res.sseqid + " cov: " + str(coverage)
	        elif len(v) == 1:
	            res = v[0]
	            hit_cov = int(res.length)/int(res.qlen)
	            if check_hit_position(res, coverage,
	                 max_query_start = max_query_start, max_subject_start = max_subject_start, min_coverage = min_coverage):
	                    coverage = hit_cov
	                    annotation = res.sseqid + " cov: " + str(coverage)
	        else:
	            print("No entries for " + record_dict[record].id + "!")


	        record_dict[record].description = record_dict[record].id +" Annotation: "+ annotation
	        SeqIO.write(record_dict[record], output_handle, "fasta")
        


if __name__ == "__main__":
	main()
