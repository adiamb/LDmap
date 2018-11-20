"""
scripts to parse meta results from GWAS and match the LD structure to 1000 genomes and find genotyped snps across cohorts.
@author: Aditya Ambati ambati@stanford.edu, Mignot Lab, Stanford University
"""
import sys
import argparse
import datetime
import os
import numpy as np
dttime=datetime.datetime.now().strftime ("%Y%m%d")

parser = argparse.ArgumentParser(description='A class method to map LD structure of top snps with 1000genomes and further indicate if top snps from GWAS meta were genotyped or not')
parser.add_argument('-meta', required=True, help='Meta results file with column order of <chr,rsid,pos,allele_A,allele_B,P_value> delimited by single space sorted by position')
parser.add_argument('-ld', required=True, help='ld file from impute files use the gtool to convert to plink and get LD window')
parser.add_argument('-refld', required=True, help='ld file from 1000 genomes ,use same command as the above to get LD window')
parser.add_argument('-snplist', required=True, help='SNPLIST of the orginal unimputed genotypes to check whether top variants are imputed or genotyped')
parser.add_argument('-chrnum', required=True, help='Chromosome number just a integer')
parser.add_argument('-outfile', required=True, help='full path to an outfile NO suffix please such as (.txt or .csv)')
### parse the arguments
args=parser.parse_args()
ldfile = args.ld
refldfile = args.refld
snpfile = args.snplist
meta = args.meta
chrnum =args.chrnum
outfile = args.outfile

####### define the LDmap as an object ########
class LD_map(object):
	"""A class of an LD object:

	Attributes:
		name: A string representing the LD object.
		ldfile: Filelist of LD files from impute files.
		reffile:Filelist of LD files from 1000 genomes
		snpfile: SNPLIST of the orginal unimputed genotypes
	"""
	instance_count =0

	def __init__(self, name, meta, ldfile, reffile, snpfile, outfile):
		"""intitate the sseq object class."""
		self.name = name
		self.ldfile = ldfile
		self.reffile =reffile
		self.snpfile = snpfile
		self.meta = meta
		self.chrnum = int(chrnum)
		self.outfile = outfile
		self.instance_count += 1

	def get_attr(self):
		print(' NAME :- {} \n SOURCE LD FILE :- {} \n REF LD FILE :- {} \n META RESULTS :- {} \n CHR NUMBER :- {} \n OUTFILE :- {} \n INSTANCE COUNT :- {}'.format(self.name, self.ldfile, self.reffile, self.meta, self.chrnum, self.outfile, self.instance_count))

	@staticmethod
	def process_ldlist(file, threshold):
		LD_snps = {}	
		with open(file) as ld_in:
			for n, line in enumerate(ld_in):
				if n > 0:
					line_parse = line.strip().split(' ')
					print line_parse
					make_key = line_parse[2]
					make_val = ','.join(line_parse[5:])
					if float(line_parse[6]) >= threshold:
						if make_key in LD_snps:
							get_snp = LD_snps.get(make_key)
							LD_snps[make_key] = get_snp + '$'+make_val
						else:
							LD_snps[make_key] = make_val
		return LD_snps

	#@staticmethod
	def process_snplist(self):
		genotyped_dic={}
		with open(self.snpfile) as snplist:
			for file in snplist:
				file = file.strip()
				if file:
					mak_name = file.replace('.snplist', '')
					with open(file) as snplist_in:
						for line in snplist_in:
							line_parse= line.strip()
							if line_parse in genotyped_dic:
								get_gen = genotyped_dic.get(line_parse)
								genotyped_dic[line_parse] = get_gen + '&' + mak_name
							else:
								genotyped_dic[line_parse] = mak_name
		return genotyped_dic

	@staticmethod
	def lead_snp(parse_loci, genotyped_dic):
		#global rsid_track
		if len(parse_loci) >= 0:
			pval_track = 1
			rsid_track = ''
			for n, item in enumerate(parse_loci):
				parse_item = item.split(' ')
				if parse_item[1] in genotyped_dic:
					if float(parse_item[5]) < pval_track:
						pval_track = float(parse_item[5])
						rsid_track = parse_item[1]
			if rsid_track:
				return rsid_track
		# else:
		# 	print 'input meta file has inconsistent column number '

	def make_loci(self, distance_kb):
		pos_track=[]
		pos_d =[]
		loci_n = 0
		make_loci={}
		header=''
		print('the LD distance has been selected as {} KB'.format(distance_kb))
		with open(self.meta) as top_snp_in:
			for n, line in enumerate(top_snp_in):
				if n > 0:
					line_parse = line.split(' ')
					if len(line_parse) == 6:
						assert len(line_parse) == 6
						if int(line_parse[0]) == self.chrnum:
							assert int(line_parse[0]) == self.chrnum
							pos_track.append(int(line_parse[2]))
							if n > 1:
								pos_diff = pos_track[-2]-int(line_parse[2])
								print pos_diff
								if np.absolute(pos_diff) < distance_kb:#250000 :## 100kb
									if loci_n in make_loci:
										get_loci = make_loci.get(loci_n)
										make_loci[loci_n] = get_loci +'%'+ line.strip()#','.join(line_parse[:3])
								else:
									loci_n += 1
									make_loci[loci_n] = line.strip()#','.join(line_parse[:3])
									print line
							else:
								loci_n += 1
								make_loci[loci_n] = line.strip()
						else:
							print('Have you input the right chromosome number at the command line ? INPUT at command is {} while the meta file has {}'.format(self.chrnum, line_parse[0]))
					else:
						print('input meta file {} has inconsistent column number at row number {} '.format(self.meta, n))
				else:
					header=line.strip().split(' ')
					#header=header.split(' ')
		print('IDENTIFIED {} LOCI FROM {} '.format(loci_n, self.chrnum))
		return make_loci

	def process_loci(self):
		loci=self.make_loci(distance_kb=250000)
		genotyped_dic = self.process_snplist()
		LD_snps = self.process_ldlist(file=self.ldfile, threshold=0.2)
		LD_1kg_snps = self.process_ldlist(file=self.reffile, threshold=0.2)
		outfile = open(self.outfile+'_ANNOT_'+dttime+'.txt', 'w')
		outfile.write('chr rsid pos allele_A allele_B P_value LOCINO LEADSNP-LD 1000G_LD GENOTYPED'+'\n')

		for key, loci in loci.iteritems():
			parse_loci = loci.split('%')
			rsid_track=self.lead_snp(parse_loci, genotyped_dic)
			print 'FROM LOCI no ', key, len(parse_loci), ' LEAD SNP IDENTIFIED AS ', rsid_track 
			for item in parse_loci:
				parse_item = item.split(' ')
				outfile.write(' '.join(parse_item)+' '+str(key)+' ')
				if rsid_track:
					if rsid_track in LD_snps:
						get_ld_snps =  {k:v for k,v in (x.split(',') for x in LD_snps.get(rsid_track).split('$')) }
						if parse_item[1] in get_ld_snps:
							print ' '.join(parse_item[:6]), get_ld_snps.get(parse_item[1]), 'LEAD SNP ', rsid_track
							outfile.write(rsid_track+'-'+get_ld_snps.get(parse_item[1]))
						else:
							outfile.write('-')
					else:
						outfile.write('-')
					if rsid_track in LD_1kg_snps:
						get_1k_ld_snps = {k:v for k,v in (x.split(',') for x in LD_1kg_snps.get(rsid_track).split('$') ) }
						if parse_item[1] in get_1k_ld_snps:
							print ' '.join(parse_item[:6]), get_1k_ld_snps.get(parse_item[1])
							outfile.write(' '+rsid_track+'-'+get_1k_ld_snps.get(parse_item[1]))
						else:
							outfile.write(' '+'-')
					else:
						outfile.write(' '+'-')
				else: ## if lead snps cannot be identified then write out
					outfile.write('-'+' '+'-')

				if parse_item[1] in genotyped_dic:
					outfile.write(' '+genotyped_dic.get(parse_item[1]))
				else:
					outfile.write(' '+'IMPUTED')
				outfile.write('\n')
		outfile.close()

if __name__ == '__main__':
	LD_object = LD_map(name=str(chrnum)+'LD', meta=meta, ldfile=ldfile, reffile=refldfile, snpfile=snpfile, outfile=outfile)
	LD_object.get_attr()
	LD_object.process_loci()




# ldfile = args.ld
# refldfile = arg.refld
# snpfile = arg.snplist
# meta = arg.meta
# chrnum =arg.chrnum
# outfile = arg.outfile




	
