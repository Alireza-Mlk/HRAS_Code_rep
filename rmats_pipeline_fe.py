##################################################

#
# rMats pipeline to detect existing .bam files fol-
# lowing the convention %SAMPLE_NAME%-T%SAMPLE_TIME%
# -R%SAMPLE_REPLICATE%.star.bam
# The pipeline aims to find and grop same sample ti-
# mes together and compare all of them one to one without duplicates through rMats with predefined (should make a .yml file) settings
#
import sys
import os
import re
from ruffus import *
from ruffus.combinatorics import *
import subprocess
from cgatcore import pipeline as P
from cgatcore import experiment as E
import glob
from cgatcore import iotools


PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


#PARAMS = P.get_parameters()
PARAMS["cwd"] = os.getcwd()
#PARAMS["gtf_dir"] = "/shared/sudlab1/General/annotations/hg38_noalt_ensembl93/ensembl.dir/geneset_all.gtf"
#PARAMS["rmats_threads"] = int(2)
#PARAMS["rmats_mem"] = int(8)
cwd = os.getcwd()


## Group Bam files based on T and collect their paths into a text file ##

#rx = r"/(?P<s_name>.+)-T(?P<s_time>.+)-(?P<s_rep>.+)\.star\.bam"
rx2 = r"/T(?P<TT>.+)\.txt"
@collate("*.bam", formatter(PARAMS["data_rx"]), "./T{s_time[0]}.txt")
def create_txt(input_file_names, output_file):
	with open(output_file, "w") as oo:
                input_file_names = [os.path.join(cwd, fn) for fn in
                                    input_file_names]
                fns = ",".join(input_file_names)
                oo.write(fns)


## Use the txt file to mix and match and run rMats ##
@follows(create_txt, mkdir("./temp.dir"))
@permutations(create_txt, formatter(), 2,
              "./{basename[0][0]}_{basename[1][0]}_od/summary.txt", "{basename[0][0]}", "{basename[1][0]}", "{basename[0][0]}_{basename[1][0]}")
def run_rmats(input_files, output_files, name1, name2, comb_name):
                statement = ''' ~/micromamba/envs/rMats_env/bin/rmats.py 
                --gtf %(gtf_dir)s 
                --tmp ./temp.dir/%(comb_name)s.temp.dir 
                --od %(comb_name)s_od                                                                                                                                        
                --variable-read-length                                                                                                                                                            
                --readLength 150                                                                                                                                                                  
                --b1 %(cwd)s/%(name1)s.txt                                                                                                                                                        
                --b2 %(cwd)s/%(name2)s.txt                                                                                                                                                       
                --anchorLength 1                                                                                                                                                                  
                --nthread %(rmats_threads)s                                                                                                                                                               
                --task both                                                                                                                                                                       
                --novelSS '''
                P.run(statement, job_threads=PARAMS["rmats_threads"], job_memory=PARAMS["rmats_mem"])
## Sorting the summaries innto the same folder##
@follows(run_rmats, mkdir("all_summaries"))
@transform("./*_od/*.txt", formatter(), "./all_summaries/{subdir[0][0]}_{basename[0]}.txt", "{subdir[0][0]}", "{subdir[0][1]}", "{path[0]}", "{basename[0]}")
def all_sum(infile, outfile, spath0, spath1, path0, inname):
       statement = ''' cp %(path0)s/%(inname)s.txt %(cwd)s/all_summaries/%(spath0)s_%(inname)s.txt '''
       P.run(statement)
@follows(all_sum)
@collate("./all_summaries/*.JCEC.txt",
	# formatter("(.+)_([^_]+).MATS.JCEC.txt"),
        #  r"./{2[0]}_Events.tsv", "{basename[0]}" )
        regex("(.+)_([^_]+).MATS.JCEC.txt"),
	r"\2_Events.tsv")
def combine_files(infiles_a, outfiles_a):

    #header = open(infile[0]).readline()
    
   # template_statement = '''awk -v FS="\t" '{n=split(FILENAME,a,"_"); 
              #                               n=split(a[4],f,"."); 
             #                                print $0, a[1], a[2], f[1]}' %(inf)s
            #              | grep -v "IJC_SAMPLE_1" 
           #               | sed -e 's/\,/-/g' >> %(outfile)s '''
                          
    
    #statement = '''  awk '{if (NR == 0) 
     #                        {n=split(FILENAME,a,"_"); 
      #                        print $0, "TA", "TB", "EventType"}
       #                    else {n=split(FILENAME,a,"_"); 
        #                         n=split(a[4],f,"."); 
         #                        print $0, a[1], a[2], f[1]}}' %s 
          #            | sed -e 's/\,/-/g' >> %s ''' % (infiles[0], outfile)
                      
    #for inf in infiles:
    #	statement += " && " + template_statement % locals()
	statement = r''' LSOA=%(LSOO)s &&
        if ((LSOA==0)); then
                awk -v LST=0 -v OFS="\t" '{if (LST == 0) {n=split("%(inname)s",a,"_"); print $0,	"TA",	"TB",	"EventType"; LST++} else {n=split("%(inname)s",a,"_"); n=split(a[4],f,"."); print $0,	a[1],	a[2],	f[1]}}' %(inf)s | sed -e 's/\,/-/g' >> %(outfiles_a)s ;
        else
        awk -v FS="\t" -v OFS="\t" '{n=split("%(inname)s",a,"_"); n=split(a[4],f,"."); print $0,	a[1],	a[2],	f[1]}' %(inf)s | grep -v "IJC_SAMPLE_1" | sed -e 's/\,/-/g' >> %(outfiles_a)s ;
        fi '''

	LSOO = -1
	for inf in infiles_a:
		inname = os.path.basename(inf)
		LSOO += 1
		P.run(statement)
## generic tasks ##

@follows(combine_files)
def full():
        pass







########################################
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
  P.main(sys.argv)

############################



