import argparse
import os, sys
import shutil
import multiprocessing.dummy as mp
parser = argparse.ArgumentParser()
parser.description='This script is used to run GWAS analysis with GLM parallely.'
parser.add_argument('-in_gen',type = str,help=('input genotype file name with bfile (bed, bim, and fam) format, required'),required=True)
parser.add_argument('-pheno',type = str,help=('input file, phenotype residual, adjusted by covariate, required'),required=True)
parser.add_argument('-parallel',type = str,default = 30, help=('parallel computing, default is 30, optional'),required=False)
parser.add_argument('-pvalue',type = str,default = 5e-7, help=('p value filter, optional'),required=False)
parser.add_argument('-odir',type = str,default = './', help=('output dir'),required=False)
args = parser.parse_args()

# creat dir
if not os.path.exists(args.odir):
    os.mkdir(args.odir)

phenotype_name = args.pheno.split('/')[-1]
phenotype_dir = os.path.dirname(os.path.abspath(args.pheno))
phenotype_full_path = os.path.abspath(args.pheno)
gen_name = args.in_gen.split('/')[-1]
gen_dir = os.path.dirname(os.path.abspath(args.in_gen))
gen_full_path = os.path.abspath(args.in_gen)



def check(genotype,phenotype):
    subject_from_genotype =[]
    genotype_file_name = genotype+'.fam'
    infile = open (genotype_file_name,'r')
    for myline in infile:
        myitem = myline.strip().split()
        if myitem[1] not in subject_from_genotype:
            subject_from_genotype.append(myitem[1])
        else:
            print (myitem[1] +' exists repeatly in genotype.'+'\n')
            exit()
    infile.close()
    subject_from_phenotype =[]

    phenotype_dic = {}
    infile = open (phenotype,'r')

    new_phenotype_name = os.getcwd() +"/"+ args.odir +'/'+".".join(args.pheno.split('/')[-1].split('.')[:-1])+'_order.txt'

    outfile = open (new_phenotype_name,'w')
    head = infile.readline()
    outfile.write(head)
    for myline in infile:
        myline = myline.strip()
        myitem = myline.split('\t')
        length = len(myitem)-1
        if myitem[0] not in subject_from_phenotype:
            subject_from_phenotype.append(myitem[0])
            phenotype_dic[myitem[0]] = myline
        else:
            print (myitem[0] +' exists repeatly in phenotype.'+'\n')
    if subject_from_genotype == subject_from_phenotype:
        print ('The subject order are same.')
        full_file = os.path.abspath(args.pheno)
        shutil.copy(full_file, new_phenotype_name)
    elif subject_from_genotype != subject_from_phenotype or len(subject_from_genotype) != len(subject_from_phenotype):
        print ('The subjects\' order or the number of subjects are not same!!!'+'\n'+'Ordering the subject in phenotype file.'+'\n')
        for item in subject_from_genotype:
            if item in phenotype_dic:
                outfile.write(phenotype_dic[item]+'\n')
            else:
                outfile.write(item+('\t'+'Na')*length+'\n')
    infile.close()
    outfile.close()

def split_phenotype(pheno):
    content = []
    phen_file = open (pheno,'r')
    dir_split = os.getcwd()+'/'+args.odir+'/phenotype_split/'
    if not os.path.exists (dir_split):     
        os.makedirs(dir_split)
    head = phen_file.readline()
    name = head.strip().split("\t")
    for myline in phen_file:
        myline = myline.strip()
        content.append(myline)
    for i in range(1,len(name)):
        file_name = dir_split+name[i]+".txt"
        outfile = open (file_name,"w")
        outfile.write("FID"+"\t"+"IID"+"\t"+name[i]+"\n")
        for myitem in content:
            item = myitem.strip().split("\t")
            outfile.write('0'+"\t"+item[0]+"\t"+item[i]+"\n")
        outfile.close()
    phen_file.close()

def glm(myfile):
    do_glm = sys.path[0]+'/bin/plink2 --bfile ' +args.in_gen +' --pheno ' +myfile+ ' --pfilter '+str(args.pvalue)+' --glm hide-covar --out ' +args.odir +'/'+   '.'.join(myfile.split('/')[-1].split('.')[:-1]) +' --adjust '+"\n"  #May 26,2022 remove --variance-standardize  add --pfilter
    os.system (do_glm)

#checking subject between genotype and phenotype and ordering
check(args.in_gen,args.pheno)

#split phenotype into single file
files = os.listdir(os.getcwd()+'/'+args.odir)
for myitem in files:

    if 'order.txt' in myitem:
        new_phenotype_name = os.getcwd()+'/'+args.odir+'/'+myitem


split_phenotype(new_phenotype_name)

#perform glm analysis
files = os.listdir(os.getcwd()+'/'+args.odir+'/phenotype_split/')

pheno_files = []
for myitem in files:
    pheno_files.append(os.getcwd()+'/'+args.odir+'/phenotype_split/'+myitem)
pheno_files.sort()

if __name__=="__main__":
    p=mp.Pool(int(args.parallel))
    p.map(glm,pheno_files)
    p.close()
    p.join()

#May 26, 2022  get full gwas result for the first trait. The markers' p-value might be NA, which will be filtered out in the next steps.
first_trait_glm = sys.path[0]+'/bin/plink2 --bfile ' +args.in_gen +' --pheno ' +pheno_files[0]+ ' --glm hide-covar --out ' +args.odir +'/'+   '.'.join(pheno_files[0].split('/')[-1].split('.')[:-1]) +' --adjust '+"\n"
os.system(first_trait_glm)
