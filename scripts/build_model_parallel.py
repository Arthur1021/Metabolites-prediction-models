import argparse
import os, sys
import shutil
import multiprocessing.dummy as mp
parser = argparse.ArgumentParser()
parser.description='This script is used to build model and perform fusion analysis.'
parser.add_argument('-in_gen',type = str,help=('input genotype file name with bfile (bed, bim, and fam) format,used for generating 1000g LD list file, required'),required=True)
parser.add_argument('-pheno',type = str,help=('input file, phenotype file, required'),required=True)
parser.add_argument('-bfdir',type = str,help=('input dir, extract up and down 100k bfile, used for build model, required'),required=True)
parser.add_argument('-parallel',type = str,default = 30, help=('parallel computing, default is 30, optional'),required=False)
parser.add_argument('-odir',type = str,default = './', help=('output dir, default is ./, do not change it'),required=False)
args = parser.parse_args()

#creat dirs
if not os.path.exists(args.odir):
    os.mkdir(args.odir)

pheno_dir = args.odir+'/phenotype_split/'
if not os.path.exists(pheno_dir):
    os.mkdir(pheno_dir)

os.system("mkdir -p phenotype_split")
os.system("mkdir -p TWAS_fusion/tmp/")
os.system("mkdir -p TWAS_fusion/out/")
os.system("mkdir -p output/TWAS_fusion/tmp/")
os.system("mkdir -p output/TWAS_fusion/out")




#split each phenotype into single file
phen_file = open(args.pheno,"r")
content = []
head = phen_file.readline()
name = head.strip().split("\t")
for myline in phen_file:
    myline = myline.strip()
    content.append(myline)

for i in range(1,len(name)):
    file_name = args.odir+"/phenotype_split/"+name[i]+".txt"
    outfile = open (file_name,"w")
    outfile.write("FID"+"\t"+"IID"+"\t"+name[i]+"\n")
    for myitem in content:
        item = myitem.strip().split("\t")
        outfile.write('0'+"\t"+item[0]+"\t"+item[i]+"\n")
    outfile.close()

phen_file.close()

path = args.bfdir
files = os.listdir( path )
full_path = os.path.abspath(path)
fbile = []
for myfile in files:
    if "unambig.Z.bed" in myfile:
        fbile.append(myfile[:-4])
print (len(fbile))
fbile.sort()

def model_TWAS_fusion(myfile):
    print ("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/bin/FUSION.compute_weights.modified.R --bfile " +full_path+"/"+ myfile + " --tmp TWAS_fusion/tmp/" + myfile.split(".")[0]+ ".all.100kb --out TWAS_fusion/out/" + myfile.split(".")[0] +".all.100kb --verbose 2 --pheno phenotype_split/" + myfile.split(".")[0]+ ".txt --noclean --PATH_gcta /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/bin/gcta_nr_robust --PATH_gemma /mnt/lvm_vol_1/dghoneim/EGA/gemma-0.98.1-linux-static")
    os.system("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/bin/FUSION.compute_weights.modified.R --bfile " +full_path+"/"+ myfile + " --tmp TWAS_fusion/tmp/" + myfile.split(".")[0]+ ".all.100kb --out TWAS_fusion/out/" + myfile.split(".")[0] +".all.100kb --verbose 2 --pheno phenotype_split/" + myfile.split(".")[0]+ ".txt --noclean --PATH_gcta /mnt/lvm_vol_1/sliu/work/pipeline/TWAS_fusion/bin/gcta_nr_robust --PATH_gemma /mnt/lvm_vol_1/dghoneim/EGA/gemma-0.98.1-linux-static")
    print ("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/read.wgt.RDat.R TWAS_fusion/out/"+ myfile.split(".")[0]+".all.100kb.wgt.RDat TWAS_fusion/out/" + myfile.split(".")[0] +".all.100kb.wgt.RDat.ld.list")
    os.system("/mnt/lvm_vol_1/dghoneim/EGA/R-3.6.1/bin/Rscript /mnt/lvm_vol_2/sliu/pipeline/TWAS_fusion/bin/read.wgt.RDat.R TWAS_fusion/out/"+ myfile.split(".")[0]+".all.100kb.wgt.RDat TWAS_fusion/out/" + myfile.split(".")[0] +".all.100kb.wgt.RDat.ld.list")
    print ("plink --bfile " + args.in_gen + " --allow-extra-chr --extract TWAS_fusion/out/"+myfile.split(".")[0]+ ".all.100kb.wgt.RDat.ld.list --make-bed --out " +myfile.split(".")[0] +".ld.1000g.Z")
    os.system("plink --bfile " + args.in_gen + " --allow-extra-chr --extract TWAS_fusion/out/"+myfile.split(".")[0]+ ".all.100kb.wgt.RDat.ld.list --make-bed --out " +myfile.split(".")[0] +".ld.1000g.Z")
    print ("echo \"WGT     ID      CHR     P0      P1\" > " +myfile.split(".")[0]+".pos")
    os.system("echo \"WGT     ID      CHR     P0      P1\" > " +myfile.split(".")[0]+".pos")
    print ("echo \""+myfile.split(".")[0]+".all.100kb.wgt.RDat       "+myfile.split(".")[0]+"    Z       1       9999999999\" >> " +myfile.split(".")[0] +".pos")
    os.system("echo \""+myfile.split(".")[0]+".all.100kb.wgt.RDat       "+myfile.split(".")[0]+"    Z       1       9999999999\" >> " +myfile.split(".")[0] +".pos")



if __name__=="__main__":
    p=mp.Pool(int(args.parallel))
    p.map(model_TWAS_fusion,fbile)
    p.close()
    p.join()
