#!/usr/bin/env python
import os,getopt,subprocess,sys,re,warnings,random

def usage():
	sys.exit('\nUsage:\npython3 RealSVsim.py -b SNP_INDEL_SV.bed(optional) -g genome.fa\n\n	-----------optional-----------\n		-s snp number\n		-d deletion number\n		-i insertion number\n		-u duplication number\n		-v inversion number\n		-p translocation pair\n		-n translocation number\n')

if sys.version_info[0] < 3:
	print('\nPlease Use Python3 !!\n')
	usage()

try:
	opts, args = getopt.getopt(sys.argv[1:], "b:g:d:s:i:u:v:p:n:h",["help"])
except getopt.GetoptError as err:
	print(err)
	usage()

[bed,genome,snp_num,del_num,ins_num,dup_num,inv_num,trans_pair,trans_num]=['']*2+[0]*7
if not opts:
	usage()
for o, a in opts:
	if o =="-b":
		bed=os.path.abspath(a)
	elif o =="-g":
		genome=os.path.abspath(a)
	elif o =="-s":
		snp_num=int(a)
	elif o =="-d":
		del_num=int(a)
	elif o =="-i":
		ins_num=int(a)
	elif o =="-u":
		dup_num=int(a)
	elif o =="-v":
		inv_num=int(a)
	elif o =="-p":
		trans_pair=int(a)
	elif o =="-n":
		trans_num=int(a)
	elif o in ("-h","--help"):
		usage()

if not genome:
	usage()	
if not bed:
	if sum([snp_num,del_num,ins_num,dup_num,inv_num,trans_pair,trans_num])==0:
		usage()

blacklist_file=genome.split('.')[0]+'_masked.bed'

def Masked_genome(genome):
	out=open(blacklist_file,'w')
	g= open(genome,'r')
	info={}
	chromosome=[]
	cont=g.read().split('>')
	del cont[0]
	for i in range(len(cont)):
		m=re.search(r'^([^\n]*)\n([\s\S]*)',cont[i],re.M)
		if m:
			chromosome.append(m.group(1).split(" ")[0])
			info[m.group(1).split(" ")[0]]=m.group(2).replace('\n','').lower()
		else:
			exit("could not find ref sequence for \n"+cont[i]+"\n!\n")
	nucleotide=['A','T','C','G','a','t','c','g']
	for n in range(len(chromosome)):
		start,end=-1,-1
		for m in range(len(info[chromosome[n]])):
			if info[chromosome[n]][m] not in nucleotide:
				if start==-1:
					start=m
					end=m+1
					continue
				if m-end>1000:
					out.write('\t'.join([chromosome[n],str(start),str(end),str(end-start)])+'\n')
					start=m
					end=m+1
				else:
					end=m+1
		if start!=-1:
			out.write('\t'.join([chromosome[n],str(start),str(end),str(end-start)])+'\n')

def Simulator(bed, genome, blacklist_file):
	out_hap1=open(genome+'.hap1.fa','w')
	out_hap2=open(genome+'.hap2.fa','w')
	vcf_out=open(genome+'.hap.vcf','w')
	
	hap1,hap2,info,blacklist={},{},{},{}
	chromosome_total,chromosome_sizes=[],{}
	SNP,DEL,INS,DUP,INV=[],[],[],[],[]
	
	
	g= open(genome,'r')
	cont=g.read().split('>')
	del cont[0]
	for i in range(len(cont)):
		m=re.search(r'^([^\n]*)\n([\s\S]*)',cont[i],re.M)
		if m:
			chromosome_name=m.group(1).split(" ")[0]
			chromosome_sequence=m.group(2).replace('\n','').lower()
			info[chromosome_name]=chromosome_sequence
			hap1[chromosome_name]=[]
			hap2[chromosome_name]=[]
			if len(chromosome_sequence)>=40000000:
				chromosome_total.append(chromosome_name)
				chromosome_sizes[chromosome_name]=len(chromosome_sequence)
		else:
			exit("could not find ref sequence for \n"+cont[i]+"\n!\n")
	
	with open(blacklist_file) as f:
		while True:
			l=f.readline().rstrip()
			if not l: break
			line=l.split("\t")
			if not line[0] in blacklist:
				blacklist[line[0]]=[[line[1],line[2]]]
			else:
				blacklist[line[0]].append([line[1],line[2]])
	
	def pick_site(chromosome,chromosome_total,chromosome_sizes,blacklist):
		while True:
			inside,outside=0,0
			if not chromosome in chromosome_total:
				item=random.randint(0,len(chromosome_total)-1)
				chromosome=chromosome_total[item]
			length=chromosome_sizes[chromosome]
			position=random.randint(1,int(length))
			if chromosome in list(blacklist.keys()):
				for i in range(len(blacklist[chromosome])):
					if position >=int(blacklist[chromosome][i][0])-1000 and position <=int(blacklist[chromosome][i][1])+1000:
						inside=1
						break
					elif position <int(blacklist[chromosome][i][0])-1000:
						outside=1
						blacklist[chromosome].insert(i,[str(position),str(position+1)])
						break
					elif i==len(blacklist[chromosome])-1 and position >int(blacklist[chromosome][i][1])+1000:
						outside=1
						blacklist[chromosome].append([str(position),str(position+1)])
			else:
				outside=1
				blacklist[chromosome]=[[str(position),str(position+1)]]
			if outside==1:
				return([chromosome,str(position),blacklist])
				break

	if bed:
		with open(bed) as f:
			while True:
				l=f.readline().rstrip()
				if not l: break
				if l[0]=='#': continue
				line=l.split("\t")
				if not line[0] in blacklist:
					blacklist[line[0]]=[[line[1],line[2]]]
				else:
					blacklist[line[0]].append([line[1],line[2]])
				if line[3]=="SNP":
					SNP.append(l)
				elif line[3]=="DEL":
					DEL.append(l)
				elif line[3]=="DUP":
					DUP.append(l)
				elif line[3]=="INS":
					INS.append(l)
				elif line[3]=="INV":
					INV.append(l)
	for i in list(blacklist.keys()):
		 blacklist[i]=sorted(blacklist[i],key=lambda x: (int(x[0]), int(x[1])))
	
	total=[[],[],[],[],[]]
	for x in range(len(total)):
		if x==0:
			if snp_num==0:
				number=len(SNP)
			else:
				number=snp_num
			svtype="SNP"
			while len(total[x])<number:
				if len(SNP)>=1:
					line=SNP.pop(0).split("\t")
				else:
					chromosome,position,blacklist=pick_site("chr0",chromosome_total,chromosome_sizes,blacklist)
					line=[chromosome,position,"-","-",0,"-"]
				ref=info[line[0]][int(line[1]):int(line[1])+1].upper()
				if not line[-1][0].upper in ['A','T','C','G']:
					a=['A','T','C','G']
					a.remove(ref)
					allele=random.choices(a,k=1)
				else:
					allele=line[-1].upper()
				site=str(int(line[1])+1)
				total[x].append(line[0]+"\t"+site+'\t'+ref+'\t'+allele)
		if x==1:
			if del_num==0:
				number=len(DEL)
			else:
				number=del_num
			svtype="DEL"
			while len(total[x])<number:
				if len(DEL)>=1:
					line=DEL.pop(0).split("\t")
				else:
					chromosome,position,blacklist=pick_site("chr0",chromosome_total,chromosome_sizes,blacklist)
					length=random.randint(50,10000)
					line=[chromosome,position,"-","-",length,"-"]
				ref=info[line[0]][int(line[1]):int(line[1])+int(line[4])+1].upper()
				allele=info[line[0]][int(line[1])].upper()
				site=str(int(line[1]))
				total[x].append(line[0]+"\t"+site+'\t'+ref+'\t'+allele)
		elif x==2:
			if ins_num==0:
				number=len(INS)
			else:
				number=ins_num
			svtype="INS"
			while len(total[x])<number:
				if len(INS)>=1:
					line=INS.pop(0).split("\t")
				else:
					chromosome,position,blacklist=pick_site("chr0",chromosome_total,chromosome_sizes,blacklist)
					length=random.randint(50,10000)
					line=[chromosome,position,"-","-",length,"-"]
				ref=info[line[0]][int(line[1])+1].upper()
				if not line[-1][0].upper in ['A','T','C','G']:
					allele=info[line[0]][int(line[1])+1].upper()+''.join(random.choices('ATCG',k=int(line[4])))
				else:
					allele=info[line[0]][int(line[1])+1].upper()+line[-1].upper()
				site=str(int(line[1])+1)
				total[x].append(line[0]+"\t"+site+'\t'+ref+'\t'+allele)
		elif x==3:
			if dup_num==0:
				number=len(DUP)
			else:
				number=dup_num
			svtype="DUP"
			while len(total[x])<number:
				if len(DUP)>=1:
					line=DUP.pop(0).split("\t")
				else:
					chromosome,position,blacklist=pick_site("chr0",chromosome_total,chromosome_sizes,blacklist)
					length=random.randint(50,10000)
					line=[chromosome,position,"-","-",length,"-"]
				ref=info[line[0]][int(line[1])].upper()
				allele=info[line[0]][int(line[1]):int(line[1])+int(line[4])+1].upper()
				site=str(int(line[1]))
				total[x].append(line[0]+"\t"+site+'\t'+ref+'\t'+allele)
		elif x==4:
			if inv_num==0:
				number=len(INV)
			else:
				number=inv_num
			svtype="INV"
			while len(total[x])<number:
				if len(INV)>=1:
					line=INV.pop(0).split("\t")
				else:
					chromosome,position,blacklist=pick_site("chr0",chromosome_total,chromosome_sizes,blacklist)
					length=random.randint(50,10000)
					line=[chromosome,position,"-","-",length,"-"]
				sequence=info[line[0]][int(line[1])+1:int(line[1])+int(line[4])+1].upper()
				inversion=list(sequence)
				for i in range(len(inversion)):
					if inversion[i]=='A':
						inversion[i]='T'
					elif inversion[i]=='T':
						inversion[i]='A'
					elif inversion[i]=='G':
						inversion[i]='C'
					elif inversion[i]=='C':
						inversion[i]='G'
				ref=info[line[0]][int(line[1])].upper()+sequence
				allele=info[line[0]][int(line[1])].upper()+''.join(inversion)[::-1]
				site=str(int(line[1]))
				total[x].append(line[0]+"\t"+site+'\t'+ref+'\t'+allele)
		sub_group=random.sample(total[x],number)
		for i in range(len(sub_group)):
			name,site,ref,allele=sub_group[i].split('\t')
			genotype=random.randint(1,3)
			if genotype == 1:
				vcf_out.write("\t".join([name,site,".",ref,allele,".",".",svtype,"GT","1|1"])+"\n")
				hap1[name].append(site+'\t'+ref+'\t'+allele)
				hap2[name].append(site+'\t'+ref+'\t'+allele)
			elif genotype == 2:
				vcf_out.write("\t".join([name,site,".",ref,allele,".",".",svtype,"GT","1|0"])+"\n")
				hap1[name].append(site+'\t'+ref+'\t'+allele)
			elif genotype == 3:
				vcf_out.write("\t".join([name,site,".",ref,allele,".",".",svtype,"GT","0|1"])+"\n")
				hap2[name].append(site+'\t'+ref+'\t'+allele)
	
	for i in range(2):
		translocation={}
		new_chromosomes=chromosome_total[:]
		new_chromosomes.remove("chrY")
		for n in range(trans_pair):
			item1=random.randint(0,len(new_chromosomes)-1)
			chromosome1=new_chromosomes[item1]
			del new_chromosomes[item1]
			item2=random.randint(0,len(new_chromosomes)-1)
			chromosome2=new_chromosomes[item2]
			del new_chromosomes[item2]
			event=[]
			for x in range(trans_num):
				chromosome,position,blacklist=pick_site(chromosome1,chromosome_total,chromosome_sizes,blacklist)
				event.append([chromosome,position])
				chromosome,position,blacklist=pick_site(chromosome2,chromosome_total,chromosome_sizes,blacklist)
				event.append([chromosome,position])
			event_sorted=sorted(event,key=lambda x: (x[0], int(x[1])))
			translocation[event_sorted[0][0]]=event_sorted
			if i ==0:
				haplotype="1|0"
			elif i==1:
				haplotype="0|1"
			vcf_out.write("\t".join([event_sorted[0][0],event_sorted[0][1],".",".",event_sorted[3][0]+":"+event_sorted[3][1],".",".","TRANS","GT",haplotype])+"\n")
			vcf_out.write("\t".join([event_sorted[1][0],event_sorted[1][1],".",".",event_sorted[4][0]+":"+event_sorted[4][1],".",".","TRANS","GT",haplotype])+"\n")
			vcf_out.write("\t".join([event_sorted[2][0],event_sorted[2][1],".",".",event_sorted[5][0]+":"+event_sorted[5][1],".",".","TRANS","GT",haplotype])+"\n")
			vcf_out.write("\t".join([event_sorted[3][0],event_sorted[3][1],".",".",event_sorted[0][0]+":"+event_sorted[0][1],".",".","TRANS","GT",haplotype])+"\n")
			vcf_out.write("\t".join([event_sorted[4][0],event_sorted[4][1],".",".",event_sorted[1][0]+":"+event_sorted[1][1],".",".","TRANS","GT",haplotype])+"\n")
			vcf_out.write("\t".join([event_sorted[5][0],event_sorted[5][1],".",".",event_sorted[2][0]+":"+event_sorted[2][1],".",".","TRANS","GT",haplotype])+"\n")
		if i ==0:
			translocation_hap1=translocation
		elif i==1:
			translocation_hap2=translocation
	
	for y in range(2):
		if y ==0:
			hap=hap1
			translocation_hap=translocation_hap1
		elif y ==1:
			hap=hap2
			translocation_hap=translocation_hap2
		new_info={}
		for name in list(hap.keys()):
			seq=info[name]
			new=list(seq)
			for n in range(len(hap[name])):
				site,left,right=hap[name][n].split('\t')
				i=int(site)-1
				if new[i]=='':
					continue
				for x in range(len(left)):
					new[i+x]=''
				new[i]=right
			new_info[name]=new
		for chromosome1 in list(translocation_hap.keys()):
			sequence1=new_info[chromosome1][:]
			chromosome2=translocation_hap[chromosome1][3][0]
			sequence2=new_info[chromosome2][:]
			new_info[chromosome1][int(translocation_hap[chromosome1][2][1]):]=sequence2[int(translocation_hap[chromosome1][5][1]):]
			new_info[chromosome1][int(translocation_hap[chromosome1][0][1]):int(translocation_hap[chromosome1][1][1])]=sequence2[int(translocation_hap[chromosome1][3][1]):int(translocation_hap[chromosome1][4][1])]
			new_info[chromosome2][int(translocation_hap[chromosome1][5][1]):]=sequence1[int(translocation_hap[chromosome1][2][1]):]
			new_info[chromosome2][int(translocation_hap[chromosome1][3][1]):int(translocation_hap[chromosome1][4][1])]=sequence1[int(translocation_hap[chromosome1][0][1]):int(translocation_hap[chromosome1][1][1])]
		for chromosome in list(new_info.keys()):
			if y ==0:
				out_hap1.write('>'+chromosome+'\n'+''.join(new_info[chromosome])+'\n')
			elif y ==1:
				out_hap2.write('>'+chromosome+'\n'+''.join(new_info[chromosome])+'\n')

if __name__=="__main__":
	print("Get NNNN regions in the genome...\n")
	Masked_genome(genome)
	print("Simulating genome...\n")
	Simulator(bed, genome, blacklist_file)
	
