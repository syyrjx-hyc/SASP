import os,subprocess,csv,pickle,sys
from Bio.Seq import Seq
import numpy as np
from xgboost import XGBClassifier

def GetParams(string,list):
	i = list.index(string)
	return list[i+1]

names = ["AMC","AMP","AZI","FOX","TIO","CRO","CHL","CIP","GEN","KAN","NAL","STR","FIS","TET","SXT"]
kinds = ["ARP","MIC"]
base2num = {"A":1,"C":2,"G":3,"T":4}
blast_dir = "/home/dell/Downloads/ncbi-blast-2.9.0+/bin"

def find_mid(line,start,end):
	if start in line and end in line:
		mid_string = line.split(start)[1].split(end)[0]
	return mid_string

def CoreProfile(name,kind,query,out_dir):
	core_seq_file = os.path.join(os.getcwd(),"profile",name+"_"+kind+"_core.fna")
	core_profile = os.path.join(os.getcwd(),"profile",name+"_"+kind+"_core.profile")
	snp_pos = {}
	with open(core_profile,"r") as f:
		for line in f.readlines():
			items = line.strip().split("	")
			snp_pos[items[0]] = {}
			for item in items[1:]:
				snp_pos[items[0]][item.split("|")[0]] = item.split("|")[1]
	
	cmd = [os.path.join(blast_dir,"makeblastdb"),\
	"-in",core_seq_file,"-out",os.path.join(out_dir,"core_db"),"-dbtype","nucl"]
	subprocess.check_output(cmd)

	cmd = [os.path.join(blast_dir,"blastn"),\
	"-query",query,"-db",os.path.join(out_dir,"core_db")\
	,"-out",os.path.join(out_dir,"core_align"),"-outfmt","5"]
	subprocess.check_output(cmd)

	snp_count = {}
	with open(os.path.join(out_dir,"core_align"),"r",encoding="utf-8") as f:
		content = list(f.readlines())
		gene_name = ""
		hseq = ""
		qseq = ""
		hlen = 0
		hstart = 0
		hend = 0
		for i in range(len(content)):
			if "<Hit_def>" in content[i]:
				gene_name = find_mid(content[i],"<Hit_def>","</Hit_def>")
				if gene_name not in snp_count:
					snp_count[gene_name] = {}
			if "<Hit_len>" in content[i]:
				hlen = int(find_mid(content[i],"<Hit_len>","</Hit_len>"))
			if "<Hsp_qseq>" in content[i]:
				qseq = find_mid(content[i],"<Hsp_qseq>","</Hsp_qseq>")
				hseq = find_mid(content[i+1],"<Hsp_hseq>","</Hsp_hseq>")
				total_gap = 0
			if "<Hsp_hit-from>" in content[i]:
				hstart = int(find_mid(content[i],"<Hsp_hit-from>","</Hsp_hit-from>"))-1
				hend = int(find_mid(content[i+1],"<Hsp_hit-to>","</Hsp_hit-to>"))-1
			if "<Hsp_midline>" in content[i]:
				if hstart > hend:
					hstart, hend = hend, hstart
					qseq = str(Seq(qseq).reverse_complement())
					hseq = str(Seq(hseq).reverse_complement())

				for i in range(len(hseq)):
					if hseq != "-":
						snp_count[gene_name][i-total_gap+hstart] = qseq[i]
					else:
						total_gap += 1
						x = 1
						while hseq[i-x] == "-":
							x += 1
						snp_count[gene_name][(i-total_gap+hstart)+x*0.01] = qseq[i]
	for cluster in snp_count:
		for site in snp_count[cluster]:
			try:
				pos = int(snp_pos[cluster][str(site+1)])
				base = snp_count[cluster][site]
				if base in base2num:
					profile_500[pos] = base2num[base.upper()]
				else:
					profile_500[pos] = 0
			except:
				pass

def AccProfile(name,kind,query,out_dir):

	acc_seq_file = os.path.join(os.getcwd(),"profile",name+"_"+kind+"_acc.fna")
	acc_profile = os.path.join(os.getcwd(),"profile",name+"_"+kind+"_acc.profile")
	acc_pos = {}
	with open(acc_profile,"r") as f:
		for line in f.readlines():
			acc_pos[line.split("|")[0]] = line.split("|")[-1].strip()
			
	cmd = [os.path.join(blast_dir,"makeblastdb"),\
	"-in",acc_seq_file,"-out",os.path.join(out_dir,"acc_db"),"-dbtype","nucl"]
	subprocess.check_output(cmd)

	cmd = [os.path.join(blast_dir,"blastn"),\
	"-query",query,"-db",os.path.join(out_dir,"acc_db"),\
	"-out",os.path.join(out_dir,"acc_out"),"-outfmt","6 std slen"]
	subprocess.check_output(cmd)

	acc_count = {}
	with open(os.path.join(out_dir,"acc_out"),"r") as f:
		plots = list(csv.reader(f,delimiter="	"))
		for row in plots:
			if float(row[2]) >= 80 and (float(row[3])/float(row[12])) >= 0.80:
				acc_count[row[1]] = True

	for cluster in acc_count:
		try:
			pos = int(acc_pos[cluster])
			profile_500[pos] = 5
		except:
			pass

def Prediction(name,kind,feature_list,out_dir):
	ARP = {
	"AMC":["susceptible","resistant"],
	"AMP":["susceptible","resistant"],
	"AZI":["susceptible","resistant"],
	"FOX":["susceptible","resistant"],
	"TIO":["susceptible","resistant"],
	"CRO":["susceptible","resistant"],
	"CHL":["susceptible","resistant"],
	"CIP":["susceptible","resistant"],
	"GEN":["susceptible","resistant"],
	"KAN":["susceptible","resistant"],
	"NAL":["susceptible","resistant"],
	"STR":["susceptible","resistant"],
	"FIS":["susceptible","resistant"],
	"TET":["susceptible","resistant"],
	"SXT":["susceptible","resistant"],
	}

	MIC = {
	"AMC":[1.0*(2**i) for i in range(7)],
	"AMP":[1.0*(2**i) for i in range(7)],
	"AZI":[0.25]+[1.0*(2**i) for i in range(6)],
	"FOX":[1.0*(2**i) for i in range(7)],
	"TIO":[0.25*(2**i) for i in range(7)],
	"CRO":[0.25*(2**i) for i in range(10)],
	"CHL":[2.0*(2**i) for i in range(6)],
	"CIP":[0.015625*(2**i) for i in range(9)],
	"GEN":[0.25*(2**i) for i in range(8)],
	"KAN":[8.0*(2**i) for i in range(5)],
	"NAL":[1.0*(2**i) for i in range(8)],
	"STR":[2.0*(2**i) for i in range(8)],
	"FIS":[16.0*(2**i) for i in range(6)],
	"TET":[4.0*(2**i) for i in range(5)],
	"SXT":[0.125*(2**i) for i in range(7)],
	}
	find_pheno = {"ARP":ARP,"MIC":MIC}
	x = np.array([profile_500])
	clf = pickle.load(open(os.path.join(os.getcwd(),"model","model_"+name+"_"+kind),'rb'))
	y_pred_prob = clf.predict_proba(x)
	y_pred = np.argmax(y_pred_prob)
	now_pheno = find_pheno[kind][name][y_pred]
	now_prob = y_pred_prob[0][y_pred]
	if kind == "ARP":
		f_out.write(name+"	"+kind+"	"+str(now_pheno[0].upper())\
			+"	probability	"+str(now_prob)+"\n")
		print(name,kind,now_pheno[0].upper(),"probability:",now_prob)
	else:
		now_index = y_pred
		total_prob = now_prob
		if now_index == 0:
			total_prob += y_pred_prob[0][now_index+1]
			plus_pheno = find_pheno[kind][name][now_index+1]
			pheno_range = str(now_pheno)+"-"+str(plus_pheno)
			f_out.write(name+"	"+kind+"	"+str(now_pheno)\
			+"	probability	"+str(now_prob)+"	"+pheno_range\
			+"	"+str(total_prob)+"\n")
		elif now_index == len(y_pred_prob[0])-1:
			total_prob += y_pred_prob[0][now_index-1]
			minus_pheno = find_pheno[kind][name][now_index-1]
			pheno_range = str(minus_pheno)+"-"+str(now_pheno)
			f_out.write(name+"	"+kind+"	"+str(now_pheno)\
			+"	probability	"+str(now_prob)+"	"+pheno_range\
			+"	"+str(total_prob)+"\n")
		else:
			total_prob += y_pred_prob[0][now_index-1]
			total_prob += y_pred_prob[0][now_index+1]
			minus_pheno = find_pheno[kind][name][now_index-1]
			plus_pheno = find_pheno[kind][name][now_index+1]
			pheno_range = str(minus_pheno)+"-"+str(plus_pheno)
			f_out.write(name+"	"+kind+"	"+str(now_pheno)\
			+"	probability	"+str(now_prob)+"	"+pheno_range\
			+"	"+str(total_prob)+"\n")
		print(name,kind,now_pheno,"probability:",now_prob,pheno_range,total_prob)

if __name__ == "__main__":

	cmds = sys.argv
	file_name = None
	file_dir = None
	try:
		file_name = GetParams("-inf",cmds)
	except:
		file_dir = GetParams("-ind",cmds)
	out_dir = GetParams("-out",cmds)
	anti = GetParams("-anti",cmds)
	pheno = GetParams("-pheno",cmds)

	if bool(file_name):
		f_out = open(os.path.join(out_dir,file_name.split("/")[-1])+".tsv","w")
		if anti.upper() != "A" and pheno.upper() != "A":
			profile_500 = [0 for i in range(500)]
			CoreProfile(anti,pheno,file_name,out_dir)
			AccProfile(anti,pheno,file_name,out_dir)
			Prediction(anti,pheno,profile_500,out_dir)
		elif anti.upper() == "A" and pheno.upper() != "A":
			for name in names:
				profile_500 = [0 for i in range(500)]
				CoreProfile(name,pheno,file_name,out_dir)
				AccProfile(name,pheno,file_name,out_dir)
				Prediction(name,pheno,profile_500,out_dir)
		elif anti.upper() != "A" and pheno.upper() == "A":
			for kind in kinds:
				profile_500 = [0 for i in range(500)]
				CoreProfile(anti,kind,file_name,out_dir)
				AccProfile(anti,kind,file_name,out_dir)
				Prediction(anti,kind,profile_500,out_dir)
		elif anti.upper() == "A" and pheno.upper() == "A":
			for name in names:
				for kind in kinds:
					profile_500 = [0 for i in range(500)]
					CoreProfile(name,kind,file_name,out_dir)
					AccProfile(name,kind,file_name,out_dir)
					Prediction(name,kind,profile_500,out_dir)
		f_out.close()

	elif bool(file_dir):
		files = []
		for file in os.listdir(file_dir):
			if file.endswith(".fna"):
				files.append(file)
		for file in files:
			file_name = os.path.join(file_dir,file)
			f_out = open(os.path.join(out_dir,file_name.split("/")[-1])+".tsv","w")
			if anti.upper() != "A" and pheno.upper() != "A":
				profile_500 = [0 for i in range(500)]
				CoreProfile(anti,pheno,file_name,out_dir)
				AccProfile(anti,pheno,file_name,out_dir)
				Prediction(anti,pheno,profile_500,out_dir)
			elif anti.upper() == "A" and pheno.upper() != "A":
				for name in names:
					profile_500 = [0 for i in range(500)]
					CoreProfile(name,pheno,file_name,out_dir)
					AccProfile(name,pheno,file_name,out_dir)
					Prediction(name,pheno,profile_500,out_dir)
			elif anti.upper() != "A" and pheno.upper() == "A":
				for kind in kinds:
					profile_500 = [0 for i in range(500)]
					CoreProfile(anti,kind,file_name,out_dir)
					AccProfile(anti,kind,file_name,out_dir)
					Prediction(anti,kind,profile_500,out_dir)
			elif anti.upper() == "A" and pheno.upper() == "A":
				for name in names:
					for kind in kinds:
						profile_500 = [0 for i in range(500)]
						CoreProfile(name,kind,file_name,out_dir)
						AccProfile(name,kind,file_name,out_dir)
						Prediction(name,kind,profile_500,out_dir)
			f_out.close()

	# del_file = ["acc_db","acc_out","core_db","core_out","core_align"]
	# for file in os.listdir(out_dir):
	# 	if file.split(".")[0] in del_file:
	# 		os.remove(os.path.join(out_dir,file))
