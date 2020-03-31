# kj(cambala) - kj(Эталон)

import sys
import cmath

depths = []

def read_kj(filename : str):
	kj = []
	with open(filename, 'r') as file:
		lst = file.readlines()
		for x in lst:
				x = x.replace('\n', '')
				x = x.replace('\r', '')
				x = x.replace('\t', ' ')
				if x == '' or not x[0].isdigit():
						continue
				sub_lst = x.split(' ')
				if len(sub_lst) == 1:
						x.replace('\t', '')
						kj.append(float(x))
				elif len(sub_lst) > 1:
						sub_kj = []
						#print(x)
						depths.append(int(sub_lst[0]))
						for y in sub_lst[1:]:
								if y=='' or not y[0].isdigit():
										continue
								y = y.replace('i','j')
								#print(y)
								sub_kj.append(complex(y))
						#print(sub_kj)
						kj.append(sub_kj)
	return kj

kj_reference_file_name = sys.argv[1]
kj_cambala_file_name = sys.argv[2]
print('kj_reference_file_name : ' + kj_reference_file_name)
print('kj_cambala_file_name : ' + kj_cambala_file_name)

kj_reference = read_kj(kj_reference_file_name)
kj_cambala = read_kj(kj_cambala_file_name)

#print(depths)

with open('ekj.txt', 'w') as ofile:
	for i in range(len(kj_reference)):
			if type(kj_reference[i]) is float:
					ofile.write(str(abs(kj_cambala[i] - kj_reference[i])/abs(kj_reference[i])) + '\n')
			elif type(kj_reference[i]) is list:
					ofile.write(str(depths[i]) + '\t')
					for j in range(len(kj_cambala[i])):
							tmp_str = str(abs(kj_cambala[i][j] - kj_reference[i][j])/abs(kj_reference[i][j]))
							tmp_str = tmp_str.replace('(','')
							tmp_str = tmp_str.replace(')','')
							tmp_str = tmp_str.replace('j','i')
							ofile.write(tmp_str + ' ')
					ofile.write('\n')