import sys
import cmath

depths = []

def read_phij(filename : str):
	phij = []
	with open(filename, 'r') as file:
		lst = file.readlines()
		for x in lst:
				x = x.replace('\n', '')
				x = x.replace('\r', '')
				x = x.replace('\t', ' ')
				if x == '' or not x[0].isdigit():
						continue
				sub_lst = x.split(' ')
				if len(sub_lst) > 1:
						sub_phij = []
						depths.append(int(sub_lst[0].split('.')[0]))
						for y in sub_lst[1:]:
								if y == '' or len(y) < 2: # or not y[1].isdigit()
										continue
								if not y[0].isdigit() and not y[1].isdigit():
										continue
								y = y.replace('i','j')
								sub_phij.append(float(y))
						phij.append(sub_phij)
	return phij

phij_reference_file_name = sys.argv[1]
phij_cambala_file_name = sys.argv[2]
print('phij_reference_file_name : ' + phij_reference_file_name)
print('phij_cambala_file_name : ' + phij_cambala_file_name)

phij_reference = read_phij(phij_reference_file_name)
phij_cambala = read_phij(phij_cambala_file_name)

print(phij_reference[0])

#print(depths)

# E_j = max| \phi_j^s - \phi_{j,ref}^s  | / max|\phi_{j,ref}^s|
err_phij_lst = []
with open('ephij.txt', 'w') as ofile:
	ofile.write('mode max|phi_cambala-phi_ref|/max|phi_ref| max|phi_cambala-phi_ref| max|phi_ref| \n')
	mod_numbers = len(phij_reference[0])
	print('mod_numbers : %d' % mod_numbers)
	for i in range(mod_numbers):
		reference_mods_all_depths = [float(x[i]) for x in phij_reference]
		abs_reference_mods_all_depths = [abs(x) for x in reference_mods_all_depths]
		max_reference_mod_all_depths = max(abs_reference_mods_all_depths)
		cambala_mods_all_depths = [float(x[i]) for x in phij_cambala]
		#print('cambala_mods_all_depths len : %d' % len(cambala_mods_all_depths))
		#print('reference_mods_all_depths len : %d' % len(reference_mods_all_depths))
		#print(reference_mods_all_depths)
		#print(cambala_mods_all_depths)
		subtraction_lst = []
		for j in range(len(cambala_mods_all_depths)):
			subtraction_lst.append(abs(cambala_mods_all_depths[j]-reference_mods_all_depths[j]))
		max_subtraction = max(subtraction_lst)
		print('mode %d' % i)
		print('max_subtraction : %f' % max_subtraction)
		print('max_reference_mod_all_depths : %f' % max_reference_mod_all_depths)
		err = max_subtraction / max_reference_mod_all_depths
		err_phij_lst.append(err)
		ofile.write('%d %f %f %f \n' % (i, err, max_subtraction, max_reference_mod_all_depths))
		#print('subtraction_lst len : %d' % len(subtraction_lst))
		#print(subtraction_lst)
#print(err_phij_lst)