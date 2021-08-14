import numpy as np
from pathlib import Path
import sys

m2au = 6.684587122268445e-12

raw_elements_file = sys.argv[1]

output_names = ["{}Myr_Baps_proper_elements.txt".format(num) for num in np.arange(0, 401, 10)]

psize = Path(raw_elements_file).stem[:-7]

flags = {'Semi-Major Axis (m)':False, 'Eccentricity':False, 'Inclination (rad)':False}
raw_data = {}
keys = []
with open(raw_elements_file) as f:
	for line in f:
		if 'Input' in line:
			continue

		if 'Semi-Major Axis (m)' in line:
			flags['Semi-Major Axis (m)'] = True
			flags['Eccentricity']		 = False
			flags['Inclination (rad)']	 = False

		if 'Eccentricity' in line:
			flags['Semi-Major Axis (m)'] = False
			flags['Eccentricity']		 = True
			flags['Inclination (rad)']	 = False

		if 'Inclination (rad)' in line:
			flags['Semi-Major Axis (m)'] = False
			flags['Eccentricity'] 		 = False
			flags['Inclination (rad)'] 	 = True

		contents = line.split('\t')

		if len(contents) > 1:
			if flags['Semi-Major Axis (m)'] == True and 'p0' in contents:
				for name in contents:
					key = '{}_{}'.format(psize, name[1])
					keys.append(key)
					raw_data[key] = {'a': [], 'e':[], 'i':[]}

			elif flags['Semi-Major Axis (m)'] == True:
				for j, a in enumerate(contents):
					if j >= len(keys): continue
					raw_data[keys[j]]['a'].append(float(a))
			
			elif flags['Eccentricity'] == True and 'p0' not in contents:
				for j, e in enumerate(contents):
					if j >= len(keys): continue
					raw_data[keys[j]]['e'].append(float(e))

			elif flags['Inclination (rad)'] == True and 'p0' not in contents:
				for j, i in enumerate(contents):
					if j >= len(keys): continue
					raw_data[keys[j]]['i'].append(float(i))

for j, out in enumerate(output_names):
	cleaned_data = []
	try:
		with open(out) as f:
			for line in f:
				row = list(line.split('\t'))
				for k, d in enumerate(row):
					if k ==0: continue
					row[k] = float(d)
				cleaned_data.append(row)
	except:
		pass

	year_index = j*1000
	lb = year_index - 100 if year_index > 100 else year_index
	ub = year_index + 100 if (40000 - year_index) > 100 else year_index

	for key in keys:
		mean_a = np.mean(raw_data[key]['a'][lb:ub]) * m2au
		mean_e = np.mean(raw_data[key]['e'][lb:ub])
		mean_i = np.mean(raw_data[key]['i'][lb:ub])
		p = mean_a*(1-mean_e)

		if mean_a < 0 or mean_a > 10 or p < 1.5: continue

		more_clean_data = [key, mean_a, mean_e, mean_i]
		cleaned_data.append(more_clean_data)

	with open(out, 'w') as f:
		for row in cleaned_data:
			if len(row) != 4: continue
			line = '{}\t{:.16f}\t{:.16f}\t{:.16f}\n'.format(*row)
			f.write(line)
	

	

					
			
