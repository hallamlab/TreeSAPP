import sys

fasta_in = sys.argv[1]
fasta_out = fasta_in.rsplit('.', 1)[0] + '.cleaned.fasta'
with open(fasta_in, 'r') as f_in:
	data = f_in.readlines()

with open(fasta_out, 'w') as f_out:
	for line in data:
		if ">" in line:
			split_line = line.lstrip('>').split('.', 1)
			new_line = '>eggnog|' + split_line[0] + '|' + split_line[1]
			f_out.write(new_line)
		else:
			f_out.write(line)

