import os
import re
import sys
import argparse
import tempfile
import shutil

def doloop2d_treat(result,lines,file_out_main):

	keywords = ["CD", "DC", "CF", "FFC", "FC", "dZ", "dR", "BC", "wrk1", "wrk2","my_kbl","Cr"]
	replacements = {}
	for keyword in keywords:
		pattern = rf"{keyword}\(i,([^)]*)\)"
		replacement = f"{keyword}1D"
		lines, count = re.subn(pattern, fr"{replacement}(\1)", lines)
		if count == 0:
			pattern = rf"{keyword}\(i\)"
			replacement = f"{keyword}1D"
			lines, count = re.subn(pattern, fr"{replacement}", lines)			
		replacements[replacement] = count

	file_out_main.write(f"      DO j={result.group(3)},{result.group(4)}\n")
    
	private_vars = [var for var, count in replacements.items() if count > 0]
    
	if private_vars:
		accprivate = "!$acc loop private(" + ",".join(private_vars) + ") vector\n"
		file_out_main.write(accprivate)

	file_out_main.write("        DO i={},{}\n".format(result.group(1), result.group(2)))

	pattern_doloop = re.compile(
		r"^[ \t]*do[ \t]*i[ \t]*=[ \t]*{}[ \t]*,[ \t]*{}.*".format(result.group(1), result.group(2)),
		re.IGNORECASE
	)
	pattern_enddo = re.compile(r"^[ \t]*enddo.*", re.IGNORECASE)

	in_loop = False

	for line in lines.splitlines(keepends=True):
		if not in_loop and pattern_doloop.match(line):
			in_loop = True
			continue
		if in_loop and pattern_enddo.match(line):
			in_loop = False
			continue
		file_out_main.write(line)

	file_out_main.write("      ENDDO\n")
	file_out_main.write("      ENDDO\n")

def doextend_treat(result,lines,file_out_main):
	file_out_main.write("      DO %s=%s,%s\n" % (result.group(1),result.group(2),result.group(3)))
	lines = re.sub(r"FX\(([^)]*),([^)]*)\)",r"FX_3D(\1,\2,k)",lines)
	arrays=result.group(4).split(',')
	for array in arrays:
		replace_string_in=r"%s\(([^)]*),([^)]*)\)"  % (array)
		replace_string_out=r"%s_3D(\1,\2,k)" % (array)
		lines = re.sub(replace_string_in,replace_string_out,lines)
	file_out_main.write(lines)
	file_out_main.write("      ENDDO\n")
	(lines,cd1d) = re.subn(r"CD\(i,([^)]*)\)",r"CD1D(\1)",lines)

if __name__=="__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('outfile')
	args=parser.parse_args()
	input_file = open(args.infile,'r',encoding='iso-8859-1')
	output_file=open(args.outfile,'w',encoding='iso-8859-1')

	file_content =''.join(input_file.readlines())
	# Compile regular expressions for detecting DOLOOP2D and its end
	doloop2d_pattern = re.compile(r"^[ \t]*DOLOOP2D\((.*),(.*),(.*),(.*)\).*", re.IGNORECASE)
	end_doloop2d_pattern = re.compile(r"^[ \t]*ENDDOLOOP2D.*", re.IGNORECASE)

	# Compile regular expressions for detecting DOEXTEND and its end
	doextend_pattern = re.compile(r"^[ \t]*DOEXTEND\(([^,]*),([^,]*),([^,]*),(.*)\).*", re.IGNORECASE)
	end_doextend_pattern = re.compile(r"^[ \t]*ENDDOEXTEND.*", re.IGNORECASE)

	# Assume 'file_content' contains the entire file content as a string
	lines = file_content.splitlines(keepends=True)
	line_iterator = iter(lines)

	for line in line_iterator:
		# Check if the line starts a DOLOOP2D block
		doloop2d_match = doloop2d_pattern.match(line)
		if doloop2d_match:
			block_lines = []
			# Gather all lines until the end of the DOLOOP2D block is found
			for subline in line_iterator:
				if end_doloop2d_pattern.match(subline):
					# Process the collected block with a dedicated function
					doloop2d_treat(doloop2d_match, ''.join(block_lines), output_file)
					break
				block_lines.append(subline)
			continue

		# Check if the line starts a DOEXTEND block
		doextend_match = doextend_pattern.match(line)
		if doextend_match:
			block_lines = []
			# Gather all lines until the end of the DOEXTEND block is found
			for subline in line_iterator:
				if end_doextend_pattern.match(subline):
					# Process the collected block with a dedicated function
					doextend_treat(doextend_match, ''.join(block_lines), output_file)
					break
				block_lines.append(subline)
			continue

		# Write any other line directly to the output file
		output_file.write(line)

	input_file.close()
	output_file.close()
