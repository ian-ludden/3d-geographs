import sys

"""
Convert a Voro++ output file, 
stored as a CSV but with tuples 
for faces and vertices, into 
a semicolon-separated file. 

Parameters:
===============
    in_fname - input filename (columns are comma-separated)
    out_fname - output filename (columns will be semicolon-separated)
"""
def convert_commas_to_semicolons(in_fname, out_fname):
    with open(in_fname, 'r') as in_f:
        with open(out_fname, 'w') as out_f:
            for line in in_f:
                first_three_cols = line.split('(')[0]
                ftc_semicol = first_three_cols.replace(',', ';')
                
                rest_of_line = line[len(first_three_cols):]
                rol_semicol = rest_of_line.replace('),(', ');(')
                
                out_f.write(ftc_semicol + rol_semicol)
    return

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: python convert_commas_to_semicolons.py [in_fname] [out_fname]')
        exit(1)
    
    in_fname = str(sys.argv[1])
    out_fname = str(sys.argv[2])

    convert_commas_to_semicolons(in_fname=in_fname, out_fname=out_fname)