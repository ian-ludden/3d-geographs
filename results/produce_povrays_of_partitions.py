import numpy as np
import os
import pandas as pd

colors = np.array([\
    [0.02352941, 0.60392157, 0.95294118], \
    [0.90196078, 0.85490196, 0.65098039], \
    [0.01176471, 0.2627451, 0.8745098 ], \
    [0.98823529, 0.35294118, 0.31372549], \
    [0.98039216, 0.76078431, 0.01960784], \
    [0.78039216, 0.62352941, 0.9372549 ], \
    [0.99607843, 0.25882353, 0.05882353], \
    [0.49411765, 0.11764706, 0.61176471], \
    [0.6627451, 0.3372549, 0.11764706], \
    [0.00784314, 0.57647059, 0.5254902 ], \
    [0.98431373, 0.86666667, 0.49411765], \
    [0.73333333, 0.97647059, 0.05882353], \
    [0.75686275, 0.97254902, 0.03921569], \
    [0.01960784, 0.28627451, 0.02745098], \
    [0.57254902, 0.58431373, 0.56862745], \
    [0.21960784, 0.00784314, 0.50980392], \
    [0.48235294, 0.78431373, 0.96470588], \
    [0.76078431, 0., 0.47058824], \
    [0.43137255, 0.45882353, 0.05490196], \
    [0.97647059, 0.45098039, 0.02352941], \
    [1., 0.50588235, 0.75294118], \
    [0.89803922, 0., 0.        ], \
    [0.77254902, 0.78823529, 0.78039216], \
    [0.9372549, 0.25098039, 0.14901961], \
    [0.60392157, 0.05490196, 0.91764706]])

colors_str = ['<' + ','.join(list(color.astype(str))) + '>' for color in colors]

"""
Generate POV-Ray code for visualizing the Voronoi cell 
with the given faces and vertices. 
"""
def gen_cell_str(faces_str, vertices_str, color_rgb='<red_val, green_val, blue_val>'):
    total_str = ""
    for face_str in faces_str:
        total_str += '// ' + face_str + '\n'
        total_str += 'polygon {' + '\n'

        face = face_str[1:-1].split(',')
        total_str += '\t{},'.format(len(face)) + '\n'

        vtx_str = ""
        for vi in face:
            vtx_str += "<" + vertices_str[int(vi)][1:-1] + ">, "
        vtx_str = vtx_str[:-2] # Remove last comma
        total_str += '\t' + vtx_str + '\n'

        total_str += '\tfinish { ambient 0.42 specular 0.5 }\n'
        total_str += '\tpigment { rgb ' + color_rgb + ' }\n'
        total_str += '}\n\n'
    
    return total_str


if __name__ == '__main__':
    dirlist = os.listdir('./')
    partitions = [file[:-4] for file in dirlist if (file[:9] == 'partition' and file[-4:] == '.csv')]

    for partition in partitions:
        filename_parts = partition.split('_')
        honeycomb_type = filename_parts[1]
        S = int(filename_parts[2])
        K = int(filename_parts[3][1:])

        dfpart = pd.read_csv(f"{partition}.csv")
        dfpart = dfpart.rename(columns={'cell_id': 'ID', 'part': 'Part'})
        dfpart = dfpart.set_index('ID')

        in_fname = f"../generate_input/{honeycomb_type}/{honeycomb_type}_{S}.csv"
        in_fname_sc = f"../tmp/{honeycomb_type}_{S}_sc.csv"
        if not os.path.exists(in_fname_sc):
            raise Exception(f"Missing semicolon-separated version of {in_fname}, expected at {in_fname_sc}. Please run '../generate_input/convert_all_semicolon_versions.sh'.")
        
        df = pd.read_csv(in_fname_sc, sep=';', header=1)
        df = df.set_index('ID')
        df = df.join(dfpart)

        parts = sorted(list(np.unique(df['Part'])))

        for part in parts:
            out_fname = f"partition_{honeycomb_type}_{S}_K{K}_{part}.pov"
            with open(out_fname, 'w') as f:
                f.write('#version 3.6;\n')
                f.write('// Right-handed coordinate system in which the z-axis points upwards\n')
                f.write('camera {\n')
                f.write(f"\tlocation <{3*S}, {-6.5*S}, {4*S}>\n")
                f.write('\tsky z\n')
                f.write('\tright -0.22*x*image_width/image_height\n')
                f.write('\tup 0.22*z\n')
                f.write(f"\tlook_at <{S/2}, {S/2}, {S/2}>\n")
                f.write('\t}\n')

                f.write(f"#declare red_val={colors[part - 1][0]};\n");
                f.write(f"#declare green_val={colors[part - 1][1]};\n");
                f.write(f"#declare blue_val={colors[part - 1][2]};\n");

                f.write('// White background\n')
                f.write('background{rgb 1}\n')

                f.write('// Two lights with slightly different colors\n')
                f.write('light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}\n')
                f.write('light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}\n')
                f.write('')
                
                # Write faces of each cell in part, as polygons
                for idx in df.index:
                    if df.loc[idx, 'Part'] == part:
                        cell_str = gen_cell_str(df.loc[idx, 'Faces'].split(' '), df.loc[idx, 'Vertices'].split(' '))
                        f.write(cell_str)