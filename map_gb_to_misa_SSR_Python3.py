__author__ = 'Kinggerm'

# change the output of function read_annotation_of_gb(): CDS first, then intron, then tRNA, then rRNA, neglect gene.
# add column region type

# add statistics
# implement excel output

import time
import xlwt


def write_excel(this_matrixes, sheet_names, this_dir):
    f = xlwt.Workbook()
    count_sheet = 0
    for this_matrix in this_matrixes:
        this_sheet = f.add_sheet(sheet_names[count_sheet], cell_overwrite_ok=True)
        for i in range(len(this_matrix)):
            for j in range(len(this_matrix[i])):
                this_sheet.write(i, j, this_matrix[i][j])
        count_sheet += 1
    f.save(this_dir)

def get_parentheses_pairs(tree_string, sign=('(', ')')):
    tree_line = list(tree_string)
    left_sign = []
    count_ls = 0
    left_level = []
    right_sign = []
    count_rs = 0
    right_level = []
    if tree_line.count(sign[0]) == tree_line.count(sign[1]):
        for i in range(0, len(tree_line)):
            if tree_line[i] == sign[0]:
                left_sign.append(i)
                left_level.append(count_ls-count_rs)
                count_ls += 1
                right_sign.append(int)
                continue
            if tree_line[i] == sign[1]:
                count_rs += 1
                right_level.append(count_ls-count_rs)
                for j in range(0, count_ls):
                    if left_level[count_ls-1-j] == right_level[count_rs-1]:
                        right_sign[count_ls-1-j] = i
                        break
    else:
        raise Exception('Unbalanced signs in tree description: '+''.join(map(str, tree_line)))
    return left_sign, right_sign
# get_parentheses_pairs('th((i))s')
# get_parentheses_pairs('this')
# if long string, make this dict


def read_annotation_of_gb(annotation_lines, seq_len):
    # annotation_list[0] will not be utilized
    # annotation_list[1] means the dictionary for the 1st base
    # example for annotationlist[x]: [{'transl_table': '11', 'protein_id': '"BAK69443.1"', 'db_xref': '"GI:345433621"', 'product': '"hypothetical chloroplast RF2"', 'codon_start': '1', 'type': 'CDS', 'direction': 'reverse', 'gene': '"ycf2"'}, {'type': 'gene', 'direction': 'reverse', 'gene': '"ycf2"'}]
    annotation_list = [[] for x in range(0, seq_len+1)]
    #
    regions = []
    for i in range(0, len(annotation_lines)):
        if annotation_lines[i][0] in ['gene', 'tRNA', 'rRNA', 'exon', 'intron', 'CDS'] and 'gene=' in ''.join(annotation_lines[i][2:]):
            # delete join() and order()
            while 'join' in annotation_lines[i][1]:
                join_loc = annotation_lines[i][1].index('join')
                pairs = get_parentheses_pairs(annotation_lines[i][1])
                right = pairs[1][pairs[0].index(join_loc+4)]
                annotation_lines[i][1] = annotation_lines[i][1][:right]+annotation_lines[i][1][right+1:]
                annotation_lines[i][1] = annotation_lines[i][1][:join_loc]+annotation_lines[i][1][join_loc+5:]
            while 'order' in annotation_lines[i][1]:
                order_loc = annotation_lines[i][1].index('order')
                pairs = get_parentheses_pairs(annotation_lines[i][1])
                right = pairs[1][pairs[0].index(order_loc+5)]
                annotation_lines[i][1] = annotation_lines[i][1][:right]+annotation_lines[i][1][right+1:]
                annotation_lines[i][1] = annotation_lines[i][1][:order_loc]+annotation_lines[i][1][order_loc+6:]
            # assume all parentheses left are complements
            # find out all complements
            annotation_lines[i][1] = annotation_lines[i][1].replace('complement', '')
            complements = get_parentheses_pairs(annotation_lines[i][1])
            # find out all region pairs
            ellipsis = {}
            for l in range(0, len(annotation_lines[i][1])-1):
                if annotation_lines[i][1][l:l+2] == '..':
                    ellipsis[l] = 'forward'
                    # and find out the directions according to the position of complement tag
                    for m in range(0, len(complements[0])):
                        if complements[0][m] < l < complements[1][m]:
                            ellipsis[l] = 'reverse'
                            break
            ellipsis_key = [x for x in ellipsis]
            ellipsis_key.sort()
            # read into regions
            region_parts = [[int(y.replace('>', '').replace('<', '')) for y in x.split('..')] for x in annotation_lines[i][1].replace('(', '').replace(')', '').split(',')]
            for l in range(0, len(region_parts)):
                this_dict = {'type': annotation_lines[i][0], 'direction': ellipsis[ellipsis_key[l]]}
                for j in range(2, len(annotation_lines[i])):
                    this_dict[annotation_lines[i][j].split('=')[0]] = '='.join(annotation_lines[i][j].split('=')[1:])
                regions.append([region_parts[l][0], region_parts[l][1], this_dict])
                for base in range(region_parts[l][0], region_parts[l][1]+1):
                    annotation_list[base].append(this_dict)
    # regions.sort(key=lambda x:x[0])
    # fill the gaps
    names = {}
    gene_regions = []
    gene_regions += [x for x in regions if x[2]['type'] == 'CDS']
    for x in gene_regions:
        names[x[2]['gene']] = 0
    gene_regions += [x for x in regions if x[2]['type'] == 'exon' and x[2]['gene'] not in names]
    for x in gene_regions:
        names[x[2]['gene']] = 0
    gene_regions += [x for x in regions if x[2]['type'] == 'tRNA' and x[2]['gene'] not in names]
    for x in gene_regions:
        names[x[2]['gene']] = 0
    gene_regions += [x for x in regions if x[2]['type'] == 'rRNA' and x[2]['gene'] not in names]
    for x in gene_regions:
        names[x[2]['gene']] = 0
    gene_regions += [x for x in regions if x[2]['type'] == 'gene' and x[2]['gene'] not in names]
    gene_regions.sort(key=lambda x:x[0])
    gene_regions_last = gene_regions[:]
    gene_regions_last.sort(key=lambda x:x[1])
    # fill the start with noncoding
    if gene_regions[0][0] > 1:
        this_dict = {'type': 'noncoding', 'direction': 'none', 'gene': gene_regions_last[-1][2]['gene']+'--'+gene_regions[0][2]['gene']}
        regions.append([1, gene_regions[0][0]-1, this_dict])
        for base in range(1, gene_regions[0][0]):
            annotation_list[base].append(this_dict)
    # fill the end with noncoding
    if gene_regions_last[-1][1] < seq_len:
        this_dict = {'type': 'noncoding', 'direction': 'none', 'gene': gene_regions_last[-1][2]['gene']+'--'+gene_regions[0][2]['gene']}
        regions.append([gene_regions_last[-1][1]+1, seq_len, this_dict])
        for base in range(gene_regions_last[-1][1]+1, seq_len+1):
            annotation_list[base].append(this_dict)
    del gene_regions_last
    # fill the middle with noncoding
    for i in range(0, len(gene_regions)-1):
        if gene_regions[i][1] < gene_regions[i+1][0]-1:
            this_dict = {'type': 'noncoding', 'direction': 'none', 'gene': gene_regions[i][2]['gene']+'--'+gene_regions[i+1][2]['gene']}
            regions.append([gene_regions[i][1]+1, gene_regions[i+1][0]-1, this_dict])
            for base in range(gene_regions[i][1]+1, gene_regions[i+1][0]):
                annotation_list[base].append(this_dict)
    return annotation_list


def read_gb(gb_dir):
    gb_file = [x.strip('\n') for x in open(gb_dir, 'rU').readlines()]
    i = 0
    gb_structure = {}
    while i < len(gb_file):
        gb_structure[gb_file[i].split('  ')[0]] = {'description':'  '.join([x.strip() for x in gb_file[i].split('  ')[1:] if x])}
        j = i + 1
        # special
        if gb_file[i].split('  ')[0] in ['FEATURES', 'features']:
            gb_structure[gb_file[i].split('  ')[0]]['Annotations lines'] = []
        elif gb_file[i].split('  ')[0] in ['LOCUS', 'locus']:
            locus_line = [x for x in gb_file[i].split(' ') if x]
            gb_structure[gb_file[i].split('  ')[0]]['sequence length'] = int(locus_line[locus_line.index('bp')-1])
        # batch reading
        while j < len(gb_file) and gb_file[j].startswith(' '):
            #
            blank_len = 0
            while blank_len < len(gb_file[j]):
                if gb_file[j][blank_len] != ' ':
                    break
                blank_len += 1
            if gb_file[i].split('  ')[0] not in ['FEATURES', 'features', 'ORIGIN', 'origin', 'BASE COUNT', 'base count']:
                temp = [x for x in gb_file[j].split(' ') if x]
                this_title = temp[0]
                this_content = [' '.join(temp[1:])]
                while j + 1 < len(gb_file) and gb_file[j+1].startswith(' '*(blank_len+4)):
                    this_content.append(gb_file[j+1].strip())
                    j += 1
                gb_structure[gb_file[i].split(' ')[0]][this_title] = this_content
            elif gb_file[i].split('  ')[0] in ['FEATURES', 'features']:
                temp = [x for x in gb_file[j].split(' ') if x]
                while j + 1 < len(gb_file) and gb_file[j+1].startswith(' '*(blank_len+4)):
                    if gb_file[j+1].strip().startswith('/'):
                        temp.append(gb_file[j+1].strip()[1:])
                    else:
                        temp[-1] += gb_file[j+1].strip()
                    j += 1
                gb_structure[gb_file[i].split('  ')[0]]['Annotations lines'].append(temp)
            else:
                #
                pass
            j += 1
        i = j
    return gb_structure['LOCUS']['sequence length'], gb_structure['FEATURES']['Annotations lines']

        
def main():
    # ------------------read GenBank file
    time0 = time.time()
    gb_dir = input("Please input the genbank file:").strip().strip('"').strip('\'')
    time1 = time.time()
    # annotation_lines = [type, base_range, gene, other_annotations...]
    seq_len, annotation_lines = read_gb(gb_dir)
    # the (x)th element in annotation_list represent the annotation for (x)th base
    annotation_list = read_annotation_of_gb(annotation_lines, seq_len)
    # create statistics for output
    statistics = True
    if statistics:
        statistics_title = ['region type', 'num. of SSR', 'SSR bases', 'total bases', 'SSR bases/total bases']
        statistics = {'IGS':[ 0, 0, 0, '-'], 'CDS':[0, 0, 0, '-'], 'intron':[0, 0, 0, '-'], 'tRNA':[0, 0, 0, '-'], 'rRNA':[0, 0, 0, '-'], 'gene':[0, 0, 0, '-'], 'Total':[0, 0, seq_len, '-']}
        for j in range(seq_len):
            these_types = [k['type'] for k in annotation_list[j]]
            if 'noncoding' in these_types:
                temp = annotation_list[j][these_types.index('noncoding')]['gene'].replace('"', '').split('--')
                if temp[0] == temp[1]:
                    statistics['intron'][2] += 1
                else:
                    statistics['IGS'][2] += 1
            else:
                for l in ['CDS', 'tRNA', 'rRNA']:
                    if l in these_types:
                        statistics[l][2] += 1
                        break
    # ---------------------------------------------------------------------------
    # ------------------read misa tab file
    time2 = time.time()
    misa_dir = input('Please input the *.misa file:').strip().strip('"').strip('\'')
    time3 = time.time()
    misa_tab = [[y.strip() for y in x.strip().split('\t')] for x in open(misa_dir, 'rU')]
    title_line, context = misa_tab[0], misa_tab[1:]
    which_start = title_line.index('start')
    which_end = title_line.index('end')
    title_line += ['region', 'region type']
    # ---------------------------------------------------------------------------
    # ------------------calculate the region name by lines (by context)
    SSR_sites = []
    for i in range(0, len(context)):
        # keep region (gene) names and types
        gene_name_list = []
        for j in range(int(context[i][which_start]), int(context[i][which_end])+1):
            SSR_sites.append(j)
            these_types = [k['type'] for k in annotation_list[j]]
            if 'noncoding' in these_types:
                temp = annotation_list[j][these_types.index('noncoding')]['gene'].replace('"', '').split('--')
                if temp[0] == temp[1]:
                    if statistics:
                        statistics['intron'][1] += 1
                    gene_name_list.append((temp[0]+' intron', 'intron'))
                else:
                    if statistics:
                        statistics['IGS'][1] += 1
                    gene_name_list.append((annotation_list[j][these_types.index('noncoding')]['gene'].replace('"', ''), 'IGS'))
            else:
                for l in ['CDS', 'tRNA', 'rRNA', 'exon', 'gene']:
                    if l in these_types:
                        gene_name_list.append((annotation_list[j][these_types.index(l)]['gene'].replace('"', ''), annotation_list[j][these_types.index(l)]['type'].replace('"', '')))
                        if statistics:
                            statistics[l][1] += 1
                        break
        # remove duplicate region names
        gene_name_list = list(set(gene_name_list))
        if statistics:
            for l in [x[1] for x in gene_name_list]:
                statistics[l][0] += 1
        # to avoid random sequence of name
        gene_name_list.sort()
        # add region (gene) name for this line to the end of this line.
        context[i].append(', '.join([x[0] for x in gene_name_list]))
        context[i].append(', '.join([x[1] for x in gene_name_list]))
    # ---------------------------------------------------------------------------
    # ------------------write the results
    misa_matrix1 = [title_line]+context
    # open(misa_dir+'.new.tab', 'wb').writelines([('\t'.join(x)+'\n').encode('utf-8') for x in misa_matrix1])
    misa_matrix2 = [[]]
    if statistics:
        del statistics['gene']
        misa_matrix2 = []
        statistics['Total'][0] = len(context)
        statistics['Total'][1] = len(set(SSR_sites))
        for term in statistics:
            if statistics[term][2]:
                statistics[term][3] = statistics[term][1]/statistics[term][2]
        for term in statistics:
            misa_matrix2.append([term] + statistics[term])
        misa_matrix2.sort(key=lambda x:x[1])
        misa_matrix2 = [statistics_title] + misa_matrix2
    write_excel([misa_matrix1, misa_matrix2], ('misa_result', 'statistics'), misa_dir+'.new.xls')
    # print time
    time4 = time.time()
    print('Manu-typing Cost:', time1-time0+time3-time2, 's')
    print('Calculating Cost:', time4-time0-(time1-time0)-(time3-time2), 's')

if __name__ == '__main__':
    main()

