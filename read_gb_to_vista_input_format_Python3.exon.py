__author__ = 'Kinggerm'

import time
import os
import platform

this_dir_split = '/'
if 'Win' in platform.architecture()[1]:
    this_dir_split = '\\'


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


def read_annotation_of_gb(annotation_lines, seq_len, by_site=True):
    # annotation_list[0] will not be utilized
    # annotation_list[1] means the dictionary for the 1st base
    # example for annotationlist[x]: [{'transl_table': '11', 'protein_id': '"BAK69443.1"', 'db_xref': '"GI:345433621"', 'product': '"hypothetical chloroplast RF2"', 'codon_start': '1', 'type': 'CDS', 'direction': 'reverse', 'gene': '"ycf2"'}, {'type': 'gene', 'direction': 'reverse', 'gene': '"ycf2"'}]
    if by_site:
        annotation_list = [[] for x in range(0, seq_len+1)]
    # parse lines
    regions = []
    error_in_gb = []
    for i in range(0, len(annotation_lines)):
        if annotation_lines[i][0] in ['gene', 'tRNA', 'rRNA', 'exon', 'intron', 'CDS'] and not 'gene=' in ''.join(annotation_lines[i][2:]):
            error_in_gb.append('\t'.join(annotation_lines[i]))
        elif annotation_lines[i][0] in ['gene', 'tRNA', 'rRNA', 'exon', 'intron', 'CDS'] and 'gene=' in ''.join(annotation_lines[i][2:]):
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
                    this_dict[annotation_lines[i][j].split('=')[0]] = '='.join(annotation_lines[i][j].split('=')[1:]).strip('"').strip('\'')
                regions.append([region_parts[l][0], region_parts[l][1], this_dict])
                if by_site:
                    for base in range(region_parts[l][0], region_parts[l][1]+1):
                        annotation_list[base].append(this_dict)
    if error_in_gb:
        return {'error': error_in_gb}
    else:
        # locate region that occupied by genes
        names = {}
        gene_regions = []
        gene_regions += [x for x in regions if x[2]['type'] == 'CDS']
        for x in gene_regions:
            names[x[2]['gene']] = 0
        gene_regions += [x for x in regions if x[2]['type'] == 'tRNA' and x[2]['gene'] not in names]
        for x in gene_regions:
            names[x[2]['gene']] = 0
        gene_regions += [x for x in regions if x[2]['type'] == 'rRNA' and x[2]['gene'] not in names]
        for x in gene_regions:
            names[x[2]['gene']] = 0
        gene_regions += [x for x in regions if x[2]['type'] == 'intron' and x[2]['gene'] not in names]
        gene_regions.sort(key=lambda x:x[0])
        gene_regions_last = gene_regions[:]
        gene_regions_last.sort(key=lambda x:x[1])
        # fill the start with IGS
        if gene_regions[0][0] > 1:
            this_dict = {'type': 'IGS', 'direction': 'none', 'gene': gene_regions_last[-1][2]['gene']+'--'+gene_regions[0][2]['gene']}
            regions.append([1, gene_regions[0][0]-1, this_dict])
            if by_site:
                for base in range(1, gene_regions[0][0]):
                    annotation_list[base].append(this_dict)
        # fill the end with IGS
        if gene_regions_last[-1][1] < seq_len:
            this_dict = {'type': 'IGS', 'direction': 'none', 'gene': gene_regions_last[-1][2]['gene']+'--'+gene_regions[0][2]['gene']}
            regions.append([gene_regions_last[-1][1]+1, seq_len, this_dict])
            if by_site:
                for base in range(gene_regions_last[-1][1]+1, seq_len+1):
                    annotation_list[base].append(this_dict)
        del gene_regions_last
        # fill the middle with IGS
        for i in range(0, len(gene_regions)-1):
            if gene_regions[i][1] < gene_regions[i+1][0]-1:
                this_dict = {'type': 'IGS', 'direction': 'none', 'gene': gene_regions[i][2]['gene']+'--'+gene_regions[i+1][2]['gene']}
                regions.append([gene_regions[i][1]+1, gene_regions[i+1][0]-1, this_dict])
                if by_site:
                    for base in range(gene_regions[i][1]+1, gene_regions[i+1][0]):
                        annotation_list[base].append(this_dict)
        regions.sort(key=lambda x:(x[0], x[1]))
        if by_site:
            return {'by_region': regions, 'by_site': annotation_list}
        else:
            return {'by_region': regions}


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


def write_vista_line(regions, which, out, repeat_exon=False, start_with_direction=True):
    if not start_with_direction:
        this_direction = ''
        this_name = '\texon'
    elif regions[0][2]['direction'] == 'forward':
        this_direction = '>\t'
        this_name = '\t'+regions[which][2]['gene']
    else:
        this_direction = '<\t'
        this_name = '\t'+regions[which][2]['gene']
    out.append(this_direction+str(regions[which][0])+'\t'+str(regions[which][1])+this_name)
    if repeat_exon:
        out.append(str(regions[which][0])+'\t'+str(regions[which][1])+'\texon')
    del regions[which]


def vista_formate_with_gb(this_gb_dir):
    # annotation_lines = [type, base_range, gene, other_annotations...]
    seq_len, annotation_lines = read_gb(this_gb_dir)
    # the (x)th element in annotation_list represent the annotation for (x)th base
    annotation_dict = read_annotation_of_gb(annotation_lines, seq_len, False)
    if not 'error' in annotation_dict:
        region_list = annotation_dict['by_region']
        # ---------------------------------------------------------------------------
        # ------------------delete unrelated regions
        # for i in range(len(region_list)):
        #     print(region_list[i])
        i = 0
        filter_type = {'exon', 'gene', 'tRNA', 'rRNA'}
        while i < len(region_list):
            if region_list[i][2]['type'] not in filter_type or (i > 0 and region_list[i][0] == region_list[i-1][0] and region_list[i][1] == region_list[i-1][1] and region_list[i][2]['gene'] == region_list[i-1][2]['gene']):
                del region_list[i]
            else:
                i += 1
        # ---------------------------------------------------------------------------
        # ------------------prepare the result format
        out_put = []

        while region_list:
            if len(region_list) > 1:
                if not region_list[0][2]['gene'] == region_list[1][2]['gene']:
                    write_vista_line(region_list, 0, out_put, True, True)
                else:
                    this_gene = region_list[1]
                    write_vista_line(region_list, 1, out_put, False, True)
                    write_vista_line(region_list, 0, out_put, False, False)
                    i = 0
                    while region_list and region_list[i][1] <= this_gene[1]:
                        if region_list[i][2]['gene'] == this_gene[2]['gene']:
                            write_vista_line(region_list, i, out_put, False, False)
                        else:
                            i += 1
            else:
                write_vista_line(region_list, 0, out_put)
        # ---------------------------------------------------------------------------
        # ------------------write the results
        open(this_gb_dir+'.for_vista.txt', 'wb').writelines([(x+'\n').encode('utf-8') for x in out_put])
    else:
        print('ERROR in reading '+this_gb_dir+':')
        for line in annotation_dict['error']:
            print('\t'+line)


def main():
    # ------------------read GenBank file
    time0 = time.time()
    time_manu = 0
    print("\nEnter 'q' or 'quit' or 'exit' to end this script.\n")
    try:
        while True:
            time1 = time.time()
            input_dir = input("Please input a genbank file(*.gb) or a directory:").strip()
            time_manu += time.time()-time1
            if not input_dir:
                continue
            elif os.path.isfile(input_dir):
                vista_formate_with_gb(input_dir)
            elif os.path.isdir(input_dir):
                gb_files = [input_dir.rstrip(this_dir_split)+this_dir_split+x for x in os.listdir(input_dir) if x.endswith('.gb')]
                for gb_file in gb_files:
                    vista_formate_with_gb(gb_file)
            elif input_dir in {'q', 'quit', 'exit'}:
                break
    except KeyboardInterrupt:
        time_manu += time.time()-time1
    # print time
    time2 = time.time()
    print('\nManu-typing Cost:', time_manu, 's')
    print('Calculating Cost:', time2-time0-time_manu, 's')


if __name__ == '__main__':
    main()
