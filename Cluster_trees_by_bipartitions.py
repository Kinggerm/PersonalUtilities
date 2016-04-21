import sys
import time
import os
from shlex import shlex
from types import StringType, ListType
from cStringIO import StringIO

#

# Nov 8, 2015
##########################################
def delete_prefix(tab_matrix, is_titled=False, column=0):
    import copy
    del_prefix_matrix = copy.deepcopy(tab_matrix)
    i = len(tab_matrix[is_titled][column])
    all_the_same_prefix = False
    while i > 0:
        j = is_titled + 1
        while j < len(tab_matrix):
            if tab_matrix[is_titled][column][0:i] == tab_matrix[j][column][0:i]:
                if j == len(tab_matrix)-1:
                    all_the_same_prefix = True
                    break
                j = j + 1
            else:
                break
        if all_the_same_prefix:
            break
        else:
            i = i - 1
    if all_the_same_prefix:
        for k in range(is_titled, len(tab_matrix)):
            del_prefix_matrix[k][column] = del_prefix_matrix[k][column][i:len(del_prefix_matrix[k][column])]
    return del_prefix_matrix


def delete_suffix(tab_matrix, is_titled=False, column=0):
    import copy
    del_suffix_matrix = copy.deepcopy(tab_matrix)
    i = len(del_suffix_matrix[is_titled][column])
    all_the_same_suffix = False
    while i > 0:
        j = is_titled + 1
        while j < len(del_suffix_matrix):
            if del_suffix_matrix[is_titled][column][(len(del_suffix_matrix[is_titled][column])-i):len(del_suffix_matrix[is_titled][column])] == del_suffix_matrix[j][column][(len(del_suffix_matrix[j][column])-i):len(del_suffix_matrix[j][column])]:
                if j == len(del_suffix_matrix)-1:
                    all_the_same_suffix = True
                    break
                j = j + 1
            else:
                break
        if all_the_same_suffix:
            break
        else:
            i = i - 1
    if all_the_same_suffix:
        for k in range(is_titled, len(tab_matrix)):
            del_suffix_matrix[k][column] = del_suffix_matrix[k][column][0:(len(del_suffix_matrix[k][column])-i)]
    return del_suffix_matrix


def read_tab(TabFile, Titled = False):
    file_original = open(TabFile, 'rU')
    TabRows = file_original.readlines()
    TabList = []
    for i in range(0, len(TabRows)):
        #only deal with simple tab file
        TabList.append(TabRows[i].split('\t'))
        TabList[i][-1] = TabList[i][-1].replace('\r','').replace('\n','')
        while not TabList[i][-1]:
            del TabList[i][-1]
    file_original.close()
    if not Titled:
        return TabList
    else:
        TitleLine = TabList[0]
        del TabList[0]
        return [TitleLine, TabList]


def write_tab(ListToWrite, FileName = './Extracted', WriteFormat = 'all'):
    #generate new file by order
    #if os.path.isfile(FileName):
        #FileList = os.listdir('./')
        #i=0
        #while [(FileName == FileList[i]) or (FileName+'_'+str(i+1) == FileList[i])] and (i<len(FileList)):
        #    i = i + 1
        #    if i = len(FileList):
        #          break
        #if i > 0:
        #FileName == FileName+'_'+str(i+1)
    #generate new file by time
    FileName = FileName+'_'+str(int(time.time()))
    FileGenerate = open(FileName+'.Temp', 'w+')

    #change the type of items in ListToWrite into string
    #ListToWrite = [str(ListToWrite) for ListToWrite in ListToWrite if ListToWrite]
    for i in range(0, len(ListToWrite)):
        for j in range(0, len(ListToWrite[i])):
            ListToWrite[i][j] = str(ListToWrite[i][j])

    if WriteFormat == 'last':
        for i in range(0,len(ListToWrite)):
            FileGenerate.write(ListToWrite[i][len(ListToWrite[i])-1])
            FileGenerate.write('\r\n')
    if WriteFormat == 'all':
        for i in range(0,len(ListToWrite)):
            FileGenerate.writelines('\t'.join(ListToWrite[i]))
            FileGenerate.write('\r\n')
    FileGenerate.close()
    os.rename(FileName+'.Temp', FileName+'.tab')
    print 'The file is generated:  '+FileName+'.tab'


# tree_line = ['(','name1',':','1',',','(','(','name2',',','name3',')',',','name4',')',')']
# return [[0, 5, 6], [14, 13, 10]]
def get_parentheses_pairs(tree_string, sign=['(', ')']):
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
                        right_sign[count_ls-1-j]=i
                        break
    else:
        raise StandardError('Unbalanced signs in tree description: '+''.join(map(str, tree_line)))
    return [left_sign, right_sign]


# tree_line = ['(','name1',':','1',',','(','(','name2',',','name3',')',',','name4',')',')']
# return [['name1', 'name2', 'name3', 'name4'], [1, 7, 9, 12]]
def get_taxa_from_newick(tree_line):
    taxa_names = []
    taxa_locate = []
    i = 0
    while i < len(tree_line):
        if tree_line[i] == '(':
            i = i + 1
            while tree_line[i] == '(':
                i = i + 1
            taxa_names.append(tree_line[i])
            taxa_locate.append(i)
            i = i + 1
            continue
        if tree_line[i] == ',':
            i = i + 1
            while tree_line[i] == '(':
                i = i + 1
            taxa_names.append(tree_line[i])
            taxa_locate.append(i)
            i = i + 1
            continue
        i = i + 1
    return [taxa_names, taxa_locate]


# tree_line = ['(','name1',':','0.11',',','(','(','name2',':','0.5',',','name3',')','98',':','0.12',',','name4',')',')']
# tree_line = '(name1:0.11,((name2:0.5,name3)98:0.12,name4))'
# return [[set(['name2', 'name3']), set(['name4', 'name1'])], [98, 98], [0.12, 0.12]]
# tree_line = '(Aphananthe_aspera:0.01028904103994690393,(((Chaetachme_aristata:0.00000043361344959758,Pteroceltis_tatarinowii:0.01028433168183005841)100:0.02072827102260334925,(Cannabis_sativa:0.00000043361344959758,Humulus_scandens:0.01034079306489057125)100:0.02080390161275600183)4:0.00000043361344959758,((Trema_orientalis:0.00000043361344959758,(Parasponia_rugosa:0.00000043361344959758,Gironniera_subaequalis:0.00000043361344959758)1:0.00000043361344959758)0:0.00000043361344959758,(Trema_levigata:0.00000043361344959758,Lozanella_enantiophylla:0.00000043361344959758)4:0.00000043361344959758)3:0.00000043361344959758)11:0.00000043361344959758,Celtis_biondii:0.01028904052294378133)'
# return [[set(['Trema_levigata', 'Humulus_scandens', 'Lozanella_enantiophylla', 'Gironniera_subaequalis', 'Parasponia_rugosa', 'Chaetachme_aristata', 'Pteroceltis_tatarinowii', 'Trema_orientalis', 'Cannabis_sativa']), set(['Aphananthe_aspera', 'Celtis_biondii']), set(['Pteroceltis_tatarinowii', 'Chaetachme_aristata', 'Humulus_scandens', 'Cannabis_sativa']), set(['Parasponia_rugosa', 'Celtis_biondii', 'Trema_levigata', 'Trema_orientalis', 'Aphananthe_aspera', 'Lozanella_enantiophylla', 'Gironniera_subaequalis']), set(['Pteroceltis_tatarinowii', 'Chaetachme_aristata']), set(['Trema_levigata', 'Humulus_scandens', 'Lozanella_enantiophylla', 'Gironniera_subaequalis', 'Parasponia_rugosa', 'Celtis_biondii', 'Trema_orientalis', 'Cannabis_sativa', 'Aphananthe_aspera']), set(['Cannabis_sativa', 'Humulus_scandens']), set(['Trema_levigata', 'Lozanella_enantiophylla', 'Gironniera_subaequalis', 'Parasponia_rugosa', 'Chaetachme_aristata', 'Celtis_biondii', 'Trema_orientalis', 'Pteroceltis_tatarinowii', 'Aphananthe_aspera']), set(['Trema_orientalis', 'Parasponia_rugosa', 'Lozanella_enantiophylla', 'Gironniera_subaequalis', 'Trema_levigata']), set(['Chaetachme_aristata', 'Celtis_biondii', 'Cannabis_sativa', 'Pteroceltis_tatarinowii', 'Aphananthe_aspera', 'Humulus_scandens']), set(['Trema_orientalis', 'Parasponia_rugosa', 'Gironniera_subaequalis']), set(['Trema_levigata', 'Chaetachme_aristata', 'Celtis_biondii', 'Cannabis_sativa', 'Pteroceltis_tatarinowii', 'Aphananthe_aspera', 'Humulus_scandens', 'Lozanella_enantiophylla']), set(['Parasponia_rugosa', 'Gironniera_subaequalis']), set(['Trema_levigata', 'Humulus_scandens', 'Lozanella_enantiophylla', 'Chaetachme_aristata', 'Celtis_biondii', 'Pteroceltis_tatarinowii', 'Trema_orientalis', 'Cannabis_sativa', 'Aphananthe_aspera']), set(['Trema_levigata', 'Lozanella_enantiophylla']), set(['Humulus_scandens', 'Gironniera_subaequalis', 'Parasponia_rugosa', 'Chaetachme_aristata', 'Celtis_biondii', 'Pteroceltis_tatarinowii', 'Trema_orientalis', 'Cannabis_sativa', 'Aphananthe_aspera'])], [11, 11, 4, 4, 100, 100, 100, 100, 3, 3, 0, 0, 1, 1, 4, 4], [4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07, 0.02072827102260335, 0.02072827102260335, 0.020803901612756002, 0.020803901612756002, 4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07, 4.3361344959758e-07]]
def transform_to_bipartition(tree_line):
    #[[{set1},{set2}],[bootstrap1,bootstrap2],[length1,length2]]
    bipartition_matrix = [[],[],[]]
    if not type(tree_line) is ListType:
        tree_line = parse(tree_line)
    parentheses_pairs = get_parentheses_pairs(tree_line)
    taxa_info = get_taxa_from_newick(tree_line)
    for i in range(0, len(parentheses_pairs[0])):
        #print str(parentheses_pairs[0][i])+'-'+str(parentheses_pairs[1][i])
        temp_set1 = set([])
        temp_set2 = set([])
        for j in range(0, len(taxa_info[1])):
            if taxa_info[1][j] > parentheses_pairs[0][i] and taxa_info[1][j] < parentheses_pairs[1][i]:
                #print taxa_info[0][j]+'in',
                temp_set1 = temp_set1 | {taxa_info[0][j]}
            else:
                #print taxa_info[0][j]+'out',
                temp_set2 = temp_set2 | {taxa_info[0][j]}
        #print ''
        if len(temp_set1) >= 2 and len(temp_set2) >= 2:
            bipartition_matrix[0].append(temp_set1)
            bipartition_matrix[0].append(temp_set2)
            temp_num = 0
            try:
                is_int = float(tree_line[parentheses_pairs[1][i]+1]) == int(float(tree_line[parentheses_pairs[1][i]+1]))
            except:
                bipartition_matrix[1].append(0)
                bipartition_matrix[1].append(0)
            else:
                if is_int:
                    bipartition_matrix[1].append(int(tree_line[parentheses_pairs[1][i]+1]))
                    bipartition_matrix[1].append(int(tree_line[parentheses_pairs[1][i]+1]))
                else:
                    bipartition_matrix[1].append(float(tree_line[parentheses_pairs[1][i]+1]))
                    bipartition_matrix[1].append(float(tree_line[parentheses_pairs[1][i]+1]))
                temp_num = 1
            if tree_line[parentheses_pairs[1][i]+temp_num+1] == ':':
                bipartition_matrix[2].append(float(tree_line[parentheses_pairs[1][i]+temp_num+2]))
                bipartition_matrix[2].append(float(tree_line[parentheses_pairs[1][i]+temp_num+2]))
            else:
                bipartition_matrix[2].append(0)
                bipartition_matrix[2].append(0)
    for i in range(0, len(taxa_info[0])):
        bipartition_matrix[0].append(set([taxa_info[0][i]]))
        bipartition_matrix[0].append(set(taxa_info[0])-set([taxa_info[0][i]]))
        bipartition_matrix[1].append(100)
        bipartition_matrix[1].append(100)
        if tree_line[taxa_info[1][i]+1] == ':':
            bipartition_matrix[2].append(float(tree_line[taxa_info[1][i]+2]))
            bipartition_matrix[2].append(float(tree_line[taxa_info[1][i]+2]))
        else:
            bipartition_matrix[2].append(0)
            bipartition_matrix[2].append(0)
    return [taxa_info[0], bipartition_matrix]


class Tokenizer(shlex):
    def __init__(self, infile):
        shlex.__init__(self, infile)
        self.commenters = ''
        self.wordchars = self.wordchars+'-.'
        self.quotes = "'"


def parse(src):
    if type(input) is StringType:
        src = StringIO(src)
    tokens = Tokenizer(src)
    result = []
    while 1:
        token = tokens.get_token()
        if token == '' : break
        result.append(token)
    return result


def main():

    time0 = time.time()
    print '\
------------------------------------------------------------\n\
< Cluster trees by bipartition Program >\n\n\
This is a little script written by jinjianjun@mail.kib.ac.cn in python.\nIt has to do with a Criteria.tab file under "./" and Newick format *.tre files under "./Trees/".\n\
Each line in Criteria.tab assigns a kind of bipartition to cluster (Bipartition Name\ttaxon1,taxon2).\
It can... .\nThe author is a new pythoner. Please feel free to send him suggestions.\nBTW, He is working on phylogeny and biogeography of Cannabaceae.'
    bootstrap_threshold = int(raw_input("\nPlease assign the bootstrap threshold:"))
    statistics_for_trees = [['Tree']]
    #read criteria.tab
    criteria_tab = []
    while not len(criteria_tab):
        input_dir = raw_input('\nPlease assign the Criteria.tab:\n').strip()
        try:
            if input_dir == '':
                criteria_tab = read_tab('./Criteria.tab')
                this_directory = './'
            else:
                this_directory = input_dir.split('/')
                end_name = this_directory[-1]
                del this_directory[-1]
                this_directory = '/'.join(this_directory)+'/'
                criteria_tab = read_tab(this_directory+end_name)
        except:
            continue
        else:
            break
    criteria = []
    separate_taxa_in_criteria = []
    Prefix = []
    for i in range(0, len(criteria_tab)):
        if len(criteria_tab[i]) == 2:
            criteria.append([criteria_tab[i][0], set(criteria_tab[i][1].replace('"','').replace(', ',',').split(','))])
            separate_taxa_in_criteria.append(False)
        elif len(criteria_tab[i]) == 3:
            separate_taxa_in_criteria.append(True)
            criteria.append([criteria_tab[i][0], set(criteria_tab[i][1].replace('"','').replace(', ',',').split(',')), set(criteria_tab[i][2].replace('"','').replace(', ',',').split(','))])
        Prefix.append('')
        while os.path.exists(this_directory+'ClusterTrees_with_'+criteria[i][0]+'_Threshold_'+str(bootstrap_threshold)+Prefix[i]+'/'):
            Prefix[i] += '_'
        os.mkdir(this_directory+'ClusterTrees_with_'+criteria[i][0]+'_Threshold_'+str(bootstrap_threshold)+Prefix[i]+'/')
        statistics_for_trees[0].append(criteria[i][0])
    #read and transform trees into bipartition matrix
    all_files = os.listdir(this_directory+'Trees/')
    tree_files = []
    for i in range(0, len(all_files)):
        if all_files[i].endswith('.tre'):
            tree_files.append(all_files[i])
    print_warn_for_this_criteria = [False for x in criteria]
    for i in range(0, len(tree_files)):
        read_tree = open(this_directory+'Trees/'+tree_files[i], 'rU')
        tree_line = read_tree.read().split(';')[0].strip()
        trees_bipartition_infos = transform_to_bipartition(tree_line)
        taxa_names = trees_bipartition_infos[0]
        bipartition_matrix = trees_bipartition_infos[1]
        statistics_for_trees.append([tree_files[i]])
        for j in range(0, len(criteria)):
            if not separate_taxa_in_criteria[j]:
                #check
                if not print_warn_for_this_criteria[j] and not criteria[j][1].issubset(taxa_names):
                    print_warn_for_this_criteria[j] = True
                    print 'Error in Criteria tab file: Line '+str(j+1)+': '+criteria_tab[j][1]
                #
                if (criteria[j][1] in bipartition_matrix[0]):
                    if bipartition_matrix[1][bipartition_matrix[0].index(criteria[j][1])]>=bootstrap_threshold:
                        os.system('cp '+this_directory+'Trees/'+tree_files[i]+' '+this_directory+'ClusterTrees_with_'+criteria[j][0]+'_Threshold_'+str(bootstrap_threshold)+Prefix[j]+'/')
                    statistics_for_trees[i+1].append(bipartition_matrix[1][bipartition_matrix[0].index(criteria[j][1])])
                else:
                    statistics_for_trees[i+1].append('-')
            else:
                #check
                if not print_warn_for_this_criteria[j] and not set(list(criteria[j][1])+list(criteria[j][2])).issubset(taxa_names):
                    print_warn_for_this_criteria[j] = True
                    print 'Error in Criteria tab file: Line '+str(j+1)+': '+criteria_tab[j][1]+'\t'+criteria_tab[j][2]
                #
                if (criteria[j][1] in bipartition_matrix[0]) and (criteria[j][2] in bipartition_matrix[0]) and (set(list(criteria[j][1])+list(criteria[j][2])) not in bipartition_matrix[0]):
                    conflit_node_bs =[]
                    for k in range(0, len(bipartition_matrix[0])):
                        if (len(bipartition_matrix[0][k]) < len(taxa_names)-1) and (len(bipartition_matrix[0][k]) > 1) and abs(criteria[j][1].issubset(bipartition_matrix[0][k]) - criteria[j][2].issubset(bipartition_matrix[0][k])):
                            conflit_node_bs.append(bipartition_matrix[1][k])
                    if bipartition_matrix[1][bipartition_matrix[0].index(criteria[j][1])]>=bootstrap_threshold and bipartition_matrix[1][bipartition_matrix[0].index(criteria[j][2])]>=bootstrap_threshold and max(conflit_node_bs)>=bootstrap_threshold:
                        os.system('cp '+this_directory+'Trees/'+tree_files[i]+' '+this_directory+'ClusterTrees_with_'+criteria[j][0]+'_Threshold_'+str(bootstrap_threshold)+Prefix[j]+'/')
                    statistics_for_trees[i+1].append(min(max(conflit_node_bs), bipartition_matrix[1][bipartition_matrix[0].index(criteria[j][1])], bipartition_matrix[1][bipartition_matrix[0].index(criteria[j][2])]))
                else:
                    statistics_for_trees[i+1].append('-')
    statistics_for_trees = delete_suffix(delete_prefix(statistics_for_trees, True, 0), True, 0)
    write_tab(statistics_for_trees, this_directory+'Statistics_for_trees')
    
###########################################
if __name__=='__main__':
    main()