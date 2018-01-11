###########################################################################
###########################################################################
## Copyright (C) 2018  Guignard Leo <guingardl__@__janelia.hhmi.org>     ##
##                                                                       ##
## This program is free software: you can redistribute it and/or modify  ##
## it under the terms of the GNU General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or     ##
## (at your option) any later version.                                   ##
##                                                                       ##
## This program is distributed in the hope that it will be useful,       ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
## GNU General Public License for more details.                          ##
##                                                                       ##
## You should have received a copy of the GNU General Public License     ##
## along with this program.  If not, see <http://www.gnu.org/licenses/>. ##
###########################################################################
###########################################################################

import os
import pandas as pd
from loading_data import loading_data
import numpy as np
from matplotlib import cm
from random import randint
from matplotlib import pyplot as plt

def kinetic_stage(s, g, kinetic_groups, basic_value, stages):
    v = kinetic_groups.get(g, basic_value)
    return stages[np.where(stages==s)[0] + v][0]

def short_name(n):
    return n.split('.')[0]+'.'+str(int(n.split('.')[1][:-1]))

def get_daughter_names(n):
    letter=n.split('.')[0][0]
    z_cycle=str(int(n.split('.')[0][1:])+1)
    num1='%04d'%(int(n.split('.')[1][:-1])*2-1)
    num2='%04d'%(int(n.split('.')[1][:-1])*2)
    end=n[-1]
    return letter+z_cycle+'.'+num1+end, letter+z_cycle+'.'+num2+end

def order_cells(c1, c2, names):
    n=names[c1]
    f=int(n.split('.')[1][:-1])
    n=names[c2]
    s=int(n.split('.')[1][:-1])
    if f<s:
        return c1, c2
    else:
        return c2, c1

def get_sister_name(n):
    return n.split('.')[0]+'.'+'%04d'%(int(n.split('.')[1][:-1])+1)+n[-1]

def print_summary(f_path, test_model, values, results, cells_to_find_induction, cells_in_grey):
    from datetime import date
    f = open(f_path + 'Model_output_' + date.today().isoformat() + '.txt', 'w')
    print f_path + 'Model_output_' + date.today().isoformat() + '.txt'
    f.write('''
Model output ''' + date.today().isoformat() + '''

Parameter sweep:
    see parameter_sweep.pdf for the distribution of th possible combinations of FP/FN
    see Ground_truth.csv for the ground truth cell/ligand couples.
    Cells for which to find induction:\n''' +
    str(cells_to_find_induction) + '''
    Cells in the grey zone:\n''' +
    str(cells_in_grey) + '''

The best combination of parameters is the following:

(rho, tau, integration time, gamma, kfLt, delta1, delta6, delta8, delta10): (#Missed cell-ligand couple, #over-predictions, #Miss predictions)\n'''
            )
    for i, v in enumerate(values):
        tmp = v[:2] + v[3:]
        f.write('(%.1f, %d, %d, %.2E, %d, %.2E, %.2E, %.2E, %.2E): '%tmp)
        f.write('(%d, %d, %d)\n'%results[i])
    f.write('\n\n')
    f.write(test_model.format_output())
    f.write('''

For the recapitulative figure please find summary-figure.pdf
Where:
green means known cells where no induction happen
red means known cells where induction happen
blue means grey zone
''')

def load_model_data(path, stages, kinetic_groups, basic_value, small_half_life, z_c_end=9, time_end=180):

    (pos_AP, lin_tree, fates,
     fates2, fates3, vol,
     inv_lin_tree, surf_ex,
     surfaces, bary, properties,
     ColorMap) = loading_data(path)
    names = properties['Names'][0]

    tissue_order = ['Endoderm Head', 'Endodermal Strand 1',
                    'Endodermal Strand 2', 'Epidermis Head',
                    'Epidermis Tail', 'Germline',
                    'Mesoderm Mesenchyme', 'Mesoderm Muscle 1',
                    'Mesoderm Muscle 2', 'Mesoderm Notochord 1',
                    'Mesoderm Notochord 2',
                    'Mesoderm Trunk Lateral Cell',
                    'Mesoderm Trunk Ventral Cell',
                    'Neural Plate Head Dorsal',
                    'Neural Plate Head Ventral',
                    'Neural Plate Tail Dorsal',
                    'Neural Plate Tail Lateral',
                    'Neural Plate Tail Ventral']
    CMapFates = ColorMap(tissue_order+['undeter'], 'rainbow')

    f = open(path+'/Gene_expression/Ground_truth.csv')
    data = pd.read_csv(f, sep=',')
    f.close()

    ground_truth = {}
    ground_truth_pw = {}
    for row in data.values:
        if type(row[0])!=float:
            mother_name=row[0].lower()
            if not ground_truth.has_key(mother_name):
                ground_truth[mother_name]=[]
            if not ground_truth_pw.has_key(mother_name):
                ground_truth_pw[mother_name]=set()
            if type(row[4])==str:
                ground_truth[mother_name].append(tuple(row[4].split(';')))
            if type(row[5])==str:
                ground_truth[mother_name].append(tuple(row[5].split(';')))
            if type(row[7])==str:
                ground_truth[mother_name].append(tuple(row[7].split(';')))
            if type(row[8])==str:
                ground_truth[mother_name].append(tuple(row[8].split(';')))
            if type(row[3])==str:
                ground_truth_pw[mother_name].add(row[3])
            if type(row[6])==str:
                ground_truth_pw[mother_name].add(row[6])

    f=open(path+'/Gene_expression/2016-04-21-stage6-8-10-11-12-expression-pathways.csv', 'rU')
    import csv
    csv_file=csv.reader(f, delimiter=',')

    rows=[]
    for row in csv_file:
        rows.append(row)
    f.close()

    roles=np.array(['Rec', 'Rec', 'Mod', 'Lig'])

    possibilities=set()
    gene_expression={}
    gene_role={}
    gene_effect={}
    gene_affect_neighb={}
    gene_affect_self={}
    pathways={}
    gene_presence={}
    remaining=set()
    for r in rows[1:]:
        if (not 'Cadherin' in r[3] and 
            not 'TTK' in r[1] and 
            not 'whole embryo' in r[13] and 
            not 'Maternal' in r[13] and 
            not 'No expression' in r[13] and
            not 'BMP3' in r[1] and
            not 'Noggin' in r[1] and
            not 'ci-bmp2' in r[1] and
            not '12' in r[12]):
            # if 'FZD4' in r:
            #     break
            g_n = r[3].split(',')[0]+r[4]
            if not 'Stage' in r[12]:
                s_n = 'Stage ' + r[12]
            else:
                s_n = r[12]
            stage=kinetic_stage(s_n, g_n, kinetic_groups, basic_value, stages)
            stage_next = stages[np.where(stages==stage)[0]+1][0]
            if 'ci-tolloid' in r:
                name=tuple(['tolloid', 'tolloid', stage])
                name_next=tuple(['tolloid', 'tolloid', stage_next])
            else:
                name=tuple([g_n, g_n, stage])
                name_next=tuple([g_n, g_n, stage_next])
            tags=r[13].split(',')
            for p in r[3].split(','):
                pathways.setdefault(p.strip(), set()).add(name)
            gene_effect[name]=r[4]=='+'
            if not g_n in small_half_life:
                gene_effect[name_next]=r[4]=='+'
                gene_role[name_next]=roles[np.array(r[5:9])=='1'][0]
                gene_affect_neighb[name_next]=gene_affect_neighb.get(name_next, False) or r[10]=='1'
                gene_affect_self[name_next]=gene_affect_neighb.get(name_next, False) or r[9]=='1'
                for p in r[3].split(','):
                    pathways.setdefault(p.strip(), set()).add(name_next)
            gene_role[name]=roles[np.array(r[5:9])=='1'][0]
            gene_affect_neighb[name] = gene_affect_neighb.get(name, False) or r[10]=='1'
            gene_affect_self[name]= gene_affect_self.get(name, False) or r[9]=='1'
            gene_presence[name]=s_n
            if not (gene_role[name]=='Rec'):

                if not gene_expression.has_key(name):    
                    gene_expression[name]=set()
                if not gene_expression.has_key(name_next):    
                    gene_expression[name_next] = set()
                for terr in tags:
                    if 'cell pair' in terr.lower():
                        if len(terr)>1 and terr[0]==' ':
                            possibilities.add(terr[1:].lower().split(' cell pair')[0])
                            gene_expression[name].add(terr[1:].lower().split(' cell pair')[0])
                            if not g_n in small_half_life:
                                gene_expression[name_next].add(terr[1:].lower().split(' cell pair')[0])
                        else:
                            possibilities.add(terr.lower().split(' cell pair')[0])
                            gene_expression[name].add(terr.lower().split(' cell pair')[0])
                            if not g_n in small_half_life:
                                gene_expression[name_next].add(terr.lower().split(' cell pair')[0])
                    elif 'line)' in terr:
                        possibilities.add(terr.split('(')[-1].split(')')[0].split(' line')[0].lower())
                        gene_expression[name].add(terr.split('(')[-1].split(')')[0].split(' line')[0].lower())
                        if not g_n in small_half_life:
                            gene_expression[name_next].add(terr.split('(')[-1].split(')')[0].split(' line')[0].lower())
                    elif 'lines)' in terr:
                        for t in terr.split('(')[-1].split(')')[0].split(' lines')[0].replace(' ', '').split('&'):
                            possibilities.add(t.lower())
                            gene_expression[name].add(t.lower())
                            if not g_n in small_half_life:
                                gene_expression[name_next].add(t.lower())
                    elif 'line' in terr and terr.split('line')[-1]=='':
                        possibilities.add(terr.replace(' ', '')[0])
                        gene_expression[name].add(terr.replace(' ', '')[0])
                        if not g_n in small_half_life:
                                gene_expression[name_next].add(terr.replace(' ', '')[0])
                    elif 'whole embryo' in terr or 'Maternal' in terr:# or 'No expression' in terr:
                        terr = set('A', 'B', 'a', 'b')
                        gene_expression[name].add(terr)
                        if not g_n in small_half_life:
                                gene_expression[name_next].union_update(terr)
                        possibilities.union_update(terr)
                    else:
                        # print r[:5], s_n, terr
                        remaining.add(terr)

    def get_mother_name(n):
        letter=n[0]
        z_c=str(int(n.split('.')[0][1:])-1)
        num='%04d'%np.ceil(float(n.split('.')[1][:-1])/2)
        end=n[-1]
        return letter+z_c+'.'+num+end

    def short_name(n):
        return n.split('.')[0]+'.'+str(int(n.split('.')[1][:-1]))


    cells_by_stage = {
        'Stage 5':{'a':set(['a5.3', 'a5.4']),
                   'b':set(['b5.3', 'b5.4']),
                   'A':set(['a5.1', 'a5.2']),
                   'B':set(['b5.1', 'b5.2'])},
        
        'Stage 6':{'a':set(['a6.5', 'a6.6', 'a6.7', 'a6.8']),
                   'b':set(['b6.5', 'b6.6', 'b6.7', 'b6.8']),
                   'A':set(['a6.1', 'a6.2', 'a6.3', 'a6.4']),
                   'B':set(['b6.1', 'b6.2', 'b6.3', 'b6.4'])},
        
        'Stage 8':{'a':set(['a7.9', 'a7.10', 'a7.11', 'a7.12', 'a7.13', 'a7.14', 'a7.15', 'a7.16']),
                   'b':set(['b7.9', 'b7.10', 'b7.11', 'b7.12', 'b7.13', 'b7.14', 'b7.15', 'b7.16']),
                   'A':set(['a7.1', 'a7.2', 'a7.3', 'a7.4', 'a7.5', 'a7.6', 'a7.7', 'a7.8']),
                   'B':set(['b7.1', 'b7.2', 'b7.3', 'b7.4', 'b7.5', 'b7.6', 'b7.7', 'b7.8'])},
        
        'Stage 10':{'a':set(['a8.17', 'a8.18', 'a8.19', 'a8.20', 'a8.21', 'a8.22', 'a8.23', 'a8.24',
                         'a8.25', 'a8.26', 'a8.27', 'a8.28', 'a8.29', 'a8.30', 'a8.31', 'a8.32']),
                    'b':set(['b8.17', 'b8.18', 'b8.19', 'b8.20', 'b8.21', 'b8.22', 'b8.23', 'b8.24',
                         'b8.25', 'b8.26', 'b8.27', 'b8.28', 'b8.29', 'b8.30', 'b8.31', 'b8.32']),
                    'A':set(['a7.1', 'a7.2', 'a8.5', 'a8.6', 'a8.7', 'a8.8', 'a7.5', 'a8.11', 'a8.12',
                         'a8.13', 'a8.14', 'a8.15', 'a8.16']),
                    'B':set(['b7.1', 'b7.2', 'b8.5', 'b8.6', 'b8.7', 'b8.8', 'b7.5', 'b7.6', 'b7.7',
                         'b8.15', 'b8.16'])},

        'Stage 11':{'a':set(['a8.17', 'a8.18', 'a8.19', 'a8.20', 'a8.21', 'a8.22', 'a8.23', 'a8.24',
                         'a8.25', 'a8.26', 'a8.27', 'a8.28', 'a8.29', 'a8.30', 'a8.31', 'a8.32']),
                    'b':set(['b8.17', 'b8.18', 'b8.19', 'b8.20', 'b8.21', 'b8.22', 'b8.23', 'b8.24',
                         'b8.25', 'b8.26', 'b8.27', 'b8.28', 'b8.29', 'b8.30', 'b8.31', 'b8.32']),
                    'A':set(['a7.1', 'a7.2', 'a8.5', 'a8.6', 'a8.7', 'a8.8', 'a7.5', 'a8.11', 'a8.12',
                         'a8.13', 'a8.14', 'a8.15', 'a8.16']),
                    'B':set(['b7.1', 'b7.2', 'b9.9', 'b9.10', 'b8.6', 'b9.13', 'b9.14', 
                         'b9.15', 'b9.16', 'b7.5', 'b7.6', 'b7.7', 'b8.15', 'b8.16'])},
                    
        'Stage 12':{'a':set(['a8.17', 'a8.18', 'a8.19', 'a8.20', 'a8.21', 'a8.22', 'a8.23', 'a8.24',
                         'a8.25', 'a8.26', 'a8.27', 'a8.28', 'a8.29', 'a8.30', 'a8.31', 'a8.32']),
                    'b':set(['b8.17', 'b8.18', 'b8.19', 'b8.20', 'b8.21', 'b8.22', 'b8.23', 'b8.24',
                         'b8.25', 'b8.26', 'b8.27', 'b8.28', 'b8.29', 'b8.30', 'b8.31', 'b8.32']),
                    'A':set(['a7.1', 'a7.2', 'a8.5', 'a8.6', 'a8.7', 'a8.8', 'a7.5', 'a8.11', 'a8.12',
                         'a8.13', 'a8.14', 'a8.15', 'a8.16']),
                    'B':set(['b7.1', 'b7.2', 'b9.9', 'b9.10', 'b8.6', 'b9.13', 'b9.14', 
                         'b9.15', 'b9.16', 'b7.5', 'b7.6', 'b7.7', 'b7.15', 'b8.16'])}
    }

    new_gene_expression = {}

    for g, t in gene_expression.iteritems():
        new_t = set()
        for ti in t:
            if len(ti)!=1:
                new_t.add(ti)
            else:
                new_t.update([short_name(get_mother_name(n+'*')) for n in cells_by_stage[g[-1]][ti]])
        new_gene_expression[g] = new_t

    gene_expression = new_gene_expression


    import cPickle as pkl
    f=open(path+'sim_LR.pkl')
    sim_tree=pkl.load(f)
    f.close()
    f=open(path+'sim_nv.pkl')
    sim_nv=pkl.load(f)
    f.close()

    sim_tree={(2*10**4+k[0], 2*10**4+k[1]):v for k, v in sim_tree}

    cells=[k for k in lin_tree.keys() if k/10**4==1]

    for c in cells:
        vol[c]=vol[lin_tree[c][0]]

    inv_lin_tree={ v : k for k, values in lin_tree.iteritems() for v in values }

    i=1
    for c in cells:
        if int(names[c].split('.')[1][:-1])%2==1:
            s_n=get_sister_name(names[c])
            for ci in cells:
                if names[ci]==s_n:
                    if inv_lin_tree.get(c, [])==[]:
                        lin_tree[i]=[c, ci]
                    if sim_tree.has_key((lin_tree[c][0], lin_tree[ci][0])):
                        sim_nv[i]=sim_tree[(lin_tree[c][0], lin_tree[ci][0])]
                    else:
                        sim_nv[i]=sim_tree[(lin_tree[ci][0], lin_tree[c][0])]
                    names[i]=get_mother_name(names[c])
                    vol[i]=vol[c]+vol[ci]
                    i+=1

    fates3_names={names[c][:-1]:fates3.get(c, 'undetermined') for c in lin_tree.keys()}

    interest_cells=[c for c in lin_tree.keys() if (c/10**4<time_end and len(lin_tree[c])==2 and
                                                   int(names[lin_tree[c][0]].split('.')[0][1:])<=z_c_end)]
    couples=[]
    for c in interest_cells:
        if names[c][-1]=='_':
            for ci in interest_cells:
                if names[c][:-1]==names[ci][:-1] and c!=ci:
                    couples.append((c, ci))

    known_interactions = 0
    for v in ground_truth.values():
        known_interactions+=len(v)

    return (pos_AP, lin_tree, fates,
            fates2, fates3, vol, sim_nv,
            inv_lin_tree, surf_ex, names,
            surfaces, bary, properties,
            ColorMap, gene_expression,
            ground_truth, ground_truth_pw,
            fates3_names, couples, known_interactions,
            gene_affect_neighb, gene_affect_self, 
            pathways, gene_effect)

def plot_model_results(model, init_pos, nb, ax, col, space = 2, values = None, error_bars = None, hatch = ''):
    if values is None:
        a, b, c, d, e, f = model.retrieve_the_cat()
    else:
        a, b, c, d, e, f = values
    p1 = .2 + init_pos
    p2 = p1 + space*nb + 1
    ax.bar(p1, a/float(a+b)*100, width = 1.6, color = col, alpha = .7, lw = 2, hatch = hatch)
    ax.bar(p1, b/float(a+b)*100, width = 1.6, bottom = a/float(a+b)*100, color = 'w', alpha = .7, lw = 2)
    ax.bar(p2, e/float(e+f)*100, width = 1.6, color = col, alpha = .7, lw = 2, hatch = hatch)
    ax.bar(p2, f/float(e+f)*100, width = 1.6, bottom = e/float(e+f)*100, color = 'w', alpha = .7, lw = 2)
    if not error_bars is None:
        ax.errorbar([p1, p2], [a/float(a+b)*100, e/float(e+f)*100],
                    yerr=[error_bars[0], error_bars[4]], fmt='none', ecolor='k', linewidth=3, mew=3, alpha = .7)
        
        
def plot_model_results_hor(model, init_pos, nb, ax, col, space = 2, values = None, error_bars = None, hatch = ''):
    if values is None:
        a, b, c, d, e, f = model.retrieve_the_cat()
    else:
        a, b, c, d, e, f = values
    p1 = init_pos + .2
    p2 = init_pos + .2
    ax.barh(p1, a/float(a+b)*100-1, left = 1, height = 1.6, color = cm.Paired(.4), lw = 2)
    ax.barh(p2, -e/float(e+f)*100+1, left = -1, height = 1.6, color = cm.Paired(0), lw = 2)
    if not error_bars is None:
        ax.errorbar([a/float(a+b)*100, -e/float(e+f)*100], [p1, p2],
                    xerr=[error_bars[0], error_bars[4]], fmt='none', ecolor='k', linewidth=3, mew=3, alpha = .7)


def plot_box_like_distrib(X, mapping, ax, color, marker, label, val_color=50, max_val=np.inf, div=1., reverse = False):
    Y_med, Y_25, Y_75 = [], [], []
    X_ = []
    for k in X:
        if k<max_val:
            X_.append(k)
            if reverse:
                Y_med.append(100-100*np.median(mapping[k])/div)
                Y_25.append(100-100*np.percentile(mapping[k], 25)/div)
                Y_75.append(100-100*np.percentile(mapping[k], 75)/div)
            else:
                Y_med.append(100*np.median(mapping[k])/div)
                Y_25.append(100*np.percentile(mapping[k], 25)/div)
                Y_75.append(100*np.percentile(mapping[k], 75)/div)

    ax.plot(X_, Y_med, '-o', color=cm.Paired(color), marker=marker, label=label, lw = 2, ms = 8)
    ax.fill_between(X_, Y_25, Y_75, color=cm.Paired(color), alpha=.2)
    return Y_med, Y_25, Y_75

def produce_basic_model_results(output_path, model, values, results, cells_to_find_induction, cells_in_grey, title):
    folder = output_path + 'E_distributions'
    
    if not os.path.exists(folder):
        os.makedirs(folder)

    print_summary(output_path, model, values, results, cells_to_find_induction, cells_in_grey)

    model.print_cell_graph(folder = folder)

    fig = model.full_output_figure()
    fig.axes[0].set_title(title, fontsize = 35)
    fig.tight_layout()
    plt.savefig(output_path + 'recap-figure.pdf')
    
def get_random_id(z_c):
    return str(randint(1, 2**(z_c-1)/4))

def randomize_name(n):
    z_c = int(n.split('.')[0][1:])
    l = ['a', 'b'][randint(0,1)]
    id_ = get_random_id(z_c)
    return l+str(z_c)+'.'+id_
