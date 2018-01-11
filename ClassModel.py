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

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint as integration
from natsort import natsorted


class InteractionModel(object):
    """docstring for InteractionModel"""
    from random import seed, randint, random

    def build_training_matrix(self):
        ''' Args:
            Returns:
                M: {string: {string: bool} } {cell name: {pathway: True if differential induction} }
            This function builds a dictionary necessary to build the *full_output_figure*
        '''
        M={}
        reversed_pw = {vi[1]:k for k, v in self.pathways.iteritems() for vi in v}
        for c1, c2 in self.couples:
            n = self.short_name(self.names[c1])
            M[n] = {}
            poss_stages = list(self.stages[self.get_poss_stage(self.names[self.lin_tree[c1][0]], c1)])
            if self.ground_truth.has_key(n):
                for l, s in set([(k[0], k[2]) for k in self.ground_truth[n]]):
                    if poss_stages.index(s) == 0:
                        M[n][reversed_pw[l]+' pol'] = True
                    else:
                        M[n][reversed_pw[l]+' ind'] = True
            if n in self.cells_to_find_induction:
                M[n]['all'] = True
            elif n in self.no_induction_to_find:
                M[n]['all'] = False
        self.training_M = M
        undeter = []
        deter_l = []
        deter = []
        not_induced = []
        for c_n, v in self.training_M.iteritems():
            if v=={}:
                undeter.append(c_n)
            elif not v['all']:
                not_induced.append(c_n)
            elif len(v.keys())>1:
                deter_l.append(c_n)
            else:
                deter.append(c_n)
        self.names_order = natsorted(deter_l) + natsorted(deter) + natsorted(undeter) + natsorted(not_induced)
        return M

    def build_resulting_matrix(self):
        ''' Args:
            Returns:
                M: {string: {string: bool} } {cell name: {pathway: True if differential induction} }
            This function builds a dictionary necessary to build the *full_output_figure*
        '''
        M = {}
        reversed_pw = {vi[1]:k for k, v in self.pathways.iteritems() for vi in v}
        for c1, c2 in self.couples:
            n = self.short_name(self.names[c1])
            if self.inducer.has_key(n):
                poss_stages = list(self.stages[self.get_poss_stage(self.names[self.lin_tree[c1][0]], c1)])
                ind = self.inducer[n]
                M[n] = {'all':True}
                for pw, pw, s in ind:
                    if poss_stages.index(s) == 0:
                        M[n][pw[:-1] + ' pol'] = True
                    else:
                        M[n][pw[:-1] + ' ind'] = True
        self.resulting_M = M
        return M

    def print_resulting_M(self):
        ''' Args:
            Returns:
            Historical function.
            Now it just calls two other functions
            if necessary ...
        '''
        if not hasattr(self, 'resulting_M'):
            self.build_resulting_matrix()
        if not hasattr(self, 'names_order'):
            self.build_training_matrix()

    def full_output_figure(self):
        ''' Args:
            Returns:
                fig: matplotlib figure
            This function builds the "big" figure that sums up
            the whole output of the model. It might be difficult
            to read at first but it does the job
        '''
        self.print_resulting_M()
        pw_order = {'all': 0, 'ERK pol': 2, 'ERK ind': 3,
                              'Nodal pol': 5, 'Nodal ind': 6,
                              'Notch pol': 8, 'Notch ind': 9,
                              'BMP pol': 11, 'BMP ind': 12,
                              'Wnt pol': 14, 'Wnt ind': 15}


        fig = plt.figure(figsize = (25, 10))
        ax = fig.add_subplot(111)

        tmp = self.resulting_M
        for i, c in enumerate(self.names_order):
            for pw, on in tmp.get(c, {}).iteritems():
                if on:
                    symbol = 's'
                else:
                    symbol = 'x'
                ax.plot(i, 15-pw_order[pw], symbol+'w', mec = 'k', mew = 3, ms=10, alpha = .5)

        tmp = self.training_M

        nb_on_pw = 0
        nb_on = 0
        nb_off = 0
        for i, c in enumerate(self.names_order):
            done_pw = False
            done_on = False
            for pw, on in tmp[c].iteritems():
                if on:
                    symbol = 's'
                else:
                    symbol = 'o'
                if pw == 'all' and not on:
                    for pos in pw_order.itervalues():
                        ax.plot(i, 15-pos, symbol+'r', )
                else:
                    ax.plot(i, 15-pw_order[pw], symbol+'g')
                if pw != 'all' and on and not done_pw:
                    nb_on_pw += 1
                    done_pw = True
                elif on and not done_on:
                    nb_on += 1
                    done_on = True
                elif not on:
                    nb_off += 1

        ax.plot([-2, 200], [14, 14], 'k--', lw=3, alpha = .8)
        ax.plot([nb_on_pw - .5, nb_on_pw - .5], [-2, 20], 'k--', lw=3, alpha = .8)
        ax.plot([nb_on - .5, nb_on - .5], [-2, 20], 'k--', lw=3, alpha = .8)
        ax.plot([len(self.names_order) - nb_off - .5, len(self.names_order) - nb_off - .5], [-2, 20], 'k--', lw=3, alpha = .8)

        ax.set_yticks(15 - np.array(sorted(pw_order.values() + [2.5, 5.5, 8.5, 11.5, 14.5])))
        ax.set_yticklabels(('Any\nInduction',
                            'Polarisation', 'ERK', 'Induction',
                            'Polarisation', 'Nodal', 'Induction',
                            'Polarisation', 'Notch', 'Induction',
                            'Polarisation', 'BMP', 'Induction',
                            'Polarisation', 'Wnt', 'Induction',
                            ), fontsize = 18)
        ax.set_xticks(range(len(self.names_order)))
        ax.set_xticklabels(self.names_order, rotation = 'vertical', fontsize = 18)
        ax.set_xlim(-1, len(self.names_order) + 1)

        ax.set_ylim(-.5, 15.5)
        fig.tight_layout()
        return fig

    def get_A_B(self, n):
        ''' Args:
                n: string name of a cell
            Returns:
                tmp: 'A'/'B'/'a'/'b'
            This function returns one of these three letters
            according to the name n
        '''
        z_c=int(n.split('.')[0][1:])
        tmp=n[0]
        if int(n.split('.')[1])<=2**(z_c-1)/8:
            return tmp.upper()
        else:
            return tmp

    def short_name(self, n):
        ''' Args:
                n: string name of a cell
            Returns
                n: string name of a cell
            This function formats the name of cell
            from 'a7.0012*' to 'a7.12'
        '''
        return n.split('.')[0]+'.'+str(int(n.split('.')[1][:-1]))

    def get_avg_vol(self, c):
        ''' Args:
                c: int, cell id
            Returns:
                -: float, average volume
            This function compute the average volume of
            a cell c between t(c) + self.begin and 
            t(c) + self.begin + self.end
            where t(c) is the time of c (c//10**4)
            the volume is normalised by .3**3, the resolution of the image
            The norm is hard-coded, SAD!
        '''
        ii = 0
        c_tmp = c
        volumes = []
        while ii < self.begin and len(self.lin_tree.get(c_tmp, [])) == 1:
            c_tmp = self.lin_tree[c_tmp][0]
            ii += 1

        while ii <= self.end and len(self.lin_tree.get(c_tmp, [])) == 1:
            volumes += [self.vol[c_tmp]]
            c_tmp = self.lin_tree[c_tmp][0]
            ii += 1
        return np.mean(volumes)*.3**3


    def get_min_vol(self, k, param, size=10, decal=5):
        ''' Args:
                k: int, cell id
                param: {int: float/int} {cell id: metric} dictionary of metric on cells
                size: period of time to consider in time points
                decal: number of starting points to not consider
            Returns:
                -: float
            This function computes and returns the min value between the
            two sister cell of *k* of the average in time during *size*
            times points from *decal* time points of the metric held in *param*. 
        '''
        c1=self.lin_tree[k][0]
        c2=self.lin_tree[k][1]
        out=[]
        for i in range(decal):
            c1=self.lin_tree[c1][0]
            c2=self.lin_tree[c2][0]
        for i in range(size):
            out.append([param[c1], param[c2]])
            c1=self.lin_tree[c1][0]
            c2=self.lin_tree[c2][0]
        return min(np.mean(out, axis=0))

    def get_max_vol(self, k, param, size=10, decal=5):
        ''' Args:
                k: int, cell id
                param: {int: float/int} {cell id: metric} dictionary of metric on cells
                size: period of time to consider in time points
                decal: number of starting points to not consider
            Returns:
                -: float
            This function computes and returns the max value between the
            two sister cell of *k* of the average in time during *size*
            times points from *decal* time points of the metric held in *param*. 
        '''
        c1=self.lin_tree[k][0]
        c2=self.lin_tree[k][1]
        out=[]
        for i in range(decal):
            c1=self.lin_tree[c1][0]
            c2=self.lin_tree[c2][0]
        for i in range(size):
            out.append([param[c1], param[c2]])
            c1=self.lin_tree[c1][0]
            c2=self.lin_tree[c2][0]
        return max(np.mean(out, axis=0))

    def get_poss_stage(self, n, c):
        ''' Args:
                n: string, cell name
                c: int, cell id
            Returns:
                -: [int, int] list of indices of stages to take into account
            This function follow rules to return the stages that correspond
            to a given cell and its mother based on the zygotic cycle and
            on a parameter of the model
        '''
        z_c = (int(n.split('.')[0][1:]))
        if z_c == 6:
            return [0, 1]
        elif z_c == 7:
            return [1, 2]
        elif c/10**4 < self.time_112:
            return [2, 3]
        else:
            return [3, 4]

    def get_avg_contact_surf(self, c, begin=0, end=5):
        ''' Args:
                c: int, cell id
                begin: number of time points after c to start
                end: number of time points after begin to end
            Returns:
                -: {int: float}, {neighbor id: average area}
            This function compute the average area of contact
            between a cell and its neighbors, for each neighbor,
            between t(c) + begin to t(c) + begin + end
            where t(c) is the time of c (c//10**4)
        '''
        c_tmp = c
        # move forward to begin
        for i in range(begin):
            c_tmp=self.lin_tree[c_tmp][0]

        # build the temporary dictionary
        surf_ex_tmp = {n:[s] for n, s in self.surf_ex[c_tmp].iteritems() if not n%10**4 == 1}
        i = begin + 1
        c_tmp = self.lin_tree[c_tmp][0]
        # build a corresponding dictionary to take into account dividing cells between begin and end
        n_corres = {ni:n for n in surf_ex_tmp.iterkeys() for ni in self.lin_tree[n]} #if not n%10**4 == 1 

        while i<=end and len(self.lin_tree.get(c_tmp, [])) == 1:
            # going forward in the lineage tree, retrieving the surface of contacts
            # and adding them to their corresponding orignal cells
            surf_ex_t = {}
            n_corres_tmp = {}
            for n, s in self.surf_ex[c_tmp].iteritems():
                if n%10**4 != 1:
                    surf_ex_t.setdefault(n_corres.get(n, n), []).append(s) # adding the surface to the list of surfaces
                    for ni in self.lin_tree.get(n, []):
                        n_corres_tmp[ni] = n_corres.get(n, n) # updating the corresponding dictionary
            for n, s in surf_ex_t.iteritems():
                surf_ex_tmp.setdefault(n, []).append(np.sum(s)) 
                # if two cells have been added to the same corresponding cell, it is a division
                # in this case we add the area of the two sister cells (double check here?)
            n_corres = n_corres_tmp
            c_tmp = self.lin_tree[c_tmp][0]
            i += 1
        max_len = np.max([len(s) for s in surf_ex_tmp.itervalues()])
        # Note that that I had to rewrite this function several times,
        # it is likely that there is a more elegant (and maybe correct?) way to write that
        return {n:np.sum(s)/max_len for n, s in surf_ex_tmp.iteritems()}


    def get_previous(self, n):
        ''' Args:
                n: string cell name without '_'/'*'
            Returns:
                string: name of the mother of n
        '''
        n1, n2=n.split('.')
        return n1[0]+str(int(n1[1:])-1)+'.'+str(int(np.ceil(int(n2)/2.)))

    def get_possible_tags(self, n, up=4):
        ''' Args:
                n: string cell name without '_'/'*'
                up: earliest zygotic cycle to create
            Returns:
               c_possibilities: set of possible tags/territories for the cell n 
        '''
        c_possibilities=set([n])
        tmp=n
        while int(tmp.split('.')[0][1:])>up:
            # adding to the set the name of the mother cells
            tmp = self.get_previous(tmp)
            c_possibilities.add(tmp)

        z_c = int(n.split('.')[0][1:]) # zygotic cycle
        tmp = n[0] # a/b
        if int(n.split('.')[1]) <= 2**(z_c-1)/8: # animal/vegetal side (a/A, b/B)
            c_possibilities.add(tmp.upper())
        else:
            c_possibilities.add(tmp)

        c_possibilities.add('Maternal')
        c_possibilities.add('whole embryo')
        return c_possibilities

    def get_surf_bilateral(self, c, sis, surf, m_d = 'M'):
        ''' Args:
                c: int, cell id to consider
                sis: [int, int] list of the two sister cells
                     note that one of the ids is c
                surf: {c(int): {n(int): area(float)} }; {cell id: {neighb id: area of contact} }
                m_d: string 'D' if we are looking at daughter induction
                            'M' if we are looking at mother polarization
            Returns:
                neighbs: {string: float}; {cell name: area of contact} (the name does not include the L/R differentiation)
            This function compute the surface of contact to the neighbors of a given cell c
            grouping L/R cells together and adding noise/scaling if necessary
        '''
        neighbs = {}
        # Creates an entry for the multiplicative coefficient necessary
        # to preserve symmetry and reflexivity when adding noise
        if not hasattr(self, 'multi_coef'):
            self.multi_coef = {} 
        if not self.multi_coef.has_key(self.short_name(self.names[c])):
            self.multi_coef[self.short_name(self.names[c])] = {}

        for ci, si in surf.iteritems():
            # ci: neighbor, si: area to that neighbor
            if ci%10**4!=1 and self.names.has_key(ci) and (not ci in sis or m_d == 'D'):
                # We ensure that we are not considering the surface to the exterior, id%1e4 == 1
                # and that the cell is an actual cell and that we are not considering contact 
                # between sister cells when considering polarisation of the mother
                if self.rand and not self.multi_coef[self.short_name(self.names[c])].has_key(self.short_name(self.names[ci])):
                    # in the case of noise addition
                    multi = self.random()*(2*self.ratio_noise/100.) # computation of the multiplicator coefficient
                    # note that random() returns a value in [0, 1[ following a normal distribution
                    if self.randint(0, 1):
                        multi = 1 - multi
                    else:
                        multi = 1 + multi
                    # increasing or decreasing noise
                    self.multi_coef[self.short_name(self.names[c])][self.short_name(self.names[ci])] = multi
                    if not self.multi_coef.has_key(self.short_name(self.names[ci])):
                        self.multi_coef[self.short_name(self.names[ci])] = {}
                    self.multi_coef[self.short_name(self.names[ci])][self.short_name(self.names[c])] = multi
                    # the coefficient is saved for c -> ci and for ci -> c to ensure reflexivity
                    # the name used is the name without '_'/'*' to ensure symmetry
                    
                neighbs.setdefault(self.short_name(self.names[ci]), 0)
                if not self.take_surf:
                    neighbs[self.short_name(self.names[ci])] += (si!=0) * 2200.
                    # in the case where we don't consider the surfaces (semi-binary run), it is replaced, when existing by 817.
                    # 817 hard-coded ... SAD!
                elif self.rand:
                    neighbs[self.short_name(self.names[ci])] += si * self.multi_coef[self.short_name(self.names[c])][self.short_name(self.names[ci])]
                    # when we want to apply a randomisation of the surface of contact
                else:    
                    if self.size_increase!=1:
                        neighbs[self.short_name(self.names[ci])] += si * self.size_increase
                        # when we want to apply a scaling
                    else:
                        neighbs[self.short_name(self.names[ci])] += si
        return neighbs

    def get_terr(self, c, n):
        ''' Args:
                c: int, cell id
                n: {int: float}; {neighb: area}
            Returns:
                terr_inter: {string: float}; {territory: area}
            This function creates a mapping of the territories to
            their corresponding area allowing to assess the 
            surfaces of contact back in time. For example, if c is in contact
            with A7.1 we can extrapolate the surface it had with A6.1
            by combining this surface (A7.1) to the surface with A7.2 if it exists
        '''
        terr_inter = {}
        surf_tot = np.sum(n.values())
        terr_inter={self.get_previous(self.short_name(self.names[c])):surf_tot, self.short_name(self.names[c]):surf_tot}
        for ni, s in n.iteritems():
            p_t=self.get_possible_tags(ni)
            for p_ti in p_t:
                if terr_inter.has_key(p_ti):
                    terr_inter[p_ti]+=s
                else:
                    terr_inter[p_ti]=s
        terr_inter['whole embryo']=surf_tot
        terr_inter['Maternal']=surf_tot
        return terr_inter

    def express(self, tag, gene, self_tag):
        return (tag in self.gene_expression.get(gene, []) and
                (self.gene_affect_self[gene] or not tag in self_tag))

    def get_ligand_expression(self, c, tags, self_tag, l=None):
        ''' Args:
                c: int, cell id to consider
                tags: {string: float}; {territory: area}
                self_tag: tags corresponding to the cell 
                          itself for autocrine ligands
                l: (string, )*3; (l name, l name, stage) 
                   ligand to consider, None if all ligands
            Returns:
                -: {(string, )*3: float}; {l: area of contact}
        '''
        gene_inter_s = {}
        # restriction of the ligands to look for
        if l is None:
            tmp = self.gene_expression
        else:
            tmp = {l: self.gene_expression[l]}

        tot_surf = tags['Maternal']
        for gene, t in tmp.iteritems():
            # case by case dealing of agonist/antagonist interactions
            # note that gene is the name of the ligand/antagonist in the form
            # (ligand-name, ligand-name, stage), ie: ('ERK+', 'ERK+', 'Stage 8')
            # and that t is the set of territories that express that
            # ligand/antagonist at that stage
            # Note the redundancy of ligand-name, it's historical
            if not 'ERK' in gene[0] and not 'tolloid' in gene[0] and '+' in gene[0]:
                # in the case of everything but ERK or tolloid (they are specials)
                ligand = gene # ligand names
                inhib = (gene[0][:-1] + '-', gene[0][:-1] + '-', gene[-1]) #antagonist name
                # Note that for historical reason the ligands are called genes ...
                gene_inter_s[gene] = 0.
                if self.gene_affect_self[ligand] and set(self_tag).intersection(self.gene_expression[ligand]) != set():
                    # in the case where the cell itself expresses the ligand and the ligand is autocrine
                    # the whole surface is considered from which will be removed the surfaces that express antagonists
                    gene_inter_s[gene] = tot_surf*(.3**2)
                    # Note that the .3**2 is the resolution of the image, hard coded here ... SAD!
                    if not ('BMP+' == ligand[0] and set(self_tag).intersection(self.gene_expression.get(('tolloid', 'tolloid', ligand[-1]), [])) != set()):
                        # when we are not in the case of BMP "protected" by tolloid (it avoids BMP antagonists to act)
                        for tag, s in tags.iteritems():
                            if (self.express(tag, inhib, self_tag) and not
                                ('BMP+' == ligand[0] and self.express(tag, ('tolloid', 'tolloid', ligand[-1]), self_tag))):
                                # when the antagonist is expressed in a given tag/territory the corresponding surface is removed
                                gene_inter_s[gene]-=s*(.3**2)
                else:
                    # when the cell does not express to itself
                    for tag, s in tags.iteritems():
                        if (self.express(tag, ligand, self_tag) and (not self.express(tag, inhib, self_tag) or
                                                                ('BMP+' == ligand[0] and self.express(tag, ('tolloid', 'tolloid', ligand[-1]), self_tag)))):
                            # if the ligand is expressed with no expression of the antagonist or protection of tolloid for BMP
                            # the area of contact is added to the area of contact to the given ligand
                            gene_inter_s[gene]+=s*(.3**2)


            elif 'ERK' in gene[0]:
                # in the case of ERK (FGF/Eph) it is much simpler
                for tag, s in tags.iteritems():
                    if tag in t:
                        if 'ERK-' in gene[0]:
                            if not gene_inter_s.has_key(gene):
                                gene_inter_s[gene]=0.
                            gene_inter_s[gene]+=s*(.3**2)
                        elif 'ERK+' in gene[0]:
                            if not gene_inter_s.has_key(gene):
                                gene_inter_s[gene]=0.
                            if tag in self_tag:
                                gene_inter_s[gene]+=tot_surf*(.3**2)
                            else:
                                gene_inter_s[gene]+=s*(.3**2)
                    # note that Eph (ERK-) is not autocrine but FGF (ERK+) is
        return {k: max(0, min(v, tot_surf*(.3**2))) for k, v in gene_inter_s.iteritems()}

    def Rt(self, c, sis, l, surf = None, m_d = 'M'):
        ''' Args:
                c: int id of the sister cell to consider
                sis: [int, int] list of the two sister cells
                     note that one of the ids is c
                l: (string, string, string) tuple representing the ligand:
                   (PW, PW, stage), note that for historical reasons the
                   first two stings are now equals
                surf: {c(int): {n(int): area(float)} }; {cell id: {neighb id: area of contact} }
                m_d: string 'D' if we are looking at daughter induction
                            'M' if we are looking at mother polarization
            Returns:
                -: float AiL, the surface of contact of c with ligand l
        '''
        # a precomputed areas of contact to neighbors can be given if not the ponctual areas are used
        if surf is None: 
            surf = self.surf_ex.get(c, {})

        # computes the surface of contact grouping together bilateral cells
        # and adding noise/scaling if necessary
        neighbs = self.get_surf_bilateral(c, sis, surf, m_d) 
        
        # save the new surfaces to allow the ratio before/after noise/scaling
        if not hasattr(self, 'avg_contact_surf_af_noise'):
            self.avg_contact_surf_af_noise = {}
        self.avg_contact_surf_af_noise[c] = neighbs

        # Transform the cell names in territories.
        # It is mainly historical when gene expression was in weird forms
        # A line/b line ...
        # Now it just adds "whole embryo" and "maternal" as territories
        # It also allows to "go back in time" by adding the surfaces of sister cell
        # allocating it ot the mother cell
        terr_inter = self.get_terr(c, neighbs)

        # building of the territories that represent the cell itself
        self_tag=[self.get_previous(self.short_name(self.names[c])), # mother cell
                  self.get_previous(self.get_previous(self.short_name(self.names[c]))),  # grand mother cell
                  self.short_name(self.names[c]), # cell itself
                  self.get_A_B(self.short_name(self.names[c]))] # A/B/a/b

        # Build the surfaces of contact of the cell c to a given ligand l (AiL)
        gene_inter_s = self.get_ligand_expression(c, terr_inter, self_tag, l)
        
        if not hasattr(self, 'surf_to_ligand'):
            self.surf_to_ligand = {}
        if not self.surf_to_ligand.has_key(c):
            self.surf_to_ligand[c] = {}
        self.surf_to_ligand[c][l] = gene_inter_s.get(l, 0)

        return gene_inter_s.get(l, 0)

    def analytic_Rc(self, t, Aij, kfLt):
        ''' Args:
                t: float, time (in seconds)
                Aij: float area of contact ($\mu m^2$)
                kfLt: $k_f.L^T$ value
            Returns:
                -: Rc for the given t, Aij, kfLt
            This function computes the Rc for the
             analytic solution found (see paper)
        '''
        return Aij * (1 - np.exp(-kfLt * t))

    def F(self, alphaR_ke, kfLt, t):
        return (alphaR_ke/kfLt) * (kfLt * t - ( 1 - np.exp(-kfLt * t) ))

    def analytic_E(self, A, alphaR_ke, kfLt, t):
        return (1 - np.exp( -A * self.F(alphaR_ke, kfLt, t) ))

    # def analytic_E(self, Rc, gamma):
    #     return Rc / (Rc + gamma)

    def get_E_C(self, c, sis, ligand, const = False):
        ''' Args:
                c: int id of the sister cell to consider
                sis: [int, int] list of the two sister cells
                     note that one of the ids is c
                ligand: (string, string, string) tuple representing the ligand:
                        (PW, PW, stage), note that for historical reasons the
                        first two stings are now equals
                const: Bool whether the computation is for Polarisation
                       or for sister induction
            Returns:
                E: [float, ] list of E* in time (can be ploted against t of this function for example)
                C: [float, ] list of Rc in time (can be ploted against t of this function for example)
                   note that for historical reason this variable is called C instead of Rc
                -: {c(int): float}; {cell id: surface ratio} ratio of surface in contact 
                   with ligand before and after noise/scaling, necessary for normalisation by the volume
            This function computes E* and Rc given a cell and a ligand knowing
             whether it is a polarisation or an induction
        '''
        sis = list(sis)
        kfLt = self.kfLt.get(ligand, 10**-3)
        if const:
            if not self.avg_contact_surf_const.has_key(c):
                # Computes an average of the surface of contacts
                self.avg_contact_surf_const[c] = self.get_avg_contact_surf(c, self.begin, self.begin + 5)
            A = self.Rt(c, sis, ligand, self.avg_contact_surf_const[c]) # computes A not considering sis-sis cell exchanges
            tmp_s = self.avg_contact_surf_const[c]
        else:
            if not self.avg_contact_surf.has_key(c):
                # Computes an average of the surface of contacts
                self.avg_contact_surf[c] = self.get_avg_contact_surf(c, self.begin, self.begin + self.integration_time/120)
            A = self.Rt(c, sis, ligand, self.avg_contact_surf[c], 'D') # computes A
            tmp_s = self.avg_contact_surf[c]

        # Note that for historical reasons AiL is wrongly named A

        def E_prime(E, t, kfLt, Rt): # definition dE*/(ET.dt)
            dE = kfLt * (1 - E) * self.analytic_Rc(t, Rt, kfLt)
            return dE

        # times for which to compute the E* and Rc
        t = list(np.arange(0, self.integration_time, self.dt)) + [self.integration_time]
        # we use integration function of scipy package (see paper)
        # E0 = 0
        # E = integration(E_prime, E0, t, args = (kfLt, A))
        # This is not the right computation of Rc, it assumes alpha_R = 1
        Rc = self.analytic_Rc(np.array(t), A, kfLt)
        E  = self.analytic_E(A, self.alphaR_ke, kfLt, np.array(t))

        # computation of the ratio after/before noise/scaling
        tot_surf_after_noise = sum(self.avg_contact_surf_af_noise[c].values())
        tot_surf_before_noise = sum([v for ci, v in tmp_s.iteritems() if ci%1e4 != 1])
        if not hasattr(self, 'surf_ratio'):
            self.surf_ratio = {}
        self.surf_ratio[c] = tot_surf_after_noise / tot_surf_before_noise
        return E, Rc, self.surf_ratio[c]

    def taking_care_of_ERK(self, gi, v_n):
        ''' Args:
                gi: {di: {Lj: L-Eij*} }; {daughter name: {Ligand name: E*} }
            Returns:
                gi: {di: {Lj: L-Eij*} }; {daughter name: {Ligand name: E*} }
            This function performs the necessary division for the ERK pathway
        '''
        if self.ERK_delta:
            for c_d, pw_info in gi.items():
                z_c = int(c_d.split('.')[0][1:])     
                ERK_v = {}
                Eph_v = {}
                to_pop = []
                for pw, v in pw_info.items():
                    if 'ERK+' in pw:
                        to_pop += [(pw, v)]
                    elif 'ERK-' in pw:
                        to_pop += [(pw, v)]
                for pw, v in to_pop:
                    if 'ERK+' in pw:
                        pw_info.pop(pw)
                        ERK_v[pw] = v
                    elif 'ERK-' in pw:
                        pw_info.pop(pw)
                        Eph_v[pw] = v                    
                for pw, v in ERK_v.iteritems():
                    if '6' in pw[2]:
                        delta = self.delta6
                    elif '8' in  pw[2]:
                        delta = self.delta8
                    else:
                        delta = self.delta10
                    # gi[c_d][pw] = v / (Eph_v.get(('ERK-', 'ERK-', pw[2]), 0) + delta)
                    gi[c_d][pw] = (v / (self.delta1 * Eph_v.get(('ERK-', 'ERK-', pw[2]), 0) + v + delta))# / v_n[c_d][pw]
        return gi

    def get_expression_sisters(self, sis, stage_to_care):
        ''' Args:
                sis: [string, string] list of two sister cells
                stage_to_care: [string, string] list of two stages
            Returns:
                tmp: {di: {Lj: L-Eij*} }; {daughter name: {Ligand name: E*} }
            Builds a dictionary for two sister cells accounting for the E* value
            for the ligands at the specified stages
        '''
        gene_inter_s = {}
        v_n = {}
        for c in sis:
            v_n[c] = {}
            gene_inter_s[c] = {}
            v_to_norm = self.get_avg_vol(c) # Volume with which the E* will be normed
            for ligand in self.gene_expression.iterkeys():
                if (ligand[2] in stage_to_care):
                    E, Rc, R = self.get_E_C(c, sis, ligand, const=ligand[2] == stage_to_care[0]) # get the E* values
                    # note that the const=ligand[2] == stage_to_care[0] refers to the fact that the function has to
                    # know if we are considering mother polarisation or sister induction
                    if self.take_surf:
                        v_to_norm = (v_to_norm * R**(3./2))
                    E = E[-1]
                    v_n[c][ligand] = v_to_norm
                    if not 'ERK' in ligand[0]:
                        gene_inter_s[c][ligand] = E #/ v_to_norm
                    # elif not 'ERK' in ligand[0]:
                    #     gene_inter_s[c][ligand] = E / v_to_norm
                    else: # if it is ERK (FGF/Eph) we do not norm by the volume yet ()
                        gene_inter_s[c][ligand] = E

        # Formating of the output for historical reasons
        gi = {self.names[c]: {k: v * self.signal_increase for k, v in n.iteritems()} for c, n in gene_inter_s.iteritems()}
        v_n = {self.names[c]: {k: v * self.signal_increase for k, v in n.iteritems()} for c, n in v_n.iteritems()}
        tmp = self.taking_care_of_ERK(gi, v_n) # Special treatment done for ERK pathway (cf paper)
        return tmp

    def common(self, giL, giR):
        ''' Args:
                giL: {di: {Lj: L-Eij*} }; {daughter name: {Ligand name: E*} } for the left cell
                giR: {di: {Lj: R-Eij*} }; {daughter name: {Ligand name: E*} } for the corresponding right cell
            Returns:
                common_genes: {di: {Lj: Eij*} } where Eij* = min(L-Eij*, R-Eij*)
            Builds a consensus on the E* between the left and right cell.
        '''
        common_genes={}
        for k, v in giL.iteritems():
            if k[-1]=='*':
                other='_'
            else:
                other='*'
            genes2_L=set(v.keys())
            genes2_R=set(giR[k[:-1]+other].keys())
            inter2=genes2_L.intersection(genes2_R)

            if self.take_right:
                to_be_R = giR[k[:-1]+other]
            else:
                to_be_R = v
            if self.take_left:
                to_be_L = v
            else:
                to_be_L = giR[k[:-1]+other]

            common_genes[k[:-1]]={g:np.min([to_be_L[g], to_be_R[g]]) for g in inter2}

        return common_genes

    def get_distributions(self):
        ''' Args:
                ---
            Returns:
                ---
            Builds the E* distributions and from that extracts the threshold values for every pathway
        '''
        tmp={}
        # See comments in self.run()
        for c1, c2 in self.couples:
            sis1 = tuple(self.lin_tree[c1])
            sis2 = tuple(self.lin_tree[c2])
            if self.names[c1][1]=='7' and c1/10**4<14:
                while sis1[0]/10**4<self.time_112 or sis2[0]/10**4<self.time_112:
                    sis1 = tuple(self.lin_tree[sis1[0]]) + tuple(self.lin_tree[sis1[1]])
                    sis2 = tuple(self.lin_tree[sis2[0]]) + tuple(self.lin_tree[sis2[1]])

            n = self.names[self.lin_tree[c1][0]]
            if not self.gene_repartition.has_key(self.short_name(self.names[c1])):
                stage_to_care = self.stages[self.get_poss_stage(n, c1)]
                giL = self.get_expression_sisters(sis1, stage_to_care)
                giR = self.get_expression_sisters(sis2, stage_to_care)
                giC = self.common(giL, giR)
                self.gene_repartition[self.short_name(self.names[c1])] = giC
                if not hasattr(self, 'giL'):
                    self.giL = {}
                if not hasattr(self, 'giR'):
                    self.giR = {}
                self.giL[self.short_name(self.names[c1])] = giL
                self.giR[self.short_name(self.names[c1])] = giR
            else:
                giC = self.gene_repartition[self.short_name(self.names[c1])]

            # Gathering of the non-zero values
            for g, v in giC.values()[0].iteritems():
                if v>0 : 
                    tmp.setdefault(g[0], []).append(v)
            for g, v in giC.values()[1].iteritems():
                if v>0 : 
                    tmp.setdefault(g[0], []).append(v)

        # Construction of the low thresholds
        self.gene_cat_low = {}
        for k, v in tmp.iteritems():
            self.gene_cat_low[k] = np.percentile(v, self.low_th)


    def print_E(self, f, print_all = True):
        ''' Args:
                f: file descriptor/file name
                print_all: bool, whether we want the thresholded values or not
            Returns:
                ---
            Builds a csv file containing all the E* star values in a "readable" fashion
        '''
        import pandas as pd
        pws = ['Nodal+','Notch+','BMP+', 'Wnt+', 'ERK+', 'ERK-']

        columns = list(np.array(zip([p + ' M' for p in pws], [p + ' D' for p in pws])).flatten())
        columns = ['Mother', 'Daughter'] + columns
        list_of_lists = []
        mothers = self.gene_repartition.keys()
        mothers.sort()
        for mother in mothers:
            sisters =  self.gene_repartition[mother]
            c = min([ci for ci in self.lin_tree.keys() if self.names[ci] == sisters.keys()[0]+'*'])
            stage_to_care = self.stages[self.get_poss_stage(sisters.keys()[0], c)]
            it = True
            sister = sisters.keys()
            sister.sort()
            for sis in sister:
                pw_values = sisters[sis]
                if it:
                    list_tmp = [mother, sis]
                else:
                    list_tmp = ['', sis]
                for pw in pws:
                    for s in stage_to_care:
                        v = pw_values.get((pw, pw, s), 0.)
                        if (v > self.gene_cat_low[pw]) or print_all:
                            list_tmp.append(v)
                        else:
                            list_tmp.append('0.')

                list_of_lists.append(list_tmp)
                it ^= True
            ratio_list = ['', 'Ratio']
            for i, (v1, v2) in enumerate(zip(list_of_lists[-1][2:], list_of_lists[-2][2:])):
                v1 = float(v1)
                v2 = float(v2)
                if v1!=0 and v2!=0:
                    ratio_list.append(max(v1, v2)/min(v1, v2))
                elif v1!=0 or v2!=0:
                    ratio_list.append(max(v1, v2))
                else:
                    ratio_list.append('-')
            list_of_lists.append(ratio_list)
            list_of_lists.append([''] * (len(pws)*2+2))
            csv_out = pd.DataFrame(list_of_lists, columns = columns)
            csv_out.to_csv(f)

    def format_output(self):
        ''' Args:
                ---
            Returns:
                s: string, a formating of the output of the model
        '''
        s = 'Number of over-predictions: ' + str(self.false_pos) + '\n'
        s += 'Number of missed inductions: ' + str(self.false_neg) + '\n'
        s += 'Number of cell/ligand couples missed: ' + str(self.known_interactions - self.found_pw) + '\n'
        missed_ligand = [v for v in self.missed if v[-1]=='missed ligand']
        missed_inductions = [v for v in self.missed if v[-1]=='induction not found']
        over_predictions = [v for v in self.missed if v[-1]=='induction found']
        s += 'Cell/ligand couples missed: ' + '\n'
        for ml in missed_ligand:
            s += '\tCell: ' + ml[0] + ', ligand: ' + ml[3][0] + ', (' + str(ml[3][1]) + ')'  + '\n'

        s += 'Missed induction: ' + '\n'
        for mi in missed_inductions:
            s += '\tCell: ' + mi[0] + '\n'
        s += 'Over predictions:' + '\n'
        for mi in over_predictions:
            s += '\tCell: ' + mi[0] + '\n'

        return s

    def print_cell_graph(self, fs=7, folder = 'E_distributions'):
        ''' Args:
                fs: int, fontsize of the labels
                folder: string, folder where the figures will be saved
            Returns:
                ---
            Builds and saves the E* distribution figures
            for each pathway at each stage
        '''

        # The combination dictionary already contains all the necessary informations
        for s, values in self.combination.items():
            fig = plt.figure(figsize = (10, 8))
            ax = fig.add_subplot(111)

            cell_names = values.keys()
            cell_names_sorted = natsorted(cell_names)
            x_pos = 1

            for name in cell_names_sorted:
                v1, v2  = values[name]
                if (self.ratio_th < max(v1, v2)/(min(v1, v2)+10**-10)) and (self.gene_cat_low[s[0]] < max(v1, v2)):
                    lcol = 'g'
                else:
                    lcol = 'r'

                # sets the colors that have to be put in for the two sister cells (in hexa form)
                if v1 < self.gene_cat_low[s[0]]:
                    v1_col = 'r'
                elif (v2 < v1) and self.pw_s_induced[s][name] and self.max_E_for_induction[s] <= v1:
                    v1_col = '#00FF38'
                elif (v1 < v2) and self.pw_s_induced[s][name] and v1 <= self.min_E_for_uninduction[s]:
                    v1_col = '#FF69B4'
                elif (v2 < v1) and self.pw_s_induced[s][name] and v1 <= self.min_E_for_uninduction[s]:
                    v1_col = '#654321'
                elif not self.pw_s_induced[s][name] and self.max_E_for_induction[s] <= v1:
                    v1_col = '#F9FB1F'
                elif not self.pw_s_induced[s][name] and v1 < self.min_E_for_uninduction[s]:
                    v1_col = '#FB991F'
                else:
                    v1_col = '#E3E3E3'

                if v2 < self.gene_cat_low[s[0]]:
                    v2_col = 'r'
                elif (v1 < v2) and self.pw_s_induced[s][name] and self.max_E_for_induction[s] <= v2:
                    v2_col = '#00FF38'
                elif (v2 < v1) and self.pw_s_induced[s][name] and v2 <= self.min_E_for_uninduction[s]:
                    v2_col = '#FF69B4'
                elif (v1 < v2) and self.pw_s_induced[s][name] and v2 <= self.min_E_for_uninduction[s]:
                    v2_col = '#654321'
                elif not self.pw_s_induced[s][name] and self.max_E_for_induction[s] <= v2:
                    v2_col = '#F9FB1F'
                elif not self.pw_s_induced[s][name] and v2 < self.min_E_for_uninduction[s]:
                    v2_col = '#FB991F'
                else:
                    v2_col = '#E3E3E3'



                ax.plot([x_pos, x_pos+1], [v1, v2], '-', c = lcol, lw=3, ms = 10)
                ax.plot([x_pos], [v1], 'o', c = v1_col, lw=3, ms = 10)
                ax.plot([x_pos+1], [v2], 'o', c = v2_col, lw=3, ms = 10)
                x_pos += 2

            # Finalise the figure
            ax.plot([0, x_pos], [self.gene_cat_low[s[0]]]*2, 'b-', alpha = .6, lw = 2)
            ax.plot([0, x_pos], [self.max_E_for_induction[s]]*2, 'g--', alpha = .6, lw = 2)
            ax.plot([0, x_pos], [self.min_E_for_uninduction[s]]*2, 'r--', alpha = .6, lw = 2)
            ax.set_xlim(-1, x_pos+1)
            if 'ERK' in s[0]:
                ax.set_yscale('log')
            ax.grid(axis='x')
            ax.set_title(s)
            ax.set_xticks(np.arange(1, x_pos, 2) + .5)

            ax.set_xticklabels(cell_names_sorted, rotation = 'vertical', ha = 'center', fontsize=fs)
            plt.savefig(folder + '/cell_E_distribution_' + '_'.join(s) + '.pdf')
            plt.close()

    def print_out(self, c_bg='k', title='', th_showed=True, decal = .2, output_folder = ''):
        ''' Args:
                c_bg: string/float/(int, )*3/(int, )*4, background color
                title: string, title of the figure
                decal: float, increment in Y of position between cells
                output_folder: string, path to the output folder
            Returns:
                number: int, number of cells for which an induction has to be found
                not_number: int, number of cells for which an induction has to be NOT found
                fig: produced figure
        '''

        number=0
        not_number=0
        fig=plt.figure(figsize=(10, 8))
        ax=fig.add_subplot(111)
        ind_pos = 0

        # Organise the Y position of every couple of cells.
        # They are first sorted by lin tree distance
        # Then ther Y position is a linear function of this order
        tmp = []
        for c1, c2 in self.couples:
            mean_nv=(self.sim_nv[c1]+self.sim_nv[c2])/2
            tmp += [mean_nv]
        pos = {}
        ind_pos = 0
        for i in np.argsort(tmp):
            ind_pos = (ind_pos + decal) % 1
            pos[self.couples[i][0]] = ind_pos

        # Header of the csv file
        f=open(output_folder + 'inductions_in_text.csv', 'w')
        for c, s in self.ground_truth.iteritems():
            f.write(c)
            for si in s:
                f.write(','+si[0]+','+si[-1]+'\n')
        f.write('cell name,lineage tree distance average,min volume ratio,polarisation,induction,inducers,polarisation/induction to find?\n')

        # Each couple of L/R cell is treated independantly
        for c1, c2 in self.couples:
            mean_nv = (self.sim_nv[c1]+self.sim_nv[c2])/2
            sis1 = self.lin_tree[c1]
            sis2 = self.lin_tree[c2]
            ratio1 = self.get_max_vol(c1, self.vol)/self.get_min_vol(c1, self.vol)
            ratio2 = self.get_max_vol(c2, self.vol)/self.get_min_vol(c2, self.vol)
            mean_ratio = np.mean([ratio1, ratio2])

            cell_name = self.short_name(self.names[c1])
            f.write(cell_name+',%f,%f,'%(mean_nv, mean_ratio))
            induced = self.induced[cell_name]
            if induced:
                Y = 2+pos[c1]
            else:
                Y = 0+pos[c1]
                
            stage_to_care = self.stages[self.get_poss_stage(self.names[self.lin_tree[c1][0]], c1)]
            inducer_stages = [s[2] for s in self.inducer.get(cell_name, [])]

            # defines the markers for the figure function of where the induction comes from
            if stage_to_care[0] in inducer_stages:
                if stage_to_care[1] in inducer_stages:
                    marker=(6,1,0) # 6 branches star
                    ms=20
                    f.write('yes,yes,')
                else:
                    marker='^' # upward triangle 
                    ms=10
                    f.write('yes,no,')
            elif stage_to_care[1] in inducer_stages:
                marker='v'  # downward triangle 
                ms=10
                f.write('no,yes,')
            else:
                marker='o'
                ms=10
                f.write('no,no,')
            for ind in self.inducer.get(cell_name, []):
                f.write(str(ind).replace(',', ';') + ' ')
            f.write(',')

            if cell_name in self.cells_to_find_induction: # The cell was meant to be induced
                ax.plot(mean_nv, Y, 'm', marker=marker, ms=ms, mec=c_bg, mew=1, alpha=1.)
                not_number+=1
                f.write('Yes,')
            elif (self.f_care.get(self.lin_tree[c1][0], '1')==self.f_care.get(self.lin_tree[c1][1], '2') and 
                  mean_nv<=self.sim_ratio_th and mean_ratio<=1.25 and not cell_name in self.grey_zone): # The cell was meant NOT to be induced
                ax.plot(mean_nv, Y, 'b', marker=marker, ms=ms, mec=c_bg, mew=1, alpha=1.)
                number+=1
                f.write('No,')
            else:
                ax.plot(mean_nv, Y, 'y', marker=marker, ms=ms, mec=c_bg, mew=1, alpha=.7)
                f.write('?')
            f.write('\n')

        # set up the figure
        ax.set_ylim((-1, 4))
        ax.set_xlim((-0.01, .5))
        ax.set_yticks([])

        ax.plot([-10, .51], [1.5, 1.5], c_bg+'--', lw=3, alpha=.5)
        if th_showed:
            ax.set_title((title + '\n' +
                          'Cell polarization Mother-Daughter; low th:'+str(self.low_th)+' ratio th:'+str(self.ratio_th)) + '\n'+
                           'kf.Lt: %.2E, int time: %d, delta: %.2E' % (self.kfLt.values()[0], self.integration_time, self.delta6))
        else:
            ax.set_title((title + '\n' +
                          'Cell polarization Mother-Daughter'))
        ax.set_xlabel('Daughters distance score', fontsize=30)
        ax.set_ylabel('Cell induction', fontsize=20)
        f.close()
        return number, not_number, fig

    def get_sister_names(self, n):
        ''' Args:
                n: string, cell name (no '_' or '*')
            Returns:
                tuple: (string, )*2 name of the sister cells
        '''
        l = n[0]
        z_c = str(int(n.split('.')[0][1:])+1)
        nb = int(n.split('.')[1])
        s1 = '%04d'%(2 * nb - 1)
        s2 = '%04d'%(2 * nb)
        return l + z_c + '.' + s1, l + z_c + '.' + s2

    def induction_responsible(self, cells):
        ''' Args: 
                cells: list of cell ids
            Returns: 
                decisive_cells: list of cell ids
            Given a list of cells returns the cells that are
            responsible of there differential induction if any.
        '''
        if not hasattr(self, 'inducer'):
            self.run()

        decisive_cells = set()
        for c in cells:
            name = self.short_name(self.names[c])
            stage_to_care = self.stages[self.get_poss_stage(name, c)[0]]
            for pw in self.inducer.get(name, []):
                if stage_to_care in pw:
                    decisive_cells.add(c)
                    for neighb in self.surf_ex[c].iterkeys():
                        if neighb%10**4 != 1 and self.short_name(self.names[neighb]) in self.gene_expression[pw]:
                            decisive_cells.add(neighb)
        return decisive_cells

    def retrieve_the_cat(self):
        ''' Function used to build the summary figure of the model
        '''
        cat1_ind = []
        cat1_no_ind = []
        cat2_ind = []
        cat2_no_ind = []
        for c1, c2 in self.couples:
            mean_nv=np.mean([self.sim_nv[c1], self.sim_nv[c2]])
            ratio1=self.get_max_vol(c1, self.vol)/self.get_min_vol(c1, self.vol)
            ratio2=self.get_max_vol(c2, self.vol)/self.get_min_vol(c2, self.vol)
            mean_ratio=np.mean([ratio1, ratio2])
            mother_name = self.short_name(self.names[c1])
            if not (mother_name in self.grey_zone or mean_nv > .1 or mean_ratio > 1.25):
                if self.induced[mother_name]:
                    cat1_ind += [mother_name]
                else:
                    cat1_no_ind += [mother_name]
                # if self.induced[mother_name]:
                #     cat2_ind += [mother_name]
                # else:
                #     cat2_no_ind += [mother_name]
            if mother_name in self.cells_to_find_induction:
                if self.induced.get(mother_name, False):
                    cat2_ind += [mother_name]
                else:
                    cat2_no_ind += [mother_name]

        found = 0
        not_found = 0
        for c, to_find in self.ground_truth.iteritems():
            for inducer in to_find:
                if inducer in self.inducer.get(c, []):
                    found += 1
                else:
                    not_found += 1
        return len(cat1_ind), len(cat1_no_ind), len(cat2_ind), len(cat2_no_ind), found, not_found

    def get_daughters(self, n):
        l = n[0]
        cycle = int(n.split('.')[0][1:])
        id_ = int(n.split('.')[1])
        # RL = n.split('.')[1][-1]
        return '%s%d.%04d'%(l, cycle+1, id_*2-1), '%s%d.%04d'%(l, cycle+1, id_*2)

    def run(self):
        ## Initialise a seed for run involving randomization
        self.seed(self.init_seed)
        self.found_pw = 0 # Number of cells with a differential induction among the cells in know_induction
        self.false_pos = 0 # Number of cells with a differential induction that shouldn't have any
        self.false_neg = 0 # Number of cells for which the right pw have not been found as a driver of differenti
        self.missed = []     
        self.combination = {}
        self.no_induction_to_find = []
        self.tmp_stage_to_care = {}
        for c1, c2 in self.couples: # each couple of mother cell (L/R) is treated independently
            # Avg of volume ratio between c1 daughters and c2's
            ratio1 = self.get_max_vol(c1, self.vol)/self.get_min_vol(c1, self.vol)
            ratio2 = self.get_max_vol(c2, self.vol)/self.get_min_vol(c2, self.vol)
            mean_ratio=np.mean([ratio1, ratio2])
            # Avg of lin tree distance between c1 daughters and c2's
            mean_nv = (self.sim_nv[c1]+self.sim_nv[c2])/2.
            sis1 = [tmp[1] for tmp in sorted([(self.names[c], c) for c in self.lin_tree[c1]])]#self.lin_tree[c1]
            sis2 = [tmp[1] for tmp in sorted([(self.names[c], c) for c in self.lin_tree[c2]])]#self.lin_tree[c2]
            if self.names[c1][1] == '7' and c1/10**4 < 14: # Take care of the issue generated by an ill definition of the 112 stage in ascidian embryos
                while sis1[0]/10**4<self.time_112 or sis2[0]/10**4<self.time_112:
                    sis1 = self.lin_tree[sis1[0]] + self.lin_tree[sis1[1]]
                    sis2 = self.lin_tree[sis2[0]] + self.lin_tree[sis2[1]]

            n = self.names[self.lin_tree[c1][0]]

            # get the relevant developmental stages
            stage_to_care = self.stages[self.get_poss_stage(n, c1)] # get the relevant developmental stages
            self.tmp_stage_to_care[self.get_daughters(self.short_name(self.names[c1]))[0]] = stage_to_care
            self.tmp_stage_to_care[self.get_daughters(self.short_name(self.names[c1]))[1]] = stage_to_care
            if not self.gene_repartition.has_key(self.short_name(self.names[c1])): # avoid to compute twice the same value
                giL = self.get_expression_sisters(sis1, stage_to_care) # dictionary that give for each sister cell its E* values for every ligand
                giR = self.get_expression_sisters(sis2, stage_to_care) # same but for the right cell
                giC = self.common(giL, giR) # build a consensus between L and R
                # save the computed values
                self.gene_repartition[self.short_name(self.names[c1])] = giC 
                if not hasattr(self, 'giL'):
                    self.giL = {}
                if not hasattr(self, 'giR'):
                    self.giR = {}
                self.giL[self.short_name(self.names[c1])] = giL
                self.giR[self.short_name(self.names[c1])] = giR
            else:
                giC = self.gene_repartition[self.short_name(self.names[c1])]
            
            mother_name = self.short_name(self.names[c1])
            
            pws = ['ERK+', 'BMP+', 'Nodal+', 'Notch+', 'Wnt+'] # Pathway to account for

            # Build the following dictionary:
            # {tuple (representing a pw, ie ('ERK+', 'ERK+', 'Stage 8')): {cell_name: [E*1, E*2 (for sisters 1 and 2)] } }
            # This dictionary allows to perform the decision on wheter a cell is induced or not
            # according to the rule described in the paper 
            for s in stage_to_care:
                for p in pws:
                    p_tmp = (p, p, s)
                    if not self.combination.has_key(p_tmp):
                        self.combination[p_tmp] = {}

                    sis1, sis2 = self.get_sister_names(mother_name)
                    self.combination[p_tmp][mother_name] = [giC[sis1].get(p_tmp, 0)]
                    self.combination[p_tmp][mother_name].append(giC[sis2].get(p_tmp, 0))

            # Fill the list of cells for which we are "sure" there are no induction
            if (self.f_care.get(self.lin_tree[c1][0], '1')==self.f_care.get(self.lin_tree[c1][1], '2')
                and mean_nv<=self.sim_ratio_th and mean_ratio<=1.25 and not mother_name in self.grey_zone):
                self.no_induction_to_find.append(mother_name)

        # Checking if the Emax and Emin have to be computed, necessary when you want to impose an external one
        max_E_for_induction_to_fill = self.max_E_for_induction == {}
        min_E_for_uninduction_to_fill = self.min_E_for_uninduction == {}
        
        # Initialisation of the Emax and Emin dictionary
        # Note that max_E_for_induction should actually be min_E_for_induction (and also for min_E...)
        # Note that I am a bit lazy to refactorize right now
        if max_E_for_induction_to_fill:
            self.max_E_for_induction = {k:np.inf for k in self.combination.iterkeys()}

        if min_E_for_uninduction_to_fill:
            self.min_E_for_uninduction = {k:0 for k in self.combination.iterkeys()}


        self.pw_s_induced = {k:{} for k in self.combination.iterkeys()}
        self.induced = {}
        self.inducer = {}

        # Computes the E* value under which there is no induction as described in the paper
        for s, values in self.combination.iteritems():
            cell_names = values.keys()
            for name in cell_names:
                v1, v2  = values[name]
                # the rule is that the value is the max of the min between E*1 and E*2 (E*1 (2) E* of sister 1 (2))
                # such that the ratio of E* ratio_th < E*1/E*2 
                # and low_th < E*1 or E*2 
                if (self.ratio_th < max(v1, v2)/(min(v1, v2)+10**-10)) and (self.gene_cat_low[s[0]] < max(v1, v2)):
                    self.pw_s_induced[s][name] = True
                    if min_E_for_uninduction_to_fill:
                        self.min_E_for_uninduction[s] = max(self.min_E_for_uninduction[s], min(v1, v2))
                else:
                    self.pw_s_induced[s][name] = False

        # Computes the E* value above which there is induction as described in the paper
        if max_E_for_induction_to_fill:
            for s, values in self.combination.iteritems():
                cell_names = values.keys()
                for name in cell_names:
                    v1, v2  = values[name]
                    # the rule is that the value is the min of the max between E*1 and E*2 (E*1 (2) E* of sister 1 (2))
                    # such that the ratio of E* ratio_th < E*1/E*2 
                    # and Emin < E*1 or E*2 
                    if self.pw_s_induced[s][name] and self.min_E_for_uninduction[s] < max(v1, v2):
                        if max(v1, v2) < self.max_E_for_induction[s]:
                            self.max_E_for_induction[s] = max(v1, v2)

        # knowing the Emin and Emax values computes what are the cells that are induced differentially
        self.induction_number = 0
        self.c_see = {}
        self.c_induced = {}
        for s, values in self.combination.items():
            cell_names = values.keys()
            for name in cell_names:
                v1, v2  = values[name]
                n1, n2 = self.get_daughters(name)
                # to be differentially induced a couple of sister cells has to satisfy the following rule:
                # E*1 (2) < Emin and Emax < E*2 (1)
                self.pw_s_induced[s][name] = ((min(v1, v2) <= self.min_E_for_uninduction[s])
                                              and (self.max_E_for_induction[s] <= max(v1, v2)))
                
                if v1 != 0:
                    self.c_see.setdefault(self.short_name(n1+'_'), []).append(s)
                if v2 != 0:
                    self.c_see.setdefault(self.short_name(n2+'_'), []).append(s)
                # Fill some variables that can be useful when looking at the results
                self.induced[name] = (self.induced.get(name, False) or self.pw_s_induced[s][name])
                if self.pw_s_induced[s][name]:
                    self.inducer.setdefault(name, []).append(s)
                
                if self.max_E_for_induction[s] <= v1:
                    if s[-1] == self.tmp_stage_to_care[n1][0]:
                        self.c_induced.setdefault(name, []).append([s, 0])
                    else:
                        self.c_induced.setdefault(self.short_name(n1+'_'), []).append([s, 1])

                    self.induction_number += 1
                if self.max_E_for_induction[s] <= v2:
                    if s[-1] == self.tmp_stage_to_care[n2][0] and not [s, 0] in self.c_induced.get(name, []):
                        self.c_induced.setdefault(name, []).append([s, 0])
                    else:
                        self.c_induced.setdefault(self.short_name(n2+'_'), []).append([s, 1])
                    self.induction_number += 1

        # Fill some variables about the performance of the model that can be useful when looking at the results
        self.cell_pairs_missed = {}
        for name, known_inductions in self.ground_truth.iteritems():
            for induction in known_inductions:
                if self.pw_s_induced[induction][name]:
                    self.found_pw += 1
                    self.cell_pairs_missed[(name, ) + induction] = 0
                else:
                    self.cell_pairs_missed[(name, ) + induction] = 1
                    self.missed.append((name, 0, 0, (induction[0], induction[2]), 'missed ligand'))

        for name in self.cells_to_find_induction:
            if not self.induced.get(name, True):
                self.false_neg += 1 
                self.missed += [(name, 'induction not found')]

        for name in self.no_induction_to_find:
            if self.induced[name]:
                self.false_pos += 1
                self.missed += [(name, 'induction found')]

    def __init__(self, lin_tree, surf_ex, couples, low_th, hi_th, ratio_th, 
                 kfLt, cells_to_find_induction, f_care,  names, fates3_names, 
                 sim_nv, vol, surfaces, gene_expression, gene_affect_neighb,
                 gene_affect_self, inv_lin_tree, pathways, stages,
                 gene_effect, ground_truth, ground_truth_pw, dt = 2, integration_time = 1800,
                 begin=2, end=6, take_kfLt=False, take_surf=True, known_interactions = 14,
                 rand=0, time_112=24, sim_ratio_th=.12, time_end=60, grey_zone = [], alphaR_ke = 1,
                 z_c_end=9, rand_kfLt=False, ratio_noise=.2, rand_gene_surf = False,
                 gene_cat_hi = None, gene_cat_low = None, gene_repartition = None, basic_euler = False,
                 delta = None, tmp = False, size_increase=1., take_left=True, take_right=True,
                 avg_contact_surf = None, avg_contact_surf_const = None, take_ERK = True,
                 min_E_for_uninduction = None, max_E_for_induction = None, ERK_delta = True,
                 signal_increase = 1.):
        ''' Args:
                see the comments inside this same function
            Returns:
            This function intialize the set of parameters of the model
        '''
        self.known_interactions = known_interactions # number of known interaction
        # self.avg_contact_surf_af_noise = {} # necessary
        # self.surf_to_ligand = {}
        self.dt = dt # dt to do the integration, not necessary anymore, historical
        # self.surf_ratio = {}
        if delta is None:
            delta = [1,]*4
        self.delta1, self.delta6, self.delta8, self.delta10 = delta # delta in SOS*/(SOS* + RasGDP + delta), see paper
        self.grey_zone = grey_zone # manual grey zone
        self.integration_time = integration_time # induction time for the integration
        self.ground_truth = ground_truth # manual ground truth
        self.stages = stages # possible developmental stages to consider ordered chronologicaly
        self.pathways = pathways # dict that maps a pathway to a ligand. outdated but still used by a few functions
        self.inv_lin_tree = inv_lin_tree # inverse of lineage tree {daughter: mother}
        self.gene_affect_self = gene_affect_self # dict that maps a ligand id to a boolean. True if the ligand is autocrine
        self.gene_affect_neighb = gene_affect_neighb # dict that maps a ligand id to a boolean. True if the ligand is paracrine (always true ...)
        self.gene_expression = gene_expression # dict that maps a ligand id to a set of territories (cell name, 'A', 'whole embryos', ...)
        self.sim_nv = sim_nv # dict that maps a cell id to the lineage tree distance between its two daughters
        self.vol = vol # dict that maps a cell id to its volume
        self.names = names # dict that maps a cell id to its name
        self.lin_tree = lin_tree # dict that maps a cell id to the list L of its decendents in the next time point (most of the time len(L) == 1)
        self.surf_ex = surf_ex # dict that maps a cell id c to a dict that maps a neighbor id to the surface it shares with c
        self.init_seed = self.randint(0, 10000) # seed initialisation for potential reproducibility
        self.ratio_noise = ratio_noise # float, noise to introduce
        self.couples = couples # list of pair of cell ids. The pairs are symmetrical cells
        self.size_increase = size_increase # scaling for the incrase of size in surface.
        self.kfLt = kfLt # dict that maps a ligand to its kfLt
        self.begin = begin # number of time points not to consider after division
        self.end = end # end - begin gives the number of time points to consider after begin
        self.take_surf = take_surf # boolean to decide to consider or not the surfaces
        self.rand = rand # boolean to decide to randomize or not the surfaces
        self.low_th = low_th # low threshold in %
        self.ratio_th = ratio_th # ratio threshold
        self.time_112 = time_112 # time point when starts the 112 cell stage
        self.cells_to_find_induction = cells_to_find_induction # manual input of cells for which we have to find inductions
        self.sim_ratio_th = sim_ratio_th # similiarity ratio under which a cell is considered not undergoing a cell fate restriction event
        self.f_care = f_care # dict that maps a cell id to its fate
        self.take_right = take_right # boolean to decide to consider or not the right part of the embryo
        self.take_left = take_left # boolean to decide to consider or not the left part of the embryo
        self.ERK_delta = ERK_delta # boolean to decide to consider ERK as special entity
        self.alphaR_ke = alphaR_ke # this parameter is the $alpha_R.k_e$
        self.signal_increase = signal_increase
        if gene_repartition is None:
            self.gene_repartition = {}
        else:
            self.gene_repartition = gene_repartition
        # dict that maps a cell name to a dict that maps daughter name to a dict that maps ligand to its E*
        # this dict can be given as an input to increase the parameter space search
        
        if avg_contact_surf is None:
            self.avg_contact_surf = {}
        else:
            self.avg_contact_surf = avg_contact_surf
        # dict that maps a cell id c to a dict that maps neighbor id to its average contact area to c
        # this dict can be given as an input to increase the parameter space search

        if avg_contact_surf_const is None:
            self.avg_contact_surf_const = {}
        else:
            self.avg_contact_surf_const = avg_contact_surf_const
        # dict that maps a cell id c to a dict that maps neighbor id to its average contact area to c (for a short period of time)
        # this dict can be given as an input to increase the parameter space search

        if gene_cat_low is None or gene_cat_hi is None:
            self.get_distributions()
        else:
            self.gene_cat_hi = gene_cat_hi
            self.gene_cat_low = gene_cat_low
        # dicts that map a ligand to its low and high (historical) threshold
        # these dicts can be given as an input to simulate mutation in the embryo

        if min_E_for_uninduction is None:
            self.min_E_for_uninduction = {}
        else:
            self.min_E_for_uninduction = min_E_for_uninduction
        # dicts that map a ligand to its Emin
        # these dicts can be given as an input to simulate mutation in the embryo

        if max_E_for_induction is None:
            self.max_E_for_induction = {}
        else:
            self.max_E_for_induction = max_E_for_induction
        # dicts that map a ligand to its Emax
        # these dicts can be given as an input to simulate mutation in the embryo
