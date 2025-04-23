#!/usr/bin/env python3

import argparse
import sys, re
from ete3 import Tree, TreeStyle, NodeStyle, Face, TextFace, add_face_to_node

"""
    augment_genetree takes a species/gene tree in ETE Tree() format and a list of reconciliation lines.
    
    It adds eventsRec information to each node such as:
        1. event type (speciation, duplication, transfer)
        2. leaf node
        3. transferBack tag (if the node recieves a transfer)
        4. species tree mapping
    
    It also inserts losses (where appropriate) into the tree. See insert_loss, loss_speciation,
    loss_duplication, and loss_transfer for more information.

    It returns the annotated gene tree in ETE Tree() format.
"""
def augment_genetree(species_tree,gene_tree,reclines):
    def sort_by_node(val):
        #Used as a key in later list.sort() calls
        return val[0].name

    def leaf(node,line):
        #Adds leaf event type to node
        stree_map = node.name
        if '_' in stree_map: stree_map=stree_map[:stree_map.index('_')]

        node.add_features(type='leaf',mapping=stree_map)        

    def speciation(node,line):
        #Adds speciation event type to node
        #Line is one line of the reconciliation output from Ranger (from reclines)
        stree_map = line[line.find('Mapping') + 12:len(line)]
        lca_leaves = line[line.find('[')+1:line.find(']:')].split(', ')
        #Takes data from reconciliation output line and adds to node

        node.add_features(type='speciation',mapping=stree_map,lca=lca_leaves)
 
    def duplication(node,line):
        #Adds duplication event type to node
        stree_map = line[line.find('Mapping') + 12:len(line)]
        lca_leaves = line[line.find('[')+1:line.find(']:')].split(', ')

        node.add_features(type='duplication',mapping=stree_map,lca=lca_leaves)

    def transfer(node,line):
        #Adds "branchingOut" (part 1 of a transfer event) event type to node. 
        def transferBack(node,recipient):
            #Takes a descendant of original transfer node, and adds "transferBack" event type.
            #This is the second part of a transfer event 
            global genetree

            if node.type:
            #Checks if the node already has its type specified
                node.add_features(transferBack=recipient)
            else :
                node.add_features(type='transferBack',mapping=stree_map,transferBack=recipient)
        
        global genetree
        global speciestree
        
        stree_map = line[(line.find('Mapping --> ') + 12):line.find(', R')]
        recipient = line[(line.find('Recipient -->')+14):len(line)]
        lca_leaves = line[line.find('[')+1:line.find(']:')].split(', ')
        
        node.add_features(type='branchingOut',mapping=stree_map,lca=lca_leaves)

        s_recipient = speciestree&"%s" %(recipient)
            #& is used to return first node named recipient (using formatted string) under speciestree
        g_children = node.get_children()
        ch1 = speciestree&"%s" %(g_children[0].mapping)
        ch2 = speciestree&"%s" %(g_children[1].mapping)

        """
            Based on species tree mappings, determines which edge reprents the transfer recipient.
            If g' (ch1) or g'' (ch2) directly map to the recipient, then that is the transfer edge.
            Else, whichever child maps to a descendant of the recipient will be the transfer edge.
        """
        if ch1 == s_recipient:
            transferBack(g_children[0],recipient)
        elif ch2 == s_recipient:
            transferBack(g_children[1],recipient)
        else:
            transfer_descendants = s_recipient.get_descendants()

            if ch1 in transfer_descendants:
                transferBack(g_children[0],recipient)
            elif ch2 in transfer_descendants:
                transferBack(g_children[1],recipient)
    
    def loss_speciation(node):
        global genetree
        global speciestree
        
        g = speciestree&"%s" %(node.mapping)
        children = node.get_children()
        ch1 = speciestree&"%s" %(children[0].mapping)
        ch2 = speciestree&"%s" %(children[1].mapping)
        #node and children are in gene tree, g and ch1 and ch2 are mappings

        distance1 = (g.get_distance(ch1)) - 1
        distance2 = (g.get_distance(ch2)) - 1
        path1 = []
        path2 = []
        losses1 = []
        losses2 = []

        if distance1 != 0:
            curr = ch1
            while curr != g:
                if curr.up == g: 
                    curr = curr.up
                else:
                    path1.append(curr.up)
                    curr = curr.up
            #path1 is now a list of nodes between (non-inclusive) g and ch1
        if path1:
            for p in path1:
                losses1 += [(p,c) for c in p.get_children() if c not in path1 and c != ch1 and c != ch2]
                #c is a child of a node in the path from g to ch1 that isn't part of that path
            if losses1: 
                losses1.sort(key=sort_by_node)
            prev = None
            for pair in losses1:
                p = pair[0]
                l = pair[1]             #why is this l when it was c earlier?
                if prev == None:
                    desc = children[0]
                    prev = insert_loss(node,desc,p.name,l.name)
                else: prev = insert_loss(node,prev,p.name,l.name)
        
        if distance2 != 0:
        #lines 112-119 but for path2
            curr = ch2
            while curr != g:
                if curr.up == g: 
                    curr = curr.up
                else:
                    path2.append(curr.up)
                    curr = curr.up
        if path2:
        #lines 120-132 but for path2
            for p in path2:
                losses2 += [(p,c) for c in p.get_children() if c not in path2 and c != ch1 and c != ch2]
            if losses2: 
                losses2.sort(key=sort_by_node)
            prev = None
            for pair in losses2:
                p = pair[0]
                l = pair[1]
                if prev == None:
                    desc = children[1]      #children[1] because this is for path2
                    prev = insert_loss(node,desc,p.name,l.name)
                else: prev = insert_loss(node,prev,p.name,l.name)

    def loss_duplication(node): 
        global genetree
        global speciestree

        g = speciestree&"%s" %(node.mapping)
        children = node.get_children()
        ch1 = speciestree&"%s" %(children[0].mapping)
        ch2 = speciestree&"%s" %(children[1].mapping)
        
        distance1 = 0
        distance2 = 0
        if g == ch1 and g == ch2: 
            return genetree
        elif g == ch1: distance2 = g.get_distance(ch2)
        elif g == ch2: distance1 = g.get_distance(ch1)
        else:
            distance1 = g.get_distance(ch1)
            distance2 = g.get_distance(ch2)
        #This has to be inefficient: looks to be copied from paper's definition but idk if this actually saves time
    
        path1 = []
        path2 = []
        losses1 = []
        losses2 = []

        #mostly the same as lines 111 to 157
        if distance1 != 0:
            curr = ch1
            while curr != g:
                    path1.append(curr.up)
                    curr = curr.up
                    #I think this will include g in path1: one child of g might have to be lost
        if path1:
            for p in path1:
                losses1 += [(p,c) for c in p.get_children() if c != ch1 and c != ch2 and c not in path1]
            if losses1: 
                losses1.sort(key=sort_by_node)
            prev = None
            for pair in losses1:
                p = pair[0]
                l = pair[1]
                if prev == None:
                    desc = children[0]
                    prev = insert_loss(node,desc,p.name,l.name)
                else: prev = insert_loss(node,prev,p.name,l.name)
        if distance2 != 0:
            curr = ch2
            while curr != g:
                    path2.append(curr.up)
                    curr = curr.up
        if path2:
            for p in path2:
                losses2 += [(p,c) for c in p.get_children() if c != ch1 and c != ch2 and c not in path2]
            if losses2: 
                losses2.sort(key=sort_by_node)
            prev = None
            for pair in losses2:
                p = pair[0]
                l = pair[1]
                if prev == None:
                    desc = children[1]
                    prev = insert_loss(node,desc,p.name,l.name)
                else: prev = insert_loss(node,prev,p.name,l.name)
    
    def loss_transfer(node,line):
        global speciestree
        global genetree

        donor_edge = line[(line.find('Mapping --> ') + 12):line.find(', R')]        #idk why it's "edge"
        donor_species = speciestree&"%s" %(donor_edge)
        transfer_edge =  line[(line.find('Recipient -->')+14):len(line)]
        transfer_species = speciestree&"%s" %(transfer_edge)

        children = node.get_children()
        ch1 = speciestree&"%s" %(children[0].mapping)
        ch2 = speciestree&"%s" %(children[1].mapping)

        distance1 = 0
        distance2 = 0

        path1 = []
        path2 = []
        losses1 = []
        losses2 = []

        if ch1 != donor_species or ch1 != transfer_species:
            curr = ch1
            if ch1.get_distance(donor_species) < ch1.get_distance(transfer_species):
                distance1 = ch1.get_distance(donor_species)
                while curr != donor_species:
                    path1.append(curr.up)
                    curr = curr.up
            else: 
                distance1 = ch1.get_distance(transfer_species)
                while curr != transfer_species:
                    path1.append(curr.up)
                    curr = curr.up
        #same as above
        if ch2 != donor_species or ch2 != transfer_species:
            curr = ch2
            if ch2.get_distance(donor_species) < ch2.get_distance(transfer_species):
                distance2 = ch2.get_distance(donor_species)
                while curr != donor_species:
                    path2.append(curr.up)
                    curr = curr.up
            else: 
                distance2 = ch2.get_distance(transfer_species)
                while curr != transfer_species:
                    path2.append(curr.up)
                    curr = curr.up
        
        #same as before
        if path1:
            for p in path1:
                losses1 += [(p,c) for c in p.get_children() if c != ch1 and c != ch2 and c not in path1]
            if losses1: 
                losses1.sort(key=sort_by_node)
            prev = None
            for pair in losses1:
                p = pair[0]
                l = pair[1]
                if prev == None:
                    prev = insert_loss(node,children[0],p.name,l.name)
                else: prev = insert_loss(node,prev,p.name,l.name)
        if path2:
            for p in path2:
                losses2 += [(p,c) for c in p.get_children() if c != ch1 and c != ch2 and c not in path2]
            if losses2: 
                losses2.sort(key=sort_by_node)
            prev = None
            for pair in losses2:
                p = pair[0]
                l = pair[1]
                if prev == None:
                    prev = insert_loss(node,children[1],p.name,l.name)
                else: prev = insert_loss(node,prev,p.name,l.name)
   
    def insert_loss(ancestor,descendant,mapping_node,mapping_loss):
        global genetree
        global speciestree
        global loss_count

        try:
            ancestor = genetree&"%s" %(ancestor.name)
            descendant = genetree&"%s" %(descendant.name)
            #converts to gene tree
        except:                                         #debugging
            print(306)      #ancestor.name,descendant.name
        if ancestor == descendant.up:
            """
                This section inserts a dummy node into the tree to 
                simulate a loss. The new node will have the same name 
                as its gene tree ancestor, but (possibly) a different 
                species tree mapping. If you want to change it, change
                name=ancestor.name ---> name="preferred name"
            """
            new_node = ancestor.add_child(name="INS_" + ancestor.name)

        elif ancestor == descendant:
            new_node = Tree()           #?
            
       
        else:
            """
                This section inserts a dummy node into the tree to 
                simulate a loss. The new node will have the same name 
                as its gene tree ancestor, but (possibly) a different 
                species tree mapping. If you want to change it, change
                name=ancestor.name ---> name="preferred name"
            """
            new_node = descendant.add_sister(name=ancestor.name)
        
        descendant.detach()
        new_node.add_child(descendant)

        new_node.add_features(mapping=mapping_node, type='speciation')

        loss = new_node.add_child(name="LOSS")
        loss.add_features(type="loss", mapping=mapping_loss) 


        if len(new_node.get_children()) != 2:               #debugging
            print(341)          #new_node
            exit()

        if ancestor == descendant:
            genetree = new_node

        return new_node

    global genetree
    genetree = gene_tree
    global speciestree
    speciestree = species_tree

    for line in reclines :
        
        if 'Transfer' in line :
            name = line[0:line.find(' =')]
            try:
                node = genetree&"%s" %(name)
            except:
                print(361)          #name
                exit()
            transfer(node,line)
            loss_transfer(node,line)

        elif 'Duplication' in line:
            name = line[0:line.find(' =')]
            try:
                node = genetree&"%s" %(name)
            except:
                print(371)          #name
                exit()
            duplication(node,line)
            loss_duplication(node)

        elif 'Speciation' in line :
            name = line[0:line.find(' =')]
            try: 
                node = genetree&"%s" %(name)
            except:
                print(381)          #name
                exit()
            speciation(node,line)
            loss_speciation(node)

        else :
            name = line[0:line.find(': ')]
            try:
                node = genetree&"%s" %(name)
            except:
                print(391)          #name
                exit()
            leaf(node,line)

    return genetree

def build_genetree_XML(genetree):
    def write_snippet(node,num_tabs):

        snippet = "  "*num_tabs + "<clade>\n"
        num_tabs += 1

        snippet += "  "*num_tabs + "<name>%s</name>\n" %(node.name)

        tab = "  "*num_tabs
        eventsRec = "%s<eventsRec>\n" %(tab)

        try:
            back = node.transferBack
            eventsRec += "%s<transferBack destinationSpecies=\"%s\"></transferBack>\n" %(tab,back)
        except:
            pass

        if node.type == 'speciation':
            eventsRec += "%s<speciation speciesLocation=\"%s\"></speciation>\n" %(tab,node.mapping)
        elif node.type == 'duplication':
            eventsRec += "%s<duplication speciesLocation=\"%s\"></duplication>\n" %(tab,node.mapping)
        elif node.type == 'branchingOut':
            eventsRec += "%s<branchingOut speciesLocation=\"%s\"></branchingOut>\n" %(tab,node.mapping)
        elif node.type == 'transferBack':
            eventsRec += "%s<transferBack destinationSpecies=\"%s\"></transferBack>\n%s<leaf speciesLocation=%s></leaf>\n" %(tab,node.transferBack,tab,node.mapping)
        elif node.type == 'leaf':
            eventsRec += "%s<leaf speciesLocation=\"%s\"></leaf>\n" %(tab,node.mapping)
        elif node.type =='loss':
            eventsRec += "%s<loss speciesLocation=\"%s\"></loss>\n" %(tab,node.mapping)

        eventsRec += "%s</eventsRec>\n" %(tab)
        
        snippet += eventsRec

        if node.children:
            try:
                child1, child2 = write_snippet(node.children[0],num_tabs), write_snippet(node.children[1],num_tabs)
            except IndexError:
                print('432: This child: ',node.name,'only has one child.')
                exit()
            snippet = snippet + child1 + child2
        
        num_tabs -= 1
        snippet += "  "*num_tabs + "</clade>\n"

        return snippet

    beginning = '<recGeneTree>\n<phylogeny rooted="true">\n'
    end = "</phylogeny>\n</recGeneTree>\n<recPhylo>"
    middle = write_snippet(genetree,1)
    recphylo = beginning + middle + end
    return recphylo

def build_speciestree_XML(speciestree):
    def write_snippet(node,num_tabs):
        snippet = "  "*num_tabs + "<clade>\n"
        num_tabs += 1

        snippet += "  "*num_tabs + "<name>%s</name>\n" %(node.name)

        if node.children:
            try:
                child1, child2 = write_snippet(node.children[0],num_tabs), write_snippet(node.children[1],num_tabs)
            except IndexError:
                print('458: This child: ',node.name,'only has one child.')
                exit()
            snippet = snippet + child1 + child2
        
        num_tabs -= 1
        snippet += "  "*num_tabs + "</clade>\n"

        return snippet

    beginning = "<recPhylo>\n<spTree>\n<phylogeny>\n"
    end = "</phylogeny>\n</spTree>\n"
    middle = write_snippet(speciestree,1)
    spXML = beginning + middle + end
    return spXML

def main(ranger_output,args):
    def findRec(lines) :
        start = lines.index('Reconciliation:') + 1
        end = 0
        for line in lines :
            if 'The minimum reconciliation cost' in line: end = lines.index(line) - 1
        return start,end
    def findSpTree(lines) :
        index = lines.index('Species Tree: ') + 1
        return lines[index]
    def findGeneTree(lines) :
        index = lines.index('Gene Tree: ') + 1
        return lines[index]
    def fixstring(strng):
        strng = re.sub(r'(m\d+)(H\d+_\d+)',r'\1',strng)
        return strng

    s, e = findRec(ranger_output)[0], findRec(ranger_output)[1]
    recLines = ranger_output[s:e]
    
    species_str = findSpTree(ranger_output)
    species_tree = Tree(species_str,format=8)
    gene_str = findGeneTree(ranger_output)
    gene_str = fixstring(gene_str)
    gene_tree = Tree(gene_str,format=8)
   
    complete_gene_tree = augment_genetree(species_tree,gene_tree,recLines)
    geneXML = build_genetree_XML(complete_gene_tree)
   
    if args.include_species:
        speciesXML = build_speciestree_XML(species_tree)
    if args.render:
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.layout_fn = my_layout
        complete_gene_tree.render(args.render_output,tree_style=ts)
    
    if args.output_file:
        out = open(args.output_file, "w+")
        if args.include_species: out.write(speciesXML)
        out.write(geneXML)
    else:
        if args.include_species: print(speciesXML)
        print(geneXML)

def my_layout(node):
    global args

    leaf_style = NodeStyle()
    leaf_style["size"] = 8
    leaf_style["shape"] = 'circle'
    leaf_style["fgcolor"] = '#000000'

    loss_style = NodeStyle()
    loss_style["size"] = 8
    loss_style["shape"] = 'circle'
    loss_style["fgcolor"] = '#FF0000'

    internal_style = NodeStyle()
    internal_style["size"] = 5
    internal_style["vt_line_width"] = 2
    internal_style["hz_line_width"] = 2
    internal_style["shape"] = 'square'

    if node.is_leaf():
        if node.name == 'LOSS': node.img_style = loss_style
        else: node.img_style = leaf_style
        F = TextFace(node.name,fstyle="italic",fsize=12)
        add_face_to_node(F,node,column=0,position="branch-right")
    elif args.label_with_name:
        node.img_style = internal_style
        F = TextFace(node.name,fstyle="normal",fsize=10)
        F.margin_top = F.margin_right = F.margin_left = 5
        add_face_to_node(F,node,column=0,position="branch-bottom")
    elif args.label_with_mappings:
        node.img_style = internal_style
        F = TextFace(node.mapping,fstyle="normal",fsize=10)
        F.margin_top = F.margin_right = F.margin_left = 5
        add_face_to_node(F,node,column=0,position="branch-bottom")
    else:
        node.img_style = internal_style

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A program to change the RANGER-DTL output to RecPhyloXML format.")
    group = parser.add_mutually_exclusive_group()
    
    parser.add_argument('-i','--input', action='store',dest='ranger_output_path', help='Specify path to RANGER-DTL output.')
    parser.add_argument('-o','--output',action='store',dest='output_file',help="Specify path to output file. Default is STDOUT.")
    parser.add_argument('--species', action='store_false',default=True,dest='include_species',help="If flag is used, species tree phylogeny is NOT included in final XML output. (Default is to include species tree)")
    parser.add_argument('-r','--render',action='store_true',default=False,dest='render',help='If flag is used, program will render the augmented genetree (includes losses). Must specify path with --render-output. Default render does NOT include labels for internal nodes.')
    parser.add_argument('--render-output',action='store',dest='render_output', help='Specify path to save rendered output. .svg type recommended.')
    
    group.add_argument('--label-map', action='store_true',default=False,dest='label_with_mappings',help="Label internal nodes of rendered gene tree with species tree mappings. Cannot be used with --label-name. Default is FALSE.")
    group.add_argument('--label-name', action='store_true',default=False,dest='label_with_name',help="Label internal nodes of rendered gene tree with names given by RANGER-DTL. Cannot be used with --label-map. Default is FALSE.")

    global args
    args = parser.parse_args()

    if not args.ranger_output_path:
        sys.exit("Must specify path to RANGER-DTL output file.")
    elif args.render and not args.render_output:
        sys.exit("Must specify path for render. Use --render-output.")

    with open(args.ranger_output_path,'r+') as file: 
        full_output = file.read().split('\n')
        main(full_output,args)
